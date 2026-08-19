"""HTTP transport for the local browser interface. Computes nothing.

Every analytical decision already exists: ``run_pipeline`` executes the stages, ``RunReport``
records what each one did, and the report, workbook and JSON are written by the same code
the command line uses. This module starts that work, reports its progress and hands back the
files. A defect here can lose or misrepresent a result; it cannot change one.

Two properties are deliberate and worth not undoing:

**The page is served, not opened from disk.** A double-clicked ``file://`` page has an opaque
origin, cannot be told apart from any other local page, and needs the service to allow
cross-origin requests to work at all. Serving it means the page and the service share an
origin, the token can be handed over without anyone copying it, and the checks in ``guard``
have something real to check.

**Progress carries all four stage outcomes.** A stage that did not run is not a stage that
found nothing. Collapsing them into "done" and "pending" would undo, in the one place
everybody looks, the distinction the rest of the system is built to preserve.
"""

from __future__ import annotations

import json
import threading
import traceback
import webbrowser
from dataclasses import dataclass, field
from datetime import UTC, datetime
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

from .. import __version__
from ..io import load_model
from ..models import (
    AnalysisIntent,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    PrivacySpec,
    QCPolicy,
    ReferenceLock,
    SampleManifest,
    SnifflesPolicy,
)
from ..pipeline.review import Decision, ReviewError
from ..pipeline.runner import RunConfiguration, run_pipeline
from ..review import inspect as inspect_review
from ..review import record as record_review
from ..status import scan as scan_envelopes
from .guard import (
    TOKEN_HEADER,
    GuardError,
    host_is_loopback,
    new_token,
    origin_is_loopback,
    resolve_within,
    token_matches,
    windows_to_wsl,
    wsl_to_windows,
)

PAGE = Path(__file__).with_name("ONTSeq.html")

#: Extensions the file browser offers. A picker that lists everything invites a path that
#: was never meant to be an input.
INPUT_SUFFIXES = frozenset({".bam"})


@dataclass
class ServiceConfig:
    """What the service was started with. Nothing here is taken from a request."""

    reference_lock: Path
    output_dir: Path
    allowed_roots: list[Path]
    qc_policy: Path
    sniffles_policy: Path
    host: str = "127.0.0.1"
    port: int = 8765
    threads: int = 4
    token: str = field(default_factory=new_token)


@dataclass
class RunJob:
    """One analysis, as the page sees it while it is happening."""

    run_id: str
    sample_id: str
    envelope: Path
    state: str = "running"
    detail: str = ""
    stages: list[dict[str, Any]] = field(default_factory=list)
    started_at: str = field(default_factory=lambda: datetime.now(UTC).isoformat())
    finished_at: str | None = None

    def snapshot(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "sample_id": self.sample_id,
            "state": self.state,
            "detail": self.detail,
            "stages": self.stages,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
        }


class Jobs:
    """The runs this process has started. One at a time, as the envelope lock requires."""

    def __init__(self) -> None:
        self._jobs: dict[str, RunJob] = {}
        self._lock = threading.Lock()

    def add(self, job: RunJob) -> None:
        with self._lock:
            self._jobs[job.run_id] = job

    def get(self, run_id: str) -> RunJob | None:
        with self._lock:
            return self._jobs.get(run_id)

    def running(self) -> bool:
        with self._lock:
            return any(job.state == "running" for job in self._jobs.values())


def _stage_view(record: Any) -> dict[str, Any]:
    """One stage, with its outcome kept as one of four values rather than two.

    ``not_run`` and ``no_call`` look identical to a progress bar and mean opposite things:
    nothing looked, versus something looked and found nothing.
    """
    return {
        "stage": record.stage.value,
        "title": record.title,
        # The contract's value verbatim — COMPLETED, NO_CALL, FAILED, NOT_RUN — so what the
        # page shows can be compared with provenance/run.json without translating. The page
        # lowercases it for display; nothing renames it here.
        "status": record.status.value,
        "reason": record.reason,
        "required": record.required,
        "verification": record.verification.value,
        "resumed": record.resumed,
        "duration_seconds": record.duration_seconds,
    }


def _build_manifest(payload: dict[str, Any], *, reference_id: str) -> SampleManifest:
    """Turn what the page sent into the manifest contract, refusing anything unstated."""
    raw_path = str(payload.get("bam", "")).strip()
    if not raw_path:
        raise GuardError("no BAM was selected")
    bam = windows_to_wsl(raw_path)
    index = f"{bam}.bai"

    mode = AssayMode(str(payload.get("assay", "")))
    target_bed = str(payload.get("target_bed", "")).strip() or None
    target_bed_version = str(payload.get("target_bed_version", "")).strip() or None

    return SampleManifest(
        sample_id=str(payload.get("sample_id", "")).strip(),
        run_id=str(payload.get("run_id", "")).strip(),
        input=InputSpec(kind=InputKind.ALIGNED_BAM, path=bam, index_path=index),
        assay=AssaySpec(
            mode=mode,
            genome_build=GenomeBuild(str(payload.get("genome_build", ""))),
            reference_id=reference_id,
            target_bed=windows_to_wsl(target_bed) if target_bed else None,
            target_bed_version=target_bed_version,
        ),
        analysis=AnalysisSpec(
            profile=mode.value,
            modules=[
                AnalysisModule.QC,
                AnalysisModule.SV,
                AnalysisModule.ISCN,
                AnalysisModule.REPORT,
            ],
            # Stated, never inferred. An AML workup asks a somatic question, and leaving
            # this unset would make every knowledge-base assertion report an unknown scope.
            intent=AnalysisIntent.SOMATIC,
        ),
        privacy=PrivacySpec(),
    )


def _execute(config: ServiceConfig, manifest: SampleManifest, job: RunJob) -> None:
    """Run the pipeline and mirror its report into the job. Failures land in the job."""
    try:
        run_config = RunConfiguration(
            manifest=manifest,
            reference_lock=load_model(config.reference_lock, ReferenceLock),
            output_base=config.output_dir,
            run_id=manifest.run_id,
            pipeline_version=__version__,
            git_commit="UNKNOWN",
            qc_policy=load_model(config.qc_policy, QCPolicy),
            sniffles_policy=load_model(config.sniffles_policy, SnifflesPolicy),
            threads=config.threads,
        )
        report, _bundle = run_pipeline(run_config)
        job.stages = [_stage_view(record) for record in report.stages]
        job.state = "passed" if report.passed else "failed"
        job.detail = report.verdict_reason
    except Exception as error:  # noqa: BLE001 - the page must see every failure, typed or not
        job.state = "error"
        # The class name matters to whoever reads this: the runner distinguishes a locked
        # envelope, an already-reviewed envelope and an ordinary failure, and the page
        # should not flatten them into "something went wrong".
        job.detail = f"{type(error).__name__}: {error}"
        traceback.print_exc()
    finally:
        job.finished_at = datetime.now(UTC).isoformat()


def make_handler(config: ServiceConfig, jobs: Jobs) -> type[BaseHTTPRequestHandler]:
    """Build the request handler bound to one service configuration."""

    class Handler(BaseHTTPRequestHandler):
        server_version = f"ONTSeq/{__version__}"
        protocol_version = "HTTP/1.1"

        def log_message(self, fmt: str, *args: Any) -> None:
            print(f"[{self.log_date_time_string()}] {fmt % args}")

        # -- plumbing ------------------------------------------------------------------
        def _send(self, status: HTTPStatus, body: bytes, content_type: str) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            # No caching: the page embeds a token that is only valid for this process.
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self.wfile.write(body)

        def _json(self, status: HTTPStatus, payload: dict[str, Any]) -> None:
            self._send(status, json.dumps(payload).encode("utf-8"), "application/json")

        def _refuse(self, status: HTTPStatus, reason: str) -> None:
            self._json(status, {"error": reason})

        def _authorised(self) -> bool:
            """Both checks, in this order, before anything reads a path or starts a run."""
            if not host_is_loopback(self.headers.get("Host"), port=config.port):
                self._refuse(HTTPStatus.FORBIDDEN, "request did not come from loopback")
                return False
            if not origin_is_loopback(self.headers.get("Origin"), port=config.port):
                self._refuse(HTTPStatus.FORBIDDEN, "request came from another origin")
                return False
            if not token_matches(self.headers.get(TOKEN_HEADER), config.token):
                self._refuse(HTTPStatus.UNAUTHORIZED, "missing or wrong session token")
                return False
            return True

        # -- routes --------------------------------------------------------------------
        def do_GET(self) -> None:  # noqa: N802 - name fixed by BaseHTTPRequestHandler
            route = urlparse(self.path)
            if route.path in {"/", "/index.html", "/ONTSeq.html"}:
                self._serve_page()
                return
            if not route.path.startswith("/api/"):
                self._refuse(HTTPStatus.NOT_FOUND, "no such page")
                return
            if not self._authorised():
                return
            query = parse_qs(route.query)
            if route.path == "/api/config":
                self._json(HTTPStatus.OK, self._config_view())
            elif route.path == "/api/browse":
                self._browse(query.get("path", [""])[0])
            elif route.path == "/api/findings":
                self._findings()
            elif route.path.startswith("/api/runs/"):
                self._run_status(route.path.rsplit("/", 1)[-1])
            else:
                self._refuse(HTTPStatus.NOT_FOUND, "no such endpoint")

        def do_POST(self) -> None:  # noqa: N802 - name fixed by BaseHTTPRequestHandler
            path = urlparse(self.path).path
            if path not in {"/api/runs"} and not path.startswith("/api/review/"):
                self._refuse(HTTPStatus.NOT_FOUND, "no such endpoint")
                return
            if not self._authorised():
                return
            length = int(self.headers.get("Content-Length", "0"))
            try:
                payload = json.loads(self.rfile.read(length) or b"{}")
            except ValueError:
                self._refuse(HTTPStatus.BAD_REQUEST, "request body was not JSON")
                return
            if path == "/api/runs":
                self._start_run(payload)
                return
            parts = path[len("/api/review/") :].split("/")
            if len(parts) != 2 or not all(parts):
                self._refuse(HTTPStatus.BAD_REQUEST, "expected /api/review/<run>/<sample>")
                return
            self._review(parts[0], parts[1], payload)

        # -- handlers ------------------------------------------------------------------
        def _serve_page(self) -> None:
            html = PAGE.read_text(encoding="utf-8")
            # The token reaches the page here rather than through the operator's clipboard.
            html = html.replace("__ONTSEQ_TOKEN__", config.token)
            self._send(HTTPStatus.OK, html.encode("utf-8"), "text/html; charset=utf-8")

        def _config_view(self) -> dict[str, Any]:
            return {
                "version": __version__,
                "roots": [
                    {"posix": str(root), "display": wsl_to_windows(str(root))}
                    for root in config.allowed_roots
                ],
                "output_dir": str(config.output_dir),
                "busy": jobs.running(),
                # Said once, in the place the page can show before anything is started.
                "not_wired": [
                    "cnv — no caller is wired in; the stage records NOT_RUN",
                    "basecalling — not needed, the instrument already basecalled",
                ],
            }

        def _browse(self, requested: str) -> None:
            try:
                target = (
                    resolve_within(requested, config.allowed_roots)
                    if requested
                    else Path(config.allowed_roots[0]).resolve()
                )
            except (GuardError, IndexError) as error:
                self._refuse(HTTPStatus.FORBIDDEN, str(error))
                return
            if not target.is_dir():
                self._refuse(HTTPStatus.BAD_REQUEST, f"not a directory: {target}")
                return
            entries = []
            for item in sorted(target.iterdir(), key=lambda entry: entry.name.lower()):
                if item.name.startswith("."):
                    continue
                is_dir = item.is_dir()
                if not is_dir and item.suffix.lower() not in INPUT_SUFFIXES:
                    continue
                entries.append(
                    {
                        "name": item.name,
                        "posix": str(item),
                        "display": wsl_to_windows(str(item)),
                        "directory": is_dir,
                        "size_bytes": None if is_dir else item.stat().st_size,
                        "indexed": is_dir or Path(f"{item}.bai").is_file(),
                    }
                )
            parent = str(target.parent) if target != target.parent else None
            self._json(
                HTTPStatus.OK,
                {
                    "path": str(target),
                    "display": wsl_to_windows(str(target)),
                    "parent": parent,
                    "entries": entries,
                },
            )

        def _start_run(self, payload: dict[str, Any]) -> None:
            if jobs.running():
                self._refuse(HTTPStatus.CONFLICT, "an analysis is already running")
                return
            try:
                lock = load_model(config.reference_lock, ReferenceLock)
                manifest = _build_manifest(payload, reference_id=lock.reference_id)
                resolve_within(manifest.input.path, config.allowed_roots)
            except GuardError as error:
                self._refuse(HTTPStatus.FORBIDDEN, str(error))
                return
            except (OSError, ValueError) as error:
                self._refuse(HTTPStatus.BAD_REQUEST, str(error))
                return

            job = RunJob(
                run_id=manifest.run_id,
                sample_id=manifest.sample_id,
                envelope=config.output_dir / manifest.run_id / manifest.sample_id,
            )
            jobs.add(job)
            args = (config, manifest, job)
            threading.Thread(target=_execute, args=args, daemon=True).start()
            self._json(HTTPStatus.ACCEPTED, job.snapshot())

        def _findings(self) -> None:
            """Every envelope this output directory holds, for the reviewing physician.

            Deliberately not filtered to the ones that passed. A run that failed, or one
            whose sign-off went stale because the release changed underneath it, is exactly
            what somebody needs to see — hiding it would make the list look tidier and the
            situation less true.
            """
            try:
                statuses = scan_envelopes(config.output_dir)
            except NotADirectoryError:
                self._json(HTTPStatus.OK, {"findings": []})
                return
            findings = []
            for status in statuses:
                result = status.report
                findings.append(
                    {
                        "run_id": status.run_id,
                        "sample_id": status.sample_id,
                        "state": status.state.value,
                        "detail": status.detail,
                        "review": None if status.review is None else status.review.value,
                        "review_detail": status.review_detail,
                        "unverified_stages": list(status.unverified_stages),
                        "finished_at": None
                        if result is None
                        else result.finished_at.isoformat(),
                        "not_run": []
                        if result is None
                        else [
                            record.title
                            for record in result.stages
                            if record.status.value == "NOT_RUN"
                        ],
                    }
                )
            self._json(HTTPStatus.OK, {"findings": findings})

        def _review(self, run_id: str, sample_id: str, payload: dict[str, Any]) -> None:
            """Record one judgement. The name is asserted; nothing here authenticates it.

            That is not a gap this layer introduces — the trail says the same about a
            judgement recorded from the command line. What the entry binds to is the
            checksum of the release bundle, so a judgement always names what was judged.
            """
            reviewer = str(payload.get("reviewer", "")).strip()
            if not reviewer:
                self._refuse(HTTPStatus.BAD_REQUEST, "no reviewer name was given")
                return
            envelope = config.output_dir / run_id / sample_id
            try:
                entry = record_review(
                    envelope,
                    decision=Decision(str(payload.get("decision", ""))),
                    reviewer=reviewer,
                    note=str(payload.get("note", "")).strip(),
                )
            except (NotADirectoryError, ReviewError, ValueError) as error:
                self._refuse(HTTPStatus.BAD_REQUEST, str(error))
                return
            report = inspect_review(envelope)
            self._json(
                HTTPStatus.OK,
                {
                    "decision": entry.decision.value,
                    "reviewer": entry.reviewer,
                    "state": report.state.value,
                    "detail": report.detail,
                    "release_sha256": report.release_sha256,
                },
            )

        def _run_status(self, run_id: str) -> None:
            job = jobs.get(run_id)
            if job is None:
                self._refuse(HTTPStatus.NOT_FOUND, f"no such run: {run_id}")
                return
            self._json(HTTPStatus.OK, job.snapshot())

    return Handler


def serve(config: ServiceConfig, *, open_browser: bool = True) -> None:
    """Run the service until interrupted."""
    jobs = Jobs()
    server = ThreadingHTTPServer((config.host, config.port), make_handler(config, jobs))
    url = f"http://{config.host}:{config.port}/"
    print(f"ONTSeq {__version__} — {url}", flush=True)
    print(f"  Ausgabe:   {config.output_dir}")
    for root in config.allowed_roots:
        print(f"  Freigabe:  {root}  ({wsl_to_windows(str(root))})")
    print("  Nur auf der Loopback-Schnittstelle. Beenden mit Strg+C.", flush=True)
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nbeendet")
    finally:
        server.server_close()
