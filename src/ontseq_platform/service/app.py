"""HTTP transport for the local browser and desktop interface. Computes nothing.

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
import os
import re
import sys
import threading
import traceback
import webbrowser
from dataclasses import dataclass, field
from datetime import UTC, datetime
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BufferedIOBase
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

from .. import __version__
from ..io import load_model
from ..models import (
    AmlKnowledgeLock,
    AnalysisIntent,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CuteSvPolicy,
    GenomeBuild,
    InputKind,
    InputSpec,
    IntervalResourceLock,
    PrivacySpec,
    QCPolicy,
    ReferenceLock,
    SampleManifest,
    SnifflesPolicy,
    SvConsensusPolicy,
    SvEvidencePolicy,
)
from ..pipeline.components import RunComponents
from ..pipeline.review import Decision, ReviewError
from ..pipeline.runner import RunConfiguration, run_pipeline
from ..profile_analysis import (
    AnalyzeSettings,
    ProfileRuntimeSettings,
    build_profile_run_configuration,
)
from ..resource_registry import ResourceRegistry
from ..review import inspect as inspect_review
from ..review import record as record_review
from ..status import scan as scan_envelopes
from ..target_coverage import TargetCoveragePolicy
from .guard import (
    TOKEN_HEADER,
    GuardError,
    host_is_loopback,
    new_token,
    origin_is_loopback,
    resolve_bam_index,
    resolve_envelope,
    resolve_within,
    token_matches,
    windows_to_wsl,
    wsl_to_windows,
)

RUNTIME_GIT_COMMIT = Path(sys.prefix) / "share" / "ontseq" / "git-commit.txt"

PAGE = Path(__file__).with_name("ONTSeq.html")

#: Extensions the file browser offers. A picker that lists everything invites a path that
#: was never meant to be an input.
INPUT_SUFFIXES = frozenset({".bam"})

#: Small JSON control messages only. Bound both fixed-length and chunked requests so a
#: local client cannot make the service buffer an arbitrary amount of data.
MAX_REQUEST_BODY_BYTES = 1024 * 1024
MAX_CHUNK_LINE_BYTES = 8192

#: Directories a by-name search may descend into before it gives up. A run directory holds
#: hundreds of gigabytes; a search with no bound would appear to hang, and a bound that is
#: hit must be reported rather than passed off as "not found".
SEARCH_DIRECTORY_LIMIT = 20_000

ACTIVE_GRCH38_PROFILE_IDS = ("AML_LCWGS_GRCh38", "AML_AS_111_GRCh38")


def _read_chunked_body(stream: BufferedIOBase) -> bytes:
    """Decode one bounded HTTP/1.1 chunked request body."""
    chunks: list[bytes] = []
    total = 0
    while True:
        line = stream.readline(MAX_CHUNK_LINE_BYTES + 1)
        if not line.endswith(b"\r\n") or len(line) > MAX_CHUNK_LINE_BYTES:
            raise ValueError("invalid chunk header")
        size_text = line[:-2].split(b";", 1)[0].strip()
        try:
            size = int(size_text, 16)
        except ValueError as error:
            raise ValueError("invalid chunk size") from error
        if size < 0:
            raise ValueError("invalid chunk size")
        if size == 0:
            trailer_bytes = 0
            while True:
                trailer = stream.readline(MAX_CHUNK_LINE_BYTES + 1)
                if trailer == b"\r\n":
                    return b"".join(chunks)
                if not trailer.endswith(b"\r\n") or len(trailer) > MAX_CHUNK_LINE_BYTES:
                    raise ValueError("invalid chunk trailer")
                trailer_bytes += len(trailer)
                if trailer_bytes > MAX_REQUEST_BODY_BYTES:
                    raise ValueError("chunk trailers are too large")
        total += size
        if total > MAX_REQUEST_BODY_BYTES:
            raise ValueError("request body is too large")
        chunk = stream.read(size)
        if len(chunk) != size or stream.read(2) != b"\r\n":
            raise ValueError("truncated chunked request body")
        chunks.append(chunk)


def _bam_is_indexed(path: Path) -> bool:
    try:
        resolve_bam_index(path)
    except GuardError:
        return False
    return True


@dataclass
class ServiceConfig:
    """What the service was started with. Nothing here is taken from a request."""

    reference_lock: Path | None
    output_dir: Path
    allowed_roots: list[Path]
    qc_policy: Path
    sniffles_policy: Path
    target_coverage_policy: Path
    components: RunComponents | None = None
    cutesv_policy: Path | None = None
    sv_consensus_policy: Path | None = None
    sv_evidence_policy: Path | None = None
    reference_fasta: Path | None = None
    gene_annotation: tuple[Path, IntervalResourceLock] | None = None
    cytoband_annotation: tuple[Path, IntervalResourceLock] | None = None
    sv_context_resources: tuple[tuple[Path, IntervalResourceLock], ...] = ()
    aml_knowledge: tuple[Path, AmlKnowledgeLock] | None = None
    sv_minimum_mean_depth: float = 10.0
    cutesv_executable: str = "cuteSV"
    resource_root: Path | None = None
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
    profile: str | None = None
    detected_genome_build: str | None = None

    def snapshot(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "sample_id": self.sample_id,
            "state": self.state,
            "detail": self.detail,
            "stages": self.stages,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "profile": self.profile,
            "detected_genome_build": self.detected_genome_build,
        }


class JobRejected(Exception):
    """Raised when a run cannot be registered. The message is safe to show the operator."""


class Jobs:
    """The runs this process has started. One at a time, as the envelope lock requires."""

    def __init__(self) -> None:
        self._jobs: dict[str, RunJob] = {}
        self._lock = threading.Lock()

    def claim(self, job: RunJob) -> None:
        """Register *job* only if nothing else is running, atomically.

        Asking ``running()`` and then calling ``add()`` reads correctly and is wrong: the
        service is a ``ThreadingHTTPServer``, so two POSTs to ``/api/runs`` arriving
        together can both see an idle service and both start a pipeline. The envelope lock
        does not catch that — it guards one envelope, and two runs with different run ids
        take different envelopes. What follows is two pipelines on one workstation, each
        sized for the whole machine, and a page that shows one of them.

        Registering a second run under a run id already in flight would also drop the first
        job from the table while its thread kept going, leaving a run nobody can observe.
        Both cases are refused here, where the decision is made under the lock.
        """
        with self._lock:
            for existing in self._jobs.values():
                if existing.state == "running":
                    raise JobRejected(
                        f"an analysis is already running ({existing.run_id}/"
                        f"{existing.sample_id}); wait for it to finish"
                    )
            if job.run_id in self._jobs:
                raise JobRejected(
                    f"run id {job.run_id!r} has already been used in this session; "
                    "choose a different run id"
                )
            self._jobs[job.run_id] = job

    def get(self, run_id: str) -> RunJob | None:
        with self._lock:
            return self._jobs.get(run_id)

    def running(self) -> bool:
        with self._lock:
            return any(job.state == "running" for job in self._jobs.values())


def _stage_view(record: Any) -> dict[str, Any]:
    """One stage, with its outcome kept as one of four values rather than two."""
    return {
        "stage": record.stage.value,
        "title": record.title,
        "status": record.status.value,
        "reason": record.reason,
        "required": record.required,
        "verification": record.verification.value,
        "resumed": record.resumed,
        "duration_seconds": record.duration_seconds,
    }


def _job_detail(report: Any) -> str:
    """Prefer concrete failed-stage reasons over the graph's generic failure summary."""

    if report.passed:
        return str(report.verdict_reason)
    failures = [record for record in report.stages if record.status.value == "FAILED"]
    if not failures:
        return str(report.verdict_reason)
    return " | ".join(f"{record.title}: {record.reason}" for record in failures[:3])


def _runtime_git_commit() -> str:
    """Read the exact commit embedded in the packed runtime, failing honestly if absent."""

    try:
        value = RUNTIME_GIT_COMMIT.read_text(encoding="ascii").strip().lower()
    except OSError:
        return "UNKNOWN"
    return value if re.fullmatch(r"[0-9a-f]{40}", value) else "UNKNOWN"


def _build_manifest(
    payload: dict[str, Any],
    *,
    reference_id: str,
    allowed_roots: list[Path],
) -> SampleManifest:
    """Turn what the page sent into the manifest contract, refusing anything unstated."""
    raw_path = str(payload.get("bam", "")).strip()
    if not raw_path:
        raise GuardError("no BAM was selected")
    bam = windows_to_wsl(raw_path)
    resolve_within(bam, allowed_roots)
    index = str(resolve_bam_index(bam))
    resolve_within(index, allowed_roots)

    mode = AssayMode(str(payload.get("assay", "")))
    genome_build = GenomeBuild(str(payload.get("genome_build", "")))
    target_bed = str(payload.get("target_bed", "")).strip() or None
    target_bed_version = str(payload.get("target_bed_version", "")).strip() or None
    target_bed_wsl = windows_to_wsl(target_bed) if target_bed else None
    if target_bed_wsl is not None:
        resolve_within(target_bed_wsl, allowed_roots)

    modules = [
        AnalysisModule.QC,
        AnalysisModule.CNV,
        AnalysisModule.SV,
        AnalysisModule.ISCN,
        AnalysisModule.REPORT,
    ]

    return SampleManifest(
        sample_id=str(payload.get("sample_id", "")).strip(),
        run_id=str(payload.get("run_id", "")).strip(),
        input=InputSpec(kind=InputKind.ALIGNED_BAM, path=bam, index_path=index),
        assay=AssaySpec(
            mode=mode,
            genome_build=genome_build,
            reference_id=reference_id,
            target_bed=target_bed_wsl,
            target_bed_version=target_bed_version,
        ),
        analysis=AnalysisSpec(
            profile=mode.value,
            modules=modules,
            intent=AnalysisIntent.SOMATIC,
        ),
        privacy=PrivacySpec(),
    )


def _execute(config: ServiceConfig, manifest: SampleManifest, job: RunJob) -> None:
    """Run the pipeline and mirror its report into the job. Failures land in the job."""
    try:
        if config.reference_lock is None:
            raise ValueError("legacy service runs require an explicit reference lock")
        run_config = RunConfiguration(
            manifest=manifest,
            reference_lock=load_model(config.reference_lock, ReferenceLock),
            output_base=config.output_dir,
            run_id=manifest.run_id,
            pipeline_version=__version__,
            git_commit=_runtime_git_commit(),
            qc_policy=load_model(config.qc_policy, QCPolicy),
            sniffles_policy=load_model(config.sniffles_policy, SnifflesPolicy),
            cutesv_policy=(
                load_model(config.cutesv_policy, CuteSvPolicy)
                if config.cutesv_policy is not None
                else None
            ),
            sv_consensus_policy=(
                load_model(config.sv_consensus_policy, SvConsensusPolicy)
                if config.sv_consensus_policy is not None
                else None
            ),
            sv_evidence_policy=(
                load_model(config.sv_evidence_policy, SvEvidencePolicy)
                if config.sv_evidence_policy is not None
                else None
            ),
            reference_fasta=config.reference_fasta,
            gene_annotation=config.gene_annotation,
            cytoband_annotation=config.cytoband_annotation,
            sv_context_resources=config.sv_context_resources,
            aml_knowledge=config.aml_knowledge,
            sv_minimum_mean_depth=config.sv_minimum_mean_depth,
            target_coverage_policy=load_model(config.target_coverage_policy, TargetCoveragePolicy),
            components=config.components,
            threads=config.threads,
            executables={"cutesv": config.cutesv_executable},
        )
        report, _bundle = run_pipeline(run_config)
        job.stages = [_stage_view(record) for record in report.stages]
        job.state = "passed" if report.passed else "failed"
        job.detail = _job_detail(report)
    except Exception as error:  # noqa: BLE001 - the page must see every failure, typed or not
        job.state = "error"
        job.detail = f"{type(error).__name__}: {error}"
        traceback.print_exc()
    finally:
        job.finished_at = datetime.now(UTC).isoformat()


def _execute_resolved(run_config: RunConfiguration, job: RunJob) -> None:
    """Execute a profile-resolved configuration without rebuilding any resource choice."""

    try:
        report, _bundle = run_pipeline(run_config)
        job.stages = [_stage_view(record) for record in report.stages]
        job.state = "passed" if report.passed else "failed"
        job.detail = _job_detail(report)
    except Exception as error:  # noqa: BLE001 - the page must see every failure, typed or not
        job.state = "error"
        job.detail = f"{type(error).__name__}: {error}"
        traceback.print_exc()
    finally:
        job.finished_at = datetime.now(UTC).isoformat()


def _build_profile_configuration(
    payload: dict[str, Any],
    config: ServiceConfig,
) -> RunConfiguration:
    profile_id = str(payload.get("profile", "")).strip()
    if profile_id not in ACTIVE_GRCH38_PROFILE_IDS:
        raise ValueError(f"unsupported GRCh38 analysis profile: {profile_id!r}")
    if config.resource_root is None:
        raise ValueError("profile runs require the service to start with --resource-root")
    registry = ResourceRegistry(config.resource_root, active_build=GenomeBuild.GRCH38)
    if profile_id not in registry.profiles:
        raise ValueError(f"GRCh38 analysis profile {profile_id!r} is not active")
    raw_path = str(payload.get("bam", "")).strip()
    if not raw_path:
        raise GuardError("no BAM was selected")
    bam = Path(windows_to_wsl(raw_path))
    resolve_within(bam, config.allowed_roots)
    index = resolve_bam_index(bam)
    resolve_within(index, config.allowed_roots)
    legacy_build = str(payload.get("genome_build", "")).strip()
    if legacy_build and legacy_build != GenomeBuild.GRCH38.value:
        raise ValueError("profile runs are strictly GRCh38; cross-build fallback is prohibited")
    expected_assay = (
        AssayMode.ADAPTIVE_SAMPLING.value
        if profile_id == "AML_AS_111_GRCh38"
        else AssayMode.LOW_COVERAGE_WGS.value
    )
    legacy_assay = str(payload.get("assay", "")).strip()
    if legacy_assay and legacy_assay != expected_assay:
        raise ValueError(
            f"profile {profile_id!r} requires assay {expected_assay!r}, not {legacy_assay!r}"
        )
    return build_profile_run_configuration(
        AnalyzeSettings(
            bam=bam,
            profile_id=profile_id,
            resource_root=config.resource_root,
            output_dir=config.output_dir,
            sample_id=str(payload.get("sample_id", "")).strip() or None,
            run_id=str(payload.get("run_id", "")).strip() or None,
            pipeline_version=__version__,
            git_commit=_runtime_git_commit(),
            threads=config.threads,
            executables={"cutesv": config.cutesv_executable},
            runtime_settings=ProfileRuntimeSettings(
                qc_policy=load_model(config.qc_policy, QCPolicy),
                sniffles_policy=load_model(config.sniffles_policy, SnifflesPolicy),
                cutesv_policy=(
                    load_model(config.cutesv_policy, CuteSvPolicy)
                    if config.cutesv_policy is not None
                    else None
                ),
                sv_consensus_policy=(
                    load_model(config.sv_consensus_policy, SvConsensusPolicy)
                    if config.sv_consensus_policy is not None
                    else None
                ),
                sv_evidence_policy=(
                    load_model(config.sv_evidence_policy, SvEvidencePolicy)
                    if config.sv_evidence_policy is not None
                    else None
                ),
                target_coverage_policy=load_model(
                    config.target_coverage_policy, TargetCoveragePolicy
                ),
                sv_minimum_mean_depth=config.sv_minimum_mean_depth,
                components=config.components,
            ),
            # Starting a run checks presence, declared sizes and checksum pins.  Reading
            # every byte of the multi-GB reference here would make a POST appear hung;
            # operators request that deeper audit explicitly with ``references validate``.
            verify_resource_checksums=False,
        )
    )


def _resolvable_profile_ids(config: ServiceConfig) -> list[str]:
    """Return only published profiles whose complete pinned context is locally ready."""

    if config.resource_root is None:
        return []
    try:
        registry = ResourceRegistry(config.resource_root, active_build=GenomeBuild.GRCH38)
    except (OSError, ValueError):
        return []
    ready: list[str] = []
    for profile_id in ACTIVE_GRCH38_PROFILE_IDS:
        if profile_id not in registry.profiles:
            continue
        try:
            registry.resolve_profile(profile_id, verify_files=False)
        except (KeyError, OSError, ValueError):
            continue
        ready.append(profile_id)
    return ready


def make_handler(config: ServiceConfig, jobs: Jobs) -> type[BaseHTTPRequestHandler]:
    """Build the request handler bound to one service configuration."""

    class Handler(BaseHTTPRequestHandler):
        server_version = f"ONTSeq/{__version__}"
        protocol_version = "HTTP/1.1"

        def log_message(self, fmt: str, *args: Any) -> None:
            print(f"[{self.log_date_time_string()}] {fmt % args}")

        def _send(self, status: HTTPStatus, body: bytes, content_type: str) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self.wfile.write(body)

        def _json(self, status: HTTPStatus, payload: dict[str, Any]) -> None:
            self._send(status, json.dumps(payload).encode("utf-8"), "application/json")

        def _refuse(self, status: HTTPStatus, reason: str) -> None:
            self._json(status, {"error": reason})

        def _authorised(self) -> bool:
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

        def _request_body(self) -> bytes:
            """Read a bounded fixed-length or HTTP/1.1 chunked request body."""
            transfer_encoding = self.headers.get("Transfer-Encoding")
            content_length = self.headers.get("Content-Length")
            if transfer_encoding is not None:
                if content_length is not None:
                    raise ValueError("ambiguous request framing")
                if transfer_encoding.strip().lower() != "chunked":
                    raise ValueError("unsupported transfer encoding")
                return _read_chunked_body(self.rfile)
            if content_length is None:
                return b""
            try:
                length = int(content_length)
            except ValueError as error:
                raise ValueError("invalid content length") from error
            if length < 0:
                raise ValueError("invalid content length")
            if length > MAX_REQUEST_BODY_BYTES:
                raise ValueError("request body is too large")
            body = self.rfile.read(length)
            if len(body) != length:
                raise ValueError("truncated request body")
            return body

        def do_GET(self) -> None:  # noqa: N802 - name fixed by BaseHTTPRequestHandler
            route = urlparse(self.path)
            if route.path in {"/", "/index.html", "/ONTSeq.html"}:
                # The page carries the session token, so it gets the same rebinding check
                # the API gets. Without it a page on an attacker-controlled name resolved
                # to 127.0.0.1 is same-origin with this service and can simply read the
                # token out of the response. The API's Host check still refuses the token
                # afterwards; not handing it over in the first place is cheaper than
                # relying on that being the only place it matters.
                if not host_is_loopback(self.headers.get("Host"), port=config.port):
                    self._refuse(HTTPStatus.FORBIDDEN, "request did not come from loopback")
                    return
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
            elif route.path == "/api/locate":
                self._locate(query.get("name", [""])[0])
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
            try:
                payload = json.loads(self._request_body() or b"{}")
            except ValueError as error:
                self.close_connection = True
                self._refuse(HTTPStatus.BAD_REQUEST, str(error) or "request body was not JSON")
                return
            if not isinstance(payload, dict):
                self._refuse(HTTPStatus.BAD_REQUEST, "request body was not a JSON object")
                return
            if path == "/api/runs":
                self._start_run(payload)
                return
            parts = path[len("/api/review/") :].split("/")
            if len(parts) != 2 or not all(parts):
                self._refuse(HTTPStatus.BAD_REQUEST, "expected /api/review/<run>/<sample>")
                return
            self._review(parts[0], parts[1], payload)

        def _serve_page(self) -> None:
            html = PAGE.read_text(encoding="utf-8")
            html = html.replace("__ONTSEQ_TOKEN__", config.token)
            self._send(HTTPStatus.OK, html.encode("utf-8"), "text/html; charset=utf-8")

        def _config_view(self) -> dict[str, Any]:
            not_wired = ["basecalling — not needed for an already aligned BAM"]
            return {
                "version": __version__,
                "roots": [
                    {"posix": str(root), "display": wsl_to_windows(str(root))}
                    for root in config.allowed_roots
                ],
                "output_dir": str(config.output_dir),
                "busy": jobs.running(),
                "not_wired": not_wired,
                "resource_root": None
                if config.resource_root is None
                else str(config.resource_root),
                "profiles": _resolvable_profile_ids(config),
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
            try:
                listing = sorted(target.iterdir(), key=lambda entry: entry.name.lower())
            except OSError as error:
                self._refuse(HTTPStatus.BAD_REQUEST, f"cannot list {target}: {error}")
                return
            entries = []
            skipped = 0
            for item in listing:
                if item.name.startswith("."):
                    continue
                try:
                    is_dir = item.is_dir()
                    if not is_dir and item.suffix.lower() not in INPUT_SUFFIXES:
                        continue
                    size = None if is_dir else item.stat().st_size
                except OSError:
                    # A dangling symlink or an entry this user may not stat is not a
                    # reason to fail the whole listing. Letting it raise here aborts the
                    # request without a response, and the page shows a network error for
                    # a directory whose other entries are perfectly usable.
                    skipped += 1
                    continue
                entries.append(
                    {
                        "name": item.name,
                        "posix": str(item),
                        "display": wsl_to_windows(str(item)),
                        "directory": is_dir,
                        "size_bytes": size,
                        "indexed": is_dir or _bam_is_indexed(item),
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
                    # Reported rather than swallowed: a BAM missing from a listing is
                    # something the operator has to be able to see, not guess at.
                    "unreadable_entries": skipped,
                },
            )

        def _start_run(self, payload: dict[str, Any]) -> None:
            profile_id = str(payload.get("profile", "")).strip()
            try:
                resolved_config: RunConfiguration | None = None
                if profile_id:
                    resolved_config = _build_profile_configuration(payload, config)
                    manifest = resolved_config.manifest
                else:
                    if config.reference_lock is None:
                        raise ValueError(
                            "the service has no legacy reference lock; select an installed profile"
                        )
                    lock = load_model(config.reference_lock, ReferenceLock)
                    manifest = _build_manifest(
                        payload,
                        reference_id=lock.reference_id,
                        allowed_roots=config.allowed_roots,
                    )
                    if manifest.assay.genome_build != lock.genome_build:
                        raise ValueError(
                            "selected genome build does not match the configured reference lock"
                        )
            except GuardError as error:
                self._refuse(HTTPStatus.FORBIDDEN, str(error))
                return
            except (KeyError, OSError, ValueError) as error:
                self._refuse(HTTPStatus.BAD_REQUEST, str(error))
                return

            job = RunJob(
                run_id=manifest.run_id,
                sample_id=manifest.sample_id,
                envelope=config.output_dir / manifest.run_id / manifest.sample_id,
                profile=profile_id or manifest.analysis.profile,
                detected_genome_build=manifest.assay.genome_build.value,
            )
            try:
                # Claimed before the thread starts, so the slot is taken by whichever
                # request got here first rather than by whichever thread starts fastest.
                jobs.claim(job)
            except JobRejected as error:
                self._refuse(HTTPStatus.CONFLICT, str(error))
                return
            if resolved_config is not None:
                threading.Thread(
                    target=_execute_resolved,
                    args=(resolved_config, job),
                    daemon=True,
                ).start()
            else:
                threading.Thread(target=_execute, args=(config, manifest, job), daemon=True).start()
            self._json(HTTPStatus.ACCEPTED, job.snapshot())

        def _locate(self, name: str) -> None:
            """Find a BAM by its bare name inside the allowed roots without guessing."""
            wanted = Path(name.replace("\\", "/")).name
            if not wanted or wanted.startswith("."):
                self._refuse(HTTPStatus.BAD_REQUEST, "no usable file name was given")
                return
            if Path(wanted).suffix.lower() not in INPUT_SUFFIXES:
                self._refuse(HTTPStatus.BAD_REQUEST, f"not a BAM file: {wanted}")
                return

            matches: list[dict[str, Any]] = []
            visited = 0
            exhausted = False
            for root in config.allowed_roots:
                resolved_root = Path(root).expanduser().resolve()
                for directory, _subdirs, files in os.walk(resolved_root):
                    visited += 1
                    if visited > SEARCH_DIRECTORY_LIMIT:
                        exhausted = True
                        break
                    # Construct the candidate from the filename returned by the
                    # filesystem, not from the query parameter. Apart from making the
                    # data-flow boundary explicit, resolving it again prevents a BAM
                    # symlink inside an allowed root from exposing a target outside it.
                    for entry_name in files:
                        if entry_name != wanted:
                            continue
                        candidate = Path(directory) / entry_name
                        try:
                            found = resolve_within(candidate, [resolved_root])
                            if not found.is_file():
                                continue
                            size = found.stat().st_size
                        except (GuardError, OSError):
                            # Same reasoning as the browser: an entry that cannot be
                            # resolved/stat'ed is dropped from the results, not turned
                            # into a failed search across every other root.
                            continue
                        matches.append(
                            {
                                "posix": str(found),
                                "display": wsl_to_windows(str(found)),
                                "size_bytes": size,
                                "indexed": _bam_is_indexed(found),
                            }
                        )
                        break
                if exhausted:
                    break
            self._json(
                HTTPStatus.OK,
                {
                    "name": wanted,
                    "matches": matches,
                    "search_incomplete": exhausted,
                    "roots": [wsl_to_windows(str(root)) for root in config.allowed_roots],
                },
            )

        def _findings(self) -> None:
            """Every envelope this output directory holds, including failed runs."""
            try:
                statuses = scan_envelopes(config.output_dir)
            except NotADirectoryError:
                self._json(HTTPStatus.OK, {"findings": []})
                return
            findings = []
            for status in statuses:
                result = status.report
                finished = None if result is None else result.finished_at.isoformat()
                findings.append(
                    {
                        "run_id": status.run_id,
                        "sample_id": status.sample_id,
                        "state": status.state.value,
                        "detail": status.detail,
                        "review": None if status.review is None else status.review.value,
                        "review_detail": status.review_detail,
                        "unverified_stages": list(status.unverified_stages),
                        "finished_at": finished,
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
            reviewer = str(payload.get("reviewer", "")).strip()
            if not reviewer:
                self._refuse(HTTPStatus.BAD_REQUEST, "no reviewer name was given")
                return
            try:
                # Both names come straight out of the request URL. Joining them onto the
                # output directory unchecked is how "/api/review/../other" signs off an
                # envelope this service was never pointed at.
                envelope = resolve_envelope(config.output_dir, run_id, sample_id)
            except GuardError as error:
                self._refuse(HTTPStatus.FORBIDDEN, str(error))
                return
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
