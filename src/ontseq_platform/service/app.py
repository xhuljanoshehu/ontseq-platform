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

import hashlib
import json
import os
import re
import shutil
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
from ..execution import SubprocessRunner, ToolExecutionError
from ..io import load_model
from ..methylation import MethylationPolicy, MethylationReport, modkit_version
from ..methylation_probe import probe_bam_methylation
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
    PipelineResult,
    PrivacySpec,
    QCPolicy,
    ReferenceLock,
    ResolvedResourceContext,
    SampleManifest,
    SnifflesPolicy,
    SvConsensusPolicy,
    SvEvidencePolicy,
)
from ..pipeline.components import RunComponents
from ..pipeline.lock import LOCK_FILENAME
from ..pipeline.review import Decision, ReviewError
from ..pipeline.runner import (
    METHYLATION_REPORT,
    REPORT_HTML,
    REPORT_XLSX,
    RESULT_JSON,
    RUN_REPORT,
    RunConfiguration,
    run_pipeline,
)
from ..pipeline.state import RunReport
from ..profile_analysis import (
    AnalyzeSettings,
    ProfileRuntimeSettings,
    build_profile_run_configuration,
    configuration_root,
)
from ..resource_bootstrap import GRCH37_PROFILE_IDS, PROFILE_IDS
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
from .methylation_scans import (
    AvailabilityCache,
    MethylationScans,
    ScanBusyError,
    ScanNotFoundError,
)

RUNTIME_GIT_COMMIT = Path(sys.prefix) / "share" / "ontseq" / "git-commit.txt"

PAGE = Path(__file__).with_name("ONTSeq.html")
WORKSPACE_PAGE = Path(__file__).with_name("workspace.html")
MAX_RESULT_BYTES = 64 * 1024 * 1024

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

ACTIVE_GRCH38_PROFILE_IDS = PROFILE_IDS
ACTIVE_PROFILE_BUILDS = {
    **dict.fromkeys(PROFILE_IDS, GenomeBuild.GRCH38),
    **dict.fromkeys(GRCH37_PROFILE_IDS, GenomeBuild.GRCH37),
}


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
    methylation_policy: Path | None = None
    modkit_executable: str = "modkit"
    samtools_executable: str = "samtools"
    resource_root: Path | None = None
    host: str = "127.0.0.1"
    port: int = 8765
    threads: int = 4
    token: str = field(default_factory=new_token)
    instance_id: str | None = None


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

    # An editable source tree may run with an older installed tool runtime. Its embedded
    # commit identifies that runtime, not this source code; never misattribute the run.
    if (Path(__file__).resolve().parents[3] / "pyproject.toml").is_file():
        return "LOCAL_WORKTREE"
    try:
        value = RUNTIME_GIT_COMMIT.read_text(encoding="ascii").strip().lower()
    except OSError:
        return "UNKNOWN"
    return value if re.fullmatch(r"[0-9a-f]{40}", value) else "UNKNOWN"


def _include_methylation(payload: dict[str, Any]) -> bool:
    value = payload.get("include_methylation", False)
    if not isinstance(value, bool):
        raise ValueError("include_methylation must be an explicit JSON boolean")
    return value


def _local_input_path(raw_path: str) -> str:
    translated = windows_to_wsl(raw_path)
    # The Desktop runtime runs in WSL, while local service tests/development also run
    # natively on Windows. Keep native drive paths native before the same containment gate.
    return raw_path.strip().strip('"') if os.name == "nt" else translated


def _verify_methylation_runtime(policy: MethylationPolicy, executable: str) -> None:
    try:
        check = SubprocessRunner().run([executable, "--version"], timeout_seconds=10)
        if check.returncode != 0:
            raise ValueError("modkit did not return its version")
        observed = modkit_version(f"{check.stdout}\n{check.stderr}")
    except (ToolExecutionError, ValueError) as error:
        raise ValueError(
            "Methylation requires an available modkit executable with a readable version."
        ) from error
    if observed != policy.expected_version:
        raise ValueError(
            f"Methylation requires modkit {policy.expected_version}; detected {observed}."
        )


def _methylation_policy(config: ServiceConfig, assay: AssayMode) -> MethylationPolicy:
    filename = (
        "modkit-targets.technical.yaml"
        if assay == AssayMode.ADAPTIVE_SAMPLING
        else "modkit.technical.yaml"
    )
    path = config.methylation_policy or configuration_root() / "methylation" / filename
    if not path.is_file():
        raise ValueError("Methylation requires an installed versioned modkit policy.")
    return load_model(path, MethylationPolicy)


def _methylation_availability(config: ServiceConfig) -> dict[str, Any]:
    try:
        policy = _methylation_policy(config, AssayMode.LOW_COVERAGE_WGS)
    except (OSError, ValueError):
        return {
            "methylation_available": False,
            "methylation_unavailable_reason": "A versioned methylation policy is unavailable.",
            "expected_modkit_version": None,
        }
    try:
        _verify_methylation_runtime(policy, config.modkit_executable)
    except ValueError as error:
        return {
            "methylation_available": False,
            "methylation_unavailable_reason": str(error),
            "expected_modkit_version": policy.expected_version,
        }
    return {
        "methylation_available": True,
        "methylation_unavailable_reason": None,
        "expected_modkit_version": policy.expected_version,
    }


def _methylation_availability_identity(config: ServiceConfig) -> tuple[object, ...]:
    executable = shutil.which(config.modkit_executable)
    policy = config.methylation_policy or configuration_root() / "methylation/modkit.technical.yaml"

    def fingerprint(path: Path | None) -> tuple[object, ...]:
        if path is None:
            return (None,)
        try:
            info = path.stat()
            return (
                str(path.resolve()),
                info.st_dev,
                info.st_ino,
                info.st_size,
                info.st_mtime_ns,
                info.st_ctime_ns,
            )
        except OSError:
            return (str(path), None)

    return (
        config.modkit_executable,
        fingerprint(Path(executable) if executable else None),
        fingerprint(policy),
    )


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
    bam = str(resolve_within(_local_input_path(raw_path), allowed_roots))
    resolve_within(bam, allowed_roots)
    index = str(resolve_bam_index(bam))
    resolve_within(index, allowed_roots)

    mode = AssayMode(str(payload.get("assay", "")))
    genome_build = GenomeBuild(str(payload.get("genome_build", "")))
    target_bed = str(payload.get("target_bed", "")).strip() or None
    target_bed_version = str(payload.get("target_bed_version", "")).strip() or None
    target_bed_wsl = _local_input_path(target_bed) if target_bed else None
    if target_bed_wsl is not None:
        resolve_within(target_bed_wsl, allowed_roots)

    modules = [
        AnalysisModule.QC,
        AnalysisModule.CNV,
        AnalysisModule.SV,
        AnalysisModule.ISCN,
        AnalysisModule.REPORT,
    ]
    if _include_methylation(payload):
        modules.insert(-1, AnalysisModule.METHYLATION)

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
            methylation_policy=(
                _methylation_policy(config, manifest.assay.mode)
                if AnalysisModule.METHYLATION in manifest.analysis.modules
                else None
            ),
            components=config.components,
            threads=config.threads,
            executables={
                "cutesv": config.cutesv_executable,
                "samtools": config.samtools_executable,
                "modkit": config.modkit_executable,
            },
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
    include_methylation = _include_methylation(payload)
    profile_id = str(payload.get("profile", "")).strip()
    if profile_id not in ACTIVE_PROFILE_BUILDS:
        raise ValueError(f"unsupported analysis profile: {profile_id!r}")
    expected_build = ACTIVE_PROFILE_BUILDS[profile_id]
    if config.resource_root is None:
        raise ValueError("profile runs require the service to start with --resource-root")
    registry = ResourceRegistry(config.resource_root, active_build=expected_build)
    if profile_id not in registry.profiles:
        raise ValueError(f"{expected_build.value} analysis profile {profile_id!r} is not active")
    raw_path = str(payload.get("bam", "")).strip()
    if not raw_path:
        raise GuardError("no BAM was selected")
    bam = resolve_within(_local_input_path(raw_path), config.allowed_roots)
    index = resolve_bam_index(bam)
    resolve_within(index, config.allowed_roots)
    legacy_build = str(payload.get("genome_build", "")).strip()
    if legacy_build and legacy_build != expected_build.value:
        raise ValueError(
            f"profile runs are strictly {expected_build.value}; cross-build fallback is prohibited"
        )
    expected_assay = registry.profiles[profile_id].assay_mode.value
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
            include_methylation=include_methylation,
            executables={
                "cutesv": config.cutesv_executable,
                "samtools": config.samtools_executable,
                "modkit": config.modkit_executable,
            },
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
                methylation_policy=(
                    _methylation_policy(config, registry.profiles[profile_id].assay_mode)
                    if include_methylation
                    else None
                ),
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
    ready: list[str] = []
    registries: dict[GenomeBuild, ResourceRegistry] = {}
    for profile_id, build in ACTIVE_PROFILE_BUILDS.items():
        try:
            if build not in registries:
                registries[build] = ResourceRegistry(config.resource_root, active_build=build)
            registry = registries[build]
        except (OSError, ValueError):
            continue
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

    scans = MethylationScans(
        lambda *args, **kwargs: probe_bam_methylation(*args, **kwargs),
        AvailabilityCache(
            lambda: _methylation_availability(config),
            lambda: _methylation_availability_identity(config),
        ),
        samtools=config.samtools_executable,
        validate_path=lambda path: resolve_within(path, config.allowed_roots),
        analysis_running=jobs.running,
    )

    class Handler(BaseHTTPRequestHandler):
        server_version = f"ONTSeq/{__version__}"
        protocol_version = "HTTP/1.1"
        methylation_scans = scans

        def log_message(self, fmt: str, *args: Any) -> None:
            print(f"[{self.log_date_time_string()}] {fmt % args}")

        def _send(self, status: HTTPStatus, body: bytes, content_type: str) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("X-Frame-Options", "DENY")
            self.send_header("Content-Security-Policy", "frame-ancestors 'none'")
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
            if route.path in {"/", "/index.html", "/ONTSeq.html", "/workspace"}:
                # The page carries the session token, so it gets the same rebinding check
                # the API gets. Without it a page on an attacker-controlled name resolved
                # to 127.0.0.1 is same-origin with this service and can simply read the
                # token out of the response. The API's Host check still refuses the token
                # afterwards; not handing it over in the first place is cheaper than
                # relying on that being the only place it matters.
                if not host_is_loopback(self.headers.get("Host"), port=config.port):
                    self._refuse(HTTPStatus.FORBIDDEN, "request did not come from loopback")
                    return
                self._serve_page(workspace=route.path == "/workspace")
                return
            if not route.path.startswith("/api/"):
                self._refuse(HTTPStatus.NOT_FOUND, "no such page")
                return
            if not self._authorised():
                return
            query = parse_qs(route.query)
            if route.path == "/api/config":
                self._json(HTTPStatus.OK, self._config_view())
            elif route.path in {"/api/results", "/api/artifacts"}:
                self._result_artifact(query, download=route.path == "/api/artifacts")
            elif route.path == "/api/methylation":
                self._result_artifact(query, download=False, methylation=True)
            elif route.path.startswith("/api/methylation/scans/"):
                self._methylation_scan_status(route.path[len("/api/methylation/scans/") :])
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
            if path not in {
                "/api/runs",
                "/api/methylation/probe",
                "/api/methylation/scans",
            } and not (
                path.startswith("/api/review/") or path.startswith("/api/methylation/scans/")
            ):
                self.close_connection = True
                self._refuse(HTTPStatus.NOT_FOUND, "no such endpoint")
                return
            if not self._authorised():
                self.close_connection = True
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
            if path == "/api/methylation/probe":
                self._probe_methylation(payload)
                return
            if path == "/api/methylation/scans":
                self._start_methylation_scan(payload)
                return
            if path.startswith("/api/methylation/scans/"):
                suffix = path[len("/api/methylation/scans/") :]
                if not suffix.endswith("/cancel") or payload:
                    self._refuse(
                        HTTPStatus.BAD_REQUEST,
                        "expected /api/methylation/scans/<id>/cancel with an empty object",
                    )
                    return
                self._methylation_scan_status(suffix[: -len("/cancel")], cancel=True)
                return
            parts = path[len("/api/review/") :].split("/")
            if len(parts) != 2 or not all(parts):
                self._refuse(HTTPStatus.BAD_REQUEST, "expected /api/review/<run>/<sample>")
                return
            self._review(parts[0], parts[1], payload)

        def _serve_page(self, *, workspace: bool = False) -> None:
            page = WORKSPACE_PAGE if workspace else PAGE
            if not page.is_file():
                self._refuse(HTTPStatus.NOT_FOUND, "live workspace is not packaged in this build")
                return
            html = page.read_text(encoding="utf-8")
            html = html.replace("__ONTSEQ_TOKEN__", config.token)
            self._send(HTTPStatus.OK, html.encode("utf-8"), "text/html; charset=utf-8")

        def _selected_probe_bam(self, payload: dict[str, Any]) -> Path:
            raw_path = payload.get("bam_path")
            if not isinstance(raw_path, str) or not raw_path.strip():
                raise ValueError("bam_path must identify a selected BAM")
            bam = resolve_within(_local_input_path(raw_path), config.allowed_roots)
            if bam.suffix.lower() != ".bam":
                raise ValueError("tag discovery accepts only a BAM file")
            if not bam.is_file():
                raise ValueError("the selected BAM is missing or unreadable")
            return bam

        def _probe_methylation(self, payload: dict[str, Any]) -> None:
            try:
                force_refresh = payload.get("force_refresh", False)
                if not isinstance(force_refresh, bool):
                    raise ValueError("force_refresh must be an explicit JSON boolean")
                bam = self._selected_probe_bam(payload)
                response = scans.quick(bam, force_refresh=force_refresh)
            except GuardError as error:
                self._refuse(HTTPStatus.FORBIDDEN, str(error))
                return
            except (OSError, ValueError) as error:
                self._refuse(HTTPStatus.BAD_REQUEST, str(error))
                return
            self._json(HTTPStatus.OK, response)

        def _start_methylation_scan(self, payload: dict[str, Any]) -> None:
            try:
                bam = self._selected_probe_bam(payload)
                response = scans.start(bam)
            except GuardError as error:
                self._refuse(HTTPStatus.FORBIDDEN, str(error))
                return
            except ScanBusyError:
                self._json(
                    HTTPStatus.CONFLICT,
                    {"error": "Another BAM discovery scan is in progress.", "reason_code": "busy"},
                )
                return
            except (OSError, ValueError) as error:
                self._refuse(HTTPStatus.BAD_REQUEST, str(error))
                return
            self._json(HTTPStatus.ACCEPTED, response)

        def _methylation_scan_status(self, scan_id: str, *, cancel: bool = False) -> None:
            if re.fullmatch(r"[a-f0-9]{32}", scan_id) is None:
                self._refuse(HTTPStatus.BAD_REQUEST, "invalid methylation scan identifier")
                return
            try:
                response = scans.cancel(scan_id) if cancel else scans.get(scan_id)
            except ScanNotFoundError:
                self._refuse(HTTPStatus.NOT_FOUND, "methylation scan is unavailable or expired")
                return
            self._json(HTTPStatus.OK, response)

        def _result_artifact(
            self, query: dict[str, list[str]], *, download: bool, methylation: bool = False
        ) -> None:
            """Read only named outputs of an identity-checked envelope; never an arbitrary path."""
            run_id = query.get("run_id", [""])[0]
            sample_id = query.get("sample_id", [""])[0]
            kind = query.get("kind", ["json"])[0] if download else "json"
            formats = {
                "json": (RESULT_JSON, "application/json"),
                "html": (REPORT_HTML, "text/html; charset=utf-8"),
                "xlsx": (
                    REPORT_XLSX,
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                ),
            }
            if kind not in formats:
                self._refuse(HTTPStatus.BAD_REQUEST, "unsupported artifact kind")
                return
            try:
                envelope = resolve_envelope(config.output_dir, run_id, sample_id)
                job = jobs.get(run_id)
                if job is not None and job.state == "running":
                    self._refuse(HTTPStatus.CONFLICT, "analysis is still running")
                    return
                result_path = resolve_within(
                    envelope / RESULT_JSON.format(sample=sample_id), [envelope]
                )
                with result_path.open("rb") as stream:
                    result_bytes = stream.read(MAX_RESULT_BYTES + 1)
                if len(result_bytes) > MAX_RESULT_BYTES:
                    raise ValueError("result is too large for the interactive workspace")
                result = PipelineResult.model_validate_json(result_bytes)
                if (result.manifest.run_id, result.manifest.sample_id) != (run_id, sample_id):
                    raise ValueError("result identity does not match its run envelope")
                context = result.reference_context
                if isinstance(context, ResolvedResourceContext) and (
                    context.genome_build != result.manifest.assay.genome_build
                    or context.profile_id != result.manifest.analysis.profile
                ):
                    raise ValueError("result profile/build does not match its reference context")
                if job is not None and (
                    job.sample_id != sample_id
                    or job.profile != result.manifest.analysis.profile
                    or job.detected_genome_build != result.manifest.assay.genome_build.value
                ):
                    raise ValueError("result profile/build does not match the executed job")
                relative, content_type = formats[kind]
                path = resolve_within(envelope / relative.format(sample=sample_id), [envelope])
                if path.stat().st_size > MAX_RESULT_BYTES:
                    raise ValueError("artifact is too large for the interactive workspace")
                body = path.read_bytes() if download else result.model_dump_json().encode("utf-8")
                if methylation:
                    if (envelope / LOCK_FILENAME).exists():
                        self._refuse(HTTPStatus.CONFLICT, "the run envelope is still locked")
                        return
                    relative = METHYLATION_REPORT.format(sample=sample_id)
                    path = resolve_within(envelope / relative, [envelope])
                    run_path = resolve_within(envelope / RUN_REPORT, [envelope])
                    with run_path.open("rb") as stream:
                        run_bytes = stream.read(MAX_RESULT_BYTES + 1)
                    with path.open("rb") as stream:
                        report_bytes = stream.read(MAX_RESULT_BYTES + 1)
                    if max(len(run_bytes), len(report_bytes)) > MAX_RESULT_BYTES:
                        raise ValueError("methylation artifact is too large for the workspace")
                    run = RunReport.model_validate_json(run_bytes)
                    if (run.run_id, run.sample_id, run.genome_build) != (
                        run_id,
                        sample_id,
                        result.manifest.assay.genome_build,
                    ) or run.manifest.analysis.profile != result.manifest.analysis.profile:
                        raise ValueError("methylation run identity does not match the result")
                    assembly = [stage for stage in run.stages if stage.stage.value == "assemble"]
                    if len(assembly) != 1 or assembly[0].status.value != "COMPLETED":
                        raise ValueError("the current run has no completed result assembly")
                    assembled = [
                        item
                        for item in assembly[0].outputs
                        if item.relative_path == RESULT_JSON.format(sample=sample_id)
                    ]
                    if len(assembled) != 1 or (
                        assembled[0].size_bytes != len(result_bytes)
                        or assembled[0].sha256 != hashlib.sha256(result_bytes).hexdigest()
                    ):
                        raise ValueError(
                            "result artifact does not match its current assembly checksum"
                        )
                    stages = [stage for stage in run.stages if stage.stage.value == "methylation"]
                    if len(stages) != 1 or stages[0].status.value not in {"COMPLETED", "NO_CALL"}:
                        raise ValueError("methylation did not produce a concluded artifact")
                    records = [item for item in stages[0].outputs if item.relative_path == relative]
                    if len(records) != 1 or (
                        records[0].size_bytes != len(report_bytes)
                        or records[0].sha256 != hashlib.sha256(report_bytes).hexdigest()
                    ):
                        raise ValueError("methylation artifact does not match its run checksum")
                    methylation_report = MethylationReport.model_validate_json(report_bytes)
                    if (methylation_report.sample_id, methylation_report.genome_build) != (
                        sample_id,
                        result.manifest.assay.genome_build,
                    ) or methylation_report.status != stages[0].status:
                        raise ValueError(
                            "methylation report identity/status does not match the run"
                        )
                    if AnalysisModule.METHYLATION not in result.manifest.analysis.modules:
                        raise ValueError("methylation was not requested in the run manifest")
                    outcomes = [
                        item for item in result.modules if item.module == AnalysisModule.METHYLATION
                    ]
                    if (
                        len(outcomes) != 1
                        or outcomes[0].status != methylation_report.status
                        or methylation_report.tool not in outcomes[0].tools
                        or methylation_report.tool not in stages[0].tools
                    ):
                        raise ValueError("methylation outcome/tool does not match the run result")
                    if result.provenance.reference_checksums.get("bedmethyl") != (
                        methylation_report.bedmethyl_fingerprint.sha256
                    ):
                        raise ValueError("methylation source provenance does not match the result")
                    body = json.dumps(
                        {**methylation_report.model_dump(mode="json"), "run_id": run_id}
                    ).encode("utf-8")
            except GuardError as error:
                self._refuse(HTTPStatus.FORBIDDEN, str(error))
                return
            except FileNotFoundError:
                self._refuse(HTTPStatus.NOT_FOUND, "requested result artifact is not available")
                return
            except (OSError, ValueError) as error:
                self._refuse(HTTPStatus.BAD_REQUEST, str(error))
                return
            if not download:
                self._send(HTTPStatus.OK, body, content_type)
                return
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Content-Disposition", f'attachment; filename="{path.name}"')
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Security-Policy", "sandbox")
            self.end_headers()
            self.wfile.write(body)

        def _config_view(self) -> dict[str, Any]:
            not_wired = ["basecalling — not needed for an already aligned BAM"]
            return {
                **scans.availability(),
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
                "instance_id": config.instance_id,
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
            if scans.discovery_active():
                self._json(
                    HTTPStatus.CONFLICT,
                    {
                        "error": "Wait for BAM discovery to finish or acknowledge cancellation.",
                        "reason_code": "busy",
                    },
                )
                return
            profile_id = str(payload.get("profile", "")).strip()
            try:
                _include_methylation(payload)
                resolved_config: RunConfiguration | None = None
                policy: MethylationPolicy | None = None
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
                    if AnalysisModule.METHYLATION in manifest.analysis.modules:
                        policy = _methylation_policy(config, manifest.assay.mode)
                        if policy.cpg_only and (
                            config.reference_fasta is None or not config.reference_fasta.is_file()
                        ):
                            raise ValueError("Methylation requires the locked reference FASTA.")
                        if (
                            policy.region_source.value == "target_bed"
                            and not manifest.assay.target_bed
                        ):
                            raise ValueError("Target-region methylation requires an analysis BED.")
                if AnalysisModule.METHYLATION in manifest.analysis.modules:
                    selected_policy = (
                        resolved_config.methylation_policy
                        if resolved_config is not None
                        else policy
                    )
                    if selected_policy is None:
                        raise ValueError("Methylation requires a versioned policy.")
                    _verify_methylation_runtime(selected_policy, config.modkit_executable)
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
                scans.claim_analysis(lambda: jobs.claim(job))
            except ScanBusyError:
                self._json(
                    HTTPStatus.CONFLICT,
                    {
                        "error": "Wait for BAM discovery to finish or acknowledge cancellation.",
                        "reason_code": "busy",
                    },
                )
                return
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


def close_handler_resources(handler: type[BaseHTTPRequestHandler]) -> None:
    """Stop discovery workers without touching independent analysis jobs."""
    scans = getattr(handler, "methylation_scans", None)
    if isinstance(scans, MethylationScans):
        scans.close()


def serve(config: ServiceConfig, *, open_browser: bool = True) -> None:
    """Run the service until interrupted."""
    jobs = Jobs()
    handler = make_handler(config, jobs)
    server = ThreadingHTTPServer((config.host, config.port), handler)
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
        close_handler_resources(handler)
        server.server_close()
