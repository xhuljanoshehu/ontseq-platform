from __future__ import annotations

import os
import queue
import re
import shlex
import subprocess
import threading
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import TextIO

import yaml

from ..io import load_model
from ..models import (
    AnalysisModule,
    AssayMode,
    GenomeBuild,
    ModuleRunStatus,
    PipelineResult,
    ReferenceLock,
    Verdict,
)
from .config import DesktopConfig, DesktopReferenceProfile


class DesktopBackendError(RuntimeError):
    """Raised when a desktop analysis cannot safely continue."""


class DesktopStage(StrEnum):
    INPUT = "input"
    QC = "qc"
    CNV = "cnv"
    SV = "sv"
    FUSION = "fusion"
    ISCN = "iscn"
    REPORT = "report"


class DesktopStageStatus(StrEnum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"
    NOT_RUN = "NOT_RUN"
    NO_CALL = "NO_CALL"


@dataclass(frozen=True)
class ProgressEvent:
    stage: DesktopStage
    status: DesktopStageStatus
    message: str
    percent: int


@dataclass(frozen=True)
class DesktopAnalysisRequest:
    bam_path: Path
    sample_id: str
    genome_build: GenomeBuild
    assay_mode: AssayMode
    output_dir: Path
    run_id: str | None = None


@dataclass(frozen=True)
class DesktopRunResult:
    sample_id: str
    output_dir: Path
    result_json: Path
    report_html: Path
    workbook_xlsx: Path
    run_log: Path
    pipeline_result: PipelineResult


@dataclass(frozen=True)
class DiagnosticCheck:
    name: str
    ok: bool
    detail: str


ProgressCallback = Callable[[ProgressEvent], None]
LogCallback = Callable[[str], None]

_SAMPLE_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def sanitize_sample_id(value: str) -> str:
    candidate = _SAMPLE_SAFE.sub("_", value.strip()).strip("._-")
    if len(candidate) < 3:
        candidate = f"S_{candidate or 'sample'}"
    return candidate[:64]


def locate_bam_index(bam_path: Path) -> Path | None:
    candidates = (Path(f"{bam_path}.bai"), bam_path.with_suffix(".bai"))
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def default_run_id() -> str:
    return datetime.now(UTC).strftime("RUN_%Y%m%d_%H%M%S")


class DesktopBackend:
    """Run the existing ONTSeq CLI through a local or WSL2 boundary."""

    def __init__(self, config: DesktopConfig) -> None:
        self.config = config

    def diagnose(self) -> list[DiagnosticCheck]:
        checks: list[DiagnosticCheck] = []
        checks.append(self._check_output_folder())
        if self.config.backend_mode == "wsl":
            checks.extend((self._check_wsl(), self._check_wsl_tools()))
        else:
            root = Path(self.config.wsl_project_root).expanduser()
            checks.append(DiagnosticCheck("Local project root", root.is_dir(), str(root)))
        checks.extend(self._check_reference_resources())
        return checks

    def run(
        self,
        request: DesktopAnalysisRequest,
        *,
        on_progress: ProgressCallback | None = None,
        on_log: LogCallback | None = None,
        cancel_event: threading.Event | None = None,
    ) -> DesktopRunResult:
        cancel = cancel_event or threading.Event()
        output_dir = request.output_dir.expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        log_path = output_dir / "ontseq-desktop.log"

        with log_path.open("a", encoding="utf-8") as log_handle:
            self._log(log_handle, on_log, "ONTSeq Desktop run started")
            try:
                return self._run_pipeline(
                    request,
                    output_dir=output_dir,
                    log_handle=log_handle,
                    on_progress=on_progress,
                    on_log=on_log,
                    cancel=cancel,
                )
            except Exception as exc:
                self._emit(
                    on_progress,
                    DesktopStage.REPORT,
                    DesktopStageStatus.FAIL,
                    str(exc),
                    100,
                )
                self._log(log_handle, on_log, f"FAILED: {exc}")
                raise

    def _run_pipeline(
        self,
        request: DesktopAnalysisRequest,
        *,
        output_dir: Path,
        log_handle: TextIO,
        on_progress: ProgressCallback | None,
        on_log: LogCallback | None,
        cancel: threading.Event,
    ) -> DesktopRunResult:
        bam = request.bam_path.expanduser().resolve()
        if not bam.is_file() or bam.suffix.lower() != ".bam":
            raise DesktopBackendError("Please select an existing aligned .bam file")
        bai = locate_bam_index(bam)
        if bai is None:
            raise DesktopBackendError(
                "No BAM index found. Expected <sample>.bam.bai or <sample>.bai next to the BAM."
            )

        sample_id = sanitize_sample_id(request.sample_id)
        profile = self._validate_reference_profile(request.genome_build, request.assay_mode)
        paths = self._artifact_paths(output_dir, sample_id)
        manifest = self._build_manifest(
            request,
            sample_id=sample_id,
            bam=bam,
            bai=bai,
            profile=profile,
        )
        paths["manifest"].write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")

        backend = {name: self._backend_path(path) for name, path in paths.items()}
        reference_lock = self._backend_path(Path(profile.reference_lock_path))

        self._emit(
            on_progress,
            DesktopStage.INPUT,
            DesktopStageStatus.RUNNING,
            "Input and reference provenance are being checked",
            5,
        )
        self._run_cli(
            [
                "inspect-bam",
                backend["manifest"],
                "--reference-lock",
                reference_lock,
                "--output",
                backend["intake"],
            ],
            log_handle=log_handle,
            on_log=on_log,
            cancel=cancel,
        )
        self._emit(
            on_progress,
            DesktopStage.INPUT,
            DesktopStageStatus.PASS,
            "BAM, index and reference gate passed",
            20,
        )

        self._emit(
            on_progress,
            DesktopStage.QC,
            DesktopStageStatus.RUNNING,
            "Cramino quality control is running",
            25,
        )
        self._run_cli(
            [
                "qc-cramino",
                backend["manifest"],
                "--policy",
                self.config.qc_policy_path,
                "--output",
                backend["qc"],
            ],
            log_handle=log_handle,
            on_log=on_log,
            cancel=cancel,
        )
        self._emit(
            on_progress,
            DesktopStage.QC,
            DesktopStageStatus.PASS,
            "Quality control completed",
            40,
        )

        self._emit(
            on_progress,
            DesktopStage.SV,
            DesktopStageStatus.RUNNING,
            "Sniffles2 candidate SV calling is running",
            45,
        )
        self._run_cli(
            [
                "call-sniffles",
                backend["manifest"],
                "--intake",
                backend["intake"],
                "--policy",
                self.config.sniffles_policy_path,
                "--vcf",
                backend["sniffles_vcf"],
                "--output",
                backend["sniffles"],
            ],
            log_handle=log_handle,
            on_log=on_log,
            cancel=cancel,
        )
        self._emit(
            on_progress,
            DesktopStage.SV,
            DesktopStageStatus.PASS,
            "Candidate SV evidence normalized",
            65,
        )

        self._run_cli(
            [
                "assemble-aligned-mvp",
                backend["manifest"],
                "--intake",
                backend["intake"],
                "--qc",
                backend["qc"],
                "--sniffles",
                backend["sniffles"],
                "--git-commit",
                "DESKTOP_LOCAL",
                "--output",
                backend["result_json"],
            ],
            log_handle=log_handle,
            on_log=on_log,
            cancel=cancel,
        )

        self._emit(
            on_progress,
            DesktopStage.REPORT,
            DesktopStageStatus.RUNNING,
            "HTML, Excel and JSON outputs are being generated",
            80,
        )
        self._run_cli(
            [
                "render",
                backend["result_json"],
                "--output-dir",
                self._backend_path(output_dir),
            ],
            log_handle=log_handle,
            on_log=on_log,
            cancel=cancel,
        )

        pipeline_result = load_model(paths["result_json"], PipelineResult)
        self._emit_result_module_states(pipeline_result, on_progress)
        self._require_outputs(paths)
        report_status = (
            DesktopStageStatus.WARN
            if pipeline_result.qc.verdict == Verdict.WARN or pipeline_result.warnings
            else DesktopStageStatus.PASS
        )
        self._emit(
            on_progress,
            DesktopStage.REPORT,
            report_status,
            "Analysis finished. Expert review remains required.",
            100,
        )
        self._log(log_handle, on_log, "ONTSeq Desktop run finished")
        return DesktopRunResult(
            sample_id=sample_id,
            output_dir=output_dir,
            result_json=paths["result_json"],
            report_html=paths["report_html"],
            workbook_xlsx=paths["workbook_xlsx"],
            run_log=output_dir / "ontseq-desktop.log",
            pipeline_result=pipeline_result,
        )

    def _check_output_folder(self) -> DiagnosticCheck:
        root = Path(self.config.output_root).expanduser()
        try:
            root.mkdir(parents=True, exist_ok=True)
            probe = root / ".ontseq-write-test"
            probe.write_text("ok", encoding="utf-8")
            probe.unlink()
        except OSError as exc:
            return DiagnosticCheck("Output folder", False, str(exc))
        return DiagnosticCheck("Output folder", True, str(root))

    def _check_wsl(self) -> DiagnosticCheck:
        try:
            completed = subprocess.run(
                ["wsl.exe", "--status"],
                check=False,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=15,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            return DiagnosticCheck("WSL", False, str(exc))
        detail = (completed.stdout or completed.stderr or "WSL responded").strip()
        return DiagnosticCheck("WSL", completed.returncode == 0, detail)

    def _check_wsl_tools(self) -> DiagnosticCheck:
        shell = (
            "for c in ontseq samtools cramino sniffles; do "
            'command -v "$c" >/dev/null || { echo missing:$c; exit 7; }; '
            "done; echo ready"
        )
        try:
            completed = subprocess.run(
                self._wsl_shell_command(shell),
                check=False,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=20,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            return DiagnosticCheck("ONTSeq backend tools", False, str(exc))
        detail = (completed.stdout or completed.stderr).strip()
        return DiagnosticCheck("ONTSeq backend tools", completed.returncode == 0, detail)

    def _check_reference_resources(self) -> list[DiagnosticCheck]:
        checks: list[DiagnosticCheck] = []
        for profile in self.config.reference_profiles:
            lock = Path(profile.reference_lock_path).expanduser()
            checks.append(
                DiagnosticCheck(
                    f"{profile.genome_build.value} reference lock",
                    lock.is_file(),
                    str(lock),
                )
            )
            if profile.adaptive_sampling_target_bed_path:
                bed = Path(profile.adaptive_sampling_target_bed_path).expanduser()
                checks.append(
                    DiagnosticCheck(
                        f"{profile.genome_build.value} adaptive target BED",
                        bed.is_file(),
                        str(bed),
                    )
                )
        return checks

    @staticmethod
    def _artifact_paths(output_dir: Path, sample_id: str) -> dict[str, Path]:
        return {
            "manifest": output_dir / f"{sample_id}.manifest.yaml",
            "intake": output_dir / f"{sample_id}.intake.json",
            "qc": output_dir / f"{sample_id}.qc.json",
            "sniffles_vcf": output_dir / f"{sample_id}.sniffles.vcf",
            "sniffles": output_dir / f"{sample_id}.sniffles.json",
            "result_json": output_dir / f"{sample_id}.result.json",
            "report_html": output_dir / f"{sample_id}.report.html",
            "workbook_xlsx": output_dir / f"{sample_id}.results.xlsx",
        }

    @staticmethod
    def _require_outputs(paths: dict[str, Path]) -> None:
        expected = (paths["result_json"], paths["report_html"], paths["workbook_xlsx"])
        missing = [path.name for path in expected if not path.is_file()]
        if missing:
            raise DesktopBackendError("Expected output was not created: " + ", ".join(missing))

    def _validate_reference_profile(
        self, genome_build: GenomeBuild, assay_mode: AssayMode
    ) -> DesktopReferenceProfile:
        profile = self.config.reference_for(genome_build)
        lock_path = Path(profile.reference_lock_path).expanduser().resolve()
        if not lock_path.is_file():
            raise DesktopBackendError(f"Reference lock not found: {lock_path}")
        lock = load_model(lock_path, ReferenceLock)
        if lock.genome_build != genome_build:
            raise DesktopBackendError(
                f"Reference lock build is {lock.genome_build.value}, expected {genome_build.value}"
            )
        if lock.reference_id != profile.reference_id:
            raise DesktopBackendError(
                "Configured reference_id does not match the selected reference lock"
            )
        if assay_mode == AssayMode.ADAPTIVE_SAMPLING:
            if not profile.adaptive_sampling_target_bed_path:
                raise DesktopBackendError(
                    "Adaptive Sampling requires a configured target BED and version"
                )
            bed = Path(profile.adaptive_sampling_target_bed_path).expanduser().resolve()
            if not bed.is_file():
                raise DesktopBackendError(f"Adaptive Sampling target BED not found: {bed}")
        return profile

    def _build_manifest(
        self,
        request: DesktopAnalysisRequest,
        *,
        sample_id: str,
        bam: Path,
        bai: Path,
        profile: DesktopReferenceProfile,
    ) -> dict[str, object]:
        assay: dict[str, object] = {
            "mode": request.assay_mode.value,
            "genome_build": request.genome_build.value,
            "reference_id": profile.reference_id,
        }
        if request.assay_mode == AssayMode.ADAPTIVE_SAMPLING:
            bed_path = profile.adaptive_sampling_target_bed_path
            bed_version = profile.adaptive_sampling_target_bed_version
            if bed_path is None or bed_version is None:
                raise DesktopBackendError("Adaptive Sampling target configuration is incomplete")
            assay["target_bed"] = self._backend_path(Path(bed_path))
            assay["target_bed_version"] = bed_version

        return {
            "schema_version": "0.1.0",
            "sample_id": sample_id,
            "run_id": request.run_id or default_run_id(),
            "input": {
                "kind": "aligned_bam",
                "path": self._backend_path(bam),
                "index_path": self._backend_path(bai),
            },
            "assay": assay,
            "analysis": {
                "profile": f"desktop-{request.assay_mode.value}",
                "modules": [
                    AnalysisModule.QC.value,
                    AnalysisModule.CNV.value,
                    AnalysisModule.SV.value,
                    AnalysisModule.FUSION.value,
                    AnalysisModule.ISCN.value,
                    AnalysisModule.REPORT.value,
                ],
                "parameters": {},
            },
            "privacy": {
                "pseudonymized": True,
                "contains_direct_identifiers": False,
                "cloud_upload_approved": False,
            },
        }

    def _emit_result_module_states(
        self, result: PipelineResult, on_progress: ProgressCallback | None
    ) -> None:
        stage_map = {
            AnalysisModule.CNV: DesktopStage.CNV,
            AnalysisModule.SV: DesktopStage.SV,
            AnalysisModule.FUSION: DesktopStage.FUSION,
            AnalysisModule.ISCN: DesktopStage.ISCN,
        }
        status_map = {
            ModuleRunStatus.COMPLETED: DesktopStageStatus.PASS,
            ModuleRunStatus.NOT_RUN: DesktopStageStatus.NOT_RUN,
            ModuleRunStatus.FAILED: DesktopStageStatus.FAIL,
            ModuleRunStatus.NO_CALL: DesktopStageStatus.NO_CALL,
        }
        for outcome in result.modules:
            stage = stage_map.get(outcome.module)
            if stage is not None:
                self._emit(on_progress, stage, status_map[outcome.status], outcome.reason, 90)

    def _backend_path(self, path: Path) -> str:
        local_path = path.expanduser().resolve()
        if self.config.backend_mode == "local":
            return str(local_path)
        try:
            completed = subprocess.run(
                self._wsl_exec(["wslpath", "-a", "-u", str(local_path)]),
                check=False,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=20,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise DesktopBackendError(f"Unable to translate path for WSL: {local_path}") from exc
        translated = completed.stdout.strip()
        if completed.returncode != 0 or not translated:
            detail = (completed.stderr or completed.stdout).strip()
            raise DesktopBackendError(f"WSL path translation failed for {local_path}: {detail}")
        return translated

    def _run_cli(
        self,
        args: list[str],
        *,
        log_handle: TextIO,
        on_log: LogCallback | None,
        cancel: threading.Event,
    ) -> None:
        if cancel.is_set():
            raise DesktopBackendError("Analysis cancelled")
        if self.config.backend_mode == "wsl":
            command = self._wsl_shell_command(shlex.join(["ontseq", *args]))
            cwd = None
        else:
            command = ["ontseq", *args]
            cwd = Path(self.config.wsl_project_root).expanduser()
        self._run_process(
            command,
            cwd=cwd,
            log_handle=log_handle,
            on_log=on_log,
            cancel=cancel,
        )

    def _run_process(
        self,
        command: list[str],
        *,
        cwd: Path | None,
        log_handle: TextIO,
        on_log: LogCallback | None,
        cancel: threading.Event,
    ) -> None:
        self._log(log_handle, on_log, f"RUN: {shlex.join(command)}")
        try:
            process: subprocess.Popen[str] = subprocess.Popen(
                command,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
        except OSError as exc:
            raise DesktopBackendError(f"Could not start backend command: {exc}") from exc

        output_queue: queue.Queue[str] = queue.Queue()

        def read_output() -> None:
            if process.stdout is None:
                return
            for line in process.stdout:
                output_queue.put(line.rstrip())

        reader = threading.Thread(target=read_output, daemon=True)
        reader.start()
        recent: list[str] = []
        while process.poll() is None:
            self._drain_output(output_queue, recent, log_handle, on_log)
            if cancel.wait(0.15):
                self._terminate_process(process)
                raise DesktopBackendError("Analysis cancelled")
        reader.join(timeout=2)
        self._drain_output(output_queue, recent, log_handle, on_log)
        if process.returncode != 0:
            tail = " | ".join(recent[-8:])
            detail = f": {tail}" if tail else ""
            raise DesktopBackendError(
                f"Backend command failed with exit code {process.returncode}{detail}"
            )

    def _drain_output(
        self,
        output_queue: queue.Queue[str],
        recent: list[str],
        log_handle: TextIO,
        on_log: LogCallback | None,
    ) -> None:
        while True:
            try:
                line = output_queue.get_nowait()
            except queue.Empty:
                return
            recent.append(line)
            if len(recent) > 30:
                del recent[:-30]
            self._log(log_handle, on_log, line)

    @staticmethod
    def _terminate_process(process: subprocess.Popen[str]) -> None:
        if process.poll() is not None:
            return
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                check=False,
                capture_output=True,
                timeout=10,
            )
        else:
            process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()

    def _wsl_base(self) -> list[str]:
        command = ["wsl.exe"]
        if self.config.wsl_distribution:
            command.extend(["--distribution", self.config.wsl_distribution])
        return command

    def _wsl_exec(self, argv: list[str]) -> list[str]:
        return [*self._wsl_base(), "--exec", *argv]

    def _wsl_shell_command(self, command: str) -> list[str]:
        root = self._project_root_expression(self.config.wsl_project_root)
        return self._wsl_exec(["bash", "-lc", f"cd {root} && {command}"])

    @staticmethod
    def _project_root_expression(value: str) -> str:
        if value == "~":
            return '"$HOME"'
        if value.startswith("~/"):
            return '"$HOME"/' + shlex.quote(value[2:])
        return shlex.quote(value)

    @staticmethod
    def _emit(
        callback: ProgressCallback | None,
        stage: DesktopStage,
        status: DesktopStageStatus,
        message: str,
        percent: int,
    ) -> None:
        if callback is not None:
            callback(ProgressEvent(stage, status, message, max(0, min(percent, 100))))

    @staticmethod
    def _log(log_handle: TextIO, callback: LogCallback | None, line: str) -> None:
        timestamp = datetime.now(UTC).isoformat(timespec="seconds")
        rendered = f"[{timestamp}] {line}"
        log_handle.write(rendered + "\n")
        log_handle.flush()
        if callback is not None:
            callback(rendered)
