from __future__ import annotations

import os
import queue
import re
import shlex
import subprocess
import threading
import time
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
    pass


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
    candidates = (
        Path(f"{bam_path}.bai"),
        bam_path.with_suffix(".bai"),
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def default_run_id() -> str:
    return datetime.now(UTC).strftime("RUN_%Y%m%d_%H%M%S")


class DesktopBackend:
    def __init__(self, config: DesktopConfig) -> None:
        self.config = config

    def diagnose(self) -> list[DiagnosticCheck]:
        checks: list[DiagnosticCheck] = []
        output_root = Path(self.config.output_root).expanduser()
        try:
            output_root.mkdir(parents=True, exist_ok=True)
            probe = output_root / ".ontseq-write-test"
            probe.write_text("ok", encoding="utf-8")
            probe.unlink()
            checks.append(DiagnosticCheck("Output folder", True, str(output_root)))
        except OSError as exc:
            checks.append(DiagnosticCheck("Output folder", False, str(exc)))

        if self.config.backend_mode == "wsl":
            try:
                completed = subprocess.run(
                    self._wsl_base() + ["--status"],
                    check=False,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=15,
                )
                checks.append(
                    DiagnosticCheck(
                        "WSL",
                        completed.returncode == 0,
                        (completed.stdout or completed.stderr or "WSL responded").strip(),
                    )
                )
            except (OSError, subprocess.SubprocessError) as exc:
                checks.append(DiagnosticCheck("WSL", False, str(exc)))

            tool_check = (
                "for c in ontseq samtools cramino sniffles; do "
                "command -v \"$c\" >/dev/null || { echo missing:$c; exit 7; }; done; "
                "echo ready"
            )
            try:
                completed = subprocess.run(
                    self._wsl_shell_command(tool_check),
                    check=False,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=20,
                )
                checks.append(
                    DiagnosticCheck(
                        "ONTSeq backend tools",
                        completed.returncode == 0,
                        (completed.stdout or completed.stderr).strip(),
                    )
                )
            except (OSError, subprocess.SubprocessError) as exc:
                checks.append(DiagnosticCheck("ONTSeq backend tools", False, str(exc)))
        else:
            project_root = Path(self.config.wsl_project_root).expanduser()
            checks.append(
                DiagnosticCheck(
                    "Local project root",
                    project_root.is_dir(),
                    str(project_root),
                )
            )

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
                return self._run(
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

    def _run(
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

        self._emit(
            on_progress,
            DesktopStage.INPUT,
            DesktopStageStatus.RUNNING,
            "Input and reference provenance are being checked",
            5,
        )

        manifest_path = output_dir / f"{sample_id}.manifest.yaml"
        intake_path = output_dir / f"{sample_id}.intake.json"
        qc_path = output_dir / f"{sample_id}.qc.json"
        sniffles_vcf = output_dir / f"{sample_id}.sniffles.vcf"
        sniffles_path = output_dir / f"{sample_id}.sniffles.json"
        result_json = output_dir / f"{sample_id}.result.json"

        manifest = self._build_manifest(
            request,
            sample_id=sample_id,
            bam=bam,
            bai=bai,
            profile=profile,
        )
        manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")

        backend_manifest = self._backend_path(manifest_path)
        backend_reference_lock = self._backend_path(Path(profile.reference_lock_path))
        backend_intake = self._backend_path(intake_path)
        backend_qc = self._backend_path(qc_path)
        backend_sniffles_vcf = self._backend_path(sniffles_vcf)
        backend_sniffles = self._backend_path(sniffles_path)
        backend_result = self._backend_path(result_json)
        backend_output_dir = self._backend_path(output_dir)

        self._run_cli(
            [
                "inspect-bam",
                backend_manifest,
                "--reference-lock",
                backend_reference_lock,
                "--output",
                backend_intake,
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
                backend_manifest,
                "--policy",
                self.config.qc_policy_path,
                "--output",
                backend_qc,
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
                backend_manifest,
                "--intake",
                backend_intake,
                "--policy",
                self.config.sniffles_policy_path,
                "--vcf",
                backend_sniffles_vcf,
                "--output",
                backend_sniffles,
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
                backend_manifest,
                "--intake",
                backend_intake,
                "--qc",
                backend_qc,
                "--sniffles",
                backend_sniffles,
                "--git-commit",
                "DESKTOP_LOCAL",
                "--output",
                backend_result,
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
            ["render", backend_result, "--output-dir", backend_output_dir],
            log_handle=log_handle,
            on_log=on_log,
            cancel=cancel,
        )

        pipeline_result = load_model(result_json, PipelineResult)
        self._emit_result_module_states(pipeline_result, on_progress)
        report_html = output_dir / f"{sample_id}.report.html"
        workbook_xlsx = output_dir / f"{sample_id}.results.xlsx"
        missing = [
            path
            for path in (result_json, report_html, workbook_xlsx)
            if not path.is_file()
        ]
        if missing:
            raise DesktopBackendError(
                "Expected output was not created: " + ", ".join(path.name for path in missing)
            )

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
            result_json=result_json,
            report_html=report_html,
            workbook_xlsx=workbook_xlsx,
            run_log=output_dir / "ontseq-desktop.log",
            pipeline_result=pipeline_result,
        )

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
            if (
                profile.adaptive_sampling_target_bed_path is None
                or profile.adaptive_sampling_target_bed_version is None
            ):
                raise DesktopBackendError("Adaptive Sampling target configuration is incomplete")
            assay["target_bed"] = self._backend_path(
                Path(profile.adaptive_sampling_target_bed_path)
            )
            assay["target_bed_version"] = profile.adaptive_sampling_target_bed_version

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
            if stage is None:
                continue
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
        self._log(log_handle, on_log, f"RUN: {self._redacted_command(command)}")
        try:
            process = subprocess.Popen(
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
            assert process.stdout is not None
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

    def _terminate_process(self, process: subprocess.Popen[str]) -> None:
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
        return self._wsl_base() + ["--exec", *argv]

    def _wsl_shell_command(self, command: str) -> list[str]:
        project_root = self._project_root_expression(self.config.wsl_project_root)
        shell_command = f"cd {project_root} && {command}"
        return self._wsl_exec(["bash", "-lc", shell_command])

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

    @staticmethod
    def _redacted_command(command: list[str]) -> str:
        # Commands are written only to the local run log. Keep the representation shell-safe and
        # avoid adding any data content beyond already-local filesystem paths and tool arguments.
        return shlex.join(command)
