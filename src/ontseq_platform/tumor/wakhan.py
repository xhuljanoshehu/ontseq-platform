from __future__ import annotations

import hashlib
import os
import re
import shutil
import tempfile
from pathlib import Path, PurePosixPath
from typing import Literal

from pydantic import Field, field_validator, model_validator

from ..execution import CommandRunner, SubprocessRunner
from ..models import FileFingerprint, GenomeBuild, ModuleRunStatus, StrictModel, ToolRecord
from ..multicaller_contracts import CallerInputRole
from ..reference import sha256_file
from ..tumor_inputs import TumorAuxiliaryInputArtifact, TumorInputBundle

_SHA256 = r"^[0-9a-f]{64}$"
_SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")


class WakhanPhasedCnaPolicy(StrictModel):
    """Research-only contract for Wakhan phased CNA execution."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    mode: Literal["tumor_normal", "tumor_only"]
    expected_version: str = Field(pattern=r"^\d+\.\d+(?:\.\d+)?$")
    real_tool_qualified: bool = False
    analytical_validation: Literal["not_validated"] = "not_validated"
    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    reference_sha256: str = Field(pattern=_SHA256)
    timeout_seconds: int = Field(ge=60)
    note: str = Field(min_length=12)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def pinned_runtime_version(self) -> WakhanPhasedCnaPolicy:
        if self.expected_version != "0.4.4":
            raise ValueError("Wakhan expected_version must match the pinned 0.4.4 contract")
        return self


class WakhanNativeArtifact(StrictModel):
    relative_path: str = Field(min_length=1)
    role: str = Field(min_length=1)
    media_type: str = Field(min_length=1)
    fingerprint: FileFingerprint
    artifact_kind: Literal["file", "symlink"] = "file"
    symlink_target: str | None = None
    sensitive_output: Literal[True] = True
    exportable: Literal[False] = False
    research_only: Literal[True] = True

    @field_validator("relative_path", "symlink_target")
    @classmethod
    def safe_relative_path(cls, value: str | None) -> str | None:
        if value is None:
            return None
        if "\\" in value:
            raise ValueError("Wakhan artifact paths must use POSIX separators")
        path = PurePosixPath(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Wakhan artifact path must remain inside output directory")
        return value

    @model_validator(mode="after")
    def coherent_artifact_kind(self) -> WakhanNativeArtifact:
        if self.artifact_kind == "symlink" and self.symlink_target is None:
            raise ValueError("Wakhan symlink artifact requires a relative target")
        if self.artifact_kind == "file" and self.symlink_target is not None:
            raise ValueError("Wakhan file artifact cannot carry a symlink target")
        return self


class WakhanReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str = Field(min_length=1)
    mode: Literal["tumor_normal", "tumor_only"]
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: WakhanPhasedCnaPolicy
    tool: ToolRecord
    real_tool_qualified: bool
    tumor_bam_fingerprint: FileFingerprint
    tumor_bam_index_fingerprint: FileFingerprint
    phased_vcf_fingerprint: FileFingerprint
    reference_fingerprint: FileFingerprint
    runtime_script_fingerprint: FileFingerprint
    breakpoints_fingerprint: FileFingerprint | None = None
    breakpoint_parent_lane_id: str | None = None
    breakpoint_parent_lane_sha256: str | None = Field(default=None, pattern=_SHA256)
    native_artifacts: list[WakhanNativeArtifact] = Field(min_length=1)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_report(self) -> WakhanReport:
        if self.mode != self.policy.mode:
            raise ValueError("Wakhan report mode must match policy mode")
        if self.real_tool_qualified != self.policy.real_tool_qualified:
            raise ValueError("Wakhan real-tool qualification must match policy")
        parent_pair = (
            self.breakpoint_parent_lane_id is not None,
            self.breakpoint_parent_lane_sha256 is not None,
        )
        if parent_pair[0] != parent_pair[1]:
            raise ValueError("Wakhan breakpoint parent ID and SHA-256 must be declared together")
        if self.breakpoints_fingerprint is None and parent_pair[0]:
            raise ValueError("Wakhan breakpoint parent lock requires breakpoint evidence")
        if self.breakpoints_fingerprint is not None and not parent_pair[0]:
            raise ValueError("Wakhan breakpoint evidence requires an exact parent lane lock")
        if self.status != ModuleRunStatus.COMPLETED:
            raise ValueError("Wakhan runtime report currently represents completed executions only")
        return self


def build_wakhan_argv(
    *,
    python_executable: str,
    wakhan_script: Path,
    tumor_bam: Path,
    phased_vcf: Path,
    reference_fasta: Path,
    sample_id: str,
    output_dir: Path,
    policy: WakhanPhasedCnaPolicy,
    threads: int = 8,
    breakpoints_vcf: Path | None = None,
) -> tuple[str, ...]:
    if threads < 1:
        raise ValueError("Wakhan threads must be at least 1")
    argv = [
        python_executable,
        str(wakhan_script),
        "all",
        "--threads",
        str(threads),
        "--reference",
        str(reference_fasta),
        "--target-bam",
        str(tumor_bam),
    ]
    if policy.mode == "tumor_normal":
        argv.extend(["--normal-phased-vcf", str(phased_vcf)])
    else:
        argv.extend(["--tumor-phased-vcf", str(phased_vcf)])
    argv.extend(
        [
            "--genome-name",
            sample_id,
            "--out-dir-plots",
            str(output_dir),
        ]
    )
    if breakpoints_vcf is None:
        argv.append("--change-point-detection-for-cna")
    else:
        argv.extend(["--breakpoints", str(breakpoints_vcf)])
    return tuple(argv)


def _fingerprint(path: Path) -> FileFingerprint:
    return FileFingerprint(size_bytes=path.stat().st_size, sha256=sha256_file(path))


def _registered_artifact(
    inputs: TumorInputBundle,
    role: CallerInputRole,
    *,
    label: str,
) -> TumorAuxiliaryInputArtifact:
    artifact = inputs.artifact_for(role)
    if artifact is None:
        raise ValueError(f"Wakhan requires registered {label}")
    return artifact


def _validate_file(
    path: Path,
    artifact: TumorAuxiliaryInputArtifact,
    *,
    label: str,
) -> FileFingerprint:
    if not path.is_file():
        raise ValueError(f"Wakhan requires {label}: {path}")
    fingerprint = _fingerprint(path)
    if fingerprint.sha256 != artifact.sha256:
        raise ValueError(f"Wakhan {label} SHA-256 does not match registered input")
    return fingerprint


def _bam_index(path: Path) -> Path:
    for candidate in (Path(f"{path}.bai"), path.with_suffix(".bai"), Path(f"{path}.csi")):
        if candidate.is_file():
            return candidate
    raise ValueError(f"Wakhan requires an index for BAM: {path.name}")


def _validate_common_inputs(
    *,
    tumor_bam: Path,
    phased_vcf: Path,
    reference_fasta: Path,
    sample_id: str,
    inputs: TumorInputBundle,
    policy: WakhanPhasedCnaPolicy,
) -> tuple[FileFingerprint, FileFingerprint, FileFingerprint, FileFingerprint]:
    if not _SAFE_SAMPLE_ID.fullmatch(sample_id) or sample_id in {".", ".."}:
        raise ValueError("Wakhan sample ID is unsafe")
    if inputs.analysis_sample_id != sample_id:
        raise ValueError("Wakhan tumor sample does not match registered analysis sample")
    if inputs.genome_build != policy.genome_build:
        raise ValueError("Wakhan input genome build does not match policy")
    if inputs.reference_id != policy.reference_id:
        raise ValueError("Wakhan input reference ID does not match policy")
    if inputs.reference_sha256 != policy.reference_sha256:
        raise ValueError("Wakhan input reference SHA-256 does not match policy")

    tumor = _registered_artifact(inputs, CallerInputRole.TUMOR_BAM, label="tumor BAM")
    phased = _registered_artifact(
        inputs,
        CallerInputRole.PHASED_VARIANTS,
        label="phased variants",
    )
    if tumor.source_sample_id != sample_id:
        raise ValueError("Wakhan tumor BAM source sample does not match requested sample")
    if policy.mode == "tumor_only" and phased.source_sample_id != sample_id:
        raise ValueError("Wakhan tumor-only mode requires tumor-derived phased variants")
    if policy.mode == "tumor_normal" and phased.source_sample_id == sample_id:
        raise ValueError("Wakhan tumor-normal mode requires normal-derived phased variants")

    tumor_fp = _validate_file(tumor_bam, tumor, label="tumor BAM")
    tumor_index_fp = _fingerprint(_bam_index(tumor_bam))
    phased_fp = _validate_file(phased_vcf, phased, label="phased variants")
    if not reference_fasta.is_file():
        raise ValueError("Wakhan requires the locked reference FASTA")
    if not Path(f"{reference_fasta}.fai").is_file():
        raise ValueError("Wakhan requires the reference FASTA index")
    reference_fp = _fingerprint(reference_fasta)
    if reference_fp.sha256 != policy.reference_sha256:
        raise ValueError("Wakhan reference FASTA SHA-256 does not match policy lock")
    return tumor_fp, tumor_index_fp, phased_fp, reference_fp


def _validate_breakpoints(
    *,
    breakpoints_vcf: Path | None,
    inputs: TumorInputBundle,
    expected_parent_lane_id: str | None,
    expected_parent_lane_sha256: str | None,
) -> tuple[FileFingerprint | None, str | None, str | None]:
    expected_pair = (
        expected_parent_lane_id is not None,
        expected_parent_lane_sha256 is not None,
    )
    if expected_pair[0] != expected_pair[1]:
        raise ValueError("Wakhan expected breakpoint parent lane ID and SHA-256 must be paired")

    registered = inputs.artifact_for(CallerInputRole.SEVERUS_BREAKPOINTS)
    if breakpoints_vcf is None:
        if registered is not None or expected_pair[0]:
            raise ValueError("Wakhan registered parent lane requires explicit breakpoint input")
        return None, None, None

    if registered is None:
        raise ValueError("Wakhan breakpoints require registered Severus breakpoint evidence")
    if not expected_pair[0]:
        raise ValueError("Wakhan breakpoints require an expected parent lane lock")
    assert expected_parent_lane_id is not None
    assert expected_parent_lane_sha256 is not None
    if registered.producer_lane_id != expected_parent_lane_id:
        raise ValueError("Wakhan breakpoint parent lane ID does not match registered producer")
    if registered.producer_lane_sha256 != expected_parent_lane_sha256:
        raise ValueError("Wakhan breakpoint parent lane SHA-256 does not match registered producer")

    fingerprint = _validate_file(
        breakpoints_vcf,
        registered,
        label="Severus breakpoints",
    )
    return fingerprint, expected_parent_lane_id, expected_parent_lane_sha256


def _artifact_role(relative_path: str) -> str:
    path = PurePosixPath(relative_path)
    name = path.name
    if name == "solutions_ranks.tsv":
        return "ranked_solutions"
    if name == "integer_profile.bed":
        return "integer_copy_number_profile_bed"
    if name == "integer_profile.vcf":
        return "integer_copy_number_profile_vcf"
    if name == "subclonal_profile.bed":
        return "subclonal_copy_number_profile_bed"
    if name == "subclonal_profile.vcf":
        return "subclonal_copy_number_profile_vcf"
    if name == "rephased.vcf.gz":
        return "rephased_variants"
    if "coverage_data" in path.parts:
        return "coverage_data"
    if name.startswith("solution_rank_"):
        return "solution_rank_link"
    return "caller_native"


def _media_type(relative_path: str, *, symlink: bool = False) -> str:
    if symlink:
        return "inode/symlink"
    if relative_path.endswith(".tsv"):
        return "text/tab-separated-values"
    if relative_path.endswith(".csv"):
        return "text/csv"
    if relative_path.endswith(".bed"):
        return "text/bed"
    if relative_path.endswith((".vcf", ".vcf.gz")):
        return "text/vcf"
    if relative_path.endswith(".html"):
        return "text/html"
    return "application/octet-stream"


def _symlink_fingerprint(target: str) -> FileFingerprint:
    payload = target.encode("utf-8")
    return FileFingerprint(
        size_bytes=len(payload),
        sha256=hashlib.sha256(payload).hexdigest(),
    )


def _native_artifacts(output_dir: Path) -> list[WakhanNativeArtifact]:
    root = output_dir.resolve()
    artifacts: list[WakhanNativeArtifact] = []
    for path in sorted(output_dir.rglob("*")):
        relative = path.relative_to(output_dir).as_posix()
        if path.is_symlink():
            resolved = path.resolve()
            try:
                target_relative = resolved.relative_to(root).as_posix()
            except ValueError as exc:
                raise ValueError("Wakhan symlink target escapes the output directory") from exc
            artifacts.append(
                WakhanNativeArtifact(
                    relative_path=relative,
                    role=_artifact_role(relative),
                    media_type=_media_type(relative, symlink=True),
                    fingerprint=_symlink_fingerprint(os.readlink(path)),
                    artifact_kind="symlink",
                    symlink_target=target_relative,
                )
            )
            continue
        if not path.is_file():
            continue
        artifacts.append(
            WakhanNativeArtifact(
                relative_path=relative,
                role=_artifact_role(relative),
                media_type=_media_type(relative),
                fingerprint=_fingerprint(path),
            )
        )
    if not artifacts:
        raise ValueError("Wakhan returned success without native artifacts")
    if not any(item.role == "ranked_solutions" for item in artifacts):
        raise ValueError("Wakhan output lacks solutions_ranks.tsv")
    if not any(item.role == "integer_copy_number_profile_bed" for item in artifacts):
        raise ValueError("Wakhan output lacks an integer_profile.bed solution")
    return artifacts


def run_wakhan_phased_cna(
    *,
    tumor_bam: Path,
    phased_vcf: Path,
    reference_fasta: Path,
    sample_id: str,
    output_dir: Path,
    inputs: TumorInputBundle,
    policy: WakhanPhasedCnaPolicy,
    observed_runtime_version: str,
    wakhan_script: Path,
    breakpoints_vcf: Path | None = None,
    expected_breakpoint_parent_lane_id: str | None = None,
    expected_breakpoint_parent_lane_sha256: str | None = None,
    runner: CommandRunner | None = None,
    python_executable: str = "python",
    threads: int = 8,
) -> WakhanReport:
    if output_dir.exists():
        raise ValueError("Refusing to overwrite an existing Wakhan output directory")
    if observed_runtime_version != policy.expected_version:
        raise ValueError(
            f"Wakhan version {observed_runtime_version!r} does not match "
            f"policy lock {policy.expected_version!r}"
        )
    if not wakhan_script.is_file():
        raise ValueError("Wakhan runtime script is missing")

    tumor_fp, tumor_index_fp, phased_fp, reference_fp = _validate_common_inputs(
        tumor_bam=tumor_bam,
        phased_vcf=phased_vcf,
        reference_fasta=reference_fasta,
        sample_id=sample_id,
        inputs=inputs,
        policy=policy,
    )
    breakpoints_fp, parent_id, parent_sha256 = _validate_breakpoints(
        breakpoints_vcf=breakpoints_vcf,
        inputs=inputs,
        expected_parent_lane_id=expected_breakpoint_parent_lane_id,
        expected_parent_lane_sha256=expected_breakpoint_parent_lane_sha256,
    )
    script_fp = _fingerprint(wakhan_script)

    command_runner = runner or SubprocessRunner()
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    promoted = False
    try:
        argv = build_wakhan_argv(
            python_executable=python_executable,
            wakhan_script=wakhan_script,
            tumor_bam=tumor_bam,
            phased_vcf=phased_vcf,
            reference_fasta=reference_fasta,
            sample_id=sample_id,
            output_dir=staged,
            policy=policy,
            threads=threads,
            breakpoints_vcf=breakpoints_vcf,
        )
        result = command_runner.run(argv, timeout_seconds=policy.timeout_seconds)
        if result.returncode != 0:
            diagnostic = result.stderr.strip()[-4000:]
            suffix = f": {diagnostic}" if diagnostic else ""
            raise ValueError(f"Wakhan failed with exit code {result.returncode}{suffix}")

        artifacts = _native_artifacts(staged)
        report = WakhanReport(
            sample_id=sample_id,
            mode=policy.mode,
            genome_build=policy.genome_build,
            status=ModuleRunStatus.COMPLETED,
            policy=policy,
            tool=ToolRecord(
                name="Wakhan",
                version=observed_runtime_version,
                parameters={
                    "mode": policy.mode,
                    "threads": threads,
                    "breakpoints_supplied": breakpoints_vcf is not None,
                    "change_point_detection_for_cna": breakpoints_vcf is None,
                    "use_sv_haplotypes": False,
                    "real_tool_qualified": policy.real_tool_qualified,
                    "analytical_validation": policy.analytical_validation,
                },
            ),
            real_tool_qualified=policy.real_tool_qualified,
            tumor_bam_fingerprint=tumor_fp,
            tumor_bam_index_fingerprint=tumor_index_fp,
            phased_vcf_fingerprint=phased_fp,
            reference_fingerprint=reference_fp,
            runtime_script_fingerprint=script_fp,
            breakpoints_fingerprint=breakpoints_fp,
            breakpoint_parent_lane_id=parent_id,
            breakpoint_parent_lane_sha256=parent_sha256,
            native_artifacts=artifacts,
            warnings=[policy.note],
            limitations=[
                (
                    "This adapter retains Wakhan CNA outputs but does not infer clinical "
                    "reportability or biological truth."
                ),
                (
                    "Wakhan integer_profile BED/VCF files remain caller-native evidence until "
                    "their exact schema is separately normalized and validated."
                ),
                (
                    "Severus breakpoints are optional; when absent, the adapter explicitly "
                    "requests Wakhan change-point detection for CNA segmentation."
                ),
                (
                    "The adapter does not enable --use-sv-haplotypes unless phased Severus "
                    "breakpoints are produced by a separately qualified workflow."
                ),
            ],
        )
        os.replace(staged, output_dir)
        promoted = True
        return report
    finally:
        if not promoted:
            shutil.rmtree(staged, ignore_errors=True)
