from __future__ import annotations

import gzip
import math
import os
import re
import shutil
import tempfile
from collections import Counter
from collections.abc import Iterator
from pathlib import Path, PurePosixPath
from typing import Literal

from pydantic import Field, ValidationError, field_validator, model_validator

from ..breakends import BreakendParseError, resolve_breakend
from ..execution import CommandRunner, SubprocessRunner
from ..models import (
    EventType,
    Evidence,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    StrictModel,
    ToolRecord,
)
from ..multicaller_contracts import CallerInputRole
from ..reference import sha256_file
from ..tumor_inputs import TumorAuxiliaryInputArtifact, TumorInputBundle

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")
_SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_SHA256 = r"^[0-9a-f]{64}$"
_EVENT_TYPES = {
    "DEL": EventType.DELETION,
    "DUP": EventType.DUPLICATION,
    "INV": EventType.INVERSION,
    "INS": EventType.INSERTION,
    "BND": EventType.TRANSLOCATION,
    "TRA": EventType.TRANSLOCATION,
}


class SeverusPairedPolicy(StrictModel):
    """Versioned technical contract for paired tumor-normal Severus execution."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    mode: Literal["paired"] = "paired"
    expected_version: str = Field(pattern=r"^\d+\.\d+(?:\.\d+)?$")
    real_tool_qualified: bool = False
    analytical_validation: Literal["not_validated"] = "not_validated"
    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    reference_sha256: str = Field(pattern=_SHA256)
    min_support: int = Field(ge=1)
    min_sv_size: int = Field(ge=1)
    minimum_mapping_quality: int = Field(ge=0)
    tumor_in_normal_ratio: float = Field(ge=0, le=1)
    control_coverage_threshold: int = Field(ge=0)
    timeout_seconds: int = Field(ge=60)
    output_read_ids: Literal[False] = False
    write_alignments: Literal[False] = False
    pon_mode: Literal[False] = False
    note: str = Field(min_length=12)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def pinned_runtime_version(self) -> SeverusPairedPolicy:
        if self.expected_version != "1.7":
            raise ValueError("Severus paired expected_version must match the pinned 1.7 contract")
        return self


class SeverusNativeArtifact(StrictModel):
    relative_path: str = Field(min_length=1)
    role: str = Field(min_length=1)
    media_type: str = Field(min_length=1)
    fingerprint: FileFingerprint
    research_only: Literal[True] = True

    @field_validator("relative_path")
    @classmethod
    def safe_relative_path(cls, value: str) -> str:
        if "\\" in value:
            raise ValueError("Severus artifact paths must use POSIX separators")
        path = PurePosixPath(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Severus artifact path must remain inside the output directory")
        return value


class SeverusPairedReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    tumor_sample_id: str = Field(min_length=1)
    normal_sample_id: str = Field(min_length=1)
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: SeverusPairedPolicy
    events: list[GenomicEvent] = Field(default_factory=list)
    raw_record_count: int = Field(ge=0)
    accepted_record_count: int = Field(ge=0)
    rejected_record_count: int = Field(ge=0)
    rejection_counts: dict[str, int] = Field(default_factory=dict)
    tool: ToolRecord
    tumor_bam_fingerprint: FileFingerprint
    tumor_bam_index_fingerprint: FileFingerprint
    normal_bam_fingerprint: FileFingerprint
    normal_bam_index_fingerprint: FileFingerprint
    somatic_vcf_fingerprint: FileFingerprint
    native_artifacts: list[SeverusNativeArtifact] = Field(min_length=1)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_report(self) -> SeverusPairedReport:
        if self.accepted_record_count != len(self.events):
            raise ValueError("Severus accepted count must equal normalized events")
        if self.accepted_record_count + self.rejected_record_count != self.raw_record_count:
            raise ValueError("Severus accepted and rejected counts must cover raw records")
        if sum(self.rejection_counts.values()) != self.rejected_record_count:
            raise ValueError("Severus rejection counts must sum to rejected_record_count")
        expected = ModuleRunStatus.COMPLETED if self.events else ModuleRunStatus.NO_CALL
        if self.status != expected:
            raise ValueError("Severus status is inconsistent with normalized events")
        paths = [item.relative_path for item in self.native_artifacts]
        if len(paths) != len(set(paths)):
            raise ValueError("Severus native artifact paths must be unique")
        vcf_artifacts = [
            item
            for item in self.native_artifacts
            if item.relative_path
            in {
                "somatic_SVs/severus_somatic.vcf",
                "somatic_SVs/severus_somatic.vcf.gz",
            }
        ]
        if len(vcf_artifacts) != 1:
            raise ValueError("Severus native artifacts must contain exactly one somatic VCF")
        if vcf_artifacts[0].fingerprint != self.somatic_vcf_fingerprint:
            raise ValueError("Severus somatic VCF fingerprint must match native artifact")
        return self


class _RejectedRecord(ValueError):
    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def severus_version(text: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    return text.splitlines()[0].strip()[:80] if text.strip() else "unknown"


def build_severus_paired_argv(
    *,
    severus: str,
    tumor_bam: Path,
    normal_bam: Path,
    tumor_sample_id: str,
    normal_sample_id: str,
    output_dir: Path,
    policy: SeverusPairedPolicy,
    threads: int = 8,
) -> tuple[str, ...]:
    if threads < 1:
        raise ValueError("Severus threads must be at least 1")
    return (
        severus,
        "--target-bam",
        str(tumor_bam),
        "--control-bam",
        str(normal_bam),
        "--out-dir",
        str(output_dir),
        "-t",
        str(threads),
        "--target-sample",
        tumor_sample_id,
        "--control-sample",
        normal_sample_id,
        "--min-support",
        str(policy.min_support),
        "--min-sv-size",
        str(policy.min_sv_size),
        "--min-mapq",
        str(policy.minimum_mapping_quality),
        "--TIN-ratio",
        str(policy.tumor_in_normal_ratio),
        "--control-cov-thr",
        str(policy.control_coverage_threshold),
    )


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
        raise ValueError(f"Severus paired requires registered {label}")
    return artifact


def _bam_index_path(path: Path) -> Path:
    candidates = (Path(f"{path}.bai"), path.with_suffix(".bai"), Path(f"{path}.csi"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError(f"Severus requires an index for BAM: {path.name}")


def _validate_registered_file(
    path: Path,
    artifact: TumorAuxiliaryInputArtifact,
    *,
    label: str,
) -> FileFingerprint:
    if not path.is_file():
        raise ValueError(f"Severus requires {label}: {path}")
    fingerprint = _fingerprint(path)
    if fingerprint.sha256 != artifact.sha256:
        raise ValueError(f"Severus {label} SHA-256 does not match registered input")
    return fingerprint


def _validate_input_identity(
    *,
    inputs: TumorInputBundle,
    policy: SeverusPairedPolicy,
    tumor_sample_id: str,
    normal_sample_id: str,
    tumor_bam: Path,
    normal_bam: Path,
) -> tuple[
    FileFingerprint,
    FileFingerprint,
    FileFingerprint,
    FileFingerprint,
]:
    if inputs.analysis_sample_id != tumor_sample_id:
        raise ValueError("Severus tumor sample does not match registered analysis sample")
    if inputs.genome_build != policy.genome_build:
        raise ValueError("Severus input genome build does not match policy")
    if inputs.reference_id != policy.reference_id:
        raise ValueError("Severus input reference ID does not match policy")
    if inputs.reference_sha256 != policy.reference_sha256:
        raise ValueError("Severus input reference SHA-256 does not match policy")

    tumor = _registered_artifact(inputs, CallerInputRole.TUMOR_BAM, label="tumor BAM")
    normal = _registered_artifact(
        inputs,
        CallerInputRole.MATCHED_NORMAL_BAM,
        label="matched normal BAM",
    )
    if tumor.source_sample_id != tumor_sample_id:
        raise ValueError("Severus tumor BAM source sample does not match requested tumor sample")
    if normal.source_sample_id != normal_sample_id:
        raise ValueError(
            "Severus matched normal source sample does not match requested normal sample"
        )

    tumor_fp = _validate_registered_file(tumor_bam, tumor, label="tumor BAM")
    normal_fp = _validate_registered_file(normal_bam, normal, label="matched normal BAM")
    tumor_index_fp = _fingerprint(_bam_index_path(tumor_bam))
    normal_index_fp = _fingerprint(_bam_index_path(normal_bam))
    return tumor_fp, tumor_index_fp, normal_fp, normal_index_fp


def _open_vcf(path: Path) -> Iterator[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield from handle


def _parse_info(raw: str) -> dict[str, str | bool]:
    result: dict[str, str | bool] = {}
    if raw in {"", "."}:
        return result
    for token in raw.split(";"):
        if not token:
            continue
        key, separator, value = token.partition("=")
        if not key:
            raise _RejectedRecord("malformed_info")
        result[key] = value if separator else True
    return result


def _text(info: dict[str, str | bool], key: str) -> str | None:
    value = info.get(key)
    return value if isinstance(value, str) and value not in {"", "."} else None


def _quality(raw: str) -> float | None:
    if raw in {"", "."}:
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        raise _RejectedRecord("malformed_quality") from exc
    if not math.isfinite(value) or value < 0:
        raise _RejectedRecord("malformed_quality")
    return value


def _event_loci(
    *,
    chromosome: str,
    position: int,
    alternate: str,
    event_type: EventType,
    info: dict[str, str | bool],
) -> tuple[Locus, Locus | None, int | None]:
    try:
        if event_type == EventType.TRANSLOCATION:
            try:
                resolved = resolve_breakend(
                    alternate,
                    declared_chromosome=info.get("CHR2"),
                    declared_position=info.get("END"),
                )
            except BreakendParseError as exc:
                raise _RejectedRecord(exc.reason) from exc
            return (
                Locus(chromosome=chromosome, start=position - 1, end=position),
                Locus(
                    chromosome=resolved.mate_chromosome,
                    start=resolved.mate_position_0based,
                    end=resolved.mate_position_0based + 1,
                ),
                None,
            )

        start = position - 1
        raw_end = _text(info, "END")
        if event_type == EventType.INSERTION:
            return Locus(chromosome=chromosome, start=start, end=start + 1), None, 1
        if raw_end is None:
            raise _RejectedRecord("missing_end")
        end = int(raw_end)
        if end <= start:
            raise _RejectedRecord("invalid_interval")
        return Locus(chromosome=chromosome, start=start, end=end), None, end - start
    except (ValidationError, ValueError) as exc:
        if isinstance(exc, _RejectedRecord):
            raise
        raise _RejectedRecord("invalid_locus") from exc


def _normalize_record(
    fields: list[str],
    *,
    record_number: int,
    caller_version: str,
) -> GenomicEvent:
    if len(fields) < 8:
        raise _RejectedRecord("malformed_record")
    chromosome, raw_position, native_id, _ref, alternate, raw_quality, raw_filter, raw_info = (
        fields[:8]
    )
    try:
        position = int(raw_position)
    except ValueError as exc:
        raise _RejectedRecord("malformed_position") from exc
    if position < 1:
        raise _RejectedRecord("malformed_position")

    info = _parse_info(raw_info)
    raw_type = (_text(info, "SVTYPE") or "").upper().split(":", maxsplit=1)[0]
    event_type = _EVENT_TYPES.get(raw_type)
    if event_type is None:
        raise _RejectedRecord("unsupported_svtype")
    primary, secondary, length = _event_loci(
        chromosome=chromosome,
        position=position,
        alternate=alternate,
        event_type=event_type,
        info=info,
    )
    filters = [] if raw_filter in {"", "."} else raw_filter.split(";")
    notes = [
        "Normalized from the paired Severus somatic VCF into 0-based half-open coordinates.",
        "Candidate only; analytical validation and expert review remain required.",
    ]
    if native_id not in {"", "."}:
        notes.append(f"Severus native candidate ID: {native_id}")
    for key in ("DETAILED_TYPE", "CLUSTERID", "CLUSTER_ID", "HP", "PHASESET_ID"):
        value = _text(info, key)
        if value is not None:
            notes.append(f"Severus native {key}={value}")
    return GenomicEvent(
        event_id=f"SEVERUS-{record_number:06d}",
        event_type=event_type,
        primary=primary,
        secondary=secondary,
        length_bp=length,
        evidence=[
            Evidence(
                caller="Severus",
                caller_version=caller_version,
                quality=_quality(raw_quality),
                filters=filters,
            )
        ],
        reportable=False,
        notes=notes,
    )


def _normalize_somatic_vcf(
    path: Path,
    *,
    caller_version: str,
) -> tuple[list[GenomicEvent], int, Counter[str]]:
    saw_fileformat = False
    saw_source = False
    saw_columns = False
    raw_count = 0
    events: list[GenomicEvent] = []
    rejections: Counter[str] = Counter()
    for raw_line in _open_vcf(path):
        line = raw_line.rstrip("\r\n")
        if line.startswith("##fileformat=VCF"):
            saw_fileformat = True
            continue
        if line.startswith("##source=Severus"):
            saw_source = True
            continue
        if line.startswith("#CHROM"):
            if saw_columns:
                raise ValueError("Severus output contains duplicate VCF column headers")
            columns = line.split("\t")
            expected = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
            if columns[:8] != expected:
                raise ValueError("Severus output has an unexpected VCF header")
            saw_columns = True
            continue
        if not line or line.startswith("#"):
            continue
        raw_count += 1
        try:
            events.append(
                _normalize_record(
                    line.split("\t"),
                    record_number=raw_count,
                    caller_version=caller_version,
                )
            )
        except _RejectedRecord as exc:
            rejections[exc.reason] += 1
    if not saw_fileformat or not saw_columns:
        raise ValueError("Severus output is not a complete VCF document")
    if not saw_source:
        raise ValueError("Severus output lacks the ##source=Severus declaration")
    return events, raw_count, rejections


def _artifact_role(relative_path: str) -> str:
    name = PurePosixPath(relative_path).name
    if name.startswith("severus_somatic.vcf"):
        return "somatic_vcf"
    if "breakpoint" in name and name.endswith(".csv"):
        return "breakpoint_table"
    if "cluster" in name and name.endswith(".tsv"):
        return "breakpoint_cluster"
    if relative_path.endswith(".html"):
        return "breakpoint_graph"
    if name == "read_qual.txt":
        return "read_quality_summary"
    if name.endswith(".log"):
        return "caller_log"
    return "caller_native"


def _media_type(relative_path: str) -> str:
    if relative_path.endswith((".vcf", ".vcf.gz")):
        return "text/vcf"
    if relative_path.endswith(".csv"):
        return "text/csv"
    if relative_path.endswith(".tsv"):
        return "text/tab-separated-values"
    if relative_path.endswith(".html"):
        return "text/html"
    if relative_path.endswith((".txt", ".log")):
        return "text/plain"
    return "application/octet-stream"


def _native_artifacts(output_dir: Path) -> list[SeverusNativeArtifact]:
    artifacts: list[SeverusNativeArtifact] = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.is_symlink():
            raise ValueError("Severus native output must not contain symbolic links")
        relative = path.relative_to(output_dir).as_posix()
        artifacts.append(
            SeverusNativeArtifact(
                relative_path=relative,
                role=_artifact_role(relative),
                media_type=_media_type(relative),
                fingerprint=_fingerprint(path),
            )
        )
    if not artifacts:
        raise ValueError("Severus returned success without native artifacts")
    return artifacts


def _somatic_vcf(output_dir: Path) -> Path:
    candidates = [
        candidate
        for candidate in (
            output_dir / "somatic_SVs" / "severus_somatic.vcf",
            output_dir / "somatic_SVs" / "severus_somatic.vcf.gz",
        )
        if candidate.is_file()
    ]
    if len(candidates) != 1:
        raise ValueError("Severus paired must produce exactly one somatic VCF")
    return candidates[0]


def run_severus_paired(
    *,
    tumor_bam: Path,
    normal_bam: Path,
    tumor_sample_id: str,
    normal_sample_id: str,
    output_dir: Path,
    inputs: TumorInputBundle,
    policy: SeverusPairedPolicy,
    runner: CommandRunner | None = None,
    severus: str = "severus",
    threads: int = 8,
) -> SeverusPairedReport:
    """Run paired Severus without tumor-only fallback or sensitive optional output flags."""
    if not _SAFE_SAMPLE_ID.fullmatch(tumor_sample_id) or tumor_sample_id in {".", ".."}:
        raise ValueError("Severus tumor sample ID is unsafe")
    if not _SAFE_SAMPLE_ID.fullmatch(normal_sample_id) or normal_sample_id in {".", ".."}:
        raise ValueError("Severus normal sample ID is unsafe")
    if threads < 1:
        raise ValueError("Severus threads must be at least 1")
    if output_dir.exists():
        raise ValueError("Refusing to overwrite an existing Severus output directory")

    tumor_fp, tumor_index_fp, normal_fp, normal_index_fp = _validate_input_identity(
        inputs=inputs,
        policy=policy,
        tumor_sample_id=tumor_sample_id,
        normal_sample_id=normal_sample_id,
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
    )

    command_runner = runner or SubprocessRunner()
    probe = command_runner.run([severus, "--version"], timeout_seconds=30)
    if probe.returncode != 0:
        raise ValueError(f"Severus version probe failed with exit code {probe.returncode}")
    version = severus_version(f"{probe.stdout}\n{probe.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"Severus version {version!r} does not match policy lock {policy.expected_version!r}"
        )

    parameters: dict[str, str | int | float | bool] = {
        "mode": policy.mode,
        "min_support": policy.min_support,
        "min_sv_size": policy.min_sv_size,
        "minimum_mapping_quality": policy.minimum_mapping_quality,
        "tumor_in_normal_ratio": policy.tumor_in_normal_ratio,
        "control_coverage_threshold": policy.control_coverage_threshold,
        "output_read_ids": False,
        "write_alignments": False,
        "pon_mode": False,
        "real_tool_qualified": policy.real_tool_qualified,
        "analytical_validation": policy.analytical_validation,
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    promoted = False
    try:
        argv = build_severus_paired_argv(
            severus=severus,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            tumor_sample_id=tumor_sample_id,
            normal_sample_id=normal_sample_id,
            output_dir=staged,
            policy=policy,
            threads=threads,
        )
        result = command_runner.run(argv, timeout_seconds=policy.timeout_seconds)
        if result.returncode != 0:
            diagnostic = result.stderr.strip()[-4000:]
            suffix = f": {diagnostic}" if diagnostic else ""
            raise ValueError(f"Severus failed with exit code {result.returncode}{suffix}")

        somatic_vcf = _somatic_vcf(staged)
        artifacts = _native_artifacts(staged)
        events, raw_count, rejections = _normalize_somatic_vcf(
            somatic_vcf,
            caller_version=version,
        )
        warnings = [policy.note]
        if rejections:
            warnings.append("One or more Severus records were rejected during normalization.")
        if not events:
            warnings.append(
                "NO_CALL means no paired Severus somatic record was retained; "
                "it is not a biological negative result."
            )
        report = SeverusPairedReport(
            tumor_sample_id=tumor_sample_id,
            normal_sample_id=normal_sample_id,
            genome_build=policy.genome_build,
            status=ModuleRunStatus.COMPLETED if events else ModuleRunStatus.NO_CALL,
            policy=policy,
            events=events,
            raw_record_count=raw_count,
            accepted_record_count=len(events),
            rejected_record_count=sum(rejections.values()),
            rejection_counts=dict(sorted(rejections.items())),
            tool=ToolRecord(name="Severus", version=version, parameters=parameters),
            tumor_bam_fingerprint=tumor_fp,
            tumor_bam_index_fingerprint=tumor_index_fp,
            normal_bam_fingerprint=normal_fp,
            normal_bam_index_fingerprint=normal_index_fp,
            somatic_vcf_fingerprint=_fingerprint(somatic_vcf),
            native_artifacts=artifacts,
            warnings=warnings,
            limitations=[
                (
                    "This adapter represents paired tumor-normal Severus only; "
                    "missing control input fails closed rather than switching to tumor-only."
                ),
                (
                    "Breakpoint graphs, clusters and caller-native tables are retained as "
                    "technical evidence and are not biological truth or clinical release."
                ),
                (
                    "PoN mode, read-ID output and alignment writing are disabled in this "
                    "research-only paired adapter."
                ),
            ],
        )
        os.replace(staged, output_dir)
        promoted = True
        return report
    finally:
        if not promoted:
            shutil.rmtree(staged, ignore_errors=True)
