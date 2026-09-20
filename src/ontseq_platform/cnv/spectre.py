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
from ..reference import sha256_file

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")
_COORDINATE_CONTRACT: Literal["spectre-0.2.1-mosdepth-zero-based-native-v1"] = (
    "spectre-0.2.1-mosdepth-zero-based-native-v1"
)
_EVENT_TYPES = {
    "DEL": EventType.DELETION,
    "DUP": EventType.DUPLICATION,
}
_SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_RUNTIME_PACKAGE: Literal["spectre-cnv==0.2.1"] = "spectre-cnv==0.2.1"
_RUNTIME_WHEEL_SHA256: Literal[
    "c50998798a33b22455f3681a69c5c61e82df3fe049224da6d2e54b068bc35510"
] = "c50998798a33b22455f3681a69c5c61e82df3fe049224da6d2e54b068bc35510"
_RUNTIME_SOURCE_COMMIT: Literal["b80fd32bd189c1689fc1031ff96433f2e79d1bed"] = (
    "b80fd32bd189c1689fc1031ff96433f2e79d1bed"
)


class SpectrePolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    mode: Literal["depth_only"] = "depth_only"
    expected_version: str = Field(pattern=r"^\d+\.\d+(?:\.\d+)?$")
    runtime_package: Literal["spectre-cnv==0.2.1"] = _RUNTIME_PACKAGE
    runtime_wheel_sha256: Literal[
        "c50998798a33b22455f3681a69c5c61e82df3fe049224da6d2e54b068bc35510"
    ] = _RUNTIME_WHEEL_SHA256
    runtime_source_commit: Literal["b80fd32bd189c1689fc1031ff96433f2e79d1bed"] = (
        _RUNTIME_SOURCE_COMMIT
    )
    real_tool_qualified: bool = False
    analytical_validation: Literal["not_validated"] = "not_validated"
    genome_build: GenomeBuild
    coordinate_contract: Literal["spectre-0.2.1-mosdepth-zero-based-native-v1"] = (
        _COORDINATE_CONTRACT
    )
    mosdepth_bin_size_bp: int = Field(gt=0)
    minimum_mapping_quality: int = Field(ge=0)
    minimum_cnv_length_bp: int = Field(gt=0)
    ploidy: int = Field(ge=1)
    pass_only: bool = True
    timeout_seconds: int = Field(ge=60)
    note: str = Field(min_length=12)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def version_matches_pinned_runtime(self) -> SpectrePolicy:
        if self.expected_version != "0.2.1":
            raise ValueError("Spectre expected_version must match the pinned 0.2.1 runtime")
        return self


class SpectreNativeArtifact(StrictModel):
    relative_path: str = Field(min_length=1)
    media_type: str = Field(min_length=1)
    fingerprint: FileFingerprint
    research_only: Literal[True] = True

    @field_validator("relative_path")
    @classmethod
    def safe_relative_path(cls, value: str) -> str:
        if "\\" in value:
            raise ValueError("Spectre artifact paths must use POSIX separators")
        path = PurePosixPath(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Spectre artifact path must remain inside the output directory")
        return value


class SpectreCallReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str = Field(min_length=1)
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: SpectrePolicy
    coordinate_contract: Literal["spectre-0.2.1-mosdepth-zero-based-native-v1"] = (
        _COORDINATE_CONTRACT
    )
    events: list[GenomicEvent] = Field(default_factory=list)
    raw_record_count: int = Field(ge=0)
    accepted_record_count: int = Field(ge=0)
    rejected_record_count: int = Field(ge=0)
    rejection_counts: dict[str, int] = Field(default_factory=dict)
    tool: ToolRecord
    vcf_fingerprint: FileFingerprint
    coverage_fingerprint: FileFingerprint | None = None
    coverage_index_fingerprint: FileFingerprint | None = None
    reference_fingerprint: FileFingerprint | None = None
    reference_index_fingerprint: FileFingerprint | None = None
    metadata_input_fingerprint: FileFingerprint | None = None
    blacklist_fingerprint: FileFingerprint | None = None
    native_artifacts: list[SpectreNativeArtifact] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_counts_and_status(self) -> SpectreCallReport:
        if self.accepted_record_count != len(self.events):
            raise ValueError("Spectre accepted count must equal retained normalized events")
        if self.accepted_record_count + self.rejected_record_count != self.raw_record_count:
            raise ValueError("Spectre accepted and rejected counts must cover raw records")
        if sum(self.rejection_counts.values()) != self.rejected_record_count:
            raise ValueError("Spectre rejection counts must sum to rejected_record_count")
        expected = ModuleRunStatus.COMPLETED if self.events else ModuleRunStatus.NO_CALL
        if self.status != expected:
            raise ValueError("Spectre status is inconsistent with normalized events")
        paths = [artifact.relative_path for artifact in self.native_artifacts]
        if len(paths) != len(set(paths)):
            raise ValueError("Spectre native artifact paths must be unique")
        if self.native_artifacts:
            expected_vcfs = {f"{self.sample_id}.vcf", f"{self.sample_id}.vcf.gz"}
            native_vcfs = [
                artifact
                for artifact in self.native_artifacts
                if artifact.relative_path in expected_vcfs
            ]
            if len(native_vcfs) != 1:
                raise ValueError("Spectre native artifacts must contain exactly one sample VCF")
            if native_vcfs[0].fingerprint != self.vcf_fingerprint:
                raise ValueError("Spectre VCF fingerprint must match the retained native artifact")
        return self


class _RejectedRecord(ValueError):
    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def spectre_version(text: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    return text.splitlines()[0].strip()[:80] if text.strip() else "unknown"


def build_spectre_depth_only_argv(
    *,
    spectre: str,
    coverage_path: Path,
    sample_id: str,
    output_dir: Path,
    reference_fasta: Path,
    policy: SpectrePolicy,
    threads: int = 1,
    metadata_path: Path | None = None,
    blacklist_path: Path | None = None,
) -> tuple[str, ...]:
    if threads < 1:
        raise ValueError("Spectre threads must be at least 1")
    argv = [
        spectre,
        "CNVCaller",
        "--coverage",
        str(coverage_path),
        "--sample-id",
        sample_id,
        "--output-dir",
        str(output_dir),
        "--reference",
        str(reference_fasta),
        "--ploidy",
        str(policy.ploidy),
        "--min-cnv-len",
        str(policy.minimum_cnv_length_bp),
        "--threads",
        str(threads),
    ]
    if metadata_path is not None:
        argv.extend(["--metadata", str(metadata_path)])
    if blacklist_path is not None:
        argv.extend(["--blacklist", str(blacklist_path)])
    return tuple(argv)


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
    for item in raw.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if not key:
            raise _RejectedRecord("malformed_info")
        result[key] = value if separator else True
    return result


def _required_text(info: dict[str, str | bool], key: str, *, reason: str) -> str:
    value = info.get(key)
    if not isinstance(value, str) or not value:
        raise _RejectedRecord(reason)
    return value


def _native_integer(info: dict[str, str | bool], key: str, *, reason: str) -> int:
    raw = _required_text(info, key, reason=reason)
    try:
        return int(raw)
    except ValueError as exc:
        raise _RejectedRecord(reason) from exc


def _copy_number(info: dict[str, str | bool]) -> float:
    raw = _required_text(info, "CN", reason="missing_copy_number")
    try:
        value = float(raw)
    except ValueError as exc:
        raise _RejectedRecord("malformed_copy_number") from exc
    if not math.isfinite(value):
        raise _RejectedRecord("nonfinite_copy_number")
    if value < 0:
        raise _RejectedRecord("negative_copy_number")
    return value


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


def _normalize_record(
    fields: list[str],
    *,
    record_number: int,
    policy: SpectrePolicy,
    caller_version: str,
) -> GenomicEvent:
    if len(fields) != 10:
        raise _RejectedRecord("malformed_record")
    chromosome, raw_start, native_id, _ref, alternate, raw_quality, raw_filter, raw_info = fields[
        :8
    ]
    try:
        start = int(raw_start)
    except ValueError as exc:
        raise _RejectedRecord("malformed_start") from exc
    if start < 0:
        raise _RejectedRecord("malformed_start")

    if raw_filter == "":
        raise _RejectedRecord("malformed_filter")
    filters = raw_filter.split(";") if raw_filter != "." else []
    if policy.pass_only and raw_filter not in {".", "PASS"}:
        raise _RejectedRecord("filter_not_pass")

    info = _parse_info(raw_info)
    raw_type = _required_text(info, "SVTYPE", reason="missing_svtype").upper()
    event_type = _EVENT_TYPES.get(raw_type)
    if event_type is None:
        raise _RejectedRecord("unsupported_svtype")

    expected_alt = f"<{raw_type}>"
    if alternate not in {expected_alt, "."}:
        raise _RejectedRecord("svtype_alt_mismatch")

    end = _native_integer(info, "END", reason="missing_or_invalid_end")
    if end <= start:
        raise _RejectedRecord("invalid_native_interval")
    native_length = _native_integer(info, "SVLEN", reason="missing_or_invalid_svlen")
    if native_length != end - start or native_length < 1:
        raise _RejectedRecord("native_length_mismatch")
    if native_length < policy.minimum_cnv_length_bp:
        raise _RejectedRecord("cnv_length_below_policy")

    cn = _copy_number(info)
    quality = _quality(raw_quality)
    try:
        locus = Locus(chromosome=chromosome, start=start, end=end)
    except ValidationError as exc:
        raise _RejectedRecord("invalid_primary_locus") from exc

    notes = [
        (
            "Spectre 0.2.1 native coordinates are preserved under "
            f"{policy.coordinate_contract}; POS is not shifted as generic VCF 1-based input."
        ),
        "Candidate only; analytical validation and expert review remain required.",
    ]
    if native_id not in {"", "."}:
        notes.append(f"Spectre native candidate ID: {native_id}")
    return GenomicEvent(
        event_id=f"SPECTRE-{record_number:06d}",
        event_type=event_type,
        primary=locus,
        length_bp=native_length,
        copy_number=cn,
        evidence=[
            Evidence(
                caller="Spectre",
                caller_version=caller_version,
                quality=quality,
                filters=filters,
            )
        ],
        reportable=False,
        notes=notes,
    )


def normalize_spectre_vcf(
    path: Path,
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    policy: SpectrePolicy,
    tool: ToolRecord,
    native_artifacts: list[SpectreNativeArtifact] | None = None,
    coverage_fingerprint: FileFingerprint | None = None,
    coverage_index_fingerprint: FileFingerprint | None = None,
    reference_fingerprint: FileFingerprint | None = None,
    reference_index_fingerprint: FileFingerprint | None = None,
    metadata_input_fingerprint: FileFingerprint | None = None,
    blacklist_fingerprint: FileFingerprint | None = None,
) -> SpectreCallReport:
    if not path.is_file():
        raise ValueError("Spectre VCF is missing or unreadable")
    if genome_build != policy.genome_build:
        raise ValueError("Spectre policy and requested genome build do not match")
    if tool.version != policy.expected_version:
        raise ValueError(
            f"Spectre version {tool.version!r} does not match policy lock "
            f"{policy.expected_version!r}"
        )

    saw_fileformat = False
    saw_spectre_source = False
    saw_columns = False
    raw_count = 0
    events: list[GenomicEvent] = []
    rejections: Counter[str] = Counter()
    for raw_line in _open_vcf(path):
        line = raw_line.rstrip("\r\n")
        if line.startswith("##fileformat=VCF"):
            saw_fileformat = True
            continue
        if line == "##source=Spectre":
            saw_spectre_source = True
            continue
        if line.startswith("#CHROM"):
            if saw_columns:
                raise ValueError("Spectre output contains duplicate VCF column headers")
            columns = line.split("\t")
            expected = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            if len(columns) != 10 or columns[:9] != expected:
                raise ValueError("Spectre output has an unexpected single-sample VCF header")
            if columns[9] != sample_id:
                raise ValueError("Spectre VCF sample header does not match requested sample")
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
                    policy=policy,
                    caller_version=tool.version,
                )
            )
        except _RejectedRecord as exc:
            rejections[exc.reason] += 1

    if not saw_fileformat or not saw_columns:
        raise ValueError("Spectre output is not a complete VCF document")
    if not saw_spectre_source:
        raise ValueError("Spectre output lacks the pinned ##source=Spectre declaration")

    warnings = [policy.note]
    if rejections:
        warnings.append("One or more Spectre records were rejected by normalization policy.")
    if not events:
        warnings.append(
            "NO_CALL means no Spectre record passed this technical policy; "
            "it is not a biological negative result."
        )

    return SpectreCallReport(
        sample_id=sample_id,
        genome_build=genome_build,
        status=ModuleRunStatus.COMPLETED if events else ModuleRunStatus.NO_CALL,
        policy=policy,
        events=events,
        raw_record_count=raw_count,
        accepted_record_count=len(events),
        rejected_record_count=sum(rejections.values()),
        rejection_counts=dict(sorted(rejections.items())),
        tool=tool,
        vcf_fingerprint=FileFingerprint(
            size_bytes=path.stat().st_size,
            sha256=sha256_file(path),
        ),
        coverage_fingerprint=coverage_fingerprint,
        coverage_index_fingerprint=coverage_index_fingerprint,
        reference_fingerprint=reference_fingerprint,
        reference_index_fingerprint=reference_index_fingerprint,
        metadata_input_fingerprint=metadata_input_fingerprint,
        blacklist_fingerprint=blacklist_fingerprint,
        native_artifacts=list(native_artifacts or []),
        warnings=warnings,
        limitations=[
            (
                "Spectre 0.2.1 writes Mosdepth-derived zero-based candidate positions directly "
                "into its VCF fields; ONTSeq preserves that versioned native coordinate contract."
            ),
            (
                "Depth-only mode does not consume SNV, Sniffles/SNFJ or population evidence; "
                "supported modes require separately registered lanes."
            ),
        ],
    )


def _fingerprint(path: Path) -> FileFingerprint:
    return FileFingerprint(size_bytes=path.stat().st_size, sha256=sha256_file(path))


def _require_file(path: Path, *, label: str) -> FileFingerprint:
    if not path.is_file():
        raise ValueError(f"Spectre requires {label}: {path}")
    return _fingerprint(path)


def _artifact_media_type(relative_path: str) -> str:
    if relative_path.endswith((".vcf", ".vcf.gz")):
        return "text/vcf"
    if relative_path.endswith(".bed.gz"):
        return "application/gzip"
    if relative_path.endswith(".png"):
        return "image/png"
    if relative_path.endswith(".mdr"):
        return "text/plain"
    return "application/octet-stream"


def _native_artifacts(output_dir: Path) -> list[SpectreNativeArtifact]:
    artifacts: list[SpectreNativeArtifact] = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.is_symlink():
            raise ValueError("Spectre native output must not contain symbolic links")
        relative_path = path.relative_to(output_dir).as_posix()
        artifacts.append(
            SpectreNativeArtifact(
                relative_path=relative_path,
                media_type=_artifact_media_type(relative_path),
                fingerprint=_fingerprint(path),
            )
        )
    if not artifacts:
        raise ValueError("Spectre returned success without any native artifacts")
    return artifacts


def run_spectre_depth_only(
    *,
    coverage_path: Path,
    sample_id: str,
    output_dir: Path,
    reference_fasta: Path,
    genome_build: GenomeBuild,
    policy: SpectrePolicy,
    runner: CommandRunner | None = None,
    spectre: str = "spectre",
    threads: int = 1,
    metadata_path: Path | None = None,
    blacklist_path: Path | None = None,
) -> SpectreCallReport:
    """Run pinned Spectre in independent depth-only mode and retain every native artifact."""
    if genome_build != policy.genome_build:
        raise ValueError("Spectre policy and requested genome build do not match")
    if not _SAFE_SAMPLE_ID.fullmatch(sample_id) or sample_id in {".", ".."}:
        raise ValueError("Spectre sample ID is unsafe for caller-native output paths")
    if threads < 1:
        raise ValueError("Spectre threads must be at least 1")
    if output_dir.exists():
        raise ValueError("Refusing to overwrite an existing Spectre output directory")

    coverage_fingerprint = _require_file(coverage_path, label="Mosdepth regions BED.gz")
    coverage_index = Path(f"{coverage_path}.csi")
    coverage_index_fingerprint = _require_file(
        coverage_index,
        label="Mosdepth regions BED.gz CSI index",
    )
    reference_fingerprint = _require_file(reference_fasta, label="reference FASTA")
    reference_index = Path(f"{reference_fasta}.fai")
    reference_index_fingerprint = _require_file(reference_index, label="reference FASTA index")
    metadata_input_fingerprint = (
        _require_file(metadata_path, label="Spectre metadata")
        if metadata_path is not None
        else None
    )
    blacklist_fingerprint = (
        _require_file(blacklist_path, label="Spectre blacklist")
        if blacklist_path is not None
        else None
    )

    command_runner = runner or SubprocessRunner()
    probe = command_runner.run([spectre, "version"], timeout_seconds=30)
    if probe.returncode != 0:
        raise ValueError(f"Spectre version probe failed with exit code {probe.returncode}")
    version = spectre_version(f"{probe.stdout}\n{probe.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"Spectre version {version!r} does not match policy lock {policy.expected_version!r}"
        )

    parameters: dict[str, str | int | bool | None] = {
        "mode": policy.mode,
        "threads": threads,
        "ploidy": policy.ploidy,
        "minimum_cnv_length_bp": policy.minimum_cnv_length_bp,
        "mosdepth_bin_size_bp": policy.mosdepth_bin_size_bp,
        "minimum_mapping_quality": policy.minimum_mapping_quality,
        "coordinate_contract": policy.coordinate_contract,
        "expected_version": policy.expected_version,
        "runtime_package": policy.runtime_package,
        "runtime_wheel_sha256": policy.runtime_wheel_sha256,
        "runtime_source_commit": policy.runtime_source_commit,
        "real_tool_qualified": policy.real_tool_qualified,
        "analytical_validation": policy.analytical_validation,
        "metadata_input": metadata_path is not None,
        "blacklist_input": blacklist_path is not None,
        "snfj_input": False,
        "snv_input": False,
        "population_mode": False,
        "cancer_mode": False,
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    promoted = False
    try:
        argv = build_spectre_depth_only_argv(
            spectre=spectre,
            coverage_path=coverage_path,
            sample_id=sample_id,
            output_dir=staged,
            reference_fasta=reference_fasta,
            policy=policy,
            threads=threads,
            metadata_path=metadata_path,
            blacklist_path=blacklist_path,
        )
        result = command_runner.run(argv, timeout_seconds=policy.timeout_seconds)
        if result.returncode != 0:
            diagnostic = result.stderr.strip()[-4000:]
            suffix = f": {diagnostic}" if diagnostic else ""
            raise ValueError(f"Spectre failed with exit code {result.returncode}{suffix}")

        vcf_candidates = [
            candidate
            for candidate in (staged / f"{sample_id}.vcf", staged / f"{sample_id}.vcf.gz")
            if candidate.is_file()
        ]
        if len(vcf_candidates) != 1:
            raise ValueError("Spectre must produce exactly one native .vcf or .vcf.gz output")
        artifacts = _native_artifacts(staged)
        report = normalize_spectre_vcf(
            vcf_candidates[0],
            sample_id=sample_id,
            genome_build=genome_build,
            policy=policy,
            tool=ToolRecord(name="Spectre", version=version, parameters=parameters),
            native_artifacts=artifacts,
            coverage_fingerprint=coverage_fingerprint,
            coverage_index_fingerprint=coverage_index_fingerprint,
            reference_fingerprint=reference_fingerprint,
            reference_index_fingerprint=reference_index_fingerprint,
            metadata_input_fingerprint=metadata_input_fingerprint,
            blacklist_fingerprint=blacklist_fingerprint,
        )
        os.replace(staged, output_dir)
        promoted = True
        return report
    finally:
        if not promoted:
            shutil.rmtree(staged, ignore_errors=True)
