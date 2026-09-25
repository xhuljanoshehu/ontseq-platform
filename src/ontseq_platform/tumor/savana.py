from __future__ import annotations

import csv
import math
import os
import re
import shutil
import tempfile
from collections import Counter
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
_SHA256 = r"^[0-9a-f]{64}$"
_SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")


class SavanaBasePolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    expected_version: str = Field(pattern=r"^\d+\.\d+(?:\.\d+)?$")
    mode: Literal["paired", "tumor_only"]
    real_tool_qualified: bool = False
    analytical_validation: Literal["not_validated"] = "not_validated"
    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    reference_sha256: str = Field(pattern=_SHA256)
    minimum_sv_length: int = Field(ge=1)
    minimum_mapping_quality: int = Field(ge=0)
    minimum_support: int = Field(ge=1)
    minimum_allele_fraction: float = Field(ge=0, le=1)
    cn_bin_size_kbp: int = Field(ge=1)
    timeout_seconds: int = Field(ge=60)
    note: str = Field(min_length=12)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def pinned_runtime_version(self) -> SavanaBasePolicy:
        if self.expected_version != "1.3.8":
            raise ValueError("SAVANA expected_version must match the pinned 1.3.8 contract")
        return self


class SavanaPairedPolicy(SavanaBasePolicy):
    mode: Literal["paired"] = "paired"


class SavanaTumorOnlyPolicy(SavanaBasePolicy):
    mode: Literal["tumor_only"] = "tumor_only"


class SavanaFit(StrictModel):
    purity: float = Field(ge=0, le=1)
    ploidy: float = Field(gt=0)
    distance: float = Field(ge=0)
    rank: int = Field(ge=1)


class SavanaCopyNumberSegment(StrictModel):
    chromosome: str = Field(pattern=r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    segment_id: str = Field(min_length=1)
    bin_count: int | None = Field(default=None, ge=0)
    sum_of_bin_lengths: int | None = Field(default=None, ge=0)
    weight: float | None = Field(default=None, ge=0)
    copy_number: float = Field(ge=0)
    minor_allele_copy_number: float | None = Field(default=None, ge=0)
    mean_baf: float | None = Field(default=None, ge=0, le=1)
    heterozygous_snp_count: int | None = Field(default=None, ge=0)

    @model_validator(mode="after")
    def valid_interval(self) -> SavanaCopyNumberSegment:
        if self.end <= self.start:
            raise ValueError("SAVANA copy-number segment end must be greater than start")
        return self


class SavanaNativeArtifact(StrictModel):
    relative_path: str = Field(min_length=1)
    role: str = Field(min_length=1)
    media_type: str = Field(min_length=1)
    fingerprint: FileFingerprint
    sensitive_output: bool = False
    exportable: bool = True
    research_only: Literal[True] = True

    @field_validator("relative_path")
    @classmethod
    def safe_relative_path(cls, value: str) -> str:
        if "\\" in value:
            raise ValueError("SAVANA artifact paths must use POSIX separators")
        path = PurePosixPath(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("SAVANA artifact path must remain inside output directory")
        return value

    @model_validator(mode="after")
    def sensitive_artifacts_are_not_exportable(self) -> SavanaNativeArtifact:
        if self.sensitive_output and self.exportable:
            raise ValueError("Sensitive SAVANA native artifacts cannot be exportable")
        return self


class SavanaReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str = Field(min_length=1)
    normal_sample_id: str | None = None
    mode: Literal["paired", "tumor_only"]
    genome_build: GenomeBuild
    status: ModuleRunStatus
    sv_status: ModuleRunStatus
    cna_status: ModuleRunStatus
    policy: SavanaPairedPolicy | SavanaTumorOnlyPolicy
    events: list[GenomicEvent]
    raw_sv_record_count: int = Field(ge=0)
    accepted_sv_record_count: int = Field(ge=0)
    rejected_sv_record_count: int = Field(ge=0)
    rejection_counts: dict[str, int] = Field(default_factory=dict)
    selected_fit: SavanaFit | None = None
    ranked_fits: list[SavanaFit] = Field(default_factory=list)
    copy_number_segments: list[SavanaCopyNumberSegment] = Field(default_factory=list)
    tool: ToolRecord
    tumor_bam_fingerprint: FileFingerprint
    tumor_bam_index_fingerprint: FileFingerprint
    normal_bam_fingerprint: FileFingerprint | None = None
    normal_bam_index_fingerprint: FileFingerprint | None = None
    snp_vcf_fingerprint: FileFingerprint
    reference_fingerprint: FileFingerprint
    native_artifacts: list[SavanaNativeArtifact]
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_report(self) -> SavanaReport:
        if self.mode != self.policy.mode:
            raise ValueError("SAVANA report mode must match policy mode")
        if self.mode == "paired":
            if self.normal_sample_id is None:
                raise ValueError("Paired SAVANA report requires normal sample ID")
            if self.normal_bam_fingerprint is None or self.normal_bam_index_fingerprint is None:
                raise ValueError("Paired SAVANA report requires normal BAM fingerprints")
        elif self.normal_sample_id is not None:
            raise ValueError("Tumor-only SAVANA report cannot carry a normal sample ID")

        if self.raw_sv_record_count != (
            self.accepted_sv_record_count + self.rejected_sv_record_count
        ):
            raise ValueError("SAVANA SV record counts are inconsistent")
        if self.accepted_sv_record_count != len(self.events):
            raise ValueError("SAVANA accepted SV count must match normalized events")
        if self.sv_status == ModuleRunStatus.COMPLETED and not self.events:
            raise ValueError("Completed SAVANA SV status requires at least one event")
        if self.sv_status == ModuleRunStatus.NO_CALL and self.events:
            raise ValueError("SAVANA SV NO_CALL cannot carry normalized events")
        if self.cna_status == ModuleRunStatus.COMPLETED and (
            self.selected_fit is None or not self.copy_number_segments
        ):
            raise ValueError("Completed SAVANA CNA status requires fit and segments")
        if self.status == ModuleRunStatus.NO_CALL and (
            self.sv_status == ModuleRunStatus.COMPLETED
            or self.cna_status == ModuleRunStatus.COMPLETED
        ):
            raise ValueError("Overall SAVANA NO_CALL conflicts with completed submodule")
        return self


class _RejectedRecord(ValueError):
    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def savana_version(text: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    return text.splitlines()[0].strip()[:80] if text.strip() else "unknown"


def build_savana_paired_argv(
    *,
    savana: str,
    tumor_bam: Path,
    normal_bam: Path,
    snp_vcf: Path,
    reference_fasta: Path,
    sample_id: str,
    output_dir: Path,
    policy: SavanaPairedPolicy,
    threads: int = 8,
) -> tuple[str, ...]:
    if threads < 1:
        raise ValueError("SAVANA threads must be at least 1")
    return (
        savana,
        "--tumour",
        str(tumor_bam),
        "--normal",
        str(normal_bam),
        "--outdir",
        str(output_dir),
        "--ref",
        str(reference_fasta),
        "--snp_vcf",
        str(snp_vcf),
        "--sample",
        sample_id,
        "--length",
        str(policy.minimum_sv_length),
        "--mapq",
        str(policy.minimum_mapping_quality),
        "--min_support",
        str(policy.minimum_support),
        "--min_af",
        str(policy.minimum_allele_fraction),
        "--cn_binsize",
        str(policy.cn_bin_size_kbp),
        "--threads",
        str(threads),
        "--cna_threads",
        str(threads),
    )


def build_savana_tumor_only_argv(
    *,
    savana: str,
    tumor_bam: Path,
    snp_vcf: Path,
    reference_fasta: Path,
    sample_id: str,
    output_dir: Path,
    policy: SavanaTumorOnlyPolicy,
    threads: int = 8,
) -> tuple[str, ...]:
    if threads < 1:
        raise ValueError("SAVANA threads must be at least 1")
    return (
        savana,
        "to",
        "--tumour",
        str(tumor_bam),
        "--outdir",
        str(output_dir),
        "--ref",
        str(reference_fasta),
        "--snp_vcf",
        str(snp_vcf),
        "--sample",
        sample_id,
        "--length",
        str(policy.minimum_sv_length),
        "--mapq",
        str(policy.minimum_mapping_quality),
        "--min_support",
        str(policy.minimum_support),
        "--min_af",
        str(policy.minimum_allele_fraction),
        "--cn_binsize",
        str(policy.cn_bin_size_kbp),
        "--threads",
        str(threads),
        "--cna_threads",
        str(threads),
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
        raise ValueError(f"SAVANA requires registered {label}")
    return artifact


def _bam_index_path(path: Path) -> Path:
    candidates = (Path(f"{path}.bai"), path.with_suffix(".bai"), Path(f"{path}.csi"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError(f"SAVANA requires an index for BAM: {path.name}")


def _validate_registered_file(
    path: Path,
    artifact: TumorAuxiliaryInputArtifact,
    *,
    label: str,
) -> FileFingerprint:
    if not path.is_file():
        raise ValueError(f"SAVANA requires {label}: {path}")
    fingerprint = _fingerprint(path)
    if fingerprint.sha256 != artifact.sha256:
        raise ValueError(f"SAVANA {label} SHA-256 does not match registered input")
    return fingerprint


def _validate_reference(
    reference_fasta: Path,
    inputs: TumorInputBundle,
    policy: SavanaBasePolicy,
) -> FileFingerprint:
    if inputs.genome_build != policy.genome_build:
        raise ValueError("SAVANA input genome build does not match policy")
    if inputs.reference_id != policy.reference_id:
        raise ValueError("SAVANA input reference ID does not match policy")
    if inputs.reference_sha256 != policy.reference_sha256:
        raise ValueError("SAVANA input reference SHA-256 does not match policy")
    if not reference_fasta.is_file():
        raise ValueError("SAVANA requires the locked reference FASTA")
    if not Path(f"{reference_fasta}.fai").is_file():
        raise ValueError("SAVANA requires the reference FASTA index")
    fingerprint = _fingerprint(reference_fasta)
    if fingerprint.sha256 != policy.reference_sha256:
        raise ValueError("SAVANA reference FASTA SHA-256 does not match policy lock")
    return fingerprint


def _validate_common_inputs(
    *,
    tumor_bam: Path,
    snp_vcf: Path,
    sample_id: str,
    reference_fasta: Path,
    inputs: TumorInputBundle,
    policy: SavanaBasePolicy,
) -> tuple[FileFingerprint, FileFingerprint, FileFingerprint, FileFingerprint]:
    if not _SAFE_SAMPLE_ID.fullmatch(sample_id) or sample_id in {".", ".."}:
        raise ValueError("SAVANA sample ID is unsafe")
    if inputs.analysis_sample_id != sample_id:
        raise ValueError("SAVANA tumor sample does not match registered analysis sample")

    tumor = _registered_artifact(inputs, CallerInputRole.TUMOR_BAM, label="tumor BAM")
    snps = _registered_artifact(inputs, CallerInputRole.SNP_VCF, label="SNP VCF")
    if tumor.source_sample_id != sample_id:
        raise ValueError("SAVANA tumor BAM source sample does not match requested sample")

    tumor_fp = _validate_registered_file(tumor_bam, tumor, label="tumor BAM")
    tumor_index_fp = _fingerprint(_bam_index_path(tumor_bam))
    snp_fp = _validate_registered_file(snp_vcf, snps, label="SNP VCF")
    reference_fp = _validate_reference(reference_fasta, inputs, policy)
    return tumor_fp, tumor_index_fp, snp_fp, reference_fp


def _validate_paired_inputs(
    *,
    tumor_bam: Path,
    normal_bam: Path,
    snp_vcf: Path,
    sample_id: str,
    normal_sample_id: str,
    reference_fasta: Path,
    inputs: TumorInputBundle,
    policy: SavanaPairedPolicy,
) -> tuple[
    FileFingerprint,
    FileFingerprint,
    FileFingerprint,
    FileFingerprint,
    FileFingerprint,
    FileFingerprint,
]:
    tumor_fp, tumor_index_fp, snp_fp, reference_fp = _validate_common_inputs(
        tumor_bam=tumor_bam,
        snp_vcf=snp_vcf,
        sample_id=sample_id,
        reference_fasta=reference_fasta,
        inputs=inputs,
        policy=policy,
    )
    if not _SAFE_SAMPLE_ID.fullmatch(normal_sample_id) or normal_sample_id in {".", ".."}:
        raise ValueError("SAVANA normal sample ID is unsafe")
    normal = _registered_artifact(
        inputs,
        CallerInputRole.MATCHED_NORMAL_BAM,
        label="matched normal BAM",
    )
    if normal.source_sample_id != normal_sample_id:
        raise ValueError("SAVANA matched normal source sample does not match requested normal")
    normal_fp = _validate_registered_file(normal_bam, normal, label="matched normal BAM")
    normal_index_fp = _fingerprint(_bam_index_path(normal_bam))
    return (
        tumor_fp,
        tumor_index_fp,
        normal_fp,
        normal_index_fp,
        snp_fp,
        reference_fp,
    )


def _parse_info(raw: str) -> dict[str, str | bool]:
    parsed: dict[str, str | bool] = {}
    if raw in {"", "."}:
        return parsed
    for token in raw.split(";"):
        if not token:
            continue
        key, separator, value = token.partition("=")
        if not key:
            raise _RejectedRecord("malformed_info")
        parsed[key] = value if separator else True
    return parsed


def _text(info: dict[str, str | bool], key: str) -> str | None:
    value = info.get(key)
    return value if isinstance(value, str) and value not in {"", "."} else None


def _finite_float(raw: str, *, reason: str) -> float:
    try:
        value = float(raw)
    except ValueError as exc:
        raise _RejectedRecord(reason) from exc
    if not math.isfinite(value):
        raise _RejectedRecord(reason)
    return value


def _integer(info: dict[str, str | bool], key: str, *, reason: str) -> int | None:
    raw = _text(info, key)
    if raw is None:
        return None
    try:
        value = int(raw)
    except ValueError as exc:
        raise _RejectedRecord(reason) from exc
    if value < 0:
        raise _RejectedRecord(reason)
    return value


def _event_loci(
    *,
    chromosome: str,
    position: int,
    alternate: str,
    svtype: str,
    info: dict[str, str | bool],
) -> tuple[EventType, Locus, Locus | None, int | None]:
    start = position - 1
    try:
        if svtype == "INS":
            length_raw = _text(info, "SVLEN")
            length = int(length_raw) if length_raw is not None else 1
            if length < 1:
                raise _RejectedRecord("malformed_svlen")
            return (
                EventType.INSERTION,
                Locus(chromosome=chromosome, start=start, end=start + 1),
                None,
                length,
            )
        if svtype == "BND":
            try:
                resolved = resolve_breakend(
                    alternate,
                    declared_chromosome=info.get("CHR2"),
                    declared_position=info.get("END"),
                )
            except BreakendParseError as exc:
                raise _RejectedRecord(exc.reason) from exc
            return (
                EventType.TRANSLOCATION,
                Locus(chromosome=chromosome, start=start, end=start + 1),
                Locus(
                    chromosome=resolved.mate_chromosome,
                    start=resolved.mate_position_0based,
                    end=resolved.mate_position_0based + 1,
                ),
                None,
            )
    except (ValidationError, ValueError) as exc:
        if isinstance(exc, _RejectedRecord):
            raise
        raise _RejectedRecord("invalid_locus") from exc
    raise _RejectedRecord("unsupported_svtype")


def _normalize_sv_record(
    fields: list[str],
    *,
    record_number: int,
    caller_version: str,
    policy: SavanaBasePolicy,
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
    if (_text(info, "CLASS") or "").upper() not in {"", "SOMATIC"}:
        raise _RejectedRecord("not_somatic")
    svtype = (_text(info, "SVTYPE") or "").upper()
    event_type, primary, secondary, length = _event_loci(
        chromosome=chromosome,
        position=position,
        alternate=alternate,
        svtype=svtype,
        info=info,
    )
    if length is not None and length < policy.minimum_sv_length:
        raise _RejectedRecord("sv_length_below_policy")

    support = _integer(
        info,
        "TUMOUR_READ_SUPPORT",
        reason="malformed_tumor_read_support",
    )
    if support is not None and support < policy.minimum_support:
        raise _RejectedRecord("support_below_policy")

    allele_fraction = None
    raw_af = _text(info, "TUMOUR_AF")
    if raw_af is not None:
        allele_fraction = _finite_float(raw_af, reason="malformed_tumor_af")
        if not 0 <= allele_fraction <= 1:
            raise _RejectedRecord("malformed_tumor_af")
        if allele_fraction < policy.minimum_allele_fraction:
            raise _RejectedRecord("allele_fraction_below_policy")

    quality = None
    if raw_quality not in {"", "."}:
        quality = _finite_float(raw_quality, reason="malformed_quality")
        if quality < 0:
            raise _RejectedRecord("malformed_quality")

    filters = [] if raw_filter in {"", "."} else raw_filter.split(";")
    notes = [
        "Normalized from SAVANA somatic VCF into 0-based half-open coordinates.",
        "Candidate only; analytical validation and expert review remain required.",
    ]
    if native_id not in {"", "."}:
        notes.append(f"SAVANA native candidate ID: {native_id}")
    bp_notation = _text(info, "BP_NOTATION")
    if bp_notation is not None:
        notes.append(f"SAVANA BP_NOTATION={bp_notation}")

    return GenomicEvent(
        event_id=f"SAVANA-{record_number:06d}",
        event_type=event_type,
        primary=primary,
        secondary=secondary,
        length_bp=length,
        evidence=[
            Evidence(
                caller="SAVANA",
                caller_version=caller_version,
                support_reads=support,
                variant_allele_fraction=allele_fraction,
                quality=quality,
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
    policy: SavanaBasePolicy,
) -> tuple[list[GenomicEvent], int, Counter[str]]:
    saw_fileformat = False
    saw_columns = False
    raw_count = 0
    events: list[GenomicEvent] = []
    rejections: Counter[str] = Counter()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("##fileformat=VCF"):
                saw_fileformat = True
                continue
            if line.startswith("#CHROM"):
                if saw_columns:
                    raise ValueError("SAVANA VCF contains duplicate column headers")
                columns = line.split("\t")
                expected = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                if columns[:8] != expected:
                    raise ValueError("SAVANA VCF has an unexpected header")
                saw_columns = True
                continue
            if not line or line.startswith("#"):
                continue
            raw_count += 1
            try:
                events.append(
                    _normalize_sv_record(
                        line.split("\t"),
                        record_number=raw_count,
                        caller_version=caller_version,
                        policy=policy,
                    )
                )
            except _RejectedRecord as exc:
                rejections[exc.reason] += 1
    if not saw_fileformat or not saw_columns:
        raise ValueError("SAVANA output is not a complete VCF document")
    return events, raw_count, rejections


def _column(row: dict[str, str], *names: str) -> str | None:
    for name in names:
        value = row.get(name)
        if value not in {None, ""}:
            return value
    return None


def _required_float(row: dict[str, str], *names: str) -> float:
    raw = _column(row, *names)
    if raw is None:
        raise ValueError(f"SAVANA table is missing required numeric field: {names[0]}")
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(f"SAVANA table contains invalid numeric field: {names[0]}") from exc
    if not math.isfinite(value):
        raise ValueError(f"SAVANA table contains non-finite numeric field: {names[0]}")
    return value


def _optional_float(row: dict[str, str], *names: str) -> float | None:
    raw = _column(row, *names)
    if raw is None:
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(f"SAVANA table contains invalid numeric field: {names[0]}") from exc
    if not math.isfinite(value):
        raise ValueError(f"SAVANA table contains non-finite numeric field: {names[0]}")
    return value


def _optional_int(row: dict[str, str], *names: str) -> int | None:
    raw = _column(row, *names)
    if raw is None:
        return None
    try:
        return int(raw)
    except ValueError as exc:
        raise ValueError(f"SAVANA table contains invalid integer field: {names[0]}") from exc


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _parse_fit(row: dict[str, str]) -> SavanaFit:
    purity = _required_float(row, "purity", "cellularity")
    ploidy = _required_float(row, "ploidy")
    distance = _required_float(row, "distance", "distance_score")
    rank_raw = _column(row, "rank")
    if rank_raw is None:
        rank = 1
    else:
        try:
            rank = int(rank_raw)
        except ValueError as exc:
            raise ValueError("SAVANA fit contains invalid rank") from exc
    return SavanaFit(purity=purity, ploidy=ploidy, distance=distance, rank=rank)


def _matching_outputs(
    output_dir: Path,
    *,
    exact: str,
    pattern: str,
) -> list[Path]:
    exact_path = output_dir / exact
    candidates = [exact_path] if exact_path.is_file() else sorted(output_dir.glob(pattern))
    return list(dict.fromkeys(candidates))


def _find_output(
    output_dir: Path,
    *,
    exact: str,
    pattern: str,
    label: str,
) -> Path:
    unique = _matching_outputs(output_dir, exact=exact, pattern=pattern)
    if len(unique) != 1:
        raise ValueError(f"SAVANA must produce exactly one {label}")
    return unique[0]


def _has_complete_cna_output(output_dir: Path, sample_id: str) -> bool:
    specifications = (
        (
            f"{sample_id}_ranked_solutions.tsv",
            f"{sample_id}_*_ranked_solutions.tsv",
            "ranked purity/ploidy solutions table",
        ),
        (
            f"{sample_id}_fitted_purity_ploidy.tsv",
            f"{sample_id}_*_fitted_purity_ploidy.tsv",
            "selected purity/ploidy fit table",
        ),
        (
            f"{sample_id}_segmented_absolute_copy_number.tsv",
            f"{sample_id}_*_segmented_absolute_copy_number.tsv",
            "segmented absolute copy-number table",
        ),
    )
    present: list[bool] = []
    for exact, pattern, label in specifications:
        matches = _matching_outputs(output_dir, exact=exact, pattern=pattern)
        if len(matches) > 1:
            raise ValueError(f"SAVANA must produce at most one {label}")
        present.append(bool(matches))
    if not any(present):
        return False
    if not all(present):
        raise ValueError("SAVANA CNA output is incomplete")
    return True


def _fit_outputs(output_dir: Path, sample_id: str) -> tuple[SavanaFit, list[SavanaFit]]:
    ranked_path = _find_output(
        output_dir,
        exact=f"{sample_id}_ranked_solutions.tsv",
        pattern=f"{sample_id}_*_ranked_solutions.tsv",
        label="ranked purity/ploidy solutions table",
    )
    fitted_path = _find_output(
        output_dir,
        exact=f"{sample_id}_fitted_purity_ploidy.tsv",
        pattern=f"{sample_id}_*_fitted_purity_ploidy.tsv",
        label="selected purity/ploidy fit table",
    )
    ranked_rows = _read_tsv(ranked_path)
    fitted_rows = _read_tsv(fitted_path)
    if not ranked_rows:
        raise ValueError("SAVANA ranked purity/ploidy solutions table is empty")
    if len(fitted_rows) != 1:
        raise ValueError("SAVANA selected purity/ploidy table must contain exactly one row")
    ranked = sorted((_parse_fit(row) for row in ranked_rows), key=lambda item: item.rank)
    ranks = [item.rank for item in ranked]
    if len(set(ranks)) != len(ranks):
        raise ValueError("SAVANA ranked purity/ploidy solutions table has duplicate ranks")
    if ranked[0].rank != 1:
        raise ValueError("SAVANA ranked purity/ploidy solutions table is missing rank 1")
    selected = _parse_fit(fitted_rows[0])
    top = ranked[0]
    if (selected.purity, selected.ploidy, selected.distance, selected.rank) != (
        top.purity,
        top.ploidy,
        top.distance,
        top.rank,
    ):
        raise ValueError(
            "SAVANA selected purity/ploidy fit does not match the rank-1 ranked solution"
        )
    return selected, ranked


def _copy_number_segments(output_dir: Path, sample_id: str) -> list[SavanaCopyNumberSegment]:
    path = _find_output(
        output_dir,
        exact=f"{sample_id}_segmented_absolute_copy_number.tsv",
        pattern=f"{sample_id}_*_segmented_absolute_copy_number.tsv",
        label="segmented absolute copy-number table",
    )
    result: list[SavanaCopyNumberSegment] = []
    for row in _read_tsv(path):
        chromosome = _column(row, "chromosome", "chrom")
        segment_id = _column(row, "segment_id", "segment")
        start_raw = _column(row, "start")
        end_raw = _column(row, "end")
        if chromosome is None or segment_id is None or start_raw is None or end_raw is None:
            raise ValueError("SAVANA copy-number segment is missing identity fields")
        try:
            start = int(start_raw)
            end = int(end_raw)
        except ValueError as exc:
            raise ValueError("SAVANA copy-number segment coordinates are invalid") from exc
        result.append(
            SavanaCopyNumberSegment(
                chromosome=chromosome,
                start=start,
                end=end,
                segment_id=segment_id,
                bin_count=_optional_int(row, "bin_count"),
                sum_of_bin_lengths=_optional_int(row, "sum_of_bin_lengths"),
                weight=_optional_float(row, "weight"),
                copy_number=_required_float(row, "copyNumber", "copy_number"),
                minor_allele_copy_number=_optional_float(
                    row,
                    "minorAlleleCopyNumber",
                    "minor_allele_copy_number",
                ),
                mean_baf=_optional_float(row, "meanBAF", "mean_baf"),
                heterozygous_snp_count=_optional_int(row, "no_hetSNPs", "het_snp_count"),
            )
        )
    if not result:
        raise ValueError("SAVANA segmented absolute copy-number table is empty")
    return result


def _artifact_role(relative_path: str) -> str:
    name = PurePosixPath(relative_path).name
    if name.endswith(".classified.somatic.vcf"):
        return "somatic_vcf"
    if name.endswith(".classified.somatic.bedpe"):
        return "somatic_bedpe"
    if name.endswith("_sv_breakpoints_read_support.tsv"):
        return "sv_read_support"
    if name.endswith("_ranked_solutions.tsv"):
        return "ranked_purity_ploidy_solutions"
    if name.endswith("_fitted_purity_ploidy.tsv"):
        return "selected_purity_ploidy_fit"
    if name.endswith("_segmented_absolute_copy_number.tsv"):
        return "absolute_copy_number_segments"
    if name.endswith("_raw_read_counts.tsv"):
        return "raw_copy_number_read_counts"
    if name.endswith(".inserted_sequences.fa"):
        return "inserted_sequences"
    return "caller_native"


def _media_type(relative_path: str) -> str:
    if relative_path.endswith(".vcf"):
        return "text/vcf"
    if relative_path.endswith(".bedpe"):
        return "text/tab-separated-values"
    if relative_path.endswith(".tsv"):
        return "text/tab-separated-values"
    if relative_path.endswith((".fa", ".fasta")):
        return "text/fasta"
    if relative_path.endswith((".txt", ".log")):
        return "text/plain"
    return "application/octet-stream"


def _is_sensitive(relative_path: str) -> bool:
    name = PurePosixPath(relative_path).name
    return name.endswith("_sv_breakpoints_read_support.tsv") or name.endswith(
        ".inserted_sequences.fa"
    )


def _native_artifacts(output_dir: Path) -> list[SavanaNativeArtifact]:
    result: list[SavanaNativeArtifact] = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.is_symlink():
            raise ValueError("SAVANA native output must not contain symbolic links")
        relative = path.relative_to(output_dir).as_posix()
        sensitive = _is_sensitive(relative)
        result.append(
            SavanaNativeArtifact(
                relative_path=relative,
                role=_artifact_role(relative),
                media_type=_media_type(relative),
                fingerprint=_fingerprint(path),
                sensitive_output=sensitive,
                exportable=not sensitive,
            )
        )
    if not result:
        raise ValueError("SAVANA returned success without native artifacts")
    return result


def _somatic_vcf(output_dir: Path, sample_id: str) -> Path:
    path = output_dir / f"{sample_id}.classified.somatic.vcf"
    if not path.is_file():
        raise ValueError("SAVANA must produce the classified somatic VCF")
    return path


def _probe_version(
    runner: CommandRunner,
    *,
    savana: str,
    expected_version: str,
) -> str:
    probe = runner.run([savana, "--version"], timeout_seconds=30)
    if probe.returncode != 0:
        raise ValueError(f"SAVANA version probe failed with exit code {probe.returncode}")
    version = savana_version(f"{probe.stdout}\n{probe.stderr}")
    if version != expected_version:
        raise ValueError(
            f"SAVANA version {version!r} does not match policy lock {expected_version!r}"
        )
    return version


def _parameters(policy: SavanaBasePolicy) -> dict[str, str | int | float | bool]:
    return {
        "mode": policy.mode,
        "minimum_sv_length": policy.minimum_sv_length,
        "minimum_mapping_quality": policy.minimum_mapping_quality,
        "minimum_support": policy.minimum_support,
        "minimum_allele_fraction": policy.minimum_allele_fraction,
        "cn_bin_size_kbp": policy.cn_bin_size_kbp,
        "real_tool_qualified": policy.real_tool_qualified,
        "analytical_validation": policy.analytical_validation,
    }


def _statuses(
    *,
    events: list[GenomicEvent],
    selected_fit: SavanaFit | None,
    segments: list[SavanaCopyNumberSegment],
) -> tuple[ModuleRunStatus, ModuleRunStatus, ModuleRunStatus]:
    sv_status = ModuleRunStatus.COMPLETED if events else ModuleRunStatus.NO_CALL
    cna_status = (
        ModuleRunStatus.COMPLETED
        if selected_fit is not None and segments
        else ModuleRunStatus.NO_CALL
    )
    overall = (
        ModuleRunStatus.COMPLETED
        if sv_status == ModuleRunStatus.COMPLETED or cna_status == ModuleRunStatus.COMPLETED
        else ModuleRunStatus.NO_CALL
    )
    return overall, sv_status, cna_status


def _build_report(
    *,
    staged: Path,
    sample_id: str,
    normal_sample_id: str | None,
    mode: Literal["paired", "tumor_only"],
    policy: SavanaPairedPolicy | SavanaTumorOnlyPolicy,
    version: str,
    tumor_fp: FileFingerprint,
    tumor_index_fp: FileFingerprint,
    normal_fp: FileFingerprint | None,
    normal_index_fp: FileFingerprint | None,
    snp_fp: FileFingerprint,
    reference_fp: FileFingerprint,
) -> SavanaReport:
    somatic_vcf = _somatic_vcf(staged, sample_id)
    events, raw_count, rejections = _normalize_somatic_vcf(
        somatic_vcf,
        caller_version=version,
        policy=policy,
    )

    if _has_complete_cna_output(staged, sample_id):
        selected_fit, ranked_fits = _fit_outputs(staged, sample_id)
        segments = _copy_number_segments(staged, sample_id)
    else:
        selected_fit = None
        ranked_fits = []
        segments = []

    overall, sv_status, cna_status = _statuses(
        events=events,
        selected_fit=selected_fit,
        segments=segments,
    )
    warnings = [policy.note]
    if rejections:
        warnings.append("One or more SAVANA SV records were rejected during normalization.")
    if sv_status == ModuleRunStatus.NO_CALL:
        warnings.append(
            "SAVANA SV NO_CALL means no somatic SV record was retained; "
            "it is not a biological negative result."
        )
    if cna_status == ModuleRunStatus.NO_CALL:
        warnings.append(
            "SAVANA CNA NO_CALL means no complete copy-number fit was retained; "
            "it is not a biological negative result."
        )

    return SavanaReport(
        sample_id=sample_id,
        normal_sample_id=normal_sample_id,
        mode=mode,
        genome_build=policy.genome_build,
        status=overall,
        sv_status=sv_status,
        cna_status=cna_status,
        policy=policy,
        events=events,
        raw_sv_record_count=raw_count,
        accepted_sv_record_count=len(events),
        rejected_sv_record_count=sum(rejections.values()),
        rejection_counts=dict(sorted(rejections.items())),
        selected_fit=selected_fit,
        ranked_fits=ranked_fits,
        copy_number_segments=segments,
        tool=ToolRecord(name="SAVANA", version=version, parameters=_parameters(policy)),
        tumor_bam_fingerprint=tumor_fp,
        tumor_bam_index_fingerprint=tumor_index_fp,
        normal_bam_fingerprint=normal_fp,
        normal_bam_index_fingerprint=normal_index_fp,
        snp_vcf_fingerprint=snp_fp,
        reference_fingerprint=reference_fp,
        native_artifacts=_native_artifacts(staged),
        warnings=warnings,
        limitations=[
            "SAVANA output is research evidence and does not authorize clinical release.",
            "Paired and tumor-only modes are separate analytical contracts without fallback.",
            (
                "Caller-native read identifiers and inserted sequences are retained only "
                "as sensitive artifacts."
            ),
            "Caller agreement must not be interpreted as biological truth.",
        ],
    )


def run_savana_paired(
    *,
    tumor_bam: Path,
    normal_bam: Path,
    snp_vcf: Path,
    reference_fasta: Path,
    sample_id: str,
    normal_sample_id: str,
    output_dir: Path,
    inputs: TumorInputBundle,
    policy: SavanaPairedPolicy,
    runner: CommandRunner | None = None,
    savana: str = "savana",
    threads: int = 8,
) -> SavanaReport:
    if output_dir.exists():
        raise ValueError("Refusing to overwrite an existing SAVANA output directory")
    (
        tumor_fp,
        tumor_index_fp,
        normal_fp,
        normal_index_fp,
        snp_fp,
        reference_fp,
    ) = _validate_paired_inputs(
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        snp_vcf=snp_vcf,
        sample_id=sample_id,
        normal_sample_id=normal_sample_id,
        reference_fasta=reference_fasta,
        inputs=inputs,
        policy=policy,
    )
    command_runner = runner or SubprocessRunner()
    version = _probe_version(
        command_runner,
        savana=savana,
        expected_version=policy.expected_version,
    )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    promoted = False
    try:
        argv = build_savana_paired_argv(
            savana=savana,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            snp_vcf=snp_vcf,
            reference_fasta=reference_fasta,
            sample_id=sample_id,
            output_dir=staged,
            policy=policy,
            threads=threads,
        )
        result = command_runner.run(argv, timeout_seconds=policy.timeout_seconds)
        if result.returncode != 0:
            diagnostic = result.stderr.strip()[-4000:]
            suffix = f": {diagnostic}" if diagnostic else ""
            raise ValueError(f"SAVANA failed with exit code {result.returncode}{suffix}")

        report = _build_report(
            staged=staged,
            sample_id=sample_id,
            normal_sample_id=normal_sample_id,
            mode="paired",
            policy=policy,
            version=version,
            tumor_fp=tumor_fp,
            tumor_index_fp=tumor_index_fp,
            normal_fp=normal_fp,
            normal_index_fp=normal_index_fp,
            snp_fp=snp_fp,
            reference_fp=reference_fp,
        )
        os.replace(staged, output_dir)
        promoted = True
        return report
    finally:
        if not promoted:
            shutil.rmtree(staged, ignore_errors=True)


def run_savana_tumor_only(
    *,
    tumor_bam: Path,
    snp_vcf: Path,
    reference_fasta: Path,
    sample_id: str,
    output_dir: Path,
    inputs: TumorInputBundle,
    policy: SavanaTumorOnlyPolicy,
    runner: CommandRunner | None = None,
    savana: str = "savana",
    threads: int = 8,
) -> SavanaReport:
    if output_dir.exists():
        raise ValueError("Refusing to overwrite an existing SAVANA output directory")
    tumor_fp, tumor_index_fp, snp_fp, reference_fp = _validate_common_inputs(
        tumor_bam=tumor_bam,
        snp_vcf=snp_vcf,
        sample_id=sample_id,
        reference_fasta=reference_fasta,
        inputs=inputs,
        policy=policy,
    )
    if inputs.artifact_for(CallerInputRole.MATCHED_NORMAL_BAM) is not None:
        raise ValueError("SAVANA tumor-only mode must not consume a matched normal BAM")

    command_runner = runner or SubprocessRunner()
    version = _probe_version(
        command_runner,
        savana=savana,
        expected_version=policy.expected_version,
    )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    promoted = False
    try:
        argv = build_savana_tumor_only_argv(
            savana=savana,
            tumor_bam=tumor_bam,
            snp_vcf=snp_vcf,
            reference_fasta=reference_fasta,
            sample_id=sample_id,
            output_dir=staged,
            policy=policy,
            threads=threads,
        )
        result = command_runner.run(argv, timeout_seconds=policy.timeout_seconds)
        if result.returncode != 0:
            diagnostic = result.stderr.strip()[-4000:]
            suffix = f": {diagnostic}" if diagnostic else ""
            raise ValueError(f"SAVANA failed with exit code {result.returncode}{suffix}")

        report = _build_report(
            staged=staged,
            sample_id=sample_id,
            normal_sample_id=None,
            mode="tumor_only",
            policy=policy,
            version=version,
            tumor_fp=tumor_fp,
            tumor_index_fp=tumor_index_fp,
            normal_fp=None,
            normal_index_fp=None,
            snp_fp=snp_fp,
            reference_fp=reference_fp,
        )
        os.replace(staged, output_dir)
        promoted = True
        return report
    finally:
        if not promoted:
            shutil.rmtree(staged, ignore_errors=True)
