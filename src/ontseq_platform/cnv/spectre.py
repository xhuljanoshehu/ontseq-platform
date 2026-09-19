from __future__ import annotations

import gzip
import math
import re
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from typing import Literal

from pydantic import Field, ValidationError, model_validator

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


class SpectrePolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    expected_version: str = Field(pattern=r"^\d+\.\d+(?:\.\d+)?$")
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
    if len(fields) < 8:
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

    filters = raw_filter.split(";") if raw_filter not in {"", "."} else []
    if policy.pass_only and filters != ["PASS"]:
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
    saw_columns = False
    raw_count = 0
    events: list[GenomicEvent] = []
    rejections: Counter[str] = Counter()
    for raw_line in _open_vcf(path):
        line = raw_line.rstrip("\r\n")
        if line.startswith("##fileformat=VCF"):
            saw_fileformat = True
            continue
        if line.startswith("#CHROM"):
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
