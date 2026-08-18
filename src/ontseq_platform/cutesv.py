from __future__ import annotations

import gzip
import re
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from typing import Literal

from pydantic import Field, ValidationError, model_validator

from .models import (
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
from .reference import sha256_file
from .sv_concordance import SVCallerObservation

_BND_MATE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
_SV_TYPE_MAP = {
    "DEL": EventType.DELETION,
    "DUP": EventType.DUPLICATION,
    "INV": EventType.INVERSION,
    "INS": EventType.INSERTION,
    "BND": EventType.TRANSLOCATION,
    "TRA": EventType.TRANSLOCATION,
}


class CuteSVPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    expected_version: str = "2.1.4"
    min_support: int = Field(default=1, ge=1)
    min_sv_length: int = Field(default=1, ge=1)
    pass_only: bool = True
    allowed_sv_types: set[EventType] = Field(
        default_factory=lambda: {
            EventType.DELETION,
            EventType.DUPLICATION,
            EventType.INVERSION,
            EventType.INSERTION,
            EventType.TRANSLOCATION,
        }
    )
    report_read_ids: Literal[False] = False
    status: Literal["technical_defaults_only"] = "technical_defaults_only"
    clinically_validated: Literal[False] = False
    note: str = (
        "cuteSV policy values are software-engineering defaults only and require assay-specific "
        "benchmarking before any biological or clinical interpretation."
    )


class CuteSVCallReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: CuteSVPolicy
    events: list[GenomicEvent] = Field(default_factory=list)
    raw_record_count: int = Field(ge=0)
    accepted_record_count: int = Field(ge=0)
    rejected_record_count: int = Field(ge=0)
    rejection_counts: dict[str, int] = Field(default_factory=dict)
    tool: ToolRecord
    vcf_fingerprint: FileFingerprint
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def accounting_is_consistent(self) -> CuteSVCallReport:
        if self.raw_record_count != self.accepted_record_count + self.rejected_record_count:
            raise ValueError("cuteSV record accounting is inconsistent")
        if self.accepted_record_count != len(self.events):
            raise ValueError("cuteSV accepted count does not match normalized events")
        if self.rejected_record_count != sum(self.rejection_counts.values()):
            raise ValueError("cuteSV rejection accounting is inconsistent")
        expected = ModuleRunStatus.COMPLETED if self.events else ModuleRunStatus.NO_CALL
        if self.status != expected:
            raise ValueError("cuteSV status does not match normalized event count")
        return self


class _RejectedRecord(ValueError):
    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def _open_vcf(path: Path) -> Iterator[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield from handle


def _parse_info(raw: str) -> dict[str, str | bool]:
    parsed: dict[str, str | bool] = {}
    if raw == ".":
        return parsed
    for item in raw.split(";"):
        key, separator, value = item.partition("=")
        if not key:
            raise _RejectedRecord("malformed_info")
        parsed[key] = value if separator else True
    return parsed


def _first_value(value: str | bool | None) -> str | None:
    if not isinstance(value, str):
        return None
    first = value.split(",", maxsplit=1)[0]
    return None if first in {"", "."} else first


def _required_int(value: str | bool | None, *, reason: str) -> int:
    first = _first_value(value)
    if first is None:
        raise _RejectedRecord(reason)
    try:
        return int(first)
    except ValueError as exc:
        raise _RejectedRecord(reason) from exc


def _optional_float(value: str | bool | None, *, reason: str) -> float | None:
    first = _first_value(value)
    if first is None:
        return None
    try:
        return float(first)
    except ValueError as exc:
        raise _RejectedRecord(reason) from exc


def _quality(raw: str) -> float | None:
    if raw in {"", "."}:
        return None
    try:
        result = float(raw)
    except ValueError as exc:
        raise _RejectedRecord("malformed_quality") from exc
    if result < 0:
        raise _RejectedRecord("malformed_quality")
    return result


def _simple_locus(
    chromosome: str,
    position: int,
    event_type: EventType,
    info: dict[str, str | bool],
) -> Locus:
    start = position - 1
    if event_type == EventType.INSERTION:
        end = start + 1
    else:
        raw_end = _first_value(info.get("END"))
        if raw_end is not None:
            try:
                end = int(raw_end)
            except ValueError as exc:
                raise _RejectedRecord("malformed_end") from exc
        else:
            length = abs(_required_int(info.get("SVLEN"), reason="missing_end_and_svlen"))
            end = start + max(length, 1)
    try:
        return Locus(chromosome=chromosome, start=start, end=end)
    except ValidationError as exc:
        raise _RejectedRecord("invalid_primary_locus") from exc


def _breakend_loci(chromosome: str, position: int, alternate: str) -> tuple[Locus, Locus]:
    mate = _BND_MATE.search(alternate)
    if mate is None:
        raise _RejectedRecord("missing_breakend_mate")
    secondary_chromosome = mate.group(1)
    secondary_position = int(mate.group(2))
    try:
        primary = Locus(chromosome=chromosome, start=position - 1, end=position)
        secondary = Locus(
            chromosome=secondary_chromosome,
            start=secondary_position - 1,
            end=secondary_position,
        )
    except ValidationError as exc:
        raise _RejectedRecord("invalid_breakend_locus") from exc
    return primary, secondary


def _sv_length(
    event_type: EventType,
    primary: Locus,
    info: dict[str, str | bool],
) -> int | None:
    if event_type == EventType.TRANSLOCATION:
        return None
    raw = _first_value(info.get("SVLEN"))
    if raw is not None:
        try:
            length = abs(int(raw))
        except ValueError as exc:
            raise _RejectedRecord("malformed_svlen") from exc
    elif event_type == EventType.INSERTION:
        raise _RejectedRecord("missing_svlen")
    else:
        length = primary.end - primary.start
    if length < 1:
        raise _RejectedRecord("malformed_svlen")
    return length


def _normalize_record(
    fields: list[str],
    *,
    record_number: int,
    policy: CuteSVPolicy,
    caller_version: str,
) -> GenomicEvent:
    if len(fields) < 8:
        raise _RejectedRecord("malformed_record")
    chromosome, raw_position, _record_id, _ref, alternate, raw_quality, raw_filter, raw_info = (
        fields[:8]
    )
    try:
        position = int(raw_position)
    except ValueError as exc:
        raise _RejectedRecord("malformed_position") from exc
    if position < 1:
        raise _RejectedRecord("malformed_position")

    filters = raw_filter.split(";") if raw_filter not in {"", "."} else []
    if policy.pass_only and filters != ["PASS"]:
        raise _RejectedRecord("filter_not_pass")

    info = _parse_info(raw_info)
    raw_sv_type = _first_value(info.get("SVTYPE"))
    if raw_sv_type is None:
        raise _RejectedRecord("missing_svtype")
    event_type = _SV_TYPE_MAP.get(raw_sv_type.upper().split(":", maxsplit=1)[0])
    if event_type is None:
        raise _RejectedRecord("unsupported_svtype")
    if event_type not in policy.allowed_sv_types:
        raise _RejectedRecord("svtype_not_allowed")

    support = _required_int(info.get("RE"), reason="missing_or_invalid_support")
    if support < policy.min_support:
        raise _RejectedRecord("support_below_policy")

    quality = _quality(raw_quality)
    vaf = _optional_float(info.get("AF"), reason="malformed_af")
    if vaf is not None and not 0 <= vaf <= 1:
        raise _RejectedRecord("invalid_af")

    if event_type == EventType.TRANSLOCATION:
        primary, secondary = _breakend_loci(chromosome, position, alternate)
    else:
        primary = _simple_locus(chromosome, position, event_type, info)
        secondary = None
    length_bp = _sv_length(event_type, primary, info)
    if length_bp is not None and length_bp < policy.min_sv_length:
        raise _RejectedRecord("sv_length_below_policy")

    notes = [
        "Normalized from a cuteSV VCF into 0-based, half-open coordinates.",
        "Candidate only; assay-specific analytical validation and expert review are pending.",
    ]
    if event_type == EventType.TRANSLOCATION:
        notes.append("Breakend evidence is not equivalent to an annotated gene fusion.")

    return GenomicEvent(
        event_id=f"CUTESV-{record_number:06d}",
        event_type=event_type,
        primary=primary,
        secondary=secondary,
        length_bp=length_bp,
        evidence=[
            Evidence(
                caller="cuteSV",
                caller_version=caller_version,
                support_reads=support,
                variant_allele_fraction=vaf,
                quality=quality,
                filters=filters,
            )
        ],
        confidence="unclassified",
        reportable=False,
        notes=notes,
    )


def normalize_cutesv_vcf(
    path: Path,
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    policy: CuteSVPolicy,
    tool: ToolRecord,
) -> CuteSVCallReport:
    if not path.is_file():
        raise ValueError("cuteSV VCF is missing or unreadable")
    if tool.version != policy.expected_version:
        raise ValueError(
            f"cuteSV version {tool.version!r} does not match policy lock "
            f"{policy.expected_version!r}"
        )

    raw_count = 0
    events: list[GenomicEvent] = []
    rejections: Counter[str] = Counter()
    saw_fileformat = False
    saw_columns = False
    sample_column_count: int | None = None

    for raw_line in _open_vcf(path):
        line = raw_line.rstrip("\r\n")
        if line.startswith("##fileformat=VCF"):
            saw_fileformat = True
            continue
        if line.startswith("#CHROM"):
            saw_columns = True
            sample_column_count = max(0, len(line.split("\t")) - 9)
            continue
        if not line or line.startswith("#"):
            continue
        raw_count += 1
        try:
            event = _normalize_record(
                line.split("\t"),
                record_number=raw_count,
                policy=policy,
                caller_version=tool.version,
            )
        except _RejectedRecord as exc:
            rejections[exc.reason] += 1
        else:
            events.append(event)

    if not saw_fileformat or not saw_columns:
        raise ValueError("cuteSV output is not a complete VCF document")
    if sample_column_count not in {0, 1}:
        raise ValueError("Single-sample analysis requires at most one VCF sample column")

    warnings = [policy.note]
    if rejections:
        warnings.append("One or more cuteSV records were rejected by the explicit policy.")
    if not events:
        warnings.append(
            "NO_CALL means no cuteSV record passed this technical policy; it is not a "
            "biological or clinical negative result."
        )

    return CuteSVCallReport(
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
        vcf_fingerprint=FileFingerprint(size_bytes=path.stat().st_size, sha256=sha256_file(path)),
        warnings=warnings,
        limitations=[
            "cuteSV record IDs, read names, inserted sequences and local source paths are not "
            "propagated into normalized output.",
            "BND/TRA records are normalized as DNA translocation evidence, never as confirmed "
            "gene fusions.",
            "Technical policy values are not analytically or clinically validated thresholds.",
        ],
    )


def cutesv_observations(report: CuteSVCallReport) -> list[SVCallerObservation]:
    observations: list[SVCallerObservation] = []
    for event in report.events:
        if len(event.evidence) != 1:
            raise ValueError("cuteSV normalized event must contain exactly one evidence record")
        observations.append(
            SVCallerObservation(
                observation_id=event.event_id,
                caller="cuteSV",
                caller_version=report.tool.version,
                source_event_id=event.event_id,
                event_type=event.event_type,
                primary=event.primary,
                secondary=event.secondary,
                evidence=event.evidence[0],
            )
        )
    return observations
