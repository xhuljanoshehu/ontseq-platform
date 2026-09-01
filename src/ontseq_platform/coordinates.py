"""Single coordinate and contig normalization boundary for ONTSeq resources.

Every returned interval is 0-based, half-open. Callers also receive the original values in a
structured provenance record; converting a coordinate must never erase how the source encoded it.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import Field, model_validator

from .models import CoordinateSystem, StrictModel


class CoordinateSourceFormat(StrEnum):
    BED = "bed"
    GTF = "gtf"
    GFF3 = "gff3"
    VCF = "vcf"
    ONE_BASED_REGION = "one_based_region"
    INTERNAL = "internal"


class CoordinateConversion(StrictModel):
    source_format: CoordinateSourceFormat
    source_coordinate_system: CoordinateSystem
    target_coordinate_system: CoordinateSystem = CoordinateSystem.ZERO_BASED_HALF_OPEN
    original_contig: str = Field(min_length=1)
    normalized_contig: str = Field(min_length=1)
    original_start: int = Field(ge=0)
    original_end: int | None = Field(default=None, ge=0)
    normalized_start: int = Field(ge=0)
    normalized_end: int = Field(gt=0)
    operation: str = Field(min_length=1)

    @model_validator(mode="after")
    def interval_is_nonempty(self) -> CoordinateConversion:
        if self.normalized_end <= self.normalized_start:
            raise ValueError("normalized interval end must be greater than start")
        return self


class NormalizedInterval(StrictModel):
    contig: str = Field(min_length=1)
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    conversion: CoordinateConversion

    @model_validator(mode="after")
    def values_match_conversion(self) -> NormalizedInterval:
        if self.end <= self.start:
            raise ValueError("normalized interval end must be greater than start")
        expected = (
            self.conversion.normalized_contig,
            self.conversion.normalized_start,
            self.conversion.normalized_end,
        )
        if (self.contig, self.start, self.end) != expected:
            raise ValueError("normalized interval does not match its conversion provenance")
        return self


def normalize_contig(contig: str) -> str:
    """Normalize only canonical human aliases and preserve all other contig labels.

    ``1`` through ``22`` and ``X``/``Y`` become ``chr``-prefixed. The mitochondrial aliases
    ``MT``, ``M`` and ``chrM`` become ``chrM``. ALT/decoy labels are returned unchanged so the
    normalization layer never invents an alias it cannot prove.
    """

    if not isinstance(contig, str):
        raise TypeError("contig must be a string")
    original = contig.strip()
    if not original:
        raise ValueError("contig cannot be empty")
    upper = original.upper()
    without_prefix = upper[3:] if upper.startswith("CHR") else upper
    if without_prefix in {str(number) for number in range(1, 23)} | {"X", "Y"}:
        return f"chr{without_prefix}"
    if without_prefix in {"M", "MT"}:
        return "chrM"
    return original


def normalize_contig_name(contig: str) -> str:
    """Descriptive alias kept for callers that prefer an explicit ``*_name`` API."""

    return normalize_contig(contig)


def _require_int(value: int, name: str) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer")


def to_zero_based_half_open(
    start: int,
    end: int | None,
    coordinate_system: CoordinateSystem,
) -> tuple[int, int]:
    """Convert a supported source interval into ONTSeq's internal representation."""

    _require_int(start, "start")
    if end is not None:
        _require_int(end, "end")

    if coordinate_system == CoordinateSystem.ZERO_BASED_HALF_OPEN:
        if end is None:
            raise ValueError("0-based half-open intervals require an end")
        normalized_start, normalized_end = start, end
        if normalized_start < 0:
            raise ValueError("0-based interval start cannot be negative")
    elif coordinate_system == CoordinateSystem.ONE_BASED_INCLUSIVE:
        if end is None:
            raise ValueError("1-based inclusive intervals require an end")
        if start < 1:
            raise ValueError("1-based interval start must be at least 1")
        normalized_start, normalized_end = start - 1, end
    elif coordinate_system == CoordinateSystem.ONE_BASED_POSITION:
        if start < 1:
            raise ValueError("1-based position must be at least 1")
        if end is not None and end != start:
            raise ValueError("a 1-based position cannot declare a different end")
        normalized_start, normalized_end = start - 1, start
    else:  # pragma: no cover - the enum prevents this in typed callers
        raise ValueError(f"unsupported coordinate system: {coordinate_system}")

    if normalized_end <= normalized_start:
        raise ValueError("interval end must be greater than start")
    return normalized_start, normalized_end


def normalize_interval(
    contig: str,
    start: int,
    end: int | None,
    *,
    coordinate_system: CoordinateSystem,
    source_format: CoordinateSourceFormat,
    confirmed_standard_bed: bool = False,
) -> NormalizedInterval:
    """Normalize an interval and retain the exact source values as provenance.

    A BED is accepted unchanged only when the caller explicitly confirms it follows the BED
    standard. This prevents a filename extension from silently deciding coordinate semantics.
    """

    if source_format == CoordinateSourceFormat.BED:
        if coordinate_system != CoordinateSystem.ZERO_BASED_HALF_OPEN:
            raise ValueError("BED normalization requires zero_based_half_open coordinates")
        if not confirmed_standard_bed:
            raise ValueError("BED coordinates must be confirmed as standard 0-based half-open")
    normalized_contig = normalize_contig(contig)
    normalized_start, normalized_end = to_zero_based_half_open(start, end, coordinate_system)
    operation = (
        "unchanged"
        if coordinate_system == CoordinateSystem.ZERO_BASED_HALF_OPEN
        else "start_minus_one"
    )
    conversion = CoordinateConversion(
        source_format=source_format,
        source_coordinate_system=coordinate_system,
        original_contig=contig,
        normalized_contig=normalized_contig,
        original_start=start,
        original_end=end,
        normalized_start=normalized_start,
        normalized_end=normalized_end,
        operation=operation,
    )
    return NormalizedInterval(
        contig=normalized_contig,
        start=normalized_start,
        end=normalized_end,
        conversion=conversion,
    )


def bed_to_interval(
    contig: str, start: int, end: int, *, confirmed_standard_bed: bool
) -> NormalizedInterval:
    return normalize_interval(
        contig,
        start,
        end,
        coordinate_system=CoordinateSystem.ZERO_BASED_HALF_OPEN,
        source_format=CoordinateSourceFormat.BED,
        confirmed_standard_bed=confirmed_standard_bed,
    )


def gtf_to_interval(contig: str, start: int, end: int) -> NormalizedInterval:
    return normalize_interval(
        contig,
        start,
        end,
        coordinate_system=CoordinateSystem.ONE_BASED_INCLUSIVE,
        source_format=CoordinateSourceFormat.GTF,
    )


def gff3_to_interval(contig: str, start: int, end: int) -> NormalizedInterval:
    return normalize_interval(
        contig,
        start,
        end,
        coordinate_system=CoordinateSystem.ONE_BASED_INCLUSIVE,
        source_format=CoordinateSourceFormat.GFF3,
    )


def gtf_gff3_to_interval(
    contig: str,
    start: int,
    end: int,
    *,
    source_format: CoordinateSourceFormat,
) -> NormalizedInterval:
    if source_format not in {CoordinateSourceFormat.GTF, CoordinateSourceFormat.GFF3}:
        raise ValueError("source_format must be gtf or gff3")
    return normalize_interval(
        contig,
        start,
        end,
        coordinate_system=CoordinateSystem.ONE_BASED_INCLUSIVE,
        source_format=source_format,
    )


def vcf_to_interval(contig: str, position: int) -> NormalizedInterval:
    return normalize_interval(
        contig,
        position,
        None,
        coordinate_system=CoordinateSystem.ONE_BASED_POSITION,
        source_format=CoordinateSourceFormat.VCF,
    )


def one_based_inclusive_to_interval(contig: str, start: int, end: int) -> NormalizedInterval:
    return normalize_interval(
        contig,
        start,
        end,
        coordinate_system=CoordinateSystem.ONE_BASED_INCLUSIVE,
        source_format=CoordinateSourceFormat.ONE_BASED_REGION,
    )


__all__ = [
    "CoordinateConversion",
    "CoordinateSourceFormat",
    "NormalizedInterval",
    "bed_to_interval",
    "gff3_to_interval",
    "gtf_gff3_to_interval",
    "gtf_to_interval",
    "normalize_contig",
    "normalize_contig_name",
    "normalize_interval",
    "one_based_inclusive_to_interval",
    "to_zero_based_half_open",
    "vcf_to_interval",
]
