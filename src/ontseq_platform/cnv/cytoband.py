"""Deterministic GRCh38 CNV-to-cytoband annotation.

The engine retains every raw overlap, then applies a configurable fraction-of-band
threshold.  It merges only adjacent affected bands with the same copy-number direction
and never crosses the p/q centromere boundary.  Whole-chromosome events are represented
separately from focal/arm-level band groups.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from ..coordinates import normalize_contig

CnvDirection = Literal["gain", "loss"]


class CytobandAnnotationError(ValueError):
    """CNV segments or cytobands violate the 0-based half-open contract."""


@dataclass(frozen=True)
class Cytoband:
    chromosome: str
    start: int
    end: int
    name: str
    gie_stain: str | None = None

    def __post_init__(self) -> None:
        if self.start < 0 or self.end <= self.start:
            raise CytobandAnnotationError("cytoband must be a valid 0-based half-open interval")
        if not self.name:
            raise CytobandAnnotationError("cytoband name is required")

    @property
    def arm(self) -> Literal["p", "q"] | None:
        if self.name.startswith("p"):
            return "p"
        if self.name.startswith("q"):
            return "q"
        return None

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class CnvSegment:
    event_id: str
    chromosome: str
    start: int
    end: int
    direction: CnvDirection
    whole_chromosome: bool = False

    def __post_init__(self) -> None:
        if not self.event_id:
            raise CytobandAnnotationError("CNV event_id is required")
        if self.start < 0 or self.end <= self.start:
            raise CytobandAnnotationError("CNV segment must be 0-based half-open")
        if self.direction not in {"gain", "loss"}:
            raise CytobandAnnotationError("CNV direction must be 'gain' or 'loss'")


@dataclass(frozen=True)
class RawBandOverlap:
    event_id: str
    chromosome: str
    band: str
    band_start: int
    band_end: int
    direction: CnvDirection
    overlap_bp: int
    fraction_of_band: float
    affected: bool
    whole_chromosome: bool


@dataclass(frozen=True)
class AffectedBandGroup:
    chromosome: str
    arm: Literal["p", "q"] | None
    direction: CnvDirection
    start_band: str
    end_band: str
    start: int
    end: int
    bands: tuple[str, ...]
    source_event_ids: tuple[str, ...]


@dataclass(frozen=True)
class WholeChromosomeCall:
    chromosome: str
    direction: CnvDirection
    source_event_ids: tuple[str, ...]


@dataclass(frozen=True)
class CnvCytobandResult:
    threshold: float
    raw_overlaps: tuple[RawBandOverlap, ...]
    affected_groups: tuple[AffectedBandGroup, ...]
    whole_chromosome_calls: tuple[WholeChromosomeCall, ...]


@dataclass(frozen=True)
class _AffectedBand:
    chromosome: str
    index: int
    band: Cytoband
    direction: CnvDirection
    source_event_ids: tuple[str, ...]


def _is_whole_chromosome(
    segment: CnvSegment,
    chromosome_sizes: dict[str, int],
) -> bool:
    if segment.whole_chromosome:
        return True
    chromosome = normalize_contig(segment.chromosome)
    size = chromosome_sizes.get(chromosome)
    return size is not None and segment.start == 0 and segment.end >= size


def annotate_cnv_cytobands(
    segments: list[CnvSegment],
    cytobands: list[Cytoband],
    *,
    affected_fraction: float,
    chromosome_sizes: dict[str, int] | None = None,
) -> CnvCytobandResult:
    """Annotate CNV segments while preserving raw and summarized representations."""
    if not 0 < affected_fraction <= 1:
        raise CytobandAnnotationError("affected_fraction must be greater than 0 and at most 1")
    normalized_sizes = {
        normalize_contig(chromosome): size for chromosome, size in (chromosome_sizes or {}).items()
    }
    for chromosome, size in normalized_sizes.items():
        if size <= 0:
            raise CytobandAnnotationError(f"chromosome size for {chromosome} must be positive")
    ordered_bands: dict[str, list[Cytoband]] = {}
    for band in cytobands:
        chromosome = normalize_contig(band.chromosome)
        ordered_bands.setdefault(chromosome, []).append(
            Cytoband(chromosome, band.start, band.end, band.name, band.gie_stain)
        )
    for chromosome, bands in ordered_bands.items():
        bands.sort(key=lambda band: (band.start, band.end, band.name))
        previous_end = -1
        for band in bands:
            if band.start < previous_end:
                raise CytobandAnnotationError(
                    f"overlapping cytoband definitions on {chromosome} are ambiguous"
                )
            previous_end = band.end

    raw: list[RawBandOverlap] = []
    affected_sources: dict[tuple[str, int, CnvDirection], set[str]] = {}
    whole_sources: dict[tuple[str, CnvDirection], set[str]] = {}
    for segment in segments:
        chromosome = normalize_contig(segment.chromosome)
        whole = _is_whole_chromosome(segment, normalized_sizes)
        if whole:
            whole_sources.setdefault((chromosome, segment.direction), set()).add(segment.event_id)
        for index, band in enumerate(ordered_bands.get(chromosome, [])):
            overlap_bp = max(0, min(segment.end, band.end) - max(segment.start, band.start))
            if overlap_bp == 0:
                continue
            fraction = overlap_bp / band.length
            affected = fraction >= affected_fraction
            raw.append(
                RawBandOverlap(
                    event_id=segment.event_id,
                    chromosome=chromosome,
                    band=band.name,
                    band_start=band.start,
                    band_end=band.end,
                    direction=segment.direction,
                    overlap_bp=overlap_bp,
                    fraction_of_band=fraction,
                    affected=affected,
                    whole_chromosome=whole,
                )
            )
            if affected and not whole:
                affected_sources.setdefault((chromosome, index, segment.direction), set()).add(
                    segment.event_id
                )

    affected_bands = [
        _AffectedBand(
            chromosome=chromosome,
            index=index,
            band=ordered_bands[chromosome][index],
            direction=direction,
            source_event_ids=tuple(sorted(event_ids)),
        )
        for (chromosome, index, direction), event_ids in affected_sources.items()
    ]
    affected_bands.sort(
        key=lambda item: (item.chromosome, item.direction, item.index, item.band.name)
    )

    groups: list[AffectedBandGroup] = []
    current: list[_AffectedBand] = []

    def flush() -> None:
        if not current:
            return
        first, last = current[0], current[-1]
        groups.append(
            AffectedBandGroup(
                chromosome=first.chromosome,
                arm=first.band.arm,
                direction=first.direction,
                start_band=first.band.name,
                end_band=last.band.name,
                start=first.band.start,
                end=last.band.end,
                bands=tuple(item.band.name for item in current),
                source_event_ids=tuple(
                    sorted({event_id for item in current for event_id in item.source_event_ids})
                ),
            )
        )
        current.clear()

    for item in affected_bands:
        if not current:
            current.append(item)
            continue
        previous = current[-1]
        adjacent = item.index == previous.index + 1
        compatible = (
            item.chromosome == previous.chromosome
            and item.direction == previous.direction
            and item.band.arm == previous.band.arm
            and item.band.arm is not None
        )
        if adjacent and compatible:
            current.append(item)
        else:
            flush()
            current.append(item)
    flush()

    whole_calls = tuple(
        WholeChromosomeCall(chromosome, direction, tuple(sorted(event_ids)))
        for (chromosome, direction), event_ids in sorted(whole_sources.items())
    )
    return CnvCytobandResult(
        threshold=affected_fraction,
        raw_overlaps=tuple(
            sorted(
                raw,
                key=lambda overlap: (
                    overlap.chromosome,
                    overlap.band_start,
                    overlap.direction,
                    overlap.event_id,
                ),
            )
        ),
        affected_groups=tuple(groups),
        whole_chromosome_calls=whole_calls,
    )
