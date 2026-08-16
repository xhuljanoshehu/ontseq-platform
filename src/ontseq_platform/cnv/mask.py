"""Construction of the evaluable genome.

A CNV metric is only meaningful relative to the part of the genome the assay could
actually interrogate. This module builds that region explicitly and keeps a per-reason
accounting of everything it removed, so a reader can always answer "how much of the
genome did this number apply to, and what was thrown away?".

The distinction this enables is the one the project cares about most:

- **negative** - the region was observable and no alteration was found;
- **no-call** - the region was not observable, so no statement is possible;
- **failure** - the module did not run or errored, handled by ``ModuleRunStatus``.

Folding the second into the first is the single most common way a CNV pipeline
overstates its own sensitivity.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum

from .intervals import (
    Interval,
    IntervalSet,
    canonical_contig,
    contig_lengths_to_set,
    intersect_set,
    normalize_set,
    subtract_set,
    total_length,
    union_set,
)


class ExclusionReason(StrEnum):
    """Why a region was removed from the evaluable genome.

    The vocabulary is closed so that exclusion accounting can be compared across runs
    and methods. A free-text reason would make two runs incomparable.
    """

    ASSEMBLY_GAP = "assembly_gap"
    CENTROMERE = "centromere"
    LOW_MAPPABILITY = "low_mappability"
    BLACKLIST = "blacklist"
    BELOW_COVERAGE_FLOOR = "below_coverage_floor"
    OUTSIDE_ANALYSIS_SCOPE = "outside_analysis_scope"
    CALLER_NO_CALL = "caller_no_call"
    TRUTH_NOT_INFORMATIVE = "truth_not_informative"
    CONTIG_NOT_IN_REFERENCE = "contig_not_in_reference"


@dataclass(frozen=True)
class ExclusionTrack:
    """One named set of regions removed for a single reason."""

    reason: ExclusionReason
    intervals: IntervalSet
    source: str
    #: Optional checksum of the file the intervals came from, for provenance.
    source_sha256: str | None = None


@dataclass(frozen=True)
class ObservabilityMask:
    """The evaluable genome plus a full account of what was removed and why."""

    evaluable: IntervalSet
    reference_bases: int
    evaluable_bases: int
    excluded_bases_by_reason: dict[str, int]
    #: Bases excluded by more than one track are attributed to the first track that
    #: removed them, so the per-reason values sum to the total excluded bases.
    analysis_scope_bases: int
    tracks: list[ExclusionTrack] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def excluded_bases(self) -> int:
        return self.reference_bases - self.evaluable_bases

    @property
    def evaluable_fraction(self) -> float | None:
        if self.reference_bases <= 0:
            return None
        return self.evaluable_bases / self.reference_bases


def build_mask(
    *,
    contig_lengths: Mapping[str, int],
    analysis_scope: Mapping[str, Sequence[Interval]] | None = None,
    tracks: Sequence[ExclusionTrack] = (),
) -> ObservabilityMask:
    """Build the evaluable genome from a reference, an analysis scope and exclusions.

    ``analysis_scope`` restricts evaluation before any exclusion is applied. For a
    low-coverage whole-genome lane it is ``None``, meaning the whole reference. For a
    panel or adaptive-sampling lane it is the analysis ROI, and everything outside is
    attributed to :attr:`ExclusionReason.OUTSIDE_ANALYSIS_SCOPE` rather than silently
    dropped: a caller cannot be blamed for missing an event it was never able to see.

    Overlapping exclusion tracks are attributed to the first track that removes a base,
    in the order given, so that the per-reason counts sum exactly to the total removed.
    """
    genome = contig_lengths_to_set(contig_lengths)
    reference_bases = total_length(genome)
    warnings: list[str] = []

    if analysis_scope is None:
        current = genome
        scope_bases = reference_bases
        excluded: dict[str, int] = {}
    else:
        scope = normalize_set(analysis_scope)
        unknown = sorted(set(scope) - set(genome))
        if unknown:
            warnings.append(
                "Analysis scope references contigs absent from the reference lock and "
                f"they were dropped: {', '.join(unknown)}"
            )
        clipped = intersect_set(scope, genome)
        current = clipped
        scope_bases = total_length(clipped)
        excluded = {ExclusionReason.OUTSIDE_ANALYSIS_SCOPE.value: reference_bases - scope_bases}

    applied: list[ExclusionTrack] = []
    for track in tracks:
        before = total_length(current)
        current = subtract_set(current, track.intervals)
        removed = before - total_length(current)
        if removed:
            excluded[track.reason.value] = excluded.get(track.reason.value, 0) + removed
        applied.append(track)

    evaluable_bases = total_length(current)
    if evaluable_bases == 0:
        warnings.append(
            "The evaluable genome is empty after masking; every downstream metric will "
            "be undefined rather than zero."
        )
    return ObservabilityMask(
        evaluable=current,
        reference_bases=reference_bases,
        evaluable_bases=evaluable_bases,
        excluded_bases_by_reason=excluded,
        analysis_scope_bases=scope_bases,
        tracks=applied,
        warnings=warnings,
    )


def coverage_floor_track(
    depths: Sequence[tuple[str, int, int, float]],
    *,
    minimum_depth: float,
    source: str,
    source_sha256: str | None = None,
) -> ExclusionTrack:
    """Build an exclusion track from per-interval depth below a floor.

    This is the connection point to the adaptive-sampling target-coverage QC adapter:
    intervals whose observed depth falls below the floor become no-call rather than
    being scored as copy-number neutral.

    The floor is an engineering parameter for benchmarking. It is not a validated
    adequacy threshold and must not be reused as a clinical no-call limit.
    """
    if minimum_depth < 0:
        raise ValueError("minimum depth must not be negative")
    below: dict[str, list[Interval]] = {}
    for contig, start, end, depth in depths:
        if end <= start:
            raise ValueError(f"depth interval {contig}:{start}-{end} is empty or inverted")
        if depth < minimum_depth:
            below.setdefault(canonical_contig(contig), []).append((start, end))
    return ExclusionTrack(
        reason=ExclusionReason.BELOW_COVERAGE_FLOOR,
        intervals=normalize_set(below),
        source=source,
        source_sha256=source_sha256,
    )


def no_call_track(
    intervals: Mapping[str, Sequence[Interval]],
    *,
    source: str,
    reason: ExclusionReason = ExclusionReason.CALLER_NO_CALL,
    source_sha256: str | None = None,
) -> ExclusionTrack:
    """Build an exclusion track from regions a caller or truth source declared unusable."""
    return ExclusionTrack(
        reason=reason,
        intervals=normalize_set(intervals),
        source=source,
        source_sha256=source_sha256,
    )


def merge_masks(*masks: ObservabilityMask) -> IntervalSet:
    """Intersect several masks, for example one per method in a fair comparison.

    Comparing two callers on different evaluable genomes is not a comparison. When
    several methods are benchmarked together, every method must be scored on the
    intersection of what all of them could observe, otherwise a method that declines to
    call in hard regions is rewarded for its own blind spots.
    """
    if not masks:
        return {}
    result = masks[0].evaluable
    for mask in masks[1:]:
        result = intersect_set(result, mask.evaluable)
    return result


def scope_from_intervals(
    intervals: Mapping[str, Sequence[Interval]], *, flank_bp: int = 0
) -> IntervalSet:
    """Expand an interval set by a symmetric flank, clamped at zero.

    Used to derive an evaluation scope from a target BED when the intended comparison
    includes a defined margin around each target. The flank must be recorded in
    provenance because it changes every denominator.
    """
    if flank_bp < 0:
        raise ValueError("flank must not be negative")
    if flank_bp == 0:
        return normalize_set(intervals)
    expanded: dict[str, list[Interval]] = {}
    for contig, spans in intervals.items():
        expanded[canonical_contig(contig)] = [
            (max(0, start - flank_bp), end + flank_bp) for start, end in spans
        ]
    return normalize_set(expanded)


def union_of_tracks(tracks: Sequence[ExclusionTrack]) -> IntervalSet:
    """Return the union of every track's intervals, for reporting."""
    return union_set(*[track.intervals for track in tracks]) if tracks else {}
