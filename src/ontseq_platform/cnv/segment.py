"""A transparent read-depth baseline CNV caller.

Why a baseline caller belongs in this repository
------------------------------------------------

This is deliberately *not* an attempt to compete with Spectre, ichorCNA or QDNAseq. It
exists for three reasons that a production caller cannot serve:

1. **It closes the loop.** A benchmark harness that has never scored a real call set is
   untested infrastructure. This caller lets the whole path - simulate, call, mask,
   evaluate, aggregate - run in continuous integration with no external binary, no
   reference genome and no genomic data.
2. **It is a null model.** Any candidate method must beat a plain binned read-depth
   segmenter on the same data before its added complexity is justified. Without that
   floor, "our benchmark says method X works" has no reference point.
3. **It is fully inspectable.** Every step below is a few lines of arithmetic, so when
   the harness reports something surprising, the caller can be ruled in or out quickly.

Algorithm
---------

Bin counts are normalised to the autosomal median, converted to log2 ratios, and
segmented by recursive binary segmentation using a maximum standardised-difference
statistic. This is the deterministic core idea of circular binary segmentation without
the permutation-based significance test, which would make continuous-integration results
depend on a random seed. The noise scale is estimated from the median absolute deviation
of successive differences, so a large true alteration does not inflate the very variance
estimate used to decide whether it is real.

Limitations
-----------

No GC correction, no mappability correction, no allele-fraction information, no ploidy
search, and no subclonal deconvolution. It assumes bins are comparable, which is true for
simulated data and for uniform low-coverage whole-genome data, and false for on-target
capture or adaptive-sampling enrichment. It is a research control and is never reportable.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .core import StateSegment
from .intervals import canonical_contig
from .states import CopyNumberState, state_from_copy_number, tumor_copy_number_from_mixture

METHOD_NAME = "ontseq-baseline-readdepth"
METHOD_VERSION = "0.1.0"


@dataclass(frozen=True)
class DepthBin:
    """A counting window handed to the caller."""

    contig: str
    start: int
    end: int
    count: int


@dataclass(frozen=True)
class SegmentationParameters:
    """Caller parameters. All are engineering defaults, none are validated thresholds."""

    #: Standardised difference above which a split is accepted.
    split_threshold: float = 4.0
    #: Minimum number of bins on each side of a split.
    minimum_bins_per_segment: int = 3
    #: Adjacent segments closer than this in log2 space are merged.
    merge_log2_tolerance: float = 0.15
    #: Bins with fewer counts than this are treated as no-call rather than as depletion.
    minimum_bin_count: int = 1
    baseline_ploidy: float = 2.0
    #: When known, the mixture is inverted to recover tumor copy number.
    tumor_fraction: float | None = None

    def __post_init__(self) -> None:
        if self.split_threshold <= 0:
            raise ValueError("split threshold must be positive")
        if self.minimum_bins_per_segment < 1:
            raise ValueError("minimum bins per segment must be at least 1")
        if self.merge_log2_tolerance < 0:
            raise ValueError("merge tolerance must not be negative")
        if self.tumor_fraction is not None and not 0.0 < self.tumor_fraction <= 1.0:
            raise ValueError("tumor fraction must lie in (0, 1] when provided")


@dataclass(frozen=True)
class CallResult:
    """Baseline caller output: segments plus the regions it refused to call."""

    segments: list[StateSegment]
    no_call_regions: list[tuple[str, int, int]]
    median_count: float
    warnings: list[str]


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("median of an empty sequence is undefined")
    middle = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def _noise_scale(values: Sequence[float]) -> float:
    """Estimate the per-bin noise standard deviation robustly.

    Uses the median absolute successive difference rather than the overall standard
    deviation. A genuine whole-chromosome gain shifts a large block of bins and would
    dominate a naive variance estimate, making the caller less likely to detect the very
    event that inflated it.
    """
    if len(values) < 2:
        return 0.0
    differences = [abs(b - a) for a, b in zip(values, values[1:], strict=False)]
    scale = _median(differences) / (math.sqrt(2.0) * 0.6745)
    return scale if scale > 0 else 0.0


def _best_split(values: Sequence[float], scale: float, minimum: int) -> tuple[int, float] | None:
    """Return the split index maximising the standardised mean difference."""
    n = len(values)
    if n < 2 * minimum or scale <= 0:
        return None
    prefix = [0.0]
    for value in values:
        prefix.append(prefix[-1] + value)
    total = prefix[-1]
    best_index = -1
    best_statistic = 0.0
    for k in range(minimum, n - minimum + 1):
        left_mean = prefix[k] / k
        right_mean = (total - prefix[k]) / (n - k)
        standard_error = scale * math.sqrt(1.0 / k + 1.0 / (n - k))
        if standard_error <= 0:
            continue
        statistic = abs(left_mean - right_mean) / standard_error
        if statistic > best_statistic:
            best_statistic = statistic
            best_index = k
    if best_index < 0:
        return None
    return best_index, best_statistic


def _segment_indices(
    values: Sequence[float], scale: float, parameters: SegmentationParameters
) -> list[tuple[int, int]]:
    """Recursively split a value series, returning half-open index ranges."""
    boundaries: list[tuple[int, int]] = []
    stack: list[tuple[int, int]] = [(0, len(values))]
    while stack:
        start, end = stack.pop()
        window = values[start:end]
        split = _best_split(window, scale, parameters.minimum_bins_per_segment)
        if split is not None and split[1] >= parameters.split_threshold:
            stack.append((start, start + split[0]))
            stack.append((start + split[0], end))
        else:
            boundaries.append((start, end))
    return sorted(boundaries)


def call_segments(
    bins: Sequence[DepthBin], parameters: SegmentationParameters | None = None
) -> CallResult:
    """Run the baseline caller over a set of depth bins."""
    resolved = parameters or SegmentationParameters()
    warnings: list[str] = []
    usable = [item for item in bins if item.count >= resolved.minimum_bin_count]
    no_call_regions = [
        (canonical_contig(item.contig), item.start, item.end)
        for item in bins
        if item.count < resolved.minimum_bin_count
    ]
    if not usable:
        return CallResult(
            segments=[],
            no_call_regions=no_call_regions,
            median_count=0.0,
            warnings=["No bin met the minimum count; the caller produced no segments."],
        )

    autosomal = [item for item in usable if canonical_contig(item.contig) not in {"X", "Y"}]
    reference_bins = autosomal or usable
    if not autosomal:
        warnings.append(
            "No autosomal bins were available; normalisation used all bins including "
            "sex chromosomes, which biases the baseline."
        )
    median_count = _median([float(item.count) for item in reference_bins])
    if median_count <= 0:
        return CallResult(
            segments=[],
            no_call_regions=no_call_regions,
            median_count=0.0,
            warnings=["The median bin count was zero; log ratios are undefined."],
        )

    by_contig: dict[str, list[DepthBin]] = {}
    for item in usable:
        by_contig.setdefault(canonical_contig(item.contig), []).append(item)

    all_ratios = [
        math.log2(item.count / median_count) for items in by_contig.values() for item in items
    ]
    scale = _noise_scale(all_ratios)
    if scale <= 0:
        warnings.append(
            "The estimated noise scale was zero; segmentation was skipped and the whole "
            "genome reported as one level per contig."
        )

    segments: list[StateSegment] = []
    for contig in sorted(by_contig):
        items = sorted(by_contig[contig], key=lambda item: item.start)
        ratios = [math.log2(item.count / median_count) for item in items]
        ranges = _segment_indices(ratios, scale, resolved) if scale > 0 else [(0, len(ratios))]
        raw: list[tuple[int, int, float]] = []
        for start_index, end_index in ranges:
            mean_ratio = sum(ratios[start_index:end_index]) / (end_index - start_index)
            raw.append((start_index, end_index, mean_ratio))
        merged = _merge_levels(raw, resolved.merge_log2_tolerance, ratios)
        for start_index, end_index, mean_ratio in merged:
            observed = resolved.baseline_ploidy * (2.0**mean_ratio)
            copy_number = observed
            if resolved.tumor_fraction is not None:
                copy_number = tumor_copy_number_from_mixture(
                    observed,
                    tumor_fraction=resolved.tumor_fraction,
                    normal_copy_number=resolved.baseline_ploidy,
                )
            segments.append(
                StateSegment(
                    contig=contig,
                    start=items[start_index].start,
                    end=items[end_index - 1].end,
                    state=state_from_copy_number(
                        copy_number, baseline_ploidy=resolved.baseline_ploidy
                    ),
                    copy_number=copy_number,
                )
            )
    return CallResult(
        segments=segments,
        no_call_regions=no_call_regions,
        median_count=median_count,
        warnings=warnings,
    )


def _merge_levels(
    ranges: Sequence[tuple[int, int, float]], tolerance: float, ratios: Sequence[float]
) -> list[tuple[int, int, float]]:
    """Merge adjacent index ranges whose mean log2 ratios are within ``tolerance``."""
    if not ranges:
        return []
    merged: list[tuple[int, int, float]] = [ranges[0]]
    for start_index, end_index, mean_ratio in ranges[1:]:
        previous_start, previous_end, previous_mean = merged[-1]
        if previous_end == start_index and abs(previous_mean - mean_ratio) <= tolerance:
            combined = ratios[previous_start:end_index]
            merged[-1] = (
                previous_start,
                end_index,
                sum(combined) / len(combined),
            )
        else:
            merged.append((start_index, end_index, mean_ratio))
    return merged


def neutral_background_segments(
    segments: Sequence[StateSegment],
    contig_lengths: dict[str, int],
    *,
    baseline_ploidy: float = 2.0,
) -> list[StateSegment]:
    """Fill gaps between called segments with explicit neutral segments.

    The baseline caller emits a level for every bin it used, but bins it rejected leave
    gaps. Filling them explicitly makes the call set closed-world; leaving them as
    no-call regions instead keeps it open-world. The choice belongs to the caller of this
    function because it changes what specificity means.
    """
    by_contig: dict[str, list[StateSegment]] = {}
    for segment in segments:
        by_contig.setdefault(segment.contig, []).append(segment)
    filled: list[StateSegment] = []
    for contig, length in contig_lengths.items():
        key = canonical_contig(contig)
        items = sorted(by_contig.get(key, []), key=lambda item: item.start)
        cursor = 0
        for segment in items:
            if segment.start > cursor:
                filled.append(
                    StateSegment(
                        key, cursor, segment.start, CopyNumberState.NEUTRAL, baseline_ploidy
                    )
                )
            filled.append(segment)
            cursor = segment.end
        if cursor < length:
            filled.append(
                StateSegment(key, cursor, length, CopyNumberState.NEUTRAL, baseline_ploidy)
            )
    return filled
