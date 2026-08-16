"""Segmentation-independent CNV comparison.

Why this exists
---------------

The repository's first benchmark layer matches CNV calls as *events* with a one-to-one
assignment and a reciprocal-overlap threshold. That contract is appropriate for
structural variants, where a breakend pair is the unit of truth, but it misreports copy
number for three reasons:

1. **Segmentation is not part of the biological claim.** If truth carries one 90 Mb
   deletion and a caller emits three adjacent deleted segments covering the same span,
   a one-to-one matcher scores one true positive and two false positives. The caller was
   right. Segment boundaries are an artifact of the segmentation algorithm and the bin
   size, not a claim about the genome.
2. **Event counting ignores event size.** A whole-chromosome gain and a 200 kb duplication
   contribute equally to precision and recall, even though they differ by three orders of
   magnitude in genomic content and in clinical weight.
3. **Unmatched is not the same as wrong.** Counting every unmatched call as a false
   positive assumes every base was assessable. Centromeres, assembly gaps, regions below
   the coverage floor and regions outside an adaptive-sampling target design are not
   assessable, and scoring them silently converts a known blind spot into a fabricated
   error rate.

This module therefore scores copy number the way copy number is actually defined: as a
*per-base state assignment across the genome*. Base-level agreement is invariant to how
either side chose to segment. Event-level detection is still reported, but derived from
base-level concordance with a many-to-many rule, so fragmentation cannot be penalised.

Everything here is dependency-free and operates on plain dataclasses. The pydantic
contract layer lives in :mod:`ontseq_platform.cnv.models` and converts into these types.
"""

from __future__ import annotations

from bisect import bisect_right
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum

from .intervals import Interval, IntervalSet, canonical_contig, intersect, normalize_set
from .states import (
    ConcordanceMode,
    CopyNumberState,
    concordant,
    is_altered,
)
from .stats import ProportionEstimate, mean_absolute_error, root_mean_square_error, wilson_interval


class TruthOutcome(StrEnum):
    """Per-truth-event result.

    ``NOT_ASSESSABLE`` is the reason this enum exists. An event that falls inside a
    region the assay cannot observe is neither a hit nor a miss, and folding it into
    either one corrupts sensitivity in opposite directions depending on which choice is
    made.
    """

    DETECTED = "DETECTED"
    MISSED = "MISSED"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"


class QueryOutcome(StrEnum):
    """Per-called-event result."""

    CONFIRMED = "CONFIRMED"
    UNCONFIRMED = "UNCONFIRMED"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"


class SizeClass(StrEnum):
    """Length band used to stratify detection.

    Resolution in read-depth CNV analysis is governed by bin size and coverage, so
    sensitivity is a strong function of event length. A single pooled recall number
    hides that dependence completely.
    """

    UNDER_100KB = "lt_100kb"
    KB100_TO_1MB = "100kb_1mb"
    MB1_TO_5MB = "1mb_5mb"
    MB5_TO_20MB = "5mb_20mb"
    OVER_20MB = "ge_20mb"


_SIZE_BOUNDS: tuple[tuple[int, SizeClass], ...] = (
    (100_000, SizeClass.UNDER_100KB),
    (1_000_000, SizeClass.KB100_TO_1MB),
    (5_000_000, SizeClass.MB1_TO_5MB),
    (20_000_000, SizeClass.MB5_TO_20MB),
)


def size_class(length_bp: int) -> SizeClass:
    """Return the length band of an event."""
    for bound, label in _SIZE_BOUNDS:
        if length_bp < bound:
            return label
    return SizeClass.OVER_20MB


@dataclass(frozen=True)
class StateSegment:
    """A contiguous span carrying one copy-number state.

    Boundary uncertainties travel with the segment rather than in a side table keyed by
    event identity. Events are formed by merging adjacent same-state segments, so a side
    table keyed by the merged span would silently lose the per-segment uncertainty of
    exactly the truth sets that need it most.
    """

    contig: str
    start: int
    end: int
    state: CopyNumberState
    copy_number: float | None = None
    start_uncertainty_bp: int | None = None
    end_uncertainty_bp: int | None = None

    def __post_init__(self) -> None:
        if self.end <= self.start:
            raise ValueError("segment end must be greater than start")
        if self.start < 0:
            raise ValueError("segment start must not be negative")
        for value in (self.start_uncertainty_bp, self.end_uncertainty_bp):
            if value is not None and value < 0:
                raise ValueError("boundary uncertainty must not be negative")

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class EvaluationOptions:
    """Tunable, pre-registerable comparison parameters.

    Every value here changes reported performance, so the full option set is echoed into
    the evaluation report. Thresholds must be locked before comparative results are
    inspected; see ``docs/CNV_BENCHMARKING.md``.
    """

    concordance_mode: ConcordanceMode = ConcordanceMode.DIRECTIONAL
    #: Fraction of an event's *evaluable* bases that must be concordant to count as detected.
    detection_overlap_fraction: float = 0.5
    #: Fraction of an event's total bases that must be evaluable for it to be scored at all.
    minimum_assessable_fraction: float = 0.5
    #: Absolute copy-number distance still counted as a correct quantitative estimate.
    copy_number_tolerance: float = 0.5
    #: Truth boundary uncertainty above which breakpoint accuracy is not reported.
    maximum_truth_boundary_uncertainty_bp: int = 1_000_000
    confidence_level: float = 0.95

    def __post_init__(self) -> None:
        if not 0.0 < self.detection_overlap_fraction <= 1.0:
            raise ValueError("detection_overlap_fraction must lie in (0, 1]")
        if not 0.0 < self.minimum_assessable_fraction <= 1.0:
            raise ValueError("minimum_assessable_fraction must lie in (0, 1]")
        if self.copy_number_tolerance < 0:
            raise ValueError("copy_number_tolerance must not be negative")
        if self.maximum_truth_boundary_uncertainty_bp < 0:
            raise ValueError("maximum_truth_boundary_uncertainty_bp must not be negative")


@dataclass(frozen=True)
class BoundaryUncertainty:
    """Half-width of the confidence interval around a truth event's boundaries.

    Cytogenetic truth resolves breakpoints only to a chromosome band, which can span
    several megabases; array truth resolves them to probe spacing; a simulated truth set
    knows them exactly. Recording the uncertainty is what allows breakpoint accuracy to
    be suppressed when the truth cannot support the measurement, instead of reporting a
    breakpoint error that is really just the width of a Giemsa band.
    """

    start_uncertainty_bp: int = 0
    end_uncertainty_bp: int = 0

    def __post_init__(self) -> None:
        if self.start_uncertainty_bp < 0 or self.end_uncertainty_bp < 0:
            raise ValueError("boundary uncertainty must not be negative")

    @property
    def maximum(self) -> int:
        return max(self.start_uncertainty_bp, self.end_uncertainty_bp)


@dataclass(frozen=True)
class EventRecord:
    """A maximal run of one altered state, derived from a segment list."""

    event_id: str
    contig: str
    start: int
    end: int
    state: CopyNumberState
    copy_number: float | None
    boundary: BoundaryUncertainty

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def size_class(self) -> SizeClass:
        return size_class(self.length)


@dataclass(frozen=True)
class TruthEventResult:
    """Scoring outcome for one truth event."""

    event: EventRecord
    total_bases: int
    evaluable_bases: int
    concordant_bases: int
    no_call_bases: int
    outcome: TruthOutcome
    reason: str
    concordant_fraction: float | None
    observed_start: int | None
    observed_end: int | None
    start_delta_bp: int | None
    end_delta_bp: int | None
    breakpoint_assessable: bool
    breakpoint_skip_reason: str | None


@dataclass(frozen=True)
class QueryEventResult:
    """Scoring outcome for one called event."""

    event: EventRecord
    total_bases: int
    evaluable_bases: int
    concordant_bases: int
    outcome: QueryOutcome
    reason: str
    concordant_fraction: float | None


@dataclass(frozen=True)
class BaseLevelResult:
    """Base-pair-weighted agreement over the evaluable genome."""

    evaluable_bases: int
    concordant_bases: int
    concordance: float | None
    #: ``(truth_state, query_state) -> bases``. The full matrix, not a binary summary.
    confusion: dict[tuple[CopyNumberState, CopyNumberState], int]
    truth_bases_by_state: dict[CopyNumberState, int]
    query_bases_by_state: dict[CopyNumberState, int]
    #: Per-state bp recall and precision. Descriptive only; see the note on intervals below.
    recall_by_state: dict[CopyNumberState, float | None]
    precision_by_state: dict[CopyNumberState, float | None]


@dataclass(frozen=True)
class GenomePartition:
    """How the reference genome was divided before scoring.

    These four numbers are the audit trail behind every metric. They make the difference
    between "the assay found nothing" and "the assay could not look" an explicit,
    inspectable quantity.
    """

    reference_bases: int
    evaluable_bases: int
    excluded_bases: int
    query_no_call_bases: int
    truth_silent_bases: int
    excluded_bases_by_reason: dict[str, int]


@dataclass(frozen=True)
class CopyNumberAccuracy:
    """Quantitative copy-number agreement where both sides report a number."""

    assessed_bases: int
    mean_absolute_error: float | None
    root_mean_square_error: float | None
    within_tolerance_bases: int
    within_tolerance_fraction: float | None


@dataclass(frozen=True)
class BreakpointAccuracy:
    """Breakpoint agreement, restricted to truth events that can support it."""

    assessed_events: int
    skipped_events: int
    skip_reasons: dict[str, int]
    median_absolute_start_delta_bp: int | None
    median_absolute_end_delta_bp: int | None
    maximum_absolute_delta_bp: int | None


@dataclass(frozen=True)
class StratumDetection:
    """Detection rate within one stratum, with an interval and its denominator."""

    label: str
    detected: int
    missed: int
    not_assessable: int
    detection_rate: ProportionEstimate


@dataclass(frozen=True)
class EvaluationResult:
    """The complete comparison of one call set against one truth set."""

    options: EvaluationOptions
    partition: GenomePartition
    base_level: BaseLevelResult
    truth_events: list[TruthEventResult]
    query_events: list[QueryEventResult]
    detection_rate: ProportionEstimate
    confirmation_rate: ProportionEstimate
    detection_by_size_class: list[StratumDetection]
    detection_by_state: list[StratumDetection]
    copy_number_accuracy: CopyNumberAccuracy
    breakpoint_accuracy: BreakpointAccuracy
    warnings: list[str] = field(default_factory=list)


def _validate_segments(segments: Sequence[StateSegment], label: str) -> dict[str, list[StateSegment]]:
    """Group segments by canonical contig and reject overlaps within a contig."""
    grouped: dict[str, list[StateSegment]] = {}
    for segment in segments:
        grouped.setdefault(canonical_contig(segment.contig), []).append(segment)
    for contig, items in grouped.items():
        items.sort(key=lambda item: (item.start, item.end))
        for previous, current in zip(items, items[1:], strict=False):
            if current.start < previous.end:
                raise ValueError(
                    f"{label} segments overlap on contig {contig}: "
                    f"[{previous.start}, {previous.end}) and [{current.start}, {current.end})"
                )
    return grouped


def _locate(items: Sequence[StateSegment], starts: Sequence[int], position: int) -> StateSegment | None:
    """Return the segment covering ``position``, or ``None``."""
    index = bisect_right(starts, position) - 1
    if index < 0:
        return None
    candidate = items[index]
    return candidate if candidate.start <= position < candidate.end else None


def _covered(spans: Sequence[Interval], starts: Sequence[int], position: int) -> bool:
    """Return whether ``position`` falls inside a normalized interval list."""
    index = bisect_right(starts, position) - 1
    if index < 0:
        return False
    start, end = spans[index]
    return start <= position < end


def _build_event(
    run: Sequence[StateSegment],
    contig: str,
    prefix: str,
    fallback: BoundaryUncertainty,
) -> EventRecord:
    """Fold a run of same-state segments into one event with a length-weighted copy number.

    The event's boundary uncertainty comes from the outer edges of the run: the first
    segment's start uncertainty and the last segment's end uncertainty. Interior
    boundaries vanish along with the segmentation that produced them.
    """
    start = run[0].start
    end = run[-1].end
    numbers = [item.copy_number for item in run if item.copy_number is not None]
    weights = [item.length for item in run if item.copy_number is not None]
    copy_number = (
        sum(number * weight for number, weight in zip(numbers, weights, strict=True))
        / sum(weights)
        if numbers
        else None
    )
    return EventRecord(
        event_id=f"{prefix}:{contig}:{start}-{end}:{run[0].state.value}",
        contig=contig,
        start=start,
        end=end,
        state=run[0].state,
        copy_number=copy_number,
        boundary=BoundaryUncertainty(
            start_uncertainty_bp=(
                run[0].start_uncertainty_bp
                if run[0].start_uncertainty_bp is not None
                else fallback.start_uncertainty_bp
            ),
            end_uncertainty_bp=(
                run[-1].end_uncertainty_bp
                if run[-1].end_uncertainty_bp is not None
                else fallback.end_uncertainty_bp
            ),
        ),
    )


def derive_events(
    segments: Sequence[StateSegment],
    *,
    prefix: str,
    default_boundary: BoundaryUncertainty | None = None,
) -> list[EventRecord]:
    """Collapse a segment list into maximal runs of one altered state.

    Adjacent segments carrying the same state are merged even when a caller emitted them
    separately, so that an event's identity does not depend on the caller's segmentation
    granularity. Neutral and no-call spans are not events.
    """
    grouped = _validate_segments(segments, prefix)
    fallback = default_boundary or BoundaryUncertainty()
    events: list[EventRecord] = []
    for contig in sorted(grouped):
        run: list[StateSegment] = []
        for segment in grouped[contig]:
            continues_run = (
                bool(run) and run[-1].end == segment.start and run[-1].state == segment.state
            )
            if not is_altered(segment.state) or not continues_run:
                if run:
                    events.append(_build_event(run, contig, prefix, fallback))
                    run = []
            if is_altered(segment.state):
                run.append(segment)
        if run:
            events.append(_build_event(run, contig, prefix, fallback))
    return events


def _state_at(
    grouped: Mapping[str, list[StateSegment]],
    starts: Mapping[str, list[int]],
    contig: str,
    position: int,
    background: CopyNumberState,
) -> tuple[CopyNumberState, float | None]:
    """Return the state and copy number asserted at a position, or the background."""
    items = grouped.get(contig)
    if not items:
        return background, None
    segment = _locate(items, starts[contig], position)
    if segment is None:
        return background, None
    return segment.state, segment.copy_number


def evaluate(
    *,
    truth_segments: Sequence[StateSegment],
    query_segments: Sequence[StateSegment],
    evaluable: Mapping[str, Sequence[Interval]],
    reference_bases: int,
    truth_background: CopyNumberState,
    query_background: CopyNumberState,
    options: EvaluationOptions | None = None,
    default_truth_boundary: BoundaryUncertainty | None = None,
    excluded_bases_by_reason: Mapping[str, int] | None = None,
    warnings: Iterable[str] = (),
) -> EvaluationResult:
    """Compare a call set against a truth set over an explicitly evaluable genome.

    ``truth_background`` and ``query_background`` state what a side means by *silence*.
    A SNP array asserts a neutral copy number wherever it has probes, so its background
    is ``NEUTRAL`` (closed world). A karyotype report or an alteration-only caller
    asserts nothing outside what it lists, so its background is ``NO_CALL`` (open world).
    Getting this wrong inverts specificity, which is why it has no default.
    """
    resolved = options or EvaluationOptions()
    truth_grouped = _validate_segments(truth_segments, "truth")
    query_grouped = _validate_segments(query_segments, "query")
    truth_starts = {contig: [s.start for s in items] for contig, items in truth_grouped.items()}
    query_starts = {contig: [s.start for s in items] for contig, items in query_grouped.items()}
    evaluable_set: IntervalSet = normalize_set(evaluable)
    evaluable_starts = {contig: [s for s, _ in spans] for contig, spans in evaluable_set.items()}

    confusion: dict[tuple[CopyNumberState, CopyNumberState], int] = {}
    truth_bases_by_state: dict[CopyNumberState, int] = {}
    query_bases_by_state: dict[CopyNumberState, int] = {}
    copy_number_pairs: list[tuple[float, float]] = []
    copy_number_weights: list[int] = []
    within_tolerance_bases = 0
    evaluable_bases = 0
    concordant_bases = 0
    query_no_call_bases = 0
    truth_silent_bases = 0

    contigs = sorted(set(truth_grouped) | set(query_grouped) | set(evaluable_set))
    # Per-contig accumulation of concordant bases, keyed by contig, used later to score
    # events without re-walking the partition.
    concordance_spans: dict[str, list[Interval]] = {}
    # Bases that actually entered the confusion matrix. This is narrower than the
    # observability mask: it additionally excludes bases where either side is silent.
    # Event scoring must use this set, otherwise a call in a region where an open-world
    # truth asserts nothing would be counted as a false positive rather than as
    # unassessable.
    scored_spans: dict[str, list[Interval]] = {}

    for contig in contigs:
        cuts: set[int] = set()
        for segment in truth_grouped.get(contig, []):
            cuts.add(segment.start)
            cuts.add(segment.end)
        for segment in query_grouped.get(contig, []):
            cuts.add(segment.start)
            cuts.add(segment.end)
        for start, end in evaluable_set.get(contig, []):
            cuts.add(start)
            cuts.add(end)
        ordered = sorted(cuts)
        for left, right in zip(ordered, ordered[1:], strict=False):
            if right <= left:
                continue
            width = right - left
            truth_state, truth_cn = _state_at(
                truth_grouped, truth_starts, contig, left, truth_background
            )
            query_state, query_cn = _state_at(
                query_grouped, query_starts, contig, left, query_background
            )
            inside = _covered(
                evaluable_set.get(contig, []), evaluable_starts.get(contig, []), left
            )
            if query_state == CopyNumberState.NO_CALL:
                query_no_call_bases += width
            if truth_state == CopyNumberState.NO_CALL:
                truth_silent_bases += width
            if not inside:
                continue
            # Bases where either side is silent are excluded from scoring by construction
            # of the evaluable set; this guard keeps the invariant explicit and local.
            if truth_state == CopyNumberState.NO_CALL or query_state == CopyNumberState.NO_CALL:
                continue
            evaluable_bases += width
            scored_spans.setdefault(contig, []).append((left, right))
            key = (truth_state, query_state)
            confusion[key] = confusion.get(key, 0) + width
            truth_bases_by_state[truth_state] = truth_bases_by_state.get(truth_state, 0) + width
            query_bases_by_state[query_state] = query_bases_by_state.get(query_state, 0) + width
            if concordant(truth_state, query_state, resolved.concordance_mode):
                concordant_bases += width
                concordance_spans.setdefault(contig, []).append((left, right))
            if truth_cn is not None and query_cn is not None:
                copy_number_pairs.append((truth_cn, query_cn))
                copy_number_weights.append(width)
                if abs(truth_cn - query_cn) <= resolved.copy_number_tolerance:
                    within_tolerance_bases += width

    recall_by_state = {
        state: (
            _concordant_bases_for_truth_state(confusion, state, resolved) / total
            if total
            else None
        )
        for state, total in truth_bases_by_state.items()
    }
    precision_by_state = {
        state: (
            _concordant_bases_for_query_state(confusion, state, resolved) / total if total else None
        )
        for state, total in query_bases_by_state.items()
    }

    base_level = BaseLevelResult(
        evaluable_bases=evaluable_bases,
        concordant_bases=concordant_bases,
        concordance=concordant_bases / evaluable_bases if evaluable_bases else None,
        confusion=confusion,
        truth_bases_by_state=truth_bases_by_state,
        query_bases_by_state=query_bases_by_state,
        recall_by_state=recall_by_state,
        precision_by_state=precision_by_state,
    )

    normalized_concordance = normalize_set(concordance_spans)
    scored_set = normalize_set(scored_spans)
    truth_events = derive_events(
        truth_segments,
        prefix="truth",
        default_boundary=default_truth_boundary,
    )
    query_events = derive_events(query_segments, prefix="query")

    truth_results = [
        _score_truth_event(
            event,
            scored_set,
            normalized_concordance,
            query_grouped,
            query_starts,
            resolved,
        )
        for event in truth_events
    ]
    query_results = [
        _score_query_event(event, scored_set, normalized_concordance, resolved)
        for event in query_events
    ]

    detected = sum(1 for item in truth_results if item.outcome == TruthOutcome.DETECTED)
    assessable_truth = sum(
        1 for item in truth_results if item.outcome != TruthOutcome.NOT_ASSESSABLE
    )
    confirmed = sum(1 for item in query_results if item.outcome == QueryOutcome.CONFIRMED)
    assessable_query = sum(
        1 for item in query_results if item.outcome != QueryOutcome.NOT_ASSESSABLE
    )

    collected_warnings = list(warnings)
    if evaluable_bases == 0:
        collected_warnings.append(
            "No evaluable bases remained after applying the observability mask; every "
            "metric is undefined rather than zero."
        )

    return EvaluationResult(
        options=resolved,
        partition=GenomePartition(
            reference_bases=reference_bases,
            evaluable_bases=evaluable_bases,
            excluded_bases=max(0, reference_bases - evaluable_bases),
            query_no_call_bases=query_no_call_bases,
            truth_silent_bases=truth_silent_bases,
            excluded_bases_by_reason=dict(excluded_bases_by_reason or {}),
        ),
        base_level=base_level,
        truth_events=truth_results,
        query_events=query_results,
        detection_rate=wilson_interval(
            detected, assessable_truth, confidence_level=resolved.confidence_level
        ),
        confirmation_rate=wilson_interval(
            confirmed, assessable_query, confidence_level=resolved.confidence_level
        ),
        detection_by_size_class=_stratify(
            truth_results, lambda item: item.event.size_class.value, resolved.confidence_level
        ),
        detection_by_state=_stratify(
            truth_results, lambda item: item.event.state.value, resolved.confidence_level
        ),
        copy_number_accuracy=CopyNumberAccuracy(
            assessed_bases=sum(copy_number_weights),
            mean_absolute_error=mean_absolute_error(copy_number_pairs, copy_number_weights),
            root_mean_square_error=root_mean_square_error(copy_number_pairs, copy_number_weights),
            within_tolerance_bases=within_tolerance_bases,
            within_tolerance_fraction=(
                within_tolerance_bases / sum(copy_number_weights)
                if sum(copy_number_weights)
                else None
            ),
        ),
        breakpoint_accuracy=_breakpoint_accuracy(truth_results),
        warnings=collected_warnings,
    )


def _concordant_bases_for_truth_state(
    confusion: Mapping[tuple[CopyNumberState, CopyNumberState], int],
    state: CopyNumberState,
    options: EvaluationOptions,
) -> int:
    return sum(
        bases
        for (truth_state, query_state), bases in confusion.items()
        if truth_state == state and concordant(truth_state, query_state, options.concordance_mode)
    )


def _concordant_bases_for_query_state(
    confusion: Mapping[tuple[CopyNumberState, CopyNumberState], int],
    state: CopyNumberState,
    options: EvaluationOptions,
) -> int:
    return sum(
        bases
        for (truth_state, query_state), bases in confusion.items()
        if query_state == state and concordant(truth_state, query_state, options.concordance_mode)
    )


def _span_overlap(
    contig: str, span: Interval, other: Mapping[str, Sequence[Interval]]
) -> list[Interval]:
    spans = other.get(canonical_contig(contig))
    if not spans:
        return []
    return intersect([span], spans)


def _score_truth_event(
    event: EventRecord,
    evaluable: IntervalSet,
    concordance: IntervalSet,
    query_grouped: Mapping[str, list[StateSegment]],
    query_starts: Mapping[str, list[int]],
    options: EvaluationOptions,
) -> TruthEventResult:
    span = (event.start, event.end)
    total = event.length
    evaluable_spans = _span_overlap(event.contig, span, evaluable)
    evaluable_bases = sum(end - start for start, end in evaluable_spans)
    concordant_spans = _span_overlap(event.contig, span, concordance)
    concordant_bases = sum(end - start for start, end in concordant_spans)
    no_call_bases = total - evaluable_bases

    if total == 0 or evaluable_bases / total < options.minimum_assessable_fraction:
        return TruthEventResult(
            event=event,
            total_bases=total,
            evaluable_bases=evaluable_bases,
            concordant_bases=concordant_bases,
            no_call_bases=no_call_bases,
            outcome=TruthOutcome.NOT_ASSESSABLE,
            reason=(
                f"only {evaluable_bases}/{total} bases were observable, below the "
                f"minimum assessable fraction {options.minimum_assessable_fraction}"
            ),
            concordant_fraction=None,
            observed_start=None,
            observed_end=None,
            start_delta_bp=None,
            end_delta_bp=None,
            breakpoint_assessable=False,
            breakpoint_skip_reason="event_not_assessable",
        )

    fraction = concordant_bases / evaluable_bases if evaluable_bases else 0.0
    detected = fraction >= options.detection_overlap_fraction
    observed_start: int | None = None
    observed_end: int | None = None
    start_delta: int | None = None
    end_delta: int | None = None
    skip_reason: str | None = None
    assessable_breakpoint = False

    if detected and concordant_spans:
        observed_start = concordant_spans[0][0]
        observed_end = concordant_spans[-1][1]
        observed_start, observed_end = _extend_to_called_event(
            event, observed_start, observed_end, query_grouped, query_starts, options
        )
        if event.boundary.maximum > options.maximum_truth_boundary_uncertainty_bp:
            skip_reason = "truth_boundary_resolution_insufficient"
        else:
            assessable_breakpoint = True
            start_delta = observed_start - event.start
            end_delta = observed_end - event.end
    elif detected:
        skip_reason = "no_concordant_span"
    else:
        skip_reason = "event_not_detected"

    return TruthEventResult(
        event=event,
        total_bases=total,
        evaluable_bases=evaluable_bases,
        concordant_bases=concordant_bases,
        no_call_bases=no_call_bases,
        outcome=TruthOutcome.DETECTED if detected else TruthOutcome.MISSED,
        reason=(
            f"{concordant_bases}/{evaluable_bases} observable bases were concordant "
            f"(threshold {options.detection_overlap_fraction})"
        ),
        concordant_fraction=fraction,
        observed_start=observed_start,
        observed_end=observed_end,
        start_delta_bp=start_delta,
        end_delta_bp=end_delta,
        breakpoint_assessable=assessable_breakpoint,
        breakpoint_skip_reason=skip_reason,
    )


def _extend_to_called_event(
    event: EventRecord,
    observed_start: int,
    observed_end: int,
    query_grouped: Mapping[str, list[StateSegment]],
    query_starts: Mapping[str, list[int]],
    options: EvaluationOptions,
) -> tuple[int, int]:
    """Widen the observed span to the full extent of the concordant called segments.

    Concordant spans are clipped to the truth event and to the evaluable mask, so using
    them directly would report a breakpoint error of zero whenever a caller overshoots.
    Walking outward through contiguous concordant called segments recovers the caller's
    actual boundary, which is what breakpoint accuracy is supposed to measure.
    """
    segments = query_grouped.get(canonical_contig(event.contig))
    if not segments:
        return observed_start, observed_end
    starts = query_starts[canonical_contig(event.contig)]

    start = observed_start
    index = bisect_right(starts, start) - 1
    while index >= 0:
        candidate = segments[index]
        if candidate.end < start or not concordant(
            event.state, candidate.state, options.concordance_mode
        ):
            break
        start = candidate.start
        index -= 1

    end = observed_end
    index = bisect_right(starts, max(end - 1, 0)) - 1
    while 0 <= index < len(segments):
        candidate = segments[index]
        if candidate.start > end or not concordant(
            event.state, candidate.state, options.concordance_mode
        ):
            break
        end = candidate.end
        index += 1
    return start, end


def _score_query_event(
    event: EventRecord,
    evaluable: IntervalSet,
    concordance: IntervalSet,
    options: EvaluationOptions,
) -> QueryEventResult:
    span = (event.start, event.end)
    total = event.length
    evaluable_bases = sum(
        end - start for start, end in _span_overlap(event.contig, span, evaluable)
    )
    concordant_bases = sum(
        end - start for start, end in _span_overlap(event.contig, span, concordance)
    )
    if total == 0 or evaluable_bases / total < options.minimum_assessable_fraction:
        return QueryEventResult(
            event=event,
            total_bases=total,
            evaluable_bases=evaluable_bases,
            concordant_bases=concordant_bases,
            outcome=QueryOutcome.NOT_ASSESSABLE,
            reason=(
                f"only {evaluable_bases}/{total} bases were observable; the call cannot "
                "be scored as correct or incorrect"
            ),
            concordant_fraction=None,
        )
    fraction = concordant_bases / evaluable_bases if evaluable_bases else 0.0
    confirmed = fraction >= options.detection_overlap_fraction
    return QueryEventResult(
        event=event,
        total_bases=total,
        evaluable_bases=evaluable_bases,
        concordant_bases=concordant_bases,
        outcome=QueryOutcome.CONFIRMED if confirmed else QueryOutcome.UNCONFIRMED,
        reason=(
            f"{concordant_bases}/{evaluable_bases} observable bases agreed with truth "
            f"(threshold {options.detection_overlap_fraction})"
        ),
        concordant_fraction=fraction,
    )


def _stratify(
    results: Sequence[TruthEventResult],
    key: Callable[[TruthEventResult], str],
    confidence_level: float,
) -> list[StratumDetection]:
    """Group truth-event outcomes by a label function and compute per-group rates."""
    buckets: dict[str, list[TruthEventResult]] = {}
    for item in results:
        buckets.setdefault(key(item), []).append(item)
    strata: list[StratumDetection] = []
    for label in sorted(buckets):
        items = buckets[label]
        detected = sum(1 for item in items if item.outcome == TruthOutcome.DETECTED)
        missed = sum(1 for item in items if item.outcome == TruthOutcome.MISSED)
        not_assessable = sum(
            1 for item in items if item.outcome == TruthOutcome.NOT_ASSESSABLE
        )
        strata.append(
            StratumDetection(
                label=label,
                detected=detected,
                missed=missed,
                not_assessable=not_assessable,
                detection_rate=wilson_interval(
                    detected, detected + missed, confidence_level=confidence_level
                ),
            )
        )
    return strata


def _median(values: Sequence[int]) -> int | None:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) // 2


def _breakpoint_accuracy(results: Sequence[TruthEventResult]) -> BreakpointAccuracy:
    assessed = [item for item in results if item.breakpoint_assessable]
    skipped = [item for item in results if not item.breakpoint_assessable]
    reasons: dict[str, int] = {}
    for item in skipped:
        reason = item.breakpoint_skip_reason or "unspecified"
        reasons[reason] = reasons.get(reason, 0) + 1
    start_deltas = [abs(item.start_delta_bp) for item in assessed if item.start_delta_bp is not None]
    end_deltas = [abs(item.end_delta_bp) for item in assessed if item.end_delta_bp is not None]
    combined = start_deltas + end_deltas
    return BreakpointAccuracy(
        assessed_events=len(assessed),
        skipped_events=len(skipped),
        skip_reasons=reasons,
        median_absolute_start_delta_bp=_median(start_deltas),
        median_absolute_end_delta_bp=_median(end_deltas),
        maximum_absolute_delta_bp=max(combined) if combined else None,
    )
