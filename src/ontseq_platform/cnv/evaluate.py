"""Orchestration between the CNV contract layer and the comparison core.

This module owns exactly one responsibility: turn a :class:`CnvBenchmarkCase` into an
evaluable genome plus two segment lists, hand them to :func:`ontseq_platform.cnv.core.evaluate`,
and project the result back into the contract. It contains no scoring logic of its own,
which keeps the science in one testable place.
"""

from __future__ import annotations

from collections.abc import Sequence

from .core import (
    BoundaryUncertainty,
    EvaluationOptions,
    EvaluationResult,
    StratumDetection,
    TruthEventResult,
    QueryEventResult,
)
from .core import evaluate as evaluate_core
from .intervals import Interval, IntervalSet, canonical_contig, normalize_set, subtract_set
from .mask import (
    ExclusionReason,
    ExclusionTrack,
    ObservabilityMask,
    build_mask,
    scope_from_intervals,
)
from .models import (
    BaseLevelReport,
    BreakpointAccuracyReport,
    CnvBenchmarkCase,
    CnvCallSet,
    CnvEvaluationOptions,
    CnvEvaluationReport,
    CopyNumberAccuracyReport,
    EventOutcomeReport,
    GenomePartitionReport,
    GenomicRegion,
    ProportionResult,
    StratumDetectionReport,
)
from .stats import ProportionEstimate

STANDING_LIMITATIONS = [
    "This comparison is an engineering control against a declared truth set. It is not "
    "analytical or clinical validation of any CNV method.",
    "Every reported figure depends on the observability mask, the concordance mode and "
    "the detection thresholds echoed in this report.",
    "Base-level agreement is reported without confidence intervals because base pairs "
    "within a segment are not independent observations.",
    "Absence of a call inside an excluded or no-call region is not a negative biological "
    "finding.",
]


def _regions_to_set(regions: Sequence[GenomicRegion]) -> IntervalSet:
    grouped: dict[str, list[Interval]] = {}
    for region in regions:
        grouped.setdefault(canonical_contig(region.contig), []).append((region.start, region.end))
    return normalize_set(grouped)


def _proportion(estimate: ProportionEstimate) -> ProportionResult:
    return ProportionResult(
        successes=estimate.successes,
        total=estimate.total,
        point=estimate.point,
        lower=estimate.lower,
        upper=estimate.upper,
        confidence_level=estimate.confidence_level,
    )


def _stratum(item: StratumDetection) -> StratumDetectionReport:
    return StratumDetectionReport(
        label=item.label,
        detected=item.detected,
        missed=item.missed,
        not_assessable=item.not_assessable,
        detection_rate=_proportion(item.detection_rate),
    )


def _truth_event_report(item: TruthEventResult) -> EventOutcomeReport:
    return EventOutcomeReport(
        event_id=item.event.event_id,
        contig=item.event.contig,
        start=item.event.start,
        end=item.event.end,
        state=item.event.state,
        size_class=item.event.size_class.value,
        length_bp=item.event.length,
        outcome=item.outcome.value,
        reason=item.reason,
        total_bases=item.total_bases,
        evaluable_bases=item.evaluable_bases,
        concordant_bases=item.concordant_bases,
        concordant_fraction=item.concordant_fraction,
        start_delta_bp=item.start_delta_bp,
        end_delta_bp=item.end_delta_bp,
        breakpoint_skip_reason=item.breakpoint_skip_reason,
    )


def _query_event_report(item: QueryEventResult) -> EventOutcomeReport:
    return EventOutcomeReport(
        event_id=item.event.event_id,
        contig=item.event.contig,
        start=item.event.start,
        end=item.event.end,
        state=item.event.state,
        size_class=item.event.size_class.value,
        length_bp=item.event.length,
        outcome=item.outcome.value,
        reason=item.reason,
        total_bases=item.total_bases,
        evaluable_bases=item.evaluable_bases,
        concordant_bases=item.concordant_bases,
        concordant_fraction=item.concordant_fraction,
        start_delta_bp=None,
        end_delta_bp=None,
        breakpoint_skip_reason=None,
    )


def _core_options(options: CnvEvaluationOptions) -> EvaluationOptions:
    return EvaluationOptions(
        concordance_mode=options.concordance_mode,
        detection_overlap_fraction=options.detection_overlap_fraction,
        minimum_assessable_fraction=options.minimum_assessable_fraction,
        copy_number_tolerance=options.copy_number_tolerance,
        maximum_truth_boundary_uncertainty_bp=options.maximum_truth_boundary_uncertainty_bp,
        confidence_level=options.confidence_level,
    )


def build_case_mask(case: CnvBenchmarkCase) -> ObservabilityMask:
    """Assemble the evaluable genome for a benchmark case.

    Exclusion order matters only for attribution, not for the resulting region: bases
    removed by several tracks are attributed to the first track that removed them, so
    the per-reason counts sum exactly to the total removed.
    """
    scope: IntervalSet | None = None
    if case.analysis_scope:
        scope = scope_from_intervals(
            _regions_to_set(case.analysis_scope),
            flank_bp=case.options.analysis_scope_flank_bp,
        )

    tracks: list[ExclusionTrack] = []
    if case.excluded_regions:
        tracks.append(
            ExclusionTrack(
                reason=ExclusionReason.BLACKLIST,
                intervals=_regions_to_set(case.excluded_regions),
                source=f"case:{case.case_id}:excluded_regions",
            )
        )
    if case.truth.uninformative_regions:
        tracks.append(
            ExclusionTrack(
                reason=ExclusionReason.TRUTH_NOT_INFORMATIVE,
                intervals=_regions_to_set(case.truth.uninformative_regions),
                source=f"truth:{case.truth.truth_id}:uninformative",
            )
        )
    if case.truth.informative_regions:
        # Everything outside the truth's informative footprint is unusable for scoring,
        # regardless of what the caller reported there.
        genome: dict[str, list[Interval]] = {
            canonical_contig(contig): [(0, length)]
            for contig, length in case.contig_lengths.items()
        }
        outside = subtract_set(genome, _regions_to_set(case.truth.informative_regions))
        tracks.append(
            ExclusionTrack(
                reason=ExclusionReason.TRUTH_NOT_INFORMATIVE,
                intervals=outside,
                source=f"truth:{case.truth.truth_id}:outside_informative_regions",
            )
        )
    if case.call_set.no_call_regions:
        tracks.append(
            ExclusionTrack(
                reason=ExclusionReason.CALLER_NO_CALL,
                intervals=_regions_to_set(case.call_set.no_call_regions),
                source=f"call_set:{case.call_set.call_set_id}:no_call",
            )
        )
    return build_mask(
        contig_lengths=dict(case.contig_lengths),
        analysis_scope=scope,
        tracks=tracks,
    )


def evaluate_case(case: CnvBenchmarkCase, *, evaluation_id: str | None = None) -> CnvEvaluationReport:
    """Score a benchmark case and return the full auditable report."""
    mask = build_case_mask(case)
    # Segments that declare no uncertainty of their own inherit the truth set's declared
    # resolution, which is the weakest defensible assumption about their breakpoints.
    default_boundary = BoundaryUncertainty(
        start_uncertainty_bp=case.truth.resolution_bp,
        end_uncertainty_bp=case.truth.resolution_bp,
    )
    result = evaluate_core(
        truth_segments=case.truth.state_segments(),
        query_segments=case.call_set.state_segments(),
        evaluable=mask.evaluable,
        reference_bases=mask.reference_bases,
        truth_background=case.truth.background_state,
        query_background=case.call_set.background_state,
        options=_core_options(case.options),
        default_truth_boundary=default_boundary,
        excluded_bases_by_reason=mask.excluded_bases_by_reason,
        warnings=mask.warnings,
    )
    return _to_report(case, result, mask, evaluation_id=evaluation_id)


def _to_report(
    case: CnvBenchmarkCase,
    result: EvaluationResult,
    mask: ObservabilityMask,
    *,
    evaluation_id: str | None,
) -> CnvEvaluationReport:
    warnings = list(result.warnings)
    if case.truth.resolution_bp > 0:
        small = [
            item
            for item in result.query_events
            if item.event.length < case.truth.resolution_bp
        ]
        if small:
            warnings.append(
                f"{len(small)} called event(s) are smaller than the truth set's declared "
                f"resolution of {case.truth.resolution_bp} bp. The truth cannot confirm "
                "or refute them, so they must not be read as false positives."
            )
    if case.call_set.status != "COMPLETED":
        warnings.append(
            f"The call set status is {case.call_set.status}; metrics describe whatever "
            "partial output was present."
        )

    return CnvEvaluationReport(
        evaluation_id=evaluation_id or f"{case.case_id}:{case.call_set.call_set_id}",
        sample_id=case.truth.sample_id,
        truth_id=case.truth.truth_id,
        call_set_id=case.call_set.call_set_id,
        genome_build=case.genome_build,
        method=case.call_set.method,
        method_version=case.call_set.method_version,
        truth_source=case.truth.source,
        data_basis=case.call_set.data_basis,
        strata=case.strata,
        options=case.options,
        partition=GenomePartitionReport(
            reference_bases=result.partition.reference_bases,
            evaluable_bases=result.partition.evaluable_bases,
            excluded_bases=result.partition.excluded_bases,
            query_no_call_bases=result.partition.query_no_call_bases,
            truth_silent_bases=result.partition.truth_silent_bases,
            evaluable_fraction=(
                result.partition.evaluable_bases / result.partition.reference_bases
                if result.partition.reference_bases
                else None
            ),
            excluded_bases_by_reason=dict(mask.excluded_bases_by_reason),
        ),
        base_level=BaseLevelReport(
            evaluable_bases=result.base_level.evaluable_bases,
            concordant_bases=result.base_level.concordant_bases,
            concordance=result.base_level.concordance,
            confusion={
                f"{truth.value}->{query.value}": bases
                for (truth, query), bases in sorted(
                    result.base_level.confusion.items(),
                    key=lambda item: (item[0][0].value, item[0][1].value),
                )
            },
            truth_bases_by_state={
                state.value: bases
                for state, bases in sorted(
                    result.base_level.truth_bases_by_state.items(), key=lambda i: i[0].value
                )
            },
            query_bases_by_state={
                state.value: bases
                for state, bases in sorted(
                    result.base_level.query_bases_by_state.items(), key=lambda i: i[0].value
                )
            },
            recall_by_state={
                state.value: value
                for state, value in sorted(
                    result.base_level.recall_by_state.items(), key=lambda i: i[0].value
                )
            },
            precision_by_state={
                state.value: value
                for state, value in sorted(
                    result.base_level.precision_by_state.items(), key=lambda i: i[0].value
                )
            },
        ),
        detection_rate=_proportion(result.detection_rate),
        confirmation_rate=_proportion(result.confirmation_rate),
        detection_by_size_class=[_stratum(item) for item in result.detection_by_size_class],
        detection_by_state=[_stratum(item) for item in result.detection_by_state],
        copy_number_accuracy=CopyNumberAccuracyReport(
            assessed_bases=result.copy_number_accuracy.assessed_bases,
            mean_absolute_error=result.copy_number_accuracy.mean_absolute_error,
            root_mean_square_error=result.copy_number_accuracy.root_mean_square_error,
            within_tolerance_bases=result.copy_number_accuracy.within_tolerance_bases,
            within_tolerance_fraction=result.copy_number_accuracy.within_tolerance_fraction,
        ),
        breakpoint_accuracy=BreakpointAccuracyReport(
            assessed_events=result.breakpoint_accuracy.assessed_events,
            skipped_events=result.breakpoint_accuracy.skipped_events,
            skip_reasons=dict(sorted(result.breakpoint_accuracy.skip_reasons.items())),
            median_absolute_start_delta_bp=result.breakpoint_accuracy.median_absolute_start_delta_bp,
            median_absolute_end_delta_bp=result.breakpoint_accuracy.median_absolute_end_delta_bp,
            maximum_absolute_delta_bp=result.breakpoint_accuracy.maximum_absolute_delta_bp,
        ),
        truth_events=[_truth_event_report(item) for item in result.truth_events],
        query_events=[_query_event_report(item) for item in result.query_events],
        warnings=warnings,
        limitations=list(STANDING_LIMITATIONS),
    )


def compare_methods(
    case_template: CnvBenchmarkCase, call_sets: Sequence[CnvCallSet]
) -> list[CnvEvaluationReport]:
    """Score several methods against one truth set on one shared evaluable genome.

    Each method's own no-call regions are removed from the *shared* mask before any
    method is scored. Without that step a method is rewarded for declining to call in
    difficult regions, because its blind spots shrink only its own denominator.
    """
    shared_no_call: dict[str, list[Interval]] = {}
    for call_set in call_sets:
        for region in call_set.no_call_regions:
            shared_no_call.setdefault(canonical_contig(region.contig), []).append(
                (region.start, region.end)
            )
    shared_regions = [
        GenomicRegion(contig=contig, start=start, end=end, label="shared_no_call")
        for contig, spans in sorted(normalize_set(shared_no_call).items())
        for start, end in spans
    ]
    reports: list[CnvEvaluationReport] = []
    for call_set in call_sets:
        case = case_template.model_copy(
            update={
                "call_set": call_set.model_copy(update={"no_call_regions": shared_regions}),
            }
        )
        reports.append(
            evaluate_case(case, evaluation_id=f"{case_template.case_id}:{call_set.call_set_id}")
        )
    return reports
