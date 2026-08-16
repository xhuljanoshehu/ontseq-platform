from __future__ import annotations

import unittest

from ontseq_platform.cnv.core import (
    BoundaryUncertainty,
    EvaluationOptions,
    QueryOutcome,
    StateSegment,
    TruthOutcome,
    derive_events,
    evaluate,
    size_class,
)
from ontseq_platform.cnv.states import ConcordanceMode, CopyNumberState

CHR5 = 181_538_259


def _evaluate(truth, query, evaluable=None, **kwargs):
    return evaluate(
        truth_segments=truth,
        query_segments=query,
        evaluable=evaluable if evaluable is not None else {"5": [(0, CHR5)]},
        reference_bases=kwargs.pop("reference_bases", CHR5),
        truth_background=kwargs.pop("truth_background", CopyNumberState.NEUTRAL),
        query_background=kwargs.pop("query_background", CopyNumberState.NEUTRAL),
        **kwargs,
    )


class SegmentationIndependenceTests(unittest.TestCase):
    """The central property: scoring must not depend on how either side segmented."""

    def test_fragmented_calls_are_not_penalised(self) -> None:
        truth = [StateSegment("chr5", 70_000_000, 160_000_000, CopyNumberState.LOSS, 1.0)]
        fragmented = [
            StateSegment("chr5", 70_000_000, 100_000_000, CopyNumberState.LOSS, 1.0),
            StateSegment("chr5", 100_000_000, 130_000_000, CopyNumberState.LOSS, 1.0),
            StateSegment("chr5", 130_000_000, 160_000_000, CopyNumberState.LOSS, 1.0),
        ]
        whole = [StateSegment("chr5", 70_000_000, 160_000_000, CopyNumberState.LOSS, 1.0)]

        fragmented_report = _evaluate(truth, fragmented)
        whole_report = _evaluate(truth, whole)

        self.assertEqual(fragmented_report.detection_rate.point, 1.0)
        self.assertEqual(fragmented_report.base_level.concordance, 1.0)
        # Identical scores despite completely different segmentation.
        self.assertEqual(
            fragmented_report.detection_rate.successes, whole_report.detection_rate.successes
        )
        self.assertEqual(
            fragmented_report.base_level.concordant_bases,
            whole_report.base_level.concordant_bases,
        )
        # Three adjacent same-state calls collapse into one event, not three.
        self.assertEqual(len(fragmented_report.query_events), 1)
        self.assertEqual(fragmented_report.query_events[0].outcome, QueryOutcome.CONFIRMED)

    def test_contig_prefix_mismatch_does_not_destroy_recall(self) -> None:
        report = _evaluate(
            [StateSegment("chr5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
        )
        self.assertEqual(report.detection_rate.point, 1.0)

    def test_partial_detection_below_threshold_is_a_miss(self) -> None:
        report = _evaluate(
            [StateSegment("5", 0, 30_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 0, 9_000_000, CopyNumberState.LOSS, 1.0)],
        )
        self.assertEqual(report.truth_events[0].outcome, TruthOutcome.MISSED)
        self.assertAlmostEqual(report.truth_events[0].concordant_fraction, 0.3)


class ObservabilityTests(unittest.TestCase):
    """No-call, exclusion and biological negativity must remain distinguishable."""

    def test_event_in_unobservable_region_is_not_a_false_negative(self) -> None:
        report = _evaluate(
            [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            [],
            evaluable={"5": [(9_000_000, CHR5)]},
        )
        self.assertEqual(report.truth_events[0].outcome, TruthOutcome.NOT_ASSESSABLE)
        # An unassessable event is excluded from the denominator entirely.
        self.assertEqual(report.detection_rate.total, 0)
        self.assertIsNone(report.detection_rate.point)

    def test_open_world_truth_does_not_manufacture_false_positives(self) -> None:
        report = evaluate(
            truth_segments=[StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            query_segments=[
                StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0),
                StateSegment("8", 0, 5_000_000, CopyNumberState.GAIN, 3.0),
            ],
            evaluable={"5": [(0, CHR5)], "8": [(0, 145_138_636)]},
            reference_bases=CHR5 + 145_138_636,
            truth_background=CopyNumberState.NO_CALL,
            query_background=CopyNumberState.NEUTRAL,
        )
        outcomes = {item.event.contig: item.outcome for item in report.query_events}
        self.assertEqual(outcomes["5"], QueryOutcome.CONFIRMED)
        # Truth asserts nothing on chr8, so the call there cannot be wrong.
        self.assertEqual(outcomes["8"], QueryOutcome.NOT_ASSESSABLE)
        self.assertEqual(report.confirmation_rate.point, 1.0)
        self.assertGreater(report.partition.truth_silent_bases, 0)

    def test_closed_world_truth_does_count_a_false_positive(self) -> None:
        report = evaluate(
            truth_segments=[StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            query_segments=[
                StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0),
                StateSegment("8", 0, 5_000_000, CopyNumberState.GAIN, 3.0),
            ],
            evaluable={"5": [(0, CHR5)], "8": [(0, 145_138_636)]},
            reference_bases=CHR5 + 145_138_636,
            truth_background=CopyNumberState.NEUTRAL,
            query_background=CopyNumberState.NEUTRAL,
        )
        outcomes = {item.event.contig: item.outcome for item in report.query_events}
        self.assertEqual(outcomes["8"], QueryOutcome.UNCONFIRMED)
        self.assertEqual(report.confirmation_rate.point, 0.5)

    def test_partition_reconciles_exactly(self) -> None:
        """The audit trail is only an audit trail if the numbers add up."""
        report = evaluate(
            truth_segments=[StateSegment("5", 0, 20_000_000, CopyNumberState.LOSS, 1.0)],
            query_segments=[
                StateSegment("5", 0, 5_000_000, CopyNumberState.LOSS, 1.0),
                StateSegment("5", 5_000_000, 12_000_000, CopyNumberState.NO_CALL),
            ],
            evaluable={"5": [(0, 30_000_000)]},
            reference_bases=CHR5,
            truth_background=CopyNumberState.NO_CALL,
            query_background=CopyNumberState.NEUTRAL,
        )
        partition = report.partition
        self.assertEqual(
            partition.mask_bases,
            partition.evaluable_bases
            + partition.truth_silent_bases
            + partition.query_no_call_bases,
        )
        self.assertEqual(
            partition.reference_bases, partition.mask_bases + partition.excluded_bases
        )
        # The mask allowed 30 Mb; truth is open-world so it is silent beyond 20 Mb, and
        # the caller explicitly declined between 5 Mb and 12 Mb.
        self.assertEqual(partition.mask_bases, 30_000_000)
        self.assertEqual(partition.truth_silent_bases, 10_000_000)
        self.assertEqual(partition.query_no_call_bases, 7_000_000)
        self.assertEqual(partition.evaluable_bases, 13_000_000)

    def test_no_call_bases_are_counted_only_inside_the_mask(self) -> None:
        """A caller declining outside the mask is already excluded, not double counted."""
        report = evaluate(
            truth_segments=[StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            query_segments=[StateSegment("5", 20_000_000, 30_000_000, CopyNumberState.NO_CALL)],
            evaluable={"5": [(0, 10_000_000)]},
            reference_bases=CHR5,
            truth_background=CopyNumberState.NEUTRAL,
            query_background=CopyNumberState.NEUTRAL,
        )
        self.assertEqual(report.partition.query_no_call_bases, 0)
        self.assertEqual(report.partition.mask_bases, 10_000_000)

    def test_empty_evaluable_genome_yields_undefined_not_zero(self) -> None:
        report = _evaluate(
            [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            evaluable={},
        )
        self.assertIsNone(report.base_level.concordance)
        self.assertIsNone(report.detection_rate.point)
        self.assertTrue(any("undefined" in warning for warning in report.warnings))


class BreakpointResolutionTests(unittest.TestCase):
    """Breakpoint accuracy must be withheld when the truth cannot support it."""

    def test_band_resolution_truth_suppresses_breakpoint_metrics(self) -> None:
        report = _evaluate(
            [StateSegment("5", 0, 30_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 1_000_000, 29_000_000, CopyNumberState.LOSS, 1.0)],
            default_truth_boundary=BoundaryUncertainty(5_000_000, 5_000_000),
            options=EvaluationOptions(maximum_truth_boundary_uncertainty_bp=1_000_000),
        )
        self.assertEqual(report.truth_events[0].outcome, TruthOutcome.DETECTED)
        self.assertIsNone(report.truth_events[0].start_delta_bp)
        self.assertEqual(
            report.breakpoint_accuracy.skip_reasons,
            {"truth_boundary_resolution_insufficient": 1},
        )
        self.assertEqual(report.breakpoint_accuracy.assessed_events, 0)

    def test_exact_truth_reports_breakpoint_deltas(self) -> None:
        report = _evaluate(
            [StateSegment("5", 10_000_000, 30_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 10_500_000, 30_200_000, CopyNumberState.LOSS, 1.0)],
        )
        self.assertEqual(report.breakpoint_accuracy.assessed_events, 1)
        self.assertEqual(report.truth_events[0].start_delta_bp, 500_000)
        self.assertEqual(report.truth_events[0].end_delta_bp, 200_000)

    def test_per_segment_uncertainty_survives_event_merging(self) -> None:
        """Uncertainty rides on the segment, so merging keeps the outer edges."""
        truth = [
            StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0, 8_000_000, 0),
            StateSegment("5", 10_000_000, 20_000_000, CopyNumberState.LOSS, 1.0, 0, 100),
        ]
        events = derive_events(truth, prefix="truth")
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].boundary.start_uncertainty_bp, 8_000_000)
        self.assertEqual(events[0].boundary.end_uncertainty_bp, 100)

    def test_segment_uncertainty_overrides_the_default(self) -> None:
        report = _evaluate(
            [StateSegment("5", 0, 30_000_000, CopyNumberState.LOSS, 1.0, 10, 10)],
            [StateSegment("5", 100, 30_000_000, CopyNumberState.LOSS, 1.0)],
            default_truth_boundary=BoundaryUncertainty(9_000_000, 9_000_000),
            options=EvaluationOptions(maximum_truth_boundary_uncertainty_bp=1_000_000),
        )
        # The segment declares tight boundaries, so the metric is reported despite the
        # permissive default.
        self.assertEqual(report.breakpoint_accuracy.assessed_events, 1)
        self.assertEqual(report.truth_events[0].start_delta_bp, 100)

    def test_caller_declared_uncertainty_cannot_suppress_the_metric(self) -> None:
        """A call set must not be able to excuse its own breakpoint error."""
        report = _evaluate(
            [StateSegment("5", 0, 30_000_000, CopyNumberState.LOSS, 1.0, 0, 0)],
            [
                StateSegment(
                    "5", 5_000_000, 30_000_000, CopyNumberState.LOSS, 1.0, 99_000_000, 99_000_000
                )
            ],
        )
        self.assertEqual(report.breakpoint_accuracy.assessed_events, 1)
        self.assertEqual(report.truth_events[0].start_delta_bp, 5_000_000)

    def test_overshooting_caller_is_measured_beyond_the_truth_event(self) -> None:
        """A call wider than truth must show a real error, not a clipped zero."""
        report = _evaluate(
            [StateSegment("5", 10_000_000, 20_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 5_000_000, 25_000_000, CopyNumberState.LOSS, 1.0)],
        )
        self.assertEqual(report.truth_events[0].start_delta_bp, -5_000_000)
        self.assertEqual(report.truth_events[0].end_delta_bp, 5_000_000)


class ConcordanceModeTests(unittest.TestCase):
    def test_directional_mode_accepts_a_loss_called_as_homozygous_loss(self) -> None:
        report = _evaluate(
            [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 0, 10_000_000, CopyNumberState.HOMOZYGOUS_LOSS, 0.0)],
        )
        self.assertEqual(report.detection_rate.point, 1.0)

    def test_strict_mode_rejects_the_same_pair(self) -> None:
        report = _evaluate(
            [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 0, 10_000_000, CopyNumberState.HOMOZYGOUS_LOSS, 0.0)],
            options=EvaluationOptions(concordance_mode=ConcordanceMode.STRICT),
        )
        self.assertEqual(report.detection_rate.point, 0.0)

    def test_copy_neutral_loh_is_never_collapsed_into_neutral(self) -> None:
        """A dosage-only caller must not get credit for detecting CN-LOH."""
        report = _evaluate(
            [StateSegment("5", 0, 10_000_000, CopyNumberState.COPY_NEUTRAL_LOH, 2.0)],
            [StateSegment("5", 0, 10_000_000, CopyNumberState.NEUTRAL, 2.0)],
        )
        self.assertEqual(report.truth_events[0].outcome, TruthOutcome.MISSED)


class QuantitativeTests(unittest.TestCase):
    def test_copy_number_error_is_base_weighted(self) -> None:
        report = _evaluate(
            [
                StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0),
                StateSegment("5", 10_000_000, 20_000_000, CopyNumberState.NEUTRAL, 2.0),
            ],
            [
                StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.4),
                StateSegment("5", 10_000_000, 20_000_000, CopyNumberState.NEUTRAL, 2.0),
            ],
        )
        # Half the assessed bases carry an error of 0.4, the other half zero.
        self.assertAlmostEqual(report.copy_number_accuracy.mean_absolute_error, 0.2)
        self.assertEqual(report.copy_number_accuracy.within_tolerance_fraction, 1.0)

    def test_size_class_boundaries(self) -> None:
        self.assertEqual(size_class(99_999).value, "lt_100kb")
        self.assertEqual(size_class(100_000).value, "100kb_1mb")
        self.assertEqual(size_class(5_000_000).value, "5mb_20mb")
        self.assertEqual(size_class(20_000_000).value, "ge_20mb")

    def test_detection_is_stratified_by_size(self) -> None:
        report = evaluate(
            truth_segments=[
                StateSegment("5", 0, 50_000, CopyNumberState.LOSS, 1.0),
                StateSegment("5", 10_000_000, 40_000_000, CopyNumberState.GAIN, 3.0),
            ],
            query_segments=[StateSegment("5", 10_000_000, 40_000_000, CopyNumberState.GAIN, 3.0)],
            evaluable={"5": [(0, CHR5)]},
            reference_bases=CHR5,
            truth_background=CopyNumberState.NEUTRAL,
            query_background=CopyNumberState.NEUTRAL,
        )
        by_label = {item.label: item for item in report.detection_by_size_class}
        self.assertEqual(by_label["lt_100kb"].detected, 0)
        self.assertEqual(by_label["lt_100kb"].missed, 1)
        self.assertEqual(by_label["ge_20mb"].detected, 1)


class IntervalEstimateTests(unittest.TestCase):
    def test_single_event_detection_carries_a_wide_interval(self) -> None:
        report = _evaluate(
            [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
            [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
        )
        estimate = report.detection_rate
        self.assertEqual(estimate.point, 1.0)
        # One event out of one must not be reported as certain.
        self.assertLess(estimate.lower or 1.0, 0.3)


class ValidationTests(unittest.TestCase):
    def test_overlapping_query_segments_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _evaluate(
                [StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0)],
                [
                    StateSegment("5", 0, 10_000_000, CopyNumberState.LOSS, 1.0),
                    StateSegment("5", 5_000_000, 15_000_000, CopyNumberState.GAIN, 3.0),
                ],
            )

    def test_derive_events_merges_only_contiguous_same_state_runs(self) -> None:
        events = derive_events(
            [
                StateSegment("5", 0, 100, CopyNumberState.LOSS, 1.0),
                StateSegment("5", 100, 200, CopyNumberState.LOSS, 1.0),
                StateSegment("5", 300, 400, CopyNumberState.LOSS, 1.0),
                StateSegment("5", 400, 500, CopyNumberState.GAIN, 3.0),
            ],
            prefix="truth",
        )
        spans = [(item.start, item.end, item.state.value) for item in events]
        self.assertEqual(spans, [(0, 200, "loss"), (300, 400, "loss"), (400, 500, "gain")])

    def test_event_copy_number_is_length_weighted(self) -> None:
        events = derive_events(
            [
                StateSegment("5", 0, 300, CopyNumberState.LOSS, 1.0),
                StateSegment("5", 300, 400, CopyNumberState.LOSS, 0.0),
            ],
            prefix="truth",
        )
        self.assertEqual(len(events), 1)
        self.assertAlmostEqual(events[0].copy_number, 0.75)


if __name__ == "__main__":
    unittest.main()
