from __future__ import annotations

import unittest

from ontseq_platform.cnv.evaluate import evaluate_case
from ontseq_platform.cnv.models import (
    CnvBenchmarkCase,
    CnvCallSet,
    CnvDataBasis,
    CnvSegment,
    CnvStrata,
    CnvTruthSet,
    CnvTruthSource,
)
from ontseq_platform.cnv.states import CopyNumberState
from ontseq_platform.cnv.strata import (
    aggregate,
    compare_aggregates,
    estimate_limit_of_detection,
    paired_detection_comparison,
)
from ontseq_platform.models import GenomeBuild, ModuleRunStatus

CONTIGS = {"chr5": 181_538_259}
LOSS = CnvSegment(
    contig="chr5", start=70_000_000, end=160_000_000, state=CopyNumberState.LOSS, copy_number=1.0
)


def _report(*, tumor_fraction: float, detected: bool, method: str = "baseline", sample: str):
    truth = CnvTruthSet(
        truth_id="TRUTH_001",
        sample_id=sample,
        genome_build=GenomeBuild.GRCH38,
        source=CnvTruthSource.SIMULATED,
        source_version="v1",
        background_state=CopyNumberState.NEUTRAL,
        resolution_bp=100_000,
        segments=[LOSS],
    )
    segments = (
        [LOSS]
        if detected
        else [
            CnvSegment(
                contig="chr5", start=0, end=1_000_000, state=CopyNumberState.GAIN, copy_number=3.0
            )
        ]
    )
    call_set = CnvCallSet(
        call_set_id="CALLSET_001",
        sample_id=sample,
        genome_build=GenomeBuild.GRCH38,
        method=method,
        method_version="0.1.0",
        data_basis=CnvDataBasis.SIMULATED,
        background_state=CopyNumberState.NEUTRAL,
        status=ModuleRunStatus.COMPLETED,
        segments=segments,
    )
    case = CnvBenchmarkCase(
        case_id="CASE_001",
        genome_build=GenomeBuild.GRCH38,
        contig_lengths=CONTIGS,
        truth=truth,
        call_set=call_set,
        strata=CnvStrata(tumor_fraction=tumor_fraction, mean_coverage_x=3.0),
    )
    return evaluate_case(case)


class AggregationTests(unittest.TestCase):
    def test_pooling_two_methods_is_refused(self) -> None:
        reports = [
            _report(tumor_fraction=1.0, detected=True, method="method-a", sample="SAMPLE_001"),
            _report(tumor_fraction=1.0, detected=True, method="method-b", sample="SAMPLE_002"),
        ]
        with self.assertRaises(ValueError):
            aggregate(reports, aggregate_id="AGGREGATE_001")

    def test_empty_input_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            aggregate([], aggregate_id="AGGREGATE_001")

    def test_strata_are_grouped_by_tumor_fraction(self) -> None:
        reports = [
            _report(tumor_fraction=1.0, detected=True, sample="SAMPLE_001"),
            _report(tumor_fraction=1.0, detected=True, sample="SAMPLE_002"),
            _report(tumor_fraction=0.1, detected=False, sample="SAMPLE_003"),
        ]
        result = aggregate(reports, aggregate_id="AGGREGATE_001")
        by_label = {item.label: item for item in result.by_tumor_fraction}
        self.assertEqual(by_label["1"].detected_events, 2)
        self.assertEqual(by_label["0.1"].detected_events, 0)
        self.assertEqual(by_label["0.1"].missed_events, 1)
        self.assertEqual(result.overall_detection_rate.successes, 2)
        self.assertEqual(result.overall_detection_rate.total, 3)

    def test_size_class_strata_are_pooled(self) -> None:
        reports = [_report(tumor_fraction=1.0, detected=True, sample="SAMPLE_001")]
        result = aggregate(reports, aggregate_id="AGGREGATE_001")
        labels = {item.label for item in result.by_size_class}
        self.assertIn("ge_20mb", labels)


class LimitOfDetectionTests(unittest.TestCase):
    def test_graded_series_yields_both_estimates(self) -> None:
        limit = estimate_limit_of_detection(
            [(0.05, 1, 20), (0.1, 8, 20), (0.2, 17, 20), (0.5, 20, 20)],
            predictor="tumor_fraction",
        )
        self.assertEqual(limit.levels_used, 4)
        self.assertTrue(limit.model_converged)
        self.assertIsNotNone(limit.model_based_value)

    def test_empirical_value_requires_the_lower_bound_to_meet_the_target(self) -> None:
        # 20/20 gives a lower Wilson bound around 0.83, which does not reach 0.95.
        limit = estimate_limit_of_detection(
            [(0.5, 20, 20)], predictor="tumor_fraction", target_detection_rate=0.95
        )
        self.assertIsNone(limit.empirical_value)

    def test_large_sample_can_establish_an_empirical_limit(self) -> None:
        limit = estimate_limit_of_detection(
            [(0.05, 10, 200), (0.5, 200, 200)],
            predictor="tumor_fraction",
            target_detection_rate=0.95,
        )
        self.assertEqual(limit.empirical_value, 0.5)

    def test_separated_design_withholds_the_model_value(self) -> None:
        limit = estimate_limit_of_detection(
            [(0.05, 0, 10), (0.5, 10, 10)], predictor="tumor_fraction"
        )
        self.assertIsNone(limit.model_based_value)
        self.assertFalse(limit.model_converged)
        self.assertIn("separated", limit.note)

    def test_no_assessable_events_yields_no_limit(self) -> None:
        limit = estimate_limit_of_detection([(0.1, 0, 0)], predictor="tumor_fraction")
        self.assertEqual(limit.levels_used, 0)
        self.assertIsNone(limit.empirical_value)
        self.assertIsNone(limit.model_based_value)

    def test_aggregate_emits_a_limit_when_levels_permit(self) -> None:
        detected = [
            _report(tumor_fraction=1.0, detected=True, sample=f"SAMPLE_{i}") for i in range(3)
        ]
        missed = [
            _report(tumor_fraction=0.05, detected=False, sample=f"DILUTED_{i}") for i in range(3)
        ]
        reports = detected + missed
        result = aggregate(reports, aggregate_id="AGGREGATE_001")
        predictors = {item.predictor for item in result.limits_of_detection}
        self.assertIn("tumor_fraction", predictors)


class PairedComparisonTests(unittest.TestCase):
    """Method selection needs a paired test, not two independent rates."""

    def _pair(self, a_detected: bool, b_detected: bool, sample: str):
        return (
            _report(tumor_fraction=1.0, detected=a_detected, method="method-a", sample=sample),
            _report(tumor_fraction=1.0, detected=b_detected, method="method-b", sample=sample),
        )

    def test_full_agreement_yields_no_p_value(self) -> None:
        pairs = [self._pair(True, True, f"SAMPLE_{i}") for i in range(4)]
        result = paired_detection_comparison([a for a, _ in pairs], [b for _, b in pairs])
        self.assertEqual(result.paired_events, 4)
        self.assertEqual(result.both_detected, 4)
        self.assertIsNone(result.p_value)
        self.assertEqual(result.favours, "neither")
        self.assertIn("not evidence of equivalence", result.note)

    def test_consistent_advantage_is_detected(self) -> None:
        pairs = [self._pair(True, False, f"SAMPLE_{i}") for i in range(6)]
        result = paired_detection_comparison([a for a, _ in pairs], [b for _, b in pairs])
        self.assertEqual(result.only_a_detected, 6)
        self.assertEqual(result.only_b_detected, 0)
        self.assertEqual(result.favours, "a")
        self.assertAlmostEqual(result.p_value, 0.03125)

    def test_balanced_discordance_favours_neither(self) -> None:
        pairs = [self._pair(True, False, f"SAMPLE_{i}") for i in range(3)]
        pairs += [self._pair(False, True, f"OTHER_{i}") for i in range(3)]
        result = paired_detection_comparison([a for a, _ in pairs], [b for _, b in pairs])
        self.assertEqual(result.favours, "neither")
        self.assertAlmostEqual(result.p_value, 1.0)

    def test_events_not_assessable_under_both_are_excluded(self) -> None:
        shared = self._pair(True, False, "SAMPLE_001")
        lone = _report(tumor_fraction=1.0, detected=True, method="method-a", sample="SAMPLE_999")
        result = paired_detection_comparison([shared[0], lone], [shared[1]])
        self.assertEqual(result.paired_events, 1)
        self.assertGreater(result.unpaired_events, 0)

    def test_disjoint_samples_are_reported_as_no_comparison(self) -> None:
        a = _report(tumor_fraction=1.0, detected=True, method="method-a", sample="SAMPLE_001")
        b = _report(tumor_fraction=1.0, detected=True, method="method-b", sample="SAMPLE_777")
        result = paired_detection_comparison([a], [b])
        self.assertEqual(result.paired_events, 0)
        self.assertIsNone(result.p_value)
        self.assertIn("not a tie", result.note)

    def test_mixing_methods_on_one_side_is_rejected(self) -> None:
        a = _report(tumor_fraction=1.0, detected=True, method="method-a", sample="SAMPLE_001")
        b = _report(tumor_fraction=1.0, detected=True, method="method-b", sample="SAMPLE_001")
        with self.assertRaises(ValueError):
            paired_detection_comparison([a, b], [b])

    def test_empty_side_is_rejected(self) -> None:
        a = _report(tumor_fraction=1.0, detected=True, method="method-a", sample="SAMPLE_001")
        with self.assertRaises(ValueError):
            paired_detection_comparison([a], [])


class ComparisonTests(unittest.TestCase):
    def test_comparison_returns_intervals_not_a_verdict(self) -> None:
        first = aggregate(
            [_report(tumor_fraction=1.0, detected=True, method="method-a", sample="SAMPLE_001")],
            aggregate_id="AGGREGATE_A",
        )
        second = aggregate(
            [_report(tumor_fraction=1.0, detected=False, method="method-b", sample="SAMPLE_002")],
            aggregate_id="AGGREGATE_B",
        )
        rows = compare_aggregates([first, second])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0][0], "method-a")
        # Lower bounds are surfaced so an overlapping comparison stays visible.
        self.assertIsNotNone(rows[0][3])


if __name__ == "__main__":
    unittest.main()
