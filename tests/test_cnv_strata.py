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


class InferentialClaimTests(unittest.TestCase):
    """`favours` is an inferential claim; the direction of the counts is not.

    A four-nil split reads as an obvious winner and cannot reach any conventional
    threshold, because the smallest attainable two-sided exact p-value at four discordant
    pairs is 0.125. Naming a winner there asserts something the data cannot support at any
    alpha, which is why the two are separate fields.
    """

    def _pair(self, a_detected: bool, b_detected: bool, sample: str):
        return (
            _report(tumor_fraction=1.0, detected=a_detected, method="method-a", sample=sample),
            _report(tumor_fraction=1.0, detected=b_detected, method="method-b", sample=sample),
        )

    def _compare(self, pairs, **kwargs):
        return paired_detection_comparison([a for a, _ in pairs], [b for _, b in pairs], **kwargs)

    def test_a_clean_sweep_too_small_to_test_favours_nobody(self) -> None:
        result = self._compare([self._pair(True, False, f"S_{i}") for i in range(4)])
        self.assertEqual(result.only_a_detected, 4)
        self.assertEqual(result.favours, "neither")
        self.assertTrue(result.underpowered)
        self.assertAlmostEqual(result.minimum_attainable_p_value, 0.125)

    def test_the_observed_direction_is_still_reported(self) -> None:
        """Withholding the claim must not mean withholding the description."""
        result = self._compare([self._pair(True, False, f"S_{i}") for i in range(4)])
        self.assertEqual(result.observed_direction, "a")
        self.assertIn("description and not", result.note)

    def test_the_note_says_no_result_was_attainable(self) -> None:
        result = self._compare([self._pair(True, False, f"S_{i}") for i in range(4)])
        self.assertIn("underpowered by design", result.note)

    def test_a_significant_result_does_name_a_method(self) -> None:
        result = self._compare([self._pair(True, False, f"S_{i}") for i in range(6)])
        self.assertFalse(result.underpowered)
        self.assertEqual(result.favours, "a")
        self.assertEqual(result.observed_direction, "a")

    def test_a_non_significant_but_powered_comparison_favours_neither(self) -> None:
        pairs = [self._pair(True, False, f"S_{i}") for i in range(4)]
        pairs += [self._pair(False, True, f"T_{i}") for i in range(3)]
        result = self._compare(pairs)
        self.assertFalse(result.underpowered)
        self.assertEqual(result.favours, "neither")
        self.assertEqual(result.observed_direction, "a")
        self.assertIn("not evidence of equivalence", result.note)

    def test_alpha_is_recorded_so_the_claim_can_be_audited(self) -> None:
        result = self._compare([self._pair(True, False, f"S_{i}") for i in range(6)])
        self.assertAlmostEqual(result.alpha, 0.05)

    def test_a_stricter_alpha_can_withdraw_the_claim(self) -> None:
        """The direction must follow the pre-specified rule, not the other way round."""
        pairs = [self._pair(True, False, f"S_{i}") for i in range(6)]
        self.assertEqual(self._compare(pairs, alpha=0.05).favours, "a")
        self.assertEqual(self._compare(pairs, alpha=0.01).favours, "neither")

    def test_an_impossible_alpha_is_refused(self) -> None:
        pairs = [self._pair(True, False, "S_1")]
        with self.assertRaises(ValueError):
            self._compare(pairs, alpha=0.0)


class SpecimenClusteringTests(unittest.TestCase):
    """Events inside one specimen are not independent observations.

    They share its purity, library, coverage and artefacts. Every interval here is
    computed over events, so where events cluster the intervals are narrower than the data
    support. The module cannot fix that on the caller's behalf — a specimen-level endpoint
    is a study-design decision — but it must not present the problem as absent.
    """

    def _reports(self, samples):
        return [_report(tumor_fraction=1.0, detected=True, sample=sample) for sample in samples]

    def test_one_event_per_specimen_is_not_flagged(self) -> None:
        report = aggregate(self._reports(["S_1", "S_2", "S_3"]), aggregate_id="AGG_001")
        self.assertEqual(report.clustering.specimens, 3)
        self.assertFalse(report.clustering.intervals_are_anticonservative)

    def test_repeated_specimens_are_flagged(self) -> None:
        report = aggregate(self._reports(["S_1", "S_1", "S_2"]), aggregate_id="AGG_001")
        self.assertEqual(report.clustering.specimens, 2)
        self.assertEqual(report.clustering.largest_specimen_events, 2)
        self.assertTrue(report.clustering.intervals_are_anticonservative)

    def test_the_clustering_is_stated_in_words_not_only_in_a_flag(self) -> None:
        report = aggregate(self._reports(["S_1", "S_1", "S_2"]), aggregate_id="AGG_001")
        self.assertTrue(
            any("not independent" in warning for warning in report.warnings), report.warnings
        )

    def test_a_specimen_level_rate_is_reported_beside_the_event_level_one(self) -> None:
        report = aggregate(self._reports(["S_1", "S_1", "S_2"]), aggregate_id="AGG_001")
        self.assertIsNotNone(report.specimen_level_detection_rate)
        self.assertEqual(report.specimen_level_detection_rate.total, 2)
        self.assertEqual(report.overall_detection_rate.total, 3)

    def test_the_specimen_rate_does_not_replace_the_event_rate(self) -> None:
        """They answer different questions; a reader must be able to see both."""
        report = aggregate(self._reports(["S_1", "S_1", "S_2"]), aggregate_id="AGG_001")
        self.assertIsNotNone(report.overall_detection_rate)

    def test_the_limitation_is_recorded_in_the_report(self) -> None:
        report = aggregate(self._reports(["S_1", "S_2"]), aggregate_id="AGG_001")
        self.assertTrue(
            any("event-level" in item for item in report.limitations), report.limitations
        )


class ClusteredPairedComparisonTests(unittest.TestCase):
    """McNemar treats each discordant pair as an independent coin flip."""

    def _pair(self, a_detected: bool, b_detected: bool, sample: str):
        return (
            _report(tumor_fraction=1.0, detected=a_detected, method="method-a", sample=sample),
            _report(tumor_fraction=1.0, detected=b_detected, method="method-b", sample=sample),
        )

    def test_one_specimen_per_pair_is_not_flagged(self) -> None:
        pairs = [self._pair(True, False, f"S_{i}") for i in range(6)]
        result = paired_detection_comparison([a for a, _ in pairs], [b for _, b in pairs])
        self.assertEqual(result.discordant_specimens, 6)
        self.assertFalse(result.p_value_is_anticonservative)

    def test_no_discordant_pair_leaves_the_count_at_zero(self) -> None:
        pairs = [self._pair(True, True, f"S_{i}") for i in range(3)]
        result = paired_detection_comparison([a for a, _ in pairs], [b for _, b in pairs])
        self.assertEqual(result.discordant_specimens, 0)
        self.assertFalse(result.p_value_is_anticonservative)
