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
from ontseq_platform.cnv.strata import aggregate, compare_aggregates, estimate_limit_of_detection
from ontseq_platform.models import GenomeBuild, ModuleRunStatus

CONTIGS = {"chr5": 181_538_259}
LOSS = CnvSegment(
    contig="chr5", start=70_000_000, end=160_000_000, state=CopyNumberState.LOSS, copy_number=1.0
)


def _report(*, tumor_fraction: float, detected: bool, method: str = "baseline", sample: str):
    truth = CnvTruthSet(
        truth_id="T1",
        sample_id=sample,
        genome_build=GenomeBuild.GRCH38,
        source=CnvTruthSource.SIMULATED,
        source_version="v1",
        background_state=CopyNumberState.NEUTRAL,
        resolution_bp=100_000,
        segments=[LOSS],
    )
    segments = [LOSS] if detected else [
        CnvSegment(
            contig="chr5", start=0, end=1_000_000, state=CopyNumberState.GAIN, copy_number=3.0
        )
    ]
    call_set = CnvCallSet(
        call_set_id="C1",
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
        case_id="CASE",
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
            _report(tumor_fraction=1.0, detected=True, method="a", sample="S1"),
            _report(tumor_fraction=1.0, detected=True, method="b", sample="S2"),
        ]
        with self.assertRaises(ValueError):
            aggregate(reports, aggregate_id="AGG")

    def test_empty_input_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            aggregate([], aggregate_id="AGG")

    def test_strata_are_grouped_by_tumor_fraction(self) -> None:
        reports = [
            _report(tumor_fraction=1.0, detected=True, sample="S1"),
            _report(tumor_fraction=1.0, detected=True, sample="S2"),
            _report(tumor_fraction=0.1, detected=False, sample="S3"),
        ]
        result = aggregate(reports, aggregate_id="AGG")
        by_label = {item.label: item for item in result.by_tumor_fraction}
        self.assertEqual(by_label["1"].detected_events, 2)
        self.assertEqual(by_label["0.1"].detected_events, 0)
        self.assertEqual(by_label["0.1"].missed_events, 1)
        self.assertEqual(result.overall_detection_rate.successes, 2)
        self.assertEqual(result.overall_detection_rate.total, 3)

    def test_size_class_strata_are_pooled(self) -> None:
        reports = [_report(tumor_fraction=1.0, detected=True, sample="S1")]
        result = aggregate(reports, aggregate_id="AGG")
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
        limit = estimate_limit_of_detection(
            [(0.1, 0, 0)], predictor="tumor_fraction"
        )
        self.assertEqual(limit.levels_used, 0)
        self.assertIsNone(limit.empirical_value)
        self.assertIsNone(limit.model_based_value)

    def test_aggregate_emits_a_limit_when_levels_permit(self) -> None:
        reports = [
            _report(tumor_fraction=1.0, detected=True, sample=f"S{i}") for i in range(3)
        ] + [
            _report(tumor_fraction=0.05, detected=False, sample=f"T{i}") for i in range(3)
        ]
        result = aggregate(reports, aggregate_id="AGG")
        predictors = {item.predictor for item in result.limits_of_detection}
        self.assertIn("tumor_fraction", predictors)


class ComparisonTests(unittest.TestCase):
    def test_comparison_returns_intervals_not_a_verdict(self) -> None:
        first = aggregate(
            [_report(tumor_fraction=1.0, detected=True, method="a", sample="S1")],
            aggregate_id="A",
        )
        second = aggregate(
            [_report(tumor_fraction=1.0, detected=False, method="b", sample="S2")],
            aggregate_id="B",
        )
        rows = compare_aggregates([first, second])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0][0], "a")
        # Lower bounds are surfaced so an overlapping comparison stays visible.
        self.assertIsNotNone(rows[0][3])


if __name__ == "__main__":
    unittest.main()
