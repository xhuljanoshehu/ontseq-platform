from __future__ import annotations

import unittest

from ontseq_platform.benchmark import benchmark_case
from ontseq_platform.models import (
    BenchmarkCase,
    BenchmarkKind,
    BenchmarkThresholds,
    EventType,
    GenomeBuild,
    GenomicEvent,
    Locus,
)


def _event(
    event_id: str,
    event_type: EventType,
    chromosome: str,
    start: int,
    end: int,
    *,
    secondary: Locus | None = None,
    copy_number: float | None = None,
) -> GenomicEvent:
    return GenomicEvent(
        event_id=event_id,
        event_type=event_type,
        primary=Locus(chromosome=chromosome, start=start, end=end),
        secondary=secondary,
        copy_number=copy_number,
    )


class BenchmarkTests(unittest.TestCase):
    def test_cnv_matching_is_one_to_one_and_reports_false_calls(self) -> None:
        case = BenchmarkCase(
            case_id="SYNTHETIC_CNV_001",
            kind=BenchmarkKind.CNV,
            genome_build=GenomeBuild.GRCH38,
            thresholds=BenchmarkThresholds(
                minimum_reciprocal_overlap=0.5,
                copy_number_tolerance=0.5,
            ),
            truth_events=[
                _event("truth-1", EventType.DELETION, "chr5", 100, 300, copy_number=1),
                _event("truth-2", EventType.DUPLICATION, "chr8", 500, 900, copy_number=3),
            ],
            query_events=[
                _event("query-1", EventType.DELETION, "5", 110, 290, copy_number=1),
                _event("query-fp", EventType.DUPLICATION, "chr2", 100, 200, copy_number=3),
            ],
        )
        report = benchmark_case(case)
        self.assertEqual(report.metrics.true_positive, 1)
        self.assertEqual(report.metrics.false_positive, 1)
        self.assertEqual(report.metrics.false_negative, 1)
        self.assertEqual(report.metrics.precision, 0.5)
        self.assertEqual(report.metrics.recall, 0.5)

    def test_translocation_matching_accepts_swapped_breakend_order(self) -> None:
        truth = _event(
            "truth-t",
            EventType.TRANSLOCATION,
            "chr8",
            1000,
            1001,
            secondary=Locus(chromosome="chr21", start=2000, end=2001),
        )
        query = _event(
            "query-t",
            EventType.TRANSLOCATION,
            "chr21",
            2050,
            2051,
            secondary=Locus(chromosome="chr8", start=1040, end=1041),
        )
        report = benchmark_case(
            BenchmarkCase(
                case_id="SYNTHETIC_SV_001",
                kind=BenchmarkKind.SV,
                genome_build=GenomeBuild.GRCH38,
                thresholds=BenchmarkThresholds(maximum_breakpoint_distance_bp=100),
                truth_events=[truth],
                query_events=[query],
            )
        )
        self.assertEqual(report.metrics.true_positive, 1)
        self.assertEqual(report.matches[0].maximum_breakpoint_distance_bp, 50)

    def test_matching_maximizes_true_positive_count_before_local_score(self) -> None:
        case = BenchmarkCase(
            case_id="MAX_CARDINALITY_CNV",
            kind=BenchmarkKind.CNV,
            genome_build=GenomeBuild.GRCH38,
            truth_events=[
                _event("T1", EventType.DELETION, "chr1", 0, 100),
                _event("T2", EventType.DELETION, "chr1", 80, 180),
            ],
            query_events=[
                _event("Q1", EventType.DELETION, "chr1", 50, 130),
                _event("Q2", EventType.DELETION, "chr1", 0, 50),
            ],
            thresholds=BenchmarkThresholds(minimum_reciprocal_overlap=0.5),
            research_only=True,
        )

        report = benchmark_case(case)

        self.assertEqual(report.metrics.true_positive, 2)
        self.assertEqual(report.metrics.false_positive, 0)
        self.assertEqual(report.metrics.false_negative, 0)

    def test_empty_negative_case_does_not_claim_undefined_precision(self) -> None:
        report = benchmark_case(
            BenchmarkCase(
                case_id="SYNTHETIC_NEGATIVE_001",
                kind=BenchmarkKind.CNV,
                genome_build=GenomeBuild.GRCH38,
                truth_events=[],
                query_events=[],
            )
        )
        self.assertIsNone(report.metrics.precision)
        self.assertIsNone(report.metrics.recall)
        self.assertIsNone(report.metrics.f1)


if __name__ == "__main__":
    unittest.main()
