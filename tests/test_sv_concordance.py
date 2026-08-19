from __future__ import annotations

import unittest

from ontseq_platform.models import EventType, Evidence, Locus
from ontseq_platform.sv_concordance import (
    SVCallerObservation,
    SVConcordancePolicy,
    SVConcordanceStatus,
    compare_sv_caller_observations,
)


def _observation(
    observation_id: str,
    caller: str,
    caller_version: str,
    event_type: EventType,
    primary_start: int,
    secondary_start: int | None = None,
    *,
    primary_chromosome: str = "chr1",
    secondary_chromosome: str = "chr2",
) -> SVCallerObservation:
    secondary = None
    if secondary_start is not None:
        secondary = Locus(
            chromosome=secondary_chromosome,
            start=secondary_start,
            end=secondary_start + 1,
        )
    return SVCallerObservation(
        observation_id=observation_id,
        caller=caller,
        caller_version=caller_version,
        source_event_id=f"{caller.upper()}-{observation_id}",
        event_type=event_type,
        primary=Locus(
            chromosome=primary_chromosome,
            start=primary_start,
            end=primary_start + 1,
        ),
        secondary=secondary,
        evidence=Evidence(
            caller=caller,
            caller_version=caller_version,
            support_reads=10,
            precise=True,
        ),
    )


def _policy(distance: int = 50) -> SVConcordancePolicy:
    return SVConcordancePolicy(
        maximum_breakpoint_distance_bp=distance,
        note="Synthetic software-comparison tolerance; not a clinical threshold.",
    )


class SVConcordanceTests(unittest.TestCase):
    def test_exact_multi_caller_match_is_evidence_not_truth(self) -> None:
        sniffles = _observation(
            "sniffles-1",
            "Sniffles2",
            "2.8.0",
            EventType.TRANSLOCATION,
            100_000,
            200_000,
        )
        cutesv = _observation(
            "cutesv-1",
            "cuteSV",
            "2.1.3",
            EventType.TRANSLOCATION,
            100_000,
            200_000,
        )

        report = compare_sv_caller_observations([sniffles], [cutesv], _policy())

        self.assertEqual(len(report.pairs), 1)
        pair = report.pairs[0]
        self.assertEqual(pair.status, SVConcordanceStatus.EXACT_MATCH)
        self.assertEqual(pair.maximum_breakpoint_distance_bp, 0)
        self.assertEqual(pair.evidence_semantics, "support_only_not_truth")
        self.assertTrue(pair.research_only)
        self.assertFalse(pair.reportable)
        self.assertEqual(report.conclusion_semantics, "caller_concordance_is_not_truth")

    def test_near_match_requires_explicit_tolerance(self) -> None:
        sniffles = _observation(
            "sniffles-1",
            "Sniffles2",
            "2.8.0",
            EventType.TRANSLOCATION,
            100_000,
            200_000,
        )
        cutesv = _observation(
            "cutesv-1",
            "cuteSV",
            "2.1.3",
            EventType.TRANSLOCATION,
            100_025,
            200_040,
        )

        report = compare_sv_caller_observations([sniffles], [cutesv], _policy(50))

        self.assertEqual(report.pairs[0].status, SVConcordanceStatus.NEAR_MATCH)
        self.assertEqual(report.pairs[0].maximum_breakpoint_distance_bp, 40)

    def test_outside_tolerance_remains_single_caller_evidence(self) -> None:
        sniffles = _observation(
            "sniffles-1",
            "Sniffles2",
            "2.8.0",
            EventType.TRANSLOCATION,
            100_000,
            200_000,
        )
        cutesv = _observation(
            "cutesv-1",
            "cuteSV",
            "2.1.3",
            EventType.TRANSLOCATION,
            100_100,
            200_100,
        )

        report = compare_sv_caller_observations([sniffles], [cutesv], _policy(50))

        self.assertEqual(report.pairs, [])
        self.assertEqual(report.unmatched_left_observation_ids, ["sniffles-1"])
        self.assertEqual(report.unmatched_right_observation_ids, ["cutesv-1"])

    def test_same_geometry_with_different_event_type_is_conflict(self) -> None:
        sniffles = _observation(
            "sniffles-1",
            "Sniffles2",
            "2.8.0",
            EventType.TRANSLOCATION,
            100_000,
            200_000,
        )
        cutesv = _observation(
            "cutesv-1",
            "cuteSV",
            "2.1.3",
            EventType.FUSION,
            100_000,
            200_000,
        )

        report = compare_sv_caller_observations([sniffles], [cutesv], _policy())

        self.assertEqual(report.pairs[0].status, SVConcordanceStatus.TOPOLOGY_CONFLICT)
        self.assertTrue(any("conflicting event types" in warning for warning in report.warnings))

    def test_reciprocal_breakpoint_order_is_order_invariant(self) -> None:
        sniffles = _observation(
            "sniffles-1",
            "Sniffles2",
            "2.8.0",
            EventType.TRANSLOCATION,
            100_000,
            200_000,
        )
        cutesv = _observation(
            "cutesv-1",
            "cuteSV",
            "2.1.3",
            EventType.TRANSLOCATION,
            200_000,
            100_000,
            primary_chromosome="chr2",
            secondary_chromosome="chr1",
        )

        report = compare_sv_caller_observations([sniffles], [cutesv], _policy())

        self.assertEqual(report.pairs[0].status, SVConcordanceStatus.EXACT_MATCH)
        self.assertTrue(report.pairs[0].breakpoint_order_swapped)

    def test_serialized_concordance_contains_no_raw_file_or_read_fields(self) -> None:
        sniffles = _observation(
            "sniffles-1",
            "Sniffles2",
            "2.8.0",
            EventType.TRANSLOCATION,
            100_000,
            200_000,
        )
        cutesv = _observation(
            "cutesv-1",
            "cuteSV",
            "2.1.3",
            EventType.TRANSLOCATION,
            100_000,
            200_000,
        )

        serialized = compare_sv_caller_observations(
            [sniffles],
            [cutesv],
            _policy(),
        ).model_dump_json()

        for forbidden in [
            "raw_alt",
            "read_names",
            "inserted_sequence",
            "source_path",
            ".vcf",
            ".bam",
        ]:
            self.assertNotIn(forbidden, serialized)


if __name__ == "__main__":
    unittest.main()
