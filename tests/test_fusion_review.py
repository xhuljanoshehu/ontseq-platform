from __future__ import annotations

import unittest

from ontseq_platform.fusion import (
    BreakpointGeneHit,
    FusionBreakpoint,
    FusionCandidate,
    FusionClassification,
    FusionGenePair,
    FusionInterpretationReport,
    ObservabilityStatus,
)
from ontseq_platform.fusion_review import (
    FusionReviewDisposition,
    build_fusion_reviewer_report,
)
from ontseq_platform.models import Evidence, GenomeBuild, ModuleRunStatus


def _candidate(event_id: str, candidate_id: str) -> FusionCandidate:
    primary = FusionBreakpoint(
        chromosome="chr1",
        position_0based=100_000,
        genes=[BreakpointGeneHit(gene="GENE1", strand="+", distance_bp=0)],
        observability=ObservabilityStatus.OBSERVABLE,
        observability_reason="synthetic observable region",
    )
    secondary = FusionBreakpoint(
        chromosome="chr2",
        position_0based=200_000,
        genes=[BreakpointGeneHit(gene="GENE2", strand="-", distance_bp=0)],
        observability=ObservabilityStatus.OBSERVABLE,
        observability_reason="synthetic observable region",
    )
    return FusionCandidate(
        candidate_id=candidate_id,
        source_event_id=event_id,
        primary=primary,
        secondary=secondary,
        classification=FusionClassification.GENE_GENE,
        gene_pairs=[
            FusionGenePair(
                gene_a="GENE1",
                gene_b="GENE2",
                known_pair=True,
                orientation_resolved=False,
            )
        ],
        evidence=[
            Evidence(
                caller="Sniffles2",
                caller_version="2.8.0",
                support_reads=12,
                variant_allele_fraction=0.25,
                precise=True,
            )
        ],
        limitations=["Synthetic reviewer-contract fixture only."],
    )


def _report(
    *candidates: FusionCandidate,
    status: ModuleRunStatus | None = None,
) -> FusionInterpretationReport:
    resolved_status = status or (
        ModuleRunStatus.COMPLETED if candidates else ModuleRunStatus.NO_CALL
    )
    return FusionInterpretationReport(
        sample_id="SYNTHETIC_REVIEW_001",
        genome_build=GenomeBuild.GRCH38,
        status=resolved_status,
        annotation_resource_id="synthetic-genes",
        annotation_resource_version="v1",
        annotation_source_sha256="0" * 64,
        candidates=list(candidates),
        source_translocation_count=len(candidates),
    )


class FusionReviewerContractTests(unittest.TestCase):
    def test_no_call_is_not_serialized_as_negative(self) -> None:
        reviewer = build_fusion_reviewer_report(_report())

        self.assertEqual(reviewer.disposition, FusionReviewDisposition.NO_CALL)
        self.assertEqual(reviewer.absence_interpretation, "not_established")
        self.assertEqual(reviewer.candidates, [])
        self.assertTrue(
            any("not a validated biological negative" in warning for warning in reviewer.warnings)
        )

    def test_candidate_remains_research_only_and_non_reportable(self) -> None:
        reviewer = build_fusion_reviewer_report(
            _report(_candidate("SNIFFLES2-000001", "FUSION-000001"))
        )

        self.assertEqual(reviewer.disposition, FusionReviewDisposition.CANDIDATE_EVIDENCE)
        self.assertEqual(len(reviewer.candidates), 1)
        item = reviewer.candidates[0]
        self.assertTrue(item.research_only)
        self.assertFalse(item.reportable)
        self.assertTrue(item.known_pair_present)
        self.assertTrue(item.both_breakpoints_observable)
        self.assertFalse(item.transcript_orientation_resolved)
        self.assertFalse(item.orientation_evidence_available)
        self.assertEqual(item.caller_evidence[0].caller, "Sniffles2")
        self.assertEqual(item.caller_evidence[0].support_reads, 12)

    def test_exact_duplicate_breakpoints_are_flagged_but_not_collapsed(self) -> None:
        reviewer = build_fusion_reviewer_report(
            _report(
                _candidate("SNIFFLES2-000001", "FUSION-000001"),
                _candidate("SNIFFLES2-000002", "FUSION-000002"),
            )
        )

        self.assertEqual(len(reviewer.candidates), 2)
        self.assertEqual(reviewer.redundancy_group_count, 1)
        for item in reviewer.candidates:
            self.assertTrue(item.potentially_redundant)
            self.assertEqual(
                set(item.redundancy_group_event_ids),
                {"SNIFFLES2-000001", "SNIFFLES2-000002"},
            )

    def test_reviewer_json_contains_no_raw_vcf_or_read_level_fields(self) -> None:
        reviewer = build_fusion_reviewer_report(
            _report(_candidate("SNIFFLES2-000001", "FUSION-000001"))
        )
        serialized = reviewer.model_dump_json()

        forbidden = [
            "raw_alt",
            "inserted_sequence",
            "read_names",
            "source_path",
            ".vcf",
            ".bam",
        ]
        for token in forbidden:
            self.assertNotIn(token, serialized)


if __name__ == "__main__":
    unittest.main()
