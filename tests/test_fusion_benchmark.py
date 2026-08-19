from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.breakends import (
    BreakendAltForm,
    BreakendDescriptor,
    SnifflesJunctionOrientation,
)
from ontseq_platform.fusion import (
    BreakpointGeneHit,
    FusionBreakpoint,
    FusionCandidate,
    FusionClassification,
    FusionGenePair,
    FusionInterpretationReport,
    ObservabilityStatus,
)
from ontseq_platform.fusion_benchmark import (
    PartnerExpectation,
    SyntheticFusionBenchmarkSuite,
    SyntheticFusionFixture,
    evaluate_fusion_benchmark_case,
    load_synthetic_fusion_benchmark,
)
from ontseq_platform.models import GenomeBuild, ModuleRunStatus


def _candidate(
    *,
    genes: tuple[str, str] = ("BCR", "ABL1"),
    classification: FusionClassification = FusionClassification.GENE_GENE,
    primary_observability: ObservabilityStatus = ObservabilityStatus.OBSERVABLE,
    secondary_observability: ObservabilityStatus = ObservabilityStatus.OBSERVABLE,
    with_descriptor: bool = True,
) -> FusionCandidate:
    primary = FusionBreakpoint(
        chromosome="chr1",
        position_0based=100,
        genes=[BreakpointGeneHit(gene=genes[0], strand="+", distance_bp=0)],
        observability=primary_observability,
    )
    secondary = FusionBreakpoint(
        chromosome="chr2",
        position_0based=200,
        genes=[BreakpointGeneHit(gene=genes[1], strand="-", distance_bp=0)],
        observability=secondary_observability,
    )
    descriptor = None
    if with_descriptor:
        descriptor = BreakendDescriptor(
            source_event_id="SNIFFLES2-000001",
            primary_chromosome="chr1",
            primary_position_0based=100,
            mate_chromosome="chr2",
            mate_position_0based=200,
            alt_form=BreakendAltForm.LOCAL_THEN_CLOSE,
            sniffles_junction_orientation=SnifflesJunctionOrientation.PLUS_PLUS,
        )
    return FusionCandidate(
        candidate_id="FUSION-SNIFFLES2-000001",
        source_event_id="SNIFFLES2-000001",
        primary=primary,
        secondary=secondary,
        classification=classification,
        breakend_descriptor=descriptor,
        gene_pairs=[
            FusionGenePair(
                gene_a=genes[0],
                gene_b=genes[1],
                orientation_resolved=False,
            )
        ],
        evidence=[],
        limitations=["synthetic benchmark candidate"],
    )


def _report(candidate: FusionCandidate | None) -> FusionInterpretationReport:
    candidates = [candidate] if candidate is not None else []
    return FusionInterpretationReport(
        sample_id="SYNTHETIC_BENCHMARK",
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.COMPLETED if candidates else ModuleRunStatus.NO_CALL,
        annotation_resource_id="synthetic",
        annotation_resource_version="v1",
        annotation_source_sha256="0" * 64,
        candidates=candidates,
        source_translocation_count=len(candidates),
        breakend_descriptor_count=sum(item.breakend_descriptor is not None for item in candidates),
        warnings=["synthetic benchmark only"],
    )


class FusionBenchmarkContractTests(unittest.TestCase):
    def test_repository_suite_is_explicitly_nonbiological(self) -> None:
        suite = load_synthetic_fusion_benchmark(Path("configs/sv/fusion.synthetic_benchmark.yaml"))

        self.assertEqual(suite.status, "synthetic_software_benchmark_only")
        self.assertGreaterEqual(len(suite.fixtures), 9)
        self.assertTrue(all(not fixture.clinical_truth for fixture in suite.fixtures))
        self.assertTrue(
            all(
                fixture.coordinate_semantics == "synthetic_nonbiological"
                for fixture in suite.fixtures
            )
        )
        labels = {fixture.family_label for fixture in suite.fixtures}
        self.assertTrue(any("BCR::ABL1" in label for label in labels))
        self.assertTrue(any("PML::RARA" in label for label in labels))
        self.assertTrue(any("RUNX1::RUNX1T1" in label for label in labels))
        self.assertTrue(any("CBFB::MYH11" in label for label in labels))
        self.assertTrue(any("KMT2A::*" in label for label in labels))

    def test_exact_pair_fixture_passes_expected_candidate(self) -> None:
        fixture = SyntheticFusionFixture(
            fixture_id="exact_pair",
            family_label="BCR::ABL1 synthetic",
            expected_classification=FusionClassification.GENE_GENE,
            partner_expectation=PartnerExpectation.EXACT_PAIR,
            expected_gene_pair=("BCR", "ABL1"),
            expected_primary_observability=ObservabilityStatus.OBSERVABLE,
            expected_secondary_observability=ObservabilityStatus.OBSERVABLE,
            accepted_junction_orientations=[SnifflesJunctionOrientation.PLUS_PLUS],
            note="software behavior only",
        )
        result = evaluate_fusion_benchmark_case(fixture, _report(_candidate()))

        self.assertTrue(result.passed)
        self.assertEqual(result.failures, [])

    def test_wrong_pair_and_observability_fail_with_reasons(self) -> None:
        fixture = SyntheticFusionFixture(
            fixture_id="wrong_pair_guard",
            family_label="BCR::ABL1 synthetic",
            expected_classification=FusionClassification.GENE_GENE,
            partner_expectation=PartnerExpectation.EXACT_PAIR,
            expected_gene_pair=("BCR", "ABL1"),
            expected_primary_observability=ObservabilityStatus.OBSERVABLE,
            expected_secondary_observability=ObservabilityStatus.OBSERVABLE,
            note="software behavior only",
        )
        candidate = _candidate(
            genes=("BCR", "OTHER"),
            secondary_observability=ObservabilityStatus.UNKNOWN,
        )
        result = evaluate_fusion_benchmark_case(fixture, _report(candidate))

        self.assertFalse(result.passed)
        self.assertTrue(any("expected gene pair not observed" in item for item in result.failures))
        self.assertTrue(any("secondary observability mismatch" in item for item in result.failures))

    def test_any_partner_fixture_accepts_anchor_gene_without_whitelist(self) -> None:
        fixture = SyntheticFusionFixture(
            fixture_id="kmt2a_partner",
            family_label="KMT2A::* synthetic",
            expected_classification=FusionClassification.GENE_GENE,
            partner_expectation=PartnerExpectation.ANY_PARTNER,
            anchor_gene="KMT2A",
            note="software behavior only",
        )
        result = evaluate_fusion_benchmark_case(
            fixture,
            _report(_candidate(genes=("KMT2A", "SYNTHPARTNER"))),
        )

        self.assertTrue(result.passed)

    def test_no_candidate_control_passes_no_call_report(self) -> None:
        fixture = SyntheticFusionFixture(
            fixture_id="no_candidate",
            family_label="synthetic no-source-event control",
            expected_classification=FusionClassification.UNRESOLVED,
            partner_expectation=PartnerExpectation.NONE,
            expected_candidate_count=0,
            note="software behavior only",
        )
        result = evaluate_fusion_benchmark_case(fixture, _report(None))

        self.assertTrue(result.passed)

    def test_duplicate_fixture_ids_are_rejected(self) -> None:
        fixture = SyntheticFusionFixture(
            fixture_id="duplicate",
            family_label="synthetic control",
            expected_classification=FusionClassification.UNRESOLVED,
            expected_candidate_count=0,
            note="software behavior only",
        )
        with self.assertRaisesRegex(ValueError, "fixture ids must be unique"):
            SyntheticFusionBenchmarkSuite(
                suite_id="duplicate_suite",
                description="synthetic only",
                fixtures=[fixture, fixture],
            )

    def test_invalid_partner_contract_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "exact-pair fixtures require"):
            SyntheticFusionFixture(
                fixture_id="invalid_pair",
                family_label="invalid synthetic",
                expected_classification=FusionClassification.GENE_GENE,
                partner_expectation=PartnerExpectation.EXACT_PAIR,
                note="software behavior only",
            )

    def test_loader_rejects_non_mapping_yaml(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.yaml"
            path.write_text("- not\n- a\n- mapping\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must contain a mapping"):
                load_synthetic_fusion_benchmark(path)


if __name__ == "__main__":
    unittest.main()
