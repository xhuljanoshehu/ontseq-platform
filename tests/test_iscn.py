from __future__ import annotations

import unittest

from ontseq_platform.iscn import build_iscn_proposal, iscn_module_outcome
from ontseq_platform.models import (
    EventType,
    GenomeBuild,
    GenomicEvent,
    ISCNDispositionOutcome,
    ISCNProposal,
    ISCNProposalStatus,
    ISCNResourceProvenance,
    ISCNSelectionPolicy,
    Locus,
    ModuleRunStatus,
    ReferenceDictionaryContract,
    SvObservability,
)


def _provenance(build: GenomeBuild) -> ISCNResourceProvenance:
    return ISCNResourceProvenance(
        genome_build=build,
        reference_dictionary_contract=ReferenceDictionaryContract.EXACT_FULL,
        reference_bundle_id=f"SYNTHETIC_{build.value}",
        reference_bundle_version="test-v1",
        reference_lock_sha256="a" * 64,
        annotation_cache_sha256="c" * 64,
        cytoband_release=f"synthetic-{build.value}",
        cytoband_sha256="b" * 64,
    )


class ISCNProposalTests(unittest.TestCase):
    def test_subset_renderer_is_deterministic(self) -> None:
        events = [
            GenomicEvent(
                event_id="gain8",
                event_type=EventType.CHROMOSOME_GAIN,
                primary=Locus(chromosome="chr8", start=0, end=1),
                copy_number=3.0,
                whole_chromosome_span_confirmed=True,
                reportable=True,
            ),
            GenomicEvent(
                event_id="del5q",
                event_type=EventType.DELETION,
                primary=Locus(
                    chromosome="chr5",
                    start=10,
                    end=20,
                    cytoband_start="q13",
                    cytoband_end="q34",
                ),
                copy_number=1.0,
                reportable=True,
            ),
        ]
        proposal = build_iscn_proposal(
            events,
            genome_build=GenomeBuild.GRCH38,
            resource_provenance=_provenance(GenomeBuild.GRCH38),
            cnv_source_event_ids={"del5q", "gain8"},
        )
        self.assertEqual(proposal.notation, "del(5)(q13q34),+8")
        self.assertEqual(proposal.proposal_status, ISCNProposalStatus.PARTIAL_EVENT_LEVEL)
        self.assertEqual(proposal.source_event_ids, ["del5q", "gain8"])
        self.assertNotIn("46,", proposal.notation)

    def test_unsupported_event_never_becomes_a_normal_karyotype(self) -> None:
        event = GenomicEvent(
            event_id="ins1",
            event_type=EventType.INSERTION,
            primary=Locus(chromosome="chr1", start=10, end=11),
            reportable=True,
        )
        proposal = build_iscn_proposal(
            [event],
            genome_build=GenomeBuild.GRCH38,
            resource_provenance=_provenance(GenomeBuild.GRCH38),
        )
        self.assertIsNone(proposal.notation)
        self.assertEqual(proposal.proposal_status, ISCNProposalStatus.NO_RENDERABLE_CANDIDATE)
        self.assertEqual(proposal.event_dispositions[0].event_id, "ins1")
        self.assertTrue(any("not evidence of a normal" in warning for warning in proposal.warnings))

    def test_empty_assessed_event_set_never_becomes_a_normal_karyotype(self) -> None:
        proposal = build_iscn_proposal(
            [],
            genome_build=GenomeBuild.GRCH38,
            resource_provenance=_provenance(GenomeBuild.GRCH38),
        )

        self.assertIsNone(proposal.notation)
        self.assertFalse(proposal.baseline_assessed)
        self.assertEqual(proposal.proposal_status, ISCNProposalStatus.NO_RENDERABLE_CANDIDATE)
        self.assertTrue(any("normal karyotype" in warning for warning in proposal.warnings))

    def test_threshold_based_chromosome_call_without_exact_span_is_not_rendered(self) -> None:
        cases = (
            ("ninety-percent-gain", EventType.CHROMOSOME_GAIN, "chr8", 3.0),
            ("ninety-percent-loss", EventType.CHROMOSOME_LOSS, "chr7", 1.0),
        )
        for event_id, event_type, chromosome, copy_number in cases:
            with self.subTest(event_type=event_type):
                event = GenomicEvent(
                    event_id=event_id,
                    event_type=event_type,
                    primary=Locus(chromosome=chromosome, start=0, end=90),
                    copy_number=copy_number,
                    reportable=True,
                )
                proposal = build_iscn_proposal(
                    [event],
                    genome_build=GenomeBuild.GRCH38,
                    resource_provenance=_provenance(GenomeBuild.GRCH38),
                    cnv_source_event_ids={event.event_id},
                )

                self.assertIsNone(proposal.notation)
                self.assertEqual(
                    proposal.proposal_status,
                    ISCNProposalStatus.NO_RENDERABLE_CANDIDATE,
                )
                self.assertEqual(
                    proposal.event_dispositions[0].outcome,
                    ISCNDispositionOutcome.OMITTED_UNSUPPORTED,
                )
                self.assertIn("span was not confirmed", proposal.event_dispositions[0].detail)

    def test_exact_full_chromosome_gain_and_loss_render_symmetrically(self) -> None:
        cases = (
            ("exact-gain", EventType.CHROMOSOME_GAIN, "chr8", 3.0, "+8"),
            ("exact-loss", EventType.CHROMOSOME_LOSS, "chr7", 1.0, "-7"),
        )
        for event_id, event_type, chromosome, copy_number, notation in cases:
            with self.subTest(event_type=event_type):
                event = GenomicEvent(
                    event_id=event_id,
                    event_type=event_type,
                    primary=Locus(chromosome=chromosome, start=0, end=100),
                    copy_number=copy_number,
                    whole_chromosome_span_confirmed=True,
                    reportable=True,
                )
                proposal = build_iscn_proposal(
                    [event],
                    genome_build=GenomeBuild.GRCH38,
                    resource_provenance=_provenance(GenomeBuild.GRCH38),
                    cnv_source_event_ids={event.event_id},
                )

                self.assertEqual(proposal.notation, notation)
                self.assertEqual(
                    proposal.event_dispositions[0].outcome,
                    ISCNDispositionOutcome.RENDERED,
                )

    def test_unknown_typed_cnv_source_id_is_rejected(self) -> None:
        with self.assertRaisesRegex(
            ValueError, "CNV source IDs are absent from the proposal event set: ghost-cnv"
        ):
            build_iscn_proposal(
                [],
                genome_build=GenomeBuild.GRCH38,
                resource_provenance=_provenance(GenomeBuild.GRCH38),
                cnv_source_event_ids={"ghost-cnv"},
            )

    def test_proposal_rejects_resource_provenance_from_another_build(self) -> None:
        with self.assertRaisesRegex(ValueError, "uses a different genome build"):
            build_iscn_proposal(
                [],
                genome_build=GenomeBuild.GRCH38,
                resource_provenance=_provenance(GenomeBuild.GRCH37),
            )

    def test_real_non_reportable_cnv_can_form_a_technical_fragment(self) -> None:
        event = GenomicEvent(
            event_id="cnv7",
            event_type=EventType.DELETION,
            primary=Locus(
                chromosome="chr7",
                start=1,
                end=100,
                cytoband_start="q22",
                cytoband_end="q31.1",
            ),
            copy_number=1.0,
            reportable=False,
        )
        proposal = build_iscn_proposal(
            [event],
            genome_build=GenomeBuild.GRCH37,
            selection_policy=ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1,
            resource_provenance=_provenance(GenomeBuild.GRCH37),
            policy_parameters={"cytoband_affected_fraction": 0.66},
            cnv_source_event_ids={"cnv7"},
        )

        self.assertEqual(proposal.notation, "del(7)(q22q31.1)")
        self.assertEqual(proposal.genome_build, GenomeBuild.GRCH37)
        self.assertFalse(proposal.event_fragments[0].reportable)
        self.assertEqual(iscn_module_outcome(proposal).status, ModuleRunStatus.COMPLETED)

    def test_unvalidated_sv_is_dispositioned_but_not_formalized(self) -> None:
        event = GenomicEvent(
            event_id="sv-t821",
            event_type=EventType.TRANSLOCATION,
            primary=Locus(chromosome="chr8", start=1, end=2, cytoband_start="q22"),
            secondary=Locus(chromosome="chr21", start=3, end=4, cytoband_start="q22.1"),
            confidence="high",
            known_rearrangement="RUNX1::RUNX1T1",
            reportable=True,
        )
        proposal = build_iscn_proposal(
            [event],
            genome_build=GenomeBuild.GRCH38,
            selection_policy=ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1,
            resource_provenance=_provenance(GenomeBuild.GRCH38),
        )

        self.assertIsNone(proposal.notation)
        self.assertEqual(
            proposal.event_dispositions[0].outcome,
            ISCNDispositionOutcome.EXCLUDED_BY_POLICY,
        )
        self.assertIn("SV-to-ISCN", proposal.event_dispositions[0].detail)

    def test_missing_band_and_duplicate_fragments_are_explicit(self) -> None:
        events = [
            GenomicEvent(
                event_id="missing-band",
                event_type=EventType.DELETION,
                primary=Locus(chromosome="chr5", start=1, end=2),
                copy_number=1.0,
            ),
            GenomicEvent(
                event_id="gain8-a",
                event_type=EventType.DUPLICATION,
                primary=Locus(
                    chromosome="chr8",
                    start=0,
                    end=100,
                    cytoband_start="q24",
                    cytoband_end="q24",
                ),
                copy_number=3.0,
            ),
            GenomicEvent(
                event_id="gain8-b",
                event_type=EventType.DUPLICATION,
                primary=Locus(
                    chromosome="chr8",
                    start=100,
                    end=200,
                    cytoband_start="q24",
                    cytoband_end="q24",
                ),
                copy_number=3.0,
            ),
        ]
        proposal = build_iscn_proposal(
            events,
            genome_build=GenomeBuild.GRCH38,
            selection_policy=ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1,
            resource_provenance=_provenance(GenomeBuild.GRCH38),
            cnv_source_event_ids={"missing-band", "gain8-a", "gain8-b"},
        )

        self.assertEqual(proposal.notation, "dup(8)(q24)")
        dispositions = {item.event_id: item for item in proposal.event_dispositions}
        self.assertEqual(
            dispositions["missing-band"].outcome,
            ISCNDispositionOutcome.OMITTED_UNSUPPORTED,
        )
        self.assertEqual(dispositions["gain8-b"].reason_code, "DUPLICATE_FRAGMENT")

    def test_hard_observability_and_sex_chromosome_gates_precede_selection(self) -> None:
        events = [
            GenomicEvent(
                event_id="unobservable-reportable",
                event_type=EventType.CHROMOSOME_GAIN,
                primary=Locus(chromosome="chr8", start=0, end=100),
                copy_number=3.0,
                reportable=True,
                observability=SvObservability.OUTSIDE_TARGET,
            ),
            GenomicEvent(
                event_id="x-gain",
                event_type=EventType.CHROMOSOME_GAIN,
                primary=Locus(chromosome="chrX", start=0, end=100),
                copy_number=3.0,
            ),
        ]
        proposal = build_iscn_proposal(
            events,
            genome_build=GenomeBuild.GRCH38,
            selection_policy=ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1,
            resource_provenance=_provenance(GenomeBuild.GRCH38),
            cnv_source_event_ids={"unobservable-reportable", "x-gain"},
        )

        self.assertIsNone(proposal.notation)
        details = " ".join(item.detail for item in proposal.event_dispositions)
        self.assertIn("OUTSIDE_TARGET", details)
        self.assertIn("sex-chromosome", details)

    def test_sex_chromosome_gate_also_applies_to_reportable_policy(self) -> None:
        event = GenomicEvent(
            event_id="reportable-x-gain",
            event_type=EventType.CHROMOSOME_GAIN,
            primary=Locus(chromosome="chrX", start=0, end=100),
            copy_number=3.0,
            reportable=True,
        )
        proposal = build_iscn_proposal(
            [event],
            genome_build=GenomeBuild.GRCH38,
            resource_provenance=_provenance(GenomeBuild.GRCH38),
            cnv_source_event_ids={event.event_id},
        )

        self.assertIsNone(proposal.notation)
        self.assertEqual(
            proposal.event_dispositions[0].outcome,
            ISCNDispositionOutcome.EXCLUDED_BY_POLICY,
        )
        self.assertIn("sex-chromosome", proposal.event_dispositions[0].detail)

    def test_not_requested_and_missing_resources_are_distinct(self) -> None:
        event = GenomicEvent(
            event_id="candidate",
            event_type=EventType.CHROMOSOME_GAIN,
            primary=Locus(chromosome="chr8", start=0, end=100),
            copy_number=3.0,
        )
        not_requested = build_iscn_proposal(
            [event],
            genome_build=GenomeBuild.GRCH37,
            requested=False,
            upstream_evidence_assessed=False,
        )
        not_assessed = build_iscn_proposal(
            [event],
            genome_build=GenomeBuild.GRCH37,
            upstream_evidence_assessed=False,
        )

        self.assertEqual(not_requested.proposal_status, ISCNProposalStatus.NOT_REQUESTED)
        self.assertEqual(not_assessed.proposal_status, ISCNProposalStatus.NOT_ASSESSED)
        self.assertEqual(iscn_module_outcome(not_requested).status, ModuleRunStatus.NOT_RUN)
        self.assertEqual(iscn_module_outcome(not_assessed).status, ModuleRunStatus.NOT_RUN)
        self.assertTrue(not_assessed.assessment_blockers)
        self.assertEqual(
            not_requested.event_dispositions[0].outcome,
            ISCNDispositionOutcome.NOT_REQUESTED,
        )
        self.assertEqual(
            not_assessed.event_dispositions[0].outcome,
            ISCNDispositionOutcome.ASSESSMENT_BLOCKED,
        )

    def test_each_proposal_status_allows_only_its_disposition_outcomes(self) -> None:
        unsupported = GenomicEvent(
            event_id="unsupported",
            event_type=EventType.INSERTION,
            primary=Locus(chromosome="chr6", start=10, end=11),
            reportable=True,
        )
        renderable = GenomicEvent(
            event_id="renderable",
            event_type=EventType.DELETION,
            primary=Locus(
                chromosome="chr5",
                start=10,
                end=20,
                cytoband_start="q13",
                cytoband_end="q13",
            ),
            copy_number=1.0,
            reportable=True,
        )
        provenance = _provenance(GenomeBuild.GRCH38)
        not_requested = build_iscn_proposal(
            [unsupported],
            genome_build=GenomeBuild.GRCH38,
            requested=False,
        )
        not_assessed = build_iscn_proposal(
            [unsupported],
            genome_build=GenomeBuild.GRCH38,
            upstream_evidence_assessed=False,
        )
        no_renderable = build_iscn_proposal(
            [unsupported],
            genome_build=GenomeBuild.GRCH38,
            resource_provenance=provenance,
        )
        partial = build_iscn_proposal(
            [renderable, unsupported],
            genome_build=GenomeBuild.GRCH38,
            resource_provenance=provenance,
            cnv_source_event_ids={renderable.event_id},
        )

        mutations = (
            (not_requested, "unsupported", "EXCLUDED_BY_POLICY"),
            (not_assessed, "unsupported", "EXCLUDED_BY_POLICY"),
            (no_renderable, "unsupported", "NOT_REQUESTED"),
            (partial, "unsupported", "ASSESSMENT_BLOCKED"),
        )
        for proposal, event_id, invalid_outcome in mutations:
            with self.subTest(
                proposal_status=proposal.proposal_status,
                invalid_outcome=invalid_outcome,
            ):
                payload = proposal.model_dump(mode="json")
                disposition = next(
                    item for item in payload["event_dispositions"] if item["event_id"] == event_id
                )
                disposition["outcome"] = invalid_outcome
                with self.assertRaises(ValueError):
                    ISCNProposal.model_validate(payload)


if __name__ == "__main__":
    unittest.main()
