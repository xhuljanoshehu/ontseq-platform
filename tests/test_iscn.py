# ruff: noqa: I001
from __future__ import annotations

import unittest

from ontseq_platform.iscn import build_iscn_proposal
from ontseq_platform.iscn_reference import CytobandIndex
from ontseq_platform.iscn_validation import ISCNValidationStatus, validate_subset
from ontseq_platform.models import EventType, GenomeBuild, GenomicEvent, Locus


SYNTHETIC_CYTOBANDS = """\
chr1\t0\t10\tp13\tgneg
chr1\t10\t20\tp12\tgpos50
chr1\t20\t30\tp11.1\tacen
chr1\t30\t40\tq11.1\tacen
chr1\t40\t50\tq12\tgneg
chr1\t50\t60\tq13\tgpos50
chr9\t0\t100\tp11\tgneg
chr9\t100\t200\tq34.1\tgpos50
chr22\t0\t100\tp11\tgneg
chr22\t100\t200\tq11.2\tgpos50
"""


class ISCNProposalTests(unittest.TestCase):
    def test_subset_renderer_derives_chromosome_count(self) -> None:
        events = [
            GenomicEvent(
                event_id="gain8",
                event_type=EventType.CHROMOSOME_GAIN,
                primary=Locus(chromosome="chr8", start=0, end=1),
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
                reportable=True,
            ),
        ]
        proposal = build_iscn_proposal(
            events,
            sex_chromosomes="XY",
            prefer_external_validator=False,
        )
        self.assertEqual(proposal.notation, "47,XY,del(5)(q13q34),+8")
        self.assertEqual(proposal.source_event_ids, ["del5q", "gain8"])
        self.assertTrue(any("validation PASS" in warning for warning in proposal.warnings))

    def test_explicit_chromosome_count_remains_supported(self) -> None:
        event = GenomicEvent(
            event_id="gain8",
            event_type=EventType.CHROMOSOME_GAIN,
            primary=Locus(chromosome="chr8", start=0, end=1),
            reportable=True,
        )
        proposal = build_iscn_proposal(
            [event],
            chromosome_count=48,
            sex_chromosomes="XY",
            prefer_external_validator=False,
        )
        self.assertEqual(proposal.notation, "48,XY,+8")

    def test_coordinate_to_cytoband_mapping_for_p_arm_deletion(self) -> None:
        cytobands = CytobandIndex.from_tsv_text(SYNTHETIC_CYTOBANDS, GenomeBuild.GRCH37)
        event = GenomicEvent(
            event_id="del1p",
            event_type=EventType.DELETION,
            primary=Locus(chromosome="chr1", start=1, end=20),
            reportable=True,
        )
        proposal = build_iscn_proposal(
            [event],
            cytobands=cytobands,
            prefer_external_validator=False,
        )
        self.assertEqual(proposal.notation, "46,XX,del(1)(p12p13)")

    def test_translocation_breakpoints_can_be_mapped_from_coordinates(self) -> None:
        cytobands = CytobandIndex.from_tsv_text(SYNTHETIC_CYTOBANDS, GenomeBuild.GRCH37)
        event = GenomicEvent(
            event_id="t9_22",
            event_type=EventType.TRANSLOCATION,
            primary=Locus(chromosome="chr9", start=150, end=151),
            secondary=Locus(chromosome="chr22", start=150, end=151),
            reportable=True,
        )
        proposal = build_iscn_proposal(
            [event],
            cytobands=cytobands,
            prefer_external_validator=False,
        )
        self.assertEqual(proposal.notation, "46,XX,t(9;22)(q34.1;q11.2)")

    def test_cross_centromere_simple_deletion_is_not_guessed(self) -> None:
        cytobands = CytobandIndex.from_tsv_text(SYNTHETIC_CYTOBANDS, GenomeBuild.GRCH37)
        event = GenomicEvent(
            event_id="cross_centromere",
            event_type=EventType.DELETION,
            primary=Locus(chromosome="chr1", start=20, end=41),
            reportable=True,
        )
        proposal = build_iscn_proposal(
            [event],
            cytobands=cytobands,
            prefer_external_validator=False,
        )
        self.assertEqual(proposal.notation, "46,XX")
        self.assertTrue(any("cross_centromere" in warning for warning in proposal.warnings))

    def test_unsupported_event_is_not_silently_rendered(self) -> None:
        event = GenomicEvent(
            event_id="ins1",
            event_type=EventType.INSERTION,
            primary=Locus(chromosome="chr1", start=10, end=11),
            reportable=True,
        )
        proposal = build_iscn_proposal(
            [event],
            prefer_external_validator=False,
        )
        self.assertEqual(proposal.notation, "46,XX")
        self.assertTrue(any("ins1" in warning for warning in proposal.warnings))


class CytobandIndexTests(unittest.TestCase):
    def test_centromere_bounds_are_derived_from_acen(self) -> None:
        cytobands = CytobandIndex.from_tsv_text(SYNTHETIC_CYTOBANDS, GenomeBuild.GRCH37)
        self.assertEqual(cytobands.centromere_bounds("1"), (20, 40))

    def test_half_open_boundary_maps_to_next_band(self) -> None:
        cytobands = CytobandIndex.from_tsv_text(SYNTHETIC_CYTOBANDS, GenomeBuild.GRCH37)
        band = cytobands.band_at("chr1", 10)
        self.assertIsNotNone(band)
        assert band is not None
        self.assertEqual(band.name, "p12")


class ISCNValidationTests(unittest.TestCase):
    def test_known_subset_strings_pass(self) -> None:
        for notation in (
            "46,XX",
            "47,XY,+8",
            "46,XX,del(5)(q13q34)",
            "46,XX,inv(2)(p23p13)",
            "46,XX,t(9;22)(q34;q11.2)",
        ):
            with self.subTest(notation=notation):
                self.assertEqual(validate_subset(notation).status, ISCNValidationStatus.PASS)

    def test_whitespace_fails(self) -> None:
        result = validate_subset("47, XY,+8")
        self.assertEqual(result.status, ISCNValidationStatus.FAIL)

    def test_unknown_fragment_fails(self) -> None:
        result = validate_subset("46,XX,mar")
        self.assertEqual(result.status, ISCNValidationStatus.FAIL)


if __name__ == "__main__":
    unittest.main()
