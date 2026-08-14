from __future__ import annotations

import unittest

from ontseq_platform.iscn import build_iscn_proposal
from ontseq_platform.models import EventType, GenomicEvent, Locus


class ISCNProposalTests(unittest.TestCase):
    def test_subset_renderer_is_deterministic(self) -> None:
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
        proposal = build_iscn_proposal(events, sex_chromosomes="XY")
        self.assertEqual(proposal.notation, "46,XY,del(5)(q13q34),+8")
        self.assertEqual(proposal.source_event_ids, ["del5q", "gain8"])

    def test_unsupported_event_is_not_silently_rendered(self) -> None:
        event = GenomicEvent(
            event_id="ins1",
            event_type=EventType.INSERTION,
            primary=Locus(chromosome="chr1", start=10, end=11),
            reportable=True,
        )
        proposal = build_iscn_proposal([event])
        self.assertEqual(proposal.notation, "46,XX")
        self.assertTrue(any("ins1" in warning for warning in proposal.warnings))


if __name__ == "__main__":
    unittest.main()
