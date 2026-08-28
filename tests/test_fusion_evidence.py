from __future__ import annotations

import unittest

from ontseq_platform.breakpoint_annotation import (
    AnnotatedBreakpoint,
    AnnotatedBreakpointPair,
    Breakpoint,
    TranscriptBreakpointHit,
)
from ontseq_platform.fusion_evidence import (
    FusionFrameStatus,
    fusion_evidence_from_breakpoints,
)


def _breakpoint(gene: str, transcript: str, *, region: str) -> AnnotatedBreakpoint:
    return AnnotatedBreakpoint(
        breakpoint=Breakpoint("chr1", 100),
        cytoband="q1",
        transcripts=(
            TranscriptBreakpointHit(
                gene_id=f"G_{gene}",
                gene_name=gene,
                transcript_id=transcript,
                transcript_name=transcript,
                strand="+",
                preferred=True,
                rank_tier=1,
                region=region,  # type: ignore[arg-type]
                exon_number=2 if region == "exon" else None,
                intron_number=2 if region == "intron" else None,
                cds_phase=0,
            ),
        ),
        contexts=(),
    )


class FusionEvidenceTests(unittest.TestCase):
    def test_gene_transcript_region_and_orientation_are_projected(self) -> None:
        evidence = fusion_evidence_from_breakpoints(
            AnnotatedBreakpointPair(
                primary=_breakpoint("BCR", "ENST_BCR", region="exon"),
                secondary=_breakpoint("ABL1", "ENST_ABL1", region="intron"),
            ),
            orientation="+-",
        )
        self.assertEqual(evidence.gene_a.gene, "BCR")
        self.assertEqual(evidence.gene_a.preferred_transcript, "ENST_BCR")
        self.assertEqual(evidence.gene_a.exon_number, 2)
        self.assertEqual(evidence.gene_b.gene, "ABL1")
        self.assertEqual(evidence.gene_b.intron_number, 2)
        self.assertEqual(evidence.orientation, "+-")
        self.assertEqual(evidence.frame_status, FusionFrameStatus.UNKNOWN)

    def test_ambiguous_orientation_and_missing_partner_stay_unknown(self) -> None:
        evidence = fusion_evidence_from_breakpoints(
            AnnotatedBreakpointPair(
                primary=_breakpoint("RUNX1", "ENST_RUNX1", region="transcript"),
                secondary=None,
            ),
            orientation="ambiguous",
        )
        self.assertIsNone(evidence.orientation)
        self.assertIsNone(evidence.gene_b.gene)
        self.assertEqual(evidence.gene_b.region, "unknown")
        self.assertEqual(evidence.frame_status.value, "unknown")


if __name__ == "__main__":
    unittest.main()
