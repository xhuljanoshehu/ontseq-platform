from __future__ import annotations

import unittest

from ontseq_platform.cnv.cytoband import CnvSegment, Cytoband, annotate_cnv_cytobands

BANDS = [
    Cytoband("chr1", 0, 100, "p12"),
    Cytoband("chr1", 100, 200, "p11"),
    Cytoband("chr1", 200, 300, "q11"),
    Cytoband("chr1", 300, 400, "q12"),
]


class CytobandFractionTests(unittest.TestCase):
    def test_exact_66_percent_is_affected_and_all_raw_overlaps_remain(self) -> None:
        result = annotate_cnv_cytobands(
            [
                CnvSegment("exact", "chr1", 34, 100, "loss"),
                CnvSegment("below", "chr1", 135, 200, "loss"),
            ],
            BANDS,
            affected_fraction=0.66,
        )
        overlaps = {item.event_id: item for item in result.raw_overlaps}
        self.assertAlmostEqual(overlaps["exact"].fraction_of_band, 0.66)
        self.assertTrue(overlaps["exact"].affected)
        self.assertAlmostEqual(overlaps["below"].fraction_of_band, 0.65)
        self.assertFalse(overlaps["below"].affected)
        self.assertEqual(result.affected_groups[0].bands, ("p12",))

    def test_exact_synthetic_band_boundary_is_zero_based_half_open(self) -> None:
        result = annotate_cnv_cytobands(
            [CnvSegment("right-band", "chr1", 100, 200, "gain")],
            BANDS,
            affected_fraction=0.66,
        )

        self.assertEqual(len(result.raw_overlaps), 1)
        overlap = result.raw_overlaps[0]
        self.assertEqual(overlap.band, "p11")
        self.assertEqual(overlap.overlap_bp, 100)
        self.assertEqual(overlap.fraction_of_band, 1.0)
        self.assertEqual(result.affected_groups[0].bands, ("p11",))

    def test_one_base_on_each_side_of_synthetic_boundary_remains_traceable(self) -> None:
        result = annotate_cnv_cytobands(
            [CnvSegment("boundary-span", "1", 99, 101, "loss")],
            BANDS,
            affected_fraction=0.66,
        )

        self.assertEqual(
            [(item.band, item.overlap_bp) for item in result.raw_overlaps],
            [("p12", 1), ("p11", 1)],
        )
        self.assertTrue(all(not item.affected for item in result.raw_overlaps))
        self.assertEqual(result.affected_groups, ())


class CytobandGroupingTests(unittest.TestCase):
    def test_adjacent_same_direction_bands_merge(self) -> None:
        result = annotate_cnv_cytobands(
            [CnvSegment("p-loss", "1", 0, 200, "loss")],
            BANDS,
            affected_fraction=0.66,
        )
        self.assertEqual(len(result.affected_groups), 1)
        self.assertEqual(result.affected_groups[0].bands, ("p12", "p11"))
        self.assertEqual(result.affected_groups[0].source_event_ids, ("p-loss",))

    def test_an_unaffected_band_breaks_a_group(self) -> None:
        bands = [
            Cytoband("chr1", 0, 100, "p13"),
            Cytoband("chr1", 100, 200, "p12"),
            Cytoband("chr1", 200, 300, "p11"),
        ]
        result = annotate_cnv_cytobands(
            [
                CnvSegment("left", "chr1", 0, 100, "gain"),
                CnvSegment("right", "chr1", 200, 300, "gain"),
            ],
            bands,
            affected_fraction=0.66,
        )
        self.assertEqual([group.bands for group in result.affected_groups], [("p13",), ("p11",)])

    def test_p_and_q_arms_never_merge_across_centromere(self) -> None:
        result = annotate_cnv_cytobands(
            [CnvSegment("arm-span", "chr1", 0, 400, "gain")],
            BANDS,
            affected_fraction=0.66,
        )
        self.assertEqual(
            [group.bands for group in result.affected_groups],
            [("p12", "p11"), ("q11", "q12")],
        )
        self.assertEqual([group.arm for group in result.affected_groups], ["p", "q"])

    def test_opposite_directions_never_merge(self) -> None:
        result = annotate_cnv_cytobands(
            [
                CnvSegment("loss", "chr1", 0, 100, "loss"),
                CnvSegment("gain", "chr1", 100, 200, "gain"),
            ],
            BANDS,
            affected_fraction=0.66,
        )
        self.assertEqual(len(result.affected_groups), 2)
        self.assertEqual({group.direction for group in result.affected_groups}, {"loss", "gain"})


class WholeChromosomeTests(unittest.TestCase):
    def test_whole_chromosome_event_is_separate_from_band_groups(self) -> None:
        result = annotate_cnv_cytobands(
            [CnvSegment("whole", "chr1", 0, 400, "gain")],
            BANDS,
            affected_fraction=0.66,
            chromosome_sizes={"1": 400},
        )
        self.assertEqual(result.affected_groups, ())
        self.assertEqual(len(result.whole_chromosome_calls), 1)
        self.assertEqual(result.whole_chromosome_calls[0].chromosome, "chr1")
        self.assertEqual(result.whole_chromosome_calls[0].source_event_ids, ("whole",))
        self.assertEqual(len(result.raw_overlaps), 4)
        self.assertTrue(all(overlap.whole_chromosome for overlap in result.raw_overlaps))

    def test_focal_segment_is_not_promoted_to_whole_chromosome(self) -> None:
        result = annotate_cnv_cytobands(
            [CnvSegment("focal", "chr1", 0, 399, "gain")],
            BANDS,
            affected_fraction=0.66,
            chromosome_sizes={"chr1": 400},
        )
        self.assertEqual(result.whole_chromosome_calls, ())
        self.assertTrue(result.affected_groups)


if __name__ == "__main__":
    unittest.main()
