from __future__ import annotations

import unittest

from ontseq_platform.cnv.cytoband_regions import (
    CytobandEvent,
    CytobandRegionError,
    normalize_cytoband_regions,
)
from ontseq_platform.cnv.cytobands import Cytoband, CytobandTable
from ontseq_platform.cnv.states import CopyNumberState
from ontseq_platform.models import GenomeBuild


def _table() -> CytobandTable:
    return CytobandTable(
        resource_id="synthetic-cytoband-regions",
        genome_build=GenomeBuild.GRCH37,
        source="synthetic",
        bands=[
            Cytoband(contig="chr5", start=0, end=100, name="p15.2", stain="gneg"),
            Cytoband(contig="chr5", start=100, end=200, name="p15.1", stain="gpos25"),
            Cytoband(contig="chr5", start=200, end=300, name="p14", stain="gneg"),
            Cytoband(contig="chr5", start=300, end=400, name="p13.3", stain="gpos25"),
            Cytoband(contig="chr5", start=400, end=500, name="p13.2", stain="gneg"),
            Cytoband(contig="chr5", start=500, end=600, name="p13.1", stain="gpos50"),
            Cytoband(contig="chr5", start=600, end=700, name="q13.1", stain="gneg"),
            Cytoband(contig="chr5", start=700, end=800, name="q13.2", stain="gpos25"),
            Cytoband(contig="chr5", start=800, end=900, name="q13.3", stain="gneg"),
            Cytoband(contig="chr5", start=900, end=1000, name="q14.1", stain="gpos50"),
            Cytoband(contig="chr7", start=0, end=100, name="q31.21", stain="gneg"),
            Cytoband(contig="chr7", start=100, end=200, name="q31.22", stain="gpos25"),
            Cytoband(contig="chr7", start=200, end=300, name="q31.23", stain="gneg"),
            Cytoband(contig="chrX", start=0, end=100, name="p22.2", stain="gneg"),
            Cytoband(contig="chrX", start=100, end=200, name="p22.1", stain="gpos25"),
        ],
    )


class CytobandRegionTests(unittest.TestCase):
    def test_complete_subband_set_compacts_to_parent(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("5", "q13.1", CopyNumberState.LOSS),
                CytobandEvent("5", "q13.2", CopyNumberState.LOSS),
                CytobandEvent("5", "q13.3", CopyNumberState.LOSS),
            ],
            _table(),
        )
        self.assertEqual(len(regions), 1)
        region = regions[0]
        self.assertEqual((region.start, region.end), (600, 900))
        self.assertEqual(region.display_first_band, "q13")
        self.assertEqual(region.display_last_band, "q13")
        self.assertTrue(region.is_single_band_label)
        self.assertEqual(region.source_bands, ("q13.1", "q13.2", "q13.3"))

    def test_partial_subband_set_does_not_overcompact(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("5", "q13.2", CopyNumberState.LOSS),
                CytobandEvent("5", "q13.3", CopyNumberState.LOSS),
            ],
            _table(),
        )
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].display_first_band, "q13.2")
        self.assertEqual(regions[0].display_last_band, "q13.3")

    def test_two_decimal_levels_compact_one_level_when_complete(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("7", "q31.21", CopyNumberState.GAIN),
                CytobandEvent("7", "q31.22", CopyNumberState.GAIN),
                CytobandEvent("7", "q31.23", CopyNumberState.GAIN),
            ],
            _table(),
        )
        self.assertEqual(regions[0].display_first_band, "q31.2")
        self.assertEqual(regions[0].display_last_band, "q31.2")

    def test_adjacent_bands_with_different_states_do_not_merge(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("5", "q13.1", CopyNumberState.LOSS),
                CytobandEvent("5", "q13.2", CopyNumberState.GAIN),
            ],
            _table(),
        )
        self.assertEqual(len(regions), 2)

    def test_centromere_boundary_is_never_crossed(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("5", "p13.1", CopyNumberState.LOSS),
                CytobandEvent("5", "q13.1", CopyNumberState.LOSS),
            ],
            _table(),
        )
        self.assertEqual(len(regions), 2)
        self.assertEqual({region.arm for region in regions}, {"p", "q"})

    def test_p_arm_presentation_order_does_not_change_coordinates(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("5", "p13.3", CopyNumberState.LOSS),
                CytobandEvent("5", "p13.2", CopyNumberState.LOSS),
                CytobandEvent("5", "p13.1", CopyNumberState.LOSS),
            ],
            _table(),
        )
        self.assertEqual(len(regions), 1)
        region = regions[0]
        self.assertEqual((region.start, region.end), (300, 600))
        self.assertEqual(region.display_first_band, "p13")
        self.assertEqual(region.display_last_band, "p13")

    def test_p_arm_multi_region_labels_are_proximal_to_distal_for_display(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("5", "p15.1", CopyNumberState.GAIN),
                CytobandEvent("5", "p14", CopyNumberState.GAIN),
                CytobandEvent("5", "p13.3", CopyNumberState.GAIN),
            ],
            _table(),
        )
        self.assertEqual(len(regions), 1)
        region = regions[0]
        self.assertEqual((region.start, region.end), (100, 400))
        self.assertEqual(region.display_first_band, "p13.3")
        self.assertEqual(region.display_last_band, "p15.1")

    def test_non_adjacent_bands_remain_separate(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("5", "q13.1", CopyNumberState.LOSS),
                CytobandEvent("5", "q13.3", CopyNumberState.LOSS),
            ],
            _table(),
        )
        self.assertEqual(len(regions), 2)

    def test_parent_child_overlap_is_rejected(self) -> None:
        with self.assertRaises(CytobandRegionError) as context:
            normalize_cytoband_regions(
                [
                    CytobandEvent("5", "q13", CopyNumberState.LOSS),
                    CytobandEvent("5", "q13.2", CopyNumberState.LOSS),
                ],
                _table(),
            )
        self.assertIn("overlap", str(context.exception))

    def test_unknown_band_is_rejected(self) -> None:
        with self.assertRaises(CytobandRegionError):
            normalize_cytoband_regions(
                [CytobandEvent("5", "q99", CopyNumberState.LOSS)],
                _table(),
            )

    def test_copy_neutral_loh_is_outside_dosage_region_scope(self) -> None:
        with self.assertRaises(CytobandRegionError):
            normalize_cytoband_regions(
                [CytobandEvent("5", "q13.1", CopyNumberState.COPY_NEUTRAL_LOH)],
                _table(),
            )

    def test_source_ids_are_retained_but_paths_are_rejected(self) -> None:
        regions = normalize_cytoband_regions(
            [
                CytobandEvent("X", "p22.2", CopyNumberState.GAIN, source_id="event-1"),
                CytobandEvent("X", "p22.1", CopyNumberState.GAIN, source_id="event-2"),
            ],
            _table(),
        )
        self.assertEqual(regions[0].contig, "X")
        self.assertEqual(regions[0].source_ids, ("event-1", "event-2"))
        with self.assertRaises(CytobandRegionError):
            normalize_cytoband_regions(
                [
                    CytobandEvent(
                        "X",
                        "p22.2",
                        CopyNumberState.GAIN,
                        source_id="C:/patient/result.csv",
                    )
                ],
                _table(),
            )

    def test_empty_input_is_an_empty_normalization_not_a_biological_negative(self) -> None:
        self.assertEqual(normalize_cytoband_regions([], _table()), [])


if __name__ == "__main__":
    unittest.main()
