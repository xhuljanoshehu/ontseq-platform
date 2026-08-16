from __future__ import annotations

import unittest

from ontseq_platform.cnv.mask import (
    ExclusionReason,
    ExclusionTrack,
    build_mask,
    coverage_floor_track,
    merge_masks,
    no_call_track,
    scope_from_intervals,
)

CONTIGS = {"chr1": 1_000_000, "chr2": 1_000_000}


class MaskConstructionTests(unittest.TestCase):
    def test_unmasked_genome_is_fully_evaluable(self) -> None:
        mask = build_mask(contig_lengths=CONTIGS)
        self.assertEqual(mask.reference_bases, 2_000_000)
        self.assertEqual(mask.evaluable_bases, 2_000_000)
        self.assertEqual(mask.excluded_bases, 0)
        self.assertEqual(mask.evaluable_fraction, 1.0)

    def test_analysis_scope_attributes_the_remainder_explicitly(self) -> None:
        mask = build_mask(contig_lengths=CONTIGS, analysis_scope={"chr1": [(0, 400_000)]})
        self.assertEqual(mask.evaluable_bases, 400_000)
        self.assertEqual(
            mask.excluded_bases_by_reason[ExclusionReason.OUTSIDE_ANALYSIS_SCOPE.value],
            1_600_000,
        )

    def test_exclusion_reasons_sum_to_the_total_removed(self) -> None:
        mask = build_mask(
            contig_lengths=CONTIGS,
            tracks=[
                ExclusionTrack(
                    reason=ExclusionReason.CENTROMERE,
                    intervals={"1": [(100_000, 200_000)]},
                    source="test:centromere",
                ),
                ExclusionTrack(
                    reason=ExclusionReason.BLACKLIST,
                    intervals={"1": [(150_000, 300_000)]},
                    source="test:blacklist",
                ),
            ],
        )
        self.assertEqual(sum(mask.excluded_bases_by_reason.values()), mask.excluded_bases)
        # Overlapping bases go to the first track that removed them.
        self.assertEqual(mask.excluded_bases_by_reason["centromere"], 100_000)
        self.assertEqual(mask.excluded_bases_by_reason["blacklist"], 100_000)

    def test_fully_masked_genome_warns_instead_of_reporting_zero_silently(self) -> None:
        mask = build_mask(
            contig_lengths={"chr1": 1000},
            tracks=[
                ExclusionTrack(
                    reason=ExclusionReason.BLACKLIST,
                    intervals={"1": [(0, 1000)]},
                    source="test",
                )
            ],
        )
        self.assertEqual(mask.evaluable_bases, 0)
        self.assertTrue(any("empty" in warning for warning in mask.warnings))

    def test_scope_contig_absent_from_reference_is_reported(self) -> None:
        mask = build_mask(contig_lengths={"chr1": 1000}, analysis_scope={"chr9": [(0, 500)]})
        self.assertEqual(mask.evaluable_bases, 0)
        self.assertTrue(any("chr9" in warning or "9" in warning for warning in mask.warnings))


class TrackConstructionTests(unittest.TestCase):
    def test_coverage_floor_excludes_only_regions_below_the_floor(self) -> None:
        track = coverage_floor_track(
            [("chr1", 0, 1000, 0.5), ("chr1", 1000, 2000, 12.0), ("chr2", 0, 500, 2.0)],
            minimum_depth=5.0,
            source="test:mosdepth",
        )
        self.assertEqual(track.reason, ExclusionReason.BELOW_COVERAGE_FLOOR)
        self.assertEqual(track.intervals, {"1": [(0, 1000)], "2": [(0, 500)]})

    def test_coverage_floor_rejects_inverted_intervals(self) -> None:
        with self.assertRaises(ValueError):
            coverage_floor_track([("chr1", 100, 100, 1.0)], minimum_depth=5.0, source="test")

    def test_no_call_track_defaults_to_caller_no_call(self) -> None:
        track = no_call_track({"chr1": [(0, 100)]}, source="test")
        self.assertEqual(track.reason, ExclusionReason.CALLER_NO_CALL)

    def test_scope_flank_expands_symmetrically_and_clamps_at_zero(self) -> None:
        self.assertEqual(
            scope_from_intervals({"chr1": [(100, 200)]}, flank_bp=500), {"1": [(0, 700)]}
        )

    def test_merging_masks_takes_the_intersection(self) -> None:
        first = build_mask(contig_lengths=CONTIGS, analysis_scope={"chr1": [(0, 500_000)]})
        second = build_mask(contig_lengths=CONTIGS, analysis_scope={"chr1": [(250_000, 1_000_000)]})
        self.assertEqual(merge_masks(first, second), {"1": [(250_000, 500_000)]})


if __name__ == "__main__":
    unittest.main()
