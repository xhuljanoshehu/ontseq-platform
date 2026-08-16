from __future__ import annotations

import unittest

from ontseq_platform.cnv.intervals import (
    canonical_contig,
    contig_lengths_to_set,
    intersect,
    intersect_set,
    normalize,
    normalize_set,
    overlap_length,
    subtract,
    subtract_set,
    total_length,
    union_set,
)


class NormalizationTests(unittest.TestCase):
    def test_merges_overlapping_and_touching_intervals(self) -> None:
        self.assertEqual(normalize([(10, 20), (15, 30), (30, 40)]), [(10, 40)])

    def test_drops_empty_and_inverted_spans(self) -> None:
        self.assertEqual(normalize([(10, 10), (30, 20), (5, 8)]), [(5, 8)])

    def test_empty_input_yields_empty_output(self) -> None:
        self.assertEqual(normalize([]), [])

    def test_contig_prefix_is_stripped(self) -> None:
        self.assertEqual(canonical_contig("chr7"), "7")
        self.assertEqual(canonical_contig("7"), "7")
        self.assertEqual(canonical_contig("chrX"), "X")

    def test_set_normalization_unifies_prefixed_contigs(self) -> None:
        result = union_set({"chr1": [(0, 100)]}, {"1": [(100, 200)]})
        self.assertEqual(result, {"1": [(0, 200)]})


class AlgebraTests(unittest.TestCase):
    def test_intersect(self) -> None:
        self.assertEqual(intersect([(0, 100), (200, 300)], [(50, 250)]), [(50, 100), (200, 250)])

    def test_intersect_without_overlap(self) -> None:
        self.assertEqual(intersect([(0, 10)], [(20, 30)]), [])

    def test_subtract_punches_a_hole(self) -> None:
        self.assertEqual(subtract([(0, 100)], [(40, 60)]), [(0, 40), (60, 100)])

    def test_subtract_removes_fully_covered_span(self) -> None:
        self.assertEqual(subtract([(10, 20)], [(0, 100)]), [])

    def test_subtract_with_nothing_to_remove(self) -> None:
        self.assertEqual(subtract([(0, 100)], []), [(0, 100)])

    def test_subtract_handles_multiple_removals_in_one_span(self) -> None:
        self.assertEqual(
            subtract([(0, 100)], [(10, 20), (30, 40), (90, 200)]),
            [(0, 10), (20, 30), (40, 90)],
        )

    def test_set_operations_are_contig_aware(self) -> None:
        left = {"1": [(0, 100)], "2": [(0, 100)]}
        right = {"1": [(50, 150)]}
        self.assertEqual(intersect_set(left, right), {"1": [(50, 100)]})
        self.assertEqual(subtract_set(left, right), {"1": [(0, 50)], "2": [(0, 100)]})

    def test_total_length_does_not_silently_merge(self) -> None:
        # Deliberately naive: double counting must stay visible to the caller.
        self.assertEqual(total_length({"1": [(0, 100), (50, 150)]}), 200)
        self.assertEqual(total_length(normalize_set({"1": [(0, 100), (50, 150)]})), 150)

    def test_overlap_length_for_a_single_span(self) -> None:
        self.assertEqual(overlap_length((0, 100), "chr1", {"1": [(50, 200)]}), 50)
        self.assertEqual(overlap_length((0, 100), "chr9", {"1": [(50, 200)]}), 0)

    def test_contig_lengths_to_set(self) -> None:
        self.assertEqual(contig_lengths_to_set({"chr1": 100, "chr2": 0}), {"1": [(0, 100)]})


if __name__ == "__main__":
    unittest.main()
