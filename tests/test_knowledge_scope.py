from __future__ import annotations

import unittest

from ontseq_platform.knowledge.scope import (
    Interval,
    MatchType,
    Origin,
    ScopeAlignment,
    align,
    canonical_contig,
    classify_match,
    reciprocal_overlap,
    review_stars,
)


class ContigTests(unittest.TestCase):
    """A prefix mismatch produces zero annotations with no error, which is the worst case.

    It looks exactly like a sample nobody knows anything about, so it must be impossible.
    """

    def test_a_chr_prefix_is_ignored(self) -> None:
        self.assertEqual(canonical_contig("chr7"), canonical_contig("7"))

    def test_sex_chromosomes_normalise_to_upper_case(self) -> None:
        self.assertEqual(canonical_contig("chrx"), "X")
        self.assertEqual(canonical_contig("Y"), "Y")

    def test_an_autosome_keeps_its_number(self) -> None:
        self.assertEqual(canonical_contig("chr22"), "22")


class AlignmentTests(unittest.TestCase):
    def test_matching_origins_align(self) -> None:
        alignment, _ = align(Origin.SOMATIC, Origin.SOMATIC)
        self.assertIs(alignment, ScopeAlignment.ALIGNED)

    def test_a_germline_record_under_a_somatic_question_is_mismatched(self) -> None:
        """The central case: ClinVar is germline, AML is somatic."""
        alignment, note = align(Origin.GERMLINE, Origin.SOMATIC)
        self.assertIs(alignment, ScopeAlignment.MISMATCHED)
        self.assertIn("does not answer the question", note)

    def test_a_mismatch_is_not_described_as_wrong(self) -> None:
        """A germline pathogenic finding under a somatic assay is a secondary finding."""
        _, note = align(Origin.GERMLINE, Origin.SOMATIC)
        self.assertIn("secondary finding", note)

    def test_an_undeclared_source_origin_is_unknown_not_agreement(self) -> None:
        alignment, note = align(Origin.UNKNOWN, Origin.SOMATIC)
        self.assertIs(alignment, ScopeAlignment.UNKNOWN)
        self.assertIn("the source record", note)

    def test_an_undeclared_assay_intent_is_unknown_and_says_which_side(self) -> None:
        alignment, note = align(Origin.GERMLINE, Origin.UNKNOWN)
        self.assertIs(alignment, ScopeAlignment.UNKNOWN)
        self.assertIn("the assay", note)


class ReviewStarTests(unittest.TestCase):
    def test_an_expert_panel_outranks_a_single_submitter(self) -> None:
        self.assertGreater(
            review_stars("reviewed by expert panel"),
            review_stars("criteria provided, single submitter"),
        )

    def test_a_practice_guideline_is_the_top_rating(self) -> None:
        self.assertEqual(review_stars("practice guideline"), 4)

    def test_no_assertion_criteria_is_zero_stars(self) -> None:
        self.assertEqual(review_stars("no assertion criteria provided"), 0)

    def test_an_unknown_status_is_none_not_zero(self) -> None:
        """None means this code is stale; zero means NCBI itself rated it lowest."""
        self.assertIsNone(review_stars("some future ncbi wording"))

    def test_matching_ignores_case_and_padding(self) -> None:
        self.assertEqual(review_stars("  Practice Guideline "), 4)


class OverlapTests(unittest.TestCase):
    def test_different_contigs_never_overlap(self) -> None:
        left = Interval("7", 0, 1000)
        right = Interval("8", 0, 1000)
        self.assertEqual(reciprocal_overlap(left, right), 0.0)

    def test_the_smaller_fraction_is_reported(self) -> None:
        """Otherwise a huge finding matches a small record at 1.0 by containing it."""
        finding = Interval("7", 0, 100_000_000)
        record = Interval("7", 0, 1_000_000)
        self.assertAlmostEqual(reciprocal_overlap(finding, record), 0.01)

    def test_identical_intervals_overlap_completely(self) -> None:
        same = reciprocal_overlap(Interval("7", 10, 20), Interval("7", 10, 20))
        self.assertAlmostEqual(same, 1.0)

    def test_touching_intervals_do_not_overlap(self) -> None:
        self.assertEqual(reciprocal_overlap(Interval("7", 0, 10), Interval("7", 10, 20)), 0.0)

    def test_a_zero_length_interval_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            Interval("7", 10, 10)

    def test_a_negative_coordinate_is_refused(self) -> None:
        """ClinVar's -1 placeholder passes an end-after-start check: -1 > -2."""
        with self.assertRaises(ValueError):
            Interval("7", -2, -1)


class MatchTests(unittest.TestCase):
    def _classify(self, finding: Interval, record: Interval, **kwargs):
        options = {"minimum_reciprocal_overlap": 0.5, "exact_tolerance_bp": 10_000}
        options.update(kwargs)
        return classify_match(finding, record, **options)

    def test_near_identical_breakpoints_are_exact(self) -> None:
        """Array and read-depth breakpoints are never bit-identical."""
        result = self._classify(
            Interval("7", 1_000_000, 5_000_000), Interval("7", 1_005_000, 4_995_000)
        )
        assert result is not None
        self.assertIs(result[0], MatchType.EXACT)

    def test_a_small_record_inside_a_large_finding_is_named_as_such(self) -> None:
        result = self._classify(
            Interval("7", 0, 100_000_000), Interval("7", 50_000_000, 53_000_000)
        )
        assert result is not None
        self.assertIs(result[0], MatchType.RECORD_WITHIN_FINDING)

    def test_a_small_finding_inside_a_large_record_is_named_as_such(self) -> None:
        result = self._classify(
            Interval("7", 50_000_000, 53_000_000), Interval("7", 0, 100_000_000)
        )
        assert result is not None
        self.assertIs(result[0], MatchType.FINDING_WITHIN_RECORD)

    def test_containment_is_reported_even_at_low_overlap(self) -> None:
        """A 3 Mb record inside a 90 Mb arm loss says little, and the type says so."""
        result = self._classify(Interval("7", 0, 90_000_000), Interval("7", 10_000_000, 13_000_000))
        assert result is not None
        self.assertIs(result[0], MatchType.RECORD_WITHIN_FINDING)
        self.assertLess(result[1], 0.1)

    def test_partial_overlap_below_the_threshold_does_not_match(self) -> None:
        result = self._classify(Interval("7", 0, 10_000_000), Interval("7", 9_000_000, 30_000_000))
        self.assertIsNone(result)

    def test_partial_overlap_above_the_threshold_matches(self) -> None:
        result = self._classify(Interval("7", 0, 10_000_000), Interval("7", 2_000_000, 11_000_000))
        assert result is not None
        self.assertIs(result[0], MatchType.OVERLAP)

    def test_no_overlap_is_no_match(self) -> None:
        self.assertIsNone(self._classify(Interval("7", 0, 1000), Interval("7", 5000, 6000)))

    def test_a_contig_mismatch_is_no_match(self) -> None:
        self.assertIsNone(self._classify(Interval("7", 0, 1000), Interval("8", 0, 1000)))

    def test_a_chr_prefix_does_not_prevent_a_match(self) -> None:
        result = self._classify(Interval("chr7", 0, 1000), Interval("7", 0, 1000))
        assert result is not None
        self.assertIs(result[0], MatchType.EXACT)

    def test_an_impossible_threshold_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            self._classify(
                Interval("7", 0, 1000), Interval("7", 0, 1000), minimum_reciprocal_overlap=0.0
            )


if __name__ == "__main__":
    unittest.main()
