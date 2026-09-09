from __future__ import annotations

import unittest
from pathlib import Path

from ontseq_platform.quantitation import (
    COPY_NEUTRAL_DIPLOID,
    DEPTH_TOO_LOW,
    NO_VARIANT_READS,
    NOT_COPY_NEUTRAL,
    VAF_EXCEEDS_HETEROZYGOUS_MODEL,
    AlleleObservation,
    QuantitationError,
    cancer_cell_fraction,
    copy_number_from_ratio,
    expected_vaf,
    format_report,
    minimum_detectable_cancer_cell_fraction,
    minimum_detectable_vaf,
    minimum_variant_reads,
    tumour_fraction_from_clonal_snv,
    wilson_interval,
)

NEUTRAL = COPY_NEUTRAL_DIPLOID


class AlleleObservationTests(unittest.TestCase):
    def test_vaf(self) -> None:
        self.assertAlmostEqual(AlleleObservation(20, 80).vaf, 0.25)

    def test_rejects_more_variant_than_total_reads(self) -> None:
        with self.assertRaises(QuantitationError):
            AlleleObservation(9, 8)

    def test_rejects_zero_depth(self) -> None:
        with self.assertRaises(QuantitationError):
            AlleleObservation(0, 0)

    def test_rejects_negative_variant_reads(self) -> None:
        with self.assertRaises(QuantitationError):
            AlleleObservation(-1, 10)


class WilsonIntervalTests(unittest.TestCase):
    def test_interval_contains_the_point_estimate(self) -> None:
        observation = AlleleObservation(20, 80)
        band = wilson_interval(observation)
        self.assertLess(band.low, observation.vaf)
        self.assertGreater(band.high, observation.vaf)

    def test_interval_stays_inside_zero_to_one_with_no_variant_reads(self) -> None:
        """The reason this is Wilson and not the normal interval.

        At zero variant reads the normal interval collapses to a point at zero, claiming
        certainty from the one observation that carries least.
        """
        band = wilson_interval(AlleleObservation(0, 30))
        self.assertAlmostEqual(band.low, 0.0)
        self.assertGreater(band.high, 0.0)
        self.assertLess(band.high, 1.0)

    def test_deeper_coverage_narrows_the_interval(self) -> None:
        shallow = wilson_interval(AlleleObservation(5, 20))
        deep = wilson_interval(AlleleObservation(50, 200))
        self.assertLess(deep.width, shallow.width)

    def test_rejects_confidence_outside_the_open_unit_interval(self) -> None:
        with self.assertRaises(QuantitationError):
            wilson_interval(AlleleObservation(1, 10), confidence=1.0)


class TumourFractionTests(unittest.TestCase):
    def test_vaf_of_one_quarter_is_a_half_tumour(self) -> None:
        """The identity the whole module rests on: f = 2 x VAF."""
        estimate = tumour_fraction_from_clonal_snv(
            AlleleObservation(20, 80), copy_number_state=NEUTRAL
        )
        self.assertTrue(estimate.determinable)
        assert estimate.point is not None
        self.assertAlmostEqual(estimate.point, 0.5)

    def test_the_interval_is_reported_and_brackets_the_point(self) -> None:
        estimate = tumour_fraction_from_clonal_snv(
            AlleleObservation(20, 80), copy_number_state=NEUTRAL
        )
        assert estimate.interval is not None and estimate.point is not None
        self.assertLess(estimate.interval.low, estimate.point)
        self.assertGreater(estimate.interval.high, estimate.point)
        self.assertLessEqual(estimate.interval.high, 1.0)

    def test_refuses_outside_a_copy_neutral_region(self) -> None:
        """The central guard: under altered copy number the identity silently breaks."""
        estimate = tumour_fraction_from_clonal_snv(
            AlleleObservation(20, 80), copy_number_state="single_copy_deletion"
        )
        self.assertFalse(estimate.determinable)
        self.assertEqual(estimate.status, NOT_COPY_NEUTRAL)
        self.assertIsNone(estimate.point)
        self.assertIn("single_copy_deletion", estimate.reason())

    def test_refuses_below_the_requested_depth(self) -> None:
        estimate = tumour_fraction_from_clonal_snv(
            AlleleObservation(3, 9), copy_number_state=NEUTRAL
        )
        self.assertEqual(estimate.status, DEPTH_TOO_LOW)
        self.assertIsNone(estimate.point)

    def test_no_variant_reads_is_not_a_tumour_fraction_of_zero(self) -> None:
        estimate = tumour_fraction_from_clonal_snv(
            AlleleObservation(0, 80), copy_number_state=NEUTRAL
        )
        self.assertEqual(estimate.status, NO_VARIANT_READS)
        self.assertIsNone(estimate.point)
        self.assertIn("not evidence", estimate.reason())

    def test_refuses_a_vaf_that_implies_more_than_a_whole_sample(self) -> None:
        """VAF 0.5 is a germline heterozygous SNP, not a tumour fraction of 1."""
        estimate = tumour_fraction_from_clonal_snv(
            AlleleObservation(60, 80), copy_number_state=NEUTRAL
        )
        self.assertEqual(estimate.status, VAF_EXCEEDS_HETEROZYGOUS_MODEL)
        self.assertIsNone(estimate.point)


class CancerCellFractionTests(unittest.TestCase):
    def test_clonal_variant_gives_a_fraction_of_one(self) -> None:
        result = cancer_cell_fraction(
            AlleleObservation(20, 80), tumour_fraction=0.5, copy_number_state=NEUTRAL
        )
        self.assertTrue(result.determinable)
        assert result.point is not None
        self.assertAlmostEqual(result.point, 1.0)

    def test_half_the_allele_fraction_is_half_the_tumour(self) -> None:
        result = cancer_cell_fraction(
            AlleleObservation(10, 80), tumour_fraction=0.5, copy_number_state=NEUTRAL
        )
        assert result.point is not None
        self.assertAlmostEqual(result.point, 0.5)
        self.assertIs(result.subclonal_at(0.9), True)

    def test_undeterminable_fraction_does_not_read_as_clonal(self) -> None:
        """``None`` rather than ``False``: "cannot tell" must not collapse into "clonal"."""
        result = cancer_cell_fraction(
            AlleleObservation(10, 80), tumour_fraction=0.5, copy_number_state="amplification"
        )
        self.assertFalse(result.determinable)
        self.assertIsNone(result.subclonal_at(0.9))

    def test_rejects_an_impossible_tumour_fraction(self) -> None:
        with self.assertRaises(QuantitationError):
            cancer_cell_fraction(
                AlleleObservation(10, 80), tumour_fraction=0.0, copy_number_state=NEUTRAL
            )


class ForwardModelTests(unittest.TestCase):
    def test_forward_and_inverse_agree(self) -> None:
        for fraction in (0.2, 0.5, 0.85, 1.0):
            for clone in (0.25, 0.5, 1.0):
                vaf = expected_vaf(tumour_fraction=fraction, cancer_cell_fraction=clone)
                observation = AlleleObservation(round(vaf * 1000), 1000)
                recovered = cancer_cell_fraction(
                    observation, tumour_fraction=fraction, copy_number_state=NEUTRAL
                )
                assert recovered.point is not None
                self.assertAlmostEqual(recovered.point, clone, places=2)

    def test_a_clonal_variant_in_pure_tumour_sits_at_one_half(self) -> None:
        self.assertAlmostEqual(expected_vaf(tumour_fraction=1.0, cancer_cell_fraction=1.0), 0.5)


class CopyNumberTests(unittest.TestCase):
    def test_dilution_masks_a_single_copy_deletion(self) -> None:
        """The reason a fixed ratio threshold is unsafe.

        At half tumour content a real one-copy deletion shows a ratio of 0.75, which a
        threshold set for pure tumour would read as normal.
        """
        self.assertAlmostEqual(copy_number_from_ratio(ratio=0.75, tumour_fraction=0.5), 1.0)

    def test_pure_tumour_reproduces_the_textbook_ratios(self) -> None:
        self.assertAlmostEqual(copy_number_from_ratio(ratio=0.5, tumour_fraction=1.0), 1.0)
        self.assertAlmostEqual(copy_number_from_ratio(ratio=1.0, tumour_fraction=1.0), 2.0)
        self.assertAlmostEqual(copy_number_from_ratio(ratio=1.5, tumour_fraction=1.0), 3.0)

    def test_a_neutral_ratio_is_two_copies_at_any_purity(self) -> None:
        for fraction in (0.1, 0.5, 1.0):
            self.assertAlmostEqual(copy_number_from_ratio(ratio=1.0, tumour_fraction=fraction), 2.0)

    def test_rejects_a_negative_ratio(self) -> None:
        with self.assertRaises(QuantitationError):
            copy_number_from_ratio(ratio=-0.1, tumour_fraction=0.5)


class DetectionLimitTests(unittest.TestCase):
    def test_deeper_coverage_lowers_the_detectable_allele_fraction(self) -> None:
        shallow = minimum_detectable_vaf(30)
        deep = minimum_detectable_vaf(200)
        assert shallow is not None and deep is not None
        self.assertLess(deep, shallow)

    def test_the_error_floor_needs_more_reads_as_depth_grows(self) -> None:
        self.assertLessEqual(minimum_variant_reads(30), minimum_variant_reads(200))

    def test_off_target_depth_resolves_no_subclone(self) -> None:
        """The measured off-target depth of this assay, stated as a limit rather than a hope."""
        self.assertIsNone(minimum_detectable_cancer_cell_fraction(9, tumour_fraction=0.5))

    def test_on_target_depth_resolves_only_large_subclones(self) -> None:
        resolvable = minimum_detectable_cancer_cell_fraction(80, tumour_fraction=0.5)
        assert resolvable is not None
        self.assertGreater(resolvable, 0.3)
        self.assertLess(resolvable, 0.5)

    def test_a_lower_tumour_fraction_resolves_fewer_subclones(self) -> None:
        rich = minimum_detectable_cancer_cell_fraction(200, tumour_fraction=0.9)
        poor = minimum_detectable_cancer_cell_fraction(200, tumour_fraction=0.3)
        assert rich is not None and poor is not None
        self.assertLess(rich, poor)

    def test_rejects_impossible_parameters(self) -> None:
        with self.assertRaises(QuantitationError):
            minimum_variant_reads(0)
        with self.assertRaises(QuantitationError):
            minimum_variant_reads(80, error_rate=1.0)
        with self.assertRaises(QuantitationError):
            minimum_variant_reads(80, alpha=0.0)


class GeneratedReportTests(unittest.TestCase):
    """The document is generated so its numbers cannot drift from the code that computes them."""

    PATH = Path("docs/QUANTITATIVE_MODEL.md")

    def test_report_states_both_measured_depths(self) -> None:
        text = format_report()
        self.assertIn("measured on-target", text)
        self.assertIn("measured off-target", text)
        self.assertIn("GENERATED FILE", text)

    def test_report_says_the_caller_is_not_implemented(self) -> None:
        self.assertIn("no small-variant caller", format_report())

    def test_checked_in_report_is_current(self) -> None:
        self.assertTrue(self.PATH.exists(), "run python -m ontseq_platform.quantitation")
        self.assertEqual(
            self.PATH.read_text(encoding="utf-8"),
            format_report() + "\n",
            "docs/QUANTITATIVE_MODEL.md is stale; regenerate it",
        )


if __name__ == "__main__":
    unittest.main()
