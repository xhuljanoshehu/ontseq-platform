"""Seam tests for the fitted-cellularity grading in the QDNAseq/ACE lane.

The boundaries under test are clinician-specified and are not validated on this
platform. These tests only assert that the configured boundaries are applied,
recorded and escalated consistently; they make no claim about the biology.
"""

from __future__ import annotations

import unittest

from pydantic import ValidationError

from ontseq_platform.cnv.qdnaseq import (
    CnvFit,
    QDNAseqPolicy,
    _cellularity_tier,
    _cellularity_warnings,
)


def _policy(**overrides: object) -> QDNAseqPolicy:
    values: dict[str, object] = {
        "profile_id": "cnv-unit-test-v1",
        "note": "technical defaults for unit tests only",
    }
    values.update(overrides)
    return QDNAseqPolicy(**values)


def _fit(bin_size_kbp: int, cellularity: float) -> CnvFit:
    return CnvFit(
        bin_size_kbp=bin_size_kbp,
        cellularity=cellularity,
        ploidy=2.0,
        fit_error=0.25,
        candidate_count=3,
        segment_count=12,
        segment_file=f"SYNTH.{bin_size_kbp}kbp.segments.tsv",
        chromosome_file=f"SYNTH.{bin_size_kbp}kbp.chromosomes.tsv",
        fit_plot=f"SYNTH.{bin_size_kbp}kbp.ace-fit.png",
        copy_number_plot=f"SYNTH.{bin_size_kbp}kbp.copy-number.png",
        rds_file=f"SYNTH.{bin_size_kbp}kbp.segmented.rds",
    )


class CellularityBoundaryTests(unittest.TestCase):
    def test_the_boundaries_are_carried_in_the_policy(self) -> None:
        policy = _policy()
        self.assertAlmostEqual(policy.cellularity_review_fraction, 0.20)
        self.assertAlmostEqual(policy.cellularity_critical_fraction, 0.10)

    def test_the_critical_boundary_must_lie_below_the_review_boundary(self) -> None:
        with self.assertRaises(ValidationError):
            _policy(cellularity_review_fraction=0.20, cellularity_critical_fraction=0.30)

    def test_the_tier_helper_reports_the_configured_bands(self) -> None:
        policy = _policy()
        self.assertEqual(_cellularity_tier(0.05, policy), "critical")
        self.assertEqual(_cellularity_tier(0.15, policy), "review")
        self.assertEqual(_cellularity_tier(0.45, policy), "unflagged")


class CellularityGradingTests(unittest.TestCase):
    def test_a_confident_fit_is_not_flagged(self) -> None:
        primary = _fit(500, 0.42)
        self.assertEqual(_cellularity_warnings(primary, [primary], _policy()), [])

    def test_a_low_fit_is_declared_weakly_constrained(self) -> None:
        primary = _fit(500, 0.14)
        warnings = _cellularity_warnings(primary, [primary], _policy())
        self.assertEqual(len(warnings), 1)
        self.assertIn("0.140", warnings[0])
        self.assertIn("weakly constrained", warnings[0])
        self.assertIn("not a measured tumour content", warnings[0])
        self.assertNotIn("makes no claim", warnings[0])

    def test_a_very_low_fit_escalates_without_repeating_the_softer_warning(self) -> None:
        primary = _fit(500, 0.06)
        warnings = _cellularity_warnings(primary, [primary], _policy())
        self.assertEqual(len(warnings), 1)
        self.assertIn("0.060", warnings[0])
        self.assertIn("makes no claim", warnings[0])
        self.assertIn("orthogonal confirmation", warnings[0])
        self.assertNotIn("weakly constrained", warnings[0])

    def test_disagreeing_tiers_across_bin_sizes_are_reported(self) -> None:
        primary = _fit(500, 0.42)
        fits = [_fit(100, 0.06), primary, _fit(1000, 0.14)]
        warnings = _cellularity_warnings(primary, fits, _policy())
        self.assertEqual(len(warnings), 1)
        self.assertIn("different confidence tiers", warnings[0])
        self.assertIn("100 kbp: critical", warnings[0])
        self.assertIn("1000 kbp: review", warnings[0])

    def test_agreeing_tiers_do_not_add_a_disagreement_warning(self) -> None:
        primary = _fit(500, 0.42)
        fits = [_fit(100, 0.38), primary, _fit(1000, 0.55)]
        self.assertEqual(_cellularity_warnings(primary, fits, _policy()), [])


if __name__ == "__main__":
    unittest.main()
