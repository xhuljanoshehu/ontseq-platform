from __future__ import annotations

import hashlib
import unittest
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
PANEL_STEM = "aml_fusion_adaptive_sampling.grch38.buffered"
BED_PATH = REPO_ROOT / "configs" / "panels" / f"{PANEL_STEM}.bed"
LOCK_PATH = REPO_ROOT / "configs" / "panels" / f"{PANEL_STEM}.lock.yaml"
EXPECTATIONS_PATH = REPO_ROOT / "configs" / "qc" / "target_coverage_expectations.grch38.tsv"
REVIEW_SUFFIX = "_REVIEW_REQUIRED"


def _bed_rows() -> list[tuple[str, int, int, str]]:
    rows: list[tuple[str, int, int, str]] = []
    for line in BED_PATH.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        chromosome, start, end, label = line.split("\t")
        rows.append((chromosome, int(start), int(end), label))
    return rows


def _expectation_rows() -> list[tuple[str, str, int, int]]:
    rows: list[tuple[str, str, int, int]] = []
    for line in EXPECTATIONS_PATH.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("gene\t"):
            continue
        fields = line.split("\t")
        rows.append((fields[0], fields[1], int(fields[2]), int(fields[3])))
    return rows


class PanelLockContractTests(unittest.TestCase):
    """The lock is the only place the panel may make claims about itself."""

    def setUp(self) -> None:
        self.lock = yaml.safe_load(LOCK_PATH.read_text(encoding="utf-8"))
        self.bed = self.lock["bed"]
        self.rows = _bed_rows()

    def test_lock_checksum_matches_the_committed_bed(self) -> None:
        digest = hashlib.sha256(BED_PATH.read_bytes()).hexdigest()
        self.assertEqual(self.bed["sha256"], digest)

    def test_interval_counts_and_bases_match_the_bed(self) -> None:
        self.assertEqual(len(self.rows), self.bed["target_count"])
        labels = [row[3] for row in self.rows]
        self.assertEqual(len(set(labels)), self.bed["unique_target_labels"])
        self.assertEqual(len(labels), len(set(labels)))
        bases = sum(row[2] - row[1] for row in self.rows)
        self.assertEqual(bases, self.bed["interval_bases"])

    def test_every_interval_is_forward_oriented(self) -> None:
        for chromosome, start, end, label in self.rows:
            self.assertLess(start, end, f"{label} on {chromosome} is not forward oriented")

    def test_bed_chromosomes_match_the_lock(self) -> None:
        self.assertEqual({row[0] for row in self.rows}, set(self.bed["chromosomes"]))

    def test_the_panel_does_not_claim_validated_genes(self) -> None:
        self.assertEqual(self.bed["target_type"], "target_intervals")
        self.assertIsNone(self.bed["validated_gene_count"])
        self.assertEqual(self.lock["status"], "derived_unconfirmed")
        self.assertEqual(self.lock["role"], "selection_panel_buffered")

    def test_promotion_blockers_are_recorded(self) -> None:
        blockers = self.lock["promotion_blockers"]
        self.assertTrue(blockers)
        self.assertTrue(all(isinstance(item, str) and item for item in blockers))

    def test_open_questions_stay_visible_in_the_bed(self) -> None:
        labels = {row[3] for row in self.rows}
        for question in self.lock["open_questions"]:
            gene = question["target"]
            self.assertIn(f"{gene}{REVIEW_SUFFIX}", labels)
            self.assertNotIn(gene, labels)


class CoverageExpectationContractTests(unittest.TestCase):
    """The descriptive expectations must not disagree with the panel they describe."""

    def setUp(self) -> None:
        self.lock = yaml.safe_load(LOCK_PATH.read_text(encoding="utf-8"))
        self.intervals = {row[3]: (row[0], row[1], row[2]) for row in _bed_rows()}
        self.rows = _expectation_rows()

    def test_expectations_are_not_empty(self) -> None:
        self.assertTrue(self.rows)

    def test_every_expectation_row_matches_a_panel_interval(self) -> None:
        for label, chromosome, start, end in self.rows:
            self.assertIn(label, self.intervals)
            self.assertEqual(self.intervals[label], (chromosome, start, end))

    def test_open_questions_stay_visible_in_the_expectations(self) -> None:
        labels = {row[0] for row in self.rows}
        for question in self.lock["open_questions"]:
            gene = question["target"]
            self.assertNotIn(gene, labels)
            self.assertIn(f"{gene}{REVIEW_SUFFIX}", labels)


if __name__ == "__main__":
    unittest.main()
