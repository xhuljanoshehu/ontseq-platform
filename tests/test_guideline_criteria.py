from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.guideline_criteria import (
    COMPUTABLE,
    DRAFT,
    NEEDS_SMALL_VARIANTS,
    VERIFIED,
    GuidelineCriteriaError,
    load_for_review,
    load_reportable_criteria,
    risk_group_determinable,
)

SHIPPED = Path("configs/knowledge_bundles/GUIDELINE_CRITERIA_DRAFT_v0/guideline_criteria.v0.1.json")


def _record(record_id: str, *, verification: str, assay_status: str = COMPUTABLE) -> dict:
    return {
        "record_id": record_id,
        "category": "eln2022_adverse",
        "display_name": f"criterion {record_id}",
        "pattern_type": "fusion_pair",
        "detectable_by": ["sv_breakpoint"],
        "assay_status": assay_status,
        "verification": verification,
        "guideline_reference": "Table 5" if verification == VERIFIED else None,
        "reviewer_note": "check the wording",
        "caveat": "not a rule set",
    }


def _bundle(directory: Path, records: list[dict], *, reviewer: str | None = None) -> Path:
    path = directory / "criteria.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": "0.1.0",
                "bundle_id": "TEST_v0",
                "provenance": {"reviewer": reviewer, "review_date": None},
                "records": records,
            }
        ),
        encoding="utf-8",
    )
    return path


class DraftsMustNotReachAReport(unittest.TestCase):
    """The point of the module: an unverified draft is refused, not quietly filtered."""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.base = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_a_draft_record_is_refused_and_named(self) -> None:
        path = _bundle(self.base, [_record("ELN-X", verification=DRAFT)], reviewer="Dr. Muster")
        with self.assertRaises(GuidelineCriteriaError) as raised:
            load_reportable_criteria(path)
        self.assertIn("ELN-X", str(raised.exception))

    def test_one_draft_among_verified_records_still_refuses_the_whole_bundle(self) -> None:
        """A silently shortened criteria table produces a silently wrong classification."""
        path = _bundle(
            self.base,
            [_record("OK-1", verification=VERIFIED), _record("BAD-1", verification=DRAFT)],
            reviewer="Dr. Muster",
        )
        with self.assertRaises(GuidelineCriteriaError):
            load_reportable_criteria(path)

    def test_verified_records_need_a_named_reviewer(self) -> None:
        path = _bundle(self.base, [_record("OK-1", verification=VERIFIED)], reviewer=None)
        with self.assertRaises(GuidelineCriteriaError):
            load_reportable_criteria(path)

    def test_a_fully_verified_bundle_loads(self) -> None:
        path = _bundle(self.base, [_record("OK-1", verification=VERIFIED)], reviewer="Dr. Muster")
        bundle = load_reportable_criteria(path)
        self.assertEqual(len(bundle.criteria), 1)
        self.assertEqual(bundle.reviewer, "Dr. Muster")

    def test_review_tooling_can_still_read_the_draft(self) -> None:
        path = _bundle(self.base, [_record("ELN-X", verification=DRAFT)], reviewer=None)
        self.assertEqual(len(load_for_review(path).criteria), 1)


class RiskGroupIsNotGuessedFromWhatWeHappenToSee(unittest.TestCase):
    def test_a_criterion_the_assay_cannot_evaluate_blocks_the_risk_group(self) -> None:
        bundle = load_for_review(SHIPPED)
        verdict = risk_group_determinable(bundle.criteria)
        self.assertFalse(verdict.determinable)
        self.assertIn("Not determinable", verdict.reason())

    def test_the_reason_names_what_is_missing(self) -> None:
        bundle = load_for_review(SHIPPED)
        verdict = risk_group_determinable(bundle.criteria)
        blocking = {item.record_id for item in verdict.blocking}
        self.assertIn("ELN2022-ADV-TP53", blocking)
        self.assertIn("ELN2022-ADV-MDS-RELATED-GENES", blocking)

    def test_it_is_determinable_only_when_every_criterion_is_evaluable(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = _bundle(
                Path(raw),
                [_record("A", verification=VERIFIED), _record("B", verification=VERIFIED)],
                reviewer="Dr. Muster",
            )
            self.assertTrue(risk_group_determinable(load_for_review(path).criteria).determinable)


class TheShippedDraftSaysWhatItIs(unittest.TestCase):
    """Guards against the draft ever being marked verified without a human doing it."""

    def test_every_shipped_record_is_still_an_unverified_draft(self) -> None:
        bundle = load_for_review(SHIPPED)
        self.assertEqual(len(bundle.by_verification(DRAFT)), len(bundle.criteria))
        self.assertIsNone(bundle.reviewer)

    def test_the_shipped_draft_cannot_be_loaded_for_reporting(self) -> None:
        with self.assertRaises(GuidelineCriteriaError):
            load_reportable_criteria(SHIPPED)

    def test_small_variant_criteria_are_marked_as_not_implemented(self) -> None:
        bundle = load_for_review(SHIPPED)
        npm1 = next(item for item in bundle.criteria if "NPM1" in item.display_name)
        self.assertEqual(npm1.assay_status, NEEDS_SMALL_VARIANTS)

    def test_every_record_carries_a_caveat(self) -> None:
        for item in load_for_review(SHIPPED).criteria:
            self.assertTrue(item.caveat, item.record_id)


if __name__ == "__main__":
    unittest.main()
