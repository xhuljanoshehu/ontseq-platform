from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml

from ontseq_platform.pipeline.envelope import sha256_file
from ontseq_platform.pipeline.panel_lock import (
    PanelLockError,
    check_declared_role,
    load_panel_lock,
    panel_usage_warnings,
    target_labels,
    verify_panel_bed,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
PANEL_STEM = "aml_fusion_adaptive_sampling.grch38.buffered"
REAL_BED = REPO_ROOT / "configs" / "panels" / f"{PANEL_STEM}.bed"
REAL_LOCK = REPO_ROOT / "configs" / "panels" / f"{PANEL_STEM}.lock.yaml"

BED_TEXT = "# synthetic fixture\nchr1\t100\t200\tAAA\nchr2\t300\t400\tBBB_REVIEW_REQUIRED\n"


def _write_panel(
    directory: Path,
    *,
    status: str = "derived_unconfirmed",
    role: str = "selection_panel_buffered",
    bed_text: str = BED_TEXT,
    blockers: list[str] | None = None,
    validated_gene_count: int | None = None,
    target_count: int | None = None,
) -> tuple[Path, Path]:
    bed_path = directory / "panel.bed"
    bed_path.write_text(bed_text, encoding="utf-8")
    labels = target_labels(bed_path)
    document = {
        "panel_version": "SYNTHETIC_V1",
        "status": status,
        "genome_build": "GRCh38",
        "role": role,
        "bed": {
            "path": bed_path.name,
            "sha256": sha256_file(bed_path),
            "target_type": "target_intervals",
            "target_count": len(labels) if target_count is None else target_count,
            "unique_target_labels": len(set(labels)),
            "validated_gene_count": validated_gene_count,
        },
        "open_questions": [{"target": "BBB", "detail": "label and coordinates disagree"}],
        "promotion_blockers": ["The original design has not been byte-compared."]
        if blockers is None
        else blockers,
    }
    lock_path = directory / "panel.lock.yaml"
    lock_path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return lock_path, bed_path


class LoadPanelLockTests(unittest.TestCase):
    def test_fields_are_read_and_absence_is_preserved(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(Path(raw))
            lock = load_panel_lock(lock_path)
        self.assertEqual(lock.panel_version, "SYNTHETIC_V1")
        self.assertEqual(lock.status, "derived_unconfirmed")
        self.assertEqual(lock.target_count, 2)
        self.assertEqual(lock.unique_target_labels, 2)
        self.assertIsNone(lock.validated_gene_count)
        self.assertEqual(lock.open_question_targets, ("BBB",))
        self.assertFalse(lock.confirmed)

    def test_a_confirmed_status_with_a_blocker_is_still_not_confirmed(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(Path(raw), status="confirmed")
            lock = load_panel_lock(lock_path)
        self.assertEqual(lock.status, "confirmed")
        self.assertFalse(lock.confirmed)

    def test_a_confirmed_status_without_blockers_is_confirmed(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(Path(raw), status="confirmed", blockers=[])
            lock = load_panel_lock(lock_path)
        self.assertTrue(lock.confirmed)

    def test_a_missing_identity_field_is_refused(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(Path(raw))
            document = yaml.safe_load(lock_path.read_text(encoding="utf-8"))
            del document["role"]
            lock_path.write_text(yaml.safe_dump(document), encoding="utf-8")
            with self.assertRaises(PanelLockError) as raised:
                load_panel_lock(lock_path)
        self.assertIn("'role'", str(raised.exception))

    def test_a_malformed_digest_is_refused(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(Path(raw))
            document = yaml.safe_load(lock_path.read_text(encoding="utf-8"))
            document["bed"]["sha256"] = "NOTAHASH"
            lock_path.write_text(yaml.safe_dump(document), encoding="utf-8")
            with self.assertRaises(PanelLockError):
                load_panel_lock(lock_path)


class VerifyPanelBedTests(unittest.TestCase):
    def test_matching_bed_is_accepted(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, bed_path = _write_panel(Path(raw))
            verify_panel_bed(load_panel_lock(lock_path), bed_path)

    def test_edited_bed_is_refused(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, bed_path = _write_panel(Path(raw))
            lock = load_panel_lock(lock_path)
            bed_path.write_text(BED_TEXT.replace("100", "101"), encoding="utf-8")
            with self.assertRaises(PanelLockError) as raised:
                verify_panel_bed(lock, bed_path)
        self.assertIn("provenance cannot be established", str(raised.exception))

    def test_missing_bed_is_refused(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, bed_path = _write_panel(Path(raw))
            lock = load_panel_lock(lock_path)
            bed_path.unlink()
            with self.assertRaises(PanelLockError):
                verify_panel_bed(lock, bed_path)

    def test_a_count_that_disagrees_with_the_bed_is_refused(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, bed_path = _write_panel(Path(raw), target_count=99)
            with self.assertRaises(PanelLockError) as raised:
                verify_panel_bed(load_panel_lock(lock_path), bed_path)
        self.assertIn("records 99", str(raised.exception))

    def test_an_unlabelled_interval_is_refused(self) -> None:
        with TemporaryDirectory() as raw:
            directory = Path(raw)
            lock_path, bed_path = _write_panel(directory)
            lock = load_panel_lock(lock_path)
            bed_path.write_text("chr1\t1\t2\n", encoding="utf-8")
            with self.assertRaises(PanelLockError):
                verify_panel_bed(lock, bed_path)


class DeclaredRoleTests(unittest.TestCase):
    def test_matching_role_passes(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(Path(raw))
            check_declared_role(load_panel_lock(lock_path), "selection_panel_buffered")

    def test_contradicting_role_fails_closed(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(Path(raw))
            lock = load_panel_lock(lock_path)
        with self.assertRaises(PanelLockError) as raised:
            check_declared_role(lock, "analysis_roi_unbuffered")
        message = str(raised.exception)
        self.assertIn("selection_panel_buffered", message)
        self.assertIn("analysis_roi_unbuffered", message)


class UsageWarningTests(unittest.TestCase):
    def test_an_unconfirmed_panel_explains_itself(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, bed_path = _write_panel(Path(raw))
            lock = load_panel_lock(lock_path)
            warnings = panel_usage_warnings(lock, labels=target_labels(bed_path))
        joined = "\n".join(warnings)
        self.assertIn("derived_unconfirmed", joined)
        self.assertIn("promotion blocker", joined.lower())
        self.assertIn("claims no validated gene count", joined)
        self.assertIn("buffered", joined)
        self.assertIn("BBB_REVIEW_REQUIRED", joined)

    def test_a_confirmed_panel_stays_quiet_about_status(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(
                Path(raw),
                status="confirmed",
                role="analysis_roi_unbuffered",
                blockers=[],
                validated_gene_count=2,
            )
            warnings = panel_usage_warnings(load_panel_lock(lock_path))
        self.assertFalse(any("recorded as" in warning for warning in warnings))
        self.assertFalse(any("design is buffered" in warning for warning in warnings))
        self.assertTrue(any("BBB_REVIEW_REQUIRED" in warning for warning in warnings))

    def test_a_review_label_absent_from_the_lock_is_reported(self) -> None:
        with TemporaryDirectory() as raw:
            lock_path, _ = _write_panel(Path(raw))
            lock = load_panel_lock(lock_path)
            warnings = panel_usage_warnings(lock, labels=["CCC_REVIEW_REQUIRED"])
        self.assertTrue(any("does not list as an open question" in item for item in warnings))


class CommittedPanelTests(unittest.TestCase):
    """The module must accept the panel this repository actually ships."""

    def test_the_committed_panel_satisfies_its_own_lock(self) -> None:
        lock = load_panel_lock(REAL_LOCK)
        verify_panel_bed(lock, REAL_BED)
        self.assertEqual(lock.genome_build, "GRCh38")
        self.assertFalse(lock.confirmed)
        self.assertIsNone(lock.validated_gene_count)
        warnings = panel_usage_warnings(lock, labels=target_labels(REAL_BED))
        self.assertTrue(warnings)
        self.assertTrue(any("REVIEW_REQUIRED" in item for item in warnings))


if __name__ == "__main__":
    unittest.main()
