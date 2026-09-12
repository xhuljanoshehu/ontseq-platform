"""Tests for the panel-lock seam in the target coverage adapter.

The contract functions themselves are covered in test_panel_lock.py. What is checked here is
the seam: that the adapter consults a lock at all, that an unusable lock stops the run rather
than producing a plausible depth number, and that a design without provenance is reported
instead of being passed over in silence.
"""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ontseq_platform.models import TargetBedRole
from ontseq_platform.target_coverage import _verify_panel_contract

BED_TEXT = "# RESEARCH USE ONLY\nchr1\t100\t200\tGENE_A\nchr2\t300\t400\tGENE_B\n"


class PanelContractSeamTests(unittest.TestCase):
    def _bed(self, directory: str) -> Path:
        bed_path = Path(directory) / "panel.bed"
        bed_path.write_text(BED_TEXT, encoding="utf-8")
        return bed_path

    def test_missing_lock_is_reported_and_not_refused(self) -> None:
        """A synthetic fixture has no provenance; that is allowed but must be said out loud."""
        with TemporaryDirectory() as directory:
            notes = _verify_panel_contract(
                self._bed(directory), TargetBedRole.ANALYSIS_ROI_UNBUFFERED
            )
        self.assertEqual(len(notes), 1)
        self.assertIn("panel.bed", notes[0])
        self.assertIn("carries no panel lock", notes[0])

    def test_lock_beside_the_bed_is_actually_consulted(self) -> None:
        """An unusable lock must stop the run, which also proves the lock is read at all."""
        with TemporaryDirectory() as directory:
            bed_path = self._bed(directory)
            lock_path = Path(directory) / "panel.lock.yaml"
            lock_path.write_text("status: confirmed\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                _verify_panel_contract(bed_path, TargetBedRole.ANALYSIS_ROI_UNBUFFERED)

    def test_unparseable_lock_is_not_treated_as_absent(self) -> None:
        """Broken YAML must fail closed rather than fall back to the no-lock note."""
        with TemporaryDirectory() as directory:
            bed_path = self._bed(directory)
            lock_path = Path(directory) / "panel.lock.yaml"
            lock_path.write_text("bed: [unbalanced\n", encoding="utf-8")
            with self.assertRaises(Exception) as caught:
                _verify_panel_contract(bed_path, TargetBedRole.ANALYSIS_ROI_UNBUFFERED)
        self.assertNotIsInstance(caught.exception, AssertionError)


if __name__ == "__main__":
    unittest.main()
