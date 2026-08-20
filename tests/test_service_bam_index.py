from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.service.guard import GuardError, resolve_bam_index


class BamIndexResolutionTests(unittest.TestCase):
    def test_bam_dot_bai_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            index = Path(f"{bam}.bai")
            bam.write_bytes(b"BAM")
            index.write_bytes(b"BAI")
            self.assertEqual(resolve_bam_index(bam), index)

    def test_short_bai_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            index = bam.with_suffix(".bai")
            bam.write_bytes(b"BAM")
            index.write_bytes(b"BAI")
            self.assertEqual(resolve_bam_index(bam), index)

    def test_bam_dot_bai_has_precedence_when_both_exist(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            preferred = Path(f"{bam}.bai")
            alternative = bam.with_suffix(".bai")
            for path in (bam, preferred, alternative):
                path.write_bytes(b"x")
            self.assertEqual(resolve_bam_index(bam), preferred)

    def test_missing_index_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            bam.write_bytes(b"BAM")
            with self.assertRaises(GuardError) as caught:
                resolve_bam_index(bam)
            message = str(caught.exception)
            self.assertIn("sample.bam.bai", message)
            self.assertIn("sample.bai", message)

    def test_unrelated_bai_is_not_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "sample.bam"
            bam.write_bytes(b"BAM")
            (root / "other.bai").write_bytes(b"BAI")
            with self.assertRaises(GuardError):
                resolve_bam_index(bam)


if __name__ == "__main__":
    unittest.main()
