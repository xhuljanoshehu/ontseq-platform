from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from pydantic import ValidationError

from ontseq_platform.models import GenomeBuild
from ontseq_platform.reference import reference_lock_from_fai


class ReferenceLockTests(unittest.TestCase):
    def test_lock_records_contigs_and_exact_fai_checksum(self) -> None:
        content = "chr1\t1000\t0\t80\t81\nchr2\t2000\t1013\t80\t81\n"
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "synthetic.fa.fai"
            fai.write_text(content, encoding="utf-8")

            lock = reference_lock_from_fai(
                fai,
                reference_id="SYNTHETIC_REF",
                genome_build=GenomeBuild.GRCH38,
            )

        self.assertEqual([item.name for item in lock.contigs], ["chr1", "chr2"])
        self.assertEqual([item.length for item in lock.contigs], [1000, 2000])
        self.assertEqual(lock.source_fai_sha256, hashlib.sha256(content.encode()).hexdigest())

    def test_duplicate_contigs_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "duplicate.fa.fai"
            fai.write_text("chr1\t1000\nchr1\t1000\n", encoding="utf-8")

            with self.assertRaises(ValidationError):
                reference_lock_from_fai(
                    fai,
                    reference_id="SYNTHETIC_REF",
                    genome_build=GenomeBuild.GRCH38,
                )

    def test_malformed_fai_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "malformed.fa.fai"
            fai.write_text("chr1\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "expected at least 2 fields"):
                reference_lock_from_fai(
                    fai,
                    reference_id="SYNTHETIC_REF",
                    genome_build=GenomeBuild.GRCH38,
                )


if __name__ == "__main__":
    unittest.main()
