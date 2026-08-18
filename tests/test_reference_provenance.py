from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.models import GenomeBuild
from ontseq_platform.reference import contig_signature, sha256_file
from ontseq_platform.reference_provenance import (
    build_sequence_reference_lock,
    verify_sequence_reference_lock,
)


def _write_reference(root: Path) -> tuple[Path, Path]:
    fasta = root / "synthetic.fa"
    fai = root / "synthetic.fa.fai"
    fasta.write_text(">chr1\nAAAA\n>chr2\nCCCC\n", encoding="utf-8")
    fai.write_text("chr1\t4\t6\t4\t5\nchr2\t4\t17\t4\t5\n", encoding="utf-8")
    return fasta, fai


class SequenceReferenceLockTests(unittest.TestCase):
    def test_lock_records_sequence_fai_and_contig_fingerprints_without_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, fai = _write_reference(root)

            lock = build_sequence_reference_lock(
                fasta,
                fai,
                reference_id="SYNTHETIC_REF",
                genome_build=GenomeBuild.GRCH38,
            )

            self.assertEqual(lock.source_fasta_sha256, sha256_file(fasta))
            self.assertEqual(lock.source_fai_sha256, sha256_file(fai))
            self.assertEqual(
                lock.contig_signature_sha256,
                contig_signature([("chr1", 4), ("chr2", 4)]),
            )
            serialized = lock.model_dump_json()
            self.assertNotIn(str(fasta), serialized)
            self.assertNotIn(str(fai), serialized)
            self.assertNotIn("AAAA", serialized)
            self.assertNotIn("CCCC", serialized)

    def test_same_contigs_but_changed_fasta_sequence_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, fai = _write_reference(root)
            lock = build_sequence_reference_lock(
                fasta,
                fai,
                reference_id="SYNTHETIC_REF",
                genome_build=GenomeBuild.GRCH38,
            )

            fasta.write_text(">chr1\nTTTT\n>chr2\nCCCC\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "FASTA SHA256"):
                verify_sequence_reference_lock(
                    lock,
                    fasta_path=fasta,
                    fai_path=fai,
                    reference_id="SYNTHETIC_REF",
                    genome_build=GenomeBuild.GRCH38,
                )

    def test_changed_fai_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, fai = _write_reference(root)
            lock = build_sequence_reference_lock(
                fasta,
                fai,
                reference_id="SYNTHETIC_REF",
                genome_build=GenomeBuild.GRCH38,
            )

            fai.write_text("chr1\t4\t6\t4\t5\nchr2\t5\t17\t4\t5\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "FAI SHA256"):
                verify_sequence_reference_lock(
                    lock,
                    fasta_path=fasta,
                    fai_path=fai,
                    reference_id="SYNTHETIC_REF",
                    genome_build=GenomeBuild.GRCH38,
                )

    def test_declared_reference_identity_mismatch_fails_before_hash_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, fai = _write_reference(root)
            lock = build_sequence_reference_lock(
                fasta,
                fai,
                reference_id="SYNTHETIC_REF",
                genome_build=GenomeBuild.GRCH38,
            )

            with self.assertRaisesRegex(ValueError, "Reference ID"):
                verify_sequence_reference_lock(
                    lock,
                    fasta_path=fasta,
                    fai_path=fai,
                    reference_id="OTHER_REF",
                    genome_build=GenomeBuild.GRCH38,
                )

            with self.assertRaisesRegex(ValueError, "Genome build"):
                verify_sequence_reference_lock(
                    lock,
                    fasta_path=fasta,
                    fai_path=fai,
                    reference_id="SYNTHETIC_REF",
                    genome_build=GenomeBuild.GRCH37,
                )


if __name__ == "__main__":
    unittest.main()
