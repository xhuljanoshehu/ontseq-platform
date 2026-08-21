from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from pydantic import ValidationError

from ontseq_platform.models import GenomeBuild
from ontseq_platform.reference import (
    reference_lock_from_fai,
    reference_lock_signature,
    validate_canonical_reference,
)

GRCH37_LENGTHS = (
    249250621,
    243199373,
    198022430,
    191154276,
    180915260,
    171115067,
    159138663,
    146364022,
    141213431,
    135534747,
    135006516,
    133851895,
    115169878,
    107349540,
    102531392,
    90354753,
    81195210,
    78077248,
    59128983,
    63025520,
    48129895,
    51304566,
    155270560,
    59373566,
)
GRCH38_LENGTHS = (
    248956422,
    242193529,
    198295559,
    190214555,
    181538259,
    170805979,
    159345973,
    145138636,
    138394717,
    133797422,
    135086622,
    133275309,
    114364328,
    107043718,
    101991189,
    90338345,
    83257441,
    80373285,
    58617616,
    64444167,
    46709983,
    50818468,
    156040895,
    57227415,
)


def _canonical_fai(lengths: tuple[int, ...], *, prefix: str = "") -> str:
    labels = [*(str(number) for number in range(1, 23)), "X", "Y"]
    return "".join(
        f"{prefix}{label}\t{length}\t0\t80\t81\n"
        for label, length in zip(labels, lengths, strict=True)
    )


class ReferenceLockTests(unittest.TestCase):
    def test_lock_records_contigs_and_exact_fai_checksum(self) -> None:
        content = "chr1\t1000\t0\t80\t81\nchr2\t2000\t1013\t80\t81\n"
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "synthetic.fa.fai"
            fai.write_bytes(content.encode())

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

    def test_complete_chr_prefixed_grch38_is_accepted_for_desktop(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "grch38.fa.fai"
            content = _canonical_fai(GRCH38_LENGTHS, prefix="chr") + "chrM\t16569\n"
            fai.write_bytes(content.encode())
            lock = reference_lock_from_fai(
                fai,
                reference_id="GRCh38_LOCAL_TEST",
                genome_build=GenomeBuild.GRCH38,
                require_canonical_assembly=True,
            )

        summary = validate_canonical_reference(
            ((item.name, item.length) for item in lock.contigs), lock.genome_build
        )
        self.assertEqual(summary.naming_style, "chr-prefixed")
        self.assertEqual(summary.contig_count, 25)

    def test_complete_unprefixed_grch37_allows_decoy_contigs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "hs37d5.fa.fai"
            content = _canonical_fai(GRCH37_LENGTHS) + "MT\t16569\n\nhs37d5\t35477943\n"
            fai.write_bytes(content.encode())
            lock = reference_lock_from_fai(
                fai,
                reference_id="GRCh37_LOCAL_TEST",
                genome_build=GenomeBuild.GRCH37,
                require_canonical_assembly=True,
            )

        self.assertEqual(len(lock.contigs), 26)

    def test_region_only_reference_cannot_be_labeled_full_grch37(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "GRCh37_chr1.fasta.fai"
            fai.write_bytes(b"1\t249250621\t0\t80\t81\n")
            with self.assertRaisesRegex(ValueError, "Partial or wrong-build"):
                reference_lock_from_fai(
                    fai,
                    reference_id="GRCh37_LOCAL_TEST",
                    genome_build=GenomeBuild.GRCH37,
                    require_canonical_assembly=True,
                )

    def test_mixed_canonical_naming_styles_are_rejected(self) -> None:
        mixed = _canonical_fai(GRCH38_LENGTHS, prefix="chr") + _canonical_fai(GRCH37_LENGTHS)
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "mixed.fa.fai"
            fai.write_bytes(mixed.encode())
            with self.assertRaisesRegex(ValueError, "mixes canonical naming styles"):
                reference_lock_from_fai(
                    fai,
                    reference_id="GRCh38_MIXED_TEST",
                    genome_build=GenomeBuild.GRCH38,
                    require_canonical_assembly=True,
                )

    def test_wrong_named_build_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "24 length mismatches"):
            validate_canonical_reference(
                (
                    (name, length)
                    for name, length in zip(
                        [*(str(number) for number in range(1, 23)), "X", "Y"],
                        GRCH37_LENGTHS,
                        strict=True,
                    )
                ),
                GenomeBuild.GRCH38,
            )

    def test_lock_signature_covers_allow_extra_policy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fai = Path(temporary) / "synthetic.fa.fai"
            fai.write_bytes(b"chr1\t1000\n")
            strict_lock = reference_lock_from_fai(
                fai,
                reference_id="SYNTHETIC_REF",
                genome_build=GenomeBuild.GRCH38,
            )
            permissive_lock = strict_lock.model_copy(update={"allow_extra_contigs": True})

        self.assertNotEqual(
            reference_lock_signature(strict_lock),
            reference_lock_signature(permissive_lock),
        )


if __name__ == "__main__":
    unittest.main()
