from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ontseq_platform.models import (
    GenomeBuild,
    ReferenceContig,
    ReferenceDictionaryContract,
    ReferenceLock,
)
from ontseq_platform.reference import (
    grch38_canonical_25_contigs,
    reference_lock_for_dictionary_contract,
    reference_lock_from_fai,
    reference_lock_signature,
    validate_canonical_reference,
    validate_grch38_canonical_25,
)

GRCH38_NUCLEAR: tuple[int, ...] = (
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

GRCH37_NUCLEAR: tuple[int, ...] = (
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

LABELS = (*(str(number) for number in range(1, 23)), "X", "Y")


def _profile(lengths: tuple[int, ...], prefix: str) -> list[tuple[str, int]]:
    pairs = zip(LABELS, lengths, strict=True)
    return [(f"{prefix}{label}", length) for label, length in pairs]


def _write_fai(directory: Path, contigs: list[tuple[str, int]]) -> Path:
    path = directory / "reference.fa.fai"
    rows: list[str] = []
    offset = 0
    for name, length in contigs:
        rows.append(f"{name}\t{length}\t{offset}\t60\t61")
        offset += length
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


class CanonicalReferenceTests(unittest.TestCase):
    """A named build must be a complete build, whatever the sequence-name style."""

    def test_complete_chr_prefixed_grch38_is_accepted(self) -> None:
        contigs = [*_profile(GRCH38_NUCLEAR, "chr"), ("chrM", 16569)]
        summary = validate_canonical_reference(contigs, GenomeBuild.GRCH38)
        self.assertEqual(summary.naming_style, "chr-prefixed")
        self.assertEqual(summary.contig_count, 25)
        self.assertEqual(summary.total_reference_bases, sum(item[1] for item in contigs))

    def test_complete_unprefixed_grch37_with_a_decoy_is_accepted(self) -> None:
        contigs = [*_profile(GRCH37_NUCLEAR, ""), ("hs37d5", 35477943)]
        summary = validate_canonical_reference(contigs, GenomeBuild.GRCH37)
        self.assertEqual(summary.naming_style, "unprefixed")
        self.assertEqual(summary.contig_count, 25)

    def test_region_restricted_reference_cannot_claim_a_named_build(self) -> None:
        contigs = [("chr5", 181538259), ("chr22", 50818468)]
        with self.assertRaises(ValueError) as raised:
            validate_canonical_reference(contigs, GenomeBuild.GRCH38)
        message = str(raised.exception)
        self.assertIn("canonical assembly validation failed", message)
        self.assertIn("22 missing", message)

    def test_wrong_build_lengths_are_rejected(self) -> None:
        with self.assertRaises(ValueError) as raised:
            validate_canonical_reference(_profile(GRCH37_NUCLEAR, "chr"), GenomeBuild.GRCH38)
        self.assertIn("24 length mismatches", str(raised.exception))

    def test_mixed_naming_styles_are_rejected(self) -> None:
        contigs = [*_profile(GRCH38_NUCLEAR, "chr"), *_profile(GRCH38_NUCLEAR, "")]
        with self.assertRaises(ValueError) as raised:
            validate_canonical_reference(contigs, GenomeBuild.GRCH38)
        self.assertIn("mixes canonical naming styles", str(raised.exception))

    def test_duplicate_contig_names_are_rejected(self) -> None:
        contigs = [("chr1", 248956422), ("chr1", 248956422)]
        with self.assertRaises(ValueError) as raised:
            validate_canonical_reference(contigs, GenomeBuild.GRCH38)
        self.assertIn("duplicate contig names", str(raised.exception))

    def test_exact_grch38_canonical_25_contract_includes_mitochondrial_reference(self) -> None:
        contigs = grch38_canonical_25_contigs()

        summary = validate_grch38_canonical_25(contigs)

        self.assertEqual(summary.contig_count, 25)
        self.assertEqual(summary.total_reference_bases, 3_088_286_401)
        self.assertEqual(contigs[-1], ("chrM", 16569))

    def test_exact_grch38_canonical_25_rejects_every_dictionary_variant(self) -> None:
        canonical = grch38_canonical_25_contigs()
        variants = (
            canonical[:-1],
            (*canonical, ("chrUn_TEST", 100)),
            (*canonical[:-1], ("chrM", 16570)),
            (canonical[1], canonical[0], *canonical[2:]),
            tuple((name.removeprefix("chr"), length) for name, length in canonical),
        )
        for contigs in variants:
            with (
                self.subTest(contigs=contigs[-1][0]),
                self.assertRaisesRegex(ValueError, "Canonical-25"),
            ):
                validate_grch38_canonical_25(contigs)

    def test_canonical_25_run_lock_is_an_explicit_subset_of_the_pinned_grch38_lock(self) -> None:
        source_contigs = (*grch38_canonical_25_contigs(), ("GL000008.2", 209709))
        source_lock = ReferenceLock(
            reference_id="GRCh38_FULL_TEST",
            genome_build=GenomeBuild.GRCH38,
            contigs=[ReferenceContig(name=name, length=length) for name, length in source_contigs],
            allow_extra_contigs=True,
            source_fai_sha256="a" * 64,
        )

        selected = reference_lock_for_dictionary_contract(
            source_lock,
            ReferenceDictionaryContract.GRCH38_CANONICAL_25,
        )

        self.assertEqual(
            tuple((item.name, item.length) for item in selected.contigs),
            grch38_canonical_25_contigs(),
        )
        self.assertFalse(selected.allow_extra_contigs)
        self.assertEqual(selected.source_fai_sha256, source_lock.source_fai_sha256)


class ReferenceLockFromFaiTests(unittest.TestCase):
    def test_lock_records_the_index_checksum_and_every_contig(self) -> None:
        with TemporaryDirectory() as raw:
            fai = _write_fai(Path(raw), _profile(GRCH38_NUCLEAR, "chr"))
            lock = reference_lock_from_fai(
                fai,
                reference_id="synthetic-grch38-canonical",
                genome_build=GenomeBuild.GRCH38,
                require_canonical_assembly=True,
            )
        self.assertEqual(len(lock.contigs), 24)
        self.assertEqual(len(lock.source_fai_sha256 or ""), 64)

    def test_partial_reference_is_refused_when_a_named_build_is_required(self) -> None:
        with TemporaryDirectory() as raw:
            fai = _write_fai(Path(raw), [("chr5", 181538259), ("chr22", 50818468)])
            with self.assertRaises(ValueError):
                reference_lock_from_fai(
                    fai,
                    reference_id="synthetic-two-contig",
                    genome_build=GenomeBuild.GRCH38,
                    require_canonical_assembly=True,
                )

    def test_signature_is_stable_and_covers_the_extra_contig_policy(self) -> None:
        with TemporaryDirectory() as raw:
            fai = _write_fai(Path(raw), _profile(GRCH38_NUCLEAR, "chr"))
            strict = reference_lock_from_fai(
                fai,
                reference_id="synthetic-grch38-canonical",
                genome_build=GenomeBuild.GRCH38,
            )
            permissive = reference_lock_from_fai(
                fai,
                reference_id="synthetic-grch38-canonical",
                genome_build=GenomeBuild.GRCH38,
                allow_extra_contigs=True,
            )
        self.assertEqual(reference_lock_signature(strict), reference_lock_signature(strict))
        self.assertNotEqual(reference_lock_signature(strict), reference_lock_signature(permissive))


if __name__ == "__main__":
    unittest.main()
