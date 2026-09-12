from __future__ import annotations

import unittest

from ontseq_platform.bam_intake import _reference_check, parse_sam_header
from ontseq_platform.models import CheckStatus, GenomeBuild, ReferenceContig, ReferenceLock

LOCKED: list[tuple[str, int]] = [
    ("chr1", 248956422),
    ("chr2", 242193529),
    ("chrM", 16569),
]


def _lock(contigs: list[tuple[str, int]], *, allow_extra_contigs: bool = False) -> ReferenceLock:
    return ReferenceLock(
        reference_id="synthetic-three-contig",
        genome_build=GenomeBuild.GRCH38,
        contigs=[ReferenceContig(name=name, length=length) for name, length in contigs],
        allow_extra_contigs=allow_extra_contigs,
        source_fai_sha256="a" * 64,
    )


def _header(contigs: list[tuple[str, int]]) -> object:
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines.extend(f"@SQ\tSN:{name}\tLN:{length}" for name, length in contigs)
    lines.append("@RG\tID:rg1\tSM:synthetic")
    lines.append("@PG\tID:minimap2\tPN:minimap2\tVN:2.28")
    return parse_sam_header("\n".join(lines))


class SamHeaderParsingTests(unittest.TestCase):
    def test_header_fields_are_normalized(self) -> None:
        header = parse_sam_header(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:248956422\n"
            "@RG\tID:rg1\tSM:synthetic\n"
            "@RG\tID:rg2\n"
            "@PG\tID:minimap2\n"
        )
        self.assertEqual(header.sort_order, "coordinate")
        self.assertEqual(header.contigs, (("chr1", 248956422),))
        self.assertEqual(header.read_group_count, 2)
        self.assertEqual(header.sample_tag_count, 1)
        self.assertEqual(header.program_count, 1)

    def test_duplicate_sequence_names_are_rejected(self) -> None:
        with self.assertRaises(ValueError) as raised:
            parse_sam_header("@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr1\tLN:10\n")
        self.assertIn("duplicate @SQ sequence names", str(raised.exception))

    def test_non_numeric_sequence_length_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            parse_sam_header("@SQ\tSN:chr1\tLN:not-a-number\n")


class ReferenceDictionaryCheckTests(unittest.TestCase):
    def test_identical_dictionary_passes(self) -> None:
        status, message, details = _reference_check(_header(LOCKED), _lock(LOCKED))
        self.assertEqual(status, CheckStatus.PASS)
        self.assertIn("matches the reference lock", message)
        self.assertEqual(details["missing_contigs"], 0)
        self.assertEqual(
            details["expected_dictionary_sha256"], details["observed_dictionary_sha256"]
        )

    def test_missing_contig_fails_closed(self) -> None:
        status, message, details = _reference_check(_header(LOCKED[:2]), _lock(LOCKED))
        self.assertEqual(status, CheckStatus.FAIL)
        self.assertIn("Reference dictionary mismatch", message)
        self.assertEqual(details["missing_contigs"], 1)

    def test_length_mismatch_fails_closed(self) -> None:
        observed = [("chr1", 248956422), ("chr2", 1000), ("chrM", 16569)]
        status, _, details = _reference_check(_header(observed), _lock(LOCKED))
        self.assertEqual(status, CheckStatus.FAIL)
        self.assertEqual(details["length_mismatches"], 1)

    def test_extra_contig_fails_unless_the_policy_allows_it(self) -> None:
        observed = [*LOCKED, ("chrUn_GL000195v1", 182896)]
        strict, _, details = _reference_check(_header(observed), _lock(LOCKED))
        self.assertEqual(strict, CheckStatus.FAIL)
        self.assertEqual(details["extra_contigs"], 1)
        permissive, message, _ = _reference_check(
            _header(observed), _lock(LOCKED, allow_extra_contigs=True)
        )
        self.assertEqual(permissive, CheckStatus.WARN)
        self.assertIn("explicit extra-contig policy", message)

    def test_reordered_dictionary_fails_closed(self) -> None:
        observed = [LOCKED[1], LOCKED[0], LOCKED[2]]
        status, message, details = _reference_check(_header(observed), _lock(LOCKED))
        self.assertEqual(status, CheckStatus.FAIL)
        self.assertTrue(details["contig_order_mismatch"])
        self.assertIn("contig order differs", message)

    def test_diagnostics_never_copy_sequence_names(self) -> None:
        observed = [*LOCKED, ("chrUn_SAMPLE_2024_0815", 1000)]
        _, message, details = _reference_check(_header(observed), _lock(LOCKED))
        self.assertNotIn("SAMPLE_2024_0815", message)
        self.assertNotIn("SAMPLE_2024_0815", repr(details))


if __name__ == "__main__":
    unittest.main()
