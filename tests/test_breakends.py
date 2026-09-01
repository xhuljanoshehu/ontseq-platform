from __future__ import annotations

import unittest

from ontseq_platform.breakends import BreakendParseError, resolve_breakend


class BreakendParserTests(unittest.TestCase):
    def test_all_four_vcf_breakend_alt_forms_resolve_the_mate_without_persisting_alt(self) -> None:
        alternates = (
            "AC[chr2:5000[",
            "GT]chr2:5000]",
            "[chr2:5000[CA",
            "]chr2:5000]TG",
        )
        for alternate in alternates:
            with self.subTest(alternate=alternate):
                resolved = resolve_breakend(
                    alternate,
                    declared_chromosome=None,
                    declared_position=None,
                )
                self.assertEqual(resolved.mate_chromosome, "chr2")
                self.assertEqual(resolved.mate_position_0based, 4999)
                for source_alt in alternates:
                    self.assertNotIn(source_alt, repr(resolved))

    def test_alt_is_exactly_one_ascii_allele_with_ascii_bases_and_position(self) -> None:
        invalid = (
            "<BND>",
            "N",
            "BAD[chr2:5000[",
            " N[chr2:5000[",
            "N[chr2:5000[ ",
            "N[chr 2:5000[",
            "N[chr,2:5000[",
            "N[chr2:5000[,N]chr3:6000]",
            "N[chr2:\u0665\u0660\u0660\u0660[",
            "N[chr2:5000[N",
            "[chr2:5000[",
            "N[chr2:0[",
        )
        for alternate in invalid:
            with self.subTest(alternate=alternate), self.assertRaises(BreakendParseError):
                resolve_breakend(
                    alternate,
                    declared_chromosome=None,
                    declared_position=None,
                )

    def test_official_ascii_contig_special_characters_are_accepted(self) -> None:
        resolved = resolve_breakend(
            "N[ctg!alt:part:5000[",
            declared_chromosome=None,
            declared_position=None,
        )
        self.assertEqual(resolved.mate_chromosome, "ctg!alt:part")
        self.assertEqual(resolved.mate_position_0based, 4999)

    def test_supported_symbolic_alts_require_scalar_chr2_and_end(self) -> None:
        for alternate in ("<BND>", "<TRA>"):
            with self.subTest(alternate=alternate):
                resolved = resolve_breakend(
                    alternate,
                    declared_chromosome="chr2",
                    declared_position="5000",
                )
                self.assertEqual(resolved.mate_chromosome, "chr2")
                self.assertEqual(resolved.mate_position_0based, 4999)

        for chromosome, position in ((None, "5000"), ("chr2", None), ("chr2,chr3", "5000")):
            with (
                self.subTest(chromosome=chromosome, position=position),
                self.assertRaises(BreakendParseError),
            ):
                resolve_breakend(
                    "<BND>",
                    declared_chromosome=chromosome,
                    declared_position=position,
                )

    def test_other_non_bracket_alts_do_not_use_chr2_end_fallback(self) -> None:
        for alternate in (".", "N", "<DEL>"):
            with (
                self.subTest(alternate=alternate),
                self.assertRaises(BreakendParseError) as raised,
            ):
                resolve_breakend(
                    alternate,
                    declared_chromosome="chr2",
                    declared_position="5000",
                )
            self.assertEqual(raised.exception.reason, "unsupported_breakend_alt")

    def test_bracket_alt_info_mate_must_be_equal_after_safe_normalization(self) -> None:
        resolved = resolve_breakend(
            "N[chr2:5000[",
            declared_chromosome="2",
            declared_position="05000",
        )
        self.assertEqual(resolved.mate_chromosome, "chr2")
        self.assertEqual(resolved.mate_position_0based, 4999)

        for chromosome, position in (("chr3", "5000"), ("chr2", "5001")):
            with (
                self.subTest(chromosome=chromosome, position=position),
                self.assertRaises(BreakendParseError) as raised,
            ):
                resolve_breakend(
                    "N[chr2:5000[",
                    declared_chromosome=chromosome,
                    declared_position=position,
                )
            self.assertEqual(raised.exception.reason, "conflicting_breakend_mate")

    def test_bracket_alt_validates_chr2_even_without_end(self) -> None:
        for chromosome in ("chr2", "2"):
            with self.subTest(chromosome=chromosome):
                resolved = resolve_breakend(
                    "N[chr2:5000[",
                    declared_chromosome=chromosome,
                    declared_position=None,
                )
                self.assertEqual(resolved.mate_chromosome, "chr2")
                self.assertEqual(resolved.mate_position_0based, 4999)

        with self.assertRaises(BreakendParseError) as conflicting:
            resolve_breakend(
                "N[chr2:5000[",
                declared_chromosome="chr3",
                declared_position=None,
            )
        self.assertEqual(conflicting.exception.reason, "conflicting_breakend_mate")

        for chromosome in ("chr 2", "chr2,chr3", "", True):
            with (
                self.subTest(chromosome=chromosome),
                self.assertRaises(BreakendParseError) as malformed,
            ):
                resolve_breakend(
                    "N[chr2:5000[",
                    declared_chromosome=chromosome,
                    declared_position=None,
                )
            self.assertEqual(malformed.exception.reason, "malformed_breakend_mate")

    def test_extreme_ascii_positions_fail_with_stable_parse_errors(self) -> None:
        extreme = "9" * 4301
        with self.assertRaises(BreakendParseError) as bracketed:
            resolve_breakend(
                f"N[chr2:{extreme}[",
                declared_chromosome=None,
                declared_position=None,
            )
        self.assertEqual(bracketed.exception.reason, "malformed_breakend_alt")

        with self.assertRaises(BreakendParseError) as symbolic:
            resolve_breakend(
                "<BND>",
                declared_chromosome="chr2",
                declared_position=extreme,
            )
        self.assertEqual(symbolic.exception.reason, "malformed_breakend_mate")

    def test_end_without_chr2_is_not_treated_as_a_redundant_mate_declaration(self) -> None:
        resolved = resolve_breakend(
            "N[chr2:5000[",
            declared_chromosome=None,
            declared_position="1000",
        )
        self.assertEqual(resolved.mate_chromosome, "chr2")
        self.assertEqual(resolved.mate_position_0based, 4999)


if __name__ == "__main__":
    unittest.main()
