from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.breakends import (
    BreakendAltForm,
    breakend_descriptors_from_sniffles_vcf,
    parse_breakend_alt,
)


class BreakendParserTests(unittest.TestCase):
    def test_all_four_vcf_breakend_forms_are_preserved(self) -> None:
        cases = {
            "N[chr2:5000[": BreakendAltForm.LOCAL_THEN_OPEN,
            "N]chr2:5000]": BreakendAltForm.LOCAL_THEN_CLOSE,
            "[chr2:5000[N": BreakendAltForm.OPEN_THEN_LOCAL,
            "]chr2:5000]N": BreakendAltForm.CLOSE_THEN_LOCAL,
        }
        for alternate, expected in cases.items():
            with self.subTest(alternate=alternate):
                chromosome, position, alt_form = parse_breakend_alt(alternate)
                self.assertEqual(chromosome, "chr2")
                self.assertEqual(position, 4999)
                self.assertEqual(alt_form, expected)

    def test_parser_rejects_non_breakend_or_ambiguous_forms(self) -> None:
        for alternate in ["<BND>", "N", "N[chr2:5000[N", "[chr2:5000["]:
            with self.subTest(alternate=alternate):
                with self.assertRaises(ValueError):
                    parse_breakend_alt(alternate)

    def test_vcf_extraction_joins_by_sniffles_record_number_without_retaining_sequence(self) -> None:
        text = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1000\tSECRET1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1200;SUPPORT=8\n"
            "chr9\t2000\tSECRET2\tN\tAC]chr22:3000]\t60\tPASS\tSVTYPE=BND;SUPPORT=7\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(text, encoding="utf-8")
            descriptors = breakend_descriptors_from_sniffles_vcf(path)

        self.assertEqual(set(descriptors), {"SNIFFLES2-000002"})
        item = descriptors["SNIFFLES2-000002"]
        self.assertEqual(item.primary_chromosome, "chr9")
        self.assertEqual(item.primary_position_0based, 1999)
        self.assertEqual(item.mate_chromosome, "chr22")
        self.assertEqual(item.mate_position_0based, 2999)
        self.assertEqual(item.alt_form, BreakendAltForm.LOCAL_THEN_CLOSE)
        serialized = item.model_dump_json()
        self.assertNotIn("SECRET2", serialized)
        self.assertNotIn("AC]chr22:3000]", serialized)
        self.assertFalse(item.inserted_sequence_retained)


if __name__ == "__main__":
    unittest.main()
