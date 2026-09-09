from __future__ import annotations

import unittest

from ontseq_platform.coordinates import (
    CoordinateSourceFormat,
    bed_to_interval,
    gff3_to_interval,
    gtf_gff3_to_interval,
    gtf_to_interval,
    normalize_contig,
    one_based_inclusive_to_interval,
    to_zero_based_half_open,
    vcf_to_interval,
)
from ontseq_platform.models import CoordinateSystem


class CoordinateNormalizationTests(unittest.TestCase):
    def test_canonical_contig_aliases_are_centralized(self) -> None:
        expected = {
            "1": "chr1",
            "chr1": "chr1",
            "X": "chrX",
            "chrY": "chrY",
            "MT": "chrM",
            "M": "chrM",
            "chrM": "chrM",
        }
        self.assertEqual(
            {contig: normalize_contig(contig) for contig in expected},
            expected,
        )

    def test_noncanonical_contig_is_preserved(self) -> None:
        self.assertEqual(normalize_contig("KI270728.1"), "KI270728.1")

    def test_gtf_and_gff3_convert_start_minus_one(self) -> None:
        gtf = gtf_to_interval("1", 101, 200)
        gff = gff3_to_interval("chr1", 101, 200)
        self.assertEqual((gtf.contig, gtf.start, gtf.end), ("chr1", 100, 200))
        self.assertEqual((gff.start, gff.end), (100, 200))
        self.assertEqual(gtf.conversion.original_start, 101)
        self.assertEqual(gtf.conversion.original_contig, "1")
        self.assertEqual(gtf.conversion.operation, "start_minus_one")

    def test_combined_gtf_gff_api_rejects_other_formats(self) -> None:
        with self.assertRaisesRegex(ValueError, "gtf or gff3"):
            gtf_gff3_to_interval(
                "chr1",
                1,
                1,
                source_format=CoordinateSourceFormat.BED,
            )

    def test_vcf_position_becomes_single_base_interval(self) -> None:
        interval = vcf_to_interval("22", 1001)
        self.assertEqual((interval.contig, interval.start, interval.end), ("chr22", 1000, 1001))
        self.assertEqual(
            interval.conversion.source_coordinate_system,
            CoordinateSystem.ONE_BASED_POSITION,
        )

    def test_confirmed_standard_bed_is_unchanged(self) -> None:
        interval = bed_to_interval("chr2", 10, 20, confirmed_standard_bed=True)
        self.assertEqual((interval.start, interval.end), (10, 20))
        self.assertEqual(interval.conversion.operation, "unchanged")

    def test_bed_requires_explicit_standard_confirmation(self) -> None:
        with self.assertRaisesRegex(ValueError, "confirmed"):
            bed_to_interval("chr2", 10, 20, confirmed_standard_bed=False)

    def test_one_based_inclusive_panel_region_preserves_end(self) -> None:
        interval = one_based_inclusive_to_interval("5", 100, 200)
        self.assertEqual((interval.start, interval.end), (99, 200))

    def test_invalid_and_empty_intervals_fail(self) -> None:
        with self.assertRaises(ValueError):
            to_zero_based_half_open(0, 1, CoordinateSystem.ONE_BASED_INCLUSIVE)
        with self.assertRaises(ValueError):
            to_zero_based_half_open(10, 10, CoordinateSystem.ZERO_BASED_HALF_OPEN)
        with self.assertRaises(TypeError):
            to_zero_based_half_open(True, 2, CoordinateSystem.ZERO_BASED_HALF_OPEN)


if __name__ == "__main__":
    unittest.main()
