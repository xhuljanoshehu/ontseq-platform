from __future__ import annotations

import unittest

from ontseq_platform.cnv.adapters import (
    ICHORCNA_MAPPING,
    SEG_MAPPING,
    ColumnMapping,
    SegmentParseError,
    call_set_from_segment_table,
    parse_segment_table,
)
from ontseq_platform.cnv.models import CnvDataBasis
from ontseq_platform.cnv.states import CopyNumberState
from ontseq_platform.models import GenomeBuild, ModuleRunStatus

SEG_TABLE = [
    "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean",
    "S1\tchr5\t70000001\t160000000\t900\t-1.0",
    "S1\tchr8\t1\t145138636\t1450\t0.585",
]

ICHOR_TABLE = [
    "chrom\tstart\tend\tseg.median.logR\tcopy.number\tCorrected_Copy_Number\tCorrected_Call",
    "5\t70000001\t160000000\t-1.0\t1\t1\tHETD",
    "8\t1\t145138636\t0.585\t3\t3\tGAIN",
]


class GenericSegTests(unittest.TestCase):
    def test_one_based_start_is_converted_to_half_open(self) -> None:
        segments, _ = parse_segment_table(SEG_TABLE, SEG_MAPPING)
        self.assertEqual(segments[0].start, 70_000_000)
        self.assertEqual(segments[0].end, 160_000_000)

    def test_state_is_derived_from_the_log_ratio(self) -> None:
        segments, warnings = parse_segment_table(SEG_TABLE, SEG_MAPPING)
        by_contig = {segment.contig: segment for segment in segments}
        self.assertEqual(by_contig["5"].state, CopyNumberState.LOSS)
        self.assertEqual(by_contig["8"].state, CopyNumberState.GAIN)
        self.assertAlmostEqual(by_contig["8"].copy_number, 3.0, places=2)
        self.assertTrue(any("derived from the reported ratio" in item for item in warnings))

    def test_supporting_bins_are_retained(self) -> None:
        segments, _ = parse_segment_table(SEG_TABLE, SEG_MAPPING)
        self.assertEqual(segments[0].supporting_bins, 900)


class IchorCnaTests(unittest.TestCase):
    def test_categorical_call_takes_precedence_over_the_ratio(self) -> None:
        segments, _ = parse_segment_table(ICHOR_TABLE, ICHORCNA_MAPPING)
        by_contig = {segment.contig: segment for segment in segments}
        self.assertEqual(by_contig["5"].state, CopyNumberState.LOSS)
        self.assertEqual(by_contig["5"].copy_number, 1.0)
        self.assertEqual(by_contig["8"].state, CopyNumberState.GAIN)

    def test_amplification_vocabulary(self) -> None:
        segments, _ = parse_segment_table(
            [
                ICHOR_TABLE[0],
                "8\t1\t1000000\t3.0\t16\t16\tHLAMP",
            ],
            ICHORCNA_MAPPING,
        )
        self.assertEqual(segments[0].state, CopyNumberState.HIGH_AMPLIFICATION)


class FailureModeTests(unittest.TestCase):
    """Header-driven parsing must fail loudly rather than shift columns silently."""

    def test_missing_required_column_raises(self) -> None:
        with self.assertRaises(SegmentParseError) as context:
            parse_segment_table(["chrom\tstart\tseg.mean", "5\t1\t0.1"], SEG_MAPPING)
        self.assertIn("missing required column", str(context.exception))

    def test_reordered_columns_are_handled_by_name(self) -> None:
        reordered = [
            "seg.mean\tloc.end\tchrom\tloc.start\tnum.mark",
            "-1.0\t160000000\tchr5\t70000001\t900",
        ]
        segments, _ = parse_segment_table(reordered, SEG_MAPPING)
        self.assertEqual(segments[0].contig, "5")
        self.assertEqual(segments[0].start, 70_000_000)

    def test_table_without_a_quantitative_column_raises(self) -> None:
        mapping = ColumnMapping(contig="chrom", start="start", end="end", log2_ratio="absent")
        with self.assertRaises(SegmentParseError):
            parse_segment_table(["chrom\tstart\tend", "5\t1\t100"], mapping)

    def test_overlapping_segments_are_rejected(self) -> None:
        table = [
            SEG_TABLE[0],
            "S1\tchr5\t1\t1000000\t10\t-1.0",
            "S1\tchr5\t500000\t2000000\t10\t0.5",
        ]
        with self.assertRaises(SegmentParseError) as context:
            parse_segment_table(table, SEG_MAPPING)
        self.assertIn("partition", str(context.exception))

    def test_inverted_interval_is_rejected(self) -> None:
        with self.assertRaises(SegmentParseError):
            parse_segment_table([SEG_TABLE[0], "S1\tchr5\t500\t100\t10\t-1.0"], SEG_MAPPING)

    def test_non_canonical_contigs_are_skipped_with_a_warning(self) -> None:
        table = [SEG_TABLE[0], SEG_TABLE[1], "S1\tchrM\t1\t16000\t5\t0.1"]
        segments, warnings = parse_segment_table(table, SEG_MAPPING)
        self.assertEqual(len(segments), 1)
        self.assertTrue(any("chrM" in item for item in warnings))

    def test_missing_header_raises(self) -> None:
        with self.assertRaises(SegmentParseError):
            parse_segment_table([], SEG_MAPPING)

    def test_na_values_do_not_become_zero(self) -> None:
        table = [SEG_TABLE[0], "S1\tchr5\t1\t1000000\tNA\tNA"]
        with self.assertRaises(SegmentParseError):
            parse_segment_table(table, SEG_MAPPING)


class CallSetTests(unittest.TestCase):
    def test_call_set_is_never_reportable(self) -> None:
        call_set = call_set_from_segment_table(
            SEG_TABLE,
            SEG_MAPPING,
            call_set_id="CS_001",
            sample_id="SYNTHETIC_AML_001",
            genome_build=GenomeBuild.GRCH38,
            method="example-caller",
            method_version="1.0.0",
            data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
            closed_world=True,
        )
        self.assertIs(call_set.reportable, False)
        self.assertIs(call_set.research_only, True)
        self.assertEqual(call_set.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(call_set.background_state, CopyNumberState.NEUTRAL)

    def test_open_world_call_set_uses_a_no_call_background(self) -> None:
        call_set = call_set_from_segment_table(
            SEG_TABLE,
            SEG_MAPPING,
            call_set_id="CS_002",
            sample_id="SYNTHETIC_AML_001",
            genome_build=GenomeBuild.GRCH38,
            method="example-caller",
            method_version="1.0.0",
            data_basis=CnvDataBasis.ADAPTIVE_SAMPLING_OFF_TARGET,
            closed_world=False,
        )
        self.assertEqual(call_set.background_state, CopyNumberState.NO_CALL)

    def test_empty_table_becomes_no_call_not_a_negative_result(self) -> None:
        call_set = call_set_from_segment_table(
            [SEG_TABLE[0]],
            SEG_MAPPING,
            call_set_id="CS_003",
            sample_id="SYNTHETIC_AML_001",
            genome_build=GenomeBuild.GRCH38,
            method="example-caller",
            method_version="1.0.0",
            data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
            closed_world=True,
        )
        self.assertEqual(call_set.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(call_set.segments, [])


if __name__ == "__main__":
    unittest.main()
