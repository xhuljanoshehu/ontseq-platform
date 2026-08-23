from __future__ import annotations

import unittest

import tempfile
from pathlib import Path

from ontseq_platform.cnv.adapters import (
    ICHORCNA_MAPPING,
    QDNASEQ_ACE_MAPPING,
    SEG_MAPPING,
    ColumnMapping,
    SegmentParseError,
    call_set_from_qdnaseq_report,
    call_set_from_segment_table,
    parse_segment_table,
    qdnaseq_method_version,
)
from ontseq_platform.cnv.models import CnvDataBasis
from ontseq_platform.cnv.qdnaseq import CnvFit, QDNAseqCallReport
from ontseq_platform.cnv.states import CopyNumberState
from ontseq_platform.models import (
    EventType,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    ToolRecord,
)

SEG_TABLE = [
    "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean",
    "S1\tchr5\t70000001\t160000000\t900\t-1.0",
    "S1\tchr8\t1\t145138636\t1450\t0.585",
]

# The layout scripts/run_qdnaseq_ace.R writes: one-based inclusive, ACE-adjusted
# absolute copy number, and one row per collapsed run of bins.
QDNASEQ_TABLE = [
    "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tqnorm_log10",
    "chr7\t1\t60000000\t120\t1.04\t1\t-0.42",
    "chr7\t60000001\t159345973\t199\t2.02\t2\t0.01",
    "chr8\t1\t145138636\t290\t3.11\t3\t0.36",
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


def _fit(bin_size_kbp: int, *, segment_file: str) -> CnvFit:
    return CnvFit(
        bin_size_kbp=bin_size_kbp,
        cellularity=0.62,
        ploidy=2.1,
        fit_error=0.031,
        candidate_count=4,
        segment_count=3,
        segment_file=segment_file,
        chromosome_file=f"S.{bin_size_kbp}kbp.chromosomes.tsv",
        fit_plot=f"S.{bin_size_kbp}kbp.fit.png",
        copy_number_plot=f"S.{bin_size_kbp}kbp.cn.png",
        rds_file=f"S.{bin_size_kbp}kbp.segmented.rds",
    )


#: One promoted alteration, for the case where the runtime lane found something.
LOSS_EVENT = GenomicEvent(
    event_id="CNV_CHR7_LOSS",
    event_type=EventType.CHROMOSOME_LOSS,
    primary=Locus(chromosome="chr7", start=0, end=159_138_663),
    copy_number=1.04,
)


def _report(**overrides: object) -> QDNAseqCallReport:
    """A report whose segment table is a full partition but promotes no alteration.

    ``QDNAseqCallReport.status`` is derived from its *events*, so a run that segmented the
    whole genome and found nothing worth promoting is ``NO_CALL`` there. That is a different
    question from the one the call set answers, and the tests below rely on the difference.
    """
    primary = _fit(500, segment_file="S.500kbp.segments.tsv")
    defaults: dict[str, object] = {
        "sample_id": "SYNTHETIC_AML_001",
        "genome_build": GenomeBuild.GRCH38,
        "status": ModuleRunStatus.NO_CALL,
        "primary_fit": primary,
        "fits": [primary, _fit(1000, segment_file="S.1000kbp.segments.tsv")],
        "chromosome_consensus": [],
        "events": [],
        "tools": [
            ToolRecord(name="QDNAseq", version="1.38.0"),
            ToolRecord(name="ACE", version="1.20.0"),
            ToolRecord(name="R", version="4.4.1"),
        ],
        "output_files": ["S.500kbp.segments.tsv"],
    }
    defaults.update(overrides)
    return QDNAseqCallReport(**defaults)  # type: ignore[arg-type]


class QDNAseqMappingTests(unittest.TestCase):
    def test_one_based_start_is_converted_to_half_open(self) -> None:
        segments, _ = parse_segment_table(QDNASEQ_TABLE, QDNASEQ_ACE_MAPPING)
        self.assertEqual(segments[0].start, 0)
        self.assertEqual(segments[0].end, 60_000_000)
        self.assertEqual(segments[1].start, 60_000_000)

    def test_absolute_copy_number_drives_the_state_not_the_rounded_call(self) -> None:
        segments, _ = parse_segment_table(QDNASEQ_TABLE, QDNASEQ_ACE_MAPPING)
        by_key = {(segment.contig, segment.start): segment for segment in segments}
        self.assertEqual(by_key[("7", 0)].state, CopyNumberState.LOSS)
        self.assertEqual(by_key[("7", 60_000_000)].state, CopyNumberState.NEUTRAL)
        self.assertEqual(by_key[("8", 0)].state, CopyNumberState.GAIN)
        self.assertAlmostEqual(by_key[("8", 0)].copy_number or 0.0, 3.11)

    def test_bin_count_is_carried_as_supporting_bins(self) -> None:
        segments, _ = parse_segment_table(QDNASEQ_TABLE, QDNASEQ_ACE_MAPPING)
        self.assertEqual(segments[0].supporting_bins, 120)


class QDNAseqMethodVersionTests(unittest.TestCase):
    def test_both_packages_are_named(self) -> None:
        version = qdnaseq_method_version(
            [ToolRecord(name="QDNAseq", version="1.38.0"), ToolRecord(name="ACE", version="1.20.0")]
        )
        self.assertEqual(version, "QDNAseq 1.38.0+ACE 1.20.0")

    def test_a_missing_version_refuses_rather_than_guessing(self) -> None:
        with self.assertRaises(SegmentParseError) as raised:
            qdnaseq_method_version([ToolRecord(name="QDNAseq", version="1.38.0")])
        self.assertIn("ACE", str(raised.exception))


class QDNAseqCallSetTests(unittest.TestCase):
    def _write(self, directory: Path, name: str = "S.500kbp.segments.tsv") -> None:
        (directory / name).write_text("\n".join(QDNASEQ_TABLE) + "\n", encoding="utf-8")

    def test_provenance_and_fit_reach_the_call_set(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory)
            call_set = call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_001",
                data_basis=CnvDataBasis.ADAPTIVE_SAMPLING_OFF_TARGET,
                output_dir=directory,
            )
        self.assertEqual(call_set.method, "QDNAseq+ACE")
        self.assertEqual(call_set.method_version, "QDNAseq 1.38.0+ACE 1.20.0")
        self.assertEqual(call_set.data_basis, CnvDataBasis.ADAPTIVE_SAMPLING_OFF_TARGET)
        self.assertEqual(call_set.bin_size_bp, 500_000)
        self.assertAlmostEqual(call_set.estimated_tumor_fraction or 0.0, 0.62)
        self.assertAlmostEqual(call_set.estimated_ploidy or 0.0, 2.1)
        self.assertEqual(call_set.status, ModuleRunStatus.COMPLETED)

    def test_the_reports_status_does_not_decide_the_call_sets_status(self) -> None:
        """Two different questions, deliberately not collapsed into one.

        The runtime report says whether an alteration was promoted to an event. The call
        set says whether the method produced a scoreable partition. A run that segmented
        the whole genome and promoted nothing is NO_CALL in the first sense and COMPLETED
        in the second, and reading one as the other would turn "found nothing" into "could
        not look" — the exact confusion the vocabulary exists to prevent.
        """
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory)
            quiet = call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_010",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=directory,
            )
            loud = call_set_from_qdnaseq_report(
                _report(status=ModuleRunStatus.COMPLETED, events=[LOSS_EVENT]),
                call_set_id="QDNASEQ_011",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=directory,
            )
        self.assertEqual(quiet.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(loud.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(len(quiet.segments), len(loud.segments))

    def test_nothing_here_can_make_the_lane_reportable(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory)
            call_set = call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_002",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=directory,
            )
        self.assertIs(call_set.reportable, False)
        self.assertIs(call_set.research_only, True)

    def test_uncovered_contigs_become_declared_no_calls(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory)
            call_set = call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_003",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=directory,
                contig_lengths={"chr7": 159_345_973, "chr8": 145_138_636, "chrX": 156_040_895},
            )
        no_call_contigs = {region.contig for region in call_set.no_call_regions}
        # chr7 and chr8 are fully covered by the table; chrX is not covered at all and must
        # not be scored as agreement with whatever the truth set asserts there.
        self.assertEqual(no_call_contigs, {"X"})
        region = next(item for item in call_set.no_call_regions if item.contig == "X")
        self.assertEqual((region.start, region.end), (0, 156_040_895))

    def test_a_gap_inside_a_covered_contig_is_also_a_no_call(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory)
            call_set = call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_004",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=directory,
                contig_lengths={"chr7": 200_000_000},
            )
        regions = [item for item in call_set.no_call_regions if item.contig == "7"]
        self.assertEqual([(item.start, item.end) for item in regions], [(159_345_973, 200_000_000)])

    def test_calls_and_no_calls_account_for_each_contig_exactly(self) -> None:
        """The property that makes the denominator auditable, asserted rather than assumed."""
        lengths = {"chr7": 200_000_000, "chr8": 145_138_636, "chrX": 156_040_895}
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory)
            call_set = call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_009",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=directory,
                contig_lengths=lengths,
            )
        spans: dict[str, list[tuple[int, int]]] = {}
        for segment in call_set.segments:
            spans.setdefault(segment.contig, []).append((segment.start, segment.end))
        for region in call_set.no_call_regions:
            spans.setdefault(region.contig, []).append((region.start, region.end))
        for contig, length in lengths.items():
            items = sorted(spans.get(contig.removeprefix("chr"), []))
            previous_end = 0
            covered = 0
            for start, end in items:
                self.assertGreaterEqual(start, previous_end, f"{contig} spans overlap")
                covered += end - start
                previous_end = end
            self.assertEqual(covered, length, f"{contig} does not reconcile")

    def test_without_contig_lengths_the_gap_is_stated_rather_than_hidden(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory)
            call_set = call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_005",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=directory,
            )
        self.assertEqual(call_set.no_call_regions, [])
        self.assertTrue(any("could not be declared" in item for item in call_set.warnings))

    def test_a_requested_bin_size_selects_that_fit(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory, "S.1000kbp.segments.tsv")
            call_set = call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_006",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=directory,
                bin_size_kbp=1000,
            )
        self.assertEqual(call_set.bin_size_bp, 1_000_000)

    def test_an_absent_bin_size_refuses_and_names_what_exists(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            self._write(directory)
            with self.assertRaises(SegmentParseError) as raised:
                call_set_from_qdnaseq_report(
                    _report(),
                    call_set_id="QDNASEQ_007",
                    data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                    output_dir=directory,
                    bin_size_kbp=250,
                )
        self.assertIn("500", str(raised.exception))
        self.assertIn("1000", str(raised.exception))

    def test_a_missing_segment_table_refuses_rather_than_scoring_nothing(self) -> None:
        with (
            tempfile.TemporaryDirectory() as raw,
            self.assertRaises(SegmentParseError),
        ):
            call_set_from_qdnaseq_report(
                _report(),
                call_set_id="QDNASEQ_008",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                output_dir=Path(raw),
            )


if __name__ == "__main__":
    unittest.main()
