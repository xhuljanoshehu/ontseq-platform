from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from test_qdnaseq_runtime import FakeQDNAseqRunner, _lock, _policy

from ontseq_platform.cnv.qdnaseq import (
    CnvFit,
    _assessable_bin_extents,
    _events_from_primary_segments,
    run_qdnaseq_ace,
)
from ontseq_platform.models import (
    EventType,
    GenomeBuild,
    GenomicEvent,
    ReferenceContig,
    ReferenceLock,
)

ROOT = Path(__file__).resolve().parents[1]


def _fit() -> CnvFit:
    return CnvFit(
        bin_size_kbp=500,
        cellularity=1.0,
        ploidy=2.0,
        fit_error=0.0,
        candidate_count=1,
        segment_count=4,
        segment_file="segments.tsv",
        chromosome_file="chromosomes.tsv",
        fit_plot="fit.png",
        copy_number_plot="copy.png",
        rds_file="native.rds",
    )


class QDNAseqCoordinateTests(unittest.TestCase):
    def test_summary_requires_explicit_corrected_coordinate_contract(self) -> None:
        class WrongCoordinates(FakeQDNAseqRunner):
            def run(self, argv, *, timeout_seconds=300):
                result = super().run(argv, timeout_seconds=timeout_seconds)
                args = [str(item) for item in argv]
                path = Path(args[args.index("--output-dir") + 1]) / "SYNTH.qdnaseq-ace.summary.json"
                summary = json.loads(path.read_text())
                summary.pop("coordinate_system")
                path.write_text(json.dumps(summary))
                return result

        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            bam, script = root / "synthetic.bam", root / "synthetic.R"
            bam.write_bytes(b"synthetic")
            script.write_text("# synthetic")
            with self.assertRaisesRegex(ValueError, "explicit zero_based_half_open"):
                run_qdnaseq_ace(
                    bam=bam,
                    sample_id="SYNTH",
                    genome_build=GenomeBuild.GRCH37,
                    reference_lock=_lock(),
                    policy=_policy(),
                    output_dir=root / "result",
                    script=script,
                    runner=WrongCoordinates(),
                )
            self.assertFalse((root / "result").exists())

    def test_segment_bounds_and_coordinate_tags_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "segments.tsv"
            for start, end, tag in (
                (-1, 500, "zero_based_half_open"),
                (0, 1001, "zero_based_half_open"),
                (1, 500, "one_based_inclusive"),
                (0, 500, ""),
            ):
                with self.subTest(start=start, end=end, tag=tag):
                    path.write_text(
                        "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tcoordinate_system\n"
                        f"chr7\t{start}\t{end}\t1\t1\t-1\t{tag}\n"
                    )
                    with self.assertRaises(ValueError):
                        _events_from_primary_segments(
                            path,
                            sample_id="SYNTH",
                            fit=_fit(),
                            tools=[],
                            reference_lock=_lock(),
                            minimum_segment_bins=1,
                            whole_chromosome_fraction=0.9,
                            consensus={},
                        )

    def test_fractional_integer_fields_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "segments.tsv"
            path.write_text(
                "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tcoordinate_system\n"
                "chr7\t0\t500\t1.5\t1\t-1\tzero_based_half_open\n"
            )
            with self.assertRaisesRegex(ValueError, "invalid integer value for bin_count"):
                _events_from_primary_segments(
                    path,
                    sample_id="SYNTH",
                    fit=_fit(),
                    tools=[],
                    reference_lock=_lock(),
                    minimum_segment_bins=1,
                    whole_chromosome_fraction=0.9,
                    consensus={},
                )

    def test_negative_infinite_qnorm_is_preserved_as_unbounded_not_finite_quality(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "segments.tsv"
            path.write_text(
                "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tqnorm_log10\tcoordinate_system\n"
                "chr7\t0\t500\t1\t1\t-1\t-Inf\tzero_based_half_open\n"
            )
            events, warnings = _events_from_primary_segments(
                path,
                sample_id="SYNTH",
                fit=_fit(),
                tools=[],
                reference_lock=_lock(),
                minimum_segment_bins=1,
                whole_chromosome_fraction=0.9,
                consensus={},
            )
            self.assertFalse(warnings)
            self.assertEqual(len(events), 1)
            self.assertIsNone(events[0].evidence[0].quality)
            self.assertTrue(any("qnorm_log10=-Inf" in note for note in events[0].notes))

    def test_non_finite_numeric_fields_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "segments.tsv"
            for field, absolute_copy_number, call in (
                ("absolute_copy_number", "inf", "-1"),
                ("call", "1", "-inf"),
            ):
                with self.subTest(field=field):
                    path.write_text(
                        "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tcoordinate_system\n"
                        f"chr7\t0\t500\t1\t{absolute_copy_number}\t{call}\tzero_based_half_open\n"
                    )
                    with self.assertRaisesRegex(ValueError, "non-finite numeric value"):
                        _events_from_primary_segments(
                            path,
                            sample_id="SYNTH",
                            fit=_fit(),
                            tools=[],
                            reference_lock=_lock(),
                            minimum_segment_bins=1,
                            whole_chromosome_fraction=0.9,
                            consensus={},
                        )

    @unittest.skipUnless(
        shutil.which("Rscript"), "Rscript is required for actual R exporter regression"
    )
    def test_actual_r_export_preserves_first_bin_terminal_base_and_chromosome_end(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            output = Path(raw)
            completed = subprocess.run(
                [
                    "Rscript",
                    str(ROOT / "tests/qdnaseq_coordinate_export.R"),
                    str(ROOT / "scripts/run_qdnaseq_ace.R"),
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            self.assertIn("regression PASS", completed.stdout)
            lock = ReferenceLock(
                reference_id="SYNTHETIC_COORDINATES",
                genome_build=GenomeBuild.GRCH37,
                contigs=[
                    ReferenceContig(name="chr7", length=1000),
                    ReferenceContig(name="chr8", length=1000),
                ],
                source_fai_sha256="a" * 64,
            )
            events, warnings = _events_from_primary_segments(
                output / "segments.tsv",
                sample_id="SYNTH",
                fit=_fit(),
                tools=[],
                reference_lock=lock,
                minimum_segment_bins=1,
                whole_chromosome_fraction=0.9,
                consensus={},
            )
            self.assertFalse(warnings)
            self.assertEqual(
                [(e.primary.start, e.primary.end, e.length_bp) for e in events],
                [(0, 500, 500), (500, 999, 499), (999, 1000, 1), (0, 1000, 1000)],
            )
            self.assertEqual(events[-1].event_type, EventType.CHROMOSOME_GAIN)
            self.assertTrue(all(not event.reportable for event in events))


class WholeChromosomeSpanBasisTests(unittest.TestCase):
    """QDNAseq filters telomeric and blacklisted bins, so a real whole-chromosome segment
    never reaches the raw contig ends. The versioned span basis decides which evidence a
    +chr/-chr fragment requires."""

    def _lock_chr8(self) -> ReferenceLock:
        return ReferenceLock(
            reference_id="SYNTHETIC_SPAN",
            genome_build=GenomeBuild.GRCH37,
            contigs=[ReferenceContig(name="chr8", length=146_364_022)],
            source_fai_sha256="b" * 64,
        )

    def _segments(self, root: Path, start: int) -> Path:
        path = root / "segments.tsv"
        path.write_text(
            "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tcoordinate_system\n"
            f"chr8\t{start}\t146364022\t243\t9\t3\tzero_based_half_open\n",
            encoding="utf-8",
        )
        return path

    def _events(
        self,
        path: Path,
        *,
        basis: str = "exact_contig",
        extents: dict[str, tuple[int, int]] | None = None,
    ) -> list[GenomicEvent]:
        events, _ = _events_from_primary_segments(
            path,
            sample_id="SYNTH",
            fit=_fit(),
            tools=[],
            reference_lock=self._lock_chr8(),
            minimum_segment_bins=1,
            whole_chromosome_fraction=0.9,
            consensus={},
            whole_chromosome_span_basis=basis,
            assessable_extents=extents,
        )
        return events

    def test_assessable_extent_reads_the_exporter_use_flag(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            bins = Path(raw) / "bins.tsv"
            bins.write_text(
                "chromosome\tstart\tend\tuse\tcoordinate_system\n"
                "8\t0\t500000\tFALSE\tzero_based_half_open\n"
                "8\t500000\t1000000\tTRUE\tzero_based_half_open\n"
                "8\t145864022\t146364022\tTRUE\tzero_based_half_open\n",
                encoding="utf-8",
            )
            self.assertEqual(_assessable_bin_extents(bins), {"chr8": (500_000, 146_364_022)})

    def test_exporter_without_use_flag_yields_no_extent(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            bins = Path(raw) / "bins.tsv"
            bins.write_text(
                "chromosome\tstart\tend\tcoordinate_system\n"
                "8\t500000\t1000000\tzero_based_half_open\n",
                encoding="utf-8",
            )
            self.assertEqual(_assessable_bin_extents(bins), {})

    def test_filtered_first_bin_suppresses_the_span_under_the_exact_contig_basis(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            events = self._events(self._segments(Path(raw), 500_000))
            self.assertEqual(events[0].event_type, EventType.CHROMOSOME_GAIN)
            self.assertFalse(events[0].whole_chromosome_span_confirmed)

    def test_full_assessable_coverage_confirms_the_span(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            events = self._events(
                self._segments(Path(raw), 500_000),
                basis="assessable_bin_extent",
                extents={"chr8": (500_000, 146_364_022)},
            )
            self.assertEqual(events[0].event_type, EventType.CHROMOSOME_GAIN)
            self.assertTrue(events[0].whole_chromosome_span_confirmed)
            self.assertTrue(
                any("assessable QDNAseq bin extent" in note for note in events[0].notes)
            )

    def test_partial_assessable_coverage_still_suppresses_the_span(self) -> None:
        # 5 Mbp missing at the p-arm end stays above the 0.9 classification fraction, so the
        # event is still a chromosome gain, but it no longer covers the assessable extent.
        with tempfile.TemporaryDirectory() as raw:
            events = self._events(
                self._segments(Path(raw), 5_000_000),
                basis="assessable_bin_extent",
                extents={"chr8": (500_000, 146_364_022)},
            )
            self.assertEqual(events[0].event_type, EventType.CHROMOSOME_GAIN)
            self.assertFalse(events[0].whole_chromosome_span_confirmed)
            self.assertTrue(
                any("ISCN rendering is therefore suppressed" in note for note in events[0].notes)
            )

    def test_full_assessable_extent_below_fraction_stays_a_segment(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "segments.tsv"
            for copy_number, call, expected in (
                (3, 1, EventType.DUPLICATION),
                (1, -1, EventType.DELETION),
            ):
                with self.subTest(event_type=expected):
                    path.write_text(
                        "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall"
                        "\tcoordinate_system\n"
                        f"chr8\t30000000\t146364022\t200\t{copy_number}\t{call}"
                        "\tzero_based_half_open\n",
                        encoding="utf-8",
                    )
                    events = self._events(
                        path,
                        basis="assessable_bin_extent",
                        extents={"chr8": (30_000_000, 146_364_022)},
                    )
                    self.assertEqual(len(events), 1)
                    event = events[0]
                    self.assertEqual(event.event_type, expected)
                    self.assertFalse(event.whole_chromosome_span_confirmed)
                    self.assertIsNone(event.assessable_span_start)
                    self.assertIsNone(event.assessable_span_end)
                    self.assertFalse(event.reportable)
                    self.assertFalse(any("span confirmed" in note for note in event.notes))


if __name__ == "__main__":
    unittest.main()
