from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from test_qdnaseq_runtime import FakeQDNAseqRunner, _lock, _policy

from ontseq_platform.cnv.qdnaseq import CnvFit, _events_from_primary_segments, run_qdnaseq_ace
from ontseq_platform.models import EventType, GenomeBuild, ReferenceContig, ReferenceLock

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


if __name__ == "__main__":
    unittest.main()
