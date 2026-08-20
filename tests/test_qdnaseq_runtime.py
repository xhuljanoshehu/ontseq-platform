from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cnv.qdnaseq import QDNAseqPolicy, run_qdnaseq_ace
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import EventType, GenomeBuild, ReferenceContig, ReferenceLock


class FakeQDNAseqRunner:
    def __init__(self, *, fail: bool = False, event_call: float = -1.0) -> None:
        self.fail = fail
        self.event_call = event_call

    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        del timeout_seconds
        args = [str(item) for item in argv]
        if self.fail:
            return CommandResult(tuple(args), 9, "", "synthetic R failure")
        out = Path(args[args.index("--output-dir") + 1])
        sample = args[args.index("--sample-id") + 1]
        out.mkdir(parents=True, exist_ok=True)
        runs = []
        for bin_size, cellularity, ploidy in (
            (100, 0.55, 2.0),
            (500, 0.57, 2.0),
            (1000, 0.56, 2.0),
        ):
            segment_name = f"{sample}.{bin_size}kbp.segments.tsv"
            chromosome_name = f"{sample}.{bin_size}kbp.chromosomes.tsv"
            fit_name = f"{sample}.{bin_size}kbp.ace-fit.png"
            copy_name = f"{sample}.{bin_size}kbp.copy-number.png"
            rds_name = f"{sample}.{bin_size}kbp.segmented.rds"
            (out / segment_name).write_text(
                "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tqnorm_log10\n"
                f"chr7\t0\t950\t2\t1.0\t{self.event_call}\t-5.0\n",
                encoding="utf-8",
            )
            (out / chromosome_name).write_text(
                "chromosome\tcopy_number\nchr7\t1.0\n",
                encoding="utf-8",
            )
            (out / fit_name).write_bytes(b"PNG-fit")
            (out / copy_name).write_bytes(b"PNG-copy")
            (out / rds_name).write_bytes(b"RDS")
            runs.append(
                {
                    "bin_size_kbp": bin_size,
                    "cellularity": cellularity,
                    "ploidy": ploidy,
                    "fit_error": 0.1 + bin_size / 10000,
                    "candidate_count": 3,
                    "alternatives": [
                        {
                            "cellularity": cellularity,
                            "ploidy": ploidy,
                            "fit_error": 0.1,
                        }
                    ],
                    "segment_file": segment_name,
                    "chromosome_file": chromosome_name,
                    "fit_plot": fit_name,
                    "copy_number_plot": copy_name,
                    "rds_file": rds_name,
                    "segment_count": 1,
                }
            )
        consensus_name = f"{sample}.consensus.chromosomes.tsv"
        (out / consensus_name).write_text(
            "chromosome\tmedian_copy_number\trounded_copy_number\tagreeing_bins\t"
            "contributing_bins\tmin_copy_number\tmax_copy_number\n"
            "chr7\t1.0\t1\t3\t3\t1.0\t1.0\n",
            encoding="utf-8",
        )
        summary = {
            "schema_version": "0.1.0",
            "sample_id": sample,
            "genome_build": "GRCh37",
            "genome_annotation": "hg19",
            "primary_bin_size_kbp": 500,
            "ace_penalty": 0.6,
            "ploidy_search": {"min": 1.5, "max": 4.5, "step": 0.05},
            "consensus_strategy": "median_rounded_across_bins",
            "package_versions": {
                "R": "4.4.0",
                "QDNAseq": "1.42.0",
                "ACE": "1.24.0",
                "DNAcopy": "1.80.0",
                "QDNAseq_hg19": "1.36.0",
                "QDNAseq_hg38": None,
            },
            "runs": runs,
            "primary": runs[1],
            "consensus_file": consensus_name,
        }
        (out / f"{sample}.qdnaseq-ace.summary.json").write_text(
            json.dumps(summary),
            encoding="utf-8",
        )
        return CommandResult(tuple(args), 0, "", "")

    def run_to_file(self, argv, output_path, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        del argv, output_path, timeout_seconds
        raise AssertionError("QDNAseq adapter should not use run_to_file")


def _lock() -> ReferenceLock:
    return ReferenceLock(
        reference_id="hg19-test",
        genome_build=GenomeBuild.GRCH37,
        contigs=[ReferenceContig(name="chr7", length=1000)],
        source_fai_sha256="a" * 64,
    )


def _policy() -> QDNAseqPolicy:
    return QDNAseqPolicy(
        profile_id="test",
        expected_qdnaseq_version="1.42.0",
        expected_ace_version="1.24.0",
        note="test policy",
    )


class QDNAseqRuntimeTests(unittest.TestCase):
    def test_promotes_complete_result_and_normalizes_event(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "sample.bam"
            bam.write_bytes(b"bam")
            script = root / "runner.R"
            script.write_text("# synthetic", encoding="utf-8")
            output = root / "cnv"

            report = run_qdnaseq_ace(
                bam=bam,
                sample_id="SAMPLE_001",
                genome_build=GenomeBuild.GRCH37,
                reference_lock=_lock(),
                policy=_policy(),
                output_dir=output,
                script=script,
                runner=FakeQDNAseqRunner(),
                threads=2,
            )

            self.assertTrue(output.is_dir())
            self.assertEqual(report.primary_fit.bin_size_kbp, 500)
            self.assertAlmostEqual(report.primary_fit.cellularity, 0.57)
            self.assertAlmostEqual(report.primary_fit.ploidy, 2.0)
            self.assertEqual(len(report.events), 1)
            self.assertEqual(report.events[0].event_type, EventType.CHROMOSOME_LOSS)
            self.assertAlmostEqual(report.events[0].copy_number or -1, 1.0)
            self.assertEqual(report.chromosome_consensus[0].agreeing_bins, 3)
            self.assertTrue(any(name.endswith(".rds") for name in report.output_files))

    def test_returns_no_call_when_ace_call_is_neutral(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "sample.bam"
            bam.write_bytes(b"bam")
            script = root / "runner.R"
            script.write_text("# synthetic", encoding="utf-8")

            report = run_qdnaseq_ace(
                bam=bam,
                sample_id="SAMPLE_002",
                genome_build=GenomeBuild.GRCH37,
                reference_lock=_lock(),
                policy=_policy(),
                output_dir=root / "cnv",
                script=script,
                runner=FakeQDNAseqRunner(event_call=0.0),
            )

            self.assertEqual(report.events, [])
            self.assertEqual(report.status.value, "NO_CALL")

    def test_failure_never_promotes_partial_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "sample.bam"
            bam.write_bytes(b"bam")
            script = root / "runner.R"
            script.write_text("# synthetic", encoding="utf-8")
            output = root / "cnv"

            with self.assertRaisesRegex(ValueError, "exited with code 9"):
                run_qdnaseq_ace(
                    bam=bam,
                    sample_id="SAMPLE_003",
                    genome_build=GenomeBuild.GRCH37,
                    reference_lock=_lock(),
                    policy=_policy(),
                    output_dir=output,
                    script=script,
                    runner=FakeQDNAseqRunner(fail=True),
                )

            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".cnv.*")), [])


if __name__ == "__main__":
    unittest.main()
