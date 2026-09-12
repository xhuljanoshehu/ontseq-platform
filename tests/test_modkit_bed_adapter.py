"""Protect the modkit include-BED handoff independently of the real-binary lane."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.execution import CommandResult
from ontseq_platform.methylation import MethylationPolicy, run_methylation
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    SampleManifest,
    Verdict,
)
from ontseq_platform.reference import sha256_file


class IncludeBedTests(unittest.TestCase):
    def test_labels_do_not_reach_modkit_as_unsupported_four_or_five_column_bed(self) -> None:
        for tail, expected in [
            ("", "1\t0\t30\n"),
            ("\tGENE", "1\t0\t30\n"),
            ("\tGENE\t100", "1\t0\t30\n"),
            ("\tGENE\t100\t-", "1\t0\t30\tGENE\t100\t-\n"),
        ]:
            with self.subTest(tail=tail), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                bed = root / "original.bed"
                bed.write_text("# synthetic\n1\t0\t30" + tail + "\n", encoding="utf-8")
                original = bed.read_bytes()
                fasta = root / "reference.fa"
                fasta.write_text(">1\nACGT\n", encoding="utf-8")

                class Runner:
                    def __init__(self):
                        self.calls: list[bytes] = []

                    def run(self, argv, *, timeout_seconds=300):
                        if argv[1] == "--version":
                            return CommandResult(tuple(argv), 0, "mod_kit 0.6.4", "")
                        if argv[1] == "view":
                            # This fixture is about the include-BED projection. It contains
                            # one MM-tagged read but no independent same-base cytosine MM
                            # groups, so the two samtools safety questions must not share
                            # one canned answer.
                            if any("C[+-]" in str(item) for item in argv):
                                return CommandResult(tuple(argv), 0, "0", "")
                            return CommandResult(tuple(argv), 0, "1", "")
                        include = Path(argv[argv.index("--include-bed") + 1])
                        self.calls.append(include.read_bytes())
                        # Known-confidence synthetic pileup: 15 modified out of 20 calls.
                        Path(argv[3]).write_text(
                            "1\t1\t2\tm\t20\t+\t1\t2\t0,0,0\t20\t75\t15\t5\t0\t0\t0\t0\t0\n",
                            encoding="utf-8",
                        )
                        return CommandResult(tuple(argv), 0, "", "")

                manifest = SampleManifest(
                    sample_id="SYNTHETIC_BED",
                    run_id="SYNTHETIC_RUN",
                    input=InputSpec(
                        kind=InputKind.ALIGNED_BAM,
                        path=str(root / "synthetic.bam"),
                        index_path=str(root / "synthetic.bam.bai"),
                    ),
                    assay=AssaySpec(
                        mode=AssayMode.ADAPTIVE_SAMPLING,
                        genome_build=GenomeBuild.GRCH38,
                        reference_id="synthetic",
                        target_bed=str(bed),
                        target_bed_version="SYNTHETIC_V1",
                    ),
                    analysis=AnalysisSpec(
                        profile="synthetic", modules=[AnalysisModule.METHYLATION]
                    ),
                )
                intake = AlignedBamIntakeReport(
                    sample_id=manifest.sample_id,
                    reference_id="synthetic",
                    genome_build=GenomeBuild.GRCH38,
                    checks=[],
                    verdict=Verdict.PASS,
                )
                policy = MethylationPolicy(
                    profile_id="synthetic",
                    status="technical_defaults_only",
                    cpg_only=False,
                    combine_strands=False,
                    region_source="target_bed",
                    note="Synthetic include-BED handoff regression",
                )
                runner = Runner()
                report = run_methylation(
                    manifest,
                    intake,
                    policy,
                    output_dir=root / "out",
                    reference_fasta=fasta,
                    runner=runner,
                )
                self.assertEqual(runner.calls, [expected.encode()])
                self.assertEqual(bed.read_bytes(), original)
                self.assertEqual(report.target_bed_fingerprint.sha256, sha256_file(bed))
                self.assertEqual(report.regions[0].region_id, "GENE" if tail else "1:0-30")
                self.assertEqual(report.regions[0].mean_modified_fraction, 0.75)
                self.assertEqual(report.tool.parameters["include_bed_format"], "bed3-or-bed6-v1")
