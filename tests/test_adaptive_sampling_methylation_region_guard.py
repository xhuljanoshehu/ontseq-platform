"""Regression specification for issue #87.

Adaptive Sampling methylation must not silently reuse chromosome-wide aggregation semantics.
The supported data basis is the declared target BED; incompatible direct/custom policies
must fail before any external tool is invoked.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

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


class _NoToolRunner:
    def run(self, argv, *, timeout_seconds=300):
        raise AssertionError(f"external tool must not be invoked for invalid assay/policy: {argv}")


class AdaptiveSamplingMethylationRegionGuardTests(unittest.TestCase):
    def test_chromosome_region_policy_fails_before_modkit_or_samtools(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "sample.bam"
            bam.write_bytes(b"synthetic-not-a-real-bam")
            (root / "sample.bam.bai").write_bytes(b"")
            fasta = root / "reference.fa"
            fasta.write_text(">1\nACGTCG\n", encoding="utf-8")
            bed = root / "targets.bed"
            bed.write_text("1\t0\t6\tTARGET\n", encoding="utf-8")

            # identify_modkit_binary() runs before the command runner on the current code path.
            # A tiny local executable makes this test independent of whether CI has modkit on PATH.
            modkit = root / "modkit"
            modkit.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            modkit.chmod(0o755)

            manifest = SampleManifest(
                sample_id="SYNTHETIC_AS_METH",
                run_id="SYNTHETIC_RUN",
                input=InputSpec(
                    kind=InputKind.ALIGNED_BAM,
                    path=str(bam),
                    index_path=str(bam) + ".bai",
                ),
                assay=AssaySpec(
                    mode=AssayMode.ADAPTIVE_SAMPLING,
                    genome_build=GenomeBuild.GRCH38,
                    reference_id="synthetic-grch38",
                    target_bed=str(bed),
                    target_bed_version="SYNTHETIC_V1",
                ),
                analysis=AnalysisSpec(
                    profile="synthetic", modules=[AnalysisModule.METHYLATION]
                ),
            )
            intake = AlignedBamIntakeReport(
                sample_id=manifest.sample_id,
                reference_id=manifest.assay.reference_id,
                genome_build=manifest.assay.genome_build,
                checks=[],
                verdict=Verdict.PASS,
            )
            policy = MethylationPolicy(
                profile_id="synthetic-invalid-as-chromosome",
                status="technical_defaults_only",
                cpg_only=True,
                combine_strands=True,
                region_source="chromosome",
                note="Issue #87 regression: Adaptive Sampling cannot imply genome-wide methylation.",
            )

            with self.assertRaises(ValueError) as raised:
                run_methylation(
                    manifest,
                    intake,
                    policy,
                    output_dir=root / "out",
                    reference_fasta=fasta,
                    runner=_NoToolRunner(),
                    modkit=str(modkit),
                )

        message = str(raised.exception)
        self.assertIn("Adaptive Sampling", message)
        self.assertIn("region_source=target_bed", message)
