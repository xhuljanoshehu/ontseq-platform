from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from unittest.mock import patch

import pytest

from ontseq_platform.cutesv import run_cutesv
from ontseq_platform.cutesv_build import executable_identity
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CuteSvPolicy,
    GenomeBuild,
    InputKind,
    InputSpec,
    SampleManifest,
    Verdict,
)


class CountingRunner:
    def __init__(self) -> None:
        self.calls = 0
        self.analysis_calls = 0

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        self.calls += 1
        normalized = tuple(str(item) for item in argv)
        if "--version" in normalized:
            return CommandResult(normalized, 0, "cuteSV 2.1.3\n", "")
        self.analysis_calls += 1
        return CommandResult(normalized, 0, "", "")


def _inputs(root: Path) -> tuple[SampleManifest, AlignedBamIntakeReport, Path]:
    bam = root / "synthetic.bam"
    bai = root / "synthetic.bam.bai"
    reference = root / "synthetic.fasta"
    bam.write_bytes(b"synthetic")
    bai.write_bytes(b"synthetic")
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    manifest = SampleManifest(
        sample_id="SYNTHETIC_CUTESV_QUALIFICATION",
        run_id="SYNTHETIC_CUTESV_QUALIFICATION_RUN",
        input=InputSpec(kind=InputKind.ALIGNED_BAM, path=str(bam), index_path=str(bai)),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="SYNTHETIC_REF",
        ),
        analysis=AnalysisSpec(
            profile="synthetic",
            modules=[AnalysisModule.QC, AnalysisModule.SV, AnalysisModule.REPORT],
        ),
    )
    intake = AlignedBamIntakeReport(
        sample_id=manifest.sample_id,
        reference_id=manifest.assay.reference_id,
        genome_build=manifest.assay.genome_build,
        checks=[],
        verdict=Verdict.PASS,
    )
    return manifest, intake, reference


def test_matching_version_with_unqualified_body_fails_before_any_tool_execution(
    tmp_path: Path,
) -> None:
    executable = tmp_path / "cuteSV"
    executable.write_text("#!/usr/bin/env python3\nprint('custom body')\n", encoding="utf-8")
    executable.chmod(0o700)
    identity = executable_identity(str(executable))
    assert identity["cutesv_build_id"] == "unqualified"

    manifest, intake, reference = _inputs(tmp_path)
    runner = CountingRunner()
    scratch = tmp_path / "scratch"

    with (
        patch.dict("os.environ", {"ONTSEQ_CUTESV_SCRATCH_ROOT": str(scratch)}),
        pytest.raises(ValueError, match="not qualified for the standard SV lane"),
    ):
        run_cutesv(
            manifest,
            intake,
            CuteSvPolicy(
                profile_id="synthetic-cutesv-qualification",
                status="technical_defaults_only",
                note="Synthetic qualification regression only.",
            ),
            reference_fasta=reference,
            output_vcf=tmp_path / "calls.vcf",
            runner=runner,
            cutesv=str(executable),
            expected_source_sha256=identity["cutesv_source_sha256"],
        )

    assert runner.calls == 0
    assert runner.analysis_calls == 0
    assert not (tmp_path / "calls.vcf").exists()
    assert not scratch.exists()
