from __future__ import annotations

import unittest

from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CheckStatus,
    CraminoQCReport,
    GenomeBuild,
    InputKind,
    InputSpec,
    ModuleRunStatus,
    QCMetrics,
    SampleManifest,
    ToolRecord,
    ValidationCheck,
    Verdict,
)
from ontseq_platform.mvp import assemble_aligned_bam_mvp


def _manifest() -> SampleManifest:
    return SampleManifest(
        sample_id="SYNTHETIC_001",
        run_id="SYNTHETIC_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path="/secure/SYNTHETIC_001.bam",
            index_path="/secure/SYNTHETIC_001.bam.bai",
        ),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="SYNTHETIC_REF",
        ),
        analysis=AnalysisSpec(
            profile="lcwgs",
            modules=[
                AnalysisModule.QC,
                AnalysisModule.CNV,
                AnalysisModule.SV,
                AnalysisModule.ISCN,
                AnalysisModule.REPORT,
            ],
        ),
    )


class AlignedBamMVPTests(unittest.TestCase):
    def test_unrun_callers_are_explicit_not_biological_negatives(self) -> None:
        manifest = _manifest()
        samtools = ToolRecord(name="samtools", version="1.24")
        cramino = ToolRecord(name="cramino", version="1.3.0")
        intake = AlignedBamIntakeReport(
            sample_id=manifest.sample_id,
            reference_id=manifest.assay.reference_id,
            genome_build=manifest.assay.genome_build,
            checks=[
                ValidationCheck(name="synthetic", status=CheckStatus.PASS, message="synthetic")
            ],
            verdict=Verdict.PASS,
            tool=samtools,
        )
        qc = CraminoQCReport(
            sample_id=manifest.sample_id,
            qc=QCMetrics(
                verdict=Verdict.WARN,
                metrics={"mean_coverage_x": 3.0},
                warnings=["No validated numeric gates"],
            ),
            tool=cramino,
        )
        result = assemble_aligned_bam_mvp(
            manifest,
            intake,
            qc,
            pipeline_version="0.2.0-dev",
            git_commit="SYNTHETIC",
        )

        self.assertEqual(result.events, [])
        self.assertEqual(result.iscn.notation, "NOT GENERATED")
        status = {item.module: item.status for item in result.modules}
        self.assertEqual(status[AnalysisModule.QC], ModuleRunStatus.COMPLETED)
        self.assertEqual(status[AnalysisModule.CNV], ModuleRunStatus.NOT_RUN)
        self.assertEqual(status[AnalysisModule.SV], ModuleRunStatus.NOT_RUN)
        self.assertTrue(any("not a biological negative" in item for item in result.iscn.warnings))

    def test_failed_intake_cannot_be_assembled(self) -> None:
        manifest = _manifest()
        intake = AlignedBamIntakeReport(
            sample_id=manifest.sample_id,
            reference_id=manifest.assay.reference_id,
            genome_build=manifest.assay.genome_build,
            checks=[
                ValidationCheck(name="synthetic", status=CheckStatus.FAIL, message="synthetic")
            ],
            verdict=Verdict.FAIL,
        )
        qc = CraminoQCReport(
            sample_id=manifest.sample_id,
            qc=QCMetrics(verdict=Verdict.WARN, metrics={}),
            tool=ToolRecord(name="cramino", version="synthetic"),
        )
        with self.assertRaises(ValueError):
            assemble_aligned_bam_mvp(
                manifest,
                intake,
                qc,
                pipeline_version="0.2.0-dev",
                git_commit="SYNTHETIC",
            )


if __name__ == "__main__":
    unittest.main()
