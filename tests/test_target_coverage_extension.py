from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path
from unittest import mock

from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    FileFingerprint,
    GenomeBuild,
    InputKind,
    InputSpec,
    ModuleRunStatus,
    QCPolicy,
    ReferenceContig,
    ReferenceLock,
    SampleManifest,
    ToolRecord,
)
from ontseq_platform.pipeline.envelope import RunEnvelope
from ontseq_platform.pipeline.runner import IMPLEMENTATIONS, RunConfiguration, RunContext
from ontseq_platform.pipeline.stages import SPEC_BY_STAGE, StageId, VerificationStatus
from ontseq_platform.target_coverage import (
    TargetCoveragePolicy,
    TargetCoverageRegion,
    TargetCoverageReport,
)
from ontseq_platform.target_coverage_extension import (
    TARGET_COVERAGE_REPORT,
    TARGET_COVERAGE_WORK,
    TargetCoverageExtensionSettings,
    register_target_coverage_extension,
)


class _VersionRunner:
    def __init__(self, version: str = "0.3.14") -> None:
        self.version = version
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        del timeout_seconds
        normalized = tuple(str(item) for item in argv)
        self.calls.append(normalized)
        return CommandResult(
            argv=normalized,
            returncode=0,
            stdout=f"mosdepth {self.version}\n",
            stderr="",
        )

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        del output_path, timeout_seconds
        raise AssertionError(f"unexpected streamed command: {list(argv)}")


class TargetCoverageExtensionTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.base = Path(self._temporary.name)
        self.bam = self.base / "sample.bam"
        self.bam.write_bytes(b"synthetic-bam-placeholder")
        self.bed = self.base / "targets.bed"
        self.bed.write_text("chr1\t10\t110\tROI_A\n", encoding="utf-8")
        self.policy = TargetCoveragePolicy(
            profile_id="adaptive_target_coverage_technical_v1",
            status="technical_defaults_only",
            expected_version="0.3.14",
            thresholds=[1, 10, 20, 30],
            mapq=0,
            exclude_flags=1796,
            note="Technical descriptive defaults only; not an adequacy threshold.",
        )
        register_target_coverage_extension(TargetCoverageExtensionSettings(policy=self.policy))

    def _manifest(self, mode: AssayMode) -> SampleManifest:
        assay = AssaySpec(
            mode=mode,
            genome_build=GenomeBuild.GRCH38,
            reference_id="TEST_REF",
            target_bed=str(self.bed) if mode == AssayMode.ADAPTIVE_SAMPLING else None,
            target_bed_version="synthetic-v1" if mode == AssayMode.ADAPTIVE_SAMPLING else None,
        )
        return SampleManifest(
            sample_id="SYNTHETIC_001",
            run_id="RUN_SYNTHETIC_001",
            input=InputSpec(
                kind=InputKind.ALIGNED_BAM,
                path=str(self.bam),
                index_path=str(self.base / "sample.bam.bai"),
            ),
            assay=assay,
            analysis=AnalysisSpec(
                profile=mode.value,
                modules=[AnalysisModule.QC, AnalysisModule.SV, AnalysisModule.REPORT],
            ),
        )

    def _context(self, mode: AssayMode, runner: _VersionRunner | None = None) -> RunContext:
        manifest = self._manifest(mode)
        lock = ReferenceLock(
            reference_id="TEST_REF",
            genome_build=GenomeBuild.GRCH38,
            contigs=[ReferenceContig(name="chr1", length=1000)],
            source_fai_sha256="a" * 64,
        )
        config = RunConfiguration(
            manifest=manifest,
            reference_lock=lock,
            output_base=self.base / "runs",
            run_id=manifest.run_id,
            pipeline_version="0.0.0-test",
            git_commit="0" * 40,
            qc_policy=QCPolicy(status="technical_defaults_only", note="test"),
            threads=2,
        )
        envelope = RunEnvelope.create(
            config.output_base,
            run_id=config.run_id,
            sample_id=manifest.sample_id,
        )
        return RunContext(
            config=config,
            envelope=envelope,
            runner=runner or _VersionRunner(),
            manifest=manifest,
        )

    def _report(self) -> TargetCoverageReport:
        return TargetCoverageReport(
            sample_id="SYNTHETIC_001",
            genome_build=GenomeBuild.GRCH38,
            target_bed_version="synthetic-v1",
            status=ModuleRunStatus.COMPLETED,
            policy=self.policy,
            summary_metrics={
                "region_count": 1,
                "interval_bases": 100,
                "interval_weighted_mean_depth": 18.5,
                "minimum_region_mean_depth": 18.5,
                "median_region_mean_depth": 18.5,
                "maximum_region_mean_depth": 18.5,
                "interval_bases_at_1x_fraction": 1.0,
                "interval_bases_at_10x_fraction": 0.9,
                "interval_bases_at_20x_fraction": 0.4,
                "interval_bases_at_30x_fraction": 0.1,
                "overlapping_interval_count": 0,
            },
            regions=[
                TargetCoverageRegion(
                    chromosome="chr1",
                    start=10,
                    end=110,
                    region_id="ROI_A",
                    mean_depth=18.5,
                    bases_at_threshold={"1x": 100, "10x": 90, "20x": 40, "30x": 10},
                    fraction_at_threshold={"1x": 1.0, "10x": 0.9, "20x": 0.4, "30x": 0.1},
                )
            ],
            target_bed_fingerprint=FileFingerprint(size_bytes=20, sha256="b" * 64),
            tool=ToolRecord(name="mosdepth", version="0.3.14", parameters={}),
            warnings=[self.policy.note],
            limitations=["Technical bins are descriptive only."],
        )

    def test_registration_promotes_only_engineering_verification(self) -> None:
        spec = SPEC_BY_STAGE[StageId.TARGET_COVERAGE]
        self.assertEqual(spec.verification, VerificationStatus.VERIFIED_WITH_REAL_TOOL)
        self.assertIn("observability", spec.purpose)
        self.assertIn(StageId.TARGET_COVERAGE, IMPLEMENTATIONS)

    def test_lcwgs_is_not_run_and_does_not_probe_mosdepth(self) -> None:
        runner = _VersionRunner()
        ctx = self._context(AssayMode.LOW_COVERAGE_WGS, runner)
        implementation = IMPLEMENTATIONS[StageId.TARGET_COVERAGE]
        plan = implementation.plan(ctx)
        result = implementation.execute(ctx, plan)
        self.assertEqual(plan.parameters, {"adaptive_sampling": False})
        self.assertEqual(result.status, ModuleRunStatus.NOT_RUN)
        self.assertIn("not applicable", result.reason)
        self.assertEqual(runner.calls, [])

    def test_adaptive_plan_locks_policy_tool_bed_and_bam(self) -> None:
        runner = _VersionRunner()
        ctx = self._context(AssayMode.ADAPTIVE_SAMPLING, runner)
        plan = IMPLEMENTATIONS[StageId.TARGET_COVERAGE].plan(ctx)
        self.assertEqual(plan.tool_versions, {"mosdepth": "0.3.14"})
        self.assertEqual(plan.parameters["target_bed_role"], "analysis_roi_unbuffered")
        self.assertEqual(plan.parameters["thresholds"], [1, 10, 20, 30])
        self.assertEqual(
            {name for name, _digest in plan.external_inputs},
            {"targets.bed", "sample.bam"},
        )
        self.assertEqual(runner.calls, [("mosdepth", "--version")])

    def test_version_mismatch_fails_before_analysis(self) -> None:
        ctx = self._context(AssayMode.ADAPTIVE_SAMPLING, _VersionRunner("0.3.13"))
        with self.assertRaisesRegex(ValueError, "does not match policy lock"):
            IMPLEMENTATIONS[StageId.TARGET_COVERAGE].plan(ctx)

    def test_execute_exports_only_normalized_observability_report(self) -> None:
        ctx = self._context(AssayMode.ADAPTIVE_SAMPLING)
        ctx.envelope.atomic_write_text("manifest/intake.json", "{}\n")
        implementation = IMPLEMENTATIONS[StageId.TARGET_COVERAGE]
        plan = implementation.plan(ctx)
        with (
            mock.patch(
                "ontseq_platform.target_coverage_extension.AlignedBamIntakeReport.model_validate_json",
                return_value=mock.sentinel.intake,
            ),
            mock.patch(
                "ontseq_platform.target_coverage_extension.run_target_coverage",
                return_value=self._report(),
            ) as run_mock,
        ):
            result = implementation.execute(ctx, plan)

        self.assertEqual(result.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(len(result.outputs), 1)
        artifact = result.outputs[0]
        self.assertEqual(artifact.relative_path, "qc/SYNTHETIC_001.target-coverage.json")
        self.assertTrue(artifact.exportable)
        self.assertIn("descriptive only", result.reason)
        self.assertIn("not assay-adequacy", result.reason)
        normalized = ctx.envelope.path(ctx.path(TARGET_COVERAGE_REPORT)).read_text(encoding="utf-8")
        self.assertNotIn(str(self.bam), normalized)
        self.assertNotIn(str(self.bed), normalized)
        called_output_dir = run_mock.call_args.kwargs["output_dir"]
        self.assertEqual(called_output_dir, ctx.envelope.path(ctx.path(TARGET_COVERAGE_WORK)))
        self.assertTrue(str(called_output_dir).startswith(str(ctx.envelope.root)))
        self.assertIn("work", called_output_dir.parts)


if __name__ == "__main__":
    unittest.main()
