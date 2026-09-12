"""The target-coverage stage inside the runner, and component selection around it.

Two properties are load-bearing here. An adaptive-sampling run that cannot measure its own
enrichment must stop rather than produce a report that looks complete; and a run that was
never adaptive sampling in the first place must say so as a scope statement, not as a
coverage result. The difference between those two is the whole reason the stage exists.
"""

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
    GenomeBuild,
    InputKind,
    InputSpec,
    ModuleRunStatus,
    QCPolicy,
    ReferenceContig,
    ReferenceLock,
    SampleManifest,
    TargetBedRole,
)
from ontseq_platform.pipeline.components import RunComponents
from ontseq_platform.pipeline.envelope import RunEnvelope
from ontseq_platform.pipeline.runner import (
    COMPONENTS_REPORT,
    RunConfiguration,
    RunContext,
    StageFailure,
    StageImplementation,
    StagePlan,
    StageResult,
    _target_coverage_execute,
    _target_coverage_plan,
    run_pipeline,
)
from ontseq_platform.pipeline.stages import StageId
from ontseq_platform.target_coverage import TargetCoveragePolicy


class _VersionRunner:
    """Answers a version probe and refuses anything else."""

    def __init__(self, version: str = "0.3.14") -> None:
        self.version = version
        self.calls: list[list[str]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        self.calls.append(list(argv))
        if "--version" in argv:
            return CommandResult(
                argv=list(argv), returncode=0, stdout=f"mosdepth {self.version}\n", stderr=""
            )
        raise AssertionError(f"unexpected command: {list(argv)}")

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        raise AssertionError("target coverage does not stream to a file")


def _policy() -> TargetCoveragePolicy:
    return TargetCoveragePolicy(
        profile_id="test-target-coverage",
        status="technical_defaults_only",
        note="test policy",
    )


def _reference_lock() -> ReferenceLock:
    return ReferenceLock(
        reference_id="FAKE_REFERENCE_V1",
        genome_build=GenomeBuild.GRCH38,
        contigs=[ReferenceContig(name="chr1", length=1000)],
        source_fai_sha256="a" * 64,
    )


class StageBehaviourTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.bam = self.base / "input.bam"
        self.bam.write_bytes(b"not really a bam")
        (self.base / "input.bam.bai").write_bytes(b"index")
        self.bed = self.base / "panel.bed"
        self.bed.write_text("chr1\t100\t200\tTARGET_A\n", encoding="utf-8")

    def _manifest(self, mode: AssayMode, *, with_bed: bool = True) -> SampleManifest:
        assay = AssaySpec(
            mode=mode,
            genome_build=GenomeBuild.GRCH38,
            reference_id="FAKE_REFERENCE_V1",
            target_bed=str(self.bed) if with_bed else None,
            target_bed_version="TEST_PANEL_V1" if with_bed else None,
            target_bed_role=TargetBedRole.SELECTION_PANEL_BUFFERED,
        )
        return SampleManifest(
            sample_id="TC_001",
            run_id="TC_RUN_001",
            input=InputSpec(
                kind=InputKind.ALIGNED_BAM,
                path=str(self.bam),
                index_path=str(self.bam) + ".bai",
            ),
            assay=assay,
            analysis=AnalysisSpec(profile="test", modules=[AnalysisModule.QC]),
        )

    def _context(
        self, manifest: SampleManifest, *, policy: TargetCoveragePolicy | None
    ) -> RunContext:
        config = RunConfiguration(
            manifest=manifest,
            reference_lock=_reference_lock(),
            output_base=self.base / "runs",
            run_id="TC_RUN_001",
            pipeline_version="0.0.0-test",
            git_commit="0" * 40,
            qc_policy=QCPolicy(status="technical_defaults_only", note="test"),
            target_coverage_policy=policy,
        )
        envelope = RunEnvelope.create(
            config.output_base, run_id=config.run_id, sample_id=manifest.sample_id
        )
        return RunContext(
            config=config, envelope=envelope, runner=_VersionRunner(), manifest=manifest
        )

    def test_a_non_adaptive_run_records_scope_rather_than_coverage(self) -> None:
        context = self._context(
            self._manifest(AssayMode.LOW_COVERAGE_WGS, with_bed=False), policy=_policy()
        )
        plan = _target_coverage_plan(context)
        self.assertEqual(plan.tool_versions, {})
        result = _target_coverage_execute(context, plan)
        self.assertIs(result.status, ModuleRunStatus.NOT_RUN)
        self.assertIn("does not apply", result.reason)
        self.assertIn("not a coverage finding", result.reason)

    def test_a_non_adaptive_run_never_probes_for_mosdepth(self) -> None:
        """An lcWGS run must not require a tool it will never execute."""
        context = self._context(
            self._manifest(AssayMode.LOW_COVERAGE_WGS, with_bed=False), policy=None
        )
        _target_coverage_plan(context)
        runner = context.runner
        assert isinstance(runner, _VersionRunner)
        self.assertEqual(runner.calls, [])

    def test_adaptive_sampling_without_a_policy_fails_closed(self) -> None:
        context = self._context(self._manifest(AssayMode.ADAPTIVE_SAMPLING), policy=None)
        with self.assertRaises(StageFailure) as caught:
            _target_coverage_plan(context)
        self.assertIn("Refusing to continue", str(caught.exception))

    def test_the_plan_fingerprints_the_panel_and_records_its_role(self) -> None:
        context = self._context(self._manifest(AssayMode.ADAPTIVE_SAMPLING), policy=_policy())
        plan = _target_coverage_plan(context)
        self.assertEqual(plan.tool_versions, {"mosdepth": "0.3.14"})
        self.assertEqual(plan.parameters["target_bed_version"], "TEST_PANEL_V1")
        self.assertEqual(
            plan.parameters["target_bed_role"], TargetBedRole.SELECTION_PANEL_BUFFERED.value
        )
        fingerprinted = {name for name, _ in plan.external_inputs}
        self.assertIn("panel.bed", fingerprinted)

    def test_changing_the_panel_changes_the_plan(self) -> None:
        """Resume keys on the plan, so a new panel must not be silently reused."""
        context = self._context(self._manifest(AssayMode.ADAPTIVE_SAMPLING), policy=_policy())
        before = _target_coverage_plan(context).external_inputs
        self.bed.write_text("chr1\t100\t300\tTARGET_A\n", encoding="utf-8")
        after = _target_coverage_plan(context).external_inputs
        self.assertNotEqual(before, after)

    def test_selection_and_analysis_roi_are_fingerprinted_as_distinct_inputs(self) -> None:
        context = self._context(self._manifest(AssayMode.ADAPTIVE_SAMPLING), policy=_policy())
        selection = self.base / "selection-buffered.bed"
        selection.write_text("chr1\t50\t250\tTARGET_A\n", encoding="utf-8")
        context.config.selection_target_bed = selection

        plan = _target_coverage_plan(context)

        fingerprinted = {name for name, _checksum in plan.external_inputs}
        self.assertIn("panel.bed", fingerprinted)
        self.assertIn("selection_panel_buffered", fingerprinted)


class _FakeStage:
    def __init__(self, stage: StageId, tools: dict[str, str] | None = None) -> None:
        self.stage = stage
        #: What this fake claims to have probed. Empty means the stage ran no external
        #: tool, which is exactly the case a pinned version must not fail.
        self.tools = {"sniffles": "1.0.0"} if tools is None else tools
        self.executions = 0

    def plan(self, ctx: RunContext) -> StagePlan:
        return StagePlan(parameters={}, tool_versions=dict(self.tools))

    def execute(self, ctx: RunContext, plan: StagePlan) -> StageResult:
        self.executions += 1
        return StageResult(status=ModuleRunStatus.COMPLETED, reason="fake stage finished.")

    def implementation(self) -> StageImplementation:
        return StageImplementation(self.plan, self.execute)


class _NullRunner:
    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        raise AssertionError("no fake stage runs a command")

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        raise AssertionError("no fake stage runs a command")


class SelectionInsideTheRunnerTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.stages = {
            stage: _FakeStage(stage, tools={} if stage is StageId.TARGET_COVERAGE else None)
            for stage in (
                StageId.INTAKE,
                StageId.QC,
                StageId.TARGET_COVERAGE,
                StageId.CNV,
                StageId.SV,
                StageId.ASSEMBLE,
                StageId.REPORT,
                StageId.RELEASE,
            )
        }

    def _run(self, selection: RunComponents | None):
        manifest = SampleManifest(
            sample_id="SEL_001",
            run_id="SEL_RUN_001",
            input=InputSpec(
                kind=InputKind.ALIGNED_BAM,
                path="/nowhere/input.bam",
                index_path="/nowhere/input.bam.bai",
            ),
            assay=AssaySpec(
                mode=AssayMode.LOW_COVERAGE_WGS,
                genome_build=GenomeBuild.GRCH38,
                reference_id="FAKE_REFERENCE_V1",
            ),
            analysis=AnalysisSpec(profile="test", modules=[AnalysisModule.QC]),
        )
        config = RunConfiguration(
            manifest=manifest,
            reference_lock=_reference_lock(),
            output_base=self.base / "runs",
            run_id="SEL_RUN_001",
            pipeline_version="0.0.0-test",
            git_commit="0" * 40,
            qc_policy=QCPolicy(status="technical_defaults_only", note="test"),
            components=selection,
        )
        implementations = {stage: fake.implementation() for stage, fake in self.stages.items()}
        with mock.patch.dict(
            "ontseq_platform.pipeline.runner.IMPLEMENTATIONS", implementations, clear=True
        ):
            return run_pipeline(config, runner=_NullRunner())

    @staticmethod
    def _selection(**components: object) -> RunComponents:
        return RunComponents.model_validate(
            {
                "selection_id": "runner-selection",
                "status": "technical_defaults_only",
                "components": components,
            }
        )

    def test_a_deselected_stage_is_not_run_and_names_the_selection(self) -> None:
        report, _ = self._run(self._selection(sv={"provider": "sniffles2", "enabled": False}))
        record = report.record_for(StageId.SV)
        assert record is not None
        self.assertIs(record.status, ModuleRunStatus.NOT_RUN)
        self.assertIn("runner-selection", record.reason)
        self.assertIn("not a negative finding", record.reason)
        self.assertEqual(self.stages[StageId.SV].executions, 0)

    def test_a_version_mismatch_fails_the_stage_and_names_both_versions(self) -> None:
        report, _ = self._run(self._selection(sv={"provider": "sniffles2", "version": "2.4.0"}))
        record = report.record_for(StageId.SV)
        assert record is not None
        self.assertIs(record.status, ModuleRunStatus.FAILED)
        self.assertIn("2.4.0", record.reason)
        self.assertIn("1.0.0", record.reason)
        self.assertEqual(self.stages[StageId.SV].executions, 0)

    def test_a_matching_version_runs_the_stage(self) -> None:
        report, _ = self._run(self._selection(sv={"provider": "sniffles2", "version": "1.0.0"}))
        record = report.record_for(StageId.SV)
        assert record is not None
        self.assertIs(record.status, ModuleRunStatus.COMPLETED)

    def test_the_selection_is_written_into_the_envelope_before_anything_runs(self) -> None:
        self._run(self._selection(sv={"provider": "sniffles2", "version": "1.0.0"}))
        recorded = self.base / "runs" / "SEL_RUN_001" / "SEL_001" / COMPONENTS_REPORT
        self.assertTrue(recorded.is_file())
        self.assertIn("runner-selection", recorded.read_text(encoding="utf-8"))

    def test_an_unpinned_component_warns_that_the_run_is_not_reproducible(self) -> None:
        report, _ = self._run(self._selection(sv={"provider": "sniffles2"}))
        self.assertTrue(any("does not pin a version" in warning for warning in report.warnings))

    def test_a_pin_does_not_fail_a_stage_that_ran_no_tool(self) -> None:
        """Target coverage on an lcWGS run probes nothing; a mosdepth pin must not fail it."""
        report, _ = self._run(
            self._selection(target_coverage={"provider": "mosdepth", "version": "0.3.14"})
        )
        record = report.record_for(StageId.TARGET_COVERAGE)
        assert record is not None
        self.assertIsNot(record.status, ModuleRunStatus.FAILED)

    def test_no_selection_leaves_every_stage_alone(self) -> None:
        report, _ = self._run(None)
        record = report.record_for(StageId.SV)
        assert record is not None
        self.assertIs(record.status, ModuleRunStatus.COMPLETED)


if __name__ == "__main__":
    unittest.main()
