"""Runner behaviour proven without any bioinformatics binary.

The stage implementations are swapped for fakes so the three properties that matter can be
tested in isolation from the tools: failures propagate as NOT_RUN rather than FAILED,
resume is content-addressed rather than timestamp-based, and the release bundle withholds
what must not leave the execution system.
"""

from __future__ import annotations

import hashlib
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
)
from ontseq_platform.pipeline.envelope import RunEnvelope
from ontseq_platform.pipeline.lock import LOCK_FILENAME, RunAlreadyRunning, run_lock
from ontseq_platform.pipeline.runner import (
    ALIGNED_BAI,
    ALIGNED_BAM,
    IMPLEMENTATIONS,
    RunConfiguration,
    RunContext,
    StageFailure,
    StageImplementation,
    StagePlan,
    StageResult,
    run_pipeline,
)
from ontseq_platform.pipeline.stages import StageId


class _NullRunner:
    """A command runner that must never be called; fake stages run no tools."""

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        raise AssertionError(f"no fake stage should execute a command: {list(argv)}")

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        raise AssertionError(f"no fake stage should execute a command: {list(argv)}")


class _FakeStage:
    """A stage implementation whose behaviour and call count the test controls."""

    def __init__(
        self,
        stage: StageId,
        *,
        relative_path: str | None = None,
        fail: bool = False,
        parameter: str = "v1",
        status: ModuleRunStatus = ModuleRunStatus.COMPLETED,
    ) -> None:
        self.stage = stage
        self.relative_path = relative_path
        self.fail = fail
        self.parameter = parameter
        self.status = status
        self.executions = 0

    def plan(self, ctx: RunContext) -> StagePlan:
        return StagePlan(parameters={"knob": self.parameter}, tool_versions={"fake": "1.0.0"})

    def execute(self, ctx: RunContext, plan: StagePlan) -> StageResult:
        self.executions += 1
        if self.fail:
            raise StageFailure(f"{self.stage.value} was told to fail")
        outputs = []
        if self.relative_path is not None:
            outputs.append(
                ctx.envelope.atomic_write_text(
                    ctx.path(self.relative_path), f"{self.stage.value}:{self.parameter}\n"
                )
            )
        return StageResult(
            status=self.status,
            reason=f"{self.stage.value} finished as {self.status.value}.",
            outputs=outputs,
        )

    def implementation(self) -> StageImplementation:
        return StageImplementation(self.plan, self.execute)


#: One artifact per stage, chosen so the release bundle has both exportable and withheld
#: entries: alignment/ and work/ are intermediates that must never be exported.
STAGE_ARTIFACTS: dict[StageId, str | None] = {
    StageId.INTAKE: "manifest/intake.json",
    StageId.QC: "qc/cramino.json",
    StageId.SV: "evidence/sv/{sample}.sniffles.vcf",
    StageId.ASSEMBLE: "normalized/{sample}.result.json",
    StageId.REPORT: "reports/{sample}.report.html",
    StageId.RELEASE: None,
}


def _manifest() -> SampleManifest:
    return SampleManifest(
        sample_id="FAKE_RUNNER_001",
        run_id="FAKE_RUNNER_RUN_001",
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
        analysis=AnalysisSpec(profile="fake", modules=[AnalysisModule.QC]),
    )


def _reference_lock() -> ReferenceLock:
    return ReferenceLock(
        reference_id="FAKE_REFERENCE_V1",
        genome_build=GenomeBuild.GRCH38,
        contigs=[ReferenceContig(name="chr1", length=1000)],
        source_fai_sha256="a" * 64,
    )


class RunnerCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.base = Path(self._temporary.name)
        self.stages = {
            stage: _FakeStage(stage, relative_path=path) for stage, path in STAGE_ARTIFACTS.items()
        }

    def _config(self, **overrides: object) -> RunConfiguration:
        values: dict[str, object] = {
            "manifest": _manifest(),
            "reference_lock": _reference_lock(),
            "output_base": self.base / "runs",
            "run_id": "FAKE_RUN_001",
            "pipeline_version": "0.0.0-test",
            "git_commit": "0" * 40,
            "qc_policy": QCPolicy(status="technical_defaults_only", note="test"),
        }
        values.update(overrides)
        return RunConfiguration(**values)  # type: ignore[arg-type]

    def _run(self, **overrides: object):
        implementations = {stage: fake.implementation() for stage, fake in self.stages.items()}
        with mock.patch.dict(
            "ontseq_platform.pipeline.runner.IMPLEMENTATIONS", implementations, clear=True
        ):
            return run_pipeline(self._config(**overrides), runner=_NullRunner())

    def _envelope_root(self) -> Path:
        return self.base / "runs" / "FAKE_RUN_001" / "FAKE_RUNNER_001"


class HappyPathTests(RunnerCase):
    def test_every_planned_stage_is_recorded(self) -> None:
        report, bundle = self._run()
        recorded = {record.stage for record in report.stages}
        self.assertEqual(
            recorded,
            set(STAGE_ARTIFACTS) | {StageId.TARGET_COVERAGE, StageId.CNV, StageId.METHYLATION},
        )
        self.assertTrue(report.passed)
        self.assertIsNotNone(bundle)

    def test_basecalling_and_alignment_are_absent_for_an_aligned_bam(self) -> None:
        report, _ = self._run()
        recorded = {record.stage for record in report.stages}
        self.assertNotIn(StageId.BASECALL, recorded)
        self.assertNotIn(StageId.ALIGN, recorded)

    def test_unwired_stages_are_not_run_and_say_why(self) -> None:
        report, _ = self._run()
        record = report.record_for(StageId.CNV)
        self.assertIsNotNone(record)
        self.assertEqual(record.status, ModuleRunStatus.NOT_RUN)
        self.assertIn("No adapter is wired in", record.reason)

    def test_the_run_report_is_checksummed_into_its_own_bundle(self) -> None:
        _, bundle = self._run()
        assert bundle is not None
        paths = {item.relative_path for item in bundle.artifacts}
        self.assertIn("provenance/run.json", paths)
        self.assertEqual(
            bundle.run_report_sha256,
            next(i.sha256 for i in bundle.artifacts if i.relative_path == "provenance/run.json"),
        )

    def test_checksum_manifest_matches_the_files_on_disk(self) -> None:
        _, bundle = self._run()
        assert bundle is not None
        root = self._envelope_root()
        for line in bundle.checksum_manifest().splitlines():
            digest, relative = line.split("  ", 1)
            self.assertEqual(hashlib.sha256((root / relative).read_bytes()).hexdigest(), digest)


class ExportBoundaryTests(RunnerCase):
    def test_intermediate_artifacts_are_withheld_from_the_bundle(self) -> None:
        self.stages[StageId.SV].relative_path = "alignment/{sample}.bam"
        _, bundle = self._run()
        assert bundle is not None
        self.assertIn("alignment/FAKE_RUNNER_001.bam", bundle.withheld_artifact_paths)
        exported = {item.relative_path for item in bundle.artifacts}
        self.assertNotIn("alignment/FAKE_RUNNER_001.bam", exported)


class FailurePropagationTests(RunnerCase):
    def _failed_run(self):
        self.stages[StageId.QC].fail = True
        return self._run()

    def test_the_failing_stage_is_failed(self) -> None:
        report, _ = self._failed_run()
        record = report.record_for(StageId.QC)
        self.assertEqual(record.status, ModuleRunStatus.FAILED)
        self.assertIn("told to fail", record.reason)

    def test_downstream_stages_are_not_run_rather_than_failed(self) -> None:
        """A stage that never started has not failed; conflating the two invents evidence."""
        report, _ = self._failed_run()
        for stage in (StageId.ASSEMBLE, StageId.REPORT, StageId.RELEASE):
            self.assertEqual(
                report.record_for(stage).status, ModuleRunStatus.NOT_RUN, f"{stage} status"
            )

    def test_a_stage_that_does_not_depend_on_the_failure_still_runs(self) -> None:
        """SV depends on the intake gate, not on QC, so a QC failure must not silence it."""
        report, _ = self._failed_run()
        self.assertEqual(report.record_for(StageId.SV).status, ModuleRunStatus.COMPLETED)

    def test_the_blocked_stage_names_its_cause(self) -> None:
        report, _ = self._failed_run()
        self.assertIn("qc", report.record_for(StageId.ASSEMBLE).reason)

    def test_a_failed_run_produces_no_release_bundle(self) -> None:
        report, bundle = self._failed_run()
        self.assertFalse(report.passed)
        self.assertIsNone(bundle)

    def test_downstream_stages_never_executed(self) -> None:
        self._failed_run()
        self.assertEqual(self.stages[StageId.ASSEMBLE].executions, 0)
        self.assertEqual(self.stages[StageId.REPORT].executions, 0)

    def test_a_truthful_partial_report_is_still_written(self) -> None:
        self._failed_run()
        self.assertTrue((self._envelope_root() / "provenance" / "run.json").is_file())


class ResumeTests(RunnerCase):
    def test_an_unchanged_second_run_resumes_every_stage(self) -> None:
        self._run()
        report, _ = self._run()
        for record in report.stages:
            if record.status == ModuleRunStatus.COMPLETED:
                self.assertTrue(record.resumed, f"{record.stage} did not resume")

    def test_a_resumed_stage_does_not_execute_again(self) -> None:
        self._run()
        before = {stage: fake.executions for stage, fake in self.stages.items()}
        self._run()
        for stage, fake in self.stages.items():
            self.assertEqual(fake.executions, before[stage], f"{stage} re-executed")

    def test_a_changed_parameter_forces_the_stage_to_re_execute(self) -> None:
        self._run()
        self.stages[StageId.QC].parameter = "v2"
        report, _ = self._run()
        self.assertFalse(report.record_for(StageId.QC).resumed)
        self.assertEqual(self.stages[StageId.QC].executions, 2)

    def test_a_changed_upstream_artifact_re_executes_the_dependent(self) -> None:
        """Resume must follow the content of upstream outputs, not just local parameters."""
        self._run()
        self.stages[StageId.INTAKE].parameter = "v2"
        report, _ = self._run()
        self.assertFalse(report.record_for(StageId.INTAKE).resumed)
        self.assertFalse(report.record_for(StageId.QC).resumed)

    def test_a_tampered_artifact_re_executes_the_stage(self) -> None:
        """A matching signature is not enough; the artifact must still be byte identical."""
        self._run()
        (self._envelope_root() / "qc" / "cramino.json").write_text("tampered\n", encoding="utf-8")
        report, _ = self._run()
        self.assertFalse(report.record_for(StageId.QC).resumed)
        self.assertEqual(self.stages[StageId.QC].executions, 2)

    def test_a_deleted_artifact_re_executes_the_stage(self) -> None:
        self._run()
        (self._envelope_root() / "qc" / "cramino.json").unlink()
        report, _ = self._run()
        self.assertFalse(report.record_for(StageId.QC).resumed)

    def test_force_re_executes_everything(self) -> None:
        self._run()
        report, _ = self._run(force=True)
        for record in report.stages:
            self.assertFalse(record.resumed, f"{record.stage} resumed under --force")

    def test_a_no_call_stage_resumes_like_a_completed_one(self) -> None:
        """A caller that looked and declined has concluded; re-running it changes nothing."""
        self.stages[StageId.SV].status = ModuleRunStatus.NO_CALL
        self._run()
        report, _ = self._run()
        record = report.record_for(StageId.SV)
        self.assertEqual(record.status, ModuleRunStatus.NO_CALL)
        self.assertTrue(record.resumed)
        self.assertEqual(self.stages[StageId.SV].executions, 1)

    def test_a_run_with_a_resumed_no_call_still_passes_and_releases(self) -> None:
        self.stages[StageId.SV].status = ModuleRunStatus.NO_CALL
        self._run()
        report, bundle = self._run()
        self.assertTrue(report.passed)
        self.assertIsNotNone(bundle)

    def test_a_no_call_stage_keeps_its_artifacts_through_a_resume(self) -> None:
        """Sniffles2 writes a VCF even when nothing passes its policy; it must survive."""
        self.stages[StageId.SV].status = ModuleRunStatus.NO_CALL
        self._run()
        report, _ = self._run()
        outputs = [item.relative_path for item in report.record_for(StageId.SV).outputs]
        self.assertEqual(outputs, ["evidence/sv/FAKE_RUNNER_001.sniffles.vcf"])

    def test_resuming_after_a_failure_re_attempts_the_failed_stage(self) -> None:
        self.stages[StageId.QC].fail = True
        self._run()
        self.stages[StageId.QC].fail = False
        report, bundle = self._run()
        self.assertTrue(report.passed)
        self.assertIsNotNone(bundle)
        self.assertTrue(report.record_for(StageId.INTAKE).resumed)
        self.assertFalse(report.record_for(StageId.QC).resumed)


class LockingTests(RunnerCase):
    """The run holds an exclusive lock; a watcher that double-fires must not get through."""

    def _hold(self):
        return run_lock(
            self._envelope_root(),
            run_id="OTHER_RUN",
            sample_id="FAKE_RUNNER_001",
            pipeline_version="0.0.0-test",
        )

    def test_a_run_refuses_while_the_envelope_is_locked(self) -> None:
        with self._hold(), self.assertRaises(RunAlreadyRunning):
            self._run()

    def test_a_blocked_run_executes_no_stage(self) -> None:
        """Refusing has to happen before any write, not halfway through the graph."""
        with self._hold(), self.assertRaises(RunAlreadyRunning):
            self._run()
        self.assertEqual(sum(fake.executions for fake in self.stages.values()), 0)

    def test_the_lock_is_released_when_the_run_finishes(self) -> None:
        self._run()
        self.assertFalse((self._envelope_root() / LOCK_FILENAME).exists())

    def test_a_reclaimed_stale_lock_is_recorded_in_the_run_report(self) -> None:
        """Stepping over another run's lock is not something to leave in a terminal."""
        self._envelope_root().mkdir(parents=True, exist_ok=True)
        (self._envelope_root() / LOCK_FILENAME).write_text("", encoding="utf-8")
        report, _ = self._run()
        self.assertTrue(any("stale run lock" in warning for warning in report.warnings))


class SettleTests(RunnerCase):
    """A stage's effect on the run context must survive being resumed."""

    def setUp(self) -> None:
        super().setUp()
        self.settled: list[str] = []

        def settle(ctx: RunContext, outputs) -> None:  # type: ignore[no-untyped-def]
            self.settled.append(",".join(item.relative_path for item in outputs))

        fake = self.stages[StageId.INTAKE]
        self.intake_implementation = StageImplementation(fake.plan, fake.execute, settle)
        fake.implementation = lambda: self.intake_implementation  # type: ignore[method-assign]

    def test_settle_runs_when_the_stage_executes(self) -> None:
        self._run()
        self.assertEqual(self.settled, ["manifest/intake.json"])

    def test_settle_runs_again_when_the_stage_resumes(self) -> None:
        """Otherwise a resumed alignment leaves the manifest pointing at its own input."""
        self._run()
        self._run()
        self.assertEqual(self.settled, ["manifest/intake.json", "manifest/intake.json"])

    def test_settle_is_given_the_artifacts_the_stage_recorded(self) -> None:
        self._run()
        self.assertIn("manifest/intake.json", self.settled[0])


class AlignSettleTests(unittest.TestCase):
    """The real alignment hook, tested without running minimap2."""

    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        envelope = RunEnvelope.create(
            Path(self._temporary.name), run_id="R1", sample_id="FAKE_RUNNER_001"
        )
        config = RunConfiguration(
            manifest=_manifest(),
            reference_lock=_reference_lock(),
            output_base=Path(self._temporary.name),
            run_id="R1",
            pipeline_version="0.0.0-test",
            git_commit="0" * 40,
            qc_policy=QCPolicy(status="technical_defaults_only", note="test"),
        )
        self.context = RunContext(
            config=config,
            envelope=envelope,
            runner=_NullRunner(),
            manifest=_manifest().model_copy(
                update={"input": InputSpec(kind=InputKind.UNALIGNED_BAM, path="/nowhere/reads.bam")}
            ),
        )
        self.outputs = [
            envelope.atomic_write_text(self.context.path(ALIGNED_BAM), "bam\n"),
            envelope.atomic_write_text(self.context.path(ALIGNED_BAI), "bai\n"),
        ]
        self.settle = IMPLEMENTATIONS[StageId.ALIGN].settle
        assert self.settle is not None

    def test_the_manifest_is_repointed_at_the_aligned_bam(self) -> None:
        self.settle(self.context, self.outputs)
        self.assertEqual(self.context.manifest.input.kind, InputKind.ALIGNED_BAM)
        self.assertTrue(self.context.manifest.input.path.endswith("alignment/FAKE_RUNNER_001.bam"))
        self.assertEqual(
            self.context.manifest.input.index_path, f"{self.context.manifest.input.path}.bai"
        )

    def test_the_recorded_checksum_is_reused_rather_than_recomputed(self) -> None:
        self.settle(self.context, self.outputs)
        self.assertEqual(self.context.manifest.input.sha256, self.outputs[0].sha256)

    def test_a_missing_index_is_refused(self) -> None:
        with self.assertRaises(StageFailure):
            self.settle(self.context, self.outputs[:1])

    def test_a_missing_bam_is_refused(self) -> None:
        with self.assertRaises(StageFailure):
            self.settle(self.context, [])


if __name__ == "__main__":
    unittest.main()
