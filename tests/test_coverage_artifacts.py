"""Coverage handoff: current-stage evidence inside a run, both historical names outside it.

Inside a run, SV observability and the reviewer reports may only consume the coverage the
target-coverage stage recorded for that run. Outside a run (re-rendering, the Befund view),
archived envelopes are read by file name; that reader still accepts the sample-named file the
retired runtime extension wrote and refuses an envelope holding two different answers.
"""

from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path
from unittest import mock

from ontseq_platform.coverage_artifacts import (
    ARCHIVED_COVERAGE_READER,
    COVERAGE_ARTIFACT_CONTRACT,
    load_run_coverage,
)
from ontseq_platform.demo import build_demo_result
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
from ontseq_platform.pipeline.runner import (
    LEGACY_TARGET_COVERAGE_REPORT,
    REPORT_HTML,
    RESULT_JSON,
    SELECTION_COVERAGE_REPORT,
    TARGET_COVERAGE_REPORT,
    RunConfiguration,
    RunContext,
    StageFailure,
    StagePlan,
    _report_execute,
    _sv_execute,
    load_current_coverage,
)
from ontseq_platform.pipeline.stages import StageId
from ontseq_platform.target_coverage import (
    TargetCoveragePolicy,
    TargetCoverageRegion,
    TargetCoverageReport,
)


class _NoCommands:
    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        raise AssertionError(f"unexpected command: {list(argv)}")

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        raise AssertionError(f"unexpected streamed command: {list(argv)}")


_POLICY = TargetCoveragePolicy(
    profile_id="adaptive_target_coverage_technical_v1",
    status="technical_defaults_only",
    expected_version="0.3.14",
    thresholds=[1, 10, 20, 30],
    mapq=0,
    exclude_flags=1796,
    note="Technical descriptive defaults only; not an adequacy threshold.",
)


def _report(*, region_id: str = "ROI_A", **update: object) -> TargetCoverageReport:
    report = TargetCoverageReport(
        sample_id="SYNTHETIC_001",
        genome_build=GenomeBuild.GRCH38,
        target_bed_version="synthetic-v1",
        status=ModuleRunStatus.COMPLETED,
        policy=_POLICY,
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
                region_id=region_id,
                mean_depth=18.5,
                bases_at_threshold={"1x": 100, "10x": 90, "20x": 40, "30x": 10},
                fraction_at_threshold={"1x": 1.0, "10x": 0.9, "20x": 0.4, "30x": 0.1},
            )
        ],
        target_bed_fingerprint=FileFingerprint(size_bytes=20, sha256="b" * 64),
        tool=ToolRecord(name="mosdepth", version="0.3.14", parameters={}),
        warnings=[_POLICY.note],
        limitations=["Technical bins are descriptive only."],
    )
    return report.model_copy(update=update) if update else report


class _Base(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.bam = self.base / "sample.bam"
        self.bam.write_bytes(b"synthetic-bam-placeholder")
        self.bed = self.base / "targets.bed"
        self.bed.write_text("chr1\t10\t110\tROI_A\n", encoding="utf-8")

    def _context(self, mode: AssayMode = AssayMode.ADAPTIVE_SAMPLING) -> RunContext:
        adaptive = mode == AssayMode.ADAPTIVE_SAMPLING
        manifest = SampleManifest(
            sample_id="SYNTHETIC_001",
            run_id="RUN_SYNTHETIC_001",
            input=InputSpec(
                kind=InputKind.ALIGNED_BAM,
                path=str(self.bam),
                index_path=str(self.base / "sample.bam.bai"),
            ),
            assay=AssaySpec(
                mode=mode,
                genome_build=GenomeBuild.GRCH38,
                reference_id="TEST_REF",
                target_bed=str(self.bed) if adaptive else None,
                target_bed_version="synthetic-v1" if adaptive else None,
            ),
            analysis=AnalysisSpec(
                profile=mode.value,
                modules=[AnalysisModule.QC, AnalysisModule.SV, AnalysisModule.REPORT],
            ),
        )
        config = RunConfiguration(
            manifest=manifest,
            reference_lock=ReferenceLock(
                reference_id="TEST_REF",
                genome_build=GenomeBuild.GRCH38,
                contigs=[ReferenceContig(name="chr1", length=1000)],
                source_fai_sha256="a" * 64,
            ),
            output_base=self.base / "runs",
            run_id=manifest.run_id,
            pipeline_version="0.0.0-test",
            git_commit="0" * 40,
            qc_policy=QCPolicy(status="technical_defaults_only", note="test"),
            threads=2,
        )
        envelope = RunEnvelope.create(
            config.output_base, run_id=config.run_id, sample_id=manifest.sample_id
        )
        return RunContext(config=config, envelope=envelope, runner=_NoCommands(), manifest=manifest)

    def _record(self, ctx: RunContext, relative: str, report: TargetCoverageReport) -> None:
        artifact = ctx.envelope.atomic_write_text(relative, report.model_dump_json())
        ctx.artifacts.setdefault(StageId.TARGET_COVERAGE, []).append(artifact)


class InRunHandoffTests(_Base):
    def test_the_handoff_identifiers_are_distinct_and_recorded(self) -> None:
        self.assertEqual(COVERAGE_ARTIFACT_CONTRACT, "current-stage-coverage-v2")
        self.assertEqual(ARCHIVED_COVERAGE_READER, "core-or-sample-coverage-v1")

    def test_a_file_on_disk_is_not_current_evidence(self) -> None:
        ctx = self._context()
        for name in (TARGET_COVERAGE_REPORT, ctx.path(LEGACY_TARGET_COVERAGE_REPORT)):
            ctx.envelope.atomic_write_text(name, _report().model_dump_json())
        self.assertIsNone(load_current_coverage(ctx))

    def test_the_recorded_stage_artifact_is_read(self) -> None:
        ctx = self._context()
        self._record(ctx, TARGET_COVERAGE_REPORT, _report())
        self.assertEqual(load_current_coverage(ctx), _report())
        self.assertIsNone(load_current_coverage(ctx, selection=True))

    def test_selection_coverage_is_read_only_as_selection_coverage(self) -> None:
        ctx = self._context()
        self._record(ctx, SELECTION_COVERAGE_REPORT, _report())
        self.assertIsNone(load_current_coverage(ctx))
        self.assertIsNotNone(load_current_coverage(ctx, selection=True))

    def test_a_changed_artifact_fails_its_checksum(self) -> None:
        ctx = self._context()
        self._record(ctx, TARGET_COVERAGE_REPORT, _report())
        ctx.envelope.path(TARGET_COVERAGE_REPORT).write_text(
            _report(target_bed_version="OTHER_PANEL").model_dump_json(), encoding="utf-8"
        )
        with self.assertRaisesRegex(StageFailure, "checksum"):
            load_current_coverage(ctx)

    def test_a_current_artifact_for_another_sample_or_build_is_refused(self) -> None:
        for update, reason in (
            ({"sample_id": "OTHER_SYNTHETIC"}, "different sample"),
            ({"genome_build": GenomeBuild.GRCH37}, "different genome build"),
        ):
            with self.subTest(reason=reason):
                ctx = self._context()
                self._record(ctx, TARGET_COVERAGE_REPORT, _report(**update))
                with self.assertRaisesRegex(StageFailure, reason):
                    load_current_coverage(ctx)
                ctx.envelope.path(TARGET_COVERAGE_REPORT).unlink()

    def _sv(self, ctx: RunContext):  # noqa: ANN202 - returns the stage result under test
        ctx.envelope.atomic_write_text("manifest/intake.json", "{}")
        with mock.patch(
            "ontseq_platform.pipeline.runner.AlignedBamIntakeReport.model_validate_json",
            return_value=mock.sentinel.intake,
        ):
            return _sv_execute(ctx, StagePlan(parameters={}, tool_versions={}))

    def test_sv_fails_closed_on_leftover_coverage_and_uses_the_current_report(self) -> None:
        ctx = self._context()
        ctx.envelope.atomic_write_text(
            ctx.path(LEGACY_TARGET_COVERAGE_REPORT), _report().model_dump_json()
        )
        with self.assertRaisesRegex(StageFailure, "requires a coverage report"):
            self._sv(ctx)

        self._record(ctx, TARGET_COVERAGE_REPORT, _report())
        outcome = self._sv(ctx)
        self.assertEqual(outcome.status, ModuleRunStatus.NO_CALL)
        self.assertTrue(any("consensus" in item.relative_path for item in outcome.outputs))

    def test_the_report_renders_only_current_coverage(self) -> None:
        ctx = self._context()
        result = build_demo_result()
        result.manifest = ctx.manifest
        ctx.envelope.atomic_write_text(ctx.path(RESULT_JSON), result.model_dump_json())
        ctx.envelope.atomic_write_text(
            ctx.path(LEGACY_TARGET_COVERAGE_REPORT),
            _report(region_id="LEFTOVER_ROI").model_dump_json(),
        )
        self._record(ctx, TARGET_COVERAGE_REPORT, _report())
        _report_execute(ctx, StagePlan(parameters={}, tool_versions={}))
        document = ctx.envelope.path(ctx.path(REPORT_HTML)).read_text(encoding="utf-8")
        self.assertIn("ROI_A", document)
        self.assertIn("18.5", document)
        self.assertNotIn("LEFTOVER_ROI", document)

    def test_sv_and_report_resume_track_coverage_artifacts(self) -> None:
        ctx = self._context()
        self._record(ctx, TARGET_COVERAGE_REPORT, _report())
        artifact = ctx.artifacts[StageId.TARGET_COVERAGE][0]
        for stage in (StageId.SV, StageId.REPORT):
            with self.subTest(stage=stage):
                self.assertIn(artifact, ctx.upstream(stage, InputKind.ALIGNED_BAM))


class ArchivedEnvelopeReaderTests(_Base):
    def test_either_historical_producer_name_is_accepted(self) -> None:
        ctx = self._context()
        for name in (TARGET_COVERAGE_REPORT, ctx.path(LEGACY_TARGET_COVERAGE_REPORT)):
            with self.subTest(name=name):
                artifact = ctx.envelope.atomic_write_text(name, _report().model_dump_json())
                self.assertEqual(load_run_coverage(ctx.envelope.root, ctx.manifest), _report())
                ctx.envelope.path(artifact.relative_path).unlink()
        self.assertIsNone(load_run_coverage(ctx.envelope.root, ctx.manifest))

    def test_mismatched_identity_and_conflicting_duplicates_are_refused(self) -> None:
        ctx = self._context()
        legacy = ctx.path(LEGACY_TARGET_COVERAGE_REPORT)
        for update, reason in (
            ({"sample_id": "OTHER_SYNTHETIC"}, "different sample"),
            ({"genome_build": GenomeBuild.GRCH37}, "different genome build"),
        ):
            with self.subTest(reason=reason):
                ctx.envelope.atomic_write_text(legacy, _report(**update).model_dump_json())
                with self.assertRaisesRegex(ValueError, reason):
                    load_run_coverage(ctx.envelope.root, ctx.manifest)
        ctx.envelope.atomic_write_text(legacy, _report().model_dump_json())
        ctx.envelope.atomic_write_text(
            TARGET_COVERAGE_REPORT,
            _report(target_bed_version="OTHER_PANEL").model_dump_json(),
        )
        with self.assertRaisesRegex(ValueError, "Conflicting"):
            load_run_coverage(ctx.envelope.root, ctx.manifest)

    def test_selection_coverage_is_not_substituted_for_analysis_coverage(self) -> None:
        ctx = self._context()
        ctx.envelope.atomic_write_text(SELECTION_COVERAGE_REPORT, _report().model_dump_json())
        self.assertIsNone(load_run_coverage(ctx.envelope.root, ctx.manifest))
        self.assertIsNotNone(load_run_coverage(ctx.envelope.root, ctx.manifest, selection=True))


if __name__ == "__main__":
    unittest.main()
