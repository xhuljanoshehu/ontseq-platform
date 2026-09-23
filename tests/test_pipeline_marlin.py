"""The native MARLIN lane is automatic, optional, and bound to current evidence."""

from pathlib import Path

from test_pipeline_runner import _manifest, _reference_lock

from ontseq_platform.models import AnalysisModule, QCPolicy
from ontseq_platform.pipeline.runner import RunConfiguration
from ontseq_platform.pipeline.stages import InputKindName, planned_stages


def configuration(tmp_path: Path, *, methylation: bool = True) -> RunConfiguration:
    manifest = _manifest()
    if methylation:
        manifest.analysis.modules.append(AnalysisModule.METHYLATION)
    return RunConfiguration(
        manifest=manifest,
        reference_lock=_reference_lock(),
        output_base=tmp_path,
        run_id=manifest.run_id,
        pipeline_version="0.8.2",
        git_commit="a" * 40,
        qc_policy=QCPolicy(status="technical_defaults_only", note="synthetic"),
    )


def test_selected_methylation_requests_marlin_in_persisted_scope(tmp_path: Path) -> None:
    config = configuration(tmp_path)
    assert "marlin" in [module.value for module in config.manifest.analysis.modules]
    assert config.manifest.analysis.modules.count(AnalysisModule.METHYLATION) == 1


def test_without_methylation_marlin_is_not_requested(tmp_path: Path) -> None:
    config = configuration(tmp_path, methylation=False)
    assert "marlin" not in [module.value for module in config.manifest.analysis.modules]


def test_marlin_is_scheduled_before_result_assembly() -> None:
    stages = [stage.value for stage in planned_stages(InputKindName.ALIGNED_BAM)]
    assert "marlin" in stages
    assert stages.index("methylation") < stages.index("marlin") < stages.index("assemble")


def test_missing_installation_records_actionable_not_run_without_tools(tmp_path: Path) -> None:
    from test_pipeline_runner import _NullRunner

    from ontseq_platform.models import ModuleRunStatus
    from ontseq_platform.pipeline.envelope import RunEnvelope
    from ontseq_platform.pipeline.runner import RunContext, _execute_stage
    from ontseq_platform.pipeline.stages import StageId, StageOutcome

    config = configuration(tmp_path)
    ctx = RunContext(
        config=config,
        manifest=config.manifest,
        runner=_NullRunner(),
        envelope=RunEnvelope.create(
            tmp_path, run_id=config.run_id, sample_id=config.manifest.sample_id
        ),
    )
    record = _execute_stage(
        StageId.MARLIN,
        ctx,
        InputKindName.ALIGNED_BAM,
        {StageId.INTAKE: StageOutcome.COMPLETED},
        None,
    )
    assert record.status is ModuleRunStatus.NOT_RUN
    assert "installation" in record.reason.lower()
    assert not record.outputs


def test_optional_marlin_failure_keeps_completed_genome_run_available():
    from ontseq_platform.pipeline.stages import StageId, StageOutcome, summarize

    outcomes = {
        stage: StageOutcome.COMPLETED for stage in planned_stages(InputKindName.ALIGNED_BAM)
    }
    outcomes[StageId.MARLIN] = StageOutcome.FAILED
    verdict = summarize(InputKindName.ALIGNED_BAM, outcomes)
    assert verdict.passed
    assert StageId.MARLIN in verdict.failed_stages
    assert "marlin" in verdict.describe().lower() and "FAILED" in verdict.describe()
    outcomes[StageId.QC] = StageOutcome.FAILED
    assert not summarize(InputKindName.ALIGNED_BAM, outcomes).passed


def completed_context(tmp_path):
    from test_marlin_report_presentation import native_report
    from test_pipeline_runner import _NullRunner

    from ontseq_platform.demo import build_demo_result
    from ontseq_platform.pipeline.envelope import RunEnvelope
    from ontseq_platform.pipeline.marlin import MARLIN_REPORT
    from ontseq_platform.pipeline.runner import RunContext
    from ontseq_platform.pipeline.stages import SPEC_BY_STAGE, StageId
    from ontseq_platform.pipeline.state import ArtifactRecord, StageRecord

    result = build_demo_result()
    config = configuration(tmp_path)
    config.manifest = result.manifest.model_copy(deep=True)
    config.manifest.analysis.modules.append(AnalysisModule.MARLIN)
    config.run_id = result.manifest.run_id
    ctx = RunContext(
        config=config,
        manifest=config.manifest,
        runner=_NullRunner(),
        envelope=RunEnvelope.create(
            tmp_path, run_id=config.run_id, sample_id=config.manifest.sample_id
        ),
    )
    report = native_report(result)
    artifact = ctx.envelope.atomic_write_text(ctx.path(MARLIN_REPORT), report.model_dump_json())
    ctx.artifacts[StageId.MARLIN] = [artifact]
    spec = SPEC_BY_STAGE[StageId.MARLIN]
    ctx.stage_records[StageId.MARLIN] = StageRecord(
        stage=StageId.MARLIN,
        title=spec.title,
        status=report.status,
        verification=spec.verification,
        required=False,
        reason=report.reason,
        signature="a" * 64,
        outputs=[ArtifactRecord.of(artifact)],
        tools=report.tools,
    )
    return ctx, result


def test_failed_stage_ignores_previous_prediction_and_clears_normalized_checksum(tmp_path):
    from ontseq_platform.models import ModuleRunStatus
    from ontseq_platform.pipeline.marlin import load_marlin_report, with_marlin_result
    from ontseq_platform.pipeline.stages import StageId

    ctx, result = completed_context(tmp_path)
    record = ctx.stage_records[StageId.MARLIN]
    ctx.stage_records[StageId.MARLIN] = record.model_copy(
        update={"status": ModuleRunStatus.FAILED, "reason": "Worker failed", "outputs": []}
    )
    # A previous file (or in-memory artifact list) must not become current evidence.
    result.provenance.reference_checksums["marlin_report"] = "b" * 64
    assert load_marlin_report(ctx).top_class is None
    updated = with_marlin_result(ctx, result)
    assert (
        next(item for item in updated.modules if item.module == AnalysisModule.MARLIN).status
        == ModuleRunStatus.FAILED
    )
    assert "marlin_report" not in updated.provenance.reference_checksums


def test_resumed_current_report_preserves_evidence_and_stable_reason(tmp_path):
    from ontseq_platform.pipeline.marlin import load_marlin_report, with_marlin_result
    from ontseq_platform.pipeline.stages import StageId

    ctx, result = completed_context(tmp_path)
    record = ctx.stage_records[StageId.MARLIN]
    ctx.stage_records[StageId.MARLIN] = record.model_copy(
        update={
            "resumed": True,
            "reason": "Resumed unchanged from a previous run. " + record.reason,
        }
    )
    updated = with_marlin_result(ctx, result)
    assert updated.provenance.reference_checksums["marlin_report"] == record.outputs[0].sha256
    assert updated.modules[-1].reason == load_marlin_report(ctx).reason


def test_tampered_or_missing_current_marlin_evidence_rejected(tmp_path):
    import pytest

    from ontseq_platform.pipeline.marlin import MARLIN_REPORT, load_marlin_report

    ctx, _ = completed_context(tmp_path)
    path = ctx.envelope.path(ctx.path(MARLIN_REPORT))
    path.write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        load_marlin_report(ctx)
    path.unlink()
    with pytest.raises(ValueError, match="missing"):
        load_marlin_report(ctx)


def test_marlin_resume_plan_rehashes_inputs_even_if_metadata_is_restored(tmp_path):
    import hashlib
    import os
    from types import SimpleNamespace
    from unittest.mock import patch

    from ontseq_platform.marlin_native_contracts import NativeMarlinReadiness
    from ontseq_platform.pipeline.marlin import marlin_plan

    ctx, _ = completed_context(tmp_path)
    bam = tmp_path / "synthetic.bam"
    index = tmp_path / "synthetic.bam.bai"
    reference = tmp_path / "synthetic.fa"
    reference_fai = Path(str(reference) + ".fai")
    installation = tmp_path / "installation.json"
    for path in (bam, index, reference, reference_fai, installation):
        path.write_bytes(b"first")
    ctx.manifest.input.path = str(bam)
    ctx.manifest.input.index_path = str(index)
    ctx.config.reference_fasta = reference
    ctx.config.marlin_installation = installation
    with (
        patch(
            "ontseq_platform.marlin_native.check_native_marlin_readiness",
            return_value=NativeMarlinReadiness(ready=True, reason="synthetic"),
        ),
        patch(
            "ontseq_platform.marlin_native.native_marlin_signature",
            return_value={"model": "a" * 64},
        ),
        patch(
            "ontseq_platform.pipeline.marlin.identify_modkit_binary",
            return_value=SimpleNamespace(parameters=lambda: {}),
        ),
    ):
        original = dict(marlin_plan(ctx).external_inputs)
        for label, path in [
            ("marlin_bam", bam),
            ("marlin_bam_index", index),
            ("marlin_reference", reference),
            ("marlin_reference_fai", reference_fai),
        ]:
            stamp = path.stat()
            path.write_bytes(b"other")
            os.utime(path, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
            changed = dict(marlin_plan(ctx).external_inputs)
            assert changed[label] == hashlib.sha256(b"other").hexdigest()
            assert changed[label] != original[label]


def test_execution_rejects_reference_that_changed_after_plan_before_publication(tmp_path):
    from unittest.mock import patch

    import pytest
    from test_marlin_report_presentation import native_report

    from ontseq_platform.models import AlignedBamIntakeReport
    from ontseq_platform.pipeline.marlin import MARLIN_REPORT, marlin_execute
    from ontseq_platform.pipeline.runner import INTAKE_REPORT, StageFailure, StagePlan

    ctx, result = completed_context(tmp_path)
    report = native_report(result)
    ctx.config.marlin_installation = tmp_path / "installation.json"
    ctx.config.reference_fasta = tmp_path / "synthetic.fa"
    intake = AlignedBamIntakeReport(
        sample_id=ctx.sample_id,
        reference_id=ctx.manifest.assay.reference_id,
        genome_build=ctx.manifest.assay.genome_build,
        checks=[],
        verdict="PASS",
    )
    ctx.envelope.atomic_write_text(INTAKE_REPORT, intake.model_dump_json())
    artifact = ctx.envelope.path(ctx.path(MARLIN_REPORT))
    artifact.unlink()
    for key in ("bam", "bam_index", "reference", "reference_fai"):
        planned = {
            "marlin_" + name: value.sha256 for name, value in report.input_fingerprints.items()
        }
        planned["marlin_" + key] = "e" * 64
        plan = StagePlan(
            parameters={"requested": True, "installation_signature": {"model": "a" * 64}},
            tool_versions={},
            external_inputs=tuple(planned.items()),
        )
        with (
            patch(
                "ontseq_platform.marlin_native.native_marlin_signature",
                return_value={"model": "a" * 64},
            ),
            patch("ontseq_platform.marlin_native.run_native_marlin", return_value=report),
            pytest.raises(StageFailure, match=key + " differs from its planned"),
        ):
            marlin_execute(ctx, plan)
        assert not artifact.exists()
