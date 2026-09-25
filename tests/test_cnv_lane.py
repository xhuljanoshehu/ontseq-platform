"""The QDNAseq/ACE lane as part of the declared stage graph.

Three properties are load-bearing. Copy-number evidence reaches a result only as an artifact
the current run's CNV stage recorded; the lane is configured per run, never installed into
the process; and a failed or unconfigured lane appears in the result with the reason its
stage recorded.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest
from test_methylation_pipeline_integration import SYNTHETIC_CNV_LANE, _Fixture
from test_qdnaseq_runtime import _empty_cnv_report

from ontseq_platform.cnv.lane import CNV_DIR, CNV_REPORT, LANE_ID, load_current_cnv
from ontseq_platform.models import AnalysisModule, ModuleRunStatus, PipelineResult
from ontseq_platform.pipeline.runner import (
    IMPLEMENTATIONS,
    RunConfiguration,
    RunContext,
    StageFailure,
    StageImplementation,
    StagePlan,
    StageResult,
)
from ontseq_platform.pipeline.stages import SPEC_BY_STAGE, StageId, VerificationStatus

_TABLES = ("synthetic.bins.tsv", "synthetic.segments.tsv", "synthetic.models.tsv")


def _with_cnv(config: RunConfiguration, *, lane: bool = True) -> RunConfiguration:
    modules = [*config.manifest.analysis.modules, AnalysisModule.CNV]
    manifest = config.manifest.model_copy(
        update={"analysis": config.manifest.analysis.model_copy(update={"modules": modules})}
    )
    return replace(config, manifest=manifest, cnv_lane=SYNTHETIC_CNV_LANE if lane else None)


def _producing_cnv_stage() -> StageImplementation:
    """A CNV stage that records a valid synthetic QDNAseq/ACE report and its tables."""

    def execute(ctx: RunContext, plan: StagePlan) -> StageResult:
        report = _empty_cnv_report(
            sample_id=ctx.manifest.sample_id, genome_build=ctx.manifest.assay.genome_build
        )
        outputs = [
            ctx.envelope.atomic_write_text(f"{CNV_DIR}/{name}", "chromosome\tstart\tend\n")
            for name in _TABLES
        ]
        outputs.append(
            ctx.envelope.atomic_write_text(ctx.path(CNV_REPORT), report.model_dump_json())
        )
        return StageResult(status=report.status, reason="Synthetic CNV stage", outputs=outputs)

    return StageImplementation(lambda ctx: StagePlan({"synthetic": True}, {}), execute)


def _failing_cnv_stage() -> StageImplementation:
    def execute(ctx: RunContext, plan: StagePlan) -> StageResult:
        raise StageFailure("synthetic R failure in QDNAseq")

    return StageImplementation(lambda ctx: StagePlan({"synthetic": "failing"}, {}), execute)


def _cnv_module(result: PipelineResult):  # noqa: ANN202 - the module outcome under test
    return next(item for item in result.modules if item.module is AnalysisModule.CNV)


def _run(fixture: _Fixture, config: RunConfiguration, cnv_stage: StageImplementation | None):
    stages = {StageId.CNV: cnv_stage} if cnv_stage is not None else {}
    with patch.dict(IMPLEMENTATIONS, stages):
        return fixture.run(config)


def test_the_cnv_stage_is_declared_not_registered() -> None:
    spec = SPEC_BY_STAGE[StageId.CNV]
    assert spec.verification is VerificationStatus.VERIFIED_WITH_REAL_TOOL
    assert "QDNAseq" in spec.title
    assert StageId.CNV in IMPLEMENTATIONS


def test_current_cnv_evidence_is_merged_with_its_tables(tmp_path: Path) -> None:
    fixture = _Fixture(tmp_path)
    report, _, _ = _run(fixture, _with_cnv(fixture.config), _producing_cnv_stage())
    assert report.record_for(StageId.ASSEMBLE).status is ModuleRunStatus.COMPLETED
    result = fixture.result()
    module = _cnv_module(result)
    assert module.status is ModuleRunStatus.NO_CALL
    assert module.reason.startswith("QDNAseq+ACE multi-bin CNV completed")
    assert {"cnv_bins", "cnv_segments", "ace_models"} <= {
        item.artifact_id for item in result.sidecars
    }


def test_a_failed_rerun_does_not_promote_the_previous_cnv_report(tmp_path: Path) -> None:
    fixture = _Fixture(tmp_path)
    config = _with_cnv(fixture.config)
    _run(fixture, config, _producing_cnv_stage())
    assert fixture.envelope().path(CNV_REPORT.format(sample=config.manifest.sample_id)).is_file()

    report, _, _ = _run(fixture, replace(config, force=True), _failing_cnv_stage())

    assert report.record_for(StageId.CNV).status is ModuleRunStatus.FAILED
    result = fixture.result()
    module = _cnv_module(result)
    assert module.status is ModuleRunStatus.FAILED
    assert "synthetic R failure" in module.reason
    assert not {"cnv_bins", "cnv_segments", "ace_models"} & {
        item.artifact_id for item in result.sidecars
    }
    assert any(
        blocker.reason_code == "CNV_EVIDENCE_UNAVAILABLE"
        for blocker in result.iscn.assessment_blockers
    )


def test_a_leftover_report_is_ignored_when_cnv_is_not_requested(tmp_path: Path) -> None:
    """Regression: a report on disk alone used to become the result's CNV module."""
    fixture = _Fixture(tmp_path)
    envelope = fixture.envelope()
    sample = fixture.config.manifest.sample_id
    leftover = _empty_cnv_report(
        sample_id=sample, genome_build=fixture.config.manifest.assay.genome_build
    )
    envelope.atomic_write_text(CNV_REPORT.format(sample=sample), leftover.model_dump_json())
    for name in _TABLES:
        envelope.atomic_write_text(f"{CNV_DIR}/{name}", "chromosome\tstart\tend\n")

    config = replace(fixture.config, cnv_lane=SYNTHETIC_CNV_LANE)
    report, _, _ = _run(fixture, config, None)

    assert report.record_for(StageId.CNV).status is ModuleRunStatus.NOT_RUN
    module = _cnv_module(fixture.result())
    assert module.status is ModuleRunStatus.NOT_RUN
    assert module.reason == "CNV was not requested in the sample manifest."
    assert not module.reason.startswith("QDNAseq+ACE")


def test_the_lane_belongs_to_one_run_not_to_the_process(tmp_path: Path) -> None:
    specs_before = dict(SPEC_BY_STAGE)
    implementations_before = dict(IMPLEMENTATIONS)
    (tmp_path / "first").mkdir()
    (tmp_path / "second").mkdir()
    first = _Fixture(tmp_path / "first")
    _run(first, _with_cnv(first.config), _producing_cnv_stage())

    second = _Fixture(tmp_path / "second")
    report, _, _ = _run(second, _with_cnv(second.config, lane=False), None)

    record = report.record_for(StageId.CNV)
    assert record.status is ModuleRunStatus.NOT_RUN
    assert "No copy-number lane is configured" in record.reason
    assert "manifest requests" in record.reason
    assert dict(SPEC_BY_STAGE) == specs_before
    assert dict(IMPLEMENTATIONS) == implementations_before


def test_a_changed_current_cnv_artifact_fails_its_checksum(tmp_path: Path) -> None:
    fixture = _Fixture(tmp_path)
    config = _with_cnv(fixture.config)
    report, _, _ = _run(fixture, config, _producing_cnv_stage())
    envelope = fixture.envelope()
    context = RunContext(
        config=config,
        envelope=envelope,
        runner=None,  # type: ignore[arg-type]
        manifest=config.manifest,
        artifacts={
            StageId.CNV: [item.to_artifact() for item in report.record_for(StageId.CNV).outputs]
        },
    )
    assert load_current_cnv(context) is not None
    envelope.path(CNV_REPORT.format(sample=config.manifest.sample_id)).write_text(
        "{}", encoding="utf-8"
    )
    with pytest.raises(StageFailure, match="checksum"):
        load_current_cnv(context)


def test_assembly_and_report_plans_name_the_lane(tmp_path: Path) -> None:
    fixture = _Fixture(tmp_path)
    report, _, _ = _run(fixture, _with_cnv(fixture.config), _producing_cnv_stage())
    assert report.record_for(StageId.REPORT).status is ModuleRunStatus.COMPLETED
    from ontseq_platform.pipeline.runner import _assemble_plan, _report_plan

    context = RunContext(
        config=_with_cnv(fixture.config),
        envelope=fixture.envelope(),
        runner=None,  # type: ignore[arg-type]
        manifest=_with_cnv(fixture.config).manifest,
    )
    assert _assemble_plan(context).parameters["cnv_lane"] == LANE_ID
    assert _report_plan(context).parameters["cnv_visualization"] == LANE_ID
    unconfigured = replace(context, config=fixture.config)
    assert _assemble_plan(unconfigured).parameters["cnv_lane"] is None
    assert _report_plan(unconfigured).parameters["cnv_visualization"] is None
