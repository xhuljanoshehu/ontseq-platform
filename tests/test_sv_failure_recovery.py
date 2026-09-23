"""A partial caller file cannot turn a failed SV stage into a completed module."""

from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest
from test_methylation_pipeline_integration import _Fixture

from ontseq_platform.cnv.extension import _assemble_execute as assemble_cnv
from ontseq_platform.models import (
    AnalysisModule,
    EventType,
    Evidence,
    FileFingerprint,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    SnifflesCallReport,
    SnifflesPolicy,
    ToolRecord,
)
from ontseq_platform.pipeline.runner import (
    IMPLEMENTATIONS,
    SV_REPORT,
    StageFailure,
    StageImplementation,
    StagePlan,
    StageResult,
    _assemble_execute,
)
from ontseq_platform.pipeline.stages import StageId


@pytest.mark.parametrize("assemble", [_assemble_execute, assemble_cnv])
def test_partial_sv_failure_is_failed_in_result_and_resume(tmp_path: Path, assemble) -> None:
    fixture = _Fixture(tmp_path)
    partial = SnifflesCallReport(
        sample_id=fixture.config.manifest.sample_id,
        genome_build=fixture.config.manifest.assay.genome_build,
        status=ModuleRunStatus.COMPLETED,
        policy=SnifflesPolicy(
            profile_id="synthetic", status="technical_defaults_only", note="test"
        ),
        events=[
            GenomicEvent(
                event_id="SYNTHETIC_SV",
                event_type=EventType.DELETION,
                primary=Locus(chromosome="chr1", start=1000, end=1200),
                evidence=[Evidence(caller="Sniffles2", caller_version="2.8.0", support_reads=8)],
            )
        ],
        raw_record_count=1,
        accepted_record_count=1,
        rejected_record_count=0,
        tool=ToolRecord(name="Sniffles2", version="2.8.0"),
        vcf_fingerprint=FileFingerprint(size_bytes=0, sha256="0" * 64),
    )
    attempts = 0

    def failed_sv(ctx, plan):
        nonlocal attempts
        attempts += 1
        ctx.envelope.atomic_write_text(ctx.path(SV_REPORT), partial.model_dump_json())
        raise StageFailure("cuteSV: Cannot allocate memory")

    # Exercise the runner's real status persistence and both production assemblers.
    # Only the external callers/intake are replaced with synthetic evidence.
    original_run = fixture.run

    def run_with_assembler():
        with patch("test_methylation_pipeline_integration._assemble_execute", assemble):
            return original_run()

    with patch.dict(
        IMPLEMENTATIONS,
        {
            StageId.SV: StageImplementation(lambda ctx: StagePlan({}, {}), failed_sv),
        },
    ):
        report, _, _ = run_with_assembler()
        assert report.record_for(StageId.SV).status is ModuleRunStatus.FAILED
        result = fixture.result()
        outcome = next(x for x in result.modules if x.module is AnalysisModule.SV)
        assert outcome.status is ModuleRunStatus.FAILED
        assert "cuteSV" in outcome.reason
        assert not any("SV calling was not run" in message for message in result.warnings)
        assert result.events == []
        assert (
            fixture.envelope()
            .path(SV_REPORT.format(sample=fixture.config.manifest.sample_id))
            .is_file()
        )
        assert (
            next(x for x in result.modules if x.module is AnalysisModule.FUSION).status
            is ModuleRunStatus.NOT_RUN
        )
        # Corrected status must survive content-addressed resume; methylation is reused.
        again, _, methylation_calls = run_with_assembler()
        assert attempts == 2
        assert methylation_calls == 0
        assert again.record_for(StageId.METHYLATION).resumed
        assert (
            next(x for x in fixture.result().modules if x.module is AnalysisModule.SV).status
            is ModuleRunStatus.FAILED
        )


def test_cutesv_parallelism_does_not_change_other_stage_plans(tmp_path: Path) -> None:
    from ontseq_platform.execution import CommandResult
    from ontseq_platform.models import CuteSvPolicy, SvConsensusPolicy
    from ontseq_platform.pipeline.runner import RunContext, _methylation_plan, _sv_plan

    fixture = _Fixture(tmp_path)
    manifest = fixture.config.manifest.model_copy(
        update={
            "analysis": fixture.config.manifest.analysis.model_copy(
                update={
                    "modules": [AnalysisModule.SV, AnalysisModule.METHYLATION],
                }
            ),
        }
    )
    config = replace(
        fixture.config,
        manifest=manifest,
        cutesv_threads=1,
        cutesv_policy=CuteSvPolicy(
            profile_id="synthetic", status="technical_defaults_only", note="test"
        ),
        sv_consensus_policy=SvConsensusPolicy(
            profile_id="synthetic", status="technical_defaults_only", note="test"
        ),
    )

    class Versions:
        def run(self, argv, **kwargs):
            return CommandResult(
                argv, 0, "cuteSV 2.1.3" if argv[0] == "cuteSV" else "modkit 0.6.4", ""
            )

    ctx = RunContext(config, fixture.envelope(), Versions(), manifest)
    sv_before = _sv_plan(ctx)
    meth_before = _methylation_plan(ctx)
    ctx.config = replace(config, cutesv_threads=2)
    assert _sv_plan(ctx).parameters != sv_before.parameters
    assert _methylation_plan(ctx) == meth_before
    assert ctx.config.threads == 4


@pytest.mark.parametrize("threads", [0, -1])
def test_cutesv_invalid_parallelism_is_rejected(tmp_path: Path, threads: int) -> None:
    with pytest.raises(ValueError, match="cutesv_threads"):
        replace(_Fixture(tmp_path).config, cutesv_threads=threads)


def test_cutesv_worker_setting_reaches_command_and_provenance(tmp_path: Path) -> None:
    from test_cutesv_atomic import CuteSvRunner, qualified_cutesv_mock

    from ontseq_platform.models import CuteSvPolicy, SvConsensusPolicy
    from ontseq_platform.pipeline.runner import INTAKE_REPORT, RunContext, _sv_execute, _sv_plan

    fixture = _Fixture(tmp_path)
    manifest = fixture.config.manifest.model_copy(
        update={
            "analysis": fixture.config.manifest.analysis.model_copy(
                update={"modules": [AnalysisModule.SV]}
            ),
        }
    )
    config = replace(
        fixture.config,
        manifest=manifest,
        threads=4,
        cutesv_policy=CuteSvPolicy(
            profile_id="synthetic", status="technical_defaults_only", note="test"
        ),
        sv_consensus_policy=SvConsensusPolicy(
            profile_id="synthetic", status="technical_defaults_only", note="test"
        ),
    )
    runner = CuteSvRunner()
    ctx = RunContext(config, fixture.envelope(), runner, manifest)
    ctx.envelope.atomic_write_text(INTAKE_REPORT, fixture.intake.model_dump_json())
    with qualified_cutesv_mock():
        outcome = _sv_execute(ctx, _sv_plan(ctx))
    assert outcome.status is ModuleRunStatus.COMPLETED
    assert outcome.tools[0].parameters["threads"] == 1
    assert runner.call_argv[runner.call_argv.index("--threads") + 1] == "1"
    assert config.threads == 4


@pytest.mark.parametrize("assemble", [_assemble_execute, assemble_cnv])
@pytest.mark.parametrize("status", [ModuleRunStatus.FAILED, ModuleRunStatus.NOT_RUN])
def test_failed_sv_does_not_parse_corrupt_partial_json(tmp_path: Path, assemble, status) -> None:
    fixture = _Fixture(tmp_path)

    def failed_sv(ctx, plan):
        ctx.envelope.atomic_write_text(ctx.path(SV_REPORT), "{broken partial")
        if status is ModuleRunStatus.FAILED:
            raise StageFailure("cuteSV: Cannot allocate memory")
        return StageResult(status=ModuleRunStatus.NOT_RUN, reason="SV deselected")

    with (
        patch.dict(
            IMPLEMENTATIONS,
            {
                StageId.SV: StageImplementation(lambda ctx: StagePlan({}, {}), failed_sv),
            },
        ),
        patch("test_methylation_pipeline_integration._assemble_execute", assemble),
    ):
        report, _, _ = fixture.run()
    assert report.record_for(StageId.ASSEMBLE).status is ModuleRunStatus.COMPLETED
    assert (
        next(x for x in fixture.result().modules if x.module is AnalysisModule.SV).status is status
    )
