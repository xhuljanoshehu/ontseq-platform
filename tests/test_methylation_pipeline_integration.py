"""Synthetic checks joining optional methylation to the current result/run contracts."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest

from ontseq_platform import __version__
from ontseq_platform.cnv.extension import _assemble_execute as assemble_with_cnv
from ontseq_platform.demo import build_demo_result
from ontseq_platform.execution import CommandResult
from ontseq_platform.methylation import MethylationPolicy, normalize_methylation
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    CraminoQCReport,
    ISCNSelectionPolicy,
    ModuleRunStatus,
    PipelineResult,
    QCPolicy,
    ReferenceContig,
    ReferenceLock,
    ToolRecord,
    Verdict,
)
from ontseq_platform.pipeline.envelope import RunEnvelope
from ontseq_platform.pipeline.runner import (
    IMPLEMENTATIONS,
    INTAKE_REPORT,
    METHYLATION_REPORT,
    QC_REPORT,
    RESULT_JSON,
    RunConfiguration,
    RunContext,
    StageFailure,
    StageImplementation,
    StagePlan,
    StageResult,
    _assemble_execute,
    _assemble_plan,
    _report_execute,
    _report_plan,
    load_methylation_report,
    run_pipeline,
)
from ontseq_platform.pipeline.stages import StageId


class _ModkitVersionOnly:
    def run(self, argv: list[str], *, timeout_seconds: int = 300) -> CommandResult:
        assert argv == ["modkit", "--version"]
        return CommandResult(argv=argv, returncode=0, stdout="mod_kit 0.6.4", stderr="")

    def run_to_file(
        self, argv: list[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        raise AssertionError("Synthetic integration must not invoke a biological binary")


class _Fixture:
    def __init__(self, root: Path) -> None:
        baseline = build_demo_result()
        self.bam = root / "synthetic.bam"
        self.bam.write_text("synthetic-placeholder-v1", encoding="utf-8")
        self.fasta = root / "synthetic.fa"
        self.fasta.write_text(">chr1\nACGTCG\n", encoding="utf-8")
        manifest = baseline.manifest.model_copy(
            update={
                "input": baseline.manifest.input.model_copy(update={"path": str(self.bam)}),
                "assay": baseline.manifest.assay.model_copy(
                    update={
                        "mode": AssayMode.LOW_COVERAGE_WGS,
                        "target_bed": None,
                        "target_bed_version": None,
                    }
                ),
                "analysis": AnalysisSpec(
                    profile="synthetic-integration",
                    modules=[
                        AnalysisModule.QC,
                        AnalysisModule.METHYLATION,
                        AnalysisModule.ISCN,
                        AnalysisModule.REPORT,
                    ],
                ),
            }
        )
        self.policy = MethylationPolicy(
            profile_id="synthetic-integration",
            status="technical_defaults_only",
            note="Synthetic integration fixture, not real modkit validation",
        )
        self.config = RunConfiguration(
            manifest=manifest,
            reference_lock=ReferenceLock(
                reference_id=manifest.assay.reference_id,
                genome_build=manifest.assay.genome_build,
                contigs=[ReferenceContig(name="chr1", length=6)],
                source_fai_sha256="a" * 64,
            ),
            output_base=root / "runs",
            run_id=manifest.run_id,
            pipeline_version=__version__,
            git_commit="0" * 40,
            qc_policy=QCPolicy(status="technical_defaults_only", note="synthetic"),
            methylation_policy=self.policy,
            reference_fasta=self.fasta,
            resource_context=baseline.reference_context,
        )
        self.intake = AlignedBamIntakeReport(
            sample_id=manifest.sample_id,
            reference_id=manifest.assay.reference_id,
            genome_build=manifest.assay.genome_build,
            checks=[],
            verdict=Verdict.PASS,
            tool=ToolRecord(name="samtools", version="1.24"),
        )
        self.qc = CraminoQCReport(
            sample_id=manifest.sample_id,
            qc=baseline.qc,
            tool=ToolRecord(name="cramino", version="1.3.0"),
        )
        self.bedmethyl = root / "synthetic.bedmethyl"
        self.set_count(8)

    def set_count(self, modified: int) -> None:
        self.bedmethyl.write_text(
            f"chr1\t1\t2\tm\t10\t+\t1\t2\t0,0,0\t10\t{modified * 10}\t"
            f"{modified}\t{10 - modified}\t0\t0\t0\t0\t0\n",
            encoding="utf-8",
        )

    def methylation_report(self):
        return normalize_methylation(
            sample_id=self.config.manifest.sample_id,
            genome_build=self.config.manifest.assay.genome_build,
            bedmethyl_path=self.bedmethyl,
            policy=self.policy,
            tool=ToolRecord(name="modkit", version="0.6.4"),
        )

    def stage(self, relative: str, model) -> StageImplementation:
        def execute(ctx: RunContext, plan: StagePlan) -> StageResult:
            artifact = ctx.envelope.atomic_write_text(relative, model.model_dump_json())
            return StageResult(
                status=ModuleRunStatus.COMPLETED, reason="Synthetic fixture", outputs=[artifact]
            )

        return StageImplementation(lambda ctx: StagePlan(parameters={}, tool_versions={}), execute)

    def run(self, config: RunConfiguration | None = None):
        with (
            patch.dict(
                "ontseq_platform.pipeline.runner.IMPLEMENTATIONS",
                {
                    StageId.INTAKE: self.stage(INTAKE_REPORT, self.intake),
                    StageId.QC: self.stage(QC_REPORT, self.qc),
                    StageId.ASSEMBLE: StageImplementation(_assemble_plan, _assemble_execute),
                    StageId.REPORT: StageImplementation(_report_plan, _report_execute),
                },
            ),
            patch(
                "ontseq_platform.pipeline.runner.run_methylation",
                return_value=self.methylation_report(),
            ) as adapter,
        ):
            report, release = run_pipeline(config or self.config, runner=_ModkitVersionOnly())
        return report, release, adapter.call_count

    def envelope(self) -> RunEnvelope:
        return RunEnvelope.create(
            self.config.output_base,
            run_id=self.config.run_id,
            sample_id=self.config.manifest.sample_id,
        )

    def result(self) -> PipelineResult:
        return PipelineResult.model_validate_json(
            self.envelope()
            .path(RESULT_JSON.format(sample=self.config.manifest.sample_id))
            .read_text(encoding="utf-8")
        )


def test_optional_methylation_reaches_schema_030_and_preserves_resource_and_iscn_contracts(
    tmp_path: Path,
) -> None:
    fixture = _Fixture(tmp_path)
    report, release, calls = fixture.run()
    assert report.passed, report.verdict_reason
    assert release is not None
    assert calls == 1
    assert StageId.METHYLATION not in report.unverified_stages
    result = fixture.result()
    assert result.schema_version == "0.3.0"
    assert result.reference_context == fixture.config.resource_context
    assert result.iscn.selection_policy is ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1
    assert result.iscn.notation is None
    assert result.iscn.requires_expert_review
    assert not result.iscn.clinical_release_allowed
    outcome = next(item for item in result.modules if item.module is AnalysisModule.METHYLATION)
    assert outcome.status is ModuleRunStatus.COMPLETED
    assert [tool.name for tool in outcome.tools] == ["modkit"]
    assert result.provenance.reference_checksums["bedmethyl"]
    assert any(item.relative_path.endswith(".methylation.json") for item in release.artifacts)


def test_resumed_run_does_not_reuse_methylation_after_opt_out(tmp_path: Path) -> None:
    fixture = _Fixture(tmp_path)
    fixture.run()
    manifest = fixture.config.manifest.model_copy(
        update={
            "analysis": fixture.config.manifest.analysis.model_copy(
                update={"modules": [AnalysisModule.QC, AnalysisModule.ISCN, AnalysisModule.REPORT]}
            )
        }
    )
    report, release, calls = fixture.run(replace(fixture.config, manifest=manifest))
    assert report.passed
    assert release is not None
    assert calls == 0
    assert report.record_for(StageId.METHYLATION).status is ModuleRunStatus.NOT_RUN
    assert not report.record_for(StageId.ASSEMBLE).resumed
    result = fixture.result()
    assert all(item.module is not AnalysisModule.METHYLATION for item in result.modules)
    assert "bedmethyl" not in result.provenance.reference_checksums
    assert not any(item.relative_path.endswith(".methylation.json") for item in release.artifacts)


def test_new_methylation_measurement_invalidates_assembly_resume(tmp_path: Path) -> None:
    fixture = _Fixture(tmp_path)
    first, _, _ = fixture.run()
    old_digest = fixture.result().provenance.reference_checksums["bedmethyl"]
    fixture.set_count(3)
    fixture.bam.write_text("synthetic-placeholder-v2", encoding="utf-8")
    second, _, calls = fixture.run()
    assert second.passed
    assert calls == 1
    assert not second.record_for(StageId.METHYLATION).resumed
    assert not second.record_for(StageId.ASSEMBLE).resumed
    assert (
        first.record_for(StageId.ASSEMBLE).signature
        != second.record_for(StageId.ASSEMBLE).signature
    )
    assert fixture.result().provenance.reference_checksums["bedmethyl"] != old_digest


def test_interval_method_upgrade_invalidates_old_methylation_resume(tmp_path: Path) -> None:
    fixture = _Fixture(tmp_path)
    current = IMPLEMENTATIONS[StageId.METHYLATION]

    def legacy_plan(ctx: RunContext) -> StagePlan:
        plan = current.plan(ctx)
        parameters = dict(plan.parameters)
        parameters.pop("ontseq_region_assignment", None)
        return replace(plan, parameters=parameters)

    with patch.dict(
        IMPLEMENTATIONS,
        {StageId.METHYLATION: replace(current, plan=legacy_plan)},
    ):
        first, _, _ = fixture.run()
    upgraded, _, _ = fixture.run()
    assert upgraded.passed
    assert not upgraded.record_for(StageId.METHYLATION).resumed
    assert (
        first.record_for(StageId.METHYLATION).signature
        != upgraded.record_for(StageId.METHYLATION).signature
    )
    unchanged, _, _ = fixture.run()
    assert unchanged.record_for(StageId.METHYLATION).resumed


def test_cnv_assembler_keeps_current_methylation_and_rejects_stale_artifacts(
    tmp_path: Path,
) -> None:
    fixture = _Fixture(tmp_path)
    fixture.run()
    envelope = fixture.envelope()
    relative = METHYLATION_REPORT.format(sample=fixture.config.manifest.sample_id)
    context = RunContext(
        config=fixture.config,
        envelope=envelope,
        runner=_ModkitVersionOnly(),
        manifest=fixture.config.manifest,
        artifacts={StageId.METHYLATION: [envelope.fingerprint(relative)]},
    )
    with patch("ontseq_platform.cnv.extension._load_cnv", return_value=None):
        assemble_with_cnv(context, StagePlan(parameters={}, tool_versions={}))
    result = fixture.result()
    assert result.schema_version == "0.3.0"
    assert result.reference_context == fixture.config.resource_context
    assert any(item.module is AnalysisModule.METHYLATION for item in result.modules)

    context.artifacts[StageId.METHYLATION] = []
    assert load_methylation_report(context) is None
    context.artifacts[StageId.METHYLATION] = [envelope.fingerprint(relative)]
    envelope.path(relative).write_text("changed stale report", encoding="utf-8")
    with pytest.raises(StageFailure, match="checksum"):
        load_methylation_report(context)
