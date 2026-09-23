"""Current-run integration of the native MARLIN research adapter.

Blocked stages have no output artifacts. Their typed presentation is reconstructed only
from the current run record, never from a leftover file from an earlier attempt.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import TYPE_CHECKING

from ..marlin_native_contracts import NativeMarlinReport
from ..models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    ModuleOutcome,
    ModuleRunStatus,
    PipelineResult,
)
from ..modkit_build import identify_modkit_binary
from .stages import StageId

if TYPE_CHECKING:
    from .runner import RunContext, StagePlan, StageResult

MARLIN_REPORT = "evidence/marlin/{sample}.marlin.json"
MARLIN_WORK = "work/marlin"


def requested(ctx: RunContext) -> bool:
    return AnalysisModule.MARLIN in ctx.manifest.analysis.modules


def marlin_plan(ctx: RunContext) -> StagePlan:
    from ..marlin_native import (
        check_native_marlin_readiness,
        native_marlin_signature,
        selected_marlin_bam_index,
    )
    from .runner import StageFailure, StagePlan, _stable_digest

    if not requested(ctx):
        return StagePlan(parameters={"requested": False}, tool_versions={})
    installation = ctx.config.marlin_installation
    ready = check_native_marlin_readiness(installation, ctx.manifest.assay.genome_build)
    if not ready.ready:
        return StagePlan(parameters={"requested": True, "blocked": ready.reason}, tool_versions={})
    if ctx.config.reference_fasta is None:
        return StagePlan(
            parameters={
                "requested": True,
                "blocked": "MARLIN requires the locked reference FASTA.",
            },
            tool_versions={},
        )
    assert installation is not None
    signature = native_marlin_signature(installation, ctx.manifest.assay.genome_build)
    binary = identify_modkit_binary(ctx.config.executable("modkit"))
    inputs = []
    for label, path in (
        ("marlin_bam", Path(ctx.manifest.input.path)),
        ("marlin_bam_index", selected_marlin_bam_index(ctx.manifest)),
        ("marlin_reference", ctx.config.reference_fasta),
        ("marlin_reference_fai", Path(str(ctx.config.reference_fasta) + ".fai")),
    ):
        # Resume must see current bytes, even when size and modification time were preserved.
        digest, stable = _stable_digest(path)
        if not stable:
            raise StageFailure(f"{label} changed while MARLIN inputs were fingerprinted")
        inputs.append((label, digest))
    return StagePlan(
        parameters={
            "requested": True,
            "installation_signature": signature,
            "threads": ctx.config.threads,
            **binary.parameters(),
        },
        tool_versions={"MARLIN": "1.0.0"},
        external_inputs=(
            *inputs,
            ctx.fingerprint_external_input(installation, label="marlin_installation"),
        ),
    )


def marlin_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    from ..marlin_native import native_marlin_signature, run_native_marlin
    from .runner import INTAKE_REPORT, StageFailure, StageResult

    if plan.parameters.get("requested") is False:
        return StageResult(
            status=ModuleRunStatus.NOT_RUN,
            reason="MARLIN was not requested; no methylation classification was performed.",
        )
    if plan.parameters.get("blocked"):
        return StageResult(status=ModuleRunStatus.NOT_RUN, reason=str(plan.parameters["blocked"]))
    installation = ctx.config.marlin_installation
    reference = ctx.config.reference_fasta
    assert installation is not None and reference is not None
    if (
        native_marlin_signature(installation, ctx.manifest.assay.genome_build)
        != plan.parameters["installation_signature"]
    ):
        raise StageFailure("MARLIN installation changed after stage planning")
    intake = AlignedBamIntakeReport.model_validate_json(
        ctx.envelope.path(INTAKE_REPORT).read_text()
    )
    work = ctx.envelope.path(MARLIN_WORK)
    if work.exists():
        shutil.rmtree(work)
    report = run_native_marlin(
        run_id=ctx.config.run_id,
        manifest=ctx.manifest,
        intake=intake,
        reference_fasta=reference,
        installation_path=installation,
        output_dir=work,
        runner=ctx.runner,
        modkit=ctx.config.executable("modkit"),
        samtools=ctx.config.executable("samtools"),
        threads=ctx.config.threads,
    )
    _check_report_identity(ctx, report)
    if report.status not in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL}:
        return StageResult(
            status=report.status,
            reason=report.reason,
            tools=report.tools,
            warnings=report.warnings,
            limitations=report.limitations,
        )
    if (
        native_marlin_signature(installation, ctx.manifest.assay.genome_build)
        != plan.parameters["installation_signature"]
    ):
        raise StageFailure("MARLIN installation changed during execution")
    planned_inputs = dict(plan.external_inputs)
    for key in ("bam", "bam_index", "reference", "reference_fai"):
        actual = report.input_fingerprints.get(key)
        if actual is None or actual.sha256 != planned_inputs.get("marlin_" + key):
            raise StageFailure(f"MARLIN {key} differs from its planned input fingerprint")
    artifact = ctx.envelope.atomic_write_text(
        ctx.path(MARLIN_REPORT), report.model_dump_json(indent=2) + "\n"
    )
    outputs = [artifact]
    for path in sorted(work.rglob("*")):
        if path.is_file():
            outputs.append(ctx.envelope.fingerprint(path.relative_to(ctx.envelope.root).as_posix()))
    return StageResult(
        status=report.status,
        reason=report.reason,
        outputs=outputs,
        tools=report.tools,
        warnings=report.warnings,
        limitations=report.limitations,
    )


def _check_report_identity(ctx: RunContext, report: NativeMarlinReport) -> None:
    if (report.run_id, report.sample_id, report.genome_build) != (
        ctx.config.run_id,
        ctx.sample_id,
        ctx.manifest.assay.genome_build,
    ):
        raise ValueError("MARLIN report does not belong to this run/sample/build")


def load_marlin_report(ctx: RunContext) -> NativeMarlinReport | None:
    record = ctx.stage_records.get(StageId.MARLIN)
    if record is None or not requested(ctx):
        return None
    if record.status in {ModuleRunStatus.FAILED, ModuleRunStatus.NOT_RUN}:
        return NativeMarlinReport(
            run_id=ctx.config.run_id,
            sample_id=ctx.sample_id,
            genome_build=ctx.manifest.assay.genome_build,
            status=record.status,
            reason=record.reason,
            tools=record.tools,
            warnings=record.warnings,
            limitations=record.limitations,
        )
    artifacts = ctx.artifacts.get(StageId.MARLIN, [])
    artifact = next(
        (item for item in artifacts if item.relative_path == ctx.path(MARLIN_REPORT)), None
    )
    if artifact is None or ctx.envelope.verify(artifacts):
        raise ValueError("MARLIN current report is missing or its evidence checksum changed")
    report = NativeMarlinReport.model_validate_json(
        ctx.envelope.path(artifact.relative_path).read_text()
    )
    _check_report_identity(ctx, report)
    if report.status != record.status or report.tools != record.tools:
        raise ValueError("MARLIN current report status/tools differ from its stage record")
    return report


def current_marlin_outcome(ctx: RunContext) -> ModuleOutcome | None:
    record = ctx.stage_records.get(StageId.MARLIN)
    if record is None or not requested(ctx):
        return None
    reason = record.reason
    prefix = "Resumed unchanged from a previous run. "
    while reason.startswith(prefix):
        reason = reason.removeprefix(prefix)
    return ModuleOutcome(
        module=AnalysisModule.MARLIN, status=record.status, reason=reason, tools=record.tools
    )


def with_marlin_result(ctx: RunContext, result: PipelineResult) -> PipelineResult:
    outcome = current_marlin_outcome(ctx)
    if outcome is None:
        return result
    report = load_marlin_report(ctx)
    assert report is not None
    modules = [item for item in result.modules if item.module != AnalysisModule.MARLIN]
    modules.append(outcome)
    checksums = dict(result.provenance.reference_checksums)
    checksums.pop("marlin_report", None)
    artifact = next(
        (
            item
            for item in ctx.artifacts.get(StageId.MARLIN, [])
            if item.relative_path == ctx.path(MARLIN_REPORT)
        ),
        None,
    )
    if artifact is not None and outcome.status in {
        ModuleRunStatus.COMPLETED,
        ModuleRunStatus.NO_CALL,
    }:
        checksums["marlin_report"] = artifact.sha256
    provenance = result.provenance.model_copy(
        update={
            "reference_checksums": checksums,
            "tools": [*result.provenance.tools, *outcome.tools],
        }
    )
    data = result.model_dump()
    data.update(
        modules=modules,
        provenance=provenance,
        warnings=[*result.warnings, *report.warnings, *report.limitations],
    )
    return PipelineResult.model_validate(data)
