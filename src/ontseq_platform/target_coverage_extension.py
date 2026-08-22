from __future__ import annotations

import re
import shutil
from collections.abc import MutableMapping
from dataclasses import dataclass, replace
from pathlib import Path
from typing import cast

from .models import AlignedBamIntakeReport, AssayMode, ModuleRunStatus
from .pipeline import runner as pipeline_runner
from .pipeline.envelope import sha256_file
from .pipeline.runner import StageImplementation, StagePlan, StageResult
from .pipeline.stages import SPEC_BY_STAGE, StageId, StageSpec, VerificationStatus
from .target_coverage import TargetCoveragePolicy, run_target_coverage

TARGET_COVERAGE_REPORT = "qc/{sample}.target-coverage.json"
TARGET_COVERAGE_WORK = "work/{sample}.target-coverage"
_VERSION = re.compile(r"(?<!\d)(\d+\.\d+\.\d+)(?!\d)")


@dataclass(frozen=True)
class TargetCoverageExtensionSettings:
    policy: TargetCoveragePolicy
    mosdepth: str = "mosdepth"


_SETTINGS: TargetCoverageExtensionSettings | None = None


def _settings() -> TargetCoverageExtensionSettings:
    if _SETTINGS is None:
        raise RuntimeError("target-coverage extension has not been registered")
    return _SETTINGS


def _enabled(ctx: pipeline_runner.RunContext) -> bool:
    return ctx.manifest.assay.mode == AssayMode.ADAPTIVE_SAMPLING


def _probe_mosdepth(ctx: pipeline_runner.RunContext, executable: str) -> str:
    result = ctx.runner.run([executable, "--version"], timeout_seconds=30)
    if result.returncode != 0:
        raise ValueError(f"Mosdepth version probe failed with exit code {result.returncode}")
    match = _VERSION.search(f"{result.stdout}\n{result.stderr}")
    if match is None:
        raise ValueError("could not determine the Mosdepth version")
    return match.group(1)


def _target_coverage_plan(ctx: pipeline_runner.RunContext) -> StagePlan:
    if not _enabled(ctx):
        return StagePlan(parameters={"adaptive_sampling": False}, tool_versions={})

    settings = _settings()
    target_bed_text = ctx.manifest.assay.target_bed
    target_bed_version = ctx.manifest.assay.target_bed_version
    if not target_bed_text or not target_bed_version:
        raise ValueError("Adaptive Sampling requires target BED metadata for observability")
    target_bed = Path(target_bed_text)
    bam = Path(ctx.manifest.input.path)
    if not target_bed.is_file():
        raise ValueError(f"Adaptive Sampling target BED is missing: {target_bed.name}")
    if not bam.is_file():
        raise ValueError(f"aligned BAM is missing before target coverage: {bam.name}")

    version = _probe_mosdepth(ctx, settings.mosdepth)
    if version != settings.policy.expected_version:
        raise ValueError(
            f"Mosdepth version {version!r} does not match policy lock "
            f"{settings.policy.expected_version!r}"
        )

    return StagePlan(
        parameters={
            "adaptive_sampling": True,
            "profile": settings.policy.profile_id,
            "policy_status": settings.policy.status,
            "target_bed_version": target_bed_version,
            "target_bed_role": "analysis_roi_unbuffered",
            "thresholds": settings.policy.thresholds,
            "mapq": settings.policy.mapq,
            "exclude_flags": settings.policy.exclude_flags,
            "threads": ctx.config.threads,
        },
        tool_versions={"mosdepth": version},
        external_inputs=(
            (target_bed.name, sha256_file(target_bed)),
            (bam.name, sha256_file(bam)),
        ),
    )


def _target_coverage_execute(
    ctx: pipeline_runner.RunContext,
    plan: StagePlan,
) -> StageResult:
    del plan
    if not _enabled(ctx):
        return StageResult(
            status=ModuleRunStatus.NOT_RUN,
            reason=(
                "Target coverage is not applicable to this lcWGS run; no Adaptive Sampling "
                "target design was declared."
            ),
        )

    settings = _settings()
    intake = AlignedBamIntakeReport.model_validate_json(
        ctx.envelope.path(pipeline_runner.INTAKE_REPORT).read_text(encoding="utf-8")
    )
    work_relative = ctx.path(TARGET_COVERAGE_WORK)
    work_dir = ctx.envelope.path(work_relative)
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    report = run_target_coverage(
        ctx.manifest,
        intake,
        settings.policy,
        output_dir=work_dir,
        runner=ctx.runner,
        mosdepth=settings.mosdepth,
        threads=ctx.config.threads,
    )
    report_artifact = ctx.envelope.atomic_write_text(
        ctx.path(TARGET_COVERAGE_REPORT),
        report.model_dump_json(indent=2) + "\n",
    )
    region_count = int(report.summary_metrics.get("region_count", len(report.regions)))
    bins = "/".join(f"{value}x" for value in settings.policy.thresholds)
    return StageResult(
        status=report.status,
        reason=(
            f"Normalized Adaptive Sampling observability for {region_count} target region(s). "
            f"Technical {bins} depth bins are descriptive only and are not assay-adequacy, "
            "reportability or biological-negative thresholds."
        ),
        outputs=[report_artifact],
        tools=[report.tool],
        warnings=report.warnings,
        limitations=report.limitations,
    )


def register_target_coverage_extension(settings: TargetCoverageExtensionSettings) -> None:
    """Install Mosdepth target observability into the end-to-end execution graph.

    Registration changes only the implementation/engineering-verification state of the
    pre-existing ``target_coverage`` stage. It does not add a clinical threshold or alter
    CNV/SV/fusion/ISCN interpretation.
    """

    global _SETTINGS
    _SETTINGS = settings
    specs = cast(MutableMapping[StageId, StageSpec], SPEC_BY_STAGE)
    current = specs[StageId.TARGET_COVERAGE]
    specs[StageId.TARGET_COVERAGE] = replace(
        current,
        title="Mosdepth Adaptive Sampling observability",
        verification=VerificationStatus.VERIFIED_WITH_REAL_TOOL,
        purpose=(
            "Measure per-target depth and configured technical coverage fractions for "
            "Adaptive Sampling. These values describe observability and do not establish "
            "biological negativity or clinical adequacy."
        ),
    )
    pipeline_runner.IMPLEMENTATIONS[StageId.TARGET_COVERAGE] = StageImplementation(
        _target_coverage_plan,
        _target_coverage_execute,
    )
