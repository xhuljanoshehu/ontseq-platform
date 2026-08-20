from __future__ import annotations

import base64
import html
from collections.abc import MutableMapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import cast

from openpyxl import load_workbook

from ..models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    CraminoQCReport,
    ModuleOutcome,
    ModuleRunStatus,
    PipelineResult,
    Provenance,
    SnifflesCallReport,
)
from ..mvp import assemble_aligned_bam_mvp
from ..pipeline import runner as pipeline_runner
from ..pipeline.envelope import Artifact, sha256_file
from ..pipeline.runner import StageImplementation, StagePlan, StageResult
from ..pipeline.stages import SPEC_BY_STAGE, StageId, StageSpec, VerificationStatus
from ..report import render_html
from ..workbook import render_workbook
from .qdnaseq import QDNAseqCallReport, QDNAseqPolicy, run_qdnaseq_ace

CNV_DIR = "evidence/cnv/qdnaseq"
CNV_REPORT = "evidence/cnv/{sample}.qdnaseq.json"


@dataclass(frozen=True)
class QDNAseqExtensionSettings:
    policy: QDNAseqPolicy
    rscript: str = "Rscript"
    script: Path = Path("scripts/run_qdnaseq_ace.R")


_SETTINGS: QDNAseqExtensionSettings | None = None


def _settings() -> QDNAseqExtensionSettings:
    if _SETTINGS is None:
        raise RuntimeError("QDNAseq extension has not been registered")
    return _SETTINGS


def _requested(ctx: pipeline_runner.RunContext) -> bool:
    return AnalysisModule.CNV in set(ctx.manifest.analysis.modules)


def _probe_r_packages(ctx: pipeline_runner.RunContext) -> dict[str, str]:
    settings = _settings()
    r_version = ctx.runner.run([settings.rscript, "--version"], timeout_seconds=60)
    if r_version.returncode != 0:
        raise ValueError(f"Rscript version probe failed: {r_version.stderr.strip()}")
    expression = (
        "cat(as.character(packageVersion('QDNAseq')), '\\n');"
        "cat(as.character(packageVersion('ACE')), '\\n');"
        "cat(as.character(packageVersion('DNAcopy')), '\\n')"
    )
    packages = ctx.runner.run([settings.rscript, "-e", expression], timeout_seconds=60)
    if packages.returncode != 0:
        diagnostic = packages.stderr.strip()[-2000:]
        raise ValueError(f"QDNAseq/ACE R packages are unavailable: {diagnostic}")
    lines = [line.strip() for line in packages.stdout.splitlines() if line.strip()]
    if len(lines) < 3:
        raise ValueError("could not determine QDNAseq, ACE and DNAcopy versions")
    r_text = r_version.stderr.strip() or r_version.stdout.strip()
    return {
        "Rscript": r_text,
        "QDNAseq": lines[0],
        "ACE": lines[1],
        "DNAcopy": lines[2],
    }


def _cnv_plan(ctx: pipeline_runner.RunContext) -> StagePlan:
    if not _requested(ctx):
        return StagePlan(parameters={"requested": False}, tool_versions={})
    settings = _settings()
    if not settings.script.is_file():
        raise ValueError(f"QDNAseq R runner not found: {settings.script}")
    versions = _probe_r_packages(ctx)
    bam = Path(ctx.manifest.input.path)
    return StagePlan(
        parameters={
            "requested": True,
            "profile": settings.policy.profile_id,
            "bin_sizes_kbp": settings.policy.bin_sizes_kbp,
            "primary_bin_size_kbp": settings.policy.primary_bin_size_kbp,
            "ace_penalty": settings.policy.ace_penalty,
            "ploidy_min": settings.policy.ploidy_min,
            "ploidy_max": settings.policy.ploidy_max,
            "ploidy_step": settings.policy.ploidy_step,
            "threads": ctx.config.threads,
        },
        tool_versions=versions,
        external_inputs=(
            (bam.name, sha256_file(bam)),
            (settings.script.name, sha256_file(settings.script)),
        ),
    )


def _record_cnv_outputs(
    ctx: pipeline_runner.RunContext,
    report: QDNAseqCallReport,
) -> list[Artifact]:
    artifacts: list[Artifact] = []
    relative_report = ctx.path(CNV_REPORT)
    artifacts.append(
        ctx.envelope.atomic_write_text(relative_report, report.model_dump_json(indent=2) + "\n")
    )
    for name in report.output_files:
        if name.lower().endswith(".rds"):
            continue
        artifacts.append(ctx.envelope.fingerprint(f"{CNV_DIR}/{name}"))
    return artifacts


def _cnv_execute(ctx: pipeline_runner.RunContext, plan: StagePlan) -> StageResult:
    del plan
    if not _requested(ctx):
        return StageResult(
            status=ModuleRunStatus.NOT_RUN,
            reason="CNV was not requested in the sample manifest.",
        )
    settings = _settings()
    report = run_qdnaseq_ace(
        bam=Path(ctx.manifest.input.path),
        sample_id=ctx.sample_id,
        genome_build=ctx.manifest.assay.genome_build,
        reference_lock=ctx.config.reference_lock,
        policy=settings.policy,
        output_dir=ctx.envelope.path(CNV_DIR),
        script=settings.script,
        runner=ctx.runner,
        rscript=settings.rscript,
        threads=ctx.config.threads,
    )
    artifacts = _record_cnv_outputs(ctx, report)
    bins = ", ".join(str(value) for value in settings.policy.bin_sizes_kbp)
    reason = (
        f"QDNAseq+ACE completed at {bins} kbp; primary "
        f"{report.primary_fit.bin_size_kbp} kbp fit: cellularity "
        f"{report.primary_fit.cellularity:.3f}, ploidy {report.primary_fit.ploidy:.3f}; "
        f"{len(report.events)} normalized CNV event(s)."
    )
    return StageResult(
        status=report.status,
        reason=reason,
        outputs=artifacts,
        tools=report.tools,
        warnings=report.warnings,
        limitations=report.limitations,
    )


def _load_cnv(ctx: pipeline_runner.RunContext) -> QDNAseqCallReport | None:
    path = ctx.envelope.path(ctx.path(CNV_REPORT))
    if not path.is_file():
        return None
    return QDNAseqCallReport.model_validate_json(path.read_text(encoding="utf-8"))


def _assemble_plan(ctx: pipeline_runner.RunContext) -> StagePlan:
    external: list[tuple[str, str]] = []
    for relative in (ctx.path(CNV_REPORT), ctx.path(pipeline_runner.SV_REPORT)):
        path = ctx.envelope.path(relative)
        if path.is_file():
            external.append((Path(relative).name, sha256_file(path)))
    return StagePlan(
        parameters={
            "pipeline_version": ctx.config.pipeline_version,
            "git_commit": ctx.config.git_commit,
            "cnv_extension": "qdnaseq-ace-v1",
        },
        tool_versions={},
        external_inputs=tuple(external),
    )


def _replace_module(
    modules: Sequence[ModuleOutcome],
    replacement: ModuleOutcome,
) -> list[ModuleOutcome]:
    result = [item for item in modules if item.module != replacement.module]
    result.append(replacement)
    return sorted(result, key=lambda item: item.module.value)


def _merge_provenance(base: Provenance, cnv: QDNAseqCallReport) -> Provenance:
    return base.model_copy(update={"tools": [*base.tools, *cnv.tools]})


def _assemble_execute(ctx: pipeline_runner.RunContext, plan: StagePlan) -> StageResult:
    del plan
    intake = AlignedBamIntakeReport.model_validate_json(
        ctx.envelope.path(pipeline_runner.INTAKE_REPORT).read_text(encoding="utf-8")
    )
    qc = CraminoQCReport.model_validate_json(
        ctx.envelope.path(pipeline_runner.QC_REPORT).read_text(encoding="utf-8")
    )
    sv_path = ctx.envelope.path(ctx.path(pipeline_runner.SV_REPORT))
    sniffles = (
        SnifflesCallReport.model_validate_json(sv_path.read_text(encoding="utf-8"))
        if sv_path.is_file()
        else None
    )
    result = assemble_aligned_bam_mvp(
        ctx.manifest,
        intake,
        qc,
        pipeline_version=ctx.config.pipeline_version,
        git_commit=ctx.config.git_commit,
        sniffles_report=sniffles,
    )
    cnv = _load_cnv(ctx)
    if cnv is not None:
        outcome = ModuleOutcome(
            module=AnalysisModule.CNV,
            status=cnv.status,
            reason=(
                f"QDNAseq+ACE multi-bin CNV completed; primary "
                f"{cnv.primary_fit.bin_size_kbp} kbp, cellularity "
                f"{cnv.primary_fit.cellularity:.3f}, ploidy {cnv.primary_fit.ploidy:.3f}."
            ),
            tools=cnv.tools,
        )
        result = result.model_copy(
            update={
                "events": [*cnv.events, *result.events],
                "modules": _replace_module(result.modules, outcome),
                "provenance": _merge_provenance(result.provenance, cnv),
                "warnings": [*result.warnings, *cnv.warnings, *cnv.limitations],
            }
        )
    artifact = ctx.envelope.atomic_write_text(
        ctx.path(pipeline_runner.RESULT_JSON),
        result.model_dump_json(indent=2) + "\n",
    )
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason="QC, QDNAseq/ACE CNV and available SV evidence assembled into one result.",
        outputs=[artifact],
    )


def _report_plan(ctx: pipeline_runner.RunContext) -> StagePlan:
    external: list[tuple[str, str]] = []
    cnv_path = ctx.envelope.path(ctx.path(CNV_REPORT))
    if cnv_path.is_file():
        external.append((Path(ctx.path(CNV_REPORT)).name, sha256_file(cnv_path)))
    return StagePlan(
        parameters={"formats": ["json", "html", "xlsx"], "cnv_visualization": True},
        tool_versions={},
        external_inputs=tuple(external),
    )


def _cnv_html_section(ctx: pipeline_runner.RunContext, cnv: QDNAseqCallReport) -> str:
    fit_rows = "".join(
        "<tr>"
        f"<td>{fit.bin_size_kbp}</td><td>{fit.cellularity:.3f}</td>"
        f"<td>{fit.ploidy:.3f}</td><td>{fit.fit_error:.6g}</td>"
        f"<td>{fit.segment_count}</td></tr>"
        for fit in sorted(cnv.fits, key=lambda value: value.bin_size_kbp)
    )
    chromosome_rows = "".join(
        "<tr>"
        f"<td>{html.escape(chromosome.chromosome)}</td>"
        f"<td>{chromosome.median_copy_number:.3f}</td>"
        f"<td>{chromosome.rounded_copy_number}</td>"
        f"<td>{chromosome.agreeing_bins}/{chromosome.contributing_bins}</td>"
        f"<td>{chromosome.min_copy_number:.3f}–{chromosome.max_copy_number:.3f}</td>"
        "</tr>"
        for chromosome in cnv.chromosome_consensus
    )
    images: list[str] = []
    for label, name in (
        ("ACE purity/ploidy fit landscape", cnv.primary_fit.fit_plot),
        ("Absolute copy-number profile", cnv.primary_fit.copy_number_plot),
    ):
        path = ctx.envelope.path(f"{CNV_DIR}/{name}")
        if not path.is_file():
            continue
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        images.append(
            f"<h3>{html.escape(label)}</h3>"
            f"<img alt='{html.escape(label)}' style='max-width:100%;height:auto' "
            f"src='data:image/png;base64,{encoded}'>"
        )
    return (
        "<section><h2>Copy-number analysis — QDNAseq + ACE</h2>"
        f"<p><strong>Primary:</strong> {cnv.primary_fit.bin_size_kbp} kbp; "
        f"cellularity {cnv.primary_fit.cellularity:.3f}; "
        f"ploidy {cnv.primary_fit.ploidy:.3f}; "
        f"fit error {cnv.primary_fit.fit_error:.6g}.</p>"
        "<h3>Multi-resolution fits</h3><table><thead><tr><th>Bin (kbp)</th>"
        "<th>Cellularity</th><th>Ploidy</th><th>Fit error</th><th>Segments</th>"
        f"</tr></thead><tbody>{fit_rows}</tbody></table>"
        "<h3>Chromosome-level consensus</h3><table><thead><tr><th>Chromosome</th>"
        "<th>Median CN</th><th>Rounded CN</th><th>Agreement</th><th>Range</th>"
        f"</tr></thead><tbody>{chromosome_rows}</tbody></table>" + "".join(images) + "</section>"
    )


def _enrich_workbook(
    path: Path,
    ctx: pipeline_runner.RunContext,
    cnv: QDNAseqCallReport,
) -> None:
    workbook = load_workbook(path)
    for title in ("CNV Fits", "CNV Consensus", "CNV Segments"):
        if title in workbook.sheetnames:
            del workbook[title]

    fits_sheet = workbook.create_sheet("CNV Fits")
    fits_sheet.append(
        [
            "bin_size_kbp",
            "cellularity",
            "ploidy",
            "fit_error",
            "candidate_count",
            "segment_count",
        ]
    )
    for fit in sorted(cnv.fits, key=lambda value: value.bin_size_kbp):
        fits_sheet.append(
            [
                fit.bin_size_kbp,
                fit.cellularity,
                fit.ploidy,
                fit.fit_error,
                fit.candidate_count,
                fit.segment_count,
            ]
        )

    consensus_sheet = workbook.create_sheet("CNV Consensus")
    consensus_sheet.append(
        [
            "chromosome",
            "median_copy_number",
            "rounded_copy_number",
            "agreeing_bins",
            "contributing_bins",
            "min_copy_number",
            "max_copy_number",
        ]
    )
    for chromosome in cnv.chromosome_consensus:
        consensus_sheet.append(
            [
                chromosome.chromosome,
                chromosome.median_copy_number,
                chromosome.rounded_copy_number,
                chromosome.agreeing_bins,
                chromosome.contributing_bins,
                chromosome.min_copy_number,
                chromosome.max_copy_number,
            ]
        )

    segment_sheet = workbook.create_sheet("CNV Segments")
    segment_path = ctx.envelope.path(f"{CNV_DIR}/{cnv.primary_fit.segment_file}")
    with segment_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            segment_sheet.append(line.rstrip("\n").split("\t"))
    workbook.save(path)


def _report_execute(ctx: pipeline_runner.RunContext, plan: StagePlan) -> StageResult:
    del plan
    result = PipelineResult.model_validate_json(
        ctx.envelope.path(ctx.path(pipeline_runner.RESULT_JSON)).read_text(encoding="utf-8")
    )
    html_path = ctx.envelope.path(ctx.path(pipeline_runner.REPORT_HTML))
    xlsx_path = ctx.envelope.path(ctx.path(pipeline_runner.REPORT_XLSX))
    render_html(result, html_path)
    render_workbook(result, xlsx_path)
    cnv = _load_cnv(ctx)
    if cnv is not None:
        document = html_path.read_text(encoding="utf-8")
        section = _cnv_html_section(ctx, cnv)
        marker = "<section><h2>Warnings and limitations</h2>"
        if marker in document:
            document = document.replace(marker, section + marker, 1)
        else:
            document = document.replace("</main>", section + "</main>", 1)
        html_path.write_text(document, encoding="utf-8")
        _enrich_workbook(xlsx_path, ctx, cnv)
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason="Reviewer HTML and Excel rendered with integrated QDNAseq/ACE CNV summaries.",
        outputs=[
            ctx.envelope.fingerprint(ctx.path(pipeline_runner.REPORT_HTML)),
            ctx.envelope.fingerprint(ctx.path(pipeline_runner.REPORT_XLSX)),
        ],
    )


def register_qdnaseq_extension(settings: QDNAseqExtensionSettings) -> None:
    """Install the live CNV stage into the existing execution graph for this process."""
    global _SETTINGS
    _SETTINGS = settings
    specs = cast(MutableMapping[StageId, StageSpec], SPEC_BY_STAGE)
    current = specs[StageId.CNV]
    specs[StageId.CNV] = replace(
        current,
        title="QDNAseq + ACE copy-number analysis",
        verification=VerificationStatus.UNVERIFIED_ADAPTER,
        purpose=(
            "Run multi-resolution QDNAseq read-depth correction and CBS segmentation, "
            "estimate purity/ploidy with ACE, and retain consensus plus plots."
        ),
    )
    pipeline_runner.IMPLEMENTATIONS[StageId.CNV] = StageImplementation(_cnv_plan, _cnv_execute)
    pipeline_runner.IMPLEMENTATIONS[StageId.ASSEMBLE] = StageImplementation(
        _assemble_plan,
        _assemble_execute,
    )
    pipeline_runner.IMPLEMENTATIONS[StageId.REPORT] = StageImplementation(
        _report_plan,
        _report_execute,
    )
