from __future__ import annotations

import base64
import html
import json
import sqlite3
from collections.abc import MutableMapping, Sequence
from contextlib import closing
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import cast

from openpyxl import load_workbook

from ..annotation_cache import require_annotation_cache_build, validate_annotation_cache
from ..iscn import (
    ISCN_RULE_PROFILE,
    build_iscn_proposal,
    iscn_module_outcome,
    resource_provenance_for_iscn,
)
from ..models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    CraminoQCReport,
    CuteSvCallReport,
    EventType,
    ISCNAssessmentBlocker,
    ISCNResourceProvenance,
    ISCNSelectionPolicy,
    ModuleOutcome,
    ModuleRunStatus,
    PipelineResult,
    Provenance,
    SidecarArtifact,
    SnifflesCallReport,
    SvConsensusReport,
    Verdict,
)
from ..mvp import assemble_aligned_bam_mvp
from ..pipeline import runner as pipeline_runner
from ..pipeline.envelope import Artifact, sha256_file
from ..pipeline.runner import StageImplementation, StagePlan, StageResult
from ..pipeline.stages import SPEC_BY_STAGE, StageId, StageSpec, VerificationStatus
from ..qc import read_length_histogram_from_tsv
from ..report import render_html
from ..report_plots import CnvChromosomeBar, ReadLengthBin, cnv_genome_svg
from ..sidecars import tabular_sidecar
from ..target_coverage import TargetCoverageReport
from ..workbook import render_workbook
from .cytoband import AffectedBandGroup, CnvDirection, CnvSegment, Cytoband, annotate_cnv_cytobands
from .qdnaseq import (
    WHOLE_CHROMOSOME_CONFIRMATION,
    QDNAseqCallReport,
    QDNAseqPolicy,
    run_qdnaseq_ace,
)

CNV_DIR = "evidence/cnv/qdnaseq"
CNV_REPORT = "evidence/cnv/{sample}.qdnaseq.json"
CNV_CYTOBAND_REPORT = "evidence/cnv/{sample}.cytobands.json"


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
    external_inputs = [
        ctx.fingerprint_external_input(bam),
        ctx.fingerprint_external_input(settings.script),
    ]
    if ctx.config.annotation_cache is not None:
        external_inputs.append(
            ctx.fingerprint_external_input(ctx.config.annotation_cache, label="annotation_cache")
        )
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
            "cytoband_policy_schema": settings.policy.schema_version,
            "cytoband_policy_id": settings.policy.profile_id,
            "cytoband_affected_fraction": settings.policy.cytoband_affected_fraction,
            "whole_chromosome_fraction": settings.policy.whole_chromosome_fraction,
            "whole_chromosome_span_basis": settings.policy.whole_chromosome_span_basis,
            "whole_chromosome_confirmation": WHOLE_CHROMOSOME_CONFIRMATION,
            "threads": ctx.config.threads,
        },
        tool_versions=versions,
        external_inputs=tuple(external_inputs),
    )


def _load_cytobands(annotation_cache: Path, *, expected_build: str = "GRCh38") -> list[Cytoband]:
    with closing(
        sqlite3.connect(f"file:{annotation_cache.as_posix()}?mode=ro", uri=True)
    ) as connection:
        require_annotation_cache_build(connection, expected_build=expected_build)
        rows = connection.execute(
            "SELECT chrom, start, end, name, gie_stain FROM cytobands "
            "ORDER BY chrom, start, end, name"
        ).fetchall()
    return [
        Cytoband(
            chromosome=str(chromosome),
            start=int(start),
            end=int(end),
            name=str(name),
            gie_stain=str(stain) if stain is not None else None,
        )
        for chromosome, start, end, name, stain in rows
    ]


def _annotate_cnv_cytobands(
    ctx: pipeline_runner.RunContext,
    report: QDNAseqCallReport,
) -> tuple[QDNAseqCallReport, Artifact | None]:
    cache = ctx.config.annotation_cache
    if cache is None or not report.events:
        return report, None
    directions: dict[EventType, CnvDirection] = {
        EventType.CHROMOSOME_GAIN: "gain",
        EventType.DUPLICATION: "gain",
        EventType.CHROMOSOME_LOSS: "loss",
        EventType.DELETION: "loss",
    }
    segments = [
        CnvSegment(
            event_id=event.event_id,
            chromosome=event.primary.chromosome,
            start=event.primary.start,
            end=event.primary.end,
            direction=directions[event.event_type],
            whole_chromosome=event.event_type
            in {EventType.CHROMOSOME_GAIN, EventType.CHROMOSOME_LOSS},
        )
        for event in report.events
        if event.event_type in directions
    ]
    chromosome_sizes = {item.name: item.length for item in ctx.config.reference_lock.contigs}
    annotation = annotate_cnv_cytobands(
        segments,
        _load_cytobands(cache, expected_build=ctx.config.manifest.assay.genome_build.value),
        affected_fraction=_settings().policy.cytoband_affected_fraction,
        chromosome_sizes=chromosome_sizes,
    )
    affected_groups: dict[str, list[AffectedBandGroup]] = {}
    for group in annotation.affected_groups:
        for event_id in group.source_event_ids:
            affected_groups.setdefault(event_id, []).append(group)
    events = []
    for event in report.events:
        unannotated_locus = event.primary.model_copy(
            update={"cytoband_start": None, "cytoband_end": None}
        )
        groups = sorted(
            affected_groups.get(event.event_id, []),
            key=lambda value: (value.start, value.end, value.start_band, value.end_band),
        )
        if not groups:
            events.append(event.model_copy(update={"primary": unannotated_locus}))
            continue
        if len(groups) != 1 or groups[0].source_event_ids != (event.event_id,):
            source_count = len(
                set(event_id for group in groups for event_id in group.source_event_ids)
            )
            events.append(
                event.model_copy(
                    update={
                        "primary": unannotated_locus,
                        "notes": [
                            *event.notes,
                            f"Cytoband projection produced {len(groups)} group(s) spanning "
                            f"{source_count} source event(s); no single-event ISCN band range "
                            "was assigned.",
                        ],
                    }
                )
            )
            continue
        start_band = groups[0].start_band
        end_band = groups[0].end_band
        locus = unannotated_locus.model_copy(
            update={"cytoband_start": start_band, "cytoband_end": end_band}
        )
        events.append(
            event.model_copy(
                update={
                    "primary": locus,
                    "notes": [
                        *event.notes,
                        "Cytoband affected-fraction threshold "
                        f"{annotation.threshold:g}: {start_band}-{end_band}.",
                    ],
                }
            )
        )
    payload = {
        "schema_version": "1.0.0",
        "genome_build": ctx.config.manifest.assay.genome_build.value,
        "policy": {
            "schema_version": _settings().policy.schema_version,
            "profile_id": _settings().policy.profile_id,
            "cytoband_affected_fraction": _settings().policy.cytoband_affected_fraction,
        },
        "threshold": annotation.threshold,
        "raw_overlaps": [asdict(item) for item in annotation.raw_overlaps],
        "affected_groups": [asdict(item) for item in annotation.affected_groups],
        "whole_chromosome_calls": [asdict(item) for item in annotation.whole_chromosome_calls],
    }
    artifact = ctx.envelope.atomic_write_text(
        ctx.path(CNV_CYTOBAND_REPORT), json.dumps(payload, indent=2) + "\n"
    )
    annotated_report = QDNAseqCallReport.model_validate(
        report.model_copy(update={"events": events}).model_dump(mode="python")
    )
    return annotated_report, artifact


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
    report, cytoband_artifact = _annotate_cnv_cytobands(ctx, report)
    artifacts = _record_cnv_outputs(ctx, report)
    if cytoband_artifact is not None:
        artifacts.append(cytoband_artifact)
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


def _verified_iscn_resource_provenance(
    ctx: pipeline_runner.RunContext,
) -> tuple[ISCNResourceProvenance | None, list[ISCNAssessmentBlocker]]:
    context = ctx.config.resource_context
    provenance = resource_provenance_for_iscn(context)
    cache = ctx.config.annotation_cache
    if provenance is None or context is None or cache is None:
        return None, [
            ISCNAssessmentBlocker(
                reason_code="RESOURCE_PROVENANCE_UNAVAILABLE",
                detail="resolved reference/cytoband provenance or annotation cache is unavailable",
            )
        ]
    try:
        configured_cache = cache.resolve(strict=True)
        context_cache = Path(context.resource_paths["reference.annotation_cache"]).resolve(
            strict=True
        )
        reference_lock = Path(context.resource_paths["reference.reference_lock"])
        cytobands = Path(context.resource_paths["reference.cytobands"])
        cache_summary = validate_annotation_cache(
            configured_cache,
            expected_build=ctx.manifest.assay.genome_build.value,
        )
        with closing(
            sqlite3.connect(f"file:{configured_cache.as_posix()}?mode=ro", uri=True)
        ) as connection:
            require_annotation_cache_build(
                connection,
                expected_build=ctx.manifest.assay.genome_build.value,
            )
            metadata = dict(connection.execute("SELECT key, value FROM metadata").fetchall())
        reference_lock_sha256 = sha256_file(reference_lock)
        cytoband_sha256 = sha256_file(cytobands)
    except (KeyError, OSError, sqlite3.DatabaseError, ValueError) as exc:
        return None, [
            ISCNAssessmentBlocker(
                reason_code="RESOURCE_PROVENANCE_UNREADABLE",
                detail=f"ISCN reference/cytoband/cache provenance could not be verified: {exc}",
            )
        ]
    expected = {
        "bundle_id": provenance.reference_bundle_id,
        "bundle_version": provenance.reference_bundle_version,
        "cytoband_release": provenance.cytoband_release,
        "cytoband_sha256": provenance.cytoband_sha256,
    }
    mismatches = []
    if configured_cache != context_cache:
        mismatches.append(
            "configured annotation cache is not the cache selected by the resolved resource context"
        )
    actual_resource_sha256 = {
        "reference_lock_sha256": reference_lock_sha256,
        "cytoband_sha256": cytoband_sha256,
        "annotation_cache_sha256": cache_summary.sha256,
    }
    expected_resource_sha256 = {
        "reference_lock_sha256": provenance.reference_lock_sha256,
        "cytoband_sha256": provenance.cytoband_sha256,
        "annotation_cache_sha256": provenance.annotation_cache_sha256,
    }
    mismatches.extend(
        f"{name}: expected {expected_sha!r}, observed {actual_resource_sha256[name]!r}"
        for name, expected_sha in expected_resource_sha256.items()
        if actual_resource_sha256[name] != expected_sha
    )
    mismatches.extend(
        f"{key}: expected {value!r}, observed {metadata.get(key)!r}"
        for key, value in expected.items()
        if metadata.get(key) != value
    )
    if mismatches:
        return None, [
            ISCNAssessmentBlocker(
                reason_code="RESOURCE_CACHE_PROVENANCE_MISMATCH",
                detail="; ".join(mismatches),
            )
        ]
    return provenance, []


def _assemble_plan(ctx: pipeline_runner.RunContext) -> StagePlan:
    external: list[tuple[str, str]] = []
    for relative in (
        ctx.path(CNV_REPORT),
        ctx.path(pipeline_runner.SV_REPORT),
        ctx.path(pipeline_runner.METHYLATION_REPORT),
        ctx.path(pipeline_runner.SV_CONSENSUS_REPORT),
    ):
        path = ctx.envelope.path(relative)
        if path.is_file():
            external.append((Path(relative).name, sha256_file(path)))
    annotation_cache = ctx.config.annotation_cache
    if annotation_cache is not None and annotation_cache.is_file():
        external.append(("iscn_annotation_cache", sha256_file(annotation_cache)))
    resource_context = ctx.config.resource_context
    if resource_context is not None:
        for label, key in (
            ("iscn_reference_lock", "reference.reference_lock"),
            ("iscn_cytobands", "reference.cytobands"),
        ):
            raw_path = resource_context.resource_paths.get(key)
            resource_path = Path(raw_path) if raw_path is not None else None
            if resource_path is not None and resource_path.is_file():
                external.append((label, sha256_file(resource_path)))
    iscn_reference_lock_sha256 = (
        resource_context.resource_checksums.get("reference.reference_lock", "UNAVAILABLE")
        if resource_context is not None
        else "UNAVAILABLE"
    )
    iscn_cytoband_sha256 = (
        resource_context.resource_checksums.get("reference.cytobands", "UNAVAILABLE")
        if resource_context is not None
        else "UNAVAILABLE"
    )
    iscn_annotation_cache_sha256 = (
        resource_context.resource_checksums.get("reference.annotation_cache", "UNAVAILABLE")
        if resource_context is not None
        else "UNAVAILABLE"
    )
    return StagePlan(
        parameters={
            "pipeline_version": ctx.config.pipeline_version,
            "git_commit": ctx.config.git_commit,
            "cnv_extension": "qdnaseq-ace-v1",
            "iscn_rule_profile": ISCN_RULE_PROFILE,
            "iscn_selection_policy": ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1.value,
            "iscn_exact_full_chromosome_span_required": True,
            "iscn_whole_chromosome_fraction": _settings().policy.whole_chromosome_fraction,
            "iscn_cytoband_affected_fraction": _settings().policy.cytoband_affected_fraction,
            "iscn_reference_lock_sha256": iscn_reference_lock_sha256,
            "iscn_cytoband_sha256": iscn_cytoband_sha256,
            "iscn_annotation_cache_sha256": iscn_annotation_cache_sha256,
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
    cutesv_path = ctx.envelope.path(ctx.path(pipeline_runner.CUTESV_REPORT))
    cutesv = (
        CuteSvCallReport.model_validate_json(cutesv_path.read_text(encoding="utf-8"))
        if cutesv_path.is_file()
        else None
    )
    consensus_path = ctx.envelope.path(ctx.path(pipeline_runner.SV_CONSENSUS_REPORT))
    consensus = (
        SvConsensusReport.model_validate_json(consensus_path.read_text(encoding="utf-8"))
        if consensus_path.is_file()
        else None
    )
    sidecars: list[SidecarArtifact] = []
    histogram = ctx.envelope.path(pipeline_runner.QC_READ_LENGTH_HISTOGRAM)
    if histogram.is_file():
        sidecars.append(
            tabular_sidecar(
                artifact_id="read_length_histogram",
                envelope_root=ctx.envelope.root,
                relative_path=pipeline_runner.QC_READ_LENGTH_HISTOGRAM,
            )
        )
    result = assemble_aligned_bam_mvp(
        ctx.manifest,
        intake,
        qc,
        pipeline_version=ctx.config.pipeline_version,
        git_commit=ctx.config.git_commit,
        sniffles_report=sniffles,
        methylation_report=pipeline_runner.load_methylation_report(ctx),
        cutesv_report=cutesv,
        sv_consensus_report=consensus,
        reference_context=ctx.config.resource_context,
        sidecars=sidecars,
    )
    cnv = _load_cnv(ctx)
    if cnv is not None:
        if cnv.sample_id != ctx.manifest.sample_id:
            raise ValueError("Manifest and QDNAseq/ACE artifact must refer to the same sample")
        if cnv.genome_build != ctx.manifest.assay.genome_build:
            raise ValueError("Manifest and QDNAseq/ACE artifact use different genome builds")
        primary_tables = (
            ("cnv_bins", cnv.primary_fit.bins_file),
            ("cnv_segments", cnv.primary_fit.segment_file),
            ("ace_models", cnv.primary_fit.model_file),
        )
        for artifact_id, filename in primary_tables:
            if filename is None:
                continue
            relative = f"{CNV_DIR}/{filename}"
            sidecars.append(
                tabular_sidecar(
                    artifact_id=artifact_id,
                    envelope_root=ctx.envelope.root,
                    relative_path=relative,
                )
            )
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
        combined_events = [*cnv.events, *result.events]
        iscn_provenance, iscn_blockers = _verified_iscn_resource_provenance(ctx)
        if qc.qc.verdict == Verdict.FAIL:
            iscn_blockers.append(
                ISCNAssessmentBlocker(
                    reason_code="QC_FAILED",
                    detail=(
                        "run QC failed; candidate events remain visible but ISCN assessment "
                        "is blocked"
                    ),
                )
            )
        proposal = build_iscn_proposal(
            combined_events,
            genome_build=ctx.manifest.assay.genome_build,
            selection_policy=ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1,
            requested=AnalysisModule.ISCN in set(ctx.manifest.analysis.modules),
            upstream_evidence_assessed=cnv.status
            in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL},
            resource_provenance=iscn_provenance,
            policy_parameters={
                "automatic_unvalidated_sv_to_iscn": False,
                "sex_chromosomes_assessed": False,
                "exact_full_chromosome_span_required_for_iscn": (
                    _settings().policy.whole_chromosome_span_basis == "exact_contig"
                ),
                "whole_chromosome_span_basis": _settings().policy.whole_chromosome_span_basis,
                "cnv_policy_profile_id": _settings().policy.profile_id,
                "cnv_policy_schema_version": _settings().policy.schema_version,
                "whole_chromosome_fraction": _settings().policy.whole_chromosome_fraction,
                "cytoband_affected_fraction": _settings().policy.cytoband_affected_fraction,
            },
            technical_assumptions=[
                "Whole-chromosome candidate classification uses the versioned QDNAseq "
                "segment-fraction threshold; +chr/-chr rendering additionally requires a "
                "confirmed whole-chromosome span under the versioned span basis, either an "
                "exact zero-to-contig-end segment or full coverage of that chromosome's "
                "assessable QDNAseq bins.",
                "Segmental fragments require one contiguous affected cytoband group at the "
                "versioned affected-fraction threshold.",
                "Unvalidated SV breakpoint pairs remain review evidence outside formal notation.",
            ],
            assessment_blockers=iscn_blockers,
            cnv_source_event_ids={event.event_id for event in cnv.events},
        )
        iscn_outcome = iscn_module_outcome(proposal)
        result = result.model_copy(
            update={
                "events": combined_events,
                "iscn": proposal,
                "modules": _replace_module(_replace_module(result.modules, outcome), iscn_outcome),
                "provenance": _merge_provenance(result.provenance, cnv),
                "warnings": [*result.warnings, *cnv.warnings, *cnv.limitations],
                "sidecars": sidecars,
            }
        )
    result = PipelineResult.model_validate(result.model_dump(mode="python"))
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
    bars = [
        CnvChromosomeBar(
            chromosome=chromosome.chromosome,
            median_cn=chromosome.median_copy_number,
            min_cn=chromosome.min_copy_number,
            max_cn=chromosome.max_copy_number,
            rounded_cn=chromosome.rounded_copy_number,
        )
        for chromosome in cnv.chromosome_consensus
    ]
    genome_svg = cnv_genome_svg(
        bars,
        title="Genome-wide copy-number overview",
        baseline=cnv.primary_fit.ploidy,
    )
    genome_figure = (
        "<h3>Genome-wide copy-number overview</h3><figure class='plot' "
        "style='margin:14px 0'>"
        + genome_svg
        + "<figcaption style='color:#5e687a;font-size:12px;margin-top:6px'>Median copy "
        "number per chromosome from the multi-bin consensus; whiskers show each "
        "chromosome's min–max range across bin sizes; the dashed line marks the "
        "fitted ACE ploidy. Autosomes only — no X/Y statement is made or implied. "
        "Descriptive technical evidence, not an adequacy or reportability assessment."
        "</figcaption></figure>"
        if genome_svg
        else ""
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
        + genome_figure
        + "<h3>Multi-resolution fits</h3><table><thead><tr><th>Bin (kbp)</th>"
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
    target_path = ctx.envelope.path(pipeline_runner.TARGET_COVERAGE_REPORT)
    selection_path = ctx.envelope.path(pipeline_runner.SELECTION_COVERAGE_REPORT)
    target_coverage = (
        TargetCoverageReport.model_validate_json(target_path.read_text(encoding="utf-8"))
        if target_path.is_file()
        else None
    )
    selection_coverage = (
        TargetCoverageReport.model_validate_json(selection_path.read_text(encoding="utf-8"))
        if selection_path.is_file()
        else None
    )
    histogram_path = ctx.envelope.path(pipeline_runner.QC_READ_LENGTH_HISTOGRAM)
    qc_histogram = (
        [
            ReadLengthBin(start=start, end=end, count=count, bases=bases)
            for start, end, count, bases in read_length_histogram_from_tsv(
                histogram_path.read_text(encoding="utf-8")
            )
        ]
        if histogram_path.is_file()
        else None
    )
    render_html(
        result,
        html_path,
        target_coverage=target_coverage,
        selection_coverage=selection_coverage,
        qc_histogram=qc_histogram,
        methylation_report=pipeline_runner.load_methylation_report(ctx),
    )
    render_workbook(
        result,
        xlsx_path,
        target_coverage=target_coverage,
        selection_coverage=selection_coverage,
    )
    cnv = _load_cnv(ctx)
    if cnv is not None:
        document = html_path.read_text(encoding="utf-8")
        section = _cnv_html_section(ctx, cnv)
        marker = "<!-- ONTSEQ_CNV_SECTION -->"
        if marker in document:
            document = document.replace(marker, section, 1)
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
        verification=VerificationStatus.VERIFIED_WITH_REAL_TOOL,
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
