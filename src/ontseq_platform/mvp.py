from __future__ import annotations

from .models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    CraminoQCReport,
    ISCNProposal,
    ModuleOutcome,
    ModuleRunStatus,
    PipelineResult,
    Provenance,
    ReviewStatus,
    SampleManifest,
    SnifflesCallReport,
    Verdict,
)


def assemble_aligned_bam_mvp(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    qc_report: CraminoQCReport,
    *,
    pipeline_version: str,
    git_commit: str,
    sniffles_report: SnifflesCallReport | None = None,
) -> PipelineResult:
    if manifest.sample_id != intake.sample_id or manifest.sample_id != qc_report.sample_id:
        raise ValueError("Manifest, intake and QC artifacts must refer to the same sample")
    if intake.verdict == Verdict.FAIL:
        raise ValueError("Cannot assemble an aligned-BAM result after a failed intake gate")
    if sniffles_report is not None:
        if manifest.sample_id != sniffles_report.sample_id:
            raise ValueError("Manifest and Sniffles artifact must refer to the same sample")
        if manifest.assay.genome_build != sniffles_report.genome_build:
            raise ValueError("Manifest and Sniffles artifact use different genome builds")

    requested = set(manifest.analysis.modules)
    modules: list[ModuleOutcome] = []
    for module in AnalysisModule:
        if module == AnalysisModule.QC:
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.COMPLETED,
                    reason="Aligned-BAM intake and descriptive Cramino QC completed",
                    tools=[item for item in (intake.tool, qc_report.tool) if item is not None],
                )
            )
        elif module == AnalysisModule.REPORT:
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.COMPLETED,
                    reason="Structured MVP result is ready for JSON, HTML and XLSX rendering",
                )
            )
        elif module == AnalysisModule.SV and sniffles_report is not None:
            if sniffles_report.status == ModuleRunStatus.COMPLETED:
                reason = (
                    "Sniffles2 candidate SV evidence was normalized; reportability remains "
                    "benchmark_required"
                )
            else:
                reason = (
                    "Sniffles2 produced no candidate passing the technical policy; this NO_CALL "
                    "is not a biological negative"
                )
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=sniffles_report.status,
                    reason=reason,
                    tools=[sniffles_report.tool],
                )
            )
        elif module in {AnalysisModule.CNV, AnalysisModule.SV}:
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.NOT_RUN,
                    reason="Scientific caller selection remains benchmark_required",
                )
            )
        elif module == AnalysisModule.FUSION:
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.NOT_RUN,
                    reason=(
                        "Breakend candidates require gene annotation and fusion-specific "
                        "validation; an SV call is not a fusion assertion"
                    ),
                )
            )
        elif module == AnalysisModule.ISCN:
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.NOT_RUN,
                    reason=(
                        "No validated, cytoband-normalized CNV/SV interpretation is available "
                        "for an ISCN proposal"
                    ),
                )
            )
        elif module in requested:
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.NOT_RUN,
                    reason="Requested module is outside the aligned-BAM MVP",
                )
            )

    reference_checksums: dict[str, str] = {}
    if intake.header:
        reference_checksums["bam_header_contig_signature"] = intake.header.contig_signature_sha256
    if intake.input_fingerprint and intake.input_fingerprint.sha256:
        reference_checksums["input_bam"] = intake.input_fingerprint.sha256
    if intake.index_fingerprint and intake.index_fingerprint.sha256:
        reference_checksums["input_bam_index"] = intake.index_fingerprint.sha256
    if sniffles_report is not None and sniffles_report.vcf_fingerprint.sha256:
        reference_checksums["sniffles_vcf"] = sniffles_report.vcf_fingerprint.sha256

    events = sniffles_report.events if sniffles_report is not None else []
    iscn_warnings = [
        "ISCN was not generated because validated CNV/SV interpretation is unavailable.",
        "Absence of an ISCN proposal is not a biological negative result.",
    ]
    warnings = [
        "Aligned-BAM pipeline outputs remain research-only and require expert review.",
        "CNV and fusion interpretation remain disabled until benchmark acceptance criteria pass.",
        "No output may be used for diagnosis or treatment decisions.",
    ]
    if sniffles_report is None:
        warnings.insert(1, "SV calling was not run in this artifact.")
    else:
        warnings.insert(
            1,
            "Sniffles2 events are unclassified, non-reportable candidates under an unvalidated "
            "technical policy.",
        )
        warnings.extend(sniffles_report.warnings)
        warnings.extend(sniffles_report.limitations)

    return PipelineResult(
        manifest=manifest,
        qc=qc_report.qc,
        events=events,
        iscn=ISCNProposal(
            notation="NOT GENERATED",
            review_status=ReviewStatus.DRAFT,
            warnings=iscn_warnings,
        ),
        provenance=Provenance(
            pipeline_version=pipeline_version,
            git_commit=git_commit,
            tools=[
                item
                for item in (
                    intake.tool,
                    qc_report.tool,
                    sniffles_report.tool if sniffles_report is not None else None,
                )
                if item is not None
            ],
            reference_checksums=reference_checksums,
        ),
        modules=modules,
        warnings=warnings,
        release_status=ReviewStatus.REVIEW_REQUIRED,
    )
