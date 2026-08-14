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
    Verdict,
)


def assemble_aligned_bam_mvp(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    qc_report: CraminoQCReport,
    *,
    pipeline_version: str,
    git_commit: str,
) -> PipelineResult:
    if manifest.sample_id != intake.sample_id or manifest.sample_id != qc_report.sample_id:
        raise ValueError("Manifest, intake and QC artifacts must refer to the same sample")
    if intake.verdict == Verdict.FAIL:
        raise ValueError("Cannot assemble an aligned-BAM result after a failed intake gate")

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
        elif module in {AnalysisModule.CNV, AnalysisModule.SV, AnalysisModule.FUSION}:
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.NOT_RUN,
                    reason="Scientific caller selection remains benchmark_required",
                )
            )
        elif module == AnalysisModule.ISCN:
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.NOT_RUN,
                    reason="No normalized CNV/SV evidence is available for an ISCN proposal",
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

    return PipelineResult(
        manifest=manifest,
        qc=qc_report.qc,
        events=[],
        iscn=ISCNProposal(
            notation="NOT GENERATED",
            review_status=ReviewStatus.DRAFT,
            warnings=[
                "ISCN was not generated because CNV/SV/fusion modules were not run.",
                "Absence of events in this MVP artifact is not a biological negative result.",
            ],
        ),
        provenance=Provenance(
            pipeline_version=pipeline_version,
            git_commit=git_commit,
            tools=[item for item in (intake.tool, qc_report.tool) if item is not None],
            reference_checksums=reference_checksums,
        ),
        modules=modules,
        warnings=[
            "Aligned-BAM MVP: technical intake and descriptive QC only.",
            "CNV, SV and fusion callers remain disabled until benchmark acceptance criteria pass.",
            "No output may be used for diagnosis or treatment decisions.",
        ],
        release_status=ReviewStatus.REVIEW_REQUIRED,
    )
