from __future__ import annotations

from .iscn import build_iscn_proposal, iscn_module_outcome, resource_provenance_for_iscn
from .methylation import MethylationReport
from .models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    CraminoQCReport,
    CuteSvCallReport,
    ISCNAssessmentBlocker,
    ISCNSelectionPolicy,
    ModuleOutcome,
    ModuleRunStatus,
    PipelineResult,
    Provenance,
    ResolvedResourceContext,
    ReviewStatus,
    SampleManifest,
    SidecarArtifact,
    SnifflesCallReport,
    SvConsensusReport,
    Verdict,
)
from .sv_evidence import prioritize_sv_events


def assemble_aligned_bam_mvp(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    qc_report: CraminoQCReport,
    *,
    pipeline_version: str,
    git_commit: str,
    sniffles_report: SnifflesCallReport | None = None,
    cutesv_report: CuteSvCallReport | None = None,
    sv_consensus_report: SvConsensusReport | None = None,
    methylation_report: MethylationReport | None = None,
    reference_context: ResolvedResourceContext | None = None,
    sidecars: list[SidecarArtifact] | None = None,
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
    if methylation_report is not None:
        if manifest.sample_id != methylation_report.sample_id:
            raise ValueError("Manifest and methylation artifact must refer to the same sample")
        if manifest.assay.genome_build != methylation_report.genome_build:
            raise ValueError("Manifest and methylation artifact use different genome builds")
    if cutesv_report is not None:
        if manifest.sample_id != cutesv_report.sample_id:
            raise ValueError("Manifest and cuteSV artifact must refer to the same sample")
        if manifest.assay.genome_build != cutesv_report.genome_build:
            raise ValueError("Manifest and cuteSV artifact use different genome builds")
    if sv_consensus_report is not None:
        if manifest.sample_id != sv_consensus_report.sample_id:
            raise ValueError("Manifest and SV consensus artifact must refer to the same sample")
        if manifest.assay.genome_build != sv_consensus_report.genome_build:
            raise ValueError("Manifest and SV consensus artifact use different genome builds")

    requested = set(manifest.analysis.modules)
    source_events = (
        sv_consensus_report.events
        if sv_consensus_report is not None
        else sniffles_report.events
        if sniffles_report is not None
        else cutesv_report.events
        if cutesv_report is not None
        else []
    )
    events = prioritize_sv_events(source_events)
    iscn = build_iscn_proposal(
        events,
        genome_build=manifest.assay.genome_build,
        selection_policy=ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1,
        requested=AnalysisModule.ISCN in requested,
        upstream_evidence_assessed=False,
        resource_provenance=resource_provenance_for_iscn(reference_context),
        policy_parameters={
            "automatic_unvalidated_sv_to_iscn": False,
            "sex_chromosomes_assessed": False,
            "exact_full_chromosome_span_required_for_iscn": True,
        },
        technical_assumptions=[
            "Only checksum-bound cytobands may be rendered.",
            "Unvalidated SV breakpoint pairs remain review evidence outside formal notation.",
        ],
        assessment_blockers=[
            ISCNAssessmentBlocker(
                reason_code="CNV_EVIDENCE_UNAVAILABLE",
                detail=(
                    "the base aligned-BAM assembler has no completed CNV artifact; the CNV "
                    "extension must recompute ISCN after merging copy-number evidence"
                ),
            ),
            *(
                [
                    ISCNAssessmentBlocker(
                        reason_code="QC_FAILED",
                        detail="run QC failed; ISCN proposal assessment is blocked",
                    )
                ]
                if qc_report.qc.verdict == Verdict.FAIL
                else []
            ),
        ],
    )
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
        elif module == AnalysisModule.SV and (
            sv_consensus_report is not None
            or sniffles_report is not None
            or cutesv_report is not None
        ):
            # The outcome follows the evidence that actually reached the result, not one
            # caller's report: a cuteSV-only run produced consolidated events and must not
            # be reported as though SV calling never happened.
            reached = sv_consensus_report or sniffles_report or cutesv_report
            assert reached is not None
            if reached.status == ModuleRunStatus.COMPLETED:
                reason = (
                    "SV candidates were normalized and technically prioritized into "
                    "high/moderate/low review tiers; clinical reportability remains "
                    "benchmark_required"
                )
            else:
                reason = (
                    "SV calling produced no candidate passing the technical policy; this "
                    "NO_CALL is not a biological negative"
                )
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=reached.status,
                    reason=reason,
                    tools=[
                        report.tool
                        for report in (sniffles_report, cutesv_report)
                        if report is not None
                    ],
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
        elif module == AnalysisModule.METHYLATION and methylation_report is not None:
            if methylation_report.status == ModuleRunStatus.COMPLETED:
                measured = sum(
                    item.sites_at_minimum_coverage for item in methylation_report.regions
                )
                reason = (
                    f"modkit aggregated {measured} modified-base site observation(s) into "
                    f"{len(methylation_report.regions)} region row(s); fractions are "
                    "descriptive and no methylation threshold here is validated"
                )
            else:
                reason = (
                    "No modified-base site reached the configured coverage floor; this "
                    "NO_CALL reports an unmeasurable sample, not unmethylated DNA"
                )
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=methylation_report.status,
                    reason=reason,
                    tools=[methylation_report.tool],
                )
            )
        elif module == AnalysisModule.FUSION:
            if sv_consensus_report is None:
                status = ModuleRunStatus.NOT_RUN
                reason = "No annotated SV consensus was available for fusion assessment"
            elif sv_consensus_report.status == ModuleRunStatus.NO_CALL:
                status = ModuleRunStatus.NO_CALL
                reason = (
                    "Fusion candidate assessment received no normalized SV candidate; this is "
                    "not a biological negative result"
                )
            else:
                fusion_evidence_count = sum(
                    event.fusion_evidence is not None for event in sv_consensus_report.events
                )
                knowledge_match_count = sum(
                    event.known_rearrangement is not None for event in sv_consensus_report.events
                )
                status = ModuleRunStatus.COMPLETED
                reason = (
                    "Breakpoint-level fusion candidate assessment completed: "
                    f"{fusion_evidence_count} event(s) carried fusion evidence and "
                    f"{knowledge_match_count} matched a hematology review pattern. "
                    "Candidate assessment is not analytical validation or clinical release."
                )
            modules.append(ModuleOutcome(module=module, status=status, reason=reason))
        elif module == AnalysisModule.ISCN:
            modules.append(iscn_module_outcome(iscn))
        elif module in requested:
            reason = (
                "No current methylation artifact was supplied; consult the run status for "
                "the skipped or failed stage. This is not an unmethylated result."
                if module == AnalysisModule.METHYLATION
                else "Requested module is outside the aligned-BAM MVP"
            )
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.NOT_RUN,
                    reason=reason,
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
    if methylation_report is not None and methylation_report.bedmethyl_fingerprint.sha256:
        reference_checksums["bedmethyl"] = methylation_report.bedmethyl_fingerprint.sha256
    if cutesv_report is not None and cutesv_report.vcf_fingerprint.sha256:
        reference_checksums["cutesv_vcf"] = cutesv_report.vcf_fingerprint.sha256

    warnings = [
        "Aligned-BAM pipeline outputs remain research-only and require expert review.",
        "CNV and fusion interpretation remain disabled until benchmark acceptance criteria pass.",
        "No output may be used for diagnosis or treatment decisions.",
    ]
    if sniffles_report is None and cutesv_report is None and sv_consensus_report is None:
        warnings.insert(1, "SV calling was not run in this artifact.")
    else:
        warnings.insert(
            1,
            "SV confidence tiers are automated technical prioritization only; all candidates "
            "remain non-reportable until assay-specific validation criteria pass.",
        )
        for caller_report in (sniffles_report, cutesv_report, sv_consensus_report):
            if caller_report is not None:
                warnings.extend(caller_report.warnings)
                warnings.extend(caller_report.limitations)
    if methylation_report is not None:
        warnings.extend(methylation_report.warnings)
        warnings.extend(methylation_report.limitations)

    return PipelineResult(
        manifest=manifest,
        qc=qc_report.qc,
        events=events,
        iscn=iscn,
        provenance=Provenance(
            pipeline_version=pipeline_version,
            git_commit=git_commit,
            tools=[
                item
                for item in (
                    intake.tool,
                    qc_report.tool,
                    sniffles_report.tool if sniffles_report is not None else None,
                    methylation_report.tool if methylation_report is not None else None,
                    cutesv_report.tool if cutesv_report is not None else None,
                )
                if item is not None
            ],
            reference_checksums=reference_checksums,
        ),
        **({"reference_context": reference_context} if reference_context is not None else {}),
        sidecars=sidecars or [],
        modules=modules,
        warnings=warnings,
        release_status=ReviewStatus.REVIEW_REQUIRED,
    )
