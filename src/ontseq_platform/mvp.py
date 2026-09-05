from __future__ import annotations

from .methylation import MethylationReport
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
    sv_consensus_report: SvConsensusReport | None = None,
    methylation_report: MethylationReport | None = None,
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
        elif module == AnalysisModule.SV and (
            sv_consensus_report is not None or sniffles_report is not None
        ):
            # Keyed on the evidence that actually reached the result, not on Sniffles alone.
            # A cuteSV-only run is a supported configuration: it writes a consensus and no
            # Sniffles report, and reading only the latter recorded SV as NOT_RUN while the
            # very same result carried the consensus events.
            if sv_consensus_report is not None:
                callers = ", ".join(sv_consensus_report.caller_names) or "no caller"
                if sv_consensus_report.status == ModuleRunStatus.COMPLETED:
                    reason = (
                        f"{sv_consensus_report.input_event_count} normalized call(s) from "
                        f"{callers} were consolidated into "
                        f"{sv_consensus_report.consolidated_event_count} candidate event(s) and "
                        "technically prioritized into high/moderate/low review tiers; clinical "
                        "reportability remains benchmark_required"
                    )
                else:
                    reason = (
                        f"No candidate from {callers} passed the technical policy; this NO_CALL "
                        "is not a biological negative"
                    )
                status = sv_consensus_report.status
            else:
                assert sniffles_report is not None  # the branch condition guarantees it
                if sniffles_report.status == ModuleRunStatus.COMPLETED:
                    reason = (
                        "Sniffles2 candidates were normalized and technically prioritized into "
                        "high/moderate/low review tiers; clinical reportability remains "
                        "benchmark_required"
                    )
                else:
                    reason = (
                        "Sniffles2 produced no candidate passing the technical policy; this "
                        "NO_CALL is not a biological negative"
                    )
                status = sniffles_report.status
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=status,
                    reason=reason,
                    tools=[sniffles_report.tool] if sniffles_report is not None else [],
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
        elif module == AnalysisModule.METHYLATION and module in requested:
            # The lane is inside the MVP, so "outside the aligned-BAM MVP" would be false.
            # A requested module with no report means the stage did not produce one — it
            # failed closed on a missing tag, a version lock or a policy — and the run
            # report carries which. Saying it is out of scope would hide a refusal.
            modules.append(
                ModuleOutcome(
                    module=module,
                    status=ModuleRunStatus.NOT_RUN,
                    reason=(
                        "The methylation module was requested but produced no report; see the "
                        "run report for why the stage did not complete. This is not a finding "
                        "about the sample's methylation"
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
    if methylation_report is not None and methylation_report.bedmethyl_fingerprint.sha256:
        reference_checksums["bedmethyl"] = methylation_report.bedmethyl_fingerprint.sha256

    source_events = (
        sv_consensus_report.events
        if sv_consensus_report is not None
        else sniffles_report.events
        if sniffles_report is not None
        else []
    )
    events = prioritize_sv_events(source_events)
    iscn_warnings = [
        "ISCN was not generated because validated CNV/SV interpretation is unavailable.",
        "Absence of an ISCN proposal is not a biological negative result.",
    ]
    warnings = [
        "Aligned-BAM pipeline outputs remain research-only and require expert review.",
        "CNV and fusion interpretation remain disabled until benchmark acceptance criteria pass.",
        "No output may be used for diagnosis or treatment decisions.",
    ]
    if sv_consensus_report is None and sniffles_report is None:
        warnings.insert(1, "SV calling was not run in this artifact.")
    else:
        warnings.insert(
            1,
            "SV confidence tiers are automated technical prioritization only; all candidates "
            "remain non-reportable until assay-specific validation criteria pass.",
        )
        # Both sources contribute: a cuteSV-only run has a consensus and no Sniffles report,
        # and its caveats belong in the result either way.
        if sniffles_report is not None:
            warnings.extend(sniffles_report.warnings)
            warnings.extend(sniffles_report.limitations)
        if sv_consensus_report is not None:
            warnings.extend(sv_consensus_report.warnings)
            warnings.extend(sv_consensus_report.limitations)
    if methylation_report is not None:
        warnings.extend(methylation_report.warnings)
        warnings.extend(methylation_report.limitations)

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
                    methylation_report.tool if methylation_report is not None else None,
                )
                if item is not None
            ],
            reference_checksums=reference_checksums,
        ),
        modules=modules,
        warnings=warnings,
        release_status=ReviewStatus.REVIEW_REQUIRED,
    )
