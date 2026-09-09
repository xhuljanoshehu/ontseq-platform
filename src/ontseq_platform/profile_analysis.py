"""Manifest-free, build-isolated profile resolution for ``ontseq analyze``.

This module deliberately stops at a fully resolved :class:`RunConfiguration`.  The existing
pipeline runner remains the sole orchestrator, so the convenience command cannot acquire a
second set of scientific defaults or skip the normal run envelope.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from .bam_resolution import default_run_id, resolve_bam_header, resolve_bam_input
from .io import load_model
from .methylation import MethylationPolicy
from .models import (
    AmlKnowledgeLock,
    AnalysisIntent,
    AnalysisModule,
    AnalysisProfile,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CuteSvPolicy,
    GenomeBuild,
    InputKind,
    InputSpec,
    QCPolicy,
    ReferenceDictionaryContract,
    ReferenceLock,
    SampleManifest,
    SnifflesPolicy,
    SvConsensusPolicy,
    SvEvidencePolicy,
    TargetBedRole,
)
from .pipeline.components import RunComponents
from .pipeline.runner import RunConfiguration
from .reference import (
    GRCH37_UCSC_HG19_CANONICAL_25_FAI_ROLE,
    GRCH37_UCSC_HG19_CANONICAL_25_FASTA_ROLE,
    reference_lock_for_dictionary_contract,
    reference_lock_from_fai,
)
from .resource_registry import ResourceRegistry
from .target_coverage import TargetCoveragePolicy


@dataclass(frozen=True)
class ProfileRuntimeSettings:
    """Trusted runtime policy selection applied on top of a resolved resource profile.

    Resource profiles choose reference, panel and knowledge bundles.  The service still owns
    the policy files and component selection with which it was started; keeping those choices
    in one explicit value prevents the Desktop/profile path from silently reverting to the
    command-line defaults.
    """

    qc_policy: QCPolicy
    sniffles_policy: SnifflesPolicy | None
    cutesv_policy: CuteSvPolicy | None
    sv_consensus_policy: SvConsensusPolicy | None
    sv_evidence_policy: SvEvidencePolicy | None
    target_coverage_policy: TargetCoveragePolicy | None
    sv_minimum_mean_depth: float = 10.0
    components: RunComponents | None = None
    methylation_policy: MethylationPolicy | None = None

    def __post_init__(self) -> None:
        if self.sv_minimum_mean_depth < 0:
            raise ValueError("sv_minimum_mean_depth must be non-negative")


@dataclass(frozen=True)
class AnalyzeSettings:
    bam: Path
    profile_id: str
    resource_root: Path | None
    output_dir: Path = Path("results/runs")
    configuration_root: Path | None = None
    sample_id: str | None = None
    run_id: str | None = None
    pipeline_version: str = "UNKNOWN"
    git_commit: str = "UNKNOWN"
    threads: int = 4
    force: bool = False
    verify_resource_checksums: bool = True
    runtime_settings: ProfileRuntimeSettings | None = None
    include_methylation: bool = False
    executables: Mapping[str, str] = field(
        default_factory=lambda: {
            "samtools": "samtools",
            "cramino": "cramino",
            "sniffles": "sniffles",
            "cutesv": "cuteSV",
            "minimap2": "minimap2",
            "mosdepth": "mosdepth",
            "dorado": "dorado",
            "modkit": "modkit",
        }
    )


def configuration_root(explicit: Path | None = None) -> Path:
    """Locate the versioned runtime configs in a checkout or packed Desktop runtime.

    The operator's working directory is intentionally not a candidate.  Runtime defaults are
    part of the installed ONTSeq release and must therefore resolve to either the wheel/conda
    data directory or the adjacent source checkout, never to an unrelated ``./configs`` tree.
    """

    if explicit is not None:
        root = explicit.expanduser().resolve()
        if not root.is_dir():
            raise FileNotFoundError(f"ONTSeq configuration root not found: {root}")
        return root
    package = Path(__file__).resolve()
    candidates = [
        package.parents[2] / "configs",
        Path(sys.prefix) / "share" / "ontseq" / "configs",
        *(parent / "share" / "ontseq" / "configs" for parent in package.parents),
    ]
    for candidate in candidates:
        if (candidate / "qc" / "defaults.yaml").is_file():
            return candidate.resolve()
    rendered = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"ONTSeq configuration root was not found; checked: {rendered}")


def _required_path(context_paths: Mapping[str, str], key: str) -> Path:
    try:
        return Path(context_paths[key])
    except KeyError as exc:
        raise ValueError(f"resolved profile is missing required resource role {key!r}") from exc


def registry_for_profile(resource_root: Path | None, profile_id: str) -> ResourceRegistry:
    """Select the build declared by an installed profile, never by a BAM or filename.

    Each registry still enforces its own build boundary. Trying discovery in both registries
    does not retry analysis against a different reference after a mismatch.
    """
    matches = [
        registry
        for build in (GenomeBuild.GRCH38, GenomeBuild.GRCH37)
        if profile_id in (registry := ResourceRegistry(resource_root, active_build=build)).profiles
    ]
    if len(matches) != 1:
        raise ValueError(f"analysis profile {profile_id!r} is not uniquely active")
    return matches[0]


def _manifest(
    profile: AnalysisProfile,
    *,
    bam: Path,
    index: Path,
    sample_id: str,
    run_id: str,
    reference_lock: ReferenceLock,
    analysis_roi: Path | None,
    panel_version: str | None,
    include_methylation: bool = False,
) -> SampleManifest:
    assay = AssaySpec(
        mode=profile.assay_mode,
        genome_build=profile.genome_build,
        reference_id=reference_lock.reference_id,
        target_bed=str(analysis_roi) if analysis_roi is not None else None,
        target_bed_version=panel_version,
        target_bed_role=TargetBedRole.ANALYSIS_ROI_UNBUFFERED,
    )
    return SampleManifest(
        sample_id=sample_id,
        run_id=run_id,
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path=str(bam),
            index_path=str(index),
        ),
        assay=assay,
        analysis=AnalysisSpec(
            profile=profile.profile_id,
            modules=[
                AnalysisModule.QC,
                AnalysisModule.CNV,
                AnalysisModule.SV,
                AnalysisModule.FUSION,
                AnalysisModule.ISCN,
                *([AnalysisModule.METHYLATION] if include_methylation else []),
                AnalysisModule.REPORT,
            ],
            intent=AnalysisIntent.SOMATIC,
        ),
    )


def build_profile_run_configuration(
    settings: AnalyzeSettings,
    *,
    header_text: str | None = None,
) -> RunConfiguration:
    """Resolve bundles and reject a build/dictionary mismatch before constructing a run."""

    registry = registry_for_profile(settings.resource_root, settings.profile_id)
    context = registry.resolve_profile(
        settings.profile_id,
        verify_files=settings.verify_resource_checksums,
    )
    profile = registry.profiles[settings.profile_id]
    paths = context.resource_paths
    source_reference_lock = load_model(
        _required_path(paths, "reference.reference_lock"), ReferenceLock
    )
    reference_fasta = _required_path(paths, "reference.genome_fasta")
    bam_reference_lock = source_reference_lock
    if (
        profile.reference_dictionary_contract
        == ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25
    ):
        reference_fasta = _required_path(
            paths,
            f"reference.{GRCH37_UCSC_HG19_CANONICAL_25_FASTA_ROLE}",
        )
        contract_fai = _required_path(
            paths,
            f"reference.{GRCH37_UCSC_HG19_CANONICAL_25_FAI_ROLE}",
        )
        bam_reference_lock = reference_lock_from_fai(
            contract_fai,
            reference_id=source_reference_lock.reference_id,
            genome_build=GenomeBuild.GRCH37,
            allow_extra_contigs=False,
            require_canonical_assembly=True,
        )
    reference_lock = reference_lock_for_dictionary_contract(
        bam_reference_lock,
        profile.reference_dictionary_contract,
    )
    samtools = settings.executables.get("samtools", "samtools")
    if header_text is None:
        resolved_bam = resolve_bam_input(
            settings.bam,
            bam_reference_lock,
            required_build=profile.genome_build,
            dictionary_contract=profile.reference_dictionary_contract,
            sample_id=settings.sample_id,
            samtools=samtools,
        )
    else:
        resolved_bam = resolve_bam_header(
            bam_path=settings.bam,
            header_text=header_text,
            reference_lock=bam_reference_lock,
            required_build=profile.genome_build,
            dictionary_contract=profile.reference_dictionary_contract,
            sample_id=settings.sample_id,
        )
    run_id = settings.run_id or default_run_id(resolved_bam.sample_id)
    analysis_roi: Path | None = None
    selection_bed: Path | None = None
    panel_version: str | None = None
    if profile.assay_mode == AssayMode.ADAPTIVE_SAMPLING:
        analysis_roi = _required_path(paths, "panel.analysis_roi_unbuffered")
        selection_bed = _required_path(paths, "panel.selection_panel_buffered")
        panel_version = context.panel_bundle_version

    manifest = _manifest(
        profile,
        bam=resolved_bam.bam_path,
        index=resolved_bam.index_path,
        sample_id=resolved_bam.sample_id,
        run_id=run_id,
        reference_lock=reference_lock,
        analysis_roi=analysis_roi,
        panel_version=panel_version,
        include_methylation=settings.include_methylation,
    )
    configs = configuration_root(settings.configuration_root)
    knowledge = (
        _required_path(paths, "knowledge.rearrangement_knowledge"),
        load_model(_required_path(paths, "knowledge.knowledge_lock"), AmlKnowledgeLock),
    )
    target_policy_path = configs / "qc" / "adaptive_target_coverage.technical.yaml"
    runtime = settings.runtime_settings
    qc_policy: QCPolicy
    sniffles_policy: SnifflesPolicy | None
    cutesv_policy: CuteSvPolicy | None
    sv_consensus_policy: SvConsensusPolicy | None
    sv_evidence_policy: SvEvidencePolicy | None
    target_policy: TargetCoveragePolicy | None
    components: RunComponents | None
    if runtime is None:
        qc_policy = load_model(configs / "qc" / "defaults.yaml", QCPolicy)
        sniffles_policy = load_model(
            configs / "sv" / "sniffles2.conservative.technical.yaml", SnifflesPolicy
        )
        cutesv_policy = load_model(
            configs / "sv" / "cutesv.conservative.technical.yaml", CuteSvPolicy
        )
        sv_consensus_policy = load_model(
            configs / "sv" / "sniffles2_cutesv.consensus.technical.yaml",
            SvConsensusPolicy,
        )
        sv_evidence_policy = load_model(
            configs / "sv" / "evidence-priority.technical.yaml", SvEvidencePolicy
        )
        target_policy = (
            load_model(target_policy_path, TargetCoveragePolicy)
            if profile.assay_mode == AssayMode.ADAPTIVE_SAMPLING
            else None
        )
        sv_minimum_mean_depth = 10.0
        components = None
    else:
        qc_policy = runtime.qc_policy
        sniffles_policy = runtime.sniffles_policy
        cutesv_policy = runtime.cutesv_policy
        sv_consensus_policy = runtime.sv_consensus_policy
        sv_evidence_policy = runtime.sv_evidence_policy
        target_policy = (
            runtime.target_coverage_policy
            if profile.assay_mode == AssayMode.ADAPTIVE_SAMPLING
            else None
        )
        sv_minimum_mean_depth = runtime.sv_minimum_mean_depth
        components = runtime.components
    methylation_policy = None
    if settings.include_methylation:
        methylation_policy = runtime.methylation_policy if runtime is not None else None
        if methylation_policy is None:
            filename = (
                "modkit-targets.technical.yaml"
                if profile.assay_mode == AssayMode.ADAPTIVE_SAMPLING
                else "modkit.technical.yaml"
            )
            policy_path = configs / "methylation" / filename
            if not policy_path.is_file():
                raise ValueError("Methylation requires the installed versioned modkit policy.")
            methylation_policy = load_model(policy_path, MethylationPolicy)
        if methylation_policy.cpg_only and not reference_fasta.is_file():
            raise ValueError("Methylation requires the profile's locked reference FASTA.")
        if methylation_policy.region_source.value == "target_bed" and analysis_roi is None:
            raise ValueError(
                "Target-region methylation requires the profile's analysis target BED."
            )
    context_resource_paths = {
        role: Path(paths[key])
        for role, key in {
            "blacklist": "reference.blacklist",
            "repeatmasker": "reference.repeatmasker",
            "simple_repeats": "reference.simple_repeats",
            "segmental_duplication": "reference.segmental_duplication",
            "mappability": "reference.mappability",
        }.items()
        if key in paths
    }
    return RunConfiguration(
        manifest=manifest,
        reference_lock=reference_lock,
        output_base=settings.output_dir.expanduser().resolve(),
        run_id=run_id,
        pipeline_version=settings.pipeline_version,
        git_commit=settings.git_commit,
        qc_policy=qc_policy,
        sniffles_policy=sniffles_policy,
        cutesv_policy=cutesv_policy,
        sv_consensus_policy=sv_consensus_policy,
        sv_evidence_policy=sv_evidence_policy,
        aml_knowledge=knowledge,
        sv_minimum_mean_depth=sv_minimum_mean_depth,
        target_coverage_policy=target_policy,
        methylation_policy=methylation_policy,
        reference_fasta=reference_fasta,
        annotation_cache=_required_path(paths, "reference.annotation_cache"),
        selection_target_bed=selection_bed,
        context_resource_paths=context_resource_paths,
        resource_context=context,
        components=components,
        threads=settings.threads,
        executables=settings.executables,
        force=settings.force,
    )


__all__ = [
    "AnalyzeSettings",
    "ProfileRuntimeSettings",
    "build_profile_run_configuration",
    "configuration_root",
]
