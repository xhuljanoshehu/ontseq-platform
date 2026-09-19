from __future__ import annotations

from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from .models import StrictModel

SHA256 = r"^[0-9a-f]{64}$"
ID = r"^[A-Za-z0-9][A-Za-z0-9._-]{1,127}$"


class CallerProvider(StrEnum):
    QDNASEQ_ACE = "qdnaseq_ace"
    SPECTRE = "spectre"
    SNIFFLES2 = "sniffles2"
    CUTESV = "cutesv"
    SEVERUS = "severus"
    SAVANA = "savana"
    WAKHAN = "wakhan"
    ICHORCNA = "ichorcna"


class CallerMode(StrEnum):
    QDNASEQ_ACE_MULTIBIN = "qdnaseq_ace:multibin"
    SPECTRE_DEPTH_ONLY = "spectre:depth_only"
    SPECTRE_SNIFFLES_SUPPORTED = "spectre:sniffles_supported"
    SNIFFLES2_STANDARD = "sniffles2:standard"
    SNIFFLES2_MOSAIC = "sniffles2:mosaic"
    CUTESV_STANDARD = "cutesv:standard"
    SEVERUS_PAIRED = "severus:paired"
    SEVERUS_SINGLE_SAMPLE = "severus:single_sample"
    SAVANA_PAIRED = "savana:paired"
    SAVANA_TUMOR_ONLY = "savana:tumor_only"
    WAKHAN_PHASED_CNA = "wakhan:phased_cna"
    ICHORCNA_ULP_WGS = "ichorcna:ulp_wgs"


class CallerAnalyticalDomain(StrEnum):
    GENOME_WIDE_CNV = "genome_wide_cnv"
    STRUCTURAL_VARIANT = "structural_variant"
    TUMOR_CNA = "tumor_cna"
    CFDNA_CNV = "cfdna_cnv"


class CallerAssayRegime(StrEnum):
    LCWGS = "lcwgs"
    ADAPTIVE_SAMPLING = "adaptive_sampling"
    TUMOR_NORMAL_LONG_READ = "tumor_normal_long_read"
    TUMOR_ONLY_LONG_READ = "tumor_only_long_read"
    CFDNA_ULP_WGS = "cfdna_ulp_wgs"


class CallerInputRole(StrEnum):
    ALIGNED_BAM = "aligned_bam"
    TUMOR_BAM = "tumor_bam"
    MATCHED_NORMAL_BAM = "matched_normal_bam"
    COVERAGE_PROFILE = "coverage_profile"
    SNP_VCF = "snp_vcf"
    BAF_TABLE = "baf_table"
    PHASED_VARIANTS = "phased_variants"
    PANEL_OF_NORMALS = "panel_of_normals"
    SNIFFLES_VCF = "sniffles_vcf"
    SEVERUS_BREAKPOINTS = "severus_breakpoints"


class CallerPlanningDecision(StrEnum):
    ELIGIBLE = "ELIGIBLE"
    INELIGIBLE = "INELIGIBLE"
    RUNTIME_UNAVAILABLE = "RUNTIME_UNAVAILABLE"
    NOT_REQUESTED = "NOT_REQUESTED"


class CallerInputArtifact(StrictModel):
    artifact_id: str = Field(pattern=ID)
    role: CallerInputRole
    sha256: str = Field(pattern=SHA256)
    description: str = ""
    research_only: Literal[True] = True


class CallerLaneRequest(StrictModel):
    lane_id: str = Field(pattern=ID)
    mode_id: CallerMode
    assay_regime: CallerAssayRegime
    inputs: list[CallerInputArtifact] = Field(default_factory=list)
    parent_lane_ids: list[str] = Field(default_factory=list)
    requested: bool = True
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def unique_lane_inputs(self) -> CallerLaneRequest:
        artifact_ids = [item.artifact_id for item in self.inputs]
        input_roles = [item.role for item in self.inputs]
        if len(artifact_ids) != len(set(artifact_ids)):
            raise ValueError("Caller lane input artifact IDs must be unique")
        if len(input_roles) != len(set(input_roles)):
            raise ValueError("Caller lane input roles must be unique")
        if len(self.parent_lane_ids) != len(set(self.parent_lane_ids)):
            raise ValueError("Caller lane parent IDs must be unique")
        if self.lane_id in self.parent_lane_ids:
            raise ValueError("Caller lane cannot depend on itself")
        return self


class CallerLanePlan(StrictModel):
    lane_id: str = Field(pattern=ID)
    mode_id: CallerMode
    decision: CallerPlanningDecision
    reasons: list[str] = Field(default_factory=list)
    input_artifact_ids: list[str] = Field(default_factory=list)
    parent_lane_ids: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def explicit_noneligible_reason(self) -> CallerLanePlan:
        if self.decision == CallerPlanningDecision.ELIGIBLE and self.reasons:
            raise ValueError("Eligible caller lane cannot carry blocker reasons")
        if self.decision != CallerPlanningDecision.ELIGIBLE and not self.reasons:
            raise ValueError("Non-eligible caller lane requires an explicit reason")
        return self


class CallerCatalogEntry(StrictModel):
    mode_id: CallerMode
    provider: CallerProvider
    analytical_domains: list[CallerAnalyticalDomain] = Field(min_length=1)
    compatible_regimes: list[CallerAssayRegime] = Field(min_length=1)
    required_input_roles: list[CallerInputRole] = Field(min_length=1)
    required_parent_modes: list[CallerMode] = Field(default_factory=list)
    rationale: str = Field(min_length=12)
    research_only: bool = True

    @model_validator(mode="after")
    def coherent_catalog_entry(self) -> CallerCatalogEntry:
        if not self.mode_id.value.startswith(f"{self.provider.value}:"):
            raise ValueError("Caller mode ID must use its provider prefix")
        for label, values in (
            ("analytical domains", self.analytical_domains),
            ("compatible regimes", self.compatible_regimes),
            ("required inputs", self.required_input_roles),
            ("required parent modes", self.required_parent_modes),
        ):
            if len(values) != len(set(values)):
                raise ValueError(f"Caller catalog {label} must be unique")
        if self.mode_id in self.required_parent_modes:
            raise ValueError("Caller mode cannot require itself as a parent")
        return self


def _entry(
    mode_id: CallerMode,
    provider: CallerProvider,
    *,
    domains: list[CallerAnalyticalDomain],
    regimes: list[CallerAssayRegime],
    inputs: list[CallerInputRole],
    parents: list[CallerMode] | None = None,
    rationale: str,
) -> CallerCatalogEntry:
    return CallerCatalogEntry(
        mode_id=mode_id,
        provider=provider,
        analytical_domains=domains,
        compatible_regimes=regimes,
        required_input_roles=inputs,
        required_parent_modes=parents or [],
        rationale=rationale,
    )


CALLER_CATALOG: dict[CallerMode, CallerCatalogEntry] = {
    CallerMode.QDNASEQ_ACE_MULTIBIN: _entry(
        CallerMode.QDNASEQ_ACE_MULTIBIN,
        CallerProvider.QDNASEQ_ACE,
        domains=[CallerAnalyticalDomain.GENOME_WIDE_CNV],
        regimes=[CallerAssayRegime.LCWGS],
        inputs=[CallerInputRole.ALIGNED_BAM],
        rationale="Existing multi-resolution depth and ACE fit evidence.",
    ),
    CallerMode.SPECTRE_DEPTH_ONLY: _entry(
        CallerMode.SPECTRE_DEPTH_ONLY,
        CallerProvider.SPECTRE,
        domains=[CallerAnalyticalDomain.GENOME_WIDE_CNV],
        regimes=[CallerAssayRegime.LCWGS],
        inputs=[CallerInputRole.COVERAGE_PROFILE],
        rationale="Independent depth-only Spectre comparison without caller support.",
    ),
    CallerMode.SPECTRE_SNIFFLES_SUPPORTED: _entry(
        CallerMode.SPECTRE_SNIFFLES_SUPPORTED,
        CallerProvider.SPECTRE,
        domains=[CallerAnalyticalDomain.GENOME_WIDE_CNV],
        regimes=[CallerAssayRegime.LCWGS],
        inputs=[CallerInputRole.COVERAGE_PROFILE],
        parents=[CallerMode.SNIFFLES2_STANDARD],
        rationale="Spectre mode with explicit upstream Sniffles evidence dependency.",
    ),
    CallerMode.SNIFFLES2_STANDARD: _entry(
        CallerMode.SNIFFLES2_STANDARD,
        CallerProvider.SNIFFLES2,
        domains=[CallerAnalyticalDomain.STRUCTURAL_VARIANT],
        regimes=[
            CallerAssayRegime.ADAPTIVE_SAMPLING,
            CallerAssayRegime.TUMOR_NORMAL_LONG_READ,
            CallerAssayRegime.TUMOR_ONLY_LONG_READ,
        ],
        inputs=[CallerInputRole.ALIGNED_BAM],
        rationale="Existing long-read structural-variant adapter in standard mode.",
    ),
    CallerMode.SNIFFLES2_MOSAIC: _entry(
        CallerMode.SNIFFLES2_MOSAIC,
        CallerProvider.SNIFFLES2,
        domains=[CallerAnalyticalDomain.STRUCTURAL_VARIANT],
        regimes=[
            CallerAssayRegime.ADAPTIVE_SAMPLING,
            CallerAssayRegime.TUMOR_ONLY_LONG_READ,
        ],
        inputs=[CallerInputRole.ALIGNED_BAM],
        rationale="Separate mosaic calling mode requiring independent qualification.",
    ),
    CallerMode.CUTESV_STANDARD: _entry(
        CallerMode.CUTESV_STANDARD,
        CallerProvider.CUTESV,
        domains=[CallerAnalyticalDomain.STRUCTURAL_VARIANT],
        regimes=[
            CallerAssayRegime.ADAPTIVE_SAMPLING,
            CallerAssayRegime.TUMOR_NORMAL_LONG_READ,
            CallerAssayRegime.TUMOR_ONLY_LONG_READ,
        ],
        inputs=[CallerInputRole.ALIGNED_BAM],
        rationale="Existing second long-read structural-variant software evidence.",
    ),
    CallerMode.SEVERUS_PAIRED: _entry(
        CallerMode.SEVERUS_PAIRED,
        CallerProvider.SEVERUS,
        domains=[CallerAnalyticalDomain.STRUCTURAL_VARIANT],
        regimes=[CallerAssayRegime.TUMOR_NORMAL_LONG_READ],
        inputs=[
            CallerInputRole.TUMOR_BAM,
            CallerInputRole.MATCHED_NORMAL_BAM,
        ],
        rationale="Conservative paired tumor-normal complex SV mode.",
    ),
    CallerMode.SEVERUS_SINGLE_SAMPLE: _entry(
        CallerMode.SEVERUS_SINGLE_SAMPLE,
        CallerProvider.SEVERUS,
        domains=[CallerAnalyticalDomain.STRUCTURAL_VARIANT],
        regimes=[CallerAssayRegime.TUMOR_ONLY_LONG_READ],
        inputs=[CallerInputRole.TUMOR_BAM],
        rationale="Single-sample Severus mode kept separate from somatic paired claims.",
    ),
    CallerMode.SAVANA_PAIRED: _entry(
        CallerMode.SAVANA_PAIRED,
        CallerProvider.SAVANA,
        domains=[
            CallerAnalyticalDomain.STRUCTURAL_VARIANT,
            CallerAnalyticalDomain.TUMOR_CNA,
        ],
        regimes=[CallerAssayRegime.TUMOR_NORMAL_LONG_READ],
        inputs=[
            CallerInputRole.TUMOR_BAM,
            CallerInputRole.MATCHED_NORMAL_BAM,
            CallerInputRole.SNP_VCF,
        ],
        rationale="Paired tumor-normal SAVANA SV and CNA analysis contract.",
    ),
    CallerMode.SAVANA_TUMOR_ONLY: _entry(
        CallerMode.SAVANA_TUMOR_ONLY,
        CallerProvider.SAVANA,
        domains=[
            CallerAnalyticalDomain.STRUCTURAL_VARIANT,
            CallerAnalyticalDomain.TUMOR_CNA,
        ],
        regimes=[CallerAssayRegime.TUMOR_ONLY_LONG_READ],
        inputs=[
            CallerInputRole.TUMOR_BAM,
            CallerInputRole.SNP_VCF,
        ],
        rationale="Explicit tumor-only SAVANA mode without paired-mode fallback.",
    ),
    CallerMode.WAKHAN_PHASED_CNA: _entry(
        CallerMode.WAKHAN_PHASED_CNA,
        CallerProvider.WAKHAN,
        domains=[CallerAnalyticalDomain.TUMOR_CNA],
        regimes=[
            CallerAssayRegime.TUMOR_NORMAL_LONG_READ,
            CallerAssayRegime.TUMOR_ONLY_LONG_READ,
        ],
        inputs=[
            CallerInputRole.TUMOR_BAM,
            CallerInputRole.PHASED_VARIANTS,
        ],
        rationale="Allele-aware CNA analysis requiring explicit phased evidence.",
    ),
    CallerMode.ICHORCNA_ULP_WGS: _entry(
        CallerMode.ICHORCNA_ULP_WGS,
        CallerProvider.ICHORCNA,
        domains=[CallerAnalyticalDomain.CFDNA_CNV],
        regimes=[CallerAssayRegime.CFDNA_ULP_WGS],
        inputs=[CallerInputRole.COVERAGE_PROFILE],
        rationale="Dedicated ultra-low-pass cfDNA copy-number analysis lane.",
    ),
}


def get_caller_catalog_entry(mode: CallerMode | str) -> CallerCatalogEntry:
    try:
        parsed_mode = mode if isinstance(mode, CallerMode) else CallerMode(mode)
    except ValueError as exc:
        raise ValueError(f"Unknown multi-caller provider mode: {mode!r}") from exc
    return CALLER_CATALOG[parsed_mode]
