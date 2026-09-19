"""Fail-closed research planning across complementary CNV/SV callers.

This module reads no genomic files and executes no tools. Input fingerprints and
runtime locks are declarations, not installation or analytical-validation receipts.
Every returned lane is NOT_RUN, even when its declared inputs are compatible.
"""

from __future__ import annotations

import hashlib
import json
import re
from typing import Annotated, Any, Literal, Self

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .multicaller_catalog import CALLER_CATALOG, CATALOG_VERSION

Sha256 = Annotated[str, Field(pattern=r"^[0-9a-f]{64}$")]
SampleId = Annotated[str, Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")]
CallerId = Literal[
    "qdnaseq_ace",
    "ont_spectre",
    "ichorcna",
    "sniffles2",
    "cutesv",
    "savana",
    "severus",
    "wakhan",
]
Analysis = Literal["cnv", "sv", "allelic_cn", "purity_ploidy", "complex_sv"]
Mode = Literal["single_sample", "tumour_only", "tumour_normal"]
Build = Literal["GRCh37", "GRCh38"]
Platform = Literal["ONT", "PACBIO_HIFI", "ILLUMINA"]
DataBasis = Literal[
    "lcwgs_genome_wide",
    "adaptive_sampling_on_target",
    "adaptive_sampling_off_target",
    "cfdna_ulp_wgs",
]
ArtifactKind = Literal[
    "snv_vcf",
    "phased_snv_vcf",
    "haplotagged_bam",
    "coverage_bins",
    "het_snp_counts",
    "read_count_wig",
    "gc_track",
    "mappability_track",
    "sv_breakpoints",
    "population_snps",
    "sv_pon",
    "cn_pon",
    "roi_bed",
]
_RESOURCE_KINDS = frozenset(
    {"gc_track", "mappability_track", "population_snps", "sv_pon", "cn_pon", "roi_bed"}
)
_CN_ANALYSES = frozenset({"cnv", "allelic_cn", "purity_ploidy"})


class FrozenPlanningContract(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, allow_inf_nan=False)


class PlanningArtifact(FrozenPlanningContract):
    kind: ArtifactKind
    sha256: Sha256
    genome_build: Build
    reference_sha256: Sha256
    sample_id: SampleId | None = None
    source_input_sha256: Sha256 | None = None
    bin_size_bp: int | None = Field(default=None, strict=True, gt=0)
    producer: CallerId | None = None

    @model_validator(mode="after")
    def source_identity_is_explicit(self) -> Self:
        if self.kind in _RESOURCE_KINDS:
            if self.sample_id is not None or self.source_input_sha256 is not None:
                raise ValueError("Shared resources cannot claim a specimen/input identity")
        elif self.sample_id is None or self.source_input_sha256 is None:
            raise ValueError("Sample-derived artifacts require sample and input fingerprints")
        return self


class CallerPlanningContext(FrozenPlanningContract):
    sample_id: SampleId
    genome_build: Build
    reference_sha256: Sha256
    bam_sha256: Sha256
    assessability_mask_sha256: Sha256
    platform: Platform
    data_basis: DataBasis
    coverage_x: float | None = Field(default=None, strict=True, ge=0)
    coverage_definition: str | None = Field(default=None, min_length=3)
    tumor_fraction: float | None = Field(default=None, strict=True, ge=0, le=1)
    tumor_fraction_method: str | None = Field(default=None, min_length=3)
    tumor_fraction_timepoint: str | None = Field(default=None, min_length=3)
    normal_sample_id: SampleId | None = None
    normal_bam_sha256: Sha256 | None = None
    normal_coverage_x: float | None = Field(default=None, strict=True, ge=0)
    normal_coverage_definition: str | None = Field(default=None, min_length=3)
    artifacts: tuple[PlanningArtifact, ...] = ()

    @model_validator(mode="after")
    def coherent_input_identities(self) -> Self:
        for value, definition in (
            (self.coverage_x, self.coverage_definition),
            (self.normal_coverage_x, self.normal_coverage_definition),
        ):
            if (value is None) != (definition is None):
                raise ValueError("Measured coverage and its definition must be supplied together")
        fraction_metadata = (self.tumor_fraction_method, self.tumor_fraction_timepoint)
        if self.tumor_fraction is None and any(item is not None for item in fraction_metadata):
            raise ValueError("Tumour-fraction metadata cannot imply an unmeasured fraction")
        if self.tumor_fraction is not None and any(item is None for item in fraction_metadata):
            raise ValueError("Measured tumour fraction requires method and timepoint")
        if (self.normal_sample_id is None) != (self.normal_bam_sha256 is None):
            raise ValueError("Normal sample identity and BAM fingerprint must be paired")
        if self.normal_sample_id == self.sample_id or self.normal_bam_sha256 == self.bam_sha256:
            raise ValueError("The tumour/sample cannot be its own matched normal")
        if self.normal_sample_id is None and self.normal_coverage_x is not None:
            raise ValueError("Normal coverage cannot imply an undeclared normal specimen")
        kinds = [item.kind for item in self.artifacts]
        if len(kinds) != len(set(kinds)):
            raise ValueError("Input artifact kinds must be unique")
        for item in self.artifacts:
            if (item.genome_build, item.reference_sha256) != (
                self.genome_build,
                self.reference_sha256,
            ):
                raise ValueError("Artifact genome build/reference differs from the sample")
            if item.kind in _RESOURCE_KINDS:
                continue
            if item.sample_id == self.sample_id:
                expected_input = self.bam_sha256
            elif (
                self.normal_sample_id is not None
                and item.sample_id == self.normal_sample_id
                and item.kind in {"snv_vcf", "phased_snv_vcf"}
            ):
                expected_input = self.normal_bam_sha256
            else:
                raise ValueError("Artifact specimen does not match its declared analytical role")
            if item.source_input_sha256 != expected_input:
                raise ValueError("Artifact source input fingerprint differs from the declared BAM")
        return self


class CallerStudyPolicy(FrozenPlanningContract):
    """Prospective engineering bounds; no default is a scientific recommendation."""

    caller_id: CallerId
    genome_build: Build
    platform: Platform
    data_basis: DataBasis
    reference_sha256: Sha256
    assessability_mask_sha256: Sha256
    coverage_definition: str = Field(min_length=3)
    normal_coverage_definition: str | None = Field(default=None, min_length=3)
    allowed_modes: tuple[Mode, ...] = Field(min_length=1)
    minimum_coverage_x: float = Field(strict=True, gt=0)
    minimum_normal_coverage_x: float | None = Field(default=None, strict=True, gt=0)
    bin_size_bp: int | None = Field(default=None, strict=True, gt=0)

    @model_validator(mode="after")
    def unique_modes(self) -> Self:
        if len(self.allowed_modes) != len(set(self.allowed_modes)):
            raise ValueError("Registered analysis modes must be unique")
        return self


class CallerRuntimePin(FrozenPlanningContract):
    version: str = Field(min_length=1)
    runtime_sha256: Sha256
    adapter_sha256: Sha256
    parameters_sha256: Sha256

    @field_validator("version")
    @classmethod
    def pinned_version(cls, value: str) -> str:
        if (
            value.lower() in {"latest", "main", "master", "head", "unknown", "unpinned"}
            or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._+-]*", value) is None
        ):
            raise ValueError("Caller version must be an explicit label, never a moving selector")
        return value


class CallerSelection(FrozenPlanningContract):
    caller_id: CallerId
    analyses: tuple[Analysis, ...] = Field(min_length=1)
    mode: Mode
    enabled: bool = Field(default=True, strict=True)
    policy: CallerStudyPolicy | None = None
    runtime: CallerRuntimePin | None = None
    segmentation: Literal["change_point", "sv_guided"] | None = None

    @model_validator(mode="after")
    def coherent_selection(self) -> Self:
        if len(self.analyses) != len(set(self.analyses)):
            raise ValueError("Requested analyses must be unique")
        if self.segmentation is not None and self.caller_id != "wakhan":
            raise ValueError("This segmentation selector applies only to the Wakhan route")
        return self


class MultiCallerRequest(FrozenPlanningContract):
    schema_version: Literal["0.1.0"] = "0.1.0"
    study_registration_sha256: Sha256
    context: CallerPlanningContext
    selections: tuple[CallerSelection, ...] = ()
    research_only: Literal[True] = True
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def unique_callers(self) -> Self:
        ids = [selection.caller_id for selection in self.selections]
        if len(ids) != len(set(ids)):
            raise ValueError("A caller may have only one lane in this request; use separate plans")
        return self


class PlannedCaller(FrozenPlanningContract):
    caller_id: CallerId
    analyses: tuple[Analysis, ...]
    mode: Mode | None
    disposition: Literal["NOT_SELECTED", "BLOCKED", "ADAPTER_PENDING", "PLANNED"]
    existing_adapter: bool
    input_eligible: bool
    blockers: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    dependencies: tuple[str, ...] = ()
    policy_sha256: Sha256 | None = None
    runtime_pin_sha256: Sha256 | None = None
    execution_status: Literal["NOT_RUN"] = "NOT_RUN"


def _canonical_sha256(payload: dict[str, Any]) -> str:
    rendered = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )
    return hashlib.sha256(rendered.encode("utf-8")).hexdigest()


def _catalog_sha256() -> str:
    return _canonical_sha256(
        {
            "catalog_version": CATALOG_VERSION,
            "callers": [
                {
                    "caller_id": entry.caller_id,
                    "label": entry.label,
                    "analyses": sorted(entry.analyses),
                    "signal_families": sorted(entry.signal_families),
                    "method_family": entry.method_family,
                    "existing_adapter": entry.existing_adapter,
                    "source_repository": entry.source_repository,
                }
                for _, entry in sorted(CALLER_CATALOG.items())
            ],
        }
    )


def _policy_blockers(context: CallerPlanningContext, selection: CallerSelection) -> list[str]:
    blockers: list[str] = []
    policy = selection.policy
    if selection.runtime is None:
        blockers.append("RUNTIME_LOCK_MISSING")
    if context.coverage_x is None:
        blockers.append("COVERAGE_UNKNOWN")
    if policy is None:
        blockers.append("STUDY_POLICY_MISSING")
    else:
        if (
            policy.caller_id,
            policy.genome_build,
            policy.platform,
            policy.data_basis,
            policy.reference_sha256,
            policy.assessability_mask_sha256,
        ) != (
            selection.caller_id,
            context.genome_build,
            context.platform,
            context.data_basis,
            context.reference_sha256,
            context.assessability_mask_sha256,
        ):
            blockers.append("STUDY_POLICY_IDENTITY_MISMATCH")
        if (
            context.coverage_x is not None
            and context.coverage_definition != policy.coverage_definition
        ):
            blockers.append("COVERAGE_DEFINITION_MISMATCH")
        if selection.mode not in policy.allowed_modes:
            blockers.append("MODE_NOT_REGISTERED")
        if context.coverage_x is not None and context.coverage_x < policy.minimum_coverage_x:
            blockers.append("COVERAGE_BELOW_REGISTERED_MINIMUM")
    if selection.mode == "tumour_normal":
        if context.normal_sample_id is None:
            blockers.append("MATCHED_NORMAL_MISSING")
        if context.normal_coverage_x is None:
            blockers.append("NORMAL_COVERAGE_UNKNOWN")
        if policy is not None and context.normal_coverage_x is not None:
            if policy.normal_coverage_definition is None:
                blockers.append("NORMAL_COVERAGE_POLICY_DEFINITION_MISSING")
            elif policy.normal_coverage_definition != context.normal_coverage_definition:
                blockers.append("NORMAL_COVERAGE_DEFINITION_MISMATCH")
        minimum = None if policy is None else policy.minimum_normal_coverage_x
        if minimum is None:
            blockers.append("NORMAL_COVERAGE_POLICY_MISSING")
        elif context.normal_coverage_x is not None and context.normal_coverage_x < minimum:
            blockers.append("NORMAL_COVERAGE_BELOW_REGISTERED_MINIMUM")
    return blockers


def _scope_checks(
    context: CallerPlanningContext,
    selection: CallerSelection,
    blockers: list[str],
    warnings: list[str],
) -> None:
    caller_id = selection.caller_id
    cn_requested = bool(set(selection.analyses) & _CN_ANALYSES)
    if not set(selection.analyses).issubset(CALLER_CATALOG[caller_id].analyses):
        blockers.append("ANALYSIS_NOT_SUPPORTED")
    if caller_id == "ichorcna":
        if context.data_basis != "cfdna_ulp_wgs" or context.platform != "ILLUMINA":
            blockers.append("OUTSIDE_INITIAL_CALLER_DATA_SCOPE")
    else:
        if context.platform not in {"ONT", "PACBIO_HIFI"}:
            blockers.append("OUTSIDE_INITIAL_PLATFORM_SCOPE")
        if context.data_basis == "cfdna_ulp_wgs":
            blockers.append("OUTSIDE_INITIAL_CALLER_DATA_SCOPE")
    if context.data_basis.startswith("adaptive_sampling"):
        if cn_requested:
            blockers.append("AS_GENOMEWIDE_CN_NOT_QUALIFIED")
        else:
            warnings.append("AS_BREAKPOINT_ASSESSABILITY_REQUIRES_REVIEW")
    if selection.mode == "tumour_normal" and caller_id not in {"savana", "severus", "wakhan"}:
        blockers.append("PAIRED_MODE_NOT_IMPLEMENTED_FOR_THIS_ROUTE")
    if selection.mode == "single_sample" and caller_id in {
        "savana", "severus", "wakhan", "ichorcna"
    }:
        blockers.append("EXPLICIT_TUMOUR_MODE_REQUIRED")
    if selection.mode == "tumour_only" and set(selection.analyses) & {"sv", "complex_sv"}:
        warnings.append("UNPAIRED_SV_IS_CANDIDATE_EVIDENCE")
    if selection.mode != "single_sample" and caller_id == "ont_spectre":
        warnings.append("TUMOUR_DIPLOID_BASELINE_ASSUMPTION_REQUIRES_VALIDATION")


def _artifact_checks(
    context: CallerPlanningContext,
    selection: CallerSelection,
    blockers: list[str],
    dependencies: list[str],
) -> None:
    artifacts: dict[str, PlanningArtifact] = {item.kind: item for item in context.artifacts}
    caller_id = selection.caller_id
    if caller_id == "ont_spectre":
        snvs = artifacts.get("snv_vcf")
        if snvs is None or snvs.sample_id != context.sample_id:
            blockers.append("SAMPLE_SNV_REQUIRED")
        coverage = artifacts.get("coverage_bins")
        if coverage is None:
            blockers.append("COVERAGE_BINS_REQUIRED")
        bin_size = None if selection.policy is None else selection.policy.bin_size_bp
        if bin_size is None:
            blockers.append("BIN_SIZE_UNREGISTERED")
        elif coverage is not None and coverage.bin_size_bp != bin_size:
            blockers.append("COVERAGE_BIN_SIZE_MISMATCH")
        dependencies.extend(("coverage_bins", "sample_snvs"))
    if caller_id == "savana" and set(selection.analyses) & _CN_ANALYSES:
        if not {"snv_vcf", "het_snp_counts", "population_snps"}.intersection(artifacts):
            blockers.append("SAVANA_SNP_EVIDENCE_REQUIRED")
        dependencies.extend(("savana_native_breakpoints", "snp_evidence"))
    if caller_id == "wakhan":
        if "phased_snv_vcf" not in artifacts:
            blockers.append("PHASED_SNVS_REQUIRED")
        if "haplotagged_bam" not in artifacts:
            blockers.append("HAPLOTAGGED_BAM_REQUIRED")
        if selection.segmentation is None:
            blockers.append("SEGMENTATION_UNDECLARED")
        elif selection.segmentation == "sv_guided":
            breakpoints = artifacts.get("sv_breakpoints")
            if breakpoints is None:
                blockers.append("SV_BREAKPOINTS_REQUIRED")
            elif breakpoints.producer is None:
                blockers.append("BREAKPOINT_PRODUCER_UNDECLARED")
            elif "sv" not in CALLER_CATALOG[breakpoints.producer].analyses:
                blockers.append("BREAKPOINT_PRODUCER_NOT_SV_CALLER")
            else:
                dependencies.append(f"evidence:{breakpoints.producer}")
        dependencies.extend(("phased_snvs", "haplotagged_bam"))
    if caller_id == "ichorcna":
        for kind in ("read_count_wig", "gc_track", "mappability_track"):
            if kind not in artifacts:
                blockers.append(f"{kind.upper()}_REQUIRED")
            dependencies.append(kind)


def _derive_lanes(request: MultiCallerRequest) -> tuple[PlannedCaller, ...]:
    selections: dict[str, CallerSelection] = {item.caller_id: item for item in request.selections}
    lanes: list[PlannedCaller] = []
    for caller_id, definition in sorted(CALLER_CATALOG.items()):
        selection = selections.get(caller_id)
        values: dict[str, Any] = {
            "caller_id": caller_id,
            "analyses": () if selection is None else selection.analyses,
            "mode": None if selection is None else selection.mode,
            "existing_adapter": definition.existing_adapter,
            "input_eligible": False,
            "disposition": "NOT_SELECTED",
        }
        if selection is not None and selection.enabled:
            blockers = _policy_blockers(request.context, selection)
            warnings = ["SHARED_INPUT_IS_NOT_ORTHOGONAL_TRUTH"]
            dependencies = ["intake_identity", "assay_qc", "assessability_mask"]
            _scope_checks(request.context, selection, blockers, warnings)
            _artifact_checks(request.context, selection, blockers, dependencies)
            values.update(
                {
                    "input_eligible": not blockers,
                    "disposition": (
                        "BLOCKED"
                        if blockers
                        else "PLANNED" if definition.existing_adapter else "ADAPTER_PENDING"
                    ),
                    "blockers": tuple(sorted(set(blockers))),
                    "warnings": tuple(sorted(set(warnings))),
                    "dependencies": tuple(sorted(set(dependencies))),
                    "policy_sha256": (
                        None
                        if selection.policy is None
                        else _canonical_sha256(selection.policy.model_dump(mode="json"))
                    ),
                    "runtime_pin_sha256": (
                        None
                        if selection.runtime is None
                        else _canonical_sha256(selection.runtime.model_dump(mode="json"))
                    ),
                }
            )
        lanes.append(PlannedCaller.model_validate(values))
    return tuple(lanes)


class MultiCallerPlan(FrozenPlanningContract):
    schema_version: Literal["0.1.0"] = "0.1.0"
    request: MultiCallerRequest
    catalog_sha256: Sha256
    lanes: tuple[PlannedCaller, ...]
    plan_sha256: Sha256
    comparison_policy: Literal["NO_VOTING"] = "NO_VOTING"
    execution_enabled: Literal[False] = False
    research_only: Literal[True] = True
    human_review_required: Literal[True] = True
    retain_all_evidence: Literal[True] = True
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def verify_plan(self) -> Self:
        payload = self.model_dump(mode="json", exclude={"plan_sha256"})
        if _canonical_sha256(payload) != self.plan_sha256:
            raise ValueError("Multi-caller plan does not match its content lock")
        if self.catalog_sha256 != _catalog_sha256():
            raise ValueError("Multi-caller plan references a different catalog revision")
        if self.lanes != _derive_lanes(self.request):
            raise ValueError("Multi-caller lanes do not match the declared input/policy decisions")
        return self


def plan_multicaller(request: MultiCallerRequest) -> MultiCallerPlan:
    """Return a sealed metadata-only plan, never a live-caller execution receipt."""
    payload = request.model_dump(mode="json")
    payload["context"]["artifacts"] = sorted(
        payload["context"]["artifacts"],
        key=lambda item: item["kind"],
    )
    payload["selections"] = sorted(payload["selections"], key=lambda item: item["caller_id"])
    for selection in payload["selections"]:
        selection["analyses"] = sorted(selection["analyses"])
        if selection["policy"] is not None:
            selection["policy"]["allowed_modes"] = sorted(selection["policy"]["allowed_modes"])
    canonical_request = MultiCallerRequest.model_validate(payload)
    values: dict[str, Any] = {
        "schema_version": "0.1.0",
        "request": canonical_request.model_dump(mode="json"),
        "catalog_sha256": _catalog_sha256(),
        "lanes": [lane.model_dump(mode="json") for lane in _derive_lanes(canonical_request)],
        "comparison_policy": "NO_VOTING",
        "execution_enabled": False,
        "research_only": True,
        "human_review_required": True,
        "retain_all_evidence": True,
        "sensitive_output": True,
    }
    values["plan_sha256"] = _canonical_sha256(values)
    return MultiCallerPlan.model_validate(values)
