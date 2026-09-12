from __future__ import annotations

import re
from datetime import UTC, date, datetime
from enum import StrEnum
from pathlib import PurePosixPath, PureWindowsPath
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .iscn_syntax import render_supported_cnv_fragment


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)


class InputKind(StrEnum):
    POD5 = "pod5"
    UNALIGNED_BAM = "unaligned_bam"
    ALIGNED_BAM = "aligned_bam"


class AssayMode(StrEnum):
    LOW_COVERAGE_WGS = "lcwgs"
    ADAPTIVE_SAMPLING = "adaptive_sampling"


class GenomeBuild(StrEnum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"


class ReferenceDictionaryContract(StrEnum):
    """Explicit relationship between an aligned BAM and the reference bundle dictionary."""

    EXACT_FULL = "exact_full"
    GRCH38_CANONICAL_25 = "grch38_canonical_25"
    GRCH37_UCSC_HG19_CANONICAL_25 = "grch37_ucsc_hg19_canonical_25"


class CoordinateSystem(StrEnum):
    """Coordinate systems accepted at a resource boundary.

    Analysis code uses only ``zero_based_half_open``. The other values exist so a source
    artifact can state what it actually contains before it is normalized.
    """

    ZERO_BASED_HALF_OPEN = "zero_based_half_open"
    ONE_BASED_INCLUSIVE = "one_based_inclusive"
    ONE_BASED_POSITION = "one_based_position"


class AnalysisModule(StrEnum):
    QC = "qc"
    CNV = "cnv"
    SV = "sv"
    FUSION = "fusion"
    ISCN = "iscn"
    REPORT = "report"
    SMALL_VARIANTS = "small_variants"
    METHYLATION = "methylation"


class Verdict(StrEnum):
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"
    NOT_RUN = "NOT_RUN"


class CheckStatus(StrEnum):
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class ModuleRunStatus(StrEnum):
    COMPLETED = "COMPLETED"
    NOT_RUN = "NOT_RUN"
    FAILED = "FAILED"
    NO_CALL = "NO_CALL"


class BenchmarkKind(StrEnum):
    CNV = "cnv"
    SV = "sv"


class SnifflesMode(StrEnum):
    GERMLINE = "germline"
    MOSAIC = "mosaic"


class SvObservability(StrEnum):
    OBSERVED_ADEQUATELY = "OBSERVED_ADEQUATELY"
    PARTIALLY_OBSERVED = "PARTIALLY_OBSERVED"
    INSUFFICIENT_COVERAGE = "INSUFFICIENT_COVERAGE"
    OUTSIDE_TARGET = "OUTSIDE_TARGET"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class SvValidationStatus(StrEnum):
    DETECTED = "detected"
    TECHNICALLY_SUPPORTED = "technically_supported"
    BIOLOGICALLY_PRIORITIZED = "biologically_prioritized"
    ANALYTICALLY_VALIDATED = "analytically_validated"
    REPORTABLE = "reportable"


class FusionSupportStatus(StrEnum):
    NOT_ASSESSED = "not_assessed"
    CANDIDATE = "fusion_candidate"
    SUPPORTED = "fusion_supported"
    VALIDATED = "fusion_validated"


class EventType(StrEnum):
    CHROMOSOME_GAIN = "chromosome_gain"
    CHROMOSOME_LOSS = "chromosome_loss"
    DELETION = "deletion"
    DUPLICATION = "duplication"
    INVERSION = "inversion"
    TRANSLOCATION = "translocation"
    INSERTION = "insertion"
    FUSION = "fusion"


class ReviewStatus(StrEnum):
    DRAFT = "DRAFT"
    REVIEW_REQUIRED = "REVIEW_REQUIRED"
    REVIEWED = "REVIEWED"


class ISCNProposalStatus(StrEnum):
    """How much nomenclature context the automated proposal actually contains."""

    LEGACY_UNSPECIFIED = "LEGACY_UNSPECIFIED"
    NOT_REQUESTED = "NOT_REQUESTED"
    NOT_ASSESSED = "NOT_ASSESSED"
    NO_RENDERABLE_CANDIDATE = "NO_RENDERABLE_CANDIDATE"
    PARTIAL_EVENT_LEVEL = "PARTIAL_EVENT_LEVEL"


class ISCNSelectionPolicy(StrEnum):
    """Versioned boundary selecting events for the deliberately limited renderer."""

    REPORTABLE_ONLY_V1 = "REPORTABLE_ONLY_V1"
    TECHNICAL_CANDIDATES_V1 = "TECHNICAL_CANDIDATES_V1"


class ISCNDispositionOutcome(StrEnum):
    """Auditable disposition of every event considered by the proposal renderer."""

    RENDERED = "RENDERED"
    OMITTED_UNSUPPORTED = "OMITTED_UNSUPPORTED"
    EXCLUDED_BY_POLICY = "EXCLUDED_BY_POLICY"
    NOT_REQUESTED = "NOT_REQUESTED"
    ASSESSMENT_BLOCKED = "ASSESSMENT_BLOCKED"


class InputSpec(StrictModel):
    kind: InputKind
    path: str = Field(min_length=1)
    index_path: str | None = None
    sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")

    @model_validator(mode="after")
    def aligned_bam_requires_index(self) -> InputSpec:
        if self.kind == InputKind.ALIGNED_BAM and not self.index_path:
            raise ValueError("aligned_bam requires index_path")
        return self


class TargetBedRole(StrEnum):
    """What a target BED actually is. The two are not interchangeable.

    The BED a sequencer selects on usually carries flanks so that a read starting outside a
    gene is still enriched. Coverage computed over that buffered design answers "was the
    enrichment working"; it does not answer "was the analysis region observed", because the
    flanks dilute the per-target mean. Reporting one as the other overstates or understates
    adequacy depending on the flank size, so the manifest has to say which one it points at.
    """

    ANALYSIS_ROI_UNBUFFERED = "analysis_roi_unbuffered"
    SELECTION_PANEL_BUFFERED = "selection_panel_buffered"


class AssaySpec(StrictModel):
    mode: AssayMode
    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    target_bed: str | None = None
    target_bed_version: str | None = None
    #: Defaults to the unbuffered analysis ROI because that is the stricter reading: a run
    #: that does not declare a role is treated as claiming the narrower meaning, and a
    #: buffered design has to say so explicitly.
    target_bed_role: TargetBedRole = TargetBedRole.ANALYSIS_ROI_UNBUFFERED

    @model_validator(mode="after")
    def adaptive_sampling_requires_bed(self) -> AssaySpec:
        if self.mode == AssayMode.ADAPTIVE_SAMPLING and (
            not self.target_bed or not self.target_bed_version
        ):
            raise ValueError("adaptive_sampling requires target_bed and target_bed_version")
        return self


class AnalysisIntent(StrEnum):
    """Whether this analysis is looking for acquired or inherited variation.

    Load-bearing for knowledge-base annotation: ClinVar classifies germline variation under
    ACMG rules, and an AML workup asks a somatic question. Pairing the two without saying so
    presents an inherited-disease classification as though it answered a question about a
    tumour.
    """

    SOMATIC = "somatic"
    GERMLINE = "germline"
    BOTH = "both"


class AnalysisSpec(StrictModel):
    profile: str = Field(min_length=1)
    modules: list[AnalysisModule]
    parameters: dict[str, Any] = Field(default_factory=dict)
    #: What kind of variation this analysis is asking about. Optional and **without a
    #: default**: guessing would silently decide how every knowledge-base assertion is read.
    #: Left unset, scope alignment is reported as unknown rather than assumed.
    intent: AnalysisIntent | None = None


class PrivacySpec(StrictModel):
    pseudonymized: bool = True
    contains_direct_identifiers: bool = False
    cloud_upload_approved: bool = False

    @model_validator(mode="after")
    def block_direct_identifiers(self) -> PrivacySpec:
        if self.contains_direct_identifiers:
            raise ValueError("Direct identifiers are prohibited in pipeline manifests")
        return self


class SampleManifest(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    run_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    input: InputSpec
    assay: AssaySpec
    analysis: AnalysisSpec
    privacy: PrivacySpec = Field(default_factory=PrivacySpec)


class ReferenceContig(StrictModel):
    name: str = Field(min_length=1)
    length: int = Field(gt=0)


class ReferenceLock(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    reference_id: str = Field(min_length=1)
    genome_build: GenomeBuild
    contigs: list[ReferenceContig] = Field(min_length=1)
    allow_extra_contigs: bool = False
    source_fai_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")

    @model_validator(mode="after")
    def contig_names_are_unique(self) -> ReferenceLock:
        names = [item.name for item in self.contigs]
        if len(names) != len(set(names)):
            raise ValueError("reference lock contains duplicate contig names")
        return self


_RESOURCE_ID_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$"
_RESOURCE_ROLE_PATTERN = r"^[a-z][a-z0-9_]{0,63}$"
_SHA256_PATTERN = r"^[0-9a-f]{64}$"


def _validate_relative_artifact_path(value: str) -> str:
    """Keep every manifest path inside its bundle or run envelope.

    Manifests use POSIX separators on every host. Rejecting Windows separators avoids a
    path meaning one thing under WSL and another thing in the Desktop process.
    """

    if "\\" in value:
        raise ValueError("artifact paths must use portable POSIX separators")
    posix = PurePosixPath(value)
    windows = PureWindowsPath(value)
    if value in {"", "."} or posix.is_absolute() or windows.is_absolute() or windows.drive:
        raise ValueError("artifact path must be relative")
    if ".." in posix.parts or ".." in windows.parts:
        raise ValueError("artifact path cannot traverse outside its bundle")
    return posix.as_posix()


def _valid_resource_reference(value: str) -> bool:
    """Accept a local ID or an explicit immutable ``bundle_id:resource_id`` reference."""

    parts = value.split(":")
    return len(parts) in {1, 2} and all(
        re.fullmatch(_RESOURCE_ID_PATTERN, part) is not None for part in parts
    )


class ResourceFile(StrictModel):
    """One source or derived file declared by a bundle manifest.

    Source files are checksum-pinned. A generated file may omit its checksum in a catalog
    manifest, but it then remains unresolved until an installer has materialized and pinned
    it in the activated manifest.
    """

    resource_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    role: str = Field(pattern=_RESOURCE_ROLE_PATTERN)
    path: str = Field(min_length=1)
    sha256: str | None = Field(default=None, pattern=_SHA256_PATTERN)
    size_bytes: int | None = Field(default=None, ge=0)
    source_url: str | None = Field(default=None, min_length=1)
    release: str | None = Field(default=None, min_length=1)
    source_date: date | None = None
    coordinate_system: CoordinateSystem | None = None
    generated: bool = False
    derived_from: list[str] = Field(default_factory=list)
    required: bool = True
    description: str = ""

    @field_validator("path")
    @classmethod
    def path_stays_inside_bundle(cls, value: str) -> str:
        return _validate_relative_artifact_path(value)

    @field_validator("derived_from")
    @classmethod
    def derivation_ids_are_unique(cls, value: list[str]) -> list[str]:
        if len(value) != len(set(value)):
            raise ValueError("derived_from contains duplicate resource IDs")
        for resource_id in value:
            if not _valid_resource_reference(resource_id):
                raise ValueError(f"invalid derived resource ID: {resource_id!r}")
        return value

    @model_validator(mode="after")
    def source_and_generated_files_are_explicit(self) -> ResourceFile:
        if not self.generated and self.sha256 is None:
            raise ValueError("source resource files require a SHA256 checksum")
        if self.generated and not self.derived_from:
            raise ValueError("generated resource files require derived_from provenance")
        if self.resource_id in self.derived_from:
            raise ValueError("a resource cannot be derived from itself")
        return self


class ResourceBundle(StrictModel):
    """Fields shared by all manifest-activated bundle types."""

    schema_version: Literal["1.0.0"] = "1.0.0"
    bundle_type: Literal["reference", "panel", "knowledge"]
    bundle_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    version: str = Field(min_length=1)
    genome_build: GenomeBuild | None = None
    resources: list[ResourceFile] = Field(min_length=1)
    description: str = ""

    @model_validator(mode="after")
    def resource_ids_and_roles_are_unique(self) -> ResourceBundle:
        resource_ids = [item.resource_id for item in self.resources]
        if len(resource_ids) != len(set(resource_ids)):
            raise ValueError("bundle contains duplicate resource IDs")
        roles = [item.role for item in self.resources]
        if len(roles) != len(set(roles)):
            raise ValueError("bundle contains duplicate resource roles")
        known = set(resource_ids)
        for resource in self.resources:
            local_references = {item for item in resource.derived_from if ":" not in item}
            missing = local_references.difference(known)
            if missing:
                raise ValueError(
                    f"resource {resource.resource_id!r} has unknown derived_from IDs: "
                    f"{', '.join(sorted(missing))}"
                )
        return self

    def resource(self, resource_id: str) -> ResourceFile:
        for resource in self.resources:
            if resource.resource_id == resource_id:
                return resource
        raise ValueError(f"bundle {self.bundle_id!r} does not declare resource {resource_id!r}")


class ReferenceBundle(ResourceBundle):
    bundle_type: Literal["reference"] = "reference"
    genome_build: GenomeBuild
    reference_lock_resource_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    fasta_resource_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    fai_resource_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    annotation_cache_resource_id: str | None = Field(default=None, pattern=_RESOURCE_ID_PATTERN)

    @model_validator(mode="after")
    def required_reference_resources_exist(self) -> ReferenceBundle:
        required_ids = [
            self.reference_lock_resource_id,
            self.fasta_resource_id,
            self.fai_resource_id,
        ]
        if self.annotation_cache_resource_id is not None:
            required_ids.append(self.annotation_cache_resource_id)
        if len(required_ids) != len(set(required_ids)):
            raise ValueError("reference bundle required resource IDs must be distinct")
        expected_roles = {
            self.reference_lock_resource_id: "reference_lock",
            self.fasta_resource_id: "genome_fasta",
            self.fai_resource_id: "fasta_index",
        }
        if self.annotation_cache_resource_id is not None:
            expected_roles[self.annotation_cache_resource_id] = "annotation_cache"
        for resource_id, role in expected_roles.items():
            resource = self.resource(resource_id)
            if resource.role != role:
                raise ValueError(
                    f"reference resource {resource_id!r} must have role {role!r}, "
                    f"not {resource.role!r}"
                )
        return self


class PanelCoordinateMappingLock(StrictModel):
    """Immutable provenance for an offline, build-time panel coordinate mapping.

    The lock describes an already materialized mapping result. It is never an instruction to
    fetch a chain or execute mapping during resource activation or analysis.
    """

    schema_version: Literal["1.0.0"] = "1.0.0"
    mapping_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    status: Literal["controlled_build_time_unvalidated"]
    source_genome_build: GenomeBuild
    target_genome_build: GenomeBuild
    source_panel_bundle_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    source_panel_resource_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    source_selection_sha256: str = Field(pattern=_SHA256_PATTERN)
    source_interval_count: int = Field(ge=1)
    mapping_method: Literal["ucsc_liftover_chain"]
    mapping_tool: Literal["UCSC liftOver"]
    mapping_tool_url: str = Field(min_length=1)
    mapping_tool_sha256: str = Field(pattern=_SHA256_PATTERN)
    mapping_tool_size_bytes: int = Field(ge=1)
    mapping_tool_build_id: str = Field(min_length=1)
    mapping_tool_last_modified: datetime
    chain_name: str = Field(min_length=1)
    chain_url: str = Field(min_length=1)
    chain_md5: str = Field(pattern=r"^[0-9a-f]{32}$")
    chain_sha256: str = Field(pattern=_SHA256_PATTERN)
    min_match: float = Field(gt=0, le=1)
    min_blocks: float = Field(gt=0, le=1)
    multiple: bool
    mapped_output_sha256: str = Field(pattern=_SHA256_PATTERN)
    unmapped_output_sha256: str = Field(pattern=_SHA256_PATTERN)
    mapped_interval_count: int = Field(ge=0)
    same_chromosome_mapped_count: int = Field(ge=0)
    span_identical_count: int = Field(ge=0)
    maximum_span_delta_bases: int = Field(ge=0)
    maximum_span_delta_fraction: float = Field(ge=0)
    unmapped_target_labels: list[str] = Field(default_factory=list)
    split_target_labels: list[str] = Field(default_factory=list)
    roundtrip_chain_name: str = Field(min_length=1)
    roundtrip_chain_url: str = Field(min_length=1)
    roundtrip_chain_md5: str = Field(pattern=r"^[0-9a-f]{32}$")
    roundtrip_chain_sha256: str = Field(pattern=_SHA256_PATTERN)
    roundtrip_min_match: float = Field(gt=0, le=1)
    roundtrip_min_blocks: float = Field(gt=0, le=1)
    roundtrip_multiple: bool
    roundtrip_mapped_output_sha256: str = Field(pattern=_SHA256_PATTERN)
    roundtrip_mapped_output_size_bytes: int = Field(ge=0)
    roundtrip_unmapped_output_sha256: str = Field(pattern=_SHA256_PATTERN)
    roundtrip_unmapped_output_size_bytes: int = Field(ge=0)
    roundtrip_mapped_interval_count: int = Field(ge=0)
    reciprocal_exact_interval_count: int = Field(ge=0)
    roundtrip_unmapped_target_labels: list[str] = Field(default_factory=list)
    roundtrip_split_target_labels: list[str] = Field(default_factory=list)
    roundtrip_review_required_target_labels: list[str] = Field(default_factory=list)
    final_selection_sha256: str = Field(pattern=_SHA256_PATTERN)
    runtime_mapping_allowed: Literal[False] = False
    note: str = ""

    @model_validator(mode="after")
    def mapping_counts_and_builds_are_consistent(self) -> PanelCoordinateMappingLock:
        if self.source_genome_build == self.target_genome_build:
            raise ValueError("panel coordinate mapping source and target builds must differ")
        if len(self.unmapped_target_labels) != len(set(self.unmapped_target_labels)):
            raise ValueError("panel mapping unmapped target labels must be unique")
        if len(self.split_target_labels) != len(set(self.split_target_labels)):
            raise ValueError("panel mapping split target labels must be unique")
        if not set(self.split_target_labels).issubset(self.unmapped_target_labels):
            raise ValueError("panel mapping split targets must also be declared unmapped")
        if len(self.roundtrip_review_required_target_labels) != len(
            set(self.roundtrip_review_required_target_labels)
        ):
            raise ValueError("panel mapping roundtrip review target labels must be unique")
        if len(self.roundtrip_unmapped_target_labels) != len(
            set(self.roundtrip_unmapped_target_labels)
        ):
            raise ValueError("panel mapping roundtrip unmapped target labels must be unique")
        if len(self.roundtrip_split_target_labels) != len(set(self.roundtrip_split_target_labels)):
            raise ValueError("panel mapping roundtrip split target labels must be unique")
        if not set(self.roundtrip_split_target_labels).issubset(
            self.roundtrip_unmapped_target_labels
        ):
            raise ValueError("panel mapping roundtrip split targets must be declared unmapped")
        if not set(self.roundtrip_unmapped_target_labels).issubset(
            self.roundtrip_review_required_target_labels
        ):
            raise ValueError("panel mapping roundtrip-unmapped targets must require review")
        if set(self.roundtrip_review_required_target_labels).intersection(
            self.unmapped_target_labels
        ):
            raise ValueError(
                "forward-unmapped targets cannot also require mapped-interval roundtrip review"
            )
        if self.mapped_interval_count + len(self.unmapped_target_labels) != (
            self.source_interval_count
        ):
            raise ValueError("panel mapping mapped and unmapped counts do not cover the source")
        if self.same_chromosome_mapped_count > self.mapped_interval_count:
            raise ValueError("same-chromosome mapping count exceeds mapped interval count")
        if self.span_identical_count > self.mapped_interval_count:
            raise ValueError("span-identical mapping count exceeds mapped interval count")
        if (
            self.roundtrip_mapped_interval_count + len(self.roundtrip_unmapped_target_labels)
            != self.mapped_interval_count
        ):
            raise ValueError(
                "panel mapping roundtrip mapped and unmapped counts do not cover forward mappings"
            )
        if self.reciprocal_exact_interval_count > self.roundtrip_mapped_interval_count:
            raise ValueError("reciprocal-exact count exceeds roundtrip mapped interval count")
        return self


class PanelBundle(ResourceBundle):
    bundle_type: Literal["panel"] = "panel"
    genome_build: GenomeBuild
    assay_mode: AssayMode
    reference_dictionary_contracts: list[ReferenceDictionaryContract] = Field(
        default_factory=lambda: [ReferenceDictionaryContract.EXACT_FULL],
        min_length=1,
    )
    target_gene_map_resource_id: str | None = Field(default=None, pattern=_RESOURCE_ID_PATTERN)
    coordinate_mapping_lock_resource_id: str | None = Field(
        default=None, pattern=_RESOURCE_ID_PATTERN
    )
    coordinate_mapping_source_resource_id: str | None = Field(
        default=None, pattern=_RESOURCE_ID_PATTERN
    )
    coordinate_mapping_output_resource_id: str | None = Field(
        default=None, pattern=_RESOURCE_ID_PATTERN
    )
    coordinate_mapping_unmapped_resource_id: str | None = Field(
        default=None, pattern=_RESOURCE_ID_PATTERN
    )
    coordinate_mapping_roundtrip_output_resource_id: str | None = Field(
        default=None, pattern=_RESOURCE_ID_PATTERN
    )
    coordinate_mapping_roundtrip_unmapped_resource_id: str | None = Field(
        default=None, pattern=_RESOURCE_ID_PATTERN
    )
    coordinate_mapping_origin_bundle_id: str | None = Field(
        default=None, pattern=_RESOURCE_ID_PATTERN
    )
    coordinate_mapping_origin_resource_id: str | None = Field(
        default=None, pattern=_RESOURCE_ID_PATTERN
    )
    selection_panel_resource_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    analysis_roi_resource_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    transcript_cache_resource_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    unresolved_targets: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def required_panel_resources_exist(self) -> PanelBundle:
        if self.assay_mode != AssayMode.ADAPTIVE_SAMPLING:
            raise ValueError("panel bundles require assay_mode='adaptive_sampling'")
        if len(self.reference_dictionary_contracts) != len(
            set(self.reference_dictionary_contracts)
        ):
            raise ValueError("panel reference_dictionary_contracts must be unique")
        for contract in self.reference_dictionary_contracts:
            if (
                contract == ReferenceDictionaryContract.GRCH38_CANONICAL_25
                and self.genome_build != GenomeBuild.GRCH38
            ):
                raise ValueError("grch38_canonical_25 panels require genome_build='GRCh38'")
            if (
                contract == ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25
                and self.genome_build != GenomeBuild.GRCH37
            ):
                raise ValueError(
                    "grch37_ucsc_hg19_canonical_25 panels require genome_build='GRCh37'"
                )
        expected_roles = {
            self.selection_panel_resource_id: TargetBedRole.SELECTION_PANEL_BUFFERED.value,
            self.analysis_roi_resource_id: TargetBedRole.ANALYSIS_ROI_UNBUFFERED.value,
            self.transcript_cache_resource_id: "transcript_cache",
        }
        for resource_id, role in expected_roles.items():
            resource = self.resource(resource_id)
            if resource.role != role:
                raise ValueError(
                    f"panel resource {resource_id!r} must have role {role!r}, not {resource.role!r}"
                )
            if (
                role
                in {
                    TargetBedRole.SELECTION_PANEL_BUFFERED.value,
                    TargetBedRole.ANALYSIS_ROI_UNBUFFERED.value,
                }
                and resource.coordinate_system != CoordinateSystem.ZERO_BASED_HALF_OPEN
            ):
                raise ValueError(
                    f"panel resource {resource_id!r} with role {role!r} must declare "
                    "coordinate_system='zero_based_half_open'"
                )
        if self.target_gene_map_resource_id is not None:
            gene_map = self.resource(self.target_gene_map_resource_id)
            if gene_map.role != "target_gene_map":
                raise ValueError("panel target gene map must have role 'target_gene_map'")
            if gene_map.coordinate_system is not None:
                raise ValueError("panel target gene maps must be coordinate-free")
            for resource_id in (
                self.analysis_roi_resource_id,
                self.transcript_cache_resource_id,
            ):
                if self.target_gene_map_resource_id not in self.resource(resource_id).derived_from:
                    raise ValueError(
                        "panel ROI and transcript cache must derive from the target gene map"
                    )
        mapping_ids = (
            self.coordinate_mapping_lock_resource_id,
            self.coordinate_mapping_source_resource_id,
            self.coordinate_mapping_output_resource_id,
            self.coordinate_mapping_unmapped_resource_id,
            self.coordinate_mapping_roundtrip_output_resource_id,
            self.coordinate_mapping_roundtrip_unmapped_resource_id,
        )
        mapping_origin = (
            self.coordinate_mapping_origin_bundle_id,
            self.coordinate_mapping_origin_resource_id,
        )
        if any(item is not None for item in mapping_ids):
            if any(item is None for item in mapping_ids):
                raise ValueError("panel coordinate mapping resource IDs must be declared together")
            if any(item is None for item in mapping_origin):
                raise ValueError("panel coordinate mapping origin IDs must be declared together")
            if self.target_gene_map_resource_id is None:
                raise ValueError("coordinate-mapped panels require a target gene map")
            lock_id, source_id, output_id, unmapped_id, roundtrip_id, roundtrip_unmapped_id = (
                mapping_ids
            )
            assert lock_id is not None
            assert source_id is not None
            assert output_id is not None
            assert unmapped_id is not None
            assert roundtrip_id is not None
            assert roundtrip_unmapped_id is not None
            expected_mapping_roles = {
                lock_id: "coordinate_mapping_lock",
                source_id: "coordinate_mapping_source",
                output_id: "coordinate_mapping_output",
                unmapped_id: "coordinate_mapping_unmapped",
                roundtrip_id: "coordinate_mapping_roundtrip_output",
                roundtrip_unmapped_id: "coordinate_mapping_roundtrip_unmapped",
            }
            for resource_id, role in expected_mapping_roles.items():
                if self.resource(resource_id).role != role:
                    raise ValueError(
                        f"panel coordinate mapping resource {resource_id!r} must have role {role!r}"
                    )
            output = self.resource(output_id)
            if not output.generated or source_id not in output.derived_from:
                raise ValueError("panel mapping output must derive from its source selection")
            if lock_id not in output.derived_from:
                raise ValueError("panel mapping output must derive from its mapping lock")
            unmapped = self.resource(unmapped_id)
            if (
                not unmapped.generated
                or source_id not in unmapped.derived_from
                or lock_id not in unmapped.derived_from
            ):
                raise ValueError(
                    "panel unmapped output must derive from its source selection and mapping lock"
                )
            for resource_id in (roundtrip_id, roundtrip_unmapped_id):
                resource = self.resource(resource_id)
                if (
                    not resource.generated
                    or output_id not in resource.derived_from
                    or lock_id not in resource.derived_from
                ):
                    raise ValueError(
                        "panel roundtrip outputs must derive from the forward mapping output "
                        "and mapping lock"
                    )
            for resource_id in (
                source_id,
                output_id,
                unmapped_id,
                roundtrip_id,
                roundtrip_unmapped_id,
            ):
                if (
                    self.resource(resource_id).coordinate_system
                    != CoordinateSystem.ZERO_BASED_HALF_OPEN
                ):
                    raise ValueError(
                        "panel coordinate mapping BED resources require "
                        "coordinate_system='zero_based_half_open'"
                    )
            if self.resource(lock_id).coordinate_system is not None:
                raise ValueError("panel coordinate mapping locks must be coordinate-free")
            selection = self.resource(self.selection_panel_resource_id)
            if not selection.generated or output_id not in selection.derived_from:
                raise ValueError("panel selection must derive from its mapped coordinate output")
            if self.target_gene_map_resource_id is not None and (
                self.target_gene_map_resource_id not in selection.derived_from
            ):
                raise ValueError("panel selection must derive from its target gene map")
        elif any(item is not None for item in mapping_origin):
            raise ValueError("panel coordinate mapping origin IDs require mapping resources")
        if len(self.unresolved_targets) != len(set(self.unresolved_targets)):
            raise ValueError("panel unresolved_targets must be unique")
        return self


class KnowledgeBundle(ResourceBundle):
    bundle_type: Literal["knowledge"] = "knowledge"
    coordinate_bearing: bool = True

    @model_validator(mode="after")
    def coordinate_resources_require_a_build(self) -> KnowledgeBundle:
        if self.coordinate_bearing and self.genome_build is None:
            raise ValueError("coordinate-bearing knowledge bundles require a genome build")
        if self.genome_build is None and any(
            resource.coordinate_system is not None for resource in self.resources
        ):
            raise ValueError("knowledge resources with coordinates require a genome build")
        return self


class AnalysisProfile(StrictModel):
    schema_version: Literal["1.0.0"] = "1.0.0"
    profile_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    version: str = Field(min_length=1)
    genome_build: GenomeBuild
    assay_mode: AssayMode
    reference_bundle: str = Field(pattern=_RESOURCE_ID_PATTERN)
    reference_dictionary_contract: ReferenceDictionaryContract = (
        ReferenceDictionaryContract.EXACT_FULL
    )
    knowledge_bundle: str = Field(pattern=_RESOURCE_ID_PATTERN)
    panel_bundle: str | None = Field(default=None, pattern=_RESOURCE_ID_PATTERN)
    adaptive_sampling: Literal["enabled", "disabled"]
    description: str = ""

    @model_validator(mode="after")
    def assay_and_panel_are_consistent(self) -> AnalysisProfile:
        if (
            self.reference_dictionary_contract == ReferenceDictionaryContract.GRCH38_CANONICAL_25
            and self.genome_build != GenomeBuild.GRCH38
        ):
            raise ValueError("grch38_canonical_25 is valid only for GRCh38 profiles")
        if (
            self.reference_dictionary_contract
            == ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25
            and self.genome_build != GenomeBuild.GRCH37
        ):
            raise ValueError("grch37_ucsc_hg19_canonical_25 is valid only for GRCh37 profiles")
        if self.assay_mode == AssayMode.ADAPTIVE_SAMPLING:
            if self.adaptive_sampling != "enabled" or self.panel_bundle is None:
                raise ValueError(
                    "adaptive-sampling profiles require adaptive_sampling=enabled and a panel"
                )
        elif self.adaptive_sampling != "disabled" or self.panel_bundle is not None:
            raise ValueError("lcWGS profiles require adaptive_sampling=disabled and no panel")
        return self


class PanelResolutionSummary(StrictModel):
    """Run-bound summary of how an Adaptive Sampling panel reached its active build.

    This is provenance for panel design resolution, not sample evidence.  In particular,
    coordinate mapping does not establish analytical validation, and unresolved labels do not
    describe negative findings in the analysed sample.
    """

    schema_version: Literal["1.0.0"] = "1.0.0"
    mapping_status: Literal[
        "native_build_not_required",
        "controlled_build_time_unvalidated",
    ]
    mapping_id: str | None = Field(default=None, pattern=_RESOURCE_ID_PATTERN)
    mapping_method: Literal["ucsc_liftover_chain"] | None = None
    mapping_tool: Literal["UCSC liftOver"] | None = None
    source_genome_build: GenomeBuild
    target_genome_build: GenomeBuild
    source_panel_bundle_id: str | None = Field(default=None, pattern=_RESOURCE_ID_PATTERN)
    source_panel_resource_id: str | None = Field(default=None, pattern=_RESOURCE_ID_PATTERN)
    selection_interval_count: int = Field(ge=1)
    analysis_roi_interval_count: int = Field(ge=0)
    source_interval_count: int | None = Field(default=None, ge=1)
    mapped_interval_count: int | None = Field(default=None, ge=0)
    unmapped_target_labels: list[str] = Field(default_factory=list)
    roundtrip_mapped_interval_count: int | None = Field(default=None, ge=0)
    reciprocal_exact_interval_count: int | None = Field(default=None, ge=0)
    roundtrip_review_required_target_labels: list[str] = Field(default_factory=list)
    unresolved_target_labels: list[str] = Field(default_factory=list)
    selection_panel_sha256: str = Field(pattern=_SHA256_PATTERN)
    analysis_roi_sha256: str = Field(pattern=_SHA256_PATTERN)
    coordinate_mapping_lock_sha256: str | None = Field(
        default=None,
        pattern=_SHA256_PATTERN,
    )
    roundtrip_mapping_output_sha256: str | None = Field(
        default=None,
        pattern=_SHA256_PATTERN,
    )
    roundtrip_mapping_unmapped_sha256: str | None = Field(
        default=None,
        pattern=_SHA256_PATTERN,
    )
    runtime_mapping_allowed: Literal[False] = False

    @model_validator(mode="after")
    def mapping_state_is_consistent(self) -> PanelResolutionSummary:
        label_groups = (
            self.unmapped_target_labels,
            self.roundtrip_review_required_target_labels,
            self.unresolved_target_labels,
        )
        if any(len(labels) != len(set(labels)) for labels in label_groups):
            raise ValueError("panel resolution target labels must be unique within each state")

        mapping_values = (
            self.mapping_id,
            self.mapping_method,
            self.mapping_tool,
            self.source_panel_bundle_id,
            self.source_panel_resource_id,
            self.source_interval_count,
            self.mapped_interval_count,
            self.roundtrip_mapped_interval_count,
            self.reciprocal_exact_interval_count,
            self.coordinate_mapping_lock_sha256,
            self.roundtrip_mapping_output_sha256,
            self.roundtrip_mapping_unmapped_sha256,
        )
        if self.mapping_status == "native_build_not_required":
            if self.source_genome_build != self.target_genome_build:
                raise ValueError(
                    "native panel resolution requires identical source and target builds"
                )
            if any(value is not None for value in mapping_values):
                raise ValueError("native panel resolution cannot declare coordinate-mapping fields")
            if self.unmapped_target_labels or self.roundtrip_review_required_target_labels:
                raise ValueError("native panel resolution cannot declare coordinate-mapping labels")
            if self.analysis_roi_interval_count > self.selection_interval_count:
                raise ValueError("panel analysis ROI count cannot exceed selection count")
            return self

        if self.source_genome_build == self.target_genome_build:
            raise ValueError("coordinate-mapped panel resolution requires different builds")
        if any(value is None for value in mapping_values):
            raise ValueError(
                "coordinate-mapped panel resolution requires complete mapping provenance"
            )
        assert self.source_interval_count is not None
        assert self.mapped_interval_count is not None
        assert self.roundtrip_mapped_interval_count is not None
        assert self.reciprocal_exact_interval_count is not None
        if self.mapped_interval_count + len(self.unmapped_target_labels) != (
            self.source_interval_count
        ):
            raise ValueError("panel resolution mapped and unmapped counts do not cover the source")
        if self.roundtrip_mapped_interval_count > self.mapped_interval_count:
            raise ValueError("panel resolution roundtrip count exceeds mapped interval count")
        if self.reciprocal_exact_interval_count > self.roundtrip_mapped_interval_count:
            raise ValueError("panel resolution reciprocal-exact count exceeds roundtrip count")
        if self.selection_interval_count != self.mapped_interval_count:
            raise ValueError("mapped panel selection count must equal the forward-mapped count")
        if self.analysis_roi_interval_count > self.selection_interval_count:
            raise ValueError("panel analysis ROI count cannot exceed selection count")
        return self


class ResolvedResourceContext(StrictModel):
    """Absolute, checksum-pinned resources selected for one analysis profile."""

    schema_version: Literal["1.0.0"] = "1.0.0"
    profile_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    profile_version: str = Field(min_length=1)
    genome_build: GenomeBuild
    reference_dictionary_contract: ReferenceDictionaryContract = (
        ReferenceDictionaryContract.EXACT_FULL
    )
    reference_bundle_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    reference_bundle_version: str = Field(min_length=1)
    panel_bundle_id: str | None = Field(default=None, pattern=_RESOURCE_ID_PATTERN)
    panel_bundle_version: str | None = None
    panel_resolution: PanelResolutionSummary | None = None
    knowledge_bundle_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    knowledge_bundle_version: str = Field(min_length=1)
    resource_root: str = Field(min_length=1)
    resource_paths: dict[str, str] = Field(default_factory=dict)
    resource_checksums: dict[str, str] = Field(default_factory=dict)
    resource_releases: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def paths_checksums_and_optional_panel_are_consistent(self) -> ResolvedResourceContext:
        if (
            self.reference_dictionary_contract == ReferenceDictionaryContract.GRCH38_CANONICAL_25
            and self.genome_build != GenomeBuild.GRCH38
        ):
            raise ValueError("grch38_canonical_25 is valid only for GRCh38 contexts")
        if (
            self.reference_dictionary_contract
            == ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25
            and self.genome_build != GenomeBuild.GRCH37
        ):
            raise ValueError("grch37_ucsc_hg19_canonical_25 is valid only for GRCh37 contexts")
        if (self.panel_bundle_id is None) != (self.panel_bundle_version is None):
            raise ValueError("panel bundle ID and version must either both be set or both be null")
        if self.panel_resolution is not None:
            if self.panel_bundle_id is None:
                raise ValueError("panel resolution requires a resolved panel bundle")
            if self.panel_resolution.target_genome_build != self.genome_build:
                raise ValueError("panel resolution target build must match the resolved context")
            summary_checksums = {
                "panel.selection_panel_buffered": self.panel_resolution.selection_panel_sha256,
                "panel.analysis_roi_unbuffered": self.panel_resolution.analysis_roi_sha256,
            }
            if self.panel_resolution.coordinate_mapping_lock_sha256 is not None:
                summary_checksums["panel.coordinate_mapping_lock"] = (
                    self.panel_resolution.coordinate_mapping_lock_sha256
                )
            if self.panel_resolution.roundtrip_mapping_output_sha256 is not None:
                summary_checksums["panel.coordinate_mapping_roundtrip_output"] = (
                    self.panel_resolution.roundtrip_mapping_output_sha256
                )
            if self.panel_resolution.roundtrip_mapping_unmapped_sha256 is not None:
                summary_checksums["panel.coordinate_mapping_roundtrip_unmapped"] = (
                    self.panel_resolution.roundtrip_mapping_unmapped_sha256
                )
            for key, expected in summary_checksums.items():
                if self.resource_checksums.get(key) != expected:
                    raise ValueError(
                        f"panel resolution checksum for {key!r} does not match resolved resources"
                    )
        if set(self.resource_paths) != set(self.resource_checksums):
            raise ValueError("every resolved resource path must have exactly one checksum")
        root: PurePosixPath | PureWindowsPath
        if self.resource_root.startswith("/"):
            root = PurePosixPath(self.resource_root)
        else:
            root = PureWindowsPath(self.resource_root)
        if not root.is_absolute():
            raise ValueError("resolved resource_root must be absolute")
        if ".." in root.parts:
            raise ValueError("resolved resource_root must not contain parent traversal")
        for key, value in self.resource_paths.items():
            path: PurePosixPath | PureWindowsPath
            if isinstance(root, PurePosixPath):
                path = PurePosixPath(value)
            else:
                path = PureWindowsPath(value)
            if not path.is_absolute():
                raise ValueError(f"resolved resource path {key!r} must be absolute")
            if ".." in path.parts:
                raise ValueError(f"resolved resource path {key!r} contains parent traversal")
            try:
                path.relative_to(root)
            except ValueError as exc:
                raise ValueError(f"resolved resource path {key!r} escapes resource_root") from exc
        return self


class LegacyResourceContext(StrictModel):
    """Marker used when a PipelineResult 0.1 artifact had no bundle provenance."""

    status: Literal["legacy_unspecified"] = "legacy_unspecified"


class SidecarArtifact(StrictModel):
    """A checksum-pinned large table kept outside ``PipelineResult``."""

    artifact_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    relative_path: str = Field(min_length=1)
    schema_version: str = Field(min_length=1)
    sha256: str = Field(pattern=_SHA256_PATTERN)
    row_count: int = Field(ge=0)
    size_bytes: int | None = Field(default=None, ge=0)
    media_type: str = Field(default="text/tab-separated-values", min_length=1)
    columns: list[str] = Field(default_factory=list)

    @field_validator("relative_path")
    @classmethod
    def path_stays_inside_run_envelope(cls, value: str) -> str:
        return _validate_relative_artifact_path(value)


class QCMetrics(StrictModel):
    verdict: Verdict
    metrics: dict[str, float | int | str | None]
    warnings: list[str] = Field(default_factory=list)
    failed_gates: list[str] = Field(default_factory=list)


class QCPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    status: Literal["technical_defaults_only", "validated"]
    hard_failures: list[str] = Field(default_factory=list)
    numeric_gates: dict[str, float | int | None] = Field(default_factory=dict)
    note: str


class SnifflesPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    expected_version: str = Field(default="2.8.0", pattern=r"^\d+\.\d+\.\d+$")
    mode: SnifflesMode = SnifflesMode.GERMLINE
    min_support: int = Field(default=5, ge=1)
    min_sv_length: int = Field(default=50, ge=1)
    mapq: int = Field(default=20, ge=0, le=60)
    pass_only: Literal[True] = True
    minimum_quality: float | None = Field(default=None, ge=0)
    require_precise: bool = False
    allowed_sv_types: list[EventType] = Field(
        default_factory=lambda: [
            EventType.DELETION,
            EventType.DUPLICATION,
            EventType.INVERSION,
            EventType.INSERTION,
            EventType.TRANSLOCATION,
        ],
        min_length=1,
    )
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def supported_event_types_are_unique(self) -> SnifflesPolicy:
        supported = {
            EventType.DELETION,
            EventType.DUPLICATION,
            EventType.INVERSION,
            EventType.INSERTION,
            EventType.TRANSLOCATION,
        }
        if any(item not in supported for item in self.allowed_sv_types):
            raise ValueError("Sniffles policy contains an unsupported event type")
        if len(self.allowed_sv_types) != len(set(self.allowed_sv_types)):
            raise ValueError("Sniffles policy contains duplicate event types")
        return self


class CuteSvPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    expected_version: str = Field(default="2.1.3", pattern=r"^\d+\.\d+\.\d+$")
    min_support: int = Field(default=5, ge=1)
    min_sv_length: int = Field(default=50, ge=1)
    max_cluster_bias_ins: int = Field(default=100, ge=0)
    diff_ratio_merging_ins: float = Field(default=0.3, ge=0, le=1)
    max_cluster_bias_del: int = Field(default=100, ge=0)
    diff_ratio_merging_del: float = Field(default=0.3, ge=0, le=1)
    minimum_quality: float | None = Field(default=None, ge=0)
    pass_only: Literal[True] = True
    allowed_sv_types: list[EventType] = Field(
        default_factory=lambda: [
            EventType.DELETION,
            EventType.DUPLICATION,
            EventType.INVERSION,
            EventType.INSERTION,
            EventType.TRANSLOCATION,
        ],
        min_length=1,
    )
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def supported_event_types_are_unique(self) -> CuteSvPolicy:
        supported = {
            EventType.DELETION,
            EventType.DUPLICATION,
            EventType.INVERSION,
            EventType.INSERTION,
            EventType.TRANSLOCATION,
        }
        if any(item not in supported for item in self.allowed_sv_types):
            raise ValueError("cuteSV policy contains an unsupported event type")
        if len(self.allowed_sv_types) != len(set(self.allowed_sv_types)):
            raise ValueError("cuteSV policy contains duplicate event types")
        return self


class SvConsensusPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    maximum_breakpoint_distance_bp: int = Field(default=500, ge=0)
    minimum_reciprocal_overlap: float = Field(default=0.5, ge=0, le=1)
    maximum_length_ratio_difference: float = Field(default=0.35, ge=0, le=1)
    require_orientation_when_available: bool = True
    merge_within_caller: bool = True
    note: str = Field(min_length=1)


class SvEvidencePolicy(StrictModel):
    """Transparent, versioned weights for technical SV review prioritization."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    high_score: int = 8
    moderate_score: int = 5
    support_high: int = Field(default=20, ge=1)
    support_moderate: int = Field(default=10, ge=1)
    support_minimum: int = Field(default=5, ge=1)
    vaf_high: float = Field(default=0.10, ge=0, le=1)
    vaf_minimum: float = Field(default=0.05, ge=0, le=1)
    caller_consensus_weight: int = 4
    single_caller_weight: int = 1
    support_high_weight: int = 4
    support_moderate_weight: int = 2
    support_minimum_weight: int = 1
    precise_breakpoint_weight: int = 1
    vaf_high_weight: int = 2
    vaf_minimum_weight: int = 1
    non_pass_filter_weight: int = -3
    adequate_observability_weight: int = 1
    inadequate_observability_weight: int = -1
    context_flag_weight: int = -1
    maximum_context_penalty: int = Field(default=3, ge=0)
    known_aml_pattern_weight: int = 2
    balanced_sv_weight: int = 1
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def thresholds_are_ordered(self) -> SvEvidencePolicy:
        if self.high_score <= self.moderate_score:
            raise ValueError("high_score must be greater than moderate_score")
        if not self.support_minimum <= self.support_moderate <= self.support_high:
            raise ValueError("support thresholds must be ordered minimum <= moderate <= high")
        if self.vaf_minimum > self.vaf_high:
            raise ValueError("VAF thresholds must be ordered minimum <= high")
        if self.context_flag_weight > 0:
            raise ValueError("context_flag_weight may not reward artifact context")
        return self


class Locus(StrictModel):
    chromosome: str = Field(pattern=r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    cytoband_start: str | None = None
    cytoband_end: str | None = None
    gene: str | None = None

    @model_validator(mode="after")
    def end_after_start(self) -> Locus:
        if self.end <= self.start:
            raise ValueError("locus end must be greater than start")
        return self


class Evidence(StrictModel):
    caller: str
    caller_version: str
    support_reads: int | None = Field(default=None, ge=0)
    local_coverage: float | None = Field(default=None, ge=0)
    variant_allele_fraction: float | None = Field(default=None, ge=0, le=1)
    quality: float | None = Field(default=None, ge=0)
    filters: list[str] = Field(default_factory=list)
    supporting_read_strands: str | None = Field(
        default=None,
        pattern=r"^[+-]{1,2}$",
    )
    coverage_context: list[Annotated[float, Field(ge=0)]] = Field(default_factory=list)
    mean_alignment_nm: float | None = Field(default=None, ge=0)
    position_standard_deviation: float | None = Field(default=None, ge=0)
    length_standard_deviation: float | None = Field(default=None, ge=0)
    precise: bool | None = None


class KnowledgeResourceLock(StrictModel):
    """The exact knowledge-base release an annotation came from.

    ClinVar republishes weekly. "ClinVar says Pathogenic" without saying *which* ClinVar is
    not reproducible: the same BAM can yield different reports a month apart with nothing
    recording why. Locked by checksum, exactly as the reference genome and the cytoband
    table are.
    """

    schema_version: Literal["0.1.0"] = "0.1.0"
    source_id: str = Field(min_length=1)
    #: The publisher's own release identifier, e.g. a ClinVar file date.
    release: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    genome_build: GenomeBuild
    #: The classification system this source's assertions are written in.
    assertion_vocabulary: str = Field(min_length=1)
    #: Records read, records usable, and why the rest were not.
    records_loaded: int = Field(default=0, ge=0)
    load_summary: str = ""


class IntervalResourceLock(StrictModel):
    """Version and checksum lock for a build-specific interval annotation resource."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    resource_id: str = Field(min_length=1)
    resource_type: Literal[
        "genes",
        "cytobands",
        "repeatmasker",
        "tandem_repeat",
        "segmental_duplication",
        "blacklist",
        "mappability",
        "centromere",
        "telomere",
    ]
    source: str = Field(min_length=1)
    release: str = Field(min_length=1)
    genome_build: GenomeBuild
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    source_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    coordinate_system: Literal["zero_based_half_open"] = "zero_based_half_open"
    columns: Literal["chrom_start_end_label"] = "chrom_start_end_label"
    note: str = ""


class AmlKnowledgeLock(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    resource_id: str = Field(min_length=1)
    release: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    source_ids: list[str] = Field(min_length=1)
    note: str = Field(min_length=1)


class PathologyAssociation(StrictModel):
    """One source-attributed disease association for a rearrangement review pattern.

    The association explains why a gene pair is surfaced to a reviewer.  It is deliberately
    not a diagnosis, sample-level assertion, or reportability rule.
    """

    disease_id: str = Field(pattern=r"^DOID:\d+$")
    name: str = Field(min_length=1)
    source_id: str = Field(min_length=1)
    source_record_id: str = Field(min_length=1)
    source_url: str = Field(min_length=1)
    evidence_item_count: int | None = Field(default=None, ge=0)
    assertion_count: int | None = Field(default=None, ge=0)


class AmlRearrangementRecord(StrictModel):
    record_id: str = Field(min_length=1)
    pattern_type: Literal["exact_pair", "gene_any_partner"]
    genes: list[str] = Field(min_length=1, max_length=2)
    display_name: str = Field(min_length=1)
    relevance: Literal[
        "aml_defining_pattern", "aml_relevant_pattern", "hematology_relevant_pattern"
    ]
    source_ids: list[str] = Field(min_length=1)
    pathologies: list[PathologyAssociation] = Field(default_factory=list)
    caveat: str = Field(min_length=1)

    @model_validator(mode="after")
    def pattern_has_expected_gene_count(self) -> AmlRearrangementRecord:
        expected = 2 if self.pattern_type == "exact_pair" else 1
        if len(self.genes) != expected:
            raise ValueError(f"{self.pattern_type} requires exactly {expected} gene(s)")
        if len({gene.upper() for gene in self.genes}) != len(self.genes):
            raise ValueError("AML rearrangement genes must be unique")
        pathology_ids = [item.disease_id for item in self.pathologies]
        if len(pathology_ids) != len(set(pathology_ids)):
            raise ValueError("rearrangement pathologies must have unique disease IDs")
        cited_sources = {item.source_id for item in self.pathologies}
        if not cited_sources.issubset(set(self.source_ids)):
            raise ValueError("pathology association cites a source absent from the record")
        return self


class EventAnnotation(StrictModel):
    """A knowledge-base record attached to a finding — evidence, never a verdict.

    Carries the source's assertion verbatim together with the vocabulary it belongs to, so
    ClinVar's germline *Pathogenic* cannot be read as a claim about a somatic finding. No
    field here influences ``reportable`` or ``confidence``, and none may be added that does:
    that decision needs somatic criteria this repository does not have.
    """

    source_id: str = Field(min_length=1)
    source_release: str = Field(min_length=1)
    source_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    record_id: str = Field(min_length=1)
    record_type: str = Field(min_length=1)
    #: The source's classification, in the source's own words.
    assertion: str
    #: The rule set ``assertion`` belongs to, e.g. ``acmg_germline``.
    assertion_vocabulary: str = Field(min_length=1)
    record_origin: Literal["germline", "somatic", "unknown"]
    scope_alignment: Literal["aligned", "mismatched", "unknown"]
    scope_note: str
    match_type: Literal["exact", "record_within_finding", "finding_within_record", "overlap"]
    reciprocal_overlap: float = Field(ge=0, le=1)
    review_status: str = ""
    #: NCBI's star rating; ``None`` when the review status is not in the known vocabulary.
    review_stars: int | None = Field(default=None, ge=0, le=4)
    genes: list[str] = Field(default_factory=list)
    conditions: list[str] = Field(default_factory=list)
    caveats: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def every_annotation_states_its_limits(self) -> EventAnnotation:
        """An annotation without its caveats is the failure mode this class exists to stop.

        The caveats carry the sentence saying this is a classification of a database record
        rather than a finding about the sample. Allowing an empty list would let exactly the
        reading the design forbids through the one gap left open.
        """
        if not self.caveats:
            raise ValueError(
                "an annotation must carry its caveats; without them a database "
                "classification reads as a finding about the sample"
            )
        return self


class BreakpointTranscriptAnnotation(StrictModel):
    """Transcript-level context for one independently annotated SV breakpoint."""

    gene_id: str = Field(min_length=1)
    gene_name: str = Field(min_length=1)
    transcript_id: str = Field(min_length=1)
    transcript_name: str | None = None
    strand: Literal["+", "-"]
    preferred: bool = False
    rank_tier: int = Field(ge=1)
    region: Literal["exon", "intron", "transcript"]
    exon_number: int | None = Field(default=None, ge=1)
    intron_number: int | None = Field(default=None, ge=1)
    cds_phase: int | None = Field(default=None, ge=0, le=2)


class BreakpointContextAnnotation(StrictModel):
    resource_type: str = Field(min_length=1)
    label: str = Field(min_length=1)
    value: float | None = None


class BreakpointAnnotation(StrictModel):
    """One 0-based breakpoint with cytoband, transcript and technical context."""

    label: Literal["primary", "secondary"]
    chromosome: str = Field(pattern=r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)$")
    position: int = Field(ge=0)
    cytoband: str | None = None
    transcripts: list[BreakpointTranscriptAnnotation] = Field(default_factory=list)
    contexts: list[BreakpointContextAnnotation] = Field(default_factory=list)


class FusionPartnerAnnotation(StrictModel):
    gene: str | None = None
    preferred_transcript: str | None = None
    region: Literal["exon", "intron", "transcript", "unknown"] = "unknown"
    exon_number: int | None = Field(default=None, ge=1)
    intron_number: int | None = Field(default=None, ge=1)
    strand: Literal["+", "-"] | None = None


class FusionAnnotation(StrictModel):
    gene_a: FusionPartnerAnnotation
    gene_b: FusionPartnerAnnotation
    orientation: Literal["++", "+-", "-+", "--"] | None = None
    frame_status: Literal["in_frame", "out_of_frame", "unknown"] = "unknown"


class GenomicEvent(StrictModel):
    event_id: str
    event_type: EventType
    primary: Locus
    secondary: Locus | None = None
    length_bp: int | None = Field(default=None, ge=1)
    copy_number: float | None = Field(default=None, ge=0)
    whole_chromosome_span_confirmed: bool = Field(
        default=False,
        description=(
            "True only when the source CNV span covers the whole chromosome under the versioned "
            "span basis: an exact locked-contig span, or the recorded assessable span; required "
            "before +chr/-chr ISCN rendering"
        ),
    )
    assessable_span_start: int | None = Field(
        default=None,
        ge=0,
        description=(
            "Zero-based start of the assessable region a whole-chromosome confirmation was "
            "measured against when the versioned CNV policy does not require an exact contig span"
        ),
    )
    assessable_span_end: int | None = Field(
        default=None,
        ge=1,
        description=(
            "Half-open end of that assessable region; a confirmed segment must cover it completely"
        ),
    )
    genes: list[str] = Field(default_factory=list)
    evidence: list[Evidence] = Field(default_factory=list)
    confidence: Literal["high", "moderate", "low", "unclassified"] = "unclassified"
    reportable: bool = False
    notes: list[str] = Field(default_factory=list)
    source_event_ids: list[str] = Field(default_factory=list)
    breakpoint_distance_bp: int | None = Field(default=None, ge=0)
    technical_flags: list[str] = Field(default_factory=list)
    observability: SvObservability = SvObservability.NOT_APPLICABLE
    breakpoint_mean_depths: list[float | None] = Field(default_factory=list)
    observability_target_role: TargetBedRole | None = None
    aml_relevance: str | None = None
    known_rearrangement: str | None = None
    known_pathologies: list[PathologyAssociation] = Field(default_factory=list)
    fusion_status: FusionSupportStatus = FusionSupportStatus.NOT_ASSESSED
    validation_status: SvValidationStatus = SvValidationStatus.DETECTED
    breakpoint_annotations: list[BreakpointAnnotation] = Field(default_factory=list)
    fusion_evidence: FusionAnnotation | None = None
    #: Knowledge-base records matching this event. Evidence for a reviewer; deliberately
    #: without influence on ``confidence`` or ``reportable``.
    annotations: list[EventAnnotation] = Field(default_factory=list)

    @model_validator(mode="after")
    def paired_events_require_secondary_locus(self) -> GenomicEvent:
        if self.event_type in {EventType.TRANSLOCATION, EventType.FUSION} and not self.secondary:
            raise ValueError(f"{self.event_type.value} requires secondary locus")
        labels = [item.label for item in self.breakpoint_annotations]
        if len(labels) != len(set(labels)):
            raise ValueError("breakpoint annotations must contain at most one item per label")
        pathology_ids = [item.disease_id for item in self.known_pathologies]
        if len(pathology_ids) != len(set(pathology_ids)):
            raise ValueError("event pathologies must have unique disease IDs")
        if (self.assessable_span_start is None) != (self.assessable_span_end is None):
            raise ValueError("an assessable span requires both a start and an end")
        if (
            self.assessable_span_start is not None
            and self.assessable_span_end is not None
            and self.assessable_span_end <= self.assessable_span_start
        ):
            raise ValueError("an assessable span must end after it starts")
        if self.whole_chromosome_span_confirmed:
            if self.event_type not in {
                EventType.CHROMOSOME_GAIN,
                EventType.CHROMOSOME_LOSS,
            }:
                raise ValueError(
                    "whole chromosome span confirmation is valid only for chromosome gain/loss"
                )
            # Either the segment reaches the locked contig start, or it covers every assessable
            # base the caller could measure on that chromosome. Filtered telomeric and blacklisted
            # bins are never silently treated as unchanged sequence.
            covers_assessable = (
                self.assessable_span_start is not None
                and self.assessable_span_end is not None
                and self.primary.start <= self.assessable_span_start
                and self.primary.end >= self.assessable_span_end
            )
            if self.primary.start != 0 and not covers_assessable:
                raise ValueError(
                    "a confirmed whole chromosome span must begin at coordinate zero or cover "
                    "the recorded assessable span"
                )
        return self


class BenchmarkThresholds(StrictModel):
    minimum_reciprocal_overlap: float = Field(default=0.5, ge=0, le=1)
    maximum_breakpoint_distance_bp: int = Field(default=500, ge=0)
    copy_number_tolerance: float | None = Field(default=None, ge=0)


class BenchmarkCase(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    case_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    kind: BenchmarkKind
    genome_build: GenomeBuild
    truth_events: list[GenomicEvent]
    query_events: list[GenomicEvent]
    thresholds: BenchmarkThresholds = Field(default_factory=BenchmarkThresholds)
    strata: dict[str, float | int | str | None] = Field(default_factory=dict)
    research_only: Literal[True] = True


class BenchmarkMatch(StrictModel):
    truth_event_id: str
    query_event_id: str
    score: float = Field(ge=0, le=1)
    reciprocal_overlap: float | None = Field(default=None, ge=0, le=1)
    maximum_breakpoint_distance_bp: int | None = Field(default=None, ge=0)


class BenchmarkMetrics(StrictModel):
    true_positive: int = Field(ge=0)
    false_positive: int = Field(ge=0)
    false_negative: int = Field(ge=0)
    precision: float | None = Field(default=None, ge=0, le=1)
    recall: float | None = Field(default=None, ge=0, le=1)
    f1: float | None = Field(default=None, ge=0, le=1)


class BenchmarkReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    case_id: str
    kind: BenchmarkKind
    genome_build: GenomeBuild
    thresholds: BenchmarkThresholds
    strata: dict[str, float | int | str | None]
    metrics: BenchmarkMetrics
    matches: list[BenchmarkMatch]
    unmatched_truth_event_ids: list[str]
    unmatched_query_event_ids: list[str]
    warnings: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))


class ISCNEventFragment(StrictModel):
    """One traceable, event-level nomenclature fragment in a proposal."""

    event_id: str = Field(min_length=1)
    event_type: EventType
    fragment: str = Field(min_length=1)
    confidence: Literal["high", "moderate", "low", "unclassified"]
    reportable: bool
    selection_reason: str = Field(min_length=1)


class ISCNEventDisposition(StrictModel):
    """Why one result event was rendered, omitted, or excluded."""

    event_id: str = Field(min_length=1)
    outcome: ISCNDispositionOutcome
    reason_code: str = Field(pattern=r"^[A-Z][A-Z0-9_]{2,63}$")
    detail: str = Field(min_length=1)


class ISCNResourceProvenance(StrictModel):
    """Build and checksums used for coordinate-to-cytoband projection."""

    genome_build: GenomeBuild
    reference_dictionary_contract: ReferenceDictionaryContract
    reference_bundle_id: str = Field(pattern=_RESOURCE_ID_PATTERN)
    reference_bundle_version: str = Field(min_length=1)
    reference_lock_sha256: str = Field(pattern=_SHA256_PATTERN)
    annotation_cache_resource_key: Literal["reference.annotation_cache"] = (
        "reference.annotation_cache"
    )
    annotation_cache_sha256: str = Field(pattern=_SHA256_PATTERN)
    cytoband_resource_key: Literal["reference.cytobands"] = "reference.cytobands"
    cytoband_release: str = Field(min_length=1)
    cytoband_sha256: str = Field(pattern=_SHA256_PATTERN)


class ISCNAssessmentBlocker(StrictModel):
    """A fail-closed reason why no proposal assessment was permitted."""

    reason_code: str = Field(pattern=r"^[A-Z][A-Z0-9_]{2,63}$")
    detail: str = Field(min_length=1)


class ISCNProposal(StrictModel):
    notation: str | None = None
    standard_edition: Literal["ISCN 2024"] = "ISCN 2024"
    conformance_profile: Literal[
        "legacy-unspecified",
        "subset-v0.1-unvalidated",
        "event-fragments-v0.2-unvalidated",
    ] = "legacy-unspecified"
    review_status: ReviewStatus = ReviewStatus.REVIEW_REQUIRED
    proposal_status: ISCNProposalStatus = ISCNProposalStatus.LEGACY_UNSPECIFIED
    genome_build: GenomeBuild | None = None
    selection_policy: ISCNSelectionPolicy | None = None
    resource_provenance: ISCNResourceProvenance | None = None
    baseline_assessed: Literal[False] = False
    requires_expert_review: Literal[True] = True
    clinical_release_allowed: Literal[False] = False
    policy_parameters: dict[str, str | int | float | bool] = Field(default_factory=dict)
    assessment_blockers: list[ISCNAssessmentBlocker] = Field(default_factory=list)
    event_fragments: list[ISCNEventFragment] = Field(default_factory=list)
    event_dispositions: list[ISCNEventDisposition] = Field(default_factory=list)
    considered_event_ids: list[str] = Field(default_factory=list)
    cnv_source_event_ids: list[str] = Field(default_factory=list)
    source_event_ids: list[str] = Field(default_factory=list)
    technical_assumptions: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def proposal_state_is_traceable(self) -> ISCNProposal:
        collections = {
            "considered_event_ids": self.considered_event_ids,
            "cnv_source_event_ids": self.cnv_source_event_ids,
            "source_event_ids": self.source_event_ids,
            "event fragment IDs": [item.event_id for item in self.event_fragments],
            "event disposition IDs": [item.event_id for item in self.event_dispositions],
        }
        for label, values in collections.items():
            if len(values) != len(set(values)):
                raise ValueError(f"{label} must be unique")
        if self.proposal_status == ISCNProposalStatus.LEGACY_UNSPECIFIED:
            return self
        if self.conformance_profile != "event-fragments-v0.2-unvalidated":
            raise ValueError("structured ISCN proposals require the event-fragments-v0.2 profile")
        if self.genome_build is None or self.selection_policy is None:
            raise ValueError("new ISCN proposals require genome build and selection policy")
        if (
            self.proposal_status
            in {
                ISCNProposalStatus.NO_RENDERABLE_CANDIDATE,
                ISCNProposalStatus.PARTIAL_EVENT_LEVEL,
            }
            and self.resource_provenance is None
        ):
            raise ValueError("assessed ISCN proposals require checksum-pinned resource provenance")
        if (
            self.resource_provenance is not None
            and self.resource_provenance.genome_build != self.genome_build
        ):
            raise ValueError("ISCN resource provenance uses a different genome build")
        if self.review_status != ReviewStatus.REVIEW_REQUIRED:
            raise ValueError("automated ISCN proposals always require expert review")
        if self.proposal_status == ISCNProposalStatus.NOT_ASSESSED:
            if not self.assessment_blockers:
                raise ValueError("NOT_ASSESSED ISCN proposals require an assessment blocker")
        elif self.assessment_blockers:
            raise ValueError("assessment blockers are valid only for NOT_ASSESSED proposals")
        if set(self.considered_event_ids) != {item.event_id for item in self.event_dispositions}:
            raise ValueError("every considered ISCN event requires exactly one disposition")
        disposition_outcomes = {item.outcome for item in self.event_dispositions}
        if self.proposal_status == ISCNProposalStatus.NOT_REQUESTED:
            if disposition_outcomes.difference({ISCNDispositionOutcome.NOT_REQUESTED}):
                raise ValueError("NOT_REQUESTED proposals require NOT_REQUESTED dispositions")
        elif self.proposal_status == ISCNProposalStatus.NOT_ASSESSED:
            if disposition_outcomes.difference({ISCNDispositionOutcome.ASSESSMENT_BLOCKED}):
                raise ValueError("NOT_ASSESSED proposals require ASSESSMENT_BLOCKED dispositions")
        elif disposition_outcomes.intersection(
            {
                ISCNDispositionOutcome.NOT_REQUESTED,
                ISCNDispositionOutcome.ASSESSMENT_BLOCKED,
            }
        ):
            raise ValueError("assessed proposals cannot carry unassessed event dispositions")
        if not set(self.cnv_source_event_ids).issubset(self.considered_event_ids):
            raise ValueError("ISCN CNV source IDs must refer to considered events")
        rendered_ids = {
            item.event_id
            for item in self.event_dispositions
            if item.outcome == ISCNDispositionOutcome.RENDERED
        }
        if self.source_event_ids != [item.event_id for item in self.event_fragments]:
            raise ValueError("ISCN source event IDs must match event fragment order")
        if rendered_ids != set(self.source_event_ids):
            raise ValueError("rendered ISCN dispositions must match source event IDs")
        if not rendered_ids.issubset(self.cnv_source_event_ids):
            raise ValueError("every rendered ISCN fragment must originate in the typed CNV report")
        if self.proposal_status == ISCNProposalStatus.PARTIAL_EVENT_LEVEL:
            if not self.notation or not self.event_fragments:
                raise ValueError("a partial ISCN proposal requires notation and event fragments")
            expected = ",".join(item.fragment for item in self.event_fragments)
            if self.notation != expected:
                raise ValueError("ISCN notation must equal the ordered event fragments")
            fragments = [item.fragment for item in self.event_fragments]
            if len(fragments) != len(set(fragments)):
                raise ValueError("ISCN event fragments must be unique")
            if self.selection_policy == ISCNSelectionPolicy.REPORTABLE_ONLY_V1 and any(
                not item.reportable for item in self.event_fragments
            ):
                raise ValueError("reportable-only ISCN proposals cannot render unreportable events")
        elif self.notation is not None or self.event_fragments:
            raise ValueError("a non-partial proposal cannot carry notation or event fragments")
        return self


class ToolRecord(StrictModel):
    name: str
    version: str
    parameters: dict[str, Any] = Field(default_factory=dict)
    container_digest: str | None = None


class FileFingerprint(StrictModel):
    size_bytes: int = Field(ge=0)
    sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")


class ValidationCheck(StrictModel):
    name: str = Field(min_length=1)
    status: CheckStatus
    message: str = Field(min_length=1)
    details: dict[str, int | float | str | bool | None] = Field(default_factory=dict)


class BamHeaderSummary(StrictModel):
    sort_order: str | None = None
    sequence_count: int = Field(ge=0)
    total_reference_bases: int = Field(ge=0)
    contig_signature_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    read_group_count: int = Field(ge=0)
    sample_tag_count: int = Field(ge=0)
    program_count: int = Field(ge=0)


class AlignedBamIntakeReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    reference_id: str
    genome_build: GenomeBuild
    input_fingerprint: FileFingerprint | None = None
    index_fingerprint: FileFingerprint | None = None
    header: BamHeaderSummary | None = None
    checks: list[ValidationCheck]
    verdict: Verdict
    tool: ToolRecord | None = None
    limitations: list[str] = Field(default_factory=list)
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))


class CraminoQCReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    qc: QCMetrics
    tool: ToolRecord
    limitations: list[str] = Field(default_factory=list)
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))


class SnifflesCallReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: SnifflesPolicy
    events: list[GenomicEvent]
    raw_record_count: int = Field(ge=0)
    accepted_record_count: int = Field(ge=0)
    rejected_record_count: int = Field(ge=0)
    rejection_counts: dict[str, int] = Field(default_factory=dict)
    tool: ToolRecord
    vcf_fingerprint: FileFingerprint
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def counts_and_status_are_consistent(self) -> SnifflesCallReport:
        if self.raw_record_count != self.accepted_record_count + self.rejected_record_count:
            raise ValueError("Sniffles record counts are inconsistent")
        if self.accepted_record_count != len(self.events):
            raise ValueError("Sniffles accepted count must equal normalized event count")
        if sum(self.rejection_counts.values()) != self.rejected_record_count:
            raise ValueError("Sniffles rejection reason counts are inconsistent")
        expected_status = ModuleRunStatus.COMPLETED if self.events else ModuleRunStatus.NO_CALL
        if self.status != expected_status:
            raise ValueError("Sniffles status is inconsistent with normalized events")
        return self


class CuteSvCallReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: CuteSvPolicy
    events: list[GenomicEvent]
    raw_record_count: int = Field(ge=0)
    accepted_record_count: int = Field(ge=0)
    rejected_record_count: int = Field(ge=0)
    rejection_counts: dict[str, int] = Field(default_factory=dict)
    tool: ToolRecord
    vcf_fingerprint: FileFingerprint
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def counts_and_status_are_consistent(self) -> CuteSvCallReport:
        if self.raw_record_count != self.accepted_record_count + self.rejected_record_count:
            raise ValueError("cuteSV record counts are inconsistent")
        if self.accepted_record_count != len(self.events):
            raise ValueError("cuteSV accepted count must equal normalized event count")
        if sum(self.rejection_counts.values()) != self.rejected_record_count:
            raise ValueError("cuteSV rejection reason counts are inconsistent")
        expected_status = ModuleRunStatus.COMPLETED if self.events else ModuleRunStatus.NO_CALL
        if self.status != expected_status:
            raise ValueError("cuteSV status is inconsistent with normalized events")
        return self


class SvConsensusReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: SvConsensusPolicy
    events: list[GenomicEvent]
    input_event_count: int = Field(ge=0)
    consolidated_event_count: int = Field(ge=0)
    caller_names: list[str]
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def counts_are_consistent(self) -> SvConsensusReport:
        if self.consolidated_event_count != len(self.events):
            raise ValueError("SV consensus event count is inconsistent")
        if self.consolidated_event_count > self.input_event_count:
            raise ValueError("SV consensus cannot create more events than it received")
        expected_status = ModuleRunStatus.COMPLETED if self.events else ModuleRunStatus.NO_CALL
        if self.status != expected_status:
            raise ValueError("SV consensus status is inconsistent with events")
        return self


class LocalSmokeReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    verdict: Verdict
    intake: AlignedBamIntakeReport
    qc: CraminoQCReport
    sniffles: SnifflesCallReport
    checks: list[ValidationCheck]
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))


class ModuleOutcome(StrictModel):
    module: AnalysisModule
    status: ModuleRunStatus
    reason: str
    tools: list[ToolRecord] = Field(default_factory=list)


class Provenance(StrictModel):
    pipeline_version: str
    git_commit: str
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    tools: list[ToolRecord] = Field(default_factory=list)
    reference_checksums: dict[str, str] = Field(default_factory=dict)


class PipelineResult(StrictModel):
    schema_version: Literal["0.1.0", "0.2.0", "0.3.0"] = "0.3.0"
    manifest: SampleManifest
    qc: QCMetrics
    events: list[GenomicEvent]
    iscn: ISCNProposal
    provenance: Provenance
    reference_context: ResolvedResourceContext | LegacyResourceContext = Field(
        default_factory=LegacyResourceContext
    )
    sidecars: list[SidecarArtifact] = Field(default_factory=list)
    modules: list[ModuleOutcome] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    release_status: ReviewStatus = ReviewStatus.REVIEW_REQUIRED

    @model_validator(mode="after")
    def failed_qc_cannot_be_reviewed(self) -> PipelineResult:
        if self.qc.verdict == Verdict.FAIL and self.release_status == ReviewStatus.REVIEWED:
            raise ValueError("A QC-failed result cannot be marked REVIEWED")
        if self.schema_version in {"0.1.0", "0.2.0"}:
            if self.iscn.proposal_status != ISCNProposalStatus.LEGACY_UNSPECIFIED:
                raise ValueError("legacy PipelineResult schemas require legacy ISCN semantics")
            return self
        if self.iscn.proposal_status == ISCNProposalStatus.LEGACY_UNSPECIFIED:
            raise ValueError("PipelineResult 0.3.0 requires a structured ISCN proposal state")
        if self.schema_version == "0.3.0":
            if self.iscn.genome_build != self.manifest.assay.genome_build:
                raise ValueError("ISCN proposal genome build must match the run manifest")
            resource = self.iscn.resource_provenance
            assessed = self.iscn.proposal_status in {
                ISCNProposalStatus.NO_RENDERABLE_CANDIDATE,
                ISCNProposalStatus.PARTIAL_EVENT_LEVEL,
            }
            if assessed and resource is None:
                raise ValueError("assessed ISCN proposal resource provenance is missing")
            if resource is not None:
                if not isinstance(self.reference_context, ResolvedResourceContext):
                    raise ValueError(
                        "ISCN proposal provenance requires a resolved result resource context"
                    )
                context = self.reference_context
                expected_resource_values = {
                    "genome build": (resource.genome_build, context.genome_build),
                    "reference dictionary contract": (
                        resource.reference_dictionary_contract,
                        context.reference_dictionary_contract,
                    ),
                    "reference bundle ID": (
                        resource.reference_bundle_id,
                        context.reference_bundle_id,
                    ),
                    "reference bundle version": (
                        resource.reference_bundle_version,
                        context.reference_bundle_version,
                    ),
                    "reference lock SHA256": (
                        resource.reference_lock_sha256,
                        context.resource_checksums.get("reference.reference_lock"),
                    ),
                    "annotation cache SHA256": (
                        resource.annotation_cache_sha256,
                        context.resource_checksums.get("reference.annotation_cache"),
                    ),
                    "cytoband SHA256": (
                        resource.cytoband_sha256,
                        context.resource_checksums.get("reference.cytobands"),
                    ),
                    "cytoband release": (
                        resource.cytoband_release,
                        context.resource_releases.get("reference.cytobands", "unspecified"),
                    ),
                }
                mismatches = [
                    label
                    for label, (proposal_value, context_value) in expected_resource_values.items()
                    if proposal_value != context_value
                ]
                if mismatches:
                    raise ValueError(
                        "ISCN proposal provenance does not match the result resource context: "
                        + ", ".join(mismatches)
                    )
            event_ids = [item.event_id for item in self.events]
            if len(event_ids) != len(set(event_ids)):
                raise ValueError("Pipeline result event IDs must be unique")
            if set(self.iscn.considered_event_ids) != set(event_ids):
                raise ValueError("ISCN proposal must disposition every result event")
            events_by_id = {item.event_id: item for item in self.events}
            cnv_event_types = {
                EventType.CHROMOSOME_GAIN,
                EventType.CHROMOSOME_LOSS,
                EventType.DELETION,
                EventType.DUPLICATION,
            }
            if any(
                events_by_id[event_id].event_type not in cnv_event_types
                for event_id in self.iscn.cnv_source_event_ids
            ):
                raise ValueError("ISCN CNV source IDs must refer only to typed CNV event types")
            for fragment in self.iscn.event_fragments:
                source = events_by_id[fragment.event_id]
                if (
                    fragment.event_type != source.event_type
                    or fragment.confidence != source.confidence
                    or fragment.reportable != source.reportable
                ):
                    raise ValueError("ISCN fragment provenance must match its source event")
                expected_fragment = render_supported_cnv_fragment(
                    event_type=source.event_type.value,
                    chromosome=source.primary.chromosome,
                    copy_number=source.copy_number,
                    cytoband_start=source.primary.cytoband_start,
                    cytoband_end=source.primary.cytoband_end,
                    whole_chromosome_span_confirmed=source.whole_chromosome_span_confirmed,
                )
                if fragment.fragment != expected_fragment:
                    raise ValueError(
                        "ISCN fragment notation must match the checksum-bound source event"
                    )
            iscn_outcomes = [item for item in self.modules if item.module == AnalysisModule.ISCN]
            if len(iscn_outcomes) != 1:
                raise ValueError("PipelineResult 0.3.0 requires exactly one ISCN module outcome")
            cnv_outcomes = [item for item in self.modules if item.module == AnalysisModule.CNV]
            if len(cnv_outcomes) != 1:
                raise ValueError("PipelineResult 0.3.0 requires exactly one CNV module outcome")
            expected_status = {
                ISCNProposalStatus.NOT_REQUESTED: ModuleRunStatus.NOT_RUN,
                ISCNProposalStatus.NOT_ASSESSED: ModuleRunStatus.NOT_RUN,
                ISCNProposalStatus.NO_RENDERABLE_CANDIDATE: ModuleRunStatus.NO_CALL,
                ISCNProposalStatus.PARTIAL_EVENT_LEVEL: ModuleRunStatus.COMPLETED,
            }[self.iscn.proposal_status]
            if iscn_outcomes[0].status != expected_status:
                raise ValueError("ISCN proposal and module execution statuses disagree")
            if self.iscn.proposal_status == ISCNProposalStatus.PARTIAL_EVENT_LEVEL:
                if cnv_outcomes[0].status != ModuleRunStatus.COMPLETED:
                    raise ValueError("a partial ISCN proposal requires completed CNV execution")
            elif (
                self.iscn.proposal_status == ISCNProposalStatus.NO_RENDERABLE_CANDIDATE
                and cnv_outcomes[0].status
                not in {
                    ModuleRunStatus.COMPLETED,
                    ModuleRunStatus.NO_CALL,
                }
            ):
                raise ValueError("an assessed empty ISCN proposal requires assessed CNV execution")
        return self
