from __future__ import annotations

from datetime import UTC, datetime
from enum import StrEnum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


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


class AssaySpec(StrictModel):
    mode: AssayMode
    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    target_bed: str | None = None
    target_bed_version: str | None = None

    @model_validator(mode="after")
    def adaptive_sampling_requires_bed(self) -> AssaySpec:
        if self.mode == AssayMode.ADAPTIVE_SAMPLING and (
            not self.target_bed or not self.target_bed_version
        ):
            raise ValueError("adaptive_sampling requires target_bed and target_bed_version")
        return self


class AnalysisSpec(StrictModel):
    profile: str = Field(min_length=1)
    modules: list[AnalysisModule]
    parameters: dict[str, Any] = Field(default_factory=dict)


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


class QCMetrics(StrictModel):
    verdict: Verdict
    metrics: dict[str, float | int | str | None]
    warnings: list[str] = Field(default_factory=list)
    failed_gates: list[str] = Field(default_factory=list)


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
    quality: float | None = None
    filters: list[str] = Field(default_factory=list)


class GenomicEvent(StrictModel):
    event_id: str
    event_type: EventType
    primary: Locus
    secondary: Locus | None = None
    copy_number: float | None = Field(default=None, ge=0)
    genes: list[str] = Field(default_factory=list)
    evidence: list[Evidence] = Field(default_factory=list)
    confidence: Literal["high", "moderate", "low", "unclassified"] = "unclassified"
    reportable: bool = False
    notes: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def paired_events_require_secondary_locus(self) -> GenomicEvent:
        if self.event_type in {EventType.TRANSLOCATION, EventType.FUSION} and not self.secondary:
            raise ValueError(f"{self.event_type.value} requires secondary locus")
        return self


class ISCNProposal(StrictModel):
    notation: str
    standard_edition: Literal["ISCN 2024"] = "ISCN 2024"
    conformance_profile: Literal["subset-v0.1-unvalidated"] = "subset-v0.1-unvalidated"
    review_status: ReviewStatus = ReviewStatus.REVIEW_REQUIRED
    source_event_ids: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)


class ToolRecord(StrictModel):
    name: str
    version: str
    parameters: dict[str, Any] = Field(default_factory=dict)
    container_digest: str | None = None


class Provenance(StrictModel):
    pipeline_version: str
    git_commit: str
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    tools: list[ToolRecord] = Field(default_factory=list)
    reference_checksums: dict[str, str] = Field(default_factory=dict)


class PipelineResult(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    manifest: SampleManifest
    qc: QCMetrics
    events: list[GenomicEvent]
    iscn: ISCNProposal
    provenance: Provenance
    warnings: list[str] = Field(default_factory=list)
    release_status: ReviewStatus = ReviewStatus.REVIEW_REQUIRED

    @model_validator(mode="after")
    def failed_qc_cannot_be_reviewed(self) -> PipelineResult:
        if self.qc.verdict == Verdict.FAIL and self.release_status == ReviewStatus.REVIEWED:
            raise ValueError("A QC-failed result cannot be marked REVIEWED")
        return self
