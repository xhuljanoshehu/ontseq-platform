from __future__ import annotations

from datetime import UTC, datetime
from enum import StrEnum
from typing import Annotated, Any, Literal

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


class AmlRearrangementRecord(StrictModel):
    record_id: str = Field(min_length=1)
    pattern_type: Literal["exact_pair", "gene_any_partner"]
    genes: list[str] = Field(min_length=1, max_length=2)
    display_name: str = Field(min_length=1)
    relevance: Literal["aml_defining_pattern", "aml_relevant_pattern"]
    source_ids: list[str] = Field(min_length=1)
    caveat: str = Field(min_length=1)

    @model_validator(mode="after")
    def pattern_has_expected_gene_count(self) -> AmlRearrangementRecord:
        expected = 2 if self.pattern_type == "exact_pair" else 1
        if len(self.genes) != expected:
            raise ValueError(f"{self.pattern_type} requires exactly {expected} gene(s)")
        if len({gene.upper() for gene in self.genes}) != len(self.genes):
            raise ValueError("AML rearrangement genes must be unique")
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


class GenomicEvent(StrictModel):
    event_id: str
    event_type: EventType
    primary: Locus
    secondary: Locus | None = None
    length_bp: int | None = Field(default=None, ge=1)
    copy_number: float | None = Field(default=None, ge=0)
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
    fusion_status: FusionSupportStatus = FusionSupportStatus.NOT_ASSESSED
    validation_status: SvValidationStatus = SvValidationStatus.DETECTED
    #: Knowledge-base records matching this event. Evidence for a reviewer; deliberately
    #: without influence on ``confidence`` or ``reportable``.
    annotations: list[EventAnnotation] = Field(default_factory=list)

    @model_validator(mode="after")
    def paired_events_require_secondary_locus(self) -> GenomicEvent:
        if self.event_type in {EventType.TRANSLOCATION, EventType.FUSION} and not self.secondary:
            raise ValueError(f"{self.event_type.value} requires secondary locus")
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
    schema_version: Literal["0.1.0"] = "0.1.0"
    manifest: SampleManifest
    qc: QCMetrics
    events: list[GenomicEvent]
    iscn: ISCNProposal
    provenance: Provenance
    modules: list[ModuleOutcome] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    release_status: ReviewStatus = ReviewStatus.REVIEW_REQUIRED

    @model_validator(mode="after")
    def failed_qc_cannot_be_reviewed(self) -> PipelineResult:
        if self.qc.verdict == Verdict.FAIL and self.release_status == ReviewStatus.REVIEWED:
            raise ValueError("A QC-failed result cannot be marked REVIEWED")
        return self
