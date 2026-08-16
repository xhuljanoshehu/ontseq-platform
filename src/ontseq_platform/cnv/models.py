"""Versioned contracts for CNV truth sets, call sets and evaluation reports.

The contract layer is deliberately thin. All comparison logic lives in
:mod:`ontseq_platform.cnv.core`, which has no pydantic dependency, so that the science
can be tested independently of serialization and so that a schema change cannot quietly
alter a metric.

Two fields deserve attention because they carry most of the scientific meaning:

``background_state``
    What a set means by staying silent about a region. A SNP array asserts a neutral
    copy number wherever it has probes; a cytogenetic report and an alteration-only
    caller assert nothing outside what they list. Treating an open-world source as
    closed-world manufactures false positives; the reverse hides them.

``resolution_bp``
    The smallest event the source could have detected at all. Truth sets are silent, not
    negative, below their own resolution, and an evaluation must not charge a caller with
    a false positive for finding something the truth method could never have seen.
"""

from __future__ import annotations

from datetime import UTC, datetime
from enum import StrEnum
from typing import Annotated, Literal

from pydantic import Field, model_validator

from ..models import (
    EventType,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    StrictModel,
    ToolRecord,
)
from .core import BoundaryUncertainty, StateSegment
from .states import ConcordanceMode, CopyNumberState

CONTIG_PATTERN = r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$"

#: States a truth or call set may use as its background. Restricted to the two that
#: carry a defensible meaning; any other value would be a category error.
BackgroundState = Literal[CopyNumberState.NEUTRAL, CopyNumberState.NO_CALL]


class CnvTruthSource(StrEnum):
    """Where a truth set came from.

    The source determines how much the truth can support. Cytogenetic banding resolves
    breakpoints to a band; arrays resolve them to probe spacing; a simulator knows them
    exactly. Recording the source is what lets the evaluator refuse to report a metric
    the truth cannot justify.
    """

    SIMULATED = "simulated"
    ISCN_KARYOTYPE = "iscn_karyotype"
    FISH = "fish"
    SNP_ARRAY = "snp_array"
    CGH_ARRAY = "cgh_array"
    SHORT_READ_WGS = "short_read_wgs"
    SHORT_READ_PANEL = "short_read_panel"
    ORTHOGONAL_LONG_READ = "orthogonal_long_read"
    PUBLIC_REFERENCE_MATERIAL = "public_reference_material"
    EXPERT_CONSENSUS = "expert_consensus"


class CnvDataBasis(StrEnum):
    """Which reads a copy-number estimate was derived from.

    Adaptive sampling produces two very different read populations in one run. Rejected
    reads still occupy the flow cell and form a near-uniform low-coverage whole-genome
    background, which is the population most depth-based CNV methods actually assume.
    On-target reads are deeply but non-uniformly enriched and violate that assumption.
    Mixing the two in one benchmark stratum compares incomparable things.
    """

    ADAPTIVE_SAMPLING_OFF_TARGET = "adaptive_sampling_off_target"
    ADAPTIVE_SAMPLING_ON_TARGET = "adaptive_sampling_on_target"
    ADAPTIVE_SAMPLING_COMBINED = "adaptive_sampling_combined"
    LOW_COVERAGE_WGS = "low_coverage_wgs"
    WHOLE_GENOME = "whole_genome"
    SIMULATED = "simulated"


class CnvSegment(StrictModel):
    """One contiguous copy-number call or truth interval.

    Coordinates are zero-based half-open, matching BED. Several quantitative fields are
    optional because different methods report on different scales; the evaluator uses
    ``state`` for agreement and ``copy_number`` for quantitative error, and simply
    reports fewer metrics when a field is absent rather than imputing one.
    """

    contig: str = Field(pattern=CONTIG_PATTERN)
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    state: CopyNumberState
    copy_number: float | None = Field(default=None, ge=0)
    copy_ratio: float | None = Field(default=None, ge=0)
    log2_ratio: float | None = None
    b_allele_frequency: float | None = Field(default=None, ge=0, le=1)
    cellular_fraction: float | None = Field(default=None, ge=0, le=1)
    supporting_bins: int | None = Field(default=None, ge=0)
    quality: float | None = Field(default=None, ge=0)
    start_uncertainty_bp: int = Field(default=0, ge=0)
    end_uncertainty_bp: int = Field(default=0, ge=0)
    notes: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def end_after_start(self) -> CnvSegment:
        if self.end <= self.start:
            raise ValueError("segment end must be greater than start")
        return self

    @property
    def length_bp(self) -> int:
        return self.end - self.start

    def to_state_segment(self, *, carry_uncertainty: bool = False) -> StateSegment:
        """Project onto the dependency-free comparison type.

        Boundary uncertainty is carried only for truth segments. A caller's own opinion
        about its breakpoint precision must not be allowed to suppress the breakpoint
        metric that is measuring it.
        """
        return StateSegment(
            contig=self.contig,
            start=self.start,
            end=self.end,
            state=self.state,
            copy_number=self.copy_number,
            start_uncertainty_bp=self.start_uncertainty_bp if carry_uncertainty else None,
            end_uncertainty_bp=self.end_uncertainty_bp if carry_uncertainty else None,
        )

    def boundary_uncertainty(self) -> BoundaryUncertainty:
        return BoundaryUncertainty(
            start_uncertainty_bp=self.start_uncertainty_bp,
            end_uncertainty_bp=self.end_uncertainty_bp,
        )


class GenomicRegion(StrictModel):
    """A plain interval used for scopes, exclusions and informative regions."""

    contig: str = Field(pattern=CONTIG_PATTERN)
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    label: str | None = None

    @model_validator(mode="after")
    def end_after_start(self) -> GenomicRegion:
        if self.end <= self.start:
            raise ValueError("region end must be greater than start")
        return self


class CnvTruthSet(StrictModel):
    """A reference answer against which call sets are scored."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    truth_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    sample_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    genome_build: GenomeBuild
    source: CnvTruthSource
    source_version: str = Field(min_length=1)
    background_state: BackgroundState
    #: Smallest event this source could resolve. Below it the truth is silent.
    resolution_bp: int = Field(ge=0)
    segments: list[CnvSegment] = Field(default_factory=list)
    #: Regions the source could actually interrogate. Empty means "the whole reference",
    #: which is only defensible for a closed-world genome-wide source.
    informative_regions: list[GenomicRegion] = Field(default_factory=list)
    #: Regions the source explicitly could not assess.
    uninformative_regions: list[GenomicRegion] = Field(default_factory=list)
    tumor_fraction: float | None = Field(default=None, ge=0, le=1)
    baseline_ploidy: float = Field(default=2.0, gt=0)
    fingerprint: FileFingerprint | None = None
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def open_world_truth_needs_scope_or_resolution(self) -> CnvTruthSet:
        probe_scoped_sources = {CnvTruthSource.FISH, CnvTruthSource.SHORT_READ_PANEL}
        if (
            self.background_state == CopyNumberState.NEUTRAL
            and not self.informative_regions
            and self.source in probe_scoped_sources
        ):
            raise ValueError(
                f"{self.source.value} truth cannot claim a genome-wide neutral "
                "background; declare informative_regions instead"
            )
        if self.background_state == CopyNumberState.NEUTRAL and self.resolution_bp == 0:
            raise ValueError(
                "a closed-world truth set must declare resolution_bp so that events "
                "below its detection limit are not scored as false positives"
            )
        return self

    def state_segments(self) -> list[StateSegment]:
        return [segment.to_state_segment(carry_uncertainty=True) for segment in self.segments]


class CnvCallSet(StrictModel):
    """Normalized output of one CNV method for one sample."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    call_set_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    sample_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    genome_build: GenomeBuild
    method: str = Field(min_length=1)
    method_version: str = Field(min_length=1)
    data_basis: CnvDataBasis
    background_state: BackgroundState
    status: ModuleRunStatus
    segments: list[CnvSegment] = Field(default_factory=list)
    #: Regions the method declined to call. Kept separate from segments so that a
    #: no-call is never confused with a neutral call.
    no_call_regions: list[GenomicRegion] = Field(default_factory=list)
    bin_size_bp: int | None = Field(default=None, gt=0)
    estimated_tumor_fraction: float | None = Field(default=None, ge=0, le=1)
    estimated_ploidy: float | None = Field(default=None, gt=0)
    mean_coverage_x: float | None = Field(default=None, ge=0)
    tool: ToolRecord | None = None
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    reportable: Literal[False] = False
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def status_matches_content(self) -> CnvCallSet:
        if self.status == ModuleRunStatus.COMPLETED and not self.segments:
            raise ValueError(
                "a COMPLETED call set must contain segments; an empty result is NO_CALL"
            )
        if self.status == ModuleRunStatus.NO_CALL and self.segments:
            raise ValueError("a NO_CALL call set must not contain segments")
        return self

    def state_segments(self) -> list[StateSegment]:
        return [segment.to_state_segment() for segment in self.segments]


class CnvStrata(StrictModel):
    """Typed stratification keys for cross-run aggregation.

    Typed rather than a free-form mapping so that aggregation across runs is possible at
    all. A benchmark that cannot be grouped by coverage and tumor fraction cannot answer
    the questions this project needs answered.
    """

    assay_mode: str | None = None
    data_basis: CnvDataBasis | None = None
    mean_coverage_x: float | None = Field(default=None, ge=0)
    tumor_fraction: float | None = Field(default=None, ge=0, le=1)
    bin_size_bp: int | None = Field(default=None, gt=0)
    replicate: int | None = Field(default=None, ge=0)
    sample_class: str | None = None
    notes: list[str] = Field(default_factory=list)


class CnvEvaluationOptions(StrictModel):
    """Serializable mirror of :class:`ontseq_platform.cnv.core.EvaluationOptions`."""

    concordance_mode: ConcordanceMode = ConcordanceMode.DIRECTIONAL
    detection_overlap_fraction: float = Field(default=0.5, gt=0, le=1)
    minimum_assessable_fraction: float = Field(default=0.5, gt=0, le=1)
    copy_number_tolerance: float = Field(default=0.5, ge=0)
    maximum_truth_boundary_uncertainty_bp: int = Field(default=1_000_000, ge=0)
    confidence_level: float = Field(default=0.95, gt=0, lt=1)
    minimum_depth_floor_x: float | None = Field(default=None, ge=0)
    analysis_scope_flank_bp: int = Field(default=0, ge=0)


class ProportionResult(StrictModel):
    """A proportion reported together with its interval and denominator."""

    successes: int = Field(ge=0)
    total: int = Field(ge=0)
    point: float | None = Field(default=None, ge=0, le=1)
    lower: float | None = Field(default=None, ge=0, le=1)
    upper: float | None = Field(default=None, ge=0, le=1)
    confidence_level: float = Field(gt=0, lt=1)


class GenomePartitionReport(StrictModel):
    """How the genome was divided before scoring.

    The counts reconcile exactly, and the contract enforces it, so a reader can always
    verify what fraction of the genome a metric actually applied to:

    ``reference_bases == mask_bases + excluded_bases``

    ``mask_bases == evaluable_bases + truth_silent_bases + query_no_call_bases``
    """

    reference_bases: int = Field(ge=0)
    mask_bases: int = Field(ge=0)
    evaluable_bases: int = Field(ge=0)
    excluded_bases: int = Field(ge=0)
    query_no_call_bases: int = Field(ge=0)
    truth_silent_bases: int = Field(ge=0)
    evaluable_fraction: float | None = Field(default=None, ge=0, le=1)
    excluded_bases_by_reason: dict[str, int] = Field(default_factory=dict)

    @model_validator(mode="after")
    def partition_reconciles(self) -> GenomePartitionReport:
        if self.mask_bases != (
            self.evaluable_bases + self.truth_silent_bases + self.query_no_call_bases
        ):
            raise ValueError(
                "genome partition does not reconcile: mask bases must equal evaluable "
                "plus truth-silent plus query-no-call bases"
            )
        if self.reference_bases != self.mask_bases + self.excluded_bases:
            raise ValueError(
                "genome partition does not reconcile: reference bases must equal mask "
                "plus excluded bases"
            )
        return self


class BaseLevelReport(StrictModel):
    """Base-pair-weighted agreement.

    ``recall_by_state`` and ``precision_by_state`` are descriptive point estimates
    without intervals on purpose. Base pairs within one segment are not independent
    observations, so a binomial interval computed over base counts would be absurdly
    narrow. Intervals are reported only at event level, where the unit of observation is
    defensible.
    """

    evaluable_bases: int = Field(ge=0)
    concordant_bases: int = Field(ge=0)
    concordance: float | None = Field(default=None, ge=0, le=1)
    confusion: dict[str, int] = Field(default_factory=dict)
    truth_bases_by_state: dict[str, int] = Field(default_factory=dict)
    query_bases_by_state: dict[str, int] = Field(default_factory=dict)
    recall_by_state: dict[str, float | None] = Field(default_factory=dict)
    precision_by_state: dict[str, float | None] = Field(default_factory=dict)


class EventOutcomeReport(StrictModel):
    """Per-event scoring outcome, retained in full for auditability."""

    event_id: str
    contig: str
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    state: CopyNumberState
    size_class: str
    length_bp: int = Field(gt=0)
    outcome: str
    reason: str
    total_bases: int = Field(ge=0)
    evaluable_bases: int = Field(ge=0)
    concordant_bases: int = Field(ge=0)
    concordant_fraction: float | None = Field(default=None, ge=0, le=1)
    start_delta_bp: int | None = None
    end_delta_bp: int | None = None
    breakpoint_skip_reason: str | None = None


class StratumDetectionReport(StrictModel):
    """Detection within one stratum."""

    label: str
    detected: int = Field(ge=0)
    missed: int = Field(ge=0)
    not_assessable: int = Field(ge=0)
    detection_rate: ProportionResult


class CopyNumberAccuracyReport(StrictModel):
    """Quantitative copy-number agreement."""

    assessed_bases: int = Field(ge=0)
    mean_absolute_error: float | None = Field(default=None, ge=0)
    root_mean_square_error: float | None = Field(default=None, ge=0)
    within_tolerance_bases: int = Field(ge=0)
    within_tolerance_fraction: float | None = Field(default=None, ge=0, le=1)


class BreakpointAccuracyReport(StrictModel):
    """Breakpoint agreement over truth events whose resolution supports it."""

    assessed_events: int = Field(ge=0)
    skipped_events: int = Field(ge=0)
    skip_reasons: dict[str, int] = Field(default_factory=dict)
    median_absolute_start_delta_bp: int | None = Field(default=None, ge=0)
    median_absolute_end_delta_bp: int | None = Field(default=None, ge=0)
    maximum_absolute_delta_bp: int | None = Field(default=None, ge=0)


class CnvEvaluationReport(StrictModel):
    """The complete, auditable result of comparing one call set to one truth set."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    evaluation_id: str = Field(min_length=1)
    sample_id: str
    truth_id: str
    call_set_id: str
    genome_build: GenomeBuild
    method: str
    method_version: str
    truth_source: CnvTruthSource
    data_basis: CnvDataBasis
    strata: CnvStrata = Field(default_factory=CnvStrata)
    options: CnvEvaluationOptions
    partition: GenomePartitionReport
    base_level: BaseLevelReport
    detection_rate: ProportionResult
    confirmation_rate: ProportionResult
    detection_by_size_class: list[StratumDetectionReport] = Field(default_factory=list)
    detection_by_state: list[StratumDetectionReport] = Field(default_factory=list)
    copy_number_accuracy: CopyNumberAccuracyReport
    breakpoint_accuracy: BreakpointAccuracyReport
    truth_events: list[EventOutcomeReport] = Field(default_factory=list)
    query_events: list[EventOutcomeReport] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))


_STATE_TO_EVENT_TYPE: dict[CopyNumberState, EventType] = {
    CopyNumberState.HOMOZYGOUS_LOSS: EventType.DELETION,
    CopyNumberState.LOSS: EventType.DELETION,
    CopyNumberState.GAIN: EventType.DUPLICATION,
    CopyNumberState.HIGH_AMPLIFICATION: EventType.DUPLICATION,
}


def segment_to_genomic_event(
    segment: CnvSegment,
    *,
    event_id: str,
    contig_length: int | None = None,
    whole_chromosome_fraction: float = 0.9,
) -> GenomicEvent | None:
    """Project a CNV segment onto the shared :class:`GenomicEvent` contract.

    Returns ``None`` for states the existing event vocabulary cannot express, most
    importantly copy-neutral LOH. Mapping it onto a deletion or duplication would assert
    a dosage change that does not exist, so the segment is dropped from the shared
    contract and remains available only in the CNV-native report.

    A segment covering most of a contig is promoted to a whole-chromosome event so that
    downstream ISCN logic sees the correct granularity.
    """
    event_type = _STATE_TO_EVENT_TYPE.get(segment.state)
    if event_type is None:
        return None
    if contig_length and segment.length_bp >= whole_chromosome_fraction * contig_length:
        event_type = (
            EventType.CHROMOSOME_GAIN
            if event_type == EventType.DUPLICATION
            else EventType.CHROMOSOME_LOSS
        )
    notes = list(segment.notes)
    notes.append(
        "Derived from a benchmark-stage CNV segment. Not reportable without assay validation."
    )
    return GenomicEvent(
        event_id=event_id,
        event_type=event_type,
        primary=Locus(chromosome=segment.contig, start=segment.start, end=segment.end),
        length_bp=segment.length_bp,
        copy_number=segment.copy_number,
        confidence="unclassified",
        reportable=False,
        notes=notes,
    )


class CnvBenchmarkCase(StrictModel):
    """A self-contained, reproducible comparison of one call set against one truth set."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    case_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    genome_build: GenomeBuild
    contig_lengths: dict[str, Annotated[int, Field(gt=0)]] = Field(min_length=1)
    truth: CnvTruthSet
    call_set: CnvCallSet
    options: CnvEvaluationOptions = Field(default_factory=CnvEvaluationOptions)
    strata: CnvStrata = Field(default_factory=CnvStrata)
    analysis_scope: list[GenomicRegion] = Field(default_factory=list)
    excluded_regions: list[GenomicRegion] = Field(default_factory=list)
    exclusion_reason: str = "blacklist"
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def case_is_internally_consistent(self) -> CnvBenchmarkCase:
        if self.truth.genome_build != self.genome_build:
            raise ValueError("truth genome build does not match the benchmark case")
        if self.call_set.genome_build != self.genome_build:
            raise ValueError("call set genome build does not match the benchmark case")
        if self.truth.sample_id != self.call_set.sample_id:
            raise ValueError("truth and call set refer to different samples")
        return self
