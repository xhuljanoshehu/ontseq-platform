from __future__ import annotations

import math
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from .models import BenchmarkThresholds, EventType, GenomeBuild, GenomicEvent, StrictModel

SHA256 = r"^[0-9a-f]{64}$"
CNV_EVENT_TYPES = frozenset(
    {
        EventType.CHROMOSOME_GAIN,
        EventType.CHROMOSOME_LOSS,
        EventType.DELETION,
        EventType.DUPLICATION,
    }
)


class CnvStratificationRole(StrEnum):
    PRIMARY = "primary"
    SECONDARY = "secondary"
    EXPLORATORY = "exploratory"


class CnvStratificationDimension(StrEnum):
    CALLER = "caller"
    GENOME_BUILD = "genome_build"
    DATA_BASIS = "data_basis"
    COVERAGE = "coverage"
    TUMOR_FRACTION = "tumor_fraction"
    EVENT_CLASS = "event_class"
    EVENT_SIZE = "event_size"
    BIN_SIZE = "bin_size"


class CnvDataBasis(StrEnum):
    LCWGS_GENOME_WIDE = "lcwgs_genome_wide"
    ADAPTIVE_SAMPLING_OFF_TARGET = "adaptive_sampling_off_target"
    ADAPTIVE_SAMPLING_ON_TARGET = "adaptive_sampling_on_target"


class CnvCallerLock(StrictModel):
    caller_id: str = Field(pattern=r"^[a-z][a-z0-9_]{1,63}$")
    caller_version: str = Field(min_length=1)
    adapter_policy_sha256: str = Field(pattern=SHA256)
    execution_identity_sha256: str = Field(pattern=SHA256)


class CnvStratificationPlan(StrictModel):
    primary_dimensions: list[CnvStratificationDimension] = Field(min_length=1)
    secondary_dimensions: list[CnvStratificationDimension] = Field(default_factory=list)
    exploratory_dimensions: list[CnvStratificationDimension] = Field(default_factory=list)
    coverage_cutpoints_x: list[float] = Field(default_factory=list)
    tumor_fraction_cutpoints: list[float] = Field(default_factory=list)
    event_size_cutpoints_bp: list[int] = Field(default_factory=list)

    @model_validator(mode="after")
    def coherent_strata(self) -> CnvStratificationPlan:
        groups = [
            self.primary_dimensions,
            self.secondary_dimensions,
            self.exploratory_dimensions,
        ]
        flattened = [item for group in groups for item in group]
        if len(flattened) != len(set(flattened)):
            raise ValueError("A stratification dimension may have exactly one role")
        if self.coverage_cutpoints_x != sorted(set(self.coverage_cutpoints_x)):
            raise ValueError("Coverage cutpoints must be unique and strictly sorted")
        if any(not math.isfinite(x) or x <= 0 for x in self.coverage_cutpoints_x):
            raise ValueError("Coverage cutpoints must be finite and greater than zero")
        if self.tumor_fraction_cutpoints != sorted(set(self.tumor_fraction_cutpoints)):
            raise ValueError("Tumour-fraction cutpoints must be unique and strictly sorted")
        if any(not math.isfinite(x) or x <= 0 or x >= 1 for x in self.tumor_fraction_cutpoints):
            raise ValueError("Tumour-fraction cutpoints must lie strictly between zero and one")
        if self.event_size_cutpoints_bp != sorted(set(self.event_size_cutpoints_bp)):
            raise ValueError("Event-size cutpoints must be unique and strictly sorted")
        if any(x <= 0 for x in self.event_size_cutpoints_bp):
            raise ValueError("Event-size cutpoints must be positive")
        return self


class CnvAcceptanceMetric(StrEnum):
    SENSITIVITY = "sensitivity"
    PRECISION = "precision"
    F1 = "f1"
    FALSE_POSITIVE_BURDEN = "false_positive_burden"
    NO_CALL_RATE = "no_call_rate"
    TECHNICAL_FAILURE_RATE = "technical_failure_rate"
    NOT_ASSESSABLE_RATE = "not_assessable_rate"
    SPECIFICITY = "specificity"
    COPY_NUMBER_ERROR = "copy_number_error"
    CELLULARITY_ERROR = "cellularity_error"
    PLOIDY_ERROR = "ploidy_error"
    REPRODUCIBILITY = "reproducibility"


class CnvAcceptanceQuestion(StrictModel):
    question_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    metric: CnvAcceptanceMetric
    minimum_evaluable_denominator: int = Field(ge=1)
    minimum_acceptable: float | None = None
    maximum_acceptable: float | None = None
    filters: dict[str, str | int | float | bool] = Field(default_factory=dict)

    @model_validator(mode="after")
    def bounded_question(self) -> CnvAcceptanceQuestion:
        bounds = [self.minimum_acceptable, self.maximum_acceptable]
        if all(value is None for value in bounds):
            raise ValueError("Acceptance question requires a minimum or maximum bound")
        if any(value is not None and not math.isfinite(value) for value in bounds):
            raise ValueError("Acceptance bounds must be finite")
        if (
            self.minimum_acceptable is not None
            and self.maximum_acceptable is not None
            and self.minimum_acceptable > self.maximum_acceptable
        ):
            raise ValueError("Acceptance minimum cannot exceed maximum")
        return self


class CnvValidationMatrix(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    matrix_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    callers: list[CnvCallerLock] = Field(min_length=1)
    genome_builds: list[GenomeBuild] = Field(min_length=1)
    data_bases: list[CnvDataBasis] = Field(min_length=1)
    stratification: CnvStratificationPlan
    event_classes: list[EventType] = Field(min_length=1)
    qdnaseq_bin_sizes_kbp: list[int] = Field(default_factory=list)
    matching_thresholds: BenchmarkThresholds = Field(default_factory=BenchmarkThresholds)
    acceptance: list[CnvAcceptanceQuestion] = Field(min_length=1)
    retain_all_evidence: Literal[True] = True
    rationale: str = Field(min_length=20)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_matrix(self) -> CnvValidationMatrix:
        caller_ids = [item.caller_id for item in self.callers]
        if len(caller_ids) != len(set(caller_ids)):
            raise ValueError("Caller IDs must be unique")
        if len(self.genome_builds) != len(set(self.genome_builds)):
            raise ValueError("Genome builds must be unique")
        if len(self.data_bases) != len(set(self.data_bases)):
            raise ValueError("CNV data bases must be unique")
        if any(event not in CNV_EVENT_TYPES for event in self.event_classes):
            raise ValueError("CNV validation matrix contains a non-CNV event class")
        if len(self.event_classes) != len(set(self.event_classes)):
            raise ValueError("CNV event classes must be unique")
        if self.qdnaseq_bin_sizes_kbp != sorted(set(self.qdnaseq_bin_sizes_kbp)):
            raise ValueError("QDNAseq bin sizes must be unique and sorted")
        if any(value <= 0 for value in self.qdnaseq_bin_sizes_kbp):
            raise ValueError("QDNAseq bin sizes must be positive")

        dimensions = {
            *self.stratification.primary_dimensions,
            *self.stratification.secondary_dimensions,
            *self.stratification.exploratory_dimensions,
        }
        if (
            CnvStratificationDimension.COVERAGE in dimensions
            and not self.stratification.coverage_cutpoints_x
        ):
            raise ValueError("Registered coverage stratification requires explicit cutpoints")
        if (
            CnvStratificationDimension.TUMOR_FRACTION in dimensions
            and not self.stratification.tumor_fraction_cutpoints
        ):
            raise ValueError(
                "Registered tumour-fraction stratification requires explicit cutpoints"
            )
        if (
            CnvStratificationDimension.EVENT_SIZE in dimensions
            and not self.stratification.event_size_cutpoints_bp
        ):
            raise ValueError("Registered event-size stratification requires explicit cutpoints")
        if "qdnaseq_ace" in caller_ids and not self.qdnaseq_bin_sizes_kbp:
            raise ValueError("QDNAseq+ACE validation requires explicit registered bin sizes")

        question_ids = [item.question_id for item in self.acceptance]
        if len(question_ids) != len(set(question_ids)):
            raise ValueError("Acceptance question IDs must be unique")
        return self


class CnvRepeatKind(StrEnum):
    INDEPENDENT = "independent"
    WITHIN_RUN = "within_run"
    BETWEEN_RUN = "between_run"


class CnvTruthSource(StrictModel):
    method_name: str = Field(min_length=1)
    method_version: str = Field(min_length=1)
    resource_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,127}$")
    resource_sha256: str = Field(pattern=SHA256)
    provenance_reference: str = Field(min_length=3)
    orthogonal_to_evaluated_caller: Literal[True] = True


class CnvAssessabilityMask(StrictModel):
    resource_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,127}$")
    resource_sha256: str = Field(pattern=SHA256)
    unit: Literal["base_pairs", "genomic_bins", "regions"]
    definition: str = Field(min_length=20)


class CnvNegativeUniverse(StrictModel):
    universe_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,127}$")
    unit: Literal["genomic_bins", "regions"]
    assessable_units: int = Field(ge=1)
    resource_sha256: str = Field(pattern=SHA256)
    assessability_mask_sha256: str = Field(pattern=SHA256)
    definition: str = Field(min_length=20)


class CnvValidationSpecimen(StrictModel):
    specimen_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    biological_specimen_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    repeat_kind: CnvRepeatKind
    repeat_group_id: str | None = Field(
        default=None,
        pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$",
    )
    access_basis: Literal["public", "institutionally_authorized"]
    pseudonymized: Literal[True] = True
    material_type: str = Field(min_length=3)
    genome_build: GenomeBuild
    data_basis: CnvDataBasis
    reference_id: str = Field(min_length=1)
    reference_sha256: str = Field(pattern=SHA256)
    input_sha256: str = Field(pattern=SHA256)
    coverage_x: float | None = Field(default=None, ge=0)
    coverage_definition: str | None = None
    tumor_fraction: float | None = Field(default=None, ge=0, le=1)
    tumor_fraction_method: str | None = None
    tumor_fraction_timepoint: str | None = None
    truth_sources: list[CnvTruthSource] = Field(min_length=1)
    assessability_mask: CnvAssessabilityMask
    truth_events: list[GenomicEvent] = Field(default_factory=list)
    negative_universe: CnvNegativeUniverse | None = None

    @model_validator(mode="after")
    def coherent_specimen(self) -> CnvValidationSpecimen:
        if self.coverage_x is not None and not math.isfinite(self.coverage_x):
            raise ValueError("Coverage must be finite")
        if self.coverage_x is None and self.coverage_definition is not None:
            raise ValueError("Coverage definition cannot exist without a measured coverage value")
        if self.coverage_x is not None and not self.coverage_definition:
            raise ValueError("Measured coverage requires its definition")
        if self.tumor_fraction is not None and not math.isfinite(self.tumor_fraction):
            raise ValueError("Tumour fraction must be finite")
        fraction_metadata = [self.tumor_fraction_method, self.tumor_fraction_timepoint]
        if self.tumor_fraction is None and any(value is not None for value in fraction_metadata):
            raise ValueError("Tumour-fraction metadata cannot imply a missing fraction")
        if self.tumor_fraction is not None and any(not value for value in fraction_metadata):
            raise ValueError("Measured tumour fraction requires method and timepoint")
        if any(event.event_type not in CNV_EVENT_TYPES for event in self.truth_events):
            raise ValueError("Truth events must be CNV event types")
        event_ids = [event.event_id for event in self.truth_events]
        if len(event_ids) != len(set(event_ids)):
            raise ValueError("Truth event IDs must be unique within a specimen")
        if (
            self.negative_universe is not None
            and self.negative_universe.assessability_mask_sha256
            != self.assessability_mask.resource_sha256
        ):
            raise ValueError(
                "Negative universe must be derived from the registered assessability mask"
            )
        if self.repeat_kind == CnvRepeatKind.INDEPENDENT and self.repeat_group_id is not None:
            raise ValueError("Independent specimens cannot declare a repeat group")
        if self.repeat_kind != CnvRepeatKind.INDEPENDENT and self.repeat_group_id is None:
            raise ValueError("Repeat specimens require repeat_group_id")
        return self


class CnvValidationCohort(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    specimens: list[CnvValidationSpecimen] = Field(min_length=1)
    sensitive_output: Literal[True] = True
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def unique_specimens(self) -> CnvValidationCohort:
        ids = [item.specimen_id for item in self.specimens]
        if len(ids) != len(set(ids)):
            raise ValueError("Validation specimen IDs must be unique")
        return self
