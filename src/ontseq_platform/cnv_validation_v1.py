from __future__ import annotations

import hashlib
import math
from itertools import product
from typing import Literal

from pydantic import Field, model_validator

from .cnv_validation_contracts import (
    CnvAcceptanceMetric,
    CnvAcceptanceQuestion,
    CnvAssessabilityMask,
    CnvCallerLock,
    CnvDataBasis,
    CnvRepeatKind,
    CnvStratificationDimension,
    CnvStratificationPlan,
    CnvTruthSource,
    CnvValidationCohort,
    CnvValidationMatrix,
    CnvValidationSpecimen,
)
from .cnv_validation_registration import (
    CnvValidationRegistration,
    preregister_cnv_validation,
)
from .dilution import DilutionPolicy
from .models import (
    BenchmarkThresholds,
    EventType,
    GenomeBuild,
    GenomicEvent,
    Locus,
    StrictModel,
)

_SHA256 = r"^[0-9a-f]{64}$"
_ID = r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$"
_LCWGS_V1_CALLERS = frozenset({"qdnaseq_ace", "spectre"})
_GRCH38_CHR7_LENGTH = 159_345_973
_GRCH38_CHR8_LENGTH = 145_138_636


def _sha_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _float_token(value: float) -> str:
    return format(value, "g").replace(".", "p")


def _midpoints(values: list[float]) -> list[float]:
    return [round((left + right) / 2.0, 12) for left, right in zip(values[:-1], values[1:], strict=True)]


def _integer_midpoints(values: list[int]) -> list[int]:
    result: list[int] = []
    for left, right in zip(values[:-1], values[1:], strict=True):
        total = left + right
        if total % 2:
            raise ValueError("Analytical validation event-size midpoints must be exact integers")
        result.append(total // 2)
    return result


class CnvValidationV1Design(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    design_id: str = Field(pattern=_ID)
    genome_build: Literal[GenomeBuild.GRCH38] = GenomeBuild.GRCH38
    data_basis: Literal[CnvDataBasis.LCWGS_GENOME_WIDE] = CnvDataBasis.LCWGS_GENOME_WIDE
    coverage_levels_x: list[float] = Field(min_length=2)
    tumor_fractions: list[float] = Field(min_length=2)
    subchromosomal_sizes_bp: list[int] = Field(min_length=2)
    replicates: int = Field(ge=1)
    seed: int = Field(ge=0, le=999_999)
    qdnaseq_bin_sizes_kbp: list[int] = Field(default_factory=lambda: [100, 500, 1000])
    reference_id: str = Field(min_length=3)
    reference_sha256: str = Field(pattern=_SHA256)
    truth_source_resource_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,127}$")
    truth_source_sha256: str = Field(pattern=_SHA256)
    assessability_resource_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,127}$")
    assessability_sha256: str = Field(pattern=_SHA256)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_factors(self) -> CnvValidationV1Design:
        if self.coverage_levels_x != sorted(set(self.coverage_levels_x)):
            raise ValueError("Coverage levels must be unique and ordered ascending")
        if any(not math.isfinite(value) or value <= 0 for value in self.coverage_levels_x):
            raise ValueError("Coverage levels must be finite and greater than zero")

        expected_fractions = sorted(set(self.tumor_fractions), reverse=True)
        if self.tumor_fractions != expected_fractions:
            raise ValueError("Tumour fractions must be unique and ordered highest first")
        if any(
            not math.isfinite(value) or value <= 0 or value > 1 for value in self.tumor_fractions
        ):
            raise ValueError("Tumour fractions must be finite and lie in (0, 1]")

        if self.subchromosomal_sizes_bp != sorted(set(self.subchromosomal_sizes_bp)):
            raise ValueError("Subchromosomal CNV sizes must be unique and ordered ascending")
        if any(value <= 0 for value in self.subchromosomal_sizes_bp):
            raise ValueError("Subchromosomal CNV sizes must be positive")

        if self.qdnaseq_bin_sizes_kbp != sorted(set(self.qdnaseq_bin_sizes_kbp)):
            raise ValueError("QDNAseq bin sizes must be unique and ordered")
        if any(value <= 0 for value in self.qdnaseq_bin_sizes_kbp):
            raise ValueError("QDNAseq bin sizes must be positive")
        return self


class CnvValidationV1Cell(StrictModel):
    cell_id: str = Field(pattern=_ID)
    specimen_id: str = Field(pattern=_ID)
    repeat_group_id: str = Field(pattern=_ID)
    coverage_x: float = Field(gt=0)
    tumor_fraction: float = Field(gt=0, le=1)
    event_type: EventType
    event_size_bp: int | None = Field(default=None, gt=0)
    replicate: int = Field(ge=1)
    chromosome: str
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    copy_number: float = Field(ge=0)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_event(self) -> CnvValidationV1Cell:
        whole = self.event_type in {EventType.CHROMOSOME_GAIN, EventType.CHROMOSOME_LOSS}
        if whole and self.event_size_bp is not None:
            raise ValueError("Whole-chromosome validation cells cannot declare a numeric size")
        if not whole:
            if self.event_type not in {EventType.DELETION, EventType.DUPLICATION}:
                raise ValueError("v1 subchromosomal cells support deletion or duplication")
            if self.event_size_bp is None:
                raise ValueError("Subchromosomal validation cells require event_size_bp")
            if self.end - self.start != self.event_size_bp:
                raise ValueError("Subchromosomal validation locus must exactly match event size")
        if self.end <= self.start:
            raise ValueError("Validation truth locus end must exceed start")
        return self


def default_cnv_validation_v1_design() -> CnvValidationV1Design:
    return CnvValidationV1Design(
        design_id="ontseq-analytical-validation-v1",
        coverage_levels_x=[0.5, 1.0, 2.0, 5.0, 10.0, 20.0],
        tumor_fractions=[1.0, 0.5, 0.3, 0.2, 0.1, 0.05],
        subchromosomal_sizes_bp=[
            100_000,
            250_000,
            500_000,
            1_000_000,
            5_000_000,
            10_000_000,
        ],
        replicates=3,
        seed=20_260_920,
        reference_id="synthetic-grch38-validation-v1",
        reference_sha256=_sha_text("synthetic-grch38-validation-v1-reference-contract"),
        truth_source_resource_id="analytical-validation-v1-synthetic-truth",
        truth_source_sha256=_sha_text("analytical-validation-v1-synthetic-truth-contract"),
        assessability_resource_id="analytical-validation-v1-assessability",
        assessability_sha256=_sha_text("analytical-validation-v1-assessability-contract"),
    )


def _event_template(
    event_type: EventType,
    event_size_bp: int | None,
) -> tuple[str, int, int, float]:
    if event_type == EventType.DELETION:
        if event_size_bp is None:
            raise ValueError("Deletion template requires an event size")
        start = 40_000_000
        return "chr5", start, start + event_size_bp, 1.0
    if event_type == EventType.DUPLICATION:
        if event_size_bp is None:
            raise ValueError("Duplication template requires an event size")
        start = 60_000_000
        return "chr7", start, start + event_size_bp, 3.0
    if event_type == EventType.CHROMOSOME_LOSS:
        return "chr7", 0, _GRCH38_CHR7_LENGTH, 1.0
    if event_type == EventType.CHROMOSOME_GAIN:
        return "chr8", 0, _GRCH38_CHR8_LENGTH, 3.0
    raise ValueError(f"Unsupported v1 CNV event type: {event_type.value}")


def _cell_identity(
    *,
    coverage_x: float,
    tumor_fraction: float,
    event_type: EventType,
    event_size_bp: int | None,
) -> str:
    event_token = {
        EventType.DELETION: "del",
        EventType.DUPLICATION: "dup",
        EventType.CHROMOSOME_LOSS: "chr-loss",
        EventType.CHROMOSOME_GAIN: "chr-gain",
    }[event_type]
    size_token = "whole" if event_size_bp is None else str(event_size_bp)
    return (
        f"cov{_float_token(coverage_x)}-tf{_float_token(tumor_fraction)}-"
        f"{event_token}-s{size_token}"
    )


def expand_cnv_validation_v1_cells(design: CnvValidationV1Design) -> list[CnvValidationV1Cell]:
    cells: list[CnvValidationV1Cell] = []
    subchromosomal_types = [EventType.DELETION, EventType.DUPLICATION]
    whole_chromosome_types = [EventType.CHROMOSOME_LOSS, EventType.CHROMOSOME_GAIN]

    for coverage_x, tumor_fraction, event_type, event_size_bp, replicate in product(
        design.coverage_levels_x,
        design.tumor_fractions,
        subchromosomal_types,
        design.subchromosomal_sizes_bp,
        range(1, design.replicates + 1),
    ):
        identity = _cell_identity(
            coverage_x=coverage_x,
            tumor_fraction=tumor_fraction,
            event_type=event_type,
            event_size_bp=event_size_bp,
        )
        chromosome, start, end, copy_number = _event_template(event_type, event_size_bp)
        cells.append(
            CnvValidationV1Cell(
                cell_id=f"{identity}-r{replicate}",
                specimen_id=f"AV1-{identity}-r{replicate}",
                repeat_group_id=f"AV1-{identity}",
                coverage_x=coverage_x,
                tumor_fraction=tumor_fraction,
                event_type=event_type,
                event_size_bp=event_size_bp,
                replicate=replicate,
                chromosome=chromosome,
                start=start,
                end=end,
                copy_number=copy_number,
            )
        )

    for coverage_x, tumor_fraction, event_type, replicate in product(
        design.coverage_levels_x,
        design.tumor_fractions,
        whole_chromosome_types,
        range(1, design.replicates + 1),
    ):
        identity = _cell_identity(
            coverage_x=coverage_x,
            tumor_fraction=tumor_fraction,
            event_type=event_type,
            event_size_bp=None,
        )
        chromosome, start, end, copy_number = _event_template(event_type, None)
        cells.append(
            CnvValidationV1Cell(
                cell_id=f"{identity}-r{replicate}",
                specimen_id=f"AV1-{identity}-r{replicate}",
                repeat_group_id=f"AV1-{identity}",
                coverage_x=coverage_x,
                tumor_fraction=tumor_fraction,
                event_type=event_type,
                event_size_bp=None,
                replicate=replicate,
                chromosome=chromosome,
                start=start,
                end=end,
                copy_number=copy_number,
            )
        )

    ids = [cell.cell_id for cell in cells]
    specimen_ids = [cell.specimen_id for cell in cells]
    if len(ids) != len(set(ids)) or len(specimen_ids) != len(set(specimen_ids)):
        raise ValueError("Analytical validation v1 generated duplicate cell identities")
    return cells


def materialize_synthetic_validation_v1_cohort(
    design: CnvValidationV1Design,
) -> CnvValidationCohort:
    truth_source = CnvTruthSource(
        method_name="deterministic-synthetic-CNV-template",
        method_version="v1",
        resource_id=design.truth_source_resource_id,
        resource_sha256=design.truth_source_sha256,
        provenance_reference=(
            "Repository-generated synthetic Analytical Validation Program v1 truth; "
            "not a biological specimen."
        ),
    )
    assessability = CnvAssessabilityMask(
        resource_id=design.assessability_resource_id,
        resource_sha256=design.assessability_sha256,
        unit="regions",
        definition=(
            "Synthetic locked assessability declaration for Analytical Validation Program v1 "
            "software and study-design testing."
        ),
    )

    specimens: list[CnvValidationSpecimen] = []
    for cell in expand_cnv_validation_v1_cells(design):
        event = GenomicEvent(
            event_id=f"truth-{cell.cell_id}",
            event_type=cell.event_type,
            primary=Locus(
                chromosome=cell.chromosome,
                start=cell.start,
                end=cell.end,
            ),
            copy_number=cell.copy_number,
        )
        specimens.append(
            CnvValidationSpecimen(
                specimen_id=cell.specimen_id,
                biological_specimen_id=f"BIO-{cell.repeat_group_id}",
                repeat_kind=CnvRepeatKind.WITHIN_RUN,
                repeat_group_id=cell.repeat_group_id,
                access_basis="public",
                material_type="synthetic-in-silico-CNV",
                genome_build=design.genome_build,
                data_basis=design.data_basis,
                reference_id=design.reference_id,
                reference_sha256=design.reference_sha256,
                input_sha256=_sha_text(
                    f"{design.design_id}|{cell.cell_id}|seed={design.seed}|synthetic-input"
                ),
                coverage_x=cell.coverage_x,
                coverage_definition=(
                    "Planned nominal mean autosomal coverage for synthetic validation v1."
                ),
                tumor_fraction=cell.tumor_fraction,
                tumor_fraction_method="deterministic in-silico tumour/normal read-mixture target",
                tumor_fraction_timepoint="synthetic same-source design",
                truth_sources=[truth_source],
                assessability_mask=assessability,
                truth_events=[event],
                negative_universe=None,
            )
        )
    return CnvValidationCohort(specimens=specimens)


def _default_acceptance_questions(replicates: int) -> list[CnvAcceptanceQuestion]:
    minimum_denominator = max(3, replicates)
    return [
        CnvAcceptanceQuestion(
            question_id="v1-candidate-sensitivity",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=minimum_denominator,
            minimum_acceptable=0.90,
        ),
        CnvAcceptanceQuestion(
            question_id="v1-candidate-precision",
            metric=CnvAcceptanceMetric.PRECISION,
            minimum_evaluable_denominator=minimum_denominator,
            minimum_acceptable=0.90,
        ),
        CnvAcceptanceQuestion(
            question_id="v1-candidate-fp-burden",
            metric=CnvAcceptanceMetric.FALSE_POSITIVE_BURDEN,
            minimum_evaluable_denominator=minimum_denominator,
            maximum_acceptable=1.0,
        ),
        CnvAcceptanceQuestion(
            question_id="v1-candidate-no-call-rate",
            metric=CnvAcceptanceMetric.NO_CALL_RATE,
            minimum_evaluable_denominator=minimum_denominator,
            maximum_acceptable=0.10,
        ),
        CnvAcceptanceQuestion(
            question_id="v1-candidate-technical-failure-rate",
            metric=CnvAcceptanceMetric.TECHNICAL_FAILURE_RATE,
            minimum_evaluable_denominator=minimum_denominator,
            maximum_acceptable=0.05,
        ),
        CnvAcceptanceQuestion(
            question_id="v1-candidate-copy-number-mae",
            metric=CnvAcceptanceMetric.COPY_NUMBER_ERROR,
            minimum_evaluable_denominator=minimum_denominator,
            maximum_acceptable=0.50,
        ),
    ]


def build_cnv_validation_v1_matrix(
    design: CnvValidationV1Design,
    caller_locks: list[CnvCallerLock],
) -> CnvValidationMatrix:
    if not caller_locks:
        raise ValueError("Analytical validation v1 requires at least one caller lock")
    incompatible = sorted(
        caller.caller_id for caller in caller_locks if caller.caller_id not in _LCWGS_V1_CALLERS
    )
    if incompatible:
        raise ValueError(
            "Analytical validation v1 lcWGS study accepts only compatible genome-wide callers "
            f"{sorted(_LCWGS_V1_CALLERS)}; incompatible caller lock(s): {incompatible}"
        )

    coverage_cutpoints = _midpoints(design.coverage_levels_x)
    tumor_cutpoints = _midpoints(sorted(design.tumor_fractions))
    size_cutpoints = _integer_midpoints(design.subchromosomal_sizes_bp)
    caller_ids = {caller.caller_id for caller in caller_locks}

    return CnvValidationMatrix(
        matrix_id=f"{design.design_id}-matrix",
        callers=caller_locks,
        genome_builds=[design.genome_build],
        data_bases=[design.data_basis],
        stratification=CnvStratificationPlan(
            primary_dimensions=[
                CnvStratificationDimension.COVERAGE,
                CnvStratificationDimension.TUMOR_FRACTION,
                CnvStratificationDimension.EVENT_CLASS,
                CnvStratificationDimension.EVENT_SIZE,
            ],
            secondary_dimensions=[
                CnvStratificationDimension.CALLER,
                CnvStratificationDimension.GENOME_BUILD,
                CnvStratificationDimension.DATA_BASIS,
                CnvStratificationDimension.BIN_SIZE,
            ],
            coverage_cutpoints_x=coverage_cutpoints,
            tumor_fraction_cutpoints=tumor_cutpoints,
            event_size_cutpoints_bp=size_cutpoints,
        ),
        event_classes=[
            EventType.CHROMOSOME_GAIN,
            EventType.CHROMOSOME_LOSS,
            EventType.DELETION,
            EventType.DUPLICATION,
        ],
        qdnaseq_bin_sizes_kbp=(
            design.qdnaseq_bin_sizes_kbp if "qdnaseq_ace" in caller_ids else []
        ),
        matching_thresholds=BenchmarkThresholds(
            minimum_reciprocal_overlap=0.5,
            maximum_breakpoint_distance_bp=500,
            copy_number_tolerance=0.5,
        ),
        acceptance=_default_acceptance_questions(design.replicates),
        rationale=(
            "Analytical Validation Program v1 prospective engineering candidate matrix. "
            "Thresholds are locked before outcome review and are not clinical cutoffs."
        ),
    )


def build_cnv_validation_v1_registration(
    design: CnvValidationV1Design,
    caller_locks: list[CnvCallerLock],
    *,
    registration_id: str,
    registered_at,
    code_sha256: str,
    software_version: str,
) -> CnvValidationRegistration:
    matrix = build_cnv_validation_v1_matrix(design, caller_locks)
    cohort = materialize_synthetic_validation_v1_cohort(design)
    return preregister_cnv_validation(
        matrix,
        cohort,
        registration_id=registration_id,
        registered_at=registered_at,
        code_sha256=code_sha256,
        software_version=software_version,
    )


def dilution_policy_from_validation_v1(
    design: CnvValidationV1Design,
    *,
    expected_samtools_version: str,
) -> DilutionPolicy:
    return DilutionPolicy(
        profile_id=f"{design.design_id}-dilution",
        status="technical_defaults_only",
        expected_samtools_version=expected_samtools_version,
        tumor_fractions=list(design.tumor_fractions),
        include_normal_only_control=True,
        replicates=design.replicates,
        seed=design.seed,
        total_read_target=None,
        observed_fraction_tolerance=0.02,
        note=(
            "Analytical Validation Program v1 in-silico tumour-fraction projection. "
            "This is not a wet-lab or clinical LoD."
        ),
    )


def summarize_cnv_validation_v1(
    design: CnvValidationV1Design,
    caller_locks: list[CnvCallerLock],
) -> dict[str, int | list[int] | list[float] | list[str]]:
    matrix = build_cnv_validation_v1_matrix(design, caller_locks)
    cells = expand_cnv_validation_v1_cells(design)
    caller_ids = [caller.caller_id for caller in caller_locks]
    lanes_per_cell = sum(
        len(design.qdnaseq_bin_sizes_kbp) if caller_id == "qdnaseq_ace" else 1
        for caller_id in caller_ids
    )
    return {
        "positive_truth_cells": len(cells),
        "caller_count": len(caller_locks),
        "callers": caller_ids,
        "coverage_levels_x": list(design.coverage_levels_x),
        "tumor_fractions": list(design.tumor_fractions),
        "subchromosomal_sizes_bp": list(design.subchromosomal_sizes_bp),
        "qdnaseq_bin_sizes_kbp": list(matrix.qdnaseq_bin_sizes_kbp),
        "planned_caller_resolution_lanes": len(cells) * lanes_per_cell,
    }
