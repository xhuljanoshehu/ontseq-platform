from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from itertools import combinations, product
from typing import Literal

from pydantic import Field, model_validator

from .benchmark import benchmark_case
from .cnv_validation_contracts import (
    CnvAcceptanceMetric,
    CnvAcceptanceQuestion,
    CnvDataBasis,
    CnvRepeatKind,
)
from .cnv_validation_evidence import (
    CnvContributionStatus,
    CnvEvidenceRecordKind,
    CnvFullEvidenceRecord,
    CnvNormalizedEventRecord,
    CnvRunOutcomeState,
    CnvValidationEvidenceManifest,
    canonical_evidence_sha256,
)
from .cnv_validation_registration import CnvValidationRegistration
from .models import (
    BenchmarkCase,
    BenchmarkKind,
    BenchmarkThresholds,
    EventType,
    GenomeBuild,
    GenomicEvent,
    StrictModel,
)


class CnvEvaluationLaneKey(StrictModel):
    specimen_id: str
    biological_specimen_id: str
    caller_id: str
    caller_version: str
    adapter_policy_sha256: str
    execution_identity_sha256: str
    genome_build: GenomeBuild
    data_basis: CnvDataBasis
    reference_id: str
    reference_sha256: str
    input_sha256: str
    coverage_x: float | None = None
    coverage_definition: str | None = None
    tumor_fraction: float | None = None
    tumor_fraction_method: str | None = None
    tumor_fraction_timepoint: str | None = None
    bin_size_kbp: int | None = None
    repeat_kind: CnvRepeatKind
    replicate_id: str
    repeat_group_id: str | None = None

    @property
    def lane_id(self) -> str:
        return "lane-" + canonical_evidence_sha256(self)[:24]


class CnvMatchedEventAssignment(StrictModel):
    truth_event_id: str
    normalized_event_id: str
    score: float = Field(ge=0, le=1)
    reciprocal_overlap: float | None = Field(default=None, ge=0, le=1)


class CnvLaneEventAssignment(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    lane_id: str
    lane: CnvEvaluationLaneKey
    run_summary_full_evidence_id: str
    run_outcome: CnvRunOutcomeState
    minimum_reciprocal_overlap: float = Field(ge=0, le=1)
    copy_number_tolerance: float | None = Field(default=None, ge=0)
    true_positive: int | None = Field(default=None, ge=0)
    false_positive: int | None = Field(default=None, ge=0)
    false_negative: int | None = Field(default=None, ge=0)
    matches: list[CnvMatchedEventAssignment] = Field(default_factory=list)
    unmatched_truth_event_ids: list[str] = Field(default_factory=list)
    unmatched_normalized_event_ids: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True


_LANE_FIELDS = (
    "specimen_id",
    "biological_specimen_id",
    "caller_id",
    "caller_version",
    "adapter_policy_sha256",
    "execution_identity_sha256",
    "genome_build",
    "data_basis",
    "reference_id",
    "reference_sha256",
    "input_sha256",
    "coverage_x",
    "coverage_definition",
    "tumor_fraction",
    "tumor_fraction_method",
    "tumor_fraction_timepoint",
    "bin_size_kbp",
    "repeat_kind",
    "replicate_id",
    "repeat_group_id",
)


def _lane_key(record: CnvFullEvidenceRecord | CnvNormalizedEventRecord) -> CnvEvaluationLaneKey:
    return CnvEvaluationLaneKey.model_validate(
        {field: getattr(record, field) for field in _LANE_FIELDS}
    )


def _validate_registered_identity(
    registration: CnvValidationRegistration,
    record: CnvFullEvidenceRecord | CnvNormalizedEventRecord,
) -> None:
    caller_by_id = {item.caller_id: item for item in registration.matrix.callers}
    caller = caller_by_id.get(record.caller_id)
    if caller is None:
        raise ValueError(f"Evidence caller is not registered: {record.caller_id}")
    expected_caller = {
        "caller_version": caller.caller_version,
        "adapter_policy_sha256": caller.adapter_policy_sha256,
        "execution_identity_sha256": caller.execution_identity_sha256,
    }
    for field_name, expected in expected_caller.items():
        if getattr(record, field_name) != expected:
            raise ValueError(f"Evidence caller identity mismatch for {field_name}")

    specimen_by_id = {item.specimen_id: item for item in registration.cohort.specimens}
    specimen = specimen_by_id.get(record.specimen_id)
    if specimen is None:
        raise ValueError(f"Evidence specimen is not registered: {record.specimen_id}")
    expected_specimen: dict[str, object] = {
        "biological_specimen_id": specimen.biological_specimen_id,
        "genome_build": specimen.genome_build,
        "data_basis": specimen.data_basis,
        "reference_id": specimen.reference_id,
        "reference_sha256": specimen.reference_sha256,
        "input_sha256": specimen.input_sha256,
        "coverage_x": specimen.coverage_x,
        "coverage_definition": specimen.coverage_definition,
        "tumor_fraction": specimen.tumor_fraction,
        "tumor_fraction_method": specimen.tumor_fraction_method,
        "tumor_fraction_timepoint": specimen.tumor_fraction_timepoint,
        "repeat_kind": specimen.repeat_kind,
        "repeat_group_id": specimen.repeat_group_id,
    }
    for field_name, expected_value in expected_specimen.items():
        if getattr(record, field_name) != expected_value:
            raise ValueError(f"Evidence specimen identity mismatch for {field_name}")

    if (
        record.caller_id == "qdnaseq_ace"
        and record.bin_size_kbp not in registration.matrix.qdnaseq_bin_sizes_kbp
    ):
        raise ValueError("QDNAseq+ACE evidence uses an unregistered bin size")


def _query_event(event: CnvNormalizedEventRecord) -> GenomicEvent:
    return GenomicEvent(
        event_id=event.normalized_event_id,
        event_type=event.event_type,
        primary=event.primary,
        copy_number=event.normalized_copy_number,
    )


def assign_cnv_validation_lanes(
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
) -> list[CnvLaneEventAssignment]:
    """Match primary normalized CNVs to truth once per locked analytical lane."""
    if evidence.registration_sha256 != registration.lock_sha256:
        raise ValueError("Evidence manifest is bound to a different validation registration")

    for record in evidence.full_evidence:
        _validate_registered_identity(registration, record)
    for event in evidence.normalized_events:
        _validate_registered_identity(registration, event)

    records_by_lane: dict[str, list[CnvFullEvidenceRecord]] = defaultdict(list)
    lane_keys: dict[str, CnvEvaluationLaneKey] = {}
    for record in evidence.full_evidence:
        lane = _lane_key(record)
        lane_keys[lane.lane_id] = lane
        records_by_lane[lane.lane_id].append(record)

    normalized_by_lane: dict[str, list[CnvNormalizedEventRecord]] = defaultdict(list)
    for event in evidence.normalized_events:
        lane = _lane_key(event)
        if lane.lane_id not in records_by_lane:
            raise ValueError("Normalized CNV event has no retained full-evidence lane")
        normalized_by_lane[lane.lane_id].append(event)

    specimen_by_id = {item.specimen_id: item for item in registration.cohort.specimens}
    thresholds = registration.matrix.matching_thresholds
    assignments: list[CnvLaneEventAssignment] = []

    for lane_id in sorted(records_by_lane):
        lane = lane_keys[lane_id]
        records = records_by_lane[lane_id]
        summaries = [
            record for record in records if record.record_kind == CnvEvidenceRecordKind.RUN_SUMMARY
        ]
        if len(summaries) != 1:
            raise ValueError("Each analytical CNV lane requires exactly one run-summary record")
        summary = summaries[0]
        if any(record.run_outcome != summary.run_outcome for record in records):
            raise ValueError("CNV lane contains inconsistent run-outcome states")

        if summary.run_outcome != CnvRunOutcomeState.OBSERVED:
            assignments.append(
                CnvLaneEventAssignment(
                    lane_id=lane_id,
                    lane=lane,
                    run_summary_full_evidence_id=summary.record_id,
                    run_outcome=summary.run_outcome,
                    minimum_reciprocal_overlap=thresholds.minimum_reciprocal_overlap,
                    copy_number_tolerance=thresholds.copy_number_tolerance,
                )
            )
            continue

        primary_events = [
            event
            for event in normalized_by_lane.get(lane_id, [])
            if event.contribution_status == CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS
        ]
        specimen = specimen_by_id[lane.specimen_id]
        case = BenchmarkCase(
            case_id="cnv-" + lane_id.removeprefix("lane-")[:24],
            kind=BenchmarkKind.CNV,
            genome_build=lane.genome_build,
            truth_events=specimen.truth_events,
            query_events=[_query_event(event) for event in primary_events],
            thresholds=thresholds,
        )
        report = benchmark_case(case)

        matches = [
            CnvMatchedEventAssignment(
                truth_event_id=match.truth_event_id,
                normalized_event_id=match.query_event_id,
                score=match.score,
                reciprocal_overlap=match.reciprocal_overlap,
            )
            for match in report.matches
        ]
        assignments.append(
            CnvLaneEventAssignment(
                lane_id=lane_id,
                lane=lane,
                run_summary_full_evidence_id=summary.record_id,
                run_outcome=summary.run_outcome,
                minimum_reciprocal_overlap=thresholds.minimum_reciprocal_overlap,
                copy_number_tolerance=thresholds.copy_number_tolerance,
                true_positive=report.metrics.true_positive,
                false_positive=report.metrics.false_positive,
                false_negative=report.metrics.false_negative,
                matches=matches,
                unmatched_truth_event_ids=report.unmatched_truth_event_ids,
                unmatched_normalized_event_ids=report.unmatched_query_event_ids,
            )
        )

    return assignments


_MISSING_SPECIFICITY_ASSESSMENT_REASON = (
    "Specificity requires an explicit negative-unit assessment for every observed lane."
)


class CnvMetricEvaluationState(StrEnum):
    EVALUABLE = "EVALUABLE"
    NOT_EVALUABLE = "NOT_EVALUABLE"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class CnvMetricResult(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    metric: CnvAcceptanceMetric
    stratum_key: dict[str, str | int | float | bool] = Field(min_length=1)
    state: CnvMetricEvaluationState
    value: float | None = None
    evaluable_denominator: int = Field(ge=0)
    component_counts: dict[str, int] = Field(default_factory=dict)
    summary_values: dict[str, float] = Field(default_factory=dict)
    registration_sha256: str
    evidence_manifest_sha256: str
    truth_event_ids: list[str] = Field(default_factory=list)
    normalized_event_ids: list[str] = Field(default_factory=list)
    full_evidence_ids: list[str] = Field(default_factory=list)
    lane_ids: list[str] = Field(default_factory=list)
    reason: str | None = None
    research_only: Literal[True] = True

    @property
    def metric_id(self) -> str:
        """Content-addressed metric identity including value, denominator and provenance."""
        return "metric-" + canonical_evidence_sha256(self)[:24]

    @model_validator(mode="after")
    def coherent_metric(self) -> CnvMetricResult:
        if self.state == CnvMetricEvaluationState.EVALUABLE:
            if self.value is None or not math.isfinite(self.value):
                raise ValueError("Evaluable CNV metric requires a finite value")
        elif self.value is not None:
            raise ValueError("Non-evaluable CNV metric cannot carry a value")
        if any(value < 0 for value in self.component_counts.values()):
            raise ValueError("CNV metric component counts cannot be negative")
        if any(not math.isfinite(value) for value in self.summary_values.values()):
            raise ValueError("CNV metric summary values must be finite")
        return self


class CnvNegativeUnitAssessment(StrictModel):
    """Explicit assessed negative universe for specificity; never inferred from event FPs."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    assessment_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{1,127}$")
    registration_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    evidence_manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    lane_id: str = Field(min_length=1)
    specimen_id: str = Field(min_length=1)
    universe_id: str = Field(min_length=1)
    universe_resource_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    assessability_mask_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    assessed_units: int = Field(ge=1)
    false_positive_units: int = Field(ge=0)
    full_evidence_ids: list[str] = Field(min_length=1)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_negative_unit_assessment(self) -> CnvNegativeUnitAssessment:
        if self.false_positive_units > self.assessed_units:
            raise ValueError("False-positive negative units cannot exceed assessed units")
        if len(self.full_evidence_ids) != len(set(self.full_evidence_ids)):
            raise ValueError("Negative-unit assessment evidence IDs must be unique")
        return self


@dataclass
class _RunStateAccumulator:
    executed_lane_ids: set[str] = field(default_factory=set)
    observed_lane_ids: set[str] = field(default_factory=set)
    no_call_lane_ids: set[str] = field(default_factory=set)
    failed_lane_ids: set[str] = field(default_factory=set)
    not_assessable_lane_ids: set[str] = field(default_factory=set)
    full_evidence_ids: set[str] = field(default_factory=set)


@dataclass
class _EventMetricAccumulator:
    tp_truth: int = 0
    false_negative: int = 0
    tp_query: int = 0
    false_positive: int = 0
    observed_lane_ids: set[str] = field(default_factory=set)
    truth_event_ids: set[str] = field(default_factory=set)
    normalized_event_ids: set[str] = field(default_factory=set)
    full_evidence_ids: set[str] = field(default_factory=set)
    lane_ids: set[str] = field(default_factory=set)


def _format_cutpoint(value: float | int) -> str:
    numeric = float(value)
    if numeric.is_integer():
        return str(int(numeric))
    return format(numeric, ".15g")


def cnv_numeric_band(value: float, cutpoints: Sequence[float | int]) -> str:
    """Return a deterministic [lower, upper) label using preregistered cutpoints."""
    if not math.isfinite(value):
        raise ValueError("Stratification value must be finite")
    points = [float(item) for item in cutpoints]
    if not points or points != sorted(set(points)):
        raise ValueError("Stratification cutpoints must be unique and sorted")
    if value < points[0]:
        return f"<{_format_cutpoint(points[0])}"
    for lower, upper in zip(points[:-1], points[1:], strict=True):
        if lower <= value < upper:
            return f"[{_format_cutpoint(lower)},{_format_cutpoint(upper)})"
    return f">={_format_cutpoint(points[-1])}"


def _registered_numeric_bands(cutpoints: Sequence[float | int]) -> list[str]:
    points = [float(item) for item in cutpoints]
    if not points:
        return ["all"]
    bands = [f"<{_format_cutpoint(points[0])}"]
    bands.extend(
        f"[{_format_cutpoint(lower)},{_format_cutpoint(upper)})"
        for lower, upper in zip(points[:-1], points[1:], strict=True)
    )
    bands.append(f">={_format_cutpoint(points[-1])}")
    return bands


def _event_size_band(
    event_type: EventType,
    start: int,
    end: int,
    cutpoints: Sequence[int],
) -> str:
    if event_type in {EventType.CHROMOSOME_GAIN, EventType.CHROMOSOME_LOSS}:
        return "whole_chromosome"
    return cnv_numeric_band(float(end - start), cutpoints)


def _registered_event_strata(
    registration: CnvValidationRegistration,
) -> list[dict[str, str | int | float | bool]]:
    matrix = registration.matrix
    coverage_bands = _registered_numeric_bands(matrix.stratification.coverage_cutpoints_x)
    tumor_bands = _registered_numeric_bands(matrix.stratification.tumor_fraction_cutpoints)
    subchromosomal_size_bands = _registered_numeric_bands(
        matrix.stratification.event_size_cutpoints_bp
    )
    strata: list[dict[str, str | int | float | bool]] = []
    for caller, build, data_basis, coverage, tumor, event_class in product(
        matrix.callers,
        matrix.genome_builds,
        matrix.data_bases,
        coverage_bands,
        tumor_bands,
        matrix.event_classes,
    ):
        bin_sizes: list[int | str]
        if caller.caller_id == "qdnaseq_ace":
            bin_sizes = list(matrix.qdnaseq_bin_sizes_kbp)
        else:
            bin_sizes = ["not_applicable"]
        size_bands = (
            ["whole_chromosome"]
            if event_class in {EventType.CHROMOSOME_GAIN, EventType.CHROMOSOME_LOSS}
            else subchromosomal_size_bands
        )
        for bin_size, event_size in product(bin_sizes, size_bands):
            strata.append(
                {
                    "caller": caller.caller_id,
                    "genome_build": build.value,
                    "data_basis": data_basis.value,
                    "coverage": coverage,
                    "tumor_fraction": tumor,
                    "event_class": event_class.value,
                    "event_size": event_size,
                    "bin_size": bin_size,
                }
            )
    return strata


def _lane_technical_stratum(
    lane: CnvEvaluationLaneKey,
    registration: CnvValidationRegistration,
) -> dict[str, str | int | float | bool]:
    if lane.coverage_x is None or lane.tumor_fraction is None:
        raise ValueError("Prospective event metrics require registered continuous lane strata")
    bin_size: int | str = lane.bin_size_kbp if lane.bin_size_kbp is not None else "not_applicable"
    return {
        "caller": lane.caller_id,
        "genome_build": lane.genome_build.value,
        "data_basis": lane.data_basis.value,
        "coverage": cnv_numeric_band(
            lane.coverage_x,
            registration.matrix.stratification.coverage_cutpoints_x,
        ),
        "tumor_fraction": cnv_numeric_band(
            lane.tumor_fraction,
            registration.matrix.stratification.tumor_fraction_cutpoints,
        ),
        "bin_size": bin_size,
    }


def _event_stratum(
    lane: CnvEvaluationLaneKey,
    event_type: EventType,
    start: int,
    end: int,
    registration: CnvValidationRegistration,
) -> dict[str, str | int | float | bool]:
    key = _lane_technical_stratum(lane, registration)
    key["event_class"] = event_type.value
    key["event_size"] = _event_size_band(
        event_type,
        start,
        end,
        registration.matrix.stratification.event_size_cutpoints_bp,
    )
    return key


def _stratum_address(
    value: dict[str, str | int | float | bool],
) -> tuple[tuple[str, str | int | float | bool], ...]:
    return tuple(sorted(value.items()))


def _metric_from_accumulator(
    *,
    metric: CnvAcceptanceMetric,
    stratum_key: dict[str, str | int | float | bool],
    accumulator: _EventMetricAccumulator,
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
) -> CnvMetricResult:
    sensitivity_denominator = accumulator.tp_truth + accumulator.false_negative
    precision_denominator = accumulator.tp_query + accumulator.false_positive
    observed_lanes = len(accumulator.observed_lane_ids)

    state = CnvMetricEvaluationState.EVALUABLE
    value: float | None
    denominator: int
    reason: str | None = None

    if metric == CnvAcceptanceMetric.SENSITIVITY:
        denominator = sensitivity_denominator
        value = accumulator.tp_truth / denominator if denominator else None
        if denominator == 0:
            state = CnvMetricEvaluationState.NOT_EVALUABLE
            reason = "No evaluable truth events in the registered stratum."
    elif metric == CnvAcceptanceMetric.PRECISION:
        denominator = precision_denominator
        value = accumulator.tp_query / denominator if denominator else None
        if denominator == 0:
            state = CnvMetricEvaluationState.NOT_EVALUABLE
            reason = "No evaluable query events in the registered stratum."
    elif metric == CnvAcceptanceMetric.F1:
        denominator = min(sensitivity_denominator, precision_denominator)
        if sensitivity_denominator == 0 or precision_denominator == 0:
            value = None
            state = CnvMetricEvaluationState.NOT_EVALUABLE
            reason = "F1 requires both evaluable truth-side and query-side denominators."
        else:
            sensitivity = accumulator.tp_truth / sensitivity_denominator
            precision = accumulator.tp_query / precision_denominator
            value = (
                0.0
                if sensitivity + precision == 0
                else 2 * sensitivity * precision / (sensitivity + precision)
            )
    elif metric == CnvAcceptanceMetric.FALSE_POSITIVE_BURDEN:
        denominator = observed_lanes
        value = accumulator.false_positive / denominator if denominator else None
        if denominator == 0:
            state = CnvMetricEvaluationState.NOT_EVALUABLE
            reason = "No observed caller lanes in the registered technical stratum."
    else:
        raise ValueError(f"Unsupported event metric: {metric.value}")

    return CnvMetricResult(
        metric=metric,
        stratum_key=stratum_key,
        state=state,
        value=value,
        evaluable_denominator=denominator,
        component_counts={
            "tp_truth": accumulator.tp_truth,
            "fn": accumulator.false_negative,
            "tp_query": accumulator.tp_query,
            "fp": accumulator.false_positive,
            "observed_lanes": observed_lanes,
        },
        registration_sha256=registration.lock_sha256,
        evidence_manifest_sha256=evidence.manifest_sha256,
        truth_event_ids=sorted(accumulator.truth_event_ids),
        normalized_event_ids=sorted(accumulator.normalized_event_ids),
        full_evidence_ids=sorted(accumulator.full_evidence_ids),
        lane_ids=sorted(accumulator.lane_ids),
        reason=reason,
    )


def aggregate_cnv_event_metrics(
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
    *,
    assignments: Sequence[CnvLaneEventAssignment] | None = None,
) -> list[CnvMetricResult]:
    """Aggregate locked lane assignments without re-matching inside strata."""
    assignments = (
        assign_cnv_validation_lanes(registration, evidence)
        if assignments is None
        else list(assignments)
    )
    normalized_by_id = {item.normalized_event_id: item for item in evidence.normalized_events}
    specimen_by_id = {item.specimen_id: item for item in registration.cohort.specimens}

    registered_keys = _registered_event_strata(registration)
    accumulators = {_stratum_address(key): _EventMetricAccumulator() for key in registered_keys}
    keys_by_address = {_stratum_address(key): key for key in registered_keys}
    overall_key: dict[str, str | int | float | bool] = {"scope": "overall"}
    overall = _EventMetricAccumulator()

    for assignment in assignments:
        if assignment.run_outcome != CnvRunOutcomeState.OBSERVED:
            continue

        lane = assignment.lane
        technical = _lane_technical_stratum(lane, registration)
        overall.observed_lane_ids.add(assignment.lane_id)
        overall.lane_ids.add(assignment.lane_id)
        overall.full_evidence_ids.add(assignment.run_summary_full_evidence_id)

        for address, key in keys_by_address.items():
            if all(key.get(name) == value for name, value in technical.items()):
                accumulators[address].observed_lane_ids.add(assignment.lane_id)
                accumulators[address].lane_ids.add(assignment.lane_id)
                accumulators[address].full_evidence_ids.add(assignment.run_summary_full_evidence_id)

        specimen = specimen_by_id[lane.specimen_id]
        truth_by_id = {item.event_id: item for item in specimen.truth_events}

        for match in assignment.matches:
            truth = truth_by_id[match.truth_event_id]
            query = normalized_by_id[match.normalized_event_id]
            truth_key = _event_stratum(
                lane,
                truth.event_type,
                truth.primary.start,
                truth.primary.end,
                registration,
            )
            query_key = _event_stratum(
                lane,
                query.event_type,
                query.primary.start,
                query.primary.end,
                registration,
            )

            truth_acc = accumulators[_stratum_address(truth_key)]
            truth_acc.tp_truth += 1
            truth_acc.truth_event_ids.add(truth.event_id)
            truth_acc.normalized_event_ids.add(query.normalized_event_id)
            truth_acc.full_evidence_ids.add(assignment.run_summary_full_evidence_id)
            truth_acc.full_evidence_ids.update(query.source_full_evidence_ids)
            truth_acc.lane_ids.add(assignment.lane_id)

            query_acc = accumulators[_stratum_address(query_key)]
            query_acc.tp_query += 1
            query_acc.truth_event_ids.add(truth.event_id)
            query_acc.normalized_event_ids.add(query.normalized_event_id)
            query_acc.full_evidence_ids.add(assignment.run_summary_full_evidence_id)
            query_acc.full_evidence_ids.update(query.source_full_evidence_ids)
            query_acc.lane_ids.add(assignment.lane_id)

            overall.tp_truth += 1
            overall.tp_query += 1
            overall.truth_event_ids.add(truth.event_id)
            overall.normalized_event_ids.add(query.normalized_event_id)
            overall.full_evidence_ids.update(query.source_full_evidence_ids)

        for truth_id in assignment.unmatched_truth_event_ids:
            truth = truth_by_id[truth_id]
            key = _event_stratum(
                lane,
                truth.event_type,
                truth.primary.start,
                truth.primary.end,
                registration,
            )
            accumulator = accumulators[_stratum_address(key)]
            accumulator.false_negative += 1
            accumulator.truth_event_ids.add(truth.event_id)
            accumulator.full_evidence_ids.add(assignment.run_summary_full_evidence_id)
            accumulator.lane_ids.add(assignment.lane_id)
            overall.false_negative += 1
            overall.truth_event_ids.add(truth.event_id)

        for normalized_id in assignment.unmatched_normalized_event_ids:
            query = normalized_by_id[normalized_id]
            key = _event_stratum(
                lane,
                query.event_type,
                query.primary.start,
                query.primary.end,
                registration,
            )
            accumulator = accumulators[_stratum_address(key)]
            accumulator.false_positive += 1
            accumulator.normalized_event_ids.add(query.normalized_event_id)
            accumulator.full_evidence_ids.add(assignment.run_summary_full_evidence_id)
            accumulator.full_evidence_ids.update(query.source_full_evidence_ids)
            accumulator.lane_ids.add(assignment.lane_id)
            overall.false_positive += 1
            overall.normalized_event_ids.add(query.normalized_event_id)
            overall.full_evidence_ids.update(query.source_full_evidence_ids)

    metrics_to_emit = (
        CnvAcceptanceMetric.SENSITIVITY,
        CnvAcceptanceMetric.PRECISION,
        CnvAcceptanceMetric.F1,
        CnvAcceptanceMetric.FALSE_POSITIVE_BURDEN,
    )
    results = [
        _metric_from_accumulator(
            metric=metric,
            stratum_key=overall_key,
            accumulator=overall,
            registration=registration,
            evidence=evidence,
        )
        for metric in metrics_to_emit
    ]
    for address in sorted(keys_by_address, key=str):
        key = keys_by_address[address]
        accumulator = accumulators[address]
        results.extend(
            _metric_from_accumulator(
                metric=metric,
                stratum_key=key,
                accumulator=accumulator,
                registration=registration,
                evidence=evidence,
            )
            for metric in metrics_to_emit
        )
    return results


def _registered_technical_strata(
    registration: CnvValidationRegistration,
) -> list[dict[str, str | int | float | bool]]:
    matrix = registration.matrix
    coverage_bands = _registered_numeric_bands(matrix.stratification.coverage_cutpoints_x)
    tumor_bands = _registered_numeric_bands(matrix.stratification.tumor_fraction_cutpoints)
    strata: list[dict[str, str | int | float | bool]] = []
    for caller, build, data_basis, coverage, tumor in product(
        matrix.callers,
        matrix.genome_builds,
        matrix.data_bases,
        coverage_bands,
        tumor_bands,
    ):
        bin_sizes: list[int | str]
        if caller.caller_id == "qdnaseq_ace":
            bin_sizes = list(matrix.qdnaseq_bin_sizes_kbp)
        else:
            bin_sizes = ["not_applicable"]
        for bin_size in bin_sizes:
            strata.append(
                {
                    "caller": caller.caller_id,
                    "genome_build": build.value,
                    "data_basis": data_basis.value,
                    "coverage": coverage,
                    "tumor_fraction": tumor,
                    "bin_size": bin_size,
                }
            )
    return strata


def _run_state_metric(
    *,
    metric: CnvAcceptanceMetric,
    stratum_key: dict[str, str | int | float | bool],
    accumulator: _RunStateAccumulator,
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
) -> CnvMetricResult:
    denominator = len(accumulator.executed_lane_ids)
    counts = {
        "executed_lanes": denominator,
        "observed_lanes": len(accumulator.observed_lane_ids),
        "no_call_lanes": len(accumulator.no_call_lane_ids),
        "failed_lanes": len(accumulator.failed_lane_ids),
        "not_assessable_lanes": len(accumulator.not_assessable_lane_ids),
    }
    numerator_by_metric = {
        CnvAcceptanceMetric.NO_CALL_RATE: counts["no_call_lanes"],
        CnvAcceptanceMetric.TECHNICAL_FAILURE_RATE: counts["failed_lanes"],
        CnvAcceptanceMetric.NOT_ASSESSABLE_RATE: counts["not_assessable_lanes"],
    }
    numerator = numerator_by_metric[metric]
    if denominator == 0:
        return CnvMetricResult(
            metric=metric,
            stratum_key=stratum_key,
            state=CnvMetricEvaluationState.NOT_EVALUABLE,
            evaluable_denominator=0,
            component_counts=counts,
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            full_evidence_ids=sorted(accumulator.full_evidence_ids),
            lane_ids=sorted(accumulator.executed_lane_ids),
            reason="No executed caller lanes in the registered technical stratum.",
        )
    return CnvMetricResult(
        metric=metric,
        stratum_key=stratum_key,
        state=CnvMetricEvaluationState.EVALUABLE,
        value=numerator / denominator,
        evaluable_denominator=denominator,
        component_counts=counts,
        registration_sha256=registration.lock_sha256,
        evidence_manifest_sha256=evidence.manifest_sha256,
        full_evidence_ids=sorted(accumulator.full_evidence_ids),
        lane_ids=sorted(accumulator.executed_lane_ids),
    )


def _validate_negative_assessments(
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
    assignments: list[CnvLaneEventAssignment],
    negative_assessments: Sequence[CnvNegativeUnitAssessment],
) -> dict[str, CnvNegativeUnitAssessment]:
    assignment_by_lane = {item.lane_id: item for item in assignments}
    specimen_by_id = {item.specimen_id: item for item in registration.cohort.specimens}
    full_by_id = {item.record_id: item for item in evidence.full_evidence}
    assessment_ids = [item.assessment_id for item in negative_assessments]
    if len(assessment_ids) != len(set(assessment_ids)):
        raise ValueError("Negative-unit assessment IDs must be unique")

    by_lane: dict[str, CnvNegativeUnitAssessment] = {}
    for assessment in negative_assessments:
        if assessment.registration_sha256 != registration.lock_sha256:
            raise ValueError("Negative-unit assessment registration lock mismatch")
        if assessment.evidence_manifest_sha256 != evidence.manifest_sha256:
            raise ValueError("Negative-unit assessment evidence-manifest lock mismatch")
        assignment = assignment_by_lane.get(assessment.lane_id)
        if assignment is None:
            raise ValueError("Negative-unit assessment references an unknown lane")
        if assessment.lane_id in by_lane:
            raise ValueError("At most one negative-unit assessment is allowed per lane")
        if assignment.run_outcome != CnvRunOutcomeState.OBSERVED:
            raise ValueError("Specificity assessment requires an observed caller lane")
        if assessment.specimen_id != assignment.lane.specimen_id:
            raise ValueError("Negative-unit assessment specimen does not match its lane")

        specimen = specimen_by_id[assessment.specimen_id]
        universe = specimen.negative_universe
        if universe is None:
            raise ValueError("Negative-unit assessment requires a registered negative universe")
        expected = {
            "universe_id": universe.universe_id,
            "universe_resource_sha256": universe.resource_sha256,
            "assessability_mask_sha256": specimen.assessability_mask.resource_sha256,
            "assessed_units": universe.assessable_units,
        }
        for field_name, expected_value in expected.items():
            if getattr(assessment, field_name) != expected_value:
                raise ValueError(f"Negative-unit assessment does not match registered {field_name}")

        if assignment.run_summary_full_evidence_id not in assessment.full_evidence_ids:
            raise ValueError("Negative-unit assessment must retain its run-summary evidence ID")
        for full_id in assessment.full_evidence_ids:
            record = full_by_id.get(full_id)
            if record is None:
                raise ValueError("Negative-unit assessment references missing full evidence")
            if _lane_key(record).lane_id != assessment.lane_id:
                raise ValueError("Negative-unit assessment evidence must belong to the same lane")
        by_lane[assessment.lane_id] = assessment
    return by_lane


def _specificity_metric(
    *,
    stratum_key: dict[str, str | int | float | bool],
    observed_assignments: Sequence[CnvLaneEventAssignment],
    assessments_by_lane: dict[str, CnvNegativeUnitAssessment],
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
) -> CnvMetricResult:
    lane_ids = {item.lane_id for item in observed_assignments}
    if not lane_ids:
        return CnvMetricResult(
            metric=CnvAcceptanceMetric.SPECIFICITY,
            stratum_key=stratum_key,
            state=CnvMetricEvaluationState.NOT_EVALUABLE,
            evaluable_denominator=0,
            component_counts={
                "observed_lanes": 0,
                "assessed_lanes": 0,
                "assessed_units": 0,
                "true_negative_units": 0,
                "false_positive_units": 0,
            },
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            reason="No observed caller lanes in the registered technical stratum.",
        )

    missing = lane_ids - set(assessments_by_lane)
    if missing:
        full_ids = {item.run_summary_full_evidence_id for item in observed_assignments}
        return CnvMetricResult(
            metric=CnvAcceptanceMetric.SPECIFICITY,
            stratum_key=stratum_key,
            state=CnvMetricEvaluationState.NOT_EVALUABLE,
            evaluable_denominator=0,
            component_counts={
                "observed_lanes": len(lane_ids),
                "assessed_lanes": len(lane_ids) - len(missing),
                "assessed_units": 0,
                "true_negative_units": 0,
                "false_positive_units": 0,
            },
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            full_evidence_ids=sorted(full_ids),
            lane_ids=sorted(lane_ids),
            reason=_MISSING_SPECIFICITY_ASSESSMENT_REASON,
        )

    assessments = [assessments_by_lane[lane_id] for lane_id in sorted(lane_ids)]
    assessed_units = sum(item.assessed_units for item in assessments)
    false_positive_units = sum(item.false_positive_units for item in assessments)
    true_negative_units = assessed_units - false_positive_units
    full_ids = {full_id for item in assessments for full_id in item.full_evidence_ids}
    return CnvMetricResult(
        metric=CnvAcceptanceMetric.SPECIFICITY,
        stratum_key=stratum_key,
        state=CnvMetricEvaluationState.EVALUABLE,
        value=true_negative_units / assessed_units,
        evaluable_denominator=assessed_units,
        component_counts={
            "observed_lanes": len(lane_ids),
            "assessed_lanes": len(assessments),
            "assessed_units": assessed_units,
            "true_negative_units": true_negative_units,
            "false_positive_units": false_positive_units,
        },
        registration_sha256=registration.lock_sha256,
        evidence_manifest_sha256=evidence.manifest_sha256,
        full_evidence_ids=sorted(full_ids),
        lane_ids=sorted(lane_ids),
    )


def aggregate_cnv_run_state_metrics(
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
    *,
    negative_assessments: Sequence[CnvNegativeUnitAssessment] = (),
    assignments: Sequence[CnvLaneEventAssignment] | None = None,
) -> list[CnvMetricResult]:
    """Aggregate lane outcomes and specificity without inferring negative genomic units."""
    assignments = (
        assign_cnv_validation_lanes(registration, evidence)
        if assignments is None
        else list(assignments)
    )
    assessments_by_lane = _validate_negative_assessments(
        registration,
        evidence,
        assignments,
        negative_assessments,
    )

    registered_keys = _registered_technical_strata(registration)
    accumulators = {_stratum_address(key): _RunStateAccumulator() for key in registered_keys}
    keys_by_address = {_stratum_address(key): key for key in registered_keys}
    overall = _RunStateAccumulator()
    overall_key: dict[str, str | int | float | bool] = {"scope": "overall"}
    executed_outcomes = {
        CnvRunOutcomeState.OBSERVED,
        CnvRunOutcomeState.NO_CALL,
        CnvRunOutcomeState.FAILED,
        CnvRunOutcomeState.NOT_ASSESSABLE,
    }

    for assignment in assignments:
        if assignment.run_outcome not in executed_outcomes:
            continue
        technical = _lane_technical_stratum(assignment.lane, registration)
        address = _stratum_address(technical)
        accumulator = accumulators.get(address)
        if accumulator is None:
            raise ValueError("Executed lane falls outside the registered technical strata")
        for target in (overall, accumulator):
            target.executed_lane_ids.add(assignment.lane_id)
            target.full_evidence_ids.add(assignment.run_summary_full_evidence_id)
            if assignment.run_outcome == CnvRunOutcomeState.OBSERVED:
                target.observed_lane_ids.add(assignment.lane_id)
            elif assignment.run_outcome == CnvRunOutcomeState.NO_CALL:
                target.no_call_lane_ids.add(assignment.lane_id)
            elif assignment.run_outcome == CnvRunOutcomeState.FAILED:
                target.failed_lane_ids.add(assignment.lane_id)
            elif assignment.run_outcome == CnvRunOutcomeState.NOT_ASSESSABLE:
                target.not_assessable_lane_ids.add(assignment.lane_id)

    rate_metrics = (
        CnvAcceptanceMetric.NO_CALL_RATE,
        CnvAcceptanceMetric.TECHNICAL_FAILURE_RATE,
        CnvAcceptanceMetric.NOT_ASSESSABLE_RATE,
    )
    results = [
        _run_state_metric(
            metric=metric,
            stratum_key=overall_key,
            accumulator=overall,
            registration=registration,
            evidence=evidence,
        )
        for metric in rate_metrics
    ]

    observed_overall = [
        assignment
        for assignment in assignments
        if assignment.run_outcome == CnvRunOutcomeState.OBSERVED
    ]
    results.append(
        _specificity_metric(
            stratum_key=overall_key,
            observed_assignments=observed_overall,
            assessments_by_lane=assessments_by_lane,
            registration=registration,
            evidence=evidence,
        )
    )

    for address in sorted(keys_by_address, key=str):
        key = keys_by_address[address]
        accumulator = accumulators[address]
        results.extend(
            _run_state_metric(
                metric=metric,
                stratum_key=key,
                accumulator=accumulator,
                registration=registration,
                evidence=evidence,
            )
            for metric in rate_metrics
        )
        observed = [
            assignment
            for assignment in assignments
            if assignment.run_outcome == CnvRunOutcomeState.OBSERVED
            and _stratum_address(_lane_technical_stratum(assignment.lane, registration)) == address
        ]
        results.append(
            _specificity_metric(
                stratum_key=key,
                observed_assignments=observed,
                assessments_by_lane=assessments_by_lane,
                registration=registration,
                evidence=evidence,
            )
        )
    return results


@dataclass
class _QuantitativeAccumulator:
    absolute_errors: list[float] = field(default_factory=list)
    signed_errors: list[float] = field(default_factory=list)
    truth_event_ids: set[str] = field(default_factory=set)
    normalized_event_ids: set[str] = field(default_factory=set)
    full_evidence_ids: set[str] = field(default_factory=set)
    lane_ids: set[str] = field(default_factory=set)


def _quantitative_metric(
    *,
    metric: CnvAcceptanceMetric,
    stratum_key: dict[str, str | int | float | bool],
    accumulator: _QuantitativeAccumulator,
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
    denominator_label: str,
) -> CnvMetricResult:
    denominator = len(accumulator.absolute_errors)
    if denominator == 0:
        return CnvMetricResult(
            metric=metric,
            stratum_key=stratum_key,
            state=CnvMetricEvaluationState.NOT_EVALUABLE,
            evaluable_denominator=0,
            component_counts={denominator_label: 0},
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            truth_event_ids=sorted(accumulator.truth_event_ids),
            normalized_event_ids=sorted(accumulator.normalized_event_ids),
            full_evidence_ids=sorted(accumulator.full_evidence_ids),
            lane_ids=sorted(accumulator.lane_ids),
            reason="No eligible quantitative truth/caller pairs in this stratum.",
        )
    mean_absolute_error = sum(accumulator.absolute_errors) / denominator
    signed_bias = sum(accumulator.signed_errors) / denominator
    return CnvMetricResult(
        metric=metric,
        stratum_key=stratum_key,
        state=CnvMetricEvaluationState.EVALUABLE,
        value=mean_absolute_error,
        evaluable_denominator=denominator,
        component_counts={denominator_label: denominator},
        summary_values={
            "mean_absolute_error": mean_absolute_error,
            "signed_bias": signed_bias,
        },
        registration_sha256=registration.lock_sha256,
        evidence_manifest_sha256=evidence.manifest_sha256,
        truth_event_ids=sorted(accumulator.truth_event_ids),
        normalized_event_ids=sorted(accumulator.normalized_event_ids),
        full_evidence_ids=sorted(accumulator.full_evidence_ids),
        lane_ids=sorted(accumulator.lane_ids),
    )


def aggregate_cnv_quantitative_metrics(
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
    *,
    assignments: Sequence[CnvLaneEventAssignment] | None = None,
) -> list[CnvMetricResult]:
    """Aggregate quantitative errors only where explicit orthogonal truth exists."""
    assignments = (
        assign_cnv_validation_lanes(registration, evidence)
        if assignments is None
        else list(assignments)
    )
    normalized_by_id = {item.normalized_event_id: item for item in evidence.normalized_events}
    specimen_by_id = {item.specimen_id: item for item in registration.cohort.specimens}

    technical_keys = _registered_technical_strata(registration)
    key_by_address = {_stratum_address(key): key for key in technical_keys}
    copy_acc = {_stratum_address(key): _QuantitativeAccumulator() for key in technical_keys}
    cellularity_acc = {_stratum_address(key): _QuantitativeAccumulator() for key in technical_keys}
    ploidy_acc = {_stratum_address(key): _QuantitativeAccumulator() for key in technical_keys}
    overall_copy = _QuantitativeAccumulator()
    overall_cellularity = _QuantitativeAccumulator()
    overall_ploidy = _QuantitativeAccumulator()

    fit_records_by_lane: dict[str, list[CnvFullEvidenceRecord]] = defaultdict(list)
    for record in evidence.full_evidence:
        if record.record_kind == CnvEvidenceRecordKind.CALLER_FIT:
            fit_records_by_lane[_lane_key(record).lane_id].append(record)

    for assignment in assignments:
        if assignment.run_outcome != CnvRunOutcomeState.OBSERVED:
            continue
        lane = assignment.lane
        technical_address = _stratum_address(_lane_technical_stratum(lane, registration))
        if technical_address not in key_by_address:
            raise ValueError("Observed quantitative lane falls outside registered strata")
        specimen = specimen_by_id[lane.specimen_id]
        truth_by_id = {item.event_id: item for item in specimen.truth_events}

        for match in assignment.matches:
            truth = truth_by_id[match.truth_event_id]
            query = normalized_by_id[match.normalized_event_id]
            if truth.copy_number is None or query.normalized_copy_number is None:
                continue
            signed_error = query.normalized_copy_number - truth.copy_number
            for accumulator in (overall_copy, copy_acc[technical_address]):
                accumulator.absolute_errors.append(abs(signed_error))
                accumulator.signed_errors.append(signed_error)
                accumulator.truth_event_ids.add(truth.event_id)
                accumulator.normalized_event_ids.add(query.normalized_event_id)
                accumulator.full_evidence_ids.add(assignment.run_summary_full_evidence_id)
                accumulator.full_evidence_ids.update(query.source_full_evidence_ids)
                accumulator.lane_ids.add(assignment.lane_id)

        selected_fits = [
            record
            for record in fit_records_by_lane.get(assignment.lane_id, [])
            if record.selected_fit is True
        ]
        if len(selected_fits) > 1:
            raise ValueError("Observed CNV lane has multiple selected caller fits")
        if not selected_fits or specimen.quantitative_truth is None:
            continue
        selected = selected_fits[0]
        quantitative_truth = specimen.quantitative_truth

        if quantitative_truth.cellularity is not None and selected.cellularity is not None:
            signed_error = selected.cellularity - quantitative_truth.cellularity
            for accumulator in (overall_cellularity, cellularity_acc[technical_address]):
                accumulator.absolute_errors.append(abs(signed_error))
                accumulator.signed_errors.append(signed_error)
                accumulator.full_evidence_ids.update(
                    {assignment.run_summary_full_evidence_id, selected.record_id}
                )
                accumulator.lane_ids.add(assignment.lane_id)

        if quantitative_truth.ploidy is not None and selected.ploidy is not None:
            signed_error = selected.ploidy - quantitative_truth.ploidy
            for accumulator in (overall_ploidy, ploidy_acc[technical_address]):
                accumulator.absolute_errors.append(abs(signed_error))
                accumulator.signed_errors.append(signed_error)
                accumulator.full_evidence_ids.update(
                    {assignment.run_summary_full_evidence_id, selected.record_id}
                )
                accumulator.lane_ids.add(assignment.lane_id)

    overall_key: dict[str, str | int | float | bool] = {"scope": "overall"}
    results = [
        _quantitative_metric(
            metric=CnvAcceptanceMetric.COPY_NUMBER_ERROR,
            stratum_key=overall_key,
            accumulator=overall_copy,
            registration=registration,
            evidence=evidence,
            denominator_label="eligible_pairs",
        ),
        _quantitative_metric(
            metric=CnvAcceptanceMetric.CELLULARITY_ERROR,
            stratum_key=overall_key,
            accumulator=overall_cellularity,
            registration=registration,
            evidence=evidence,
            denominator_label="eligible_lanes",
        ),
        _quantitative_metric(
            metric=CnvAcceptanceMetric.PLOIDY_ERROR,
            stratum_key=overall_key,
            accumulator=overall_ploidy,
            registration=registration,
            evidence=evidence,
            denominator_label="eligible_lanes",
        ),
    ]
    for address in sorted(key_by_address, key=str):
        key = key_by_address[address]
        results.extend(
            (
                _quantitative_metric(
                    metric=CnvAcceptanceMetric.COPY_NUMBER_ERROR,
                    stratum_key=key,
                    accumulator=copy_acc[address],
                    registration=registration,
                    evidence=evidence,
                    denominator_label="eligible_pairs",
                ),
                _quantitative_metric(
                    metric=CnvAcceptanceMetric.CELLULARITY_ERROR,
                    stratum_key=key,
                    accumulator=cellularity_acc[address],
                    registration=registration,
                    evidence=evidence,
                    denominator_label="eligible_lanes",
                ),
                _quantitative_metric(
                    metric=CnvAcceptanceMetric.PLOIDY_ERROR,
                    stratum_key=key,
                    accumulator=ploidy_acc[address],
                    registration=registration,
                    evidence=evidence,
                    denominator_label="eligible_lanes",
                ),
            )
        )
    return results


class CnvReproducibilityPair(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    pair_id: str
    repeat_group_id: str
    repeat_kind: CnvRepeatKind
    biological_specimen_id: str
    caller_id: str
    genome_build: GenomeBuild
    data_basis: CnvDataBasis
    bin_size_kbp: int | None = None
    left_lane_id: str
    right_lane_id: str
    left_run_outcome: CnvRunOutcomeState
    right_run_outcome: CnvRunOutcomeState
    left_normalized_event_ids: list[str] = Field(default_factory=list)
    right_normalized_event_ids: list[str] = Field(default_factory=list)
    matched_normalized_event_pairs: list[list[str]] = Field(default_factory=list)
    concordance: float | None = Field(default=None, ge=0, le=1)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_reproducibility_pair(self) -> CnvReproducibilityPair:
        if any(len(pair) != 2 for pair in self.matched_normalized_event_pairs):
            raise ValueError("Reproducibility match pairs must contain exactly two event IDs")
        if (
            not self.left_normalized_event_ids
            and not self.right_normalized_event_ids
            and self.concordance is not None
        ):
            raise ValueError("Two empty repeat call sets cannot have event concordance")
        return self


def _primary_events_by_lane(
    evidence: CnvValidationEvidenceManifest,
) -> dict[str, list[CnvNormalizedEventRecord]]:
    events: dict[str, list[CnvNormalizedEventRecord]] = defaultdict(list)
    for event in evidence.normalized_events:
        if event.contribution_status == CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS:
            events[_lane_key(event).lane_id].append(event)
    for values in events.values():
        values.sort(key=lambda item: item.normalized_event_id)
    return events


def aggregate_cnv_reproducibility_metrics(
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
    *,
    assignments: Sequence[CnvLaneEventAssignment] | None = None,
) -> tuple[list[CnvMetricResult], list[CnvReproducibilityPair]]:
    """Compare repeated lanes symmetrically without treating two empty call sets as perfect."""
    assignments = (
        assign_cnv_validation_lanes(registration, evidence)
        if assignments is None
        else list(assignments)
    )
    events_by_lane = _primary_events_by_lane(evidence)
    grouped: dict[
        tuple[str, str, str, str, str, int | None],
        list[CnvLaneEventAssignment],
    ] = defaultdict(list)

    for assignment in assignments:
        lane = assignment.lane
        if lane.repeat_group_id is None:
            continue
        group_key = (
            lane.repeat_group_id,
            lane.caller_id,
            lane.genome_build.value,
            lane.data_basis.value,
            lane.caller_version,
            lane.bin_size_kbp,
        )
        grouped[group_key].append(assignment)

    pairs: list[CnvReproducibilityPair] = []
    full_by_lane: dict[str, set[str]] = defaultdict(set)
    for record in evidence.full_evidence:
        full_by_lane[_lane_key(record).lane_id].add(record.record_id)

    for group_key in sorted(grouped, key=str):
        members = sorted(
            grouped[group_key],
            key=lambda item: (
                item.lane.specimen_id,
                item.lane.replicate_id,
                item.lane_id,
            ),
        )
        biological_ids = {item.lane.biological_specimen_id for item in members}
        if len(biological_ids) != 1:
            raise ValueError("CNV repeat group mixes different biological specimens")
        repeat_kinds = {item.lane.repeat_kind for item in members}
        if len(repeat_kinds) != 1 or CnvRepeatKind.INDEPENDENT in repeat_kinds:
            raise ValueError("CNV repeat group has inconsistent repeat-kind semantics")
        biological_id = next(iter(biological_ids))
        repeat_kind = next(iter(repeat_kinds))

        for left, right in combinations(members, 2):
            left_events = events_by_lane.get(left.lane_id, [])
            right_events = events_by_lane.get(right.lane_id, [])
            concordance: float | None = None
            matched_pairs: list[list[str]] = []
            if (
                left.run_outcome == CnvRunOutcomeState.OBSERVED
                and right.run_outcome == CnvRunOutcomeState.OBSERVED
                and (left_events or right_events)
            ):
                case = BenchmarkCase(
                    case_id="repeat-"
                    + canonical_evidence_sha256({"left": left.lane_id, "right": right.lane_id})[
                        :24
                    ],
                    kind=BenchmarkKind.CNV,
                    genome_build=left.lane.genome_build,
                    truth_events=[_query_event(event) for event in left_events],
                    query_events=[_query_event(event) for event in right_events],
                    thresholds=registration.matrix.matching_thresholds,
                )
                report = benchmark_case(case)
                concordance = 2 * len(report.matches) / (len(left_events) + len(right_events))
                matched_pairs = [
                    [match.truth_event_id, match.query_event_id] for match in report.matches
                ]

            pair_id = (
                "repeat-pair-"
                + canonical_evidence_sha256({"left": left.lane_id, "right": right.lane_id})[:24]
            )
            pairs.append(
                CnvReproducibilityPair(
                    pair_id=pair_id,
                    repeat_group_id=group_key[0],
                    repeat_kind=repeat_kind,
                    biological_specimen_id=biological_id,
                    caller_id=left.lane.caller_id,
                    genome_build=left.lane.genome_build,
                    data_basis=left.lane.data_basis,
                    bin_size_kbp=left.lane.bin_size_kbp,
                    left_lane_id=left.lane_id,
                    right_lane_id=right.lane_id,
                    left_run_outcome=left.run_outcome,
                    right_run_outcome=right.run_outcome,
                    left_normalized_event_ids=[item.normalized_event_id for item in left_events],
                    right_normalized_event_ids=[item.normalized_event_id for item in right_events],
                    matched_normalized_event_pairs=matched_pairs,
                    concordance=concordance,
                )
            )

    evaluable_pairs = [item for item in pairs if item.concordance is not None]
    all_normalized_ids = {
        event_id
        for item in pairs
        for event_id in item.left_normalized_event_ids + item.right_normalized_event_ids
    }
    all_lane_ids = {
        lane_id for item in pairs for lane_id in (item.left_lane_id, item.right_lane_id)
    }
    all_full_ids = {
        full_id for lane_id in all_lane_ids for full_id in full_by_lane.get(lane_id, set())
    }
    overall_key: dict[str, str | int | float | bool] = {"scope": "overall"}
    if evaluable_pairs:
        mean_concordance = sum(item.concordance or 0.0 for item in evaluable_pairs) / len(
            evaluable_pairs
        )
        overall = CnvMetricResult(
            metric=CnvAcceptanceMetric.REPRODUCIBILITY,
            stratum_key=overall_key,
            state=CnvMetricEvaluationState.EVALUABLE,
            value=mean_concordance,
            evaluable_denominator=len(evaluable_pairs),
            component_counts={
                "evaluable_pairs": len(evaluable_pairs),
                "total_pairs": len(pairs),
            },
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            normalized_event_ids=sorted(all_normalized_ids),
            full_evidence_ids=sorted(all_full_ids),
            lane_ids=sorted(all_lane_ids),
        )
    else:
        overall = CnvMetricResult(
            metric=CnvAcceptanceMetric.REPRODUCIBILITY,
            stratum_key=overall_key,
            state=CnvMetricEvaluationState.NOT_EVALUABLE,
            evaluable_denominator=0,
            component_counts={
                "evaluable_pairs": 0,
                "total_pairs": len(pairs),
            },
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            normalized_event_ids=sorted(all_normalized_ids),
            full_evidence_ids=sorted(all_full_ids),
            lane_ids=sorted(all_lane_ids),
            reason="No evaluable repeated-lane event pairs.",
        )
    return [overall], pairs


class CnvAcceptanceDecision(StrEnum):
    PASS = "PASS"
    FAIL = "FAIL"
    NO_CALL = "NO_CALL"
    NOT_EVALUABLE = "NOT_EVALUABLE"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class CnvAcceptanceResult(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    question_id: str
    metric: CnvAcceptanceMetric
    filters: dict[str, str | int | float | bool]
    decision: CnvAcceptanceDecision
    accepted: bool | None = None
    metric_id: str | None = None
    metric_value: float | None = None
    evaluable_denominator: int = Field(ge=0)
    minimum_evaluable_denominator: int = Field(ge=1)
    minimum_acceptable: float | None = None
    maximum_acceptable: float | None = None
    reason: str | None = None
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_acceptance(self) -> CnvAcceptanceResult:
        if self.metric_value is not None and not math.isfinite(self.metric_value):
            raise ValueError("Acceptance metric value must be finite")
        expected_accepted: bool | None
        if self.decision == CnvAcceptanceDecision.PASS:
            expected_accepted = True
        elif self.decision == CnvAcceptanceDecision.FAIL:
            expected_accepted = False
        else:
            expected_accepted = None
        if self.accepted != expected_accepted:
            raise ValueError("Acceptance decision and accepted flag disagree")
        if self.decision in {CnvAcceptanceDecision.PASS, CnvAcceptanceDecision.FAIL}:
            if self.metric_id is None or self.metric_value is None:
                raise ValueError("PASS/FAIL acceptance requires a concrete metric")
        return self


def _technical_acceptance_stratum(
    filters: dict[str, str | int | float | bool],
) -> dict[str, str | int | float | bool] | None:
    if filters == {"scope": "overall"}:
        return {"scope": "overall"}
    names = (
        "caller",
        "genome_build",
        "data_basis",
        "coverage",
        "tumor_fraction",
        "bin_size",
    )
    if not all(name in filters for name in names):
        return None
    return {name: filters[name] for name in names}


def _evaluate_acceptance_question(
    question: CnvAcceptanceQuestion,
    metrics: Sequence[CnvMetricResult],
) -> CnvAcceptanceResult:
    matches = [
        item
        for item in metrics
        if item.metric == question.metric and item.stratum_key == question.filters
    ]
    if len(matches) > 1:
        raise ValueError(
            f"Acceptance question {question.question_id} matches multiple metric views"
        )
    if not matches:
        return CnvAcceptanceResult(
            question_id=question.question_id,
            metric=question.metric,
            filters=question.filters,
            decision=CnvAcceptanceDecision.NOT_EVALUABLE,
            evaluable_denominator=0,
            minimum_evaluable_denominator=question.minimum_evaluable_denominator,
            minimum_acceptable=question.minimum_acceptable,
            maximum_acceptable=question.maximum_acceptable,
            reason="No metric view exactly matches the registered acceptance filters.",
        )

    metric = matches[0]
    common = {
        "question_id": question.question_id,
        "metric": question.metric,
        "filters": question.filters,
        "metric_id": metric.metric_id,
        "metric_value": metric.value,
        "evaluable_denominator": metric.evaluable_denominator,
        "minimum_evaluable_denominator": question.minimum_evaluable_denominator,
        "minimum_acceptable": question.minimum_acceptable,
        "maximum_acceptable": question.maximum_acceptable,
    }

    if metric.state == CnvMetricEvaluationState.NOT_APPLICABLE:
        return CnvAcceptanceResult(
            **common,
            decision=CnvAcceptanceDecision.NOT_APPLICABLE,
            reason=metric.reason or "Registered metric is not applicable.",
        )

    if metric.state == CnvMetricEvaluationState.NOT_EVALUABLE:
        technical = _technical_acceptance_stratum(question.filters)
        if technical is not None:
            no_call_matches = [
                item
                for item in metrics
                if item.metric == CnvAcceptanceMetric.NO_CALL_RATE
                and item.stratum_key == technical
                and item.state == CnvMetricEvaluationState.EVALUABLE
                and item.evaluable_denominator > 0
                and item.value == 1.0
            ]
            if len(no_call_matches) == 1:
                return CnvAcceptanceResult(
                    **common,
                    decision=CnvAcceptanceDecision.NO_CALL,
                    reason="All executed lanes for this technical stratum were NO_CALL.",
                )
        return CnvAcceptanceResult(
            **common,
            decision=CnvAcceptanceDecision.NOT_EVALUABLE,
            reason=metric.reason or "Metric is not evaluable.",
        )

    if metric.evaluable_denominator < question.minimum_evaluable_denominator:
        return CnvAcceptanceResult(
            **common,
            decision=CnvAcceptanceDecision.NOT_EVALUABLE,
            reason=(
                "Metric denominator is below the prospectively registered "
                "minimum_evaluable_denominator."
            ),
        )

    if metric.value is None:
        raise ValueError("Evaluable metric unexpectedly lacks a value")

    accepted = True
    if question.minimum_acceptable is not None:
        accepted = accepted and metric.value >= question.minimum_acceptable
    if question.maximum_acceptable is not None:
        accepted = accepted and metric.value <= question.maximum_acceptable
    return CnvAcceptanceResult(
        **common,
        decision=CnvAcceptanceDecision.PASS if accepted else CnvAcceptanceDecision.FAIL,
        accepted=accepted,
    )


def _verify_metric_provenance(
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
    assignments: Sequence[CnvLaneEventAssignment],
    metrics: Sequence[CnvMetricResult],
) -> None:
    expected_manifest_sha256 = canonical_evidence_sha256(
        evidence.model_dump(mode="json", exclude={"manifest_sha256"})
    )
    if evidence.manifest_sha256 != expected_manifest_sha256:
        raise ValueError("Evidence manifest changed after its Block-2 content lock")

    full_by_id = {item.record_id: item for item in evidence.full_evidence}
    normalized_by_id = {item.normalized_event_id: item for item in evidence.normalized_events}
    truth_ids = {
        event.event_id
        for specimen in registration.cohort.specimens
        for event in specimen.truth_events
    }
    lane_ids = {item.lane_id for item in assignments}

    for metric in metrics:
        if metric.registration_sha256 != registration.lock_sha256:
            raise ValueError("Metric provenance references a different registration")
        if metric.evidence_manifest_sha256 != evidence.manifest_sha256:
            raise ValueError("Metric provenance references a different evidence manifest")
        missing_truth = set(metric.truth_event_ids) - truth_ids
        if missing_truth:
            raise ValueError("Metric references unknown registered truth events")
        missing_normalized = set(metric.normalized_event_ids) - set(normalized_by_id)
        if missing_normalized:
            raise ValueError("Metric references unknown normalized CNV events")
        missing_full = set(metric.full_evidence_ids) - set(full_by_id)
        if missing_full:
            raise ValueError("Metric references unknown full-evidence records")
        missing_lanes = set(metric.lane_ids) - lane_ids
        if missing_lanes:
            raise ValueError("Metric references unknown analytical lanes")
        metric_full_ids = set(metric.full_evidence_ids)
        for normalized_id in metric.normalized_event_ids:
            sources = set(normalized_by_id[normalized_id].source_full_evidence_ids)
            if not sources.issubset(metric_full_ids):
                raise ValueError(
                    "Metric normalized-event provenance is missing its Block-2 source evidence"
                )


class CnvValidationMetricReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    report_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,127}$")
    registration_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    evidence_manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    matching_thresholds: BenchmarkThresholds
    assignments: list[CnvLaneEventAssignment]
    metrics: list[CnvMetricResult]
    acceptance_results: list[CnvAcceptanceResult]
    negative_unit_assessments: list[CnvNegativeUnitAssessment] = Field(default_factory=list)
    reproducibility_pairs: list[CnvReproducibilityPair] = Field(default_factory=list)
    report_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    retain_all_evidence: Literal[True] = True
    clinical_validity_claimed: Literal[False] = False
    sensitive_output: Literal[True] = True
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_report(self) -> CnvValidationMetricReport:
        lane_ids = [item.lane_id for item in self.assignments]
        if len(lane_ids) != len(set(lane_ids)):
            raise ValueError("Validation report lane assignments must be unique")
        metric_ids = [item.metric_id for item in self.metrics]
        if len(metric_ids) != len(set(metric_ids)):
            raise ValueError("Validation report metric views must be unique")
        question_ids = [item.question_id for item in self.acceptance_results]
        if len(question_ids) != len(set(question_ids)):
            raise ValueError("Validation report acceptance question IDs must be unique")
        metric_id_set = set(metric_ids)
        for result in self.acceptance_results:
            if result.metric_id is not None and result.metric_id not in metric_id_set:
                raise ValueError("Acceptance result references a metric outside the report")
        for metric in self.metrics:
            if metric.registration_sha256 != self.registration_sha256:
                raise ValueError("Report contains metric from another registration")
            if metric.evidence_manifest_sha256 != self.evidence_manifest_sha256:
                raise ValueError("Report contains metric from another evidence manifest")
        for assessment in self.negative_unit_assessments:
            if assessment.registration_sha256 != self.registration_sha256:
                raise ValueError("Report contains specificity assessment from another registration")
            if assessment.evidence_manifest_sha256 != self.evidence_manifest_sha256:
                raise ValueError(
                    "Report contains specificity assessment from another evidence manifest"
                )
        for assignment in self.assignments:
            if (
                assignment.minimum_reciprocal_overlap
                != self.matching_thresholds.minimum_reciprocal_overlap
                or assignment.copy_number_tolerance
                != self.matching_thresholds.copy_number_tolerance
            ):
                raise ValueError("Report assignment matching policy differs from report lock")

        expected_sha256 = canonical_evidence_sha256(
            self.model_dump(mode="json", exclude={"report_sha256"})
        )
        if self.report_sha256 != expected_sha256:
            raise ValueError("Validation report content does not match report_sha256")
        return self


def aggregate_cnv_validation(
    registration: CnvValidationRegistration,
    evidence: CnvValidationEvidenceManifest,
    *,
    report_id: str,
    negative_assessments: Sequence[CnvNegativeUnitAssessment] = (),
) -> CnvValidationMetricReport:
    """Build a deterministic RUO metric report from one locked registration/evidence pair."""
    assignments = assign_cnv_validation_lanes(registration, evidence)
    event_metrics = aggregate_cnv_event_metrics(
        registration,
        evidence,
        assignments=assignments,
    )
    run_metrics = aggregate_cnv_run_state_metrics(
        registration,
        evidence,
        negative_assessments=negative_assessments,
        assignments=assignments,
    )
    quantitative_metrics = aggregate_cnv_quantitative_metrics(
        registration,
        evidence,
        assignments=assignments,
    )
    reproducibility_metrics, reproducibility_pairs = aggregate_cnv_reproducibility_metrics(
        registration,
        evidence,
        assignments=assignments,
    )
    metrics = [
        *event_metrics,
        *run_metrics,
        *quantitative_metrics,
        *reproducibility_metrics,
    ]
    metrics.sort(
        key=lambda item: (
            item.metric.value,
            canonical_evidence_sha256(item.stratum_key),
        )
    )

    _verify_metric_provenance(registration, evidence, assignments, metrics)
    acceptance_results = [
        _evaluate_acceptance_question(question, metrics)
        for question in registration.matrix.acceptance
    ]

    payload = {
        "schema_version": "0.1.0",
        "report_id": report_id,
        "registration_sha256": registration.lock_sha256,
        "evidence_manifest_sha256": evidence.manifest_sha256,
        "matching_thresholds": registration.matrix.matching_thresholds.model_dump(mode="json"),
        "assignments": [item.model_dump(mode="json") for item in assignments],
        "metrics": [item.model_dump(mode="json") for item in metrics],
        "acceptance_results": [
            item.model_dump(mode="json") for item in acceptance_results
        ],
        "negative_unit_assessments": [
            item.model_dump(mode="json") for item in negative_assessments
        ],
        "reproducibility_pairs": [
            item.model_dump(mode="json") for item in reproducibility_pairs
        ],
        "retain_all_evidence": True,
        "clinical_validity_claimed": False,
        "sensitive_output": True,
        "research_only": True,
    }
    payload["report_sha256"] = canonical_evidence_sha256(payload)
    return CnvValidationMetricReport.model_validate(payload)
