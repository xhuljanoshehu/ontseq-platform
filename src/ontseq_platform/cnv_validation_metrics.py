from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from itertools import product
from typing import Literal

from pydantic import Field, model_validator

from .benchmark import benchmark_case
from .cnv_validation_contracts import CnvAcceptanceMetric, CnvDataBasis, CnvRepeatKind
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
    registration_sha256: str
    evidence_manifest_sha256: str
    truth_event_ids: list[str] = Field(default_factory=list)
    normalized_event_ids: list[str] = Field(default_factory=list)
    full_evidence_ids: list[str] = Field(default_factory=list)
    lane_ids: list[str] = Field(default_factory=list)
    reason: str | None = None
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_metric(self) -> CnvMetricResult:
        if self.state == CnvMetricEvaluationState.EVALUABLE:
            if self.value is None or not math.isfinite(self.value):
                raise ValueError("Evaluable CNV metric requires a finite value")
        elif self.value is not None:
            raise ValueError("Non-evaluable CNV metric cannot carry a value")
        if any(value < 0 for value in self.component_counts.values()):
            raise ValueError("CNV metric component counts cannot be negative")
        return self


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
) -> list[CnvMetricResult]:
    """Aggregate locked lane assignments without re-matching inside strata."""
    assignments = assign_cnv_validation_lanes(registration, evidence)
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
