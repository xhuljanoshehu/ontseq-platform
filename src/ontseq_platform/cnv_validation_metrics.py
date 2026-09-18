from __future__ import annotations

from collections import defaultdict
from typing import Literal

from pydantic import Field

from .benchmark import benchmark_case
from .cnv_validation_contracts import CnvDataBasis, CnvRepeatKind
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
    for field, expected in expected_caller.items():
        if getattr(record, field) != expected:
            raise ValueError(f"Evidence caller identity mismatch for {field}")

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
    for field, expected in expected_specimen.items():
        if getattr(record, field) != expected:
            raise ValueError(f"Evidence specimen identity mismatch for {field}")

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
