from __future__ import annotations

import math
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from .models import StrictModel
from .multicaller_contracts import (
    ID,
    CallerMode,
    CallerPlanningDecision,
    MultiCallerPlan,
)
from .multicaller_identity import COLLECTION_ID, MultiCallerCollectionAddress


class CallerExecutionOutcome(StrEnum):
    COMPLETED = "COMPLETED"
    NO_CALL = "NO_CALL"
    FAILED = "FAILED"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"
    NOT_RUN = "NOT_RUN"


class CallerDependencyStatus(StrEnum):
    NO_DECLARED_CALLER_DEPENDENCY = "NO_DECLARED_CALLER_DEPENDENCY"
    DIRECT_CALLER_DEPENDENCY = "DIRECT_CALLER_DEPENDENCY"


class ComparisonGroupState(StrEnum):
    AGREEING = "AGREEING"
    CONFLICTING = "CONFLICTING"
    SINGLE_SOURCE = "SINGLE_SOURCE"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"


class MultiCallerEvidenceLink(StrictModel):
    address: MultiCallerCollectionAddress
    label: str = Field(min_length=1)
    native_artifact_ids: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def unique_native_artifact_ids(self) -> MultiCallerEvidenceLink:
        if len(self.native_artifact_ids) != len(set(self.native_artifact_ids)):
            raise ValueError("Multi-caller evidence native artifact IDs must be unique")
        if any(not item.strip() for item in self.native_artifact_ids):
            raise ValueError("Multi-caller evidence native artifact IDs cannot be empty")
        return self


class CallerLaneExecutionRecord(StrictModel):
    lane_id: str = Field(pattern=ID)
    outcome: CallerExecutionOutcome
    reasons: list[str] = Field(default_factory=list)
    evidence_links: list[MultiCallerEvidenceLink] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_execution_record(self) -> CallerLaneExecutionRecord:
        if self.outcome != CallerExecutionOutcome.COMPLETED and not self.reasons:
            raise ValueError("Non-completed caller execution requires an explicit reason")
        collection_ids = [item.address.collection_id for item in self.evidence_links]
        if len(collection_ids) != len(set(collection_ids)):
            raise ValueError("Caller execution evidence links must be unique")
        return self


class MultiCallerLaneComparison(StrictModel):
    lane_id: str = Field(pattern=ID)
    mode_id: CallerMode
    planning_decision: CallerPlanningDecision
    planning_reasons: list[str] = Field(default_factory=list)
    execution_outcome: CallerExecutionOutcome
    execution_reasons: list[str] = Field(default_factory=list)
    dependency_status: CallerDependencyStatus
    parent_lane_ids: list[str] = Field(default_factory=list)
    parent_lane_sha256s: dict[str, str] = Field(default_factory=dict)
    lane_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    evidence_links: list[MultiCallerEvidenceLink] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_dependencies(self) -> MultiCallerLaneComparison:
        if set(self.parent_lane_ids) != set(self.parent_lane_sha256s):
            raise ValueError("Comparison lane parent IDs and parent hashes must match")
        expected = (
            CallerDependencyStatus.DIRECT_CALLER_DEPENDENCY
            if self.parent_lane_ids
            else CallerDependencyStatus.NO_DECLARED_CALLER_DEPENDENCY
        )
        if self.dependency_status != expected:
            raise ValueError("Comparison lane dependency status does not match parent lanes")
        for link in self.evidence_links:
            if link.address.lane_sha256 != self.lane_sha256:
                raise ValueError("Comparison evidence link lane SHA does not match planned lane SHA")
        return self


class MultiCallerComparisonMember(StrictModel):
    lane_id: str = Field(pattern=ID)
    evidence_collection_id: str = Field(pattern=COLLECTION_ID)
    observation: str = Field(min_length=1)
    numeric_value: float | None = None
    text_value: str | None = None
    unit: str | None = None
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def one_observation_value(self) -> MultiCallerComparisonMember:
        if (self.numeric_value is None) == (self.text_value is None):
            raise ValueError(
                "Comparison member requires exactly one numeric or text observation value"
            )
        if self.numeric_value is not None and not math.isfinite(self.numeric_value):
            raise ValueError("Comparison member numeric value must be finite")
        if self.text_value is not None and not self.text_value.strip():
            raise ValueError("Comparison member text value cannot be empty")
        return self


class MultiCallerComparisonGroup(StrictModel):
    group_id: str = Field(pattern=ID)
    subject: str = Field(min_length=1)
    state: ComparisonGroupState
    members: list[MultiCallerComparisonMember] = Field(min_length=1)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def unique_members(self) -> MultiCallerComparisonGroup:
        addresses = [
            (item.lane_id, item.evidence_collection_id, item.observation) for item in self.members
        ]
        if len(addresses) != len(set(addresses)):
            raise ValueError("Comparison group contains duplicate evidence observations")
        return self


class MultiCallerComparisonReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    plan_id: str = Field(pattern=ID)
    plan_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    lanes: list[MultiCallerLaneComparison] = Field(min_length=1)
    groups: list[MultiCallerComparisonGroup] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def unique_report_membership(self) -> MultiCallerComparisonReport:
        lane_ids = [item.lane_id for item in self.lanes]
        if len(lane_ids) != len(set(lane_ids)):
            raise ValueError("Comparison report lane IDs must be unique")
        group_ids = [item.group_id for item in self.groups]
        if len(group_ids) != len(set(group_ids)):
            raise ValueError("Comparison report group IDs must be unique")
        return self


def _validate_execution_against_plan(
    *,
    planning_decision: CallerPlanningDecision,
    execution: CallerLaneExecutionRecord,
) -> None:
    if (
        planning_decision != CallerPlanningDecision.ELIGIBLE
        and execution.outcome != CallerExecutionOutcome.NOT_RUN
    ):
        raise ValueError(
            "Caller lane that was not eligible for execution must remain NOT_RUN in comparison"
        )


def build_multicaller_comparison_report(
    *,
    plan: MultiCallerPlan,
    executions: list[CallerLaneExecutionRecord],
    groups: list[MultiCallerComparisonGroup] | None = None,
) -> MultiCallerComparisonReport:
    execution_ids = [item.lane_id for item in executions]
    if len(execution_ids) != len(set(execution_ids)):
        raise ValueError("Multi-caller comparison execution lane IDs must be unique")

    plan_by_id = {item.lane_id: item for item in plan.lanes}
    execution_by_id = {item.lane_id: item for item in executions}
    if set(execution_by_id) != set(plan_by_id):
        raise ValueError("Comparison executions must contain exactly every planned caller lane")

    lanes: list[MultiCallerLaneComparison] = []
    reachable_by_lane: dict[str, set[str]] = {}
    all_collection_ids: set[str] = set()
    for lane_id in plan.topological_order:
        planned = plan_by_id[lane_id]
        execution = execution_by_id[lane_id]
        _validate_execution_against_plan(
            planning_decision=planned.decision,
            execution=execution,
        )

        lane_collection_ids: set[str] = set()
        for link in execution.evidence_links:
            if link.address.lane_sha256 != planned.lane_sha256:
                raise ValueError(
                    f"Evidence link lane SHA does not match planned lane SHA for {lane_id!r}"
                )
            collection_id = link.address.collection_id
            if collection_id in all_collection_ids:
                raise ValueError("Comparison evidence collection IDs must be globally unique")
            all_collection_ids.add(collection_id)
            lane_collection_ids.add(collection_id)
        reachable_by_lane[lane_id] = lane_collection_ids

        lanes.append(
            MultiCallerLaneComparison(
                lane_id=lane_id,
                mode_id=planned.mode_id,
                planning_decision=planned.decision,
                planning_reasons=planned.reasons,
                execution_outcome=execution.outcome,
                execution_reasons=execution.reasons,
                dependency_status=(
                    CallerDependencyStatus.DIRECT_CALLER_DEPENDENCY
                    if planned.parent_lane_ids
                    else CallerDependencyStatus.NO_DECLARED_CALLER_DEPENDENCY
                ),
                parent_lane_ids=planned.parent_lane_ids,
                parent_lane_sha256s=planned.parent_lane_sha256s,
                lane_sha256=planned.lane_sha256,
                evidence_links=execution.evidence_links,
            )
        )

    comparison_groups = groups or []
    for group in comparison_groups:
        for member in group.members:
            if member.lane_id not in plan_by_id:
                raise ValueError(
                    f"Comparison group {group.group_id!r} references unknown lane "
                    f"{member.lane_id!r}"
                )
            if member.evidence_collection_id not in reachable_by_lane[member.lane_id]:
                raise ValueError(
                    f"Comparison group {group.group_id!r} references unreachable evidence "
                    f"for lane {member.lane_id!r}"
                )

    return MultiCallerComparisonReport(
        plan_id=plan.plan_id,
        plan_sha256=plan.plan_sha256,
        lanes=lanes,
        groups=comparison_groups,
    )
