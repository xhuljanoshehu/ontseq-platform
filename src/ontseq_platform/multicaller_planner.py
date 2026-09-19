from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from typing import Any

from .multicaller_contracts import (
    CallerMode,
    CallerPlanningDecision,
    MultiCallerPlan,
    SealedCallerLanePlan,
    CallerLanePlan,
    CallerLaneRequest,
    CallerPlanningDecision,
    get_caller_catalog_entry,
)


def plan_caller_lane(
    request: CallerLaneRequest,
    *,
    runtime_available: bool,
) -> CallerLanePlan:
    """Plan one caller lane without executing it or changing its requested mode."""
    if not request.requested:
        return CallerLanePlan(
            lane_id=request.lane_id,
            mode_id=request.mode_id,
            decision=CallerPlanningDecision.NOT_REQUESTED,
            reasons=["Caller lane was not requested."],
            input_artifact_ids=[item.artifact_id for item in request.inputs],
            parent_lane_ids=request.parent_lane_ids,
        )

    entry = get_caller_catalog_entry(request.mode_id)
    reasons: list[str] = []
    if request.assay_regime not in entry.compatible_regimes:
        reasons.append(
            f"Assay regime {request.assay_regime.value!r} is not compatible with "
            f"{request.mode_id.value}."
        )

    present_roles = {item.role for item in request.inputs}
    for role in entry.required_input_roles:
        if role not in present_roles:
            reasons.append(f"Missing required input: {role.value.replace('_', ' ')}.")

    if reasons:
        return CallerLanePlan(
            lane_id=request.lane_id,
            mode_id=request.mode_id,
            decision=CallerPlanningDecision.INELIGIBLE,
            reasons=reasons,
            input_artifact_ids=[item.artifact_id for item in request.inputs],
            parent_lane_ids=request.parent_lane_ids,
        )

    if not runtime_available:
        return CallerLanePlan(
            lane_id=request.lane_id,
            mode_id=request.mode_id,
            decision=CallerPlanningDecision.RUNTIME_UNAVAILABLE,
            reasons=[f"Runtime for {request.mode_id.value} is not available."],
            input_artifact_ids=[item.artifact_id for item in request.inputs],
            parent_lane_ids=request.parent_lane_ids,
        )

    return CallerLanePlan(
        lane_id=request.lane_id,
        mode_id=request.mode_id,
        decision=CallerPlanningDecision.ELIGIBLE,
        input_artifact_ids=[item.artifact_id for item in request.inputs],
        parent_lane_ids=request.parent_lane_ids,
    )


def _canonical_sha256(value: Mapping[str, Any]) -> str:
    rendered = json.dumps(
        dict(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )
    return hashlib.sha256(rendered.encode("utf-8")).hexdigest()


def _canonical_request(request: CallerLaneRequest) -> CallerLaneRequest:
    return request.model_copy(
        update={
            "inputs": sorted(request.inputs, key=lambda item: (item.role.value, item.artifact_id)),
            "parent_lane_ids": sorted(request.parent_lane_ids),
        }
    )


def _topological_order(requests: list[CallerLaneRequest]) -> list[str]:
    lane_ids = [item.lane_id for item in requests]
    if len(lane_ids) != len(set(lane_ids)):
        raise ValueError("Multi-caller lane IDs must be unique")

    by_id = {item.lane_id: item for item in requests}
    indegree = {lane_id: 0 for lane_id in lane_ids}
    children: dict[str, set[str]] = {lane_id: set() for lane_id in lane_ids}
    for request in requests:
        for parent_id in request.parent_lane_ids:
            if parent_id not in by_id:
                raise ValueError(
                    f"Caller lane {request.lane_id!r} references missing parent lane "
                    f"{parent_id!r}"
                )
            indegree[request.lane_id] += 1
            children[parent_id].add(request.lane_id)

    ready = sorted(lane_id for lane_id, degree in indegree.items() if degree == 0)
    ordered: list[str] = []
    while ready:
        current = ready.pop(0)
        ordered.append(current)
        for child_id in sorted(children[current]):
            indegree[child_id] -= 1
            if indegree[child_id] == 0:
                ready.append(child_id)
                ready.sort()

    if len(ordered) != len(requests):
        raise ValueError("Multi-caller dependency cycle detected")
    return ordered


def _validate_parent_mode_contracts(
    requests: list[CallerLaneRequest],
    by_id: dict[str, CallerLaneRequest],
) -> None:
    for request in requests:
        if not request.parent_lane_ids:
            continue
        entry = get_caller_catalog_entry(request.mode_id)
        if not entry.required_parent_modes:
            raise ValueError(
                f"{request.mode_id.value} does not accept parent caller lanes"
            )
        allowed = set(entry.required_parent_modes)
        parent_modes = {by_id[parent_id].mode_id for parent_id in request.parent_lane_ids}
        unexpected = parent_modes - allowed
        if unexpected:
            rendered = ", ".join(sorted(item.value for item in unexpected))
            raise ValueError(
                f"{request.mode_id.value} received unsupported parent mode(s): {rendered}"
            )


def _planned_lane_with_dependencies(
    request: CallerLaneRequest,
    *,
    runtime_available: bool,
    request_by_id: dict[str, CallerLaneRequest],
    sealed_by_id: dict[str, SealedCallerLanePlan],
) -> SealedCallerLanePlan:
    base = plan_caller_lane(request, runtime_available=runtime_available)
    entry = get_caller_catalog_entry(request.mode_id)
    reasons = list(base.reasons)
    decision = base.decision

    if decision == CallerPlanningDecision.ELIGIBLE:
        parent_modes = {
            request_by_id[parent_id].mode_id: parent_id
            for parent_id in request.parent_lane_ids
        }
        for required_mode in entry.required_parent_modes:
            parent_id = parent_modes.get(required_mode)
            if parent_id is None:
                reasons.append(
                    f"Missing required parent caller mode: {required_mode.value}."
                )
                continue
            parent_plan = sealed_by_id[parent_id]
            if parent_plan.decision != CallerPlanningDecision.ELIGIBLE:
                reasons.append(
                    f"Required parent lane {parent_id!r} is not eligible for execution."
                )
        if reasons:
            decision = CallerPlanningDecision.INELIGIBLE

    parent_hashes = {
        parent_id: sealed_by_id[parent_id].lane_sha256
        for parent_id in sorted(request.parent_lane_ids)
    }
    lane_payload: dict[str, Any] = {
        "request": request.model_dump(mode="json"),
        "decision": decision.value,
        "reasons": reasons,
        "parent_lane_sha256s": parent_hashes,
    }
    lane_sha256 = _canonical_sha256(lane_payload)
    return SealedCallerLanePlan(
        lane_id=base.lane_id,
        mode_id=base.mode_id,
        decision=decision,
        reasons=reasons,
        input_artifact_ids=base.input_artifact_ids,
        parent_lane_ids=base.parent_lane_ids,
        parent_lane_sha256s=parent_hashes,
        lane_sha256=lane_sha256,
    )


def seal_multicaller_plan(
    *,
    plan_id: str,
    requests: list[CallerLaneRequest],
    runtime_availability: Mapping[CallerMode, bool],
) -> MultiCallerPlan:
    """Seal a deterministic, dependency-aware plan without executing any caller."""
    canonical_requests = sorted(
        (_canonical_request(item) for item in requests),
        key=lambda item: item.lane_id,
    )
    order = _topological_order(canonical_requests)
    request_by_id = {item.lane_id: item for item in canonical_requests}
    _validate_parent_mode_contracts(canonical_requests, request_by_id)

    sealed_by_id: dict[str, SealedCallerLanePlan] = {}
    for lane_id in order:
        request = request_by_id[lane_id]
        sealed_by_id[lane_id] = _planned_lane_with_dependencies(
            request,
            runtime_available=runtime_availability.get(request.mode_id, False),
            request_by_id=request_by_id,
            sealed_by_id=sealed_by_id,
        )

    lanes = [sealed_by_id[lane_id] for lane_id in order]
    payload: dict[str, Any] = {
        "schema_version": "0.1.0",
        "plan_id": plan_id,
        "requests": [item.model_dump(mode="json") for item in canonical_requests],
        "lanes": [item.model_dump(mode="json") for item in lanes],
        "topological_order": order,
        "research_only": True,
    }
    payload["plan_sha256"] = _canonical_sha256(payload)
    return MultiCallerPlan.model_validate(payload)
