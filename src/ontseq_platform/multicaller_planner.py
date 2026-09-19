from __future__ import annotations

from .multicaller_contracts import (
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
