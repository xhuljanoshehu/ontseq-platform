from __future__ import annotations

import json
import runpy
from pathlib import Path

from ontseq_platform.multicaller_contracts import MultiCallerPlan as ExecutionPlan
from ontseq_platform.multicaller_plan import MultiCallerPlan as RoutingPlan
from ontseq_platform.multicaller_plan import MultiCallerRequest


def test_canonical_and_supplemental_routing_schemas_do_not_collide() -> None:
    render = runpy.run_path("scripts/export_schemas.py")["_render"]
    exported = render()
    contracts = {
        "multicaller-plan": ExecutionPlan,
        "multicaller-routing-plan": RoutingPlan,
        "multicaller-request": MultiCallerRequest,
    }
    for name, model in contracts.items():
        exported_schema = json.loads(exported[Path(f"schemas/{name}.schema.json")])
        assert exported_schema == model.model_json_schema()
    assert ExecutionPlan is not RoutingPlan
    assert ExecutionPlan.model_json_schema() != RoutingPlan.model_json_schema()
    assert "request" in RoutingPlan.model_json_schema()["properties"]
    assert "request" not in ExecutionPlan.model_json_schema()["properties"]
