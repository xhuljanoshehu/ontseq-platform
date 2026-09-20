"""One-time owner-authorized merge repair; executed only on isolated CI checkouts."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

BASE = "ceaab08079b8cb7c528eac6f6137ad2e4dac5f31"
TARGETS = {
    78: ("feat/marlin-dual-runtime-compatibility", "16f2ef8e9dcb63867176749131d7e0ba1266f502"),
    81: ("feat/multicaller-routing", "0cdfd6383da0e1869b4e792849a7056d2466f602"),
    82: ("feat/multicaller-real-tool-qualification", "30521bebc67f2562c26120a8f32983dfaaaad0d7"),
}
CONFLICTS = {
    78: {"pyproject.toml"},
    81: {"docs/superpowers/specs/2026-09-19-multicaller-integration-design.md", "schemas/multicaller-plan.schema.json", "scripts/export_schemas.py"},
    82: {"tests/test_severus_real_tool.py", "workflow/envs/severus_1_7.yaml"},
}


def run(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    p = subprocess.run(args, text=True, capture_output=True)
    print("$", " ".join(args), flush=True)
    if p.stdout:
        print(p.stdout, flush=True)
    if p.stderr:
        print(p.stderr, flush=True)
    if check and p.returncode:
        raise RuntimeError(f"Command failed ({p.returncode}): {args[0]}")
    return p


def content(ref: str, path: str) -> str:
    return subprocess.check_output(["git", "show", f"{ref}:{path}"], text=True)


def write(path: str, text: str) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")


def remote(branch: str) -> str:
    return subprocess.check_output(["git", "ls-remote", "origin", f"refs/heads/{branch}"], text=True).split()[0]


def prepare(pr: int) -> None:
    branch, head = TARGETS[pr]
    if remote("main") != BASE or remote(branch) != head:
        raise RuntimeError("A reviewed branch moved; refusing stale repair")
    run("git", "checkout", "--detach", head)
    result = run("git", "merge", "--no-commit", "--no-ff", BASE, check=False)
    observed = set(subprocess.check_output(["git", "diff", "--name-only", "--diff-filter=U"], text=True).splitlines())
    if result.returncode != 1 or observed != CONFLICTS[pr]:
        raise RuntimeError(f"Unexpected conflict set: {sorted(observed)}")

    if pr == 78:
        path = "pyproject.toml"
        text = Path(path).read_text()
        pattern = r"<<<<<<<[^\n]*\n(.*?)=======\n(.*?)>>>>>>>[^\n]*\n"
        matches = list(re.finditer(pattern, text, flags=re.S))
        expected = {
            '    "docs/superpowers/plans/2026-09-17-cnv-validation-block1.md",',
            '    "docs/superpowers/plans/2026-09-15-marlin-dual-runtime-compatibility.md",',
        }
        if len(matches) != 1 or set((matches[0][1] + matches[0][2]).splitlines()) != expected:
            raise RuntimeError("Unexpected pyproject conflict contents")
        write(path, re.sub(pattern, lambda m: m[1] + m[2], text, flags=re.S))
        note = "Retain both Ruff plan exclusions and the MARLIN R probe package resource. All nonconflicting dual-runtime compatibility/evidence gates are preserved."

    elif pr == 81:
        spec = "docs/superpowers/specs/2026-09-19-multicaller-integration-design.md"
        routing_spec = "docs/superpowers/specs/2026-09-19-multicaller-routing-5a-design.md"
        scope = (
            "# Supplemental declarative routing contract\n\n"
            "Merge scope (2026-09-20): this is a separate, planning-only contract. "
            "It does not replace multicaller_contracts.MultiCallerPlan, the canonical lane "
            "planner, caller runtimes, or their qualification registry. Its existing_adapter "
            "flags describe this supplemental route, not platform-wide installation status. "
            "All its lanes remain NOT_RUN and execution_enabled=false. No implicit conversion "
            "into the canonical execution plan is provided.\n\n"
        )
        write(routing_spec, scope + content(head, spec))
        write(spec, content(BASE, spec))
        write("schemas/multicaller-plan.schema.json", content(BASE, "schemas/multicaller-plan.schema.json"))
        write("scripts/export_schemas.py", content(BASE, "scripts/export_schemas.py"))
        for path in ("docs/MULTICALLER_5A_STATUS.md", "docs/superpowers/plans/2026-09-19-multicaller-routing-5a.md"):
            text = Path(path).read_text()
            text = text.replace("2026-09-19-multicaller-integration-design.md", "2026-09-19-multicaller-routing-5a-design.md")
            text = text.replace("multicaller-plan.schema.json", "multicaller-routing-plan.schema.json")
            write(path, scope + text)
        path = "tests/test_multicaller_plan.py"
        text = Path(path).read_text()
        old = '("multicaller-plan", MultiCallerPlan)'
        if old not in text:
            raise RuntimeError("Expected alternate schema assertion not found")
        write(path, text.replace(old, '("multicaller-routing-plan", MultiCallerPlan)'))
        path = "src/ontseq_platform/multicaller_plan.py"
        text = Path(path).read_text()
        text = text.replace(
            "This module reads no genomic files and executes no tools.",
            "This supplemental declaration contract is not the canonical execution planner.\n"
            "Its public schema is multicaller-routing-plan.schema.json; no runtime dispatch\n"
            "or conversion to multicaller_contracts.MultiCallerPlan is implied.\n\n"
            "This module reads no genomic files and executes no tools.",
        )
        write(path, text)
        write("tests/test_multicaller_schema_coexistence.py", '''from __future__ import annotations

import json
from pathlib import Path
import runpy

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
        assert json.loads(exported[Path(f"schemas/{name}.schema.json")]) == model.model_json_schema()
    assert ExecutionPlan is not RoutingPlan
    assert ExecutionPlan.model_json_schema() != RoutingPlan.model_json_schema()
    assert "request" in RoutingPlan.model_json_schema()["properties"]
    assert "request" not in ExecutionPlan.model_json_schema()["properties"]
''')
        note = "Preserve the canonical execution plan/schema unchanged. Retain the supplemental declarative input eligibility planner under multicaller-routing-plan.schema.json with distinct exports and scope documentation. Retain its immutable input/provenance and no-execution contracts."

    else:
        path = "tests/test_severus_real_tool.py"
        new_path = "tests/test_severus_long_read_real_tool.py"
        if Path(new_path).exists():
            raise RuntimeError("Additional Severus smoke path already exists")
        write(new_path, content(head, path))
        write(path, content(BASE, path))
        write("workflow/envs/severus_1_7.yaml", content(BASE, "workflow/envs/severus_1_7.yaml"))
        path = ".github/workflows/sv-severus.yml"
        text = Path(path).read_text().replace(
            '      - "tests/test_severus_real_tool.py"',
            '      - "tests/test_severus_real_tool.py"\n      - "tests/test_severus_long_read_real_tool.py"',
        )
        text = text.replace(
            "python -m pytest -q tests/test_severus_runtime.py tests/test_severus_real_tool.py",
            "python -m pytest -q tests/test_severus_runtime.py tests/test_severus_real_tool.py tests/test_severus_long_read_real_tool.py",
        )
        write(path, text)
        note = "Preserve the existing Severus smoke and add the independent 15-kb long-read synthetic fixture, including its nonzero control-depth threshold and BAM hash checks. Keep the bounded main pysam dependency; run both smoke tests in the supplemental workflow."

    write(f"docs/merge-repairs/PR-{pr}.md", f"# PR #{pr} merge repair\n\nOwner-authorized repair on 2026-09-20.\n\nBase: `{BASE}`\nOriginal head: `{head}`\n\n{note}\n\nVersion remains 0.8.2. Research Use Only. No analytical or clinical validation is implied. No thresholds were relaxed to achieve a passing test.\n")
    run("git", "add", "--all")
    print(json.dumps({"pr": pr, "expected_head": head, "resolved_conflicts": sorted(observed)}), flush=True)


def finalize(pr: int) -> None:
    if pr == 81:
        # The new regression must fail while the canonical exporter lacks the extra schema.
        red = run(sys.executable, "-m", "pytest", "-q", "tests/test_multicaller_schema_coexistence.py", check=False)
        if red.returncode != 1 or "multicaller-routing-plan.schema.json" not in red.stdout or "KeyError" not in red.stdout:
            raise RuntimeError("Expected schema namespace RED was not observed")
        path = "scripts/export_schemas.py"
        text = Path(path).read_text()
        anchor = "from ontseq_platform.multicaller_contracts import MultiCallerPlan\n"
        addition = (
            "from ontseq_platform.multicaller_plan import MultiCallerPlan as MultiCallerRoutingPlan\n"
            "from ontseq_platform.multicaller_plan import MultiCallerRequest\n"
        )
        if text.count(anchor) != 1:
            raise RuntimeError("Canonical export import changed")
        text = text.replace(anchor, anchor + addition)
        anchor = '                "multicaller-plan": MultiCallerPlan,\n'
        if text.count(anchor) != 1:
            raise RuntimeError("Canonical export entry changed")
        text = text.replace(anchor, anchor + '                "multicaller-routing-plan": MultiCallerRoutingPlan,\n                "multicaller-request": MultiCallerRequest,\n')
        write(path, text)
        run(sys.executable, "scripts/export_schemas.py")
        run("ruff", "check", "--fix", "scripts/export_schemas.py", "tests/test_multicaller_schema_coexistence.py")
        run("ruff", "format", "scripts/export_schemas.py", "tests/test_multicaller_schema_coexistence.py", "src/ontseq_platform/multicaller_plan.py", "tests/test_multicaller_plan.py")
        run(sys.executable, "-m", "pytest", "-q", "tests/test_multicaller_schema_coexistence.py", "tests/test_multicaller_plan.py")
        # A supplemental schema must not change any previously exported schema.
        existing = subprocess.check_output(["git", "ls-tree", "-r", "--name-only", BASE, "schemas"], text=True).splitlines()
        for path in existing:
            if Path(path).read_text() != content(BASE, path):
                raise RuntimeError(f"Existing public schema changed unexpectedly: {path}")
    run("git", "add", "--all")
    run("git", "diff", "--cached", "--check")
    run("git", "commit", "-m", f"Merge main into PR #{pr}: preserve compatible additions and resolve conflicts")
    run("git", "show", "--stat", "--oneline", "HEAD")


def promote(pr: int) -> None:
    branch, old_head = TARGETS[pr]
    if remote("main") != BASE or remote(branch) != old_head:
        raise RuntimeError("Remote changed since repair review; refusing promotion")
    run("git", "merge-base", "--is-ancestor", old_head, "HEAD")
    run("git", "merge-base", "--is-ancestor", BASE, "HEAD")
    run("git", "push", "origin", f"HEAD:refs/heads/{branch}")
    print(json.dumps({"pr": pr, "branch": branch, "tested_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["prepare", "finalize", "promote"])
    parser.add_argument("pr", type=int, choices=sorted(TARGETS))
    args = parser.parse_args()
    {"prepare": prepare, "finalize": finalize, "promote": promote}[args.stage](args.pr)
