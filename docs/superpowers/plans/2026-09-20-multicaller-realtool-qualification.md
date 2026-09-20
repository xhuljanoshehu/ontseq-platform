# Multi-Caller Real-Tool Qualification Plan

> Branch: `feat/multicaller-realtool-qualification`
> Base: `48614a1c41c5b7480f15dfe6246be53c72419499`
> Method: strict RED -> GREEN

## Task R1 — Severus 1.7 environment and RED real-tool contract

Files:
- `workflow/envs/severus_1_7.yaml`
- `tests/test_severus_real_tool.py`
- `.github/workflows/tumor-severus.yml`

RED acceptance:
- current real binary is present and version is 1.7;
- generated synthetic tumor/normal BAMs are indexed;
- real execution goes through `run_severus_paired`;
- report is COMPLETED or NO_CALL, never a biological-negative assertion;
- native VCF/artifacts are fingerprinted;
- no patient payload is used.

## Task R2 — Severus runtime compatibility

If R1 fails:
- reproduce the exact mismatch with a focused adapter test;
- fix minimally;
- rerun focused unit tests + Severus real-tool smoke;
- retain paired fail-closed semantics.

## Task R3 — Severus qualification status

After same-head green:
- mark only the pinned paired path real-tool-qualified;
- update status docs;
- do not change analytical-validation status.

## Task R4 — Wakhan version-contract RED

Establish which current package/runtime is reproducibly installable.

Acceptance:
- existing `0.5.0` policy cannot silently run another version;
- selected new version, if any, must be documented and test-pinned;
- command semantics must match current upstream usage;
- optional Severus parent remains exact-hash bound.

## Task R5 — Wakhan real-tool smoke

Only after R4:
- generated synthetic inputs;
- phased VCF and indexed BAM;
- optional breakpoints only if reproducibly generated;
- preserve real_tool_qualified=false until a successful same-head binary smoke.

## Task R6 — final verification

Run:
- repository safety
- schema freshness
- version agreement
- Ruff
- mypy
- full tests
- Severus real-tool smoke
- Wakhan real-tool smoke if qualified
- affected desktop/runtime smokes

No merge and no package version bump.
