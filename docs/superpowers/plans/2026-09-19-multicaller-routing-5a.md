# Supplemental declarative routing contract

Merge scope (2026-09-20): this is a separate, planning-only contract. It does not replace multicaller_contracts.MultiCallerPlan, the canonical lane planner, caller runtimes, or their qualification registry. Its existing_adapter flags describe this supplemental route, not platform-wide installation status. All its lanes remain NOT_RUN and execution_enabled=false. No implicit conversion into the canonical execution plan is provided.

# Multi-caller routing 5A implementation plan

> For agentic workers: execute task-by-task with test-driven development and verify
> the unchanged baseline plus the final branch. Native execution is out of scope.

**Goal:** A metadata-only planner for eight complementary CNV/SV callers that does not
misrepresent catalog entries as executed/validated integrations.
**Architecture:** Immutable catalog + typed planning contracts + deterministic pure
routing + versioned schemas. Reuse the existing pipeline after future adapter gates,
not by mutating its provider table in this block.
**Tech Stack:** Python 3.11+, Pydantic v2, pytest/unittest, canonical JSON/SHA-256.
**Spec:** `docs/superpowers/specs/2026-09-19-multicaller-routing-5a-design.md`.

## Global constraints

Research Use Only. No patient/genomic payloads. No caller executions or data downloads.
No merge/version bump. No numerical scientific defaults. No caller voting.
Preserve `feat/cnv-validation-program` and the canonical pipeline.

## Review focus

Unknown/NaN coverage must not be sufficient. The same sample cannot be its own normal.
Artifact reference/sample mismatches must fail. Disabled callers must stay visible.
An unimplemented native adapter must not become runnable because inputs are complete.

## Task 1 — catalog

Files: `src/ontseq_platform/multicaller_catalog.py`, `tests/test_multicaller_plan.py`.
- [x] Add `test_catalog_has_eight_distinct_callers` and run it RED.
- [x] Implement an immutable catalog and explicit legacy-adapter / planned distinction.
- [x] Assert exact provider identities, capabilities, and distinct signal families.
- [x] Run the test GREEN.

## Task 2 — request and routing contracts

Files: `src/ontseq_platform/multicaller_plan.py`, `tests/test_multicaller_plan.py`.
- [x] Add failing tests for each caller's missing prerequisites, assay restrictions,
  identity checks, explicit paired/unpaired intent and minimum-depth policy.
- [x] Implement `plan_multicaller(request: MultiCallerRequest) -> MultiCallerPlan`.
- [x] Keep every catalog lane visible. All lanes remain `execution_status=NOT_RUN`.
- [x] Ensure adapter availability remains explicit independently of input readiness.
- [x] Run all routing tests GREEN, then rejection and mutation tests.

## Task 3 — canonical lock and schema

Files: `src/ontseq_platform/multicaller_plan.py`, the same test file,
`schemas/multicaller-request.schema.json`, `schemas/multicaller-routing-plan.schema.json`,
`scripts/export_schemas.py`.
- [x] Add failing tests for request-order invariance, fingerprint/threshold sensitivity,
  a tampered plan, unchecked model-copy mutation and schema freshness.
- [x] Seal plans with canonical SHA-256; prohibit timestamps and execution outcomes.
- [x] Add the two models to the existing schema exporter without rewriting old schemas.
- [x] Run local focused tests (47 passed, 27 subtests passed).
- [ ] Verify repository CI at the final PR head and record exact run results in the PR.

## Completion evidence

Record actual commands, RED failures, GREEN results, source hashes and unavailable
checks in the PR. Never report full-suite or live-caller execution based on the local
focused test workspace. The remaining native-adapter blocks are not completed by 5A.
