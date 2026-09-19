# Multi-Caller Integration Implementation Plan

> **Spec:** `docs/superpowers/specs/2026-09-19-multicaller-integration-design.md`
> **Branch:** `feat/multicaller-integration`
> **Base:** `b6e0338b6f042239a1399e3438e6c2bab4e48b08`
> **Execution:** inline, strict RED → GREEN, synthetic fixtures only
> **Release boundary:** Research Use Only; Human Review Required; no merge/version bump

## Global constraints

- Do not modify the existing one-provider-per-stage `RunComponents` behavior.
- Do not add patient data or raw genomic payloads to Git.
- Do not infer biological negative from ineligible, unavailable or failed callers.
- Do not implement caller voting or winner selection.
- Preserve caller-native artifacts and alternative solutions.
- Keep provider mode semantics explicit; no silent fallback between paired and tumor-only modes.
- A new runtime is not "integrated" until a synthetic real-tool smoke test exists and passes.
- Run repository safety, version, lint, type, schema and regression checks on every completed task head.

## File map

### Shared subsystem

- `src/ontseq_platform/multicaller_contracts.py`
  - enums and typed models for provider, mode, analytical domain, input roles, lane request, lane plan and plan lock.
  - deterministic canonical SHA-256.
  - static caller catalog for all selected families and modes.
- `src/ontseq_platform/multicaller_planner.py`
  - eligibility evaluation.
  - runtime-availability handling.
  - dependency graph validation and topological ordering.
  - deterministic plan sealing.
- `schemas/multicaller-plan.schema.json`
  - exported public schema once contracts stabilize.
- `tests/test_multicaller_contracts.py`
  - caller catalog and model invariants.
- `tests/test_multicaller_planner.py`
  - fail-closed eligibility, runtime state and dependency graph tests.

### Existing caller bridges

- `src/ontseq_platform/multicaller_existing.py`
  - bridge descriptors for existing QDNAseq+ACE, Sniffles2 and cuteSV adapters.
  - no change to caller execution or normalized biological output.
- `tests/test_multicaller_existing.py`

### New runtime/adapters

- `src/ontseq_platform/cnv/spectre.py`
- `tests/test_spectre_runtime.py`
- `configs/cnv/spectre.depth-only.technical.yaml`

- `src/ontseq_platform/tumor_inputs.py`
- `tests/test_tumor_inputs.py`

- `src/ontseq_platform/sv/severus.py`
- `tests/test_severus_runtime.py`
- `configs/sv/severus.paired.technical.yaml`

- `src/ontseq_platform/tumor/savana.py`
- `tests/test_savana_runtime.py`
- `configs/tumor/savana.paired.technical.yaml`
- `configs/tumor/savana.tumor-only.technical.yaml`

- `src/ontseq_platform/tumor/wakhan.py`
- `tests/test_wakhan_runtime.py`
- `configs/tumor/wakhan.phased-cna.technical.yaml`

- `src/ontseq_platform/cnv/ichorcna.py`
- `tests/test_ichorcna_runtime.py`
- `configs/cnv/ichorcna.ulp-wgs.technical.yaml`

### Comparison/report integration

- `src/ontseq_platform/multicaller_report.py`
- `tests/test_multicaller_report.py`
- `docs/MULTICALLER_INTEGRATION_STATUS.md`

## Task 1 — 5A catalog and core contracts

**Produces:** stable typed catalog and lane request/plan models.

1. Add RED tests asserting:
   - all selected provider/mode pairs exist;
   - unknown provider/mode is rejected;
   - `savana:paired` and `savana:tumor_only` are separate catalog entries;
   - `spectre:depth_only` declares no caller parent;
   - each catalog entry declares analytical domain, compatible data basis/assay categories and required input roles;
   - catalog IDs are unique and deterministic.
2. Run only `tests/test_multicaller_contracts.py`; expect import/model failures.
3. Implement minimal contracts/catalog.
4. Re-run focused tests GREEN.
5. Run Ruff + mypy for touched modules and full existing regression suite.
6. Commit.

## Task 2 — 5A fail-closed eligibility planner

**Produces:** eligibility decisions independent of execution.

1. RED tests:
   - SAVANA paired with tumor BAM but no matched normal => `INELIGIBLE`, explicit normal reason.
   - no implicit substitution to tumor-only.
   - missing runtime for otherwise eligible lane => `RUNTIME_UNAVAILABLE`, never `NO_CALL`.
   - Adaptive Sampling BAM alone cannot satisfy lcWGS-only ichorCNA lane.
   - missing required phasing input blocks Wakhan.
2. Implement planner with explicit reasons.
3. GREEN focused tests.
4. Full checks.
5. Commit.

## Task 3 — 5A dependency DAG and deterministic plan lock

**Produces:** sealed multi-lane plan.

1. RED tests:
   - dependency cycle rejected;
   - Spectre depth-only has no Sniffles dependency;
   - Spectre supported mode requires a declared Sniffles parent;
   - changing parent lane content changes dependent plan lock;
   - unrelated lane content remains independently addressable;
   - duplicate lane IDs rejected.
2. Implement DAG validation/topological order and canonical plan SHA-256.
3. GREEN + full checks.
4. Commit.

## Task 4 — collection-level identity namespace

**Produces:** collision-proof multi-specimen evidence addresses without rewriting existing single-sample IDs.

1. RED tests with two specimens/repeats carrying identical native record names.
2. Add collection-address helper based on lane lock + native evidence identity.
3. Verify existing Block-4 evidence model is unchanged.
4. GREEN + full checks.
5. Commit.

## Task 5 — 5C existing QDNAseq/Sniffles2/cuteSV bridge

**Produces:** multi-caller lanes can point to existing adapters without changing their output.

1. RED tests that bridge metadata resolves:
   - QDNAseq+ACE multibin;
   - Sniffles2 standard;
   - cuteSV standard.
2. Verify bridge never marks caller agreement as truth.
3. Add explicit `sniffles2:mosaic` catalog entry but keep runtime availability false until separately qualified.
4. GREEN + full regression.
5. Commit.

## Task 6 — 5B Spectre depth-only runtime

**Produces:** first new independent CNV runtime.

Precondition: verify current upstream primary documentation/version/CLI before pinning.

1. Add policy model and synthetic fixtures.
2. RED parser tests for native CNV output, empty calls, malformed coordinates/non-finite values and build mismatch.
3. Implement shell-free command-vector builder and parser.
4. RED runtime test proving no Sniffles input is passed for `depth_only`.
5. Add pinned technical environment/config.
6. Add synthetic real-tool smoke workflow/test if runtime can be installed reproducibly in CI.
7. GREEN + full checks.
8. Commit.

## Task 7 — 5D tumor auxiliary input contracts

**Produces:** reusable fingerprinted normal/SNP/BAF/phasing/PoN inputs.

1. RED tests for sample/build/reference mismatch and missing fingerprints.
2. Implement typed contracts.
3. Add parent/dependency binding fields used by downstream lanes.
4. GREEN + full checks.
5. Commit.

## Task 8 — 5E Severus paired adapter

**Produces:** conservative paired tumor-normal SV lane.

Precondition: verify current primary docs and executable CLI.

1. RED tests for missing normal, build mismatch and native graph/breakpoint retention.
2. Implement paired command builder and native artifact contract.
3. Keep single-sample mode catalogued but runtime-unqualified until separately tested.
4. Add synthetic real-tool smoke if reproducible.
5. GREEN + full checks.
6. Commit.

## Task 9 — 5F SAVANA paired and tumor-only adapters

**Produces:** explicit non-fallback SAVANA modes.

Precondition: verify current primary docs and CLI.

1. RED paired-mode missing-normal test is repeated at runtime boundary.
2. RED tests for distinct paired/tumor-only plan identity.
3. Preserve SV and CNA native outputs separately.
4. Preserve selected and alternative purity/ploidy solutions when emitted.
5. Add synthetic real-tool smoke for each qualified mode.
6. GREEN + full checks.
7. Commit.

## Task 10 — 5G Wakhan phased CNA adapter

**Produces:** allele-aware CNA lane with explicit upstream dependencies.

1. RED tests for missing phased/BAF input and dependency cycles.
2. Implement native sidecar/artifact preservation.
3. If Severus breakpoints are configured, bind exact parent lane lock.
4. Add synthetic real-tool smoke if reproducible.
5. GREEN + full checks.
6. Commit.

## Task 11 — 5H ichorCNA ULP-WGS lane

**Produces:** assay-gated ichorCNA integration.

1. RED tests that incompatible Adaptive Sampling/marrow transfer is ineligible by default.
2. Implement ULP-WGS/cfDNA-compatible policy and parser.
3. Any ONT transfer mode remains a different unvalidated catalog mode rather than silently reusing ULP-WGS eligibility.
4. Add synthetic real-tool smoke if reproducible.
5. GREEN + full checks.
6. Commit.

## Task 12 — 5I dependency-aware comparison report

**Produces:** human-reviewable comparison without voting.

1. RED tests:
   - dependent evidence shown as dependent, not independent confirmation;
   - failed/unavailable/ineligible lanes remain visible;
   - no generic majority/winner field exists;
   - native evidence links remain reachable;
   - conflicts are preserved, not averaged away.
2. Implement report contract and schema.
3. GREEN + full checks.
4. Commit.

## Task 13 — docs, schema export and final branch verification

1. Export `multicaller-plan.schema.json` and report schema if public.
2. Add status doc distinguishing:
   - catalogued;
   - planner-qualified;
   - adapter-qualified;
   - real-tool-qualified;
   - analytically validated.
3. Run:
   - `make safety`
   - `make versions`
   - `make lint`
   - `make test`
   - schema freshness
   - affected real-tool workflows
   - Desktop/runtime smoke when packaging changes.
4. Review branch against spec and document any unqualified runtimes explicitly.
5. Keep PR draft; no merge/version bump.
