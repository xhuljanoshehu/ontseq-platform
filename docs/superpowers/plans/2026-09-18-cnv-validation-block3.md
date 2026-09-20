# CNV Validation Block 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the prospective CNV aggregation layer that turns the locked registration plus Block-2 evidence into deterministic event assignments, overall/stratified metrics, explicit missing/not-evaluable states, quantitative-error summaries, reproducibility summaries, acceptance decisions, and a sealed RUO validation report.

**Architecture:** Add a dedicated `cnv_validation_metrics.py` layer above `benchmark.py`, `cnv_validation_registration.py`, and `cnv_validation_evidence.py`. Matching is performed once per observed caller lane using the existing deterministic maximum-cardinality comparator; aggregation then indexes matched/unmatched truth and query events into prospectively registered strata without re-matching inside strata. Truth-indexed TP/FN and query-indexed TP/FP are kept distinct so event-size stratification cannot change matching merely because a called segment crosses a size boundary. Specificity remains fail-closed and requires explicit negative-unit assessment tied to the registered negative universe.

**Tech Stack:** Python 3.11+, Pydantic v2, existing ONTSeq `benchmark_case()`, canonical SHA-256 content locks, unittest/pytest, JSON Schema exporter, GitHub Actions CI.

**Spec:** `docs/superpowers/specs/2026-09-17-cnv-validation-program-design.md`

## Global Constraints

- Research Use Only; no analytical or clinical validity claim.
- No real patient/genomic payloads in Git; all tests synthetic.
- Reuse the existing deterministic CNV benchmark matcher; do not invent a second event matcher.
- Never infer specificity from event-level FP counts alone.
- Continuous coverage/tumour-fraction values remain in Block-2 evidence; Block 3 derives labels only.
- Matching thresholds come only from the preregistered matrix.
- No post-hoc cutpoint tuning.
- No automatic production caller selection.
- No QDNAseq/ACE, ichorCNA or Spectre execution changes in this block.
- No merge and no package version bump; keep `0.8.2`.
- PR #79 remains Draft/unmerged.
- Every aggregate must remain traceable to registration truth IDs, normalized-event IDs, full-evidence IDs and the Block-2 manifest lock.

---

## File map

- Create `src/ontseq_platform/cnv_validation_metrics.py` — deterministic lane assignment, strata, metrics, acceptance and sealed report.
- Create `tests/test_cnv_validation_metrics.py` — RED/GREEN contract tests using synthetic registration/evidence only.
- Modify `src/ontseq_platform/cnv_validation_contracts.py` — optional quantitative orthogonal truth needed for cellularity/ploidy error.
- Modify `tests/test_cnv_validation_contracts.py` — quantitative-truth linkage tests.
- Modify `scripts/export_schemas.py` — export Block-3 report schema and refreshed cohort schema.
- Create `schemas/cnv-validation-report.schema.json` — generated schema.
- Refresh `schemas/cnv-validation-cohort.schema.json` if quantitative truth extends the cohort contract.
- Create `docs/CNV_VALIDATION_BLOCK3_STATUS.md` — RUO scope and denominator semantics.
- Modify PR #79 description only after final Same-Head verification; this does not mutate repository HEAD.

### Task 1: Quantitative orthogonal-truth contract

**Interfaces**
- Produces `CnvQuantitativeTruth` and optional `quantitative_truth` on `CnvValidationSpecimen`.
- Later metrics consume `cellularity` and `ploidy` only when explicitly linked to a registered orthogonal truth source.

- [ ] **Step 1: Write failing tests** in `tests/test_cnv_validation_contracts.py`:
  - quantitative truth accepts cellularity in `[0,1]` and ploidy `>0`;
  - at least one quantitative value is required;
  - `truth_source_resource_id` must reference one of the specimen's `truth_sources`;
  - absent quantitative truth remains valid and is never imputed.

- [ ] **Step 2: Run the focused tests** and confirm RED because `CnvQuantitativeTruth`/field do not exist.

- [ ] **Step 3: Implement minimally** in `cnv_validation_contracts.py`:
  ```python
  class CnvQuantitativeTruth(StrictModel):
      truth_source_resource_id: str
      cellularity: float | None = Field(default=None, ge=0, le=1)
      ploidy: float | None = Field(default=None, gt=0)
      note: str = Field(min_length=3)

      @model_validator(mode="after")
      def at_least_one_value(self) -> "CnvQuantitativeTruth": ...
  ```
  Add `quantitative_truth: CnvQuantitativeTruth | None = None` to `CnvValidationSpecimen` and validate source-ID linkage.

- [ ] **Step 4: Run focused + existing contract/registration tests**; confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): add quantitative orthogonal truth contract`.

### Task 2: Lane identity and one-time event assignment

**Interfaces**
- Produces `CnvEvaluationLaneKey`, `CnvMatchedEventAssignment`, `CnvLaneEventAssignment`, and `assign_cnv_validation_lanes(registration, evidence)`.
- Uses `benchmark_case()` exactly once per observed lane.
- A lane is identified by specimen, caller/version/runtime policy, genome build, data basis, reference/input, coverage, tumour fraction, bin size and replicate context.

- [ ] **Step 1: Write failing tests** in `tests/test_cnv_validation_metrics.py`:
  - evidence manifest registration SHA must equal `registration.lock_sha256`;
  - caller/version/adapter/runtime identity must match the registered caller lock;
  - every analytical lane requires exactly one `RUN_SUMMARY` record;
  - normalized events must belong to exactly one lane;
  - `OBSERVED` lane reuses preregistered matching thresholds and produces deterministic matched/unmatched assignments;
  - `FAILED`, `NO_CALL`, `NOT_ASSESSABLE` lanes do not fabricate event assignments;
  - biological-negative truth evidence is not silently treated as a successful caller execution lane;
  - query `GenomicEvent` IDs are normalized-event IDs and preserve event type/locus/copy number.

- [ ] **Step 2: Verify RED** because the metrics module does not exist.

- [ ] **Step 3: Implement minimal lane models + assignment function**:
  - convert normalized events to benchmark query events;
  - call existing `benchmark_case()`;
  - store deterministic match score/reciprocal overlap and unmatched IDs;
  - strip runtime timestamps/warnings from the locked assignment object;
  - never re-match after stratum classification.

- [ ] **Step 4: Run focused tests + `tests/test_benchmark.py`** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): add locked lane event assignment`.

### Task 3: Prospective strata and denominator-safe event metrics

**Interfaces**
- Produces deterministic band functions and `CnvMetricResult`.
- Event metric aggregation keeps separate `tp_truth` and `tp_query` counts:
  - sensitivity = `tp_truth / (tp_truth + fn)`;
  - precision = `tp_query / (tp_query + fp)`;
  - F1 = harmonic mean only when both are defined;
  - FP burden = unmatched query events per observed lane in the same technical stratum.
- Whole-chromosome gain/loss use `whole_chromosome` size band; deletion/duplication use preregistered subchromosomal size bands.

- [ ] **Step 1: Write failing tests**:
  - coverage/tumour-fraction/event-size labels are derived only from registered cutpoints;
  - exact-boundary behavior is deterministic (`[lower, upper)`);
  - continuous evidence values remain unchanged;
  - a matched pair crossing an event-size boundary remains matched once; it contributes to the truth-size TP for recall and query-size TP for precision rather than being re-matched;
  - overall and fully stratified views are both emitted;
  - registered strata with denominator zero appear explicitly as `NOT_EVALUABLE`;
  - denominator counts and contributing truth/query/normalized/full-evidence IDs are present;
  - no undefined ratio is converted to 0 or 1.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement**:
  - canonical numeric-band labels using registered cutpoints;
  - registered cross-product keys for caller × build × data_basis × coverage × tumour_fraction × event_class × event_size × bin_size plus an overall event view;
  - truth-side and query-side accumulators derived from the already locked lane assignments;
  - `CnvMetricResult` with explicit `value`, `evaluable_denominator`, component counts, state/reason and provenance IDs.

- [ ] **Step 4: Run focused tests** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): aggregate prospective event metrics`.

### Task 4: Run-state rates and fail-closed specificity

**Interfaces**
- Produces `CnvNegativeUnitAssessment` plus no-call/failure/not-assessable/specificity metrics.
- Missing execution remains distinct from no-call/failure and is visible in the report.

- [ ] **Step 1: Write failing tests**:
  - no-call, failure and not-assessable are separate counts/rates;
  - rate denominator is explicit executed lane count, never silently the matrix cross-product;
  - registered strata with no executed lane remain `NOT_EVALUABLE`, not 0%;
  - specificity with no negative universe is `NOT_EVALUABLE`;
  - negative universe alone is still insufficient: explicit negative-unit assessment is required;
  - assessment must bind specimen/lane, universe ID, universe resource SHA and assessability-mask SHA;
  - `false_positive_units <= assessed_units`;
  - assessed units must equal the registered universe for the lane in this first implementation;
  - specificity = `(assessed_units - false_positive_units) / assessed_units`;
  - event-class/size-specific specificity is `NOT_APPLICABLE` unless a future registered unit classification exists.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement minimal run-state and specificity aggregation** without deriving negative units from event FP counts.

- [ ] **Step 4: Run focused tests** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): add run-state and specificity metrics`.

### Task 5: Quantitative copy-number/cellularity/ploidy error and reproducibility

**Interfaces**
- Copy-number error uses only matched event pairs where both truth and normalized copy number exist.
- Cellularity/ploidy error uses only a single selected caller fit per lane plus explicit orthogonal quantitative truth.
- Reproducibility compares repeat-group lanes within the same caller/build/data-basis/bin context using the existing matcher; pairwise F1 is treated as engineering concordance, not truth.

- [ ] **Step 1: Write failing tests**:
  - copy-number mean absolute error and signed bias use only eligible matched pairs;
  - missing quantitative copy number is excluded with explicit denominator reduction;
  - cellularity/ploidy remain `NOT_EVALUABLE` without orthogonal quantitative truth;
  - multiple selected fits in one lane fail closed;
  - one selected fit yields absolute error against explicit truth;
  - alternative fits stay in full evidence and are not silently promoted;
  - repeat groups spanning different biological specimen IDs fail closed;
  - within-run/between-run pairwise event concordance uses same registered matching thresholds;
  - two empty repeat call sets remain `NOT_EVALUABLE`, not artificial 100% reproducibility;
  - report retains pair IDs and matched normalized-event IDs.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement quantitative-error summaries and `CnvReproducibilityPair` / mean pairwise F1 metric.**

- [ ] **Step 4: Run focused tests** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): add quantitative and reproducibility metrics`.

### Task 6: Acceptance questions, sealed report and metric provenance

**Interfaces**
- Produces `CnvAcceptanceResult`, `CnvValidationMetricReport`, and `aggregate_cnv_validation(...)`.
- Acceptance decisions: `PASS`, `FAIL`, `NO_CALL`, `NOT_EVALUABLE`; optional `NOT_APPLICABLE` is allowed for metrics such as event-specific specificity.
- No acceptance value is invented when denominator is below the preregistered minimum.

- [ ] **Step 1: Write failing tests**:
  - question filters match canonical stratum keys exactly;
  - below `minimum_evaluable_denominator` => `NOT_EVALUABLE`;
  - lower/upper bounds yield PASS/FAIL deterministically;
  - missing matching stratum => explicit `NOT_EVALUABLE`;
  - every metric carries mathematical provenance IDs and Block-2 manifest SHA;
  - report binds registration lock + evidence manifest SHA + matching thresholds;
  - report SHA changes on any metric/assignment/acceptance change;
  - report contains no production-caller winner field;
  - no timestamp-dependent data enters the report lock.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement the orchestrator** that verifies Block-2 traceability first, builds lane assignments once, aggregates all metric families, evaluates questions and seals canonical report JSON.

- [ ] **Step 4: Run all Block-3 focused tests + Block-1/2 CNV tests + benchmark/dilution regression tests**.

- [ ] **Step 5: Commit** `feat(cnv): seal validation metric report`.

### Task 7: Schema, RUO documentation, diff review and final Same-Head gate

- [ ] **Step 1: Add schema-export entries** for `CnvValidationMetricReport` and refreshed `CnvValidationCohort`.

- [ ] **Step 2: Verify RED** with `python scripts/export_schemas.py --check` because generated files are absent/stale.

- [ ] **Step 3: Generate schemas reproducibly**, commit only generated schema outputs, and remove any temporary generation workflow.

- [ ] **Step 4: Add `docs/CNV_VALIDATION_BLOCK3_STATUS.md`** documenting:
  - one-time matching before stratification;
  - truth-indexed vs query-indexed TP semantics;
  - explicit denominator rules;
  - specificity negative-unit requirement;
  - missing execution vs no-call/failure;
  - quantitative/reproducibility limitations;
  - no caller selection/clinical claim.

- [ ] **Step 5: Full PR scope review**:
  - no real biological payloads;
  - no caller execution changes;
  - no production default/winner;
  - no tolerance/cutpoint mutation;
  - no temporary workflow remains;
  - package version still `0.8.2`;
  - PR remains Draft/unmerged.

- [ ] **Step 6: Fresh Same-Head verification** on one immutable final commit:
  - main CI Python 3.11/3.12/3.13 including Ruff, mypy, schema freshness, all regression tests, Snakemake DAGs and synthetic report;
  - local-real-tool-smoke;
  - MARLIN runtime smoke;
  - Desktop bundle contract;
  - Desktop CI including Linux relocation/system smoke and Windows/WPF bundle verification.

- [ ] **Step 7: Record final run IDs and Block-3 evidence in PR #79 conversation/body without mutating repository HEAD.**
