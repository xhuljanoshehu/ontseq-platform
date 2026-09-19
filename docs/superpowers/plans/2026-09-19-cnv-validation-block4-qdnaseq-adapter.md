# CNV Validation Block 4 — QDNAseq+ACE Full-Evidence Adapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Translate one already-completed QDNAseq+ACE multi-resolution result bundle into the locked CNV validation Full-Evidence/Normalized-Event contract without re-running the caller, dropping values, or privileging 500 kbp during the 100/500/1000-kbp validation comparison.

**Architecture:** Add a validation-only adapter beside the existing live QDNAseq caller. The adapter consumes an immutable `QDNAseqCallReport`, the promoted native result directory, the exact `QDNAseqPolicy`, a locked run-identity envelope and the CNV validation registration. It fingerprints every declared caller artifact, emits row-level evidence for all configured resolutions and all ACE fit candidates, derives comparable normalized events independently for every registered bin size by reusing the caller's existing normalization routine, checks the 500-kbp derivation against the caller report, and seals one Block-2 `CnvValidationEvidenceManifest`. It does not execute R/QDNAseq/ACE and does not wire a caller into the production pipeline.

**Tech Stack:** Python 3.11+, Pydantic v2, existing `ontseq_platform.cnv.qdnaseq` parser/normalizer, existing Block-1/2/3 CNV validation contracts, SHA-256 native-artifact fingerprints, unittest/pytest, JSON Schema exporter, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-09-17-cnv-validation-program-design.md`

## Global Constraints

- Research Use Only; no analytical or clinical validity claim.
- No real patient/genomic payloads in Git; tests use synthetic QDNAseq bundles.
- Do not change `scripts/run_qdnaseq_ace.R`, `run_qdnaseq_ace()`, production CNV stage wiring or caller-selection behavior in this block.
- Preserve every QDNAseq/ACE output declared by `QDNAseqCallReport.output_files` as a fingerprinted native artifact.
- Preserve 100, 500 and 1000-kbp selected fits, alternative ACE fits, bin rows, segment rows and chromosome summaries.
- Do not delete or overwrite native values because they are not used by a primary metric.
- Reuse the existing versioned QDNAseq segment-to-event normalization logic; do not create a second biological calling rule.
- Every configured bin size is an independently evaluable validation lane; production/reporting `primary_bin_size_kbp=500` must not suppress 100/1000-kbp validation events.
- A complete successful QDNAseq bundle with zero events is an `OBSERVED` empty call set for validation, not a technical `NO_CALL`.
- Caller/runtime/input/reference identity must be verified against the preregistration before any manifest is sealed.
- No automatic winning/production caller decision.
- No merge and no package version bump; keep `0.8.2`.
- Base this stacked work on verified Block-3 head `cf6a3e1487c8ad7c6ad061e941306b66f48d3a56`.
- PR #79 remains frozen, Draft and unmerged.

## Review Focus

1. **Successful zero-event caller output** — must remain an observed empty query set so truth events become FNs instead of silently disappearing from sensitivity.
2. **Primary-vs-validation semantics across 100/500/1000 kbp** — every registered resolution must generate independently evaluable primary validation events while the production 500-kbp preference remains metadata only.
3. **Alternative ACE fit duplication** — the selected fit must appear exactly once even when the R summary's alternatives list contains the selected candidate as its first row.
4. **Native artifact completeness/path safety** — traversal, absolute paths, missing files, duplicate declared paths or unreferenced declared artifacts must fail closed.
5. **Secondary-resolution normalization drift** — the same existing QDNAseq normalization routine and policy thresholds must be used at every bin size; the re-derived primary-bin events must exactly agree with `QDNAseqCallReport.events`.

---

## File map

- Create `src/ontseq_platform/cnv/qdnaseq_validation_adapter.py` — validation-only adapter policy, run identity, native-artifact fingerprinting, evidence extraction, all-bin event derivation and manifest sealing.
- Create `tests/test_qdnaseq_validation_adapter.py` — synthetic multi-resolution adapter contract tests.
- Modify `scripts/export_schemas.py` — export the versioned QDNAseq validation-adapter policy schema.
- Create `schemas/qdnaseq-validation-adapter-policy.schema.json` — generated schema.
- Create `docs/CNV_VALIDATION_BLOCK4_STATUS.md` — evidence mapping and explicit non-claims.
- Create stacked Draft PR from `feat/cnv-qdnaseq-evidence-adapter-v1` to `feat/cnv-validation-program` only after implementation/final verification.

### Task 1: Validation adapter policy and locked run identity

**Files:**
- Create: `src/ontseq_platform/cnv/qdnaseq_validation_adapter.py`
- Create: `tests/test_qdnaseq_validation_adapter.py`

**Interfaces:**
- Produces `QDNAseqValidationAdapterPolicy`, `QDNAseqValidationRunIdentity`, `canonical_qdnaseq_ace_version(report)`, and `qdnaseq_adapter_policy_sha256(policy)`.
- Later tasks consume these types before reading any result file.

- [ ] **Step 1: Write failing tests** that require:
  - adapter policy embeds the complete `QDNAseqPolicy`;
  - policy hard-locks `caller_id="qdnaseq_ace"`, all-artifact retention, all-bin retention, all-segment retention, all-fit-candidate retention and normalization ID `qdnaseq-validation-normalization-v1`;
  - canonical caller version is built from exactly one QDNAseq and one ACE `ToolRecord`, e.g. `QDNAseq=1.42.0;ACE=1.24.0`;
  - missing, duplicate or `UNKNOWN` QDNAseq/ACE versions fail closed;
  - adapter-policy SHA is canonical and changes if any embedded QDNAseq parameter or mapping rule changes;
  - observed execution/input/reference identity must exactly equal the registered caller/specimen lock;
  - sample ID/genome build in `QDNAseqCallReport` must match the registered specimen;
  - registered `adapter_policy_sha256` must equal the actual adapter policy hash.

- [ ] **Step 2: Run focused tests** and confirm RED because the adapter module does not exist.

- [ ] **Step 3: Implement minimal strict models and identity verification**. The adapter policy must carry the exact QDNAseq policy payload rather than merely copying the registration hash.

- [ ] **Step 4: Run focused tests + Block-1 registration tests** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): lock qdnaseq validation adapter identity`.

### Task 2: Complete native-artifact registry

**Files:**
- Modify: `src/ontseq_platform/cnv/qdnaseq_validation_adapter.py`
- Modify: `tests/test_qdnaseq_validation_adapter.py`

**Interfaces:**
- Produces `fingerprint_qdnaseq_artifacts(report, result_dir) -> list[CnvNativeArtifactReference]`.
- Every artifact ID is deterministic from semantic role + bin size/path; stored paths are portable relative paths only.

- [ ] **Step 1: Write failing tests**:
  - every `QDNAseqCallReport.output_files` entry becomes exactly one native artifact;
  - summary JSON, cross-bin chromosome consensus, per-bin segments/chromosomes/bins/models/plots/RDS are all retained;
  - artifact SHA-256 and size match the actual promoted files;
  - absolute paths, `..`, Windows drive paths and backslashes fail closed;
  - duplicate declared paths fail closed;
  - missing/empty declared files fail closed;
  - output files not rooted directly in the supplied promoted result directory fail closed;
  - native artifact IDs are deterministic across repeated adapter runs.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement artifact-role inference and fingerprinting** without reading or rewriting file payloads.

- [ ] **Step 4: Run focused tests** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): fingerprint complete qdnaseq result bundle`.

### Task 3: Preserve selected and alternative ACE fits at every resolution

**Files:**
- Modify: `src/ontseq_platform/cnv/qdnaseq_validation_adapter.py`
- Modify: `tests/test_qdnaseq_validation_adapter.py`

**Interfaces:**
- Produces `qdnaseq_fit_evidence(...)->list[CnvFullEvidenceRecord]`.
- Emits one selected `CALLER_FIT` per configured bin size plus every distinct alternative candidate.

- [ ] **Step 1: Write failing tests**:
  - exactly one selected fit exists for each of 100/500/1000 kbp;
  - selected fit retains cellularity, ploidy, fit error, candidate count and segment count as typed fields/measurements;
  - alternative candidates are retained with `selected_fit=False` and `SECONDARY`;
  - an alternative numerically identical to the selected fit is not duplicated;
  - distinct alternatives with equal error but different cellularity/ploidy are retained separately;
  - selected and alternative fit records link the model/summary artifacts that contain their source values;
  - all records carry the exact caller parameters, caller-parameter SHA, dependency versions and registered identity;
  - a report missing a registered bin size or containing an unregistered extra bin size fails closed.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement fit evidence** using deterministic `fit_group_id` and record IDs containing the bin resolution and candidate ordinal.

- [ ] **Step 4: Run focused tests + existing QDNAseq cellularity tests** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): retain qdnaseq ace fit candidates`.

### Task 4: Row-level bin, segment and chromosome evidence

**Files:**
- Modify: `src/ontseq_platform/cnv/qdnaseq_validation_adapter.py`
- Modify: `tests/test_qdnaseq_validation_adapter.py`

**Interfaces:**
- Produces per-bin `RUN_SUMMARY`, `CALLER_BIN`, `CALLER_SEGMENT` and `CHROMOSOME_SUMMARY` records.
- Native TSV remains canonical for every raw column; the evidence row carries exact locus, finite numeric measurements and a fingerprint link back to that TSV.

- [ ] **Step 1: Write failing tests**:
  - one observed run summary is emitted per registered bin size;
  - all rows from all 100/500/1000-kbp bin TSVs are retained, including `use=FALSE` rows;
  - `use=FALSE` rows become `EXCLUDED_FROM_PRIMARY_METRIC` with explicit reason, not deleted;
  - every segment row is retained, including neutral `call=0` segments;
  - event-producing segment rows are `USED_FOR_PRIMARY_ANALYSIS` independently within each validation lane;
  - neutral/below-threshold segments remain retained with explicit non-contributing status/reason;
  - chromosome summary rows from each resolution are retained;
  - cross-bin consensus values are retained on the primary run-summary as finite typed measurements and the native consensus TSV is linked;
  - zero-based half-open coordinate tags, non-negative starts, end>start, finite numeric fields and unique native row IDs fail closed;
  - `-Inf` qnorm sentinels remain visible through the native artifact but are never converted into an invalid infinite numeric measurement;
  - no unlisted caller-native artifact becomes an orphan in the eventual manifest.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement strict TSV row extraction**. Do not mutate TSVs. Store every finite numeric field exposed by the row in `CnvNumericMeasurement`; rely on the linked fingerprinted TSV for non-numeric native columns not represented by the generic evidence model.

- [ ] **Step 4: Run focused tests + existing QDNAseq coordinate/runtime tests** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): retain qdnaseq multibin row evidence`.

### Task 5: Comparable normalized CNV events for 100/500/1000 kbp

**Files:**
- Modify: `src/ontseq_platform/cnv/qdnaseq_validation_adapter.py`
- Modify: `tests/test_qdnaseq_validation_adapter.py`

**Interfaces:**
- Produces one set of `CnvNormalizedEventRecord` objects per registered bin-size lane.
- Reuses `_events_from_primary_segments()` and `_assessable_bin_extents()` from the existing QDNAseq implementation with the same QDNAseq policy, reference lock and tool records.

- [ ] **Step 1: Write failing tests**:
  - normalized events are generated independently for 100, 500 and 1000 kbp;
  - every bin-size lane uses the same registered `minimum_segment_bins`, whole-chromosome fraction/span basis and ACE ploidy semantics;
  - normalized IDs include bin size and are unique across resolutions;
  - each normalized event links the exact segment full-evidence row with matching locus/copy number;
  - no event can be sourced only from a fit/bin/run-summary record;
  - all three resolution event sets use `USED_FOR_PRIMARY_ANALYSIS` for the validation comparison, regardless of production `primary_bin_size_kbp=500`;
  - the re-derived 500-kbp event set must semantically equal `QDNAseqCallReport.events` (event type, locus, copy number and whole-chromosome span evidence); any drift fails closed;
  - ambiguous/missing segment-to-event source mapping fails closed;
  - a secondary-resolution event never modifies the caller's production `QDNAseqCallReport`.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement all-bin event derivation and exact-source mapping**. This is normalization only; do not invoke QDNAseq or ACE.

- [ ] **Step 4: Run focused tests + existing QDNAseq event/coordinate tests + Block-3 metric tests** and confirm GREEN.

- [ ] **Step 5: Commit** `feat(cnv): normalize qdnaseq events across bin sizes`.

### Task 6: Successful empty-call semantics and sealed evidence manifest

**Files:**
- Modify: `src/ontseq_platform/cnv/qdnaseq_validation_adapter.py`
- Modify: `tests/test_qdnaseq_validation_adapter.py`

**Interfaces:**
- Produces `adapt_qdnaseq_validation_evidence(...) -> CnvValidationEvidenceManifest`.

- [ ] **Step 1: Write failing tests**:
  - a complete `QDNAseqCallReport.status=NO_CALL` with zero events still yields `CnvRunOutcomeState.OBSERVED` per bin and zero normalized events;
  - when compared with registered positive truth, Block-3 assignment therefore produces FNs rather than excluding the lane as a technical no-call;
  - caller process failure is not fabricated from a missing report; the adapter requires a complete parseable report and fails closed;
  - every native artifact is referenced by at least one full-evidence record;
  - every normalized event source exists and matches all identity fields;
  - the final manifest has `retain_all_evidence=true` and a valid canonical manifest SHA;
  - repeated adaptation of identical inputs is byte-/hash-deterministic;
  - changing any native artifact byte changes the manifest lock;
  - changing 100/500/1000 fit/segment/bin evidence cannot silently leave the old manifest SHA valid.

- [ ] **Step 2: Verify RED**.

- [ ] **Step 3: Implement the orchestration function** that verifies identity first, fingerprints artifacts, emits all evidence classes, normalizes all resolution events and calls `seal_cnv_validation_evidence()`.

- [ ] **Step 4: Run all adapter tests + Block-1/2/3 CNV tests + QDNAseq runtime/coordinate/cellularity tests**.

- [ ] **Step 5: Commit** `feat(cnv): seal qdnaseq multibin validation evidence`.

### Task 7: Public schema, RUO documentation, stacked PR and final same-head gate

**Files:**
- Modify: `scripts/export_schemas.py`
- Create: `schemas/qdnaseq-validation-adapter-policy.schema.json`
- Create: `docs/CNV_VALIDATION_BLOCK4_STATUS.md`

- [ ] **Step 1: Export `QDNAseqValidationAdapterPolicy` schema** and verify `python scripts/export_schemas.py --check` is RED before the generated file exists.

- [ ] **Step 2: Generate the schema reproducibly**, commit only the expected schema output, and remove any temporary generator workflow.

- [ ] **Step 3: Document RUO boundaries**:
  - adapter consumes existing caller output and never executes QDNAseq/ACE;
  - all 100/500/1000 native evidence and fit alternatives remain retained;
  - successful empty calls are observed empty sets;
  - all resolutions are validation-evaluable but no production bin-size/caller winner is chosen;
  - native TSV/RDS/plots stay sensitive external artifacts referenced by fingerprint;
  - no analytical/clinical claim.

- [ ] **Step 4: Whole-branch review**:
  - no changes to `scripts/run_qdnaseq_ace.R` or `run_qdnaseq_ace()`;
  - no production pipeline wiring/default changes;
  - no patient data;
  - no post-hoc thresholds;
  - no temporary workflow;
  - version still `0.8.2`;
  - base remains verified Block-3 branch.

- [ ] **Step 5: Final immutable-head verification**:
  - main CI Python 3.11/3.12/3.13: schema freshness, Ruff, mypy, full tests, DAGs, synthetic report;
  - local-real-tool-smoke;
  - MARLIN runtime smoke;
  - Desktop bundle contract;
  - Desktop CI;
  - focused QDNAseq runtime/coordinate/cellularity/adapter tests.

- [ ] **Step 6: Create a stacked Draft PR** with base `feat/cnv-validation-program`. Record final head/run IDs and keep it unmerged pending explicit approval.
