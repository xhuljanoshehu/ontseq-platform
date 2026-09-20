# CNV Validation Block 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement an evidence-complete CNV validation layer that preserves caller-native records and artefacts, creates auditable normalised-event links without overwriting source evidence, and provides a fail-closed aggregate-to-source traceability contract for later metric aggregation.

**Architecture:** Add a dedicated `cnv_validation_evidence.py` module above the Block-1 matrix/cohort/registration contracts. Keep caller-native full evidence, normalised comparison evidence, and aggregate traceability as separate immutable layers. Seal the evidence manifest and traceability index by canonical SHA-256; validate every cross-reference fail-closed.

**Tech Stack:** Python 3.11+, Pydantic v2, existing `StrictModel`/`GenomeBuild`/`GenomicEvent`, canonical JSON SHA-256, unittest/pytest, existing schema exporter and GitHub CI.

**Spec:** `docs/superpowers/specs/2026-09-17-cnv-validation-program-design.md`

## Global Constraints

- Research Use Only; no analytical/clinical validity claim.
- No patient genomic payloads, BAM/POD5/FASTQ/VCF content, direct identifiers, or institutionally restricted truth payloads in Git.
- Repository tests use synthetic fixtures only.
- PR #79 remains Draft and unmerged; package version remains `0.8.2`.
- No caller execution changes in Block 2.
- No aggregate metric calculation in Block 2; only the traceability contract needed by Block 3.
- Prospective stratification and later aggregation must never delete, overwrite, or hide retained source evidence.
- Caller concordance is evidence, never truth.
- lcWGS, Adaptive-Sampling off-target, and Adaptive-Sampling on-target remain distinct data bases.

---

## File map

- Create `src/ontseq_platform/cnv_validation_evidence.py`: evidence enums/models, artifact references, full evidence records, normalised event records, evidence manifest sealing/validation, traceability index sealing/verification.
- Create `tests/test_cnv_validation_evidence.py`: synthetic RED/GREEN contract tests.
- Modify `scripts/export_schemas.py`: export evidence-manifest and traceability-index schemas.
- Create `schemas/cnv-validation-evidence.schema.json`: generated from Pydantic.
- Create `schemas/cnv-validation-traceability.schema.json`: generated from Pydantic.
- Create `docs/CNV_VALIDATION_BLOCK2_STATUS.md`: exact scope and evidence boundary.
- Modify `docs/CNV_VALIDATION_BLOCK1_STATUS.md`: mark Block 1 verified and point to Block 2 without changing Block-1 semantics.

### Task 1: Caller-native artifact and full-evidence contracts

**Interfaces**

Produces:
- `CnvEvidenceRecordKind`
- `CnvContributionStatus`
- `CnvRunOutcomeState`
- `CnvNumericMeasurement`
- `CnvNativeArtifactReference`
- `CnvFullEvidenceRecord`
- `canonical_evidence_sha256(value)`

Required record address:
`specimen × caller/version × build × data_basis × coverage × tumour_fraction × bin_size × record_kind/native_record_id × replicate`.

- [ ] **Step 1: RED tests**

Add synthetic tests proving:
- full records preserve continuous coverage/tumour fraction and their method/definition metadata;
- 100/500/1000-kbp records remain independently addressable;
- selected and alternative fits remain separate records;
- non-finite numeric measurements are rejected;
- caller parameter snapshots reject a mismatched SHA-256;
- artifact paths reject absolute/traversal/Windows-separator paths;
- `FAILED`, `NO_CALL`, `NOT_ASSESSABLE`, observed evidence, and an orthogonally-supported biological negative remain distinct;
- a biological-negative state requires explicit orthogonal truth resource IDs;
- excluded/secondary/exploratory records still remain valid full-evidence records.

Expected RED: import failure because `ontseq_platform.cnv_validation_evidence` does not exist.

- [ ] **Step 2: Verify RED in CI**

Expected failure must be caused by the missing evidence module, not syntax/style failures.

- [ ] **Step 3: GREEN implementation**

Implement strict models with:
- scalar identity/provenance fields;
- explicit `caller_parameters` plus validated `caller_parameters_sha256`;
- dependency versions;
- optional locus/event/copy-number/fit fields without invented defaults;
- list of caller-specific finite numeric measurements;
- `native_artifact_ids` linkage;
- fail-closed disposition/outcome consistency.

`CnvNativeArtifactReference` stores only run-relative POSIX paths plus SHA-256 and size; no host-absolute path is permitted.

- [ ] **Step 4: Verify GREEN**

Run focused tests and the full repository suite through PR CI.

### Task 2: Normalised-event layer and sealed evidence manifest

**Interfaces**

Produces:
- `CnvNormalizedEventRecord`
- `CnvValidationEvidenceManifest`
- `seal_cnv_validation_evidence(...) -> CnvValidationEvidenceManifest`

- [ ] **Step 1: RED tests**

Prove:
- a normalised event requires at least one source full-evidence ID;
- source full-evidence IDs are unique;
- source records must exist in the same manifest;
- source specimen/caller/version/build/data-basis/reference/bin-size/replicate identity must agree with the normalised record;
- terminal `FAILED`/`NO_CALL`/`NOT_ASSESSABLE` source records cannot be transformed into a positive normalised event;
- dangling artifact IDs fail closed;
- native artifacts referenced by retained evidence cannot disappear from the manifest;
- duplicate full-evidence address tuples fail closed;
- duplicate full, normalised, and artifact IDs fail closed;
- excluded source evidence can remain retained even when it contributes to no normalised event;
- tampering with any sealed manifest field invalidates `manifest_sha256`.

- [ ] **Step 2: Verify RED**

Expected failure: missing manifest/normalised-event implementation.

- [ ] **Step 3: GREEN implementation**

Manifest validators must:
- enforce unique identities;
- validate all artifact and source links;
- reject orphan caller-native artefacts;
- reject cross-specimen/cross-caller/cross-build/cross-data-basis normalisation;
- retain all full evidence regardless of contribution status;
- verify canonical manifest SHA-256 excluding only `manifest_sha256`.

`seal_cnv_validation_evidence` computes the lock from complete canonical JSON after validating unsealed content.

- [ ] **Step 4: Verify GREEN**

Focused tests and full CI.

### Task 3: Aggregate-to-source traceability contract without metric calculation

**Interfaces**

Produces:
- `CnvAggregateEvidenceTrace`
- `CnvTraceabilityIndex`
- `seal_cnv_traceability_index(...) -> CnvTraceabilityIndex`
- `verify_cnv_traceability(evidence, index) -> None`

- [ ] **Step 1: RED tests**

Prove:
- trace entries carry explicit aggregate/question IDs and stratum keys;
- numerator/denominator/excluded full-evidence IDs remain separate lists;
- normalised-event IDs can be linked independently;
- duplicate IDs inside a trace are rejected;
- overlap between numerator and excluded lists is rejected;
- dangling full-evidence or normalised-event IDs fail closed;
- `evidence_manifest_sha256` must match the exact sealed evidence manifest;
- index tampering invalidates `index_sha256`;
- a trace can walk aggregate -> normalised event -> source full evidence -> native artifact.

- [ ] **Step 2: Verify RED**

Expected failure must reflect missing traceability functionality.

- [ ] **Step 3: GREEN implementation**

Do not calculate sensitivity, specificity, TP/FP/FN, or other Block-3 metrics. Only define and validate provenance links and numerator/denominator membership declarations.

- [ ] **Step 4: Verify GREEN**

Focused tests and full CI.

### Task 4: Public schemas and documentation

- [ ] **Step 1: RED schema test/update**

Extend schema-export expectations for:
- `schemas/cnv-validation-evidence.schema.json`
- `schemas/cnv-validation-traceability.schema.json`

`python scripts/export_schemas.py --check` must fail until generated files are present.

- [ ] **Step 2: Generate exact schemas**

Generate from the final Pydantic models; do not hand-author schema JSON.

- [ ] **Step 3: Documentation**

`docs/CNV_VALIDATION_BLOCK2_STATUS.md` must state:
- full evidence is retained independent of contribution status;
- caller-native artefacts are fingerprinted and linked;
- normalised records are derived views, not replacements;
- biological negative requires orthogonal truth linkage and is distinct from no-call/failure;
- traceability is implemented but aggregate metric calculation remains Block 3;
- no real biological data were executed or committed.

- [ ] **Step 4: Verify schema freshness and repository safety**

Run schema check, safety check, Ruff, mypy, focused tests, full tests, DAGs, Desktop/MARLIN gates.

### Task 5: Final review and same-head verification

- [ ] Review the entire PR diff against the design spec with special attention to destructive filtering, hidden defaults, caller/truth conflation, cross-reference integrity, path safety, finite-number validation, and checksum semantics.
- [ ] Resolve every relevant finding before the final gate.
- [ ] Confirm no temporary workflow remains in the PR diff.
- [ ] Confirm no real sample/patient data or caller payload is committed.
- [ ] Verify the unchanged final head with CI, MARLIN runtime smoke, Desktop bundle contract, and Desktop CI.
- [ ] Record exact workflow run numbers and final commit SHA in PR #79.
- [ ] Keep PR Draft/unmerged; do not bump `0.8.2`.

## Acceptance criteria for Block 2

1. Every retained caller record is independently addressable and stays present even when excluded from primary analysis.
2. Caller-native files are referenced by immutable SHA-256 fingerprints and safe relative paths.
3. Continuous coverage/tumour-fraction values are preserved alongside categorical/technical context.
4. 100/500/1000-kbp evidence can coexist without key collision.
5. Selected and alternative fit evidence can coexist without destructive selection.
6. A normalised event always names its exact source full-evidence records.
7. Normalisation cannot silently cross specimen/caller/build/data-basis/reference boundaries.
8. Failure/no-call/not-assessable states cannot be transformed into positive normalised events.
9. Biological-negative state is distinct and requires orthogonal truth linkage.
10. Aggregate trace entries can be verified against exact full/normalised source IDs without calculating Block-3 metrics.
11. Evidence and traceability objects are content-locked and tamper-evident.
12. Versioned schemas are current and synthetic tests pass on Python 3.11/3.12/3.13.
13. No production caller is selected and no clinical-validity claim is made.
