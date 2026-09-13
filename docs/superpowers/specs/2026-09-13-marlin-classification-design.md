# ONTSeq MARLIN Classification Lane — Design

## Status

Approved design for implementation on `feat/marlin-classification-v1`.

This design records the user-approved **Approach B**: ONTSeq owns the typed contracts, provenance, feature construction, status semantics, validation and reporting boundaries; a pinned MARLIN R/Keras/TensorFlow runtime performs only the published model inference.

Tracking issue: [#76 — 0.9.0: Integrate MARLIN methylation classifier and validate against GSE280090](https://github.com/xhuljanoshehu/ontseq-platform/issues/76).

## Goal

Add a reproducible, fail-closed, Research-Use-Only MARLIN methylation-classification lane to ONTSeq that can:

1. import the published CpG-level Nanopore methylation representation used by MARLIN;
2. deterministically construct the exact model feature vector defined by the published MARLIN v1 workflow;
3. execute a locked MARLIN model runtime without allowing runtime-dependent preprocessing drift;
4. preserve all 42 class scores and expose confidence/unknown semantics without forcing a diagnosis;
5. record enough provenance to reproduce or reject a run byte-for-byte at the input/artifact/contract level;
6. validate the downstream CpG-to-classification path first against GSE280090, beginning with `GSM8587229_AL_001.txt.gz`;
7. later bridge the same classifier to ONTSeq-native `modBAM -> modkit -> CpG` output without conflating the two evidence levels.

The intended claim after this work is deliberately narrow:

> ONTSeq can reproducibly import published CpG-level Nanopore methylation data and execute a locked MARLIN downstream classification workflow with explicit provenance and fail-closed confidence semantics.

This work does **not** establish clinical validity and does **not** establish end-to-end raw-signal equivalence.

## Current ONTSeq baseline

The implementation starts from ONTSeq `main` at the 0.8.2 line. The repository already contains:

- a pinned modkit methylation lane using modkit 0.6.4 semantics;
- explicit `COMPLETED`, `NO_CALL`, `NOT_RUN` and `FAILED` handling across execution paths;
- strong Pydantic-based contracts and model validators;
- SHA-256 provenance/fingerprinting patterns;
- prospective methylation validation and holdout infrastructure;
- content-addressed/runtime-aware software locks;
- technical-vs-analytical-vs-clinical evidence separation;
- Research-Use-Only report semantics.

MARLIN must reuse these architectural principles rather than introduce an independent script-style subsystem.

## External scientific and software baseline

### Primary publication

Steinicke TL, Benfatto S et al. *Rapid epigenomic classification of acute leukemia.* Nature Genetics (2025). DOI: `10.1038/s41588-025-02321-z`.

### MARLIN code

- Repository: `https://github.com/hovestadt/MARLIN`
- Versioned code archive: `https://doi.org/10.5281/zenodo.15723932`
- GitHub code inspected for this design at commit:
  `442aa603415a54f62e7367794f9a31c6bc20fc2d`

The inspected official prediction implementation defines the v1 model-input semantics:

- input BED fields: chromosome, start, end, methylation call, probe identifier;
- features are ordered by the published `marlin_v1.features.RData` probe list;
- methylation values `>= 0.5` become `+1`;
- methylation values `< 0.5` become `-1`;
- absent/not-covered features become `0`;
- the trained model returns 42 class scores.

ONTSeq must reproduce those semantics explicitly and test them independently; it must not rely on incidental ordering or permissive behavior inside the original R scripts.

### MARLIN model

- Zenodo model record: `https://doi.org/10.5281/zenodo.15565404`
- Published model version: `1.0.0`

The model artifact is external to Git. Installation/import must compute and persist the local SHA-256; inference is refused if the artifact differs from the lock.

### MARLIN training data

- Zenodo record: `https://doi.org/10.5281/zenodo.15566584`
- Published version: `1.0.0`
- `marlin_v1.betas.RData`, 6.7 GB, Zenodo MD5 `8ace3200b17c0e7bda383bf62b6b1fb7`
- `marlin_v1.classes.RData`, 10.3 kB, Zenodo MD5 `84a8804eac2d979fd439227921c384f3`

These training files are not required for first-pass inference and therefore are **not** a dependency of the initial ONTSeq MARLIN runtime.

### External validation data

GEO series: `GSE280090` — 39 processed Nanopore methylation samples. GEO explicitly states that raw patient files were not deposited and that the processed files were generated from Dorado v0.3.2, modkit v0.1.13 and hg19.

First external validation target:

- accession: `GSM8587229`
- sample label: `AL_001`
- restricted processed file: `GSM8587229_AL_001.txt.gz`, approximately 4 MB
- assembly: hg19 / GRCh37
- GEO description: methylation status at CpG positions restricted to the reference-cohort CpGs

The companion `*_all_cpgs.txt.gz` files are not part of the first implementation. Initial scope is deliberately restricted to the published MARLIN probe-level representation.

No GSE280090 payload is committed to Git. Validation manifests may record accession, filename, expected schema, build and checksums; data remain external/local.

## Architectural decision

### Chosen approach: ONTSeq-native contracts + locked MARLIN runtime

The classifier is split into two trust domains:

1. **ONTSeq domain** — owns all data contracts, parsing, validation, feature mapping, provenance, status logic, threshold semantics, output normalization, atomic persistence and report integration.
2. **MARLIN runtime domain** — owns only loading the pinned trained model and producing the 42 model scores from an already constructed feature vector.

The runtime must not decide coordinate interpretation, feature ordering, missingness handling, confidence state, reportability or evidence level.

This preserves faithfulness to the published model while keeping ONTSeq's existing fail-closed and auditable architecture.

## Non-goals for v1

The first MARLIN integration intentionally does **not** include:

- model retraining;
- model fine-tuning;
- threshold optimization against GSE280090;
- automatic downloading of patient-derived validation files during normal CI or execution;
- automatic downloading of model artifacts during inference;
- support for arbitrary methylation BED schemas;
- support for the GSE280090 `*_all_cpgs.txt.gz` representation unless separately specified and versioned later;
- GRCh38 or T2T classification in the first validation lane;
- live/streaming classification;
- early-stop clinical decision logic;
- promotion of any MARLIN result to clinically reportable;
- claiming that precomputed GSE data validate ONTSeq's Dorado, alignment, MM/ML or modkit path.

These exclusions are deliberate to keep the first scientific gate narrow and falsifiable.

## Module boundaries

New functionality should be separated from the existing quantitative methylation module.

Recommended source layout:

```text
src/ontseq_platform/
    marlin_contracts.py
    marlin_input.py
    marlin_features.py
    marlin_runtime.py
    marlin_validation.py
    marlin_cli.py
```

Responsibilities:

### `marlin_contracts.py`

Owns immutable/strict schemas and enums only:

- source-kind enum;
- input schema identity;
- artifact lock;
- normalized CpG/probe observation;
- feature-vector summary;
- per-class score;
- classification decision;
- final MARLIN report;
- validation-manifest/result contracts.

No subprocess execution and no file-format parsing belongs here.

### `marlin_input.py`

Owns strict parsing of the published five-column MARLIN probe BED representation.

It accepts plain text or gzip input and returns normalized probe observations. It performs no model inference.

### `marlin_features.py`

Owns the deterministic published preprocessing contract:

- load the locked ordered probe identifiers;
- map input observations by probe identifier;
- transform observed beta/methylation fraction to `+1` or `-1` at exactly 0.5;
- encode absent features as `0`;
- compute observed/missing feature counts and fractions;
- emit an ordered numeric vector plus a content digest.

This module is the canonical oracle for feature construction. The R runtime must not recreate independent preprocessing logic.

### `marlin_runtime.py`

Owns:

- runtime preflight;
- verification of model/runtime artifact locks;
- creation of a minimal inference input from the canonical ONTSeq feature vector;
- execution of the pinned R/Keras/TensorFlow environment;
- strict parsing of exactly 42 finite scores;
- rejection of malformed/partial outputs;
- atomic publication of runtime output.

It does not classify confidence and does not know expected diagnoses.

### `marlin_validation.py`

Owns external downstream reproducibility evaluation:

- validation manifest for public processed CpG inputs;
- expected accession/sample metadata;
- expected build and input hash;
- expected published label or comparison target where independently available;
- repeated-run determinism checks;
- concordance/discordance/unknown accounting;
- cohort summaries;
- no post-hoc parameter changes.

It is explicitly distinct from the existing read-mixture/LoD validation subsystem because MARLIN validation measures classification reproducibility, not source-fraction recovery.

### `marlin_cli.py`

Owns the user-facing single-purpose CLI commands while keeping the main dispatcher thin.

Initial commands should be:

```text
ontseq marlin-lock
ontseq marlin-features
ontseq marlin-classify
ontseq marlin-validate
```

The final names may be wired through the repository's existing scientific CLI conventions, but there must be one canonical implementation per operation and no duplicate preprocessing path.

## Input contracts

### Source kinds

The normalized report must distinguish at minimum:

```text
PRECOMPUTED_METHYLATION
MODKIT_DERIVED
```

Only `PRECOMPUTED_METHYLATION` is enabled for the first GSE280090 validation gate.

The provenance distinction is load-bearing: a successful precomputed-data classification must never count as evidence that ONTSeq generated equivalent CpG values from a modBAM.

### Published MARLIN probe-BED adapter

Initial adapter identity:

```text
marlin_probe_bed_v1
```

Accepted file transport:

- uncompressed UTF-8/ASCII text;
- gzip-compressed text (`.gz`).

Accepted logical schema is exactly five tab-separated fields per non-empty row:

1. chromosome;
2. start;
3. end;
4. methylation fraction or `NA`;
5. MARLIN probe identifier.

Rules:

- no silently accepted extra columns;
- no inferred delimiter;
- no locale-dependent numeric parsing;
- start/end must be integers with `0 <= start < end`;
- methylation values must be finite and within `[0, 1]`, or the explicit missing token `NA`;
- probe identifiers must be non-empty;
- duplicate probe identifiers are rejected rather than resolved by first/last occurrence;
- malformed rows fail the adapter; they are not skipped;
- the declared genome build is mandatory and is never inferred from chromosome names;
- v1 external validation accepts only GRCh37/hg19 because GSE280090 was processed on hg19;
- an input SHA-256 is calculated before inference and written into the final report.

Rejecting duplicate probe IDs is intentionally stricter than R's permissive `match()` behavior. A valid published MARLIN input is expected to map one observation per probe; ambiguous input should fail closed rather than inherit incidental row ordering.

### Native modkit bridge

The future `MODKIT_DERIVED` path must not consume region-level aggregates from `methylation.py`; MARLIN needs probe-level CpG observations.

The bridge therefore requires a dedicated probe-coordinate extraction/mapping stage tied to the same build-specific MARLIN probe resource used by the artifact lock.

That bridge is a later phase of this design and becomes enabled only after precomputed-input classification is reproducible.

## Genome-build policy

The first implementation is **GRCh37/hg19 only**.

Although upstream MARLIN documentation now mentions hg38 and T2T support, ONTSeq must not generalize that support without build-specific locked probe coordinates and an explicit bridge test.

The artifact lock therefore includes a genome-build field, and input build must equal artifact-lock build.

A mismatch is a hard failure before feature construction.

Adding GRCh38 later requires a new artifact-lock identity and build-specific validation. No liftover is performed implicitly inside the classifier.

## Canonical feature construction

The feature contract must exactly reproduce the published v1 semantics while making each transformation explicit.

Given an ordered feature list `F = [f1, f2, ..., fn]` and a unique mapping from probe ID to observed methylation value:

```text
observed value >= 0.5  -> +1.0
observed value <  0.5  -> -1.0
missing / NA           ->  0.0
```

Important rules:

- feature order comes only from the locked MARLIN feature artifact;
- input row order cannot affect the vector;
- CpGs not present in the locked feature list are ignored for model input but counted separately as non-model observations;
- a locked feature present as explicit `NA` and a locked feature absent from the file both become model value `0`, but their origin is separately counted for QC/provenance;
- feature values are materialized in a deterministic numeric representation before hashing;
- vector hashing uses a documented canonical binary or canonical JSON encoding, not platform-dependent string formatting.

Feature summary must include:

- expected feature count;
- observed model-feature count;
- explicit-NA feature count;
- absent feature count;
- non-model probe count;
- observed fraction;
- feature-vector SHA-256;
- feature-artifact SHA-256;
- preprocessing-contract version.

No minimum observed-feature threshold is invented in v1 unless it is directly locked to a published MARLIN rule. Sparse inputs are passed to the model exactly as specified; confidence handling occurs on the resulting scores.

## Artifact lock

MARLIN is treated like a controlled reference/tool bundle, not like an incidental executable.

A `MarlinArtifactLock` must record at minimum:

```text
schema_version
lock_id
marlin_model_version
marlin_model_source_uri
marlin_model_sha256
marlin_code_source_uri
marlin_code_version_or_commit
marlin_code_sha256_or_manifest_sha256
feature_artifact_source
feature_artifact_sha256
class_annotation_source
class_annotation_sha256
genome_build
expected_class_count
preprocessing_contract_version
runtime_contract_version
R_version
keras_version
tensorflow_version
python_version_if_used
created_at
research_only=true
```

Requirements:

- lock creation is a separate explicit operation;
- inference never downloads or mutates the locked artifacts;
- each file is rehashed during preflight;
- mismatch is `FAILED`, not a warning;
- symlinks/reparse points must be resolved according to the repository's existing runtime-safety conventions before hashing/launch;
- lock files contain metadata/hashes only, not model bytes.

The first supported model identity is the published MARLIN model version `1.0.0` from Zenodo record `15565404`.

## Runtime design

### Boundary

The runtime receives only:

- one canonical ordered feature vector;
- the verified model path;
- verified class annotations;
- verified runtime environment.

The runtime returns only:

- exactly 42 named scores;
- runtime metadata needed for provenance;
- exit status/log metadata.

### Original-runtime fidelity

The first reference runtime should retain the MARLIN-published R/Keras/TensorFlow stack rather than reimplement the model in a new framework.

The published README lists tested versions including R 4.1.3, Keras 2.13 and TensorFlow 2.13. ONTSeq packaging may use a reproducibly installed compatible environment, but every actual version becomes part of the lock and must be validated against the fixed reference score fixture before use.

### CPU and GPU

CPU execution is sufficient for correctness validation and is the reference CI/local path where practical. GPU execution may be supported later for performance but must demonstrate score equivalence within a predeclared numerical tolerance before being treated as the same runtime profile.

Runtime profile identity must therefore include execution backend when numerical behavior can differ.

### Output validation

A runtime result is accepted only if:

- process exits successfully;
- exactly 42 class entries are present;
- class names exactly match the locked class annotation order/set;
- every score is finite;
- no duplicate class labels exist;
- scores satisfy the model's expected probability-output constraints;
- output file/artifact is complete and parses under a strict schema.

Any partial or malformed output is `FAILED`; it cannot degrade to `NO_CALL`.

### Atomic persistence

Results are written into a staging location and published atomically or with the repository's existing cross-filesystem-safe commit protocol.

This is required because ONTSeq is routinely used through Windows/WSL boundaries where a POSIX directory rename may not behave identically on a Windows-mounted filesystem.

A failed publication must never leave a directory that appears to be a completed reusable analysis.

## Result contract

Recommended normalized contract:

```text
MarlinClassificationReport
├── schema_version
├── sample_id
├── status                       # existing ModuleRunStatus
├── classification_decision      # MARLIN-specific decision
├── source_kind
├── input_schema
├── genome_build
├── artifact_lock_id
├── input_fingerprint
├── feature_summary
├── class_scores[42]
├── top_class
├── top_class_score
├── lineage_scores
├── methylation_family_scores
├── warnings
├── limitations
├── runtime_provenance
└── research_only = true
```

All 42 model scores are retained in normalized output. The report must not discard lower-ranked classes simply because one top class exists.

### Module status

Use the existing generic status contract:

- `COMPLETED` — valid model inference exists, regardless of confidence;
- `NO_CALL` — the input is technically valid but no biologically interpretable classification can be issued under an explicitly defined decision rule;
- `NOT_RUN` — classifier was not requested or required resources were intentionally unavailable before execution;
- `FAILED` — corrupted/mismatched artifacts, invalid input, runtime failure or malformed output.

For v1, a low-confidence valid model run remains technically `COMPLETED`; its MARLIN-specific decision is `LOW_CONFIDENCE` / `UNKNOWN`. This avoids conflating successful computation with classification confidence.

### MARLIN-specific decision

Initial decision enum:

```text
HIGH_CONFIDENCE
LOW_CONFIDENCE
UNKNOWN
```

The design preserves a distinct `UNKNOWN` representation for inputs outside a supported decision context. `LOW_CONFIDENCE` is the direct state for a valid prediction whose class-level confidence does not satisfy the locked high-confidence rule.

## Confidence policy

The initial high-confidence policy is locked to the published MARLIN threshold of `0.8` for the relevant prediction score aggregation used by the publication/runtime annotation logic.

Important constraints:

- the exact score level being thresholded — class, methylation-family or lineage aggregation — must be encoded in the policy and tested against the published class annotations;
- ONTSeq does not tune `0.8` against GSE280090;
- a score below threshold is not promoted to a diagnosis;
- a later threshold change requires a new policy version and separate validation, never silent configuration drift.

The implementation plan must derive the aggregation logic directly from the published class-annotation artifact and reference plotting/inference scripts before writing production decision code.

## Lineage and family aggregation

The published MARLIN class annotations include mappings used to aggregate the 42 class scores into lineage and methylation-class-family views.

ONTSeq should preserve this hierarchy as derived output:

```text
42 raw class scores
   -> grouped lineage scores
   -> grouped methylation-family scores
```

The mapping is read from the locked class-annotation artifact. ONTSeq must not maintain a second hand-written mapping table.

Raw 42-class scores remain the primary model output; group scores are deterministic derived summaries.

## Error handling and fail-closed behavior

### Invalid input

Examples:

- wrong field count;
- non-tab delimiter under the v1 adapter;
- non-integer coordinates;
- coordinate inversion;
- beta outside `[0,1]`;
- non-finite beta;
- duplicate probe ID;
- empty probe ID;
- declared build mismatch.

Result: input validation error / `FAILED`. No rows are silently dropped.

### Artifact mismatch

Examples:

- changed HDF5 model;
- changed feature file;
- changed class annotations;
- wrong runtime version under a strict profile.

Result: hard preflight failure / `FAILED`.

### Valid but sparse input

Missing model features become `0` according to the published preprocessing contract. Sparsity is surfaced quantitatively in `feature_summary`.

A valid sparse run may return low confidence; that is not a technical failure.

### No usable model features

If no locked model features are observed and all feature entries would be zero, ONTSeq must refuse to present a meaningful prediction. The implementation may execute a deterministic guard before the external model runtime and emit a `NO_CALL` with the reason `NO_MODEL_FEATURES_OBSERVED`; it must not expose a model's prior-like output as sample evidence.

This guard is an ONTSeq safety boundary and must be covered by an explicit regression test.

## CLI and operator workflow

Initial operator flow:

```text
1. Install/prepare MARLIN runtime and artifacts offline.
2. `ontseq marlin-lock ...` creates a checksummed artifact lock.
3. `ontseq marlin-features ...` optionally inspects/exports the canonical feature summary.
4. `ontseq marlin-classify ...` runs one sample.
5. `ontseq marlin-validate ...` evaluates a pre-registered public processed-data manifest.
```

The classify command requires an explicit genome build and artifact lock.

No command auto-detects build from filename or downloads a model implicitly.

## Validation strategy

Validation is staged so each claim corresponds to one evidence layer.

### Phase 0 — pure preprocessing oracle

Before invoking the real model runtime, create synthetic fixtures that independently test:

- row-order invariance;
- exact `0.5` threshold behavior;
- `<0.5` -> `-1`;
- `>=0.5` -> `+1`;
- explicit `NA` -> `0`;
- absent feature -> `0`;
- feature-order independence from input order;
- duplicate rejection;
- malformed numeric rejection;
- gzip/plain equivalence;
- vector-digest determinism.

The expected vectors must be written by hand/independent arithmetic, not generated by the implementation under test.

### Phase 1 — fixed synthetic runtime oracle

Construct a small deterministic canonical feature vector, run it through the locked real MARLIN runtime and store the resulting 42-score vector as a versioned reference fixture with full artifact/runtime identity.

Repeated runs must satisfy:

- identical class names/order;
- stable top class;
- score differences within a predeclared numerical tolerance;
- deterministic normalized JSON after excluding intentionally variable timestamps.

The tolerance is measured and then frozen before external biological validation; it is not adjusted after seeing GSE280090 results.

### Phase 2 — AL_001 external processed-data gate

Input:

`GSM8587229_AL_001.txt.gz`

Requirements:

- public accession and local SHA-256 recorded;
- explicit GRCh37/hg19 declaration;
- deterministic parse and feature digest;
- deterministic 42-score output under repeated execution;
- confidence decision follows the locked policy;
- complete provenance written to JSON/report;
- no upstream Dorado/modkit claim is made.

Expected biological comparison is recorded from the publication/source-data context before scoring is evaluated. Any expected score/tolerance entered into a validation manifest must come from an independently documented source, not from the first ONTSeq run.

### Phase 3 — GSE280090 cohort

Run all compatible restricted probe-level files under the same locked artifact/runtime/policy.

Predeclared output metrics:

- number attempted;
- number technically completed;
- number failed input/runtime;
- high-confidence count;
- low-confidence/unknown count;
- top-class concordance where an independent expected class is available;
- high-confidence concordance;
- no-call fraction;
- feature-coverage distribution;
- runtime distribution;
- deterministic-rerun failures;
- full list of discordant samples for manual review.

No model, threshold or preprocessing rule is modified after outcomes are inspected within the same registered validation run.

### Phase 4 — native ONTSeq bridge

Only after the processed-data path is stable:

```text
same specimen / equivalent source
published or locked precomputed CpG input
                versus
ONTSeq modBAM -> modkit 0.6.4 -> probe-level CpG extraction
```

Compare:

- probe IDs observed;
- methylation fractions before binarization;
- model feature vector;
- feature digest;
- 42 MARLIN scores;
- top class and confidence state.

This phase is the first evidence that ONTSeq's native methylation extraction interoperates with MARLIN. It still does not constitute clinical validation.

### Phase 5 — streaming/realtime, deferred

Real-time MARLIN inference is a separate design because early threshold crossing can be unstable over time. It requires a predeclared stability policy, cumulative evidence semantics and explicit rules for score reversals.

No realtime/early-stop behavior is added in the first implementation.

## Testing strategy

Implementation follows strict test-driven development.

### Unit tests

Add focused test modules for:

- contracts and validators;
- probe-BED parser;
- gzip equivalence;
- feature mapping and binarization;
- artifact-lock verification;
- score parser;
- class-annotation aggregation;
- confidence policy;
- no-model-features guard;
- deterministic hashing;
- serialization round-trips.

### Metamorphic/property tests

At minimum:

- input row permutation does not change the feature vector;
- gzip vs plain-text transport does not change semantic results;
- class-score serialization/deserialization preserves all values;
- changing any locked artifact byte invalidates preflight;
- changing build invalidates the run;
- repeated runs with identical locked inputs produce equivalent normalized results;
- report rendering does not mutate normalized results.

Where useful, Hypothesis should generate valid/invalid numeric boundaries rather than duplicating implementation logic.

### Integration tests

A synthetic runtime fixture exercises the actual R/Keras/TensorFlow adapter when the runtime is present.

CI must have two levels:

1. lightweight contract/preprocessing tests that are always blocking;
2. MARLIN-runtime integration job that is blocking on the designated packaged-runtime workflow but may be skipped with an explicit reason where the external model/runtime is intentionally absent.

A skipped runtime test must not be reported as validation success.

### Existing-suite protection

Run the full existing repository checks before merge:

- pytest suites;
- Ruff lint/format;
- mypy;
- schema consistency;
- version consistency;
- repository-safety checks;
- relevant Snakemake/runtime dry-runs;
- existing methylation tests;
- packaging/install smoke tests for the selected 0.9.0 bundle.

No existing biological output fixture may change as a side effect of MARLIN integration unless a separate independently justified defect fix is reviewed.

## Security and local execution

The model/runtime is third-party executable content and must be treated as such.

Requirements:

- no network access is required during inference;
- model and runtime are installed from explicit user/admin actions, never silently fetched at run time;
- hashes are verified before every run or via an existing content-addressed verified cache whose validity is rechecked according to ONTSeq runtime policy;
- subprocess working directory is isolated;
- output paths are constrained to the run workspace;
- inherited environment is minimized;
- cancellation terminates the complete process tree using the already hardened ONTSeq process-control primitives;
- runtime logs are captured without exporting patient sequence/read identifiers;
- no raw BAM/POD5/FASTQ content enters GitHub artifacts.

## Packaging

MARLIN runtime/model packaging is separate from Python source packaging.

The portable ONTSeq bundle may include a verified runtime archive only when licensing and redistribution terms permit it. Otherwise the bundle includes:

- installer/import workflow;
- expected source URI/version;
- checksum/lock creation tooling;
- clear missing-resource status.

The application must distinguish:

```text
runtime unavailable
model unavailable
artifact mismatch
input invalid
inference failed
inference completed low confidence
inference completed high confidence
```

A missing optional MARLIN runtime must not break unrelated ONTSeq analysis profiles.

## Reporting and UI boundary

The scientific core is implemented and validated before UI integration.

First mergeable scientific increment:

- contracts;
- parser;
- feature construction;
- artifact lock;
- runtime adapter;
- CLI;
- normalized JSON;
- validation runner.

UI/report additions consume only the normalized report contract.

The UI must display at minimum:

- MARLIN status;
- top class and score;
- confidence state;
- number/fraction of observed model CpGs;
- clear Research-Use-Only label;
- source kind (`PRECOMPUTED_METHYLATION` vs `MODKIT_DERIVED`);
- limitations when upstream extraction was not validated.

All 42 scores may be available in an expandable reviewer view but need not dominate the summary screen.

## Versioning

The work belongs to the 0.9.0 development line because it adds a new classification/validation subsystem and aligns with the existing independent-validation audit design.

Versioning rules:

- no `0.9.0` release is declared merely because the MARLIN lane lands;
- package/runtime/desktop/report schemas must remain version-consistent under existing guards;
- changelog wording must state that MARLIN integration is Research Use Only and not an analytically or clinically validated diagnostic release;
- the artifact-lock and MARLIN-result schemas have independent schema versions so their evolution does not rely only on package version.

## Planned repository changes

Expected new files:

```text
src/ontseq_platform/marlin_contracts.py
src/ontseq_platform/marlin_input.py
src/ontseq_platform/marlin_features.py
src/ontseq_platform/marlin_runtime.py
src/ontseq_platform/marlin_validation.py
src/ontseq_platform/marlin_cli.py

tests/test_marlin_contracts.py
tests/test_marlin_input.py
tests/test_marlin_features.py
tests/test_marlin_runtime.py
tests/test_marlin_validation.py

docs/MARLIN_CLASSIFICATION.md
configs/marlin/marlin_v1.research.yaml
```

Existing files likely touched in focused ways:

```text
src/ontseq_platform/entrypoint.py
src/ontseq_platform/cli.py or the canonical scientific-command registration surface
CHANGELOG.md
version/schema guard inputs if required
report/UI consumers only in a later reviewed increment
```

Large model/data artifacts are never committed to the source repository.

## Implementation sequence

The implementation plan should preserve the following dependency order:

1. contracts/enums with failing tests;
2. strict precomputed probe-BED parser;
3. canonical feature construction and digest;
4. artifact lock and preflight verification;
5. runtime score adapter using a synthetic fixed vector;
6. MARLIN-specific decision/aggregation semantics;
7. CLI and normalized JSON output;
8. validation manifest/runner;
9. AL_001 external gate;
10. full compatible GSE280090 cohort gate;
11. native modkit-to-MARLIN bridge;
12. UI/report consumption;
13. separate future streaming design.

Do not start at step 9 and use AL_001 as the debugging oracle for steps 1–8.

## Acceptance criteria

The first MARLIN classification increment is complete when all of the following hold:

- [ ] published MARLIN v1 preprocessing semantics are encoded in one canonical ONTSeq implementation;
- [ ] the external MARLIN model, feature and class-annotation artifacts are checksum-locked;
- [ ] inference is offline and refuses artifact mismatch;
- [ ] `marlin_probe_bed_v1` supports plain and gzip transport with strict schema validation;
- [ ] v1 requires explicit GRCh37/hg19 and rejects build mismatch;
- [ ] feature vector is independent of input row order;
- [ ] exact 0.5 threshold and missing->0 semantics are regression-tested;
- [ ] duplicate probes and malformed numerics fail closed;
- [ ] no-observed-model-feature input returns an explicit non-diagnostic state rather than a prior-like prediction;
- [ ] runtime returns and ONTSeq preserves exactly 42 finite class scores;
- [ ] lineage/family aggregation derives only from the locked class annotation artifact;
- [ ] high-confidence policy is versioned and not tuned on GSE280090 outcomes;
- [ ] repeated synthetic runtime executions meet a frozen numerical-equivalence criterion;
- [ ] `GSM8587229_AL_001.txt.gz` is processed reproducibly with recorded accession/build/hash/provenance;
- [ ] GSE280090 validation can run from an external-data manifest without committing patient-derived payloads;
- [ ] validation output separates technical completion, confidence/unknown and biological concordance;
- [ ] `PRECOMPUTED_METHYLATION` cannot be mistaken for evidence of native `MODKIT_DERIVED` validation;
- [ ] full existing ONTSeq CI/regression checks remain green;
- [ ] documentation states the exact supported scientific claim and the stronger unsupported claims;
- [ ] no clinical-reportable status is introduced.

## Risks and mitigations

### Runtime/version drift

Risk: old Keras/TensorFlow HDF5 behavior can differ across environments.

Mitigation: lock full runtime identity, verify fixed-vector score equivalence, and keep preprocessing outside the ML runtime.

### Build/coordinate drift

Risk: probe resources may exist for hg19, hg38 and T2T and can be accidentally mixed.

Mitigation: artifact lock is build-specific; v1 supports GRCh37 only; no implicit liftover.

### Overclaim from public processed data

Risk: strong GSE280090 concordance could be mislabeled as end-to-end validation.

Mitigation: source-kind provenance is mandatory and validation reports explicitly state that Dorado/alignment/MM-ML/modkit are outside this evidence layer.

### Sparse-input false certainty

Risk: a neural network always emits scores, even if evidence is minimal.

Mitigation: preserve feature-coverage metrics, published confidence semantics and explicit low-confidence/unknown state; refuse a zero-observed-feature prior-like output.

### Numerical nondeterminism

Risk: CPU/GPU/threading differences can alter low-order score bytes.

Mitigation: define semantic/numerical equivalence rather than blindly requiring byte-identical floating output; freeze tolerance before external biological validation and include runtime profile in provenance.

### Patient-data handling

Risk: validation inputs are processed human genomic/epigenomic data.

Mitigation: no payloads in Git; public accession/checksum manifests only; local paths do not enter public reports; follow existing ONTSeq data-governance rules.

## Scientific boundary after completion

Successful completion supports:

> ONTSeq reproducibly implements the published MARLIN v1 CpG-to-classification semantics against a pinned model/runtime and can evaluate public processed Nanopore methylation inputs with explicit provenance and confidence handling.

It does not support:

> ONTSeq clinically diagnoses acute leukemia from Nanopore sequencing.

It also does not yet support:

> ONTSeq's raw-signal/Dorado/modkit path is analytically equivalent to the MARLIN publication.

Those stronger statements require the later native bridge, analytical performance characterization, intended-use cohorts, orthogonal truth, governance and human review.
