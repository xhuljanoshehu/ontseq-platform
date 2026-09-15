# ONTSeq MARLIN Classification Lane — Design

## Status

Approved design for implementation on `feat/marlin-classification-v1`.

This document records the user-approved **Approach B**: ONTSeq owns all typed contracts, provenance, input parsing, feature construction, confidence semantics, validation and reporting boundaries; a pinned MARLIN R/Keras/TensorFlow runtime performs only the published neural-network inference.

Tracking issue: [#76 — 0.9.0: Integrate MARLIN methylation classifier and validate against GSE280090](https://github.com/xhuljanoshehu/ontseq-platform/issues/76).

## Goal

Add a reproducible, fail-closed, Research-Use-Only MARLIN methylation-classification lane that can:

1. import the published CpG-level Nanopore methylation representation used by MARLIN;
2. deterministically construct the exact MARLIN v1 model feature vector;
3. execute a checksum-locked MARLIN model/runtime without allowing preprocessing drift inside the ML runtime;
4. preserve all 42 model-unit scores and derive class/family/lineage summaries from the locked annotation resource;
5. expose high-confidence versus unknown classification without forcing a diagnosis;
6. record enough provenance to reproduce or reject a run at the input/artifact/contract level;
7. validate the downstream CpG-to-classification path first against GSE280090, beginning with `GSM8587229_AL_001.txt.gz`;
8. later bridge the same classifier to ONTSeq-native `modBAM -> modkit -> probe-level CpG` output without conflating the two evidence levels.

The intended claim after this work is deliberately narrow:

> ONTSeq can reproducibly import published CpG-level Nanopore methylation data and execute a locked MARLIN downstream classification workflow with explicit provenance and confidence handling.

This work does **not** establish clinical validity and does **not** establish end-to-end raw-signal equivalence.

## Current ONTSeq baseline

Implementation starts from the 0.8.2 `main` line. Existing relevant infrastructure includes:

- pinned modkit 0.6.4 methylation semantics;
- explicit `COMPLETED`, `NO_CALL`, `NOT_RUN` and `FAILED` execution states;
- Pydantic-based strict contracts and validators;
- SHA-256 provenance/fingerprinting patterns;
- prospective methylation validation/holdout infrastructure;
- software/runtime identity locks;
- fail-closed numeric and input handling;
- technical-vs-analytical-vs-clinical evidence separation;
- Research-Use-Only reporting.

MARLIN must reuse these principles rather than become an independent script-style subsystem.

## External scientific/software baseline

### Primary publication

Steinicke TL, Benfatto S et al. *Rapid epigenomic classification of acute leukemia.* Nature Genetics (2025). DOI: `10.1038/s41588-025-02321-z`.

### MARLIN code

- Repository: `https://github.com/hovestadt/MARLIN`
- Versioned code: `https://doi.org/10.5281/zenodo.15723932`
- GitHub implementation inspected for this design at commit:
  `442aa603415a54f62e7367794f9a31c6bc20fc2d`

The inspected official code establishes the core v1 semantics:

- model input contains **357,340 ordered CpG features**;
- the network has two hidden dense layers (256 and 128 nodes) and a **42-node softmax output**;
- training beta values are binarized at `0.5` to `-1/+1`;
- prediction inputs are ordered by `marlin_v1.features.RData`;
- observed values `>= 0.5` become `+1`;
- observed values `< 0.5` become `-1`;
- absent/not-covered features become `0`;
- class annotations are loaded from `marlin_v1.class_annotations.xlsx`.

ONTSeq must reproduce these semantics explicitly and test them independently; it must not depend on incidental input ordering or permissive behavior inside the original R scripts.

### MARLIN model

- Zenodo model: `https://doi.org/10.5281/zenodo.15565404`
- Published model version: `1.0.0`

The model stays external to Git. Installation/import computes a local SHA-256 and produces a lock; inference refuses any subsequent hash mismatch.

### MARLIN reference resources

At the inspected upstream commit, the MARLIN repository publishes:

- `marlin_v1.features.RData`;
- `marlin_v1.class_annotations.xlsx`;
- `marlin_v1.probes_hg19.bed.gz`;
- `marlin_v1.probes_hg38.bed.gz`;
- `marlin_v1.probes_t2t.bed.gz`.

The first ONTSeq validation lane uses only the hg19/GRCh37 resource. Presence of hg38/T2T resources upstream does not imply that ONTSeq has validated those builds.

### MARLIN training data

- Zenodo: `https://doi.org/10.5281/zenodo.15566584`
- Published version: `1.0.0`
- `marlin_v1.betas.RData`, 6.7 GB, Zenodo MD5 `8ace3200b17c0e7bda383bf62b6b1fb7`
- `marlin_v1.classes.RData`, 10.3 kB, Zenodo MD5 `84a8804eac2d979fd439227921c384f3`

Training data are not required for first-pass inference and are not an initial runtime dependency.

### External validation data

GEO series `GSE280090` contains 39 processed Nanopore methylation samples. GEO states that patient raw files were not deposited and that processed calls were produced with Dorado v0.3.2, modkit v0.1.13 and hg19.

First external target:

- accession: `GSM8587229`;
- sample label: `AL_001`;
- file: `GSM8587229_AL_001.txt.gz`;
- size: approximately 4 MB;
- assembly: hg19 / GRCh37;
- GEO description: methylation status at CpG positions restricted to the reference-cohort CpGs.

The companion `*_all_cpgs.txt.gz` files are outside v1 scope. Initial support is deliberately restricted to the MARLIN probe-level representation.

No GSE280090 payload is committed to Git. Validation manifests may store accession, expected filename/schema/build and hashes; payloads remain external/local.

## Architectural decision

### Chosen approach: ONTSeq-native contracts + locked MARLIN runtime

The system is split into two trust domains:

1. **ONTSeq domain** — contracts, parsing, build validation, feature mapping, provenance, artifact locks, status/confidence semantics, score validation, output normalization, validation and reporting.
2. **MARLIN runtime domain** — load the verified trained model and infer 42 softmax scores from an already constructed feature vector.

The runtime must not decide coordinate interpretation, feature ordering, missing-value semantics, confidence state, reportability or evidence level.

## Non-goals for v1

The first integration does not include:

- retraining or fine-tuning;
- threshold optimization against GSE280090;
- automatic model/data downloads during inference;
- arbitrary methylation BED schemas;
- `*_all_cpgs.txt.gz` input;
- GRCh38 or T2T classification;
- implicit liftover;
- live/streaming classification;
- early-stop decision logic;
- clinical reportability;
- claims that precomputed GSE data validate Dorado, alignment, MM/ML or ONTSeq modkit extraction.

These exclusions keep the first scientific gate narrow and falsifiable.

## Module boundaries

New functionality remains separate from the existing region-aggregated quantitative methylation lane.

Recommended layout:

```text
src/ontseq_platform/
    marlin_contracts.py
    marlin_input.py
    marlin_features.py
    marlin_runtime.py
    marlin_validation.py
    marlin_cli.py
```

### `marlin_contracts.py`

Owns strict schemas/enums only:

- source kind;
- input schema identity;
- artifact/runtime lock;
- normalized probe observation;
- feature summary;
- model-unit score;
- grouped class/family/lineage score;
- classification decision;
- final report;
- validation manifest/result.

No subprocess execution or input parsing belongs here.

### `marlin_input.py`

Strict parser for published five-column MARLIN probe BED. Supports plain text and gzip. Returns normalized probe observations and performs no inference.

### `marlin_features.py`

Canonical preprocessing oracle:

- load and verify the locked ordered 357,340 probe identifiers;
- map observations by probe ID;
- encode `>=0.5 -> +1`, `<0.5 -> -1`, absent/NA -> `0`;
- compute observed/NA/absent/non-model counts;
- emit the ordered vector and deterministic SHA-256.

The R runtime must not implement a second independent preprocessing path.

### `marlin_runtime.py`

Owns:

- runtime/artifact preflight;
- launch of the verified R/Keras/TensorFlow environment;
- transfer of the already-built canonical feature vector;
- strict parsing and validation of exactly 42 softmax scores;
- runtime provenance;
- cross-filesystem-safe result publication.

It does not know expected diagnoses or confidence thresholds.

### `marlin_validation.py`

Owns downstream external reproducibility evaluation:

- public processed-data manifest;
- accession/build/hash locks;
- independent expected comparison target where available;
- repeated-run determinism;
- concordant/discordant/unknown accounting;
- cohort metrics;
- prohibition of post-hoc threshold/preprocessing changes.

It is separate from ONTSeq's methylation mixture/LoD infrastructure because it evaluates classification reproducibility rather than source-fraction recovery.

### `marlin_cli.py`

Owns single-purpose CLI operations while keeping the global dispatcher thin.

Initial commands:

```text
ontseq marlin-lock
ontseq marlin-features
ontseq marlin-classify
ontseq marlin-validate
```

There must be one canonical implementation per operation and no duplicate preprocessing route.

## Input contracts

### Source kinds

Final reports distinguish:

```text
PRECOMPUTED_METHYLATION
MODKIT_DERIVED
```

Only `PRECOMPUTED_METHYLATION` is enabled for the first external gate.

This distinction is load-bearing: successful GSE classification can never count as evidence that ONTSeq itself generated equivalent CpG values from modBAM.

### Published MARLIN probe-BED adapter

Adapter identity:

```text
marlin_probe_bed_v1
```

Accepted transport:

- uncompressed UTF-8/ASCII text;
- gzip-compressed text.

Logical schema is exactly five tab-separated fields per non-empty row:

1. chromosome;
2. start;
3. end;
4. methylation fraction or literal `NA`;
5. MARLIN probe identifier.

Rules:

- extra columns are rejected;
- delimiter is not guessed;
- numeric parsing is locale-independent;
- coordinates satisfy `0 <= start < end`;
- methylation is finite within `[0,1]` or explicit `NA`;
- probe ID is non-empty;
- duplicate probe IDs are rejected;
- malformed rows fail the whole adapter and are not skipped;
- genome build is mandatory and never inferred from names;
- v1 accepts only GRCh37/hg19;
- input SHA-256 is calculated before inference and retained in the report.

Rejecting duplicate probes is intentionally stricter than R's permissive `match()` behavior. Ambiguous inputs fail closed instead of inheriting row order.

### Native modkit bridge

The future `MODKIT_DERIVED` path cannot use region-level aggregates from `methylation.py`; MARLIN requires probe-level CpGs.

The bridge therefore gets its own probe-coordinate extraction/mapping stage tied to the same build-specific probe resource recorded by the MARLIN artifact lock. It remains disabled until the precomputed path is reproducible.

## Genome-build policy

v1 is **GRCh37/hg19 only**.

Artifact lock contains the build and the hg19 probe-resource hash. Input declaration must match it. Mismatch is a hard failure before feature construction.

GRCh38/T2T later require separate build-specific artifact locks and explicit bridge validation. No classifier-internal liftover is allowed.

## Canonical feature construction

Expected feature count is exactly **357,340** for MARLIN v1.

For ordered feature list `F` and unique probe observations:

```text
observed methylation >= 0.5 -> +1.0
observed methylation <  0.5 -> -1.0
explicit NA                 ->  0.0
feature absent from input   ->  0.0
```

Rules:

- order comes only from locked `marlin_v1.features.RData`;
- input row order cannot affect output;
- probes outside the model feature list do not enter the vector but are counted;
- explicit NA and absent model features both map to zero but are counted separately;
- vector hashing uses a documented canonical representation, not platform-specific float formatting.

Feature summary contains:

- expected feature count (`357340`);
- observed model-feature count;
- explicit-NA count;
- absent-feature count;
- non-model probe count;
- observed fraction;
- vector SHA-256;
- feature-artifact SHA-256;
- preprocessing-contract version.

No new minimum-observed-feature threshold is invented for v1. Sparse valid inputs are passed to the model according to published semantics; confidence is evaluated from model output. A separate zero-evidence guard is defined below.

## Artifact lock

MARLIN is treated as a controlled reference/tool bundle.

`MarlinArtifactLock` records at minimum:

```text
schema_version
lock_id
model_version
model_source_uri
model_sha256
code_source_uri
code_version_or_commit
code_manifest_sha256
feature_source_uri
feature_sha256
class_annotation_source_uri
class_annotation_sha256
probe_resource_source_uri
probe_resource_sha256
genome_build
expected_feature_count = 357340
expected_model_unit_count = 42
preprocessing_contract_version
runtime_contract_version
R_version
keras_version
tensorflow_version
python_version_if_used
execution_backend
created_at
research_only = true
```

Requirements:

- lock creation is explicit;
- inference performs no download/mutation;
- artifact hashes are reverified before inference or through an existing verified content-addressed runtime cache;
- mismatch is hard `FAILED`;
- symlink/reparse handling follows existing runtime-safety rules;
- lock contains metadata/hashes, not model bytes.

First supported model identity is MARLIN model `1.0.0`, Zenodo record `15565404`.

## Runtime design

### Boundary

Runtime receives only:

- one canonical 357,340-value feature vector;
- verified model;
- verified class annotations;
- verified runtime.

Runtime returns only:

- exactly 42 raw model-unit scores;
- execution metadata/log identity.

### Runtime fidelity

First reference runtime retains MARLIN's R/Keras/TensorFlow model-loading path rather than reimplementing the neural network in a new framework.

Published tested versions include R 4.1.3, Keras 2.13 and TensorFlow 2.13. The actual packaged versions must be recorded in the lock and pass the fixed-vector compatibility gate.

### CPU/GPU

CPU is the correctness/reference path where practical. A GPU profile may be added later only after equivalence against the frozen reference fixture.

Execution backend is part of runtime identity.

### Score invariants

The official training code ends in a 42-node softmax layer. Therefore an accepted raw output must satisfy all of:

- exactly 42 entries;
- exact locked model-unit names/order after annotation binding;
- every value finite;
- every value in `[0,1]`;
- no duplicate model-unit labels;
- raw score sum differs from `1.0` by no more than `1e-5`;
- output parses under the strict runtime-output schema.

A violation is `FAILED`, never `NO_CALL`.

### Numerical compatibility profile

Cross-runtime score equivalence is not assumed byte-identical because TensorFlow/backend/threading can alter low-order floating values.

Before any GSE280090 biological validation, Phase 1 creates and freezes a `MarlinRuntimeCompatibilityProfile` containing:

- reference runtime lock ID;
- fixed feature-vector SHA-256;
- reference 42-score vector;
- per-score absolute tolerance;
- top-class identity requirement;
- raw-score-sum tolerance.

No external validation may run under a runtime profile that has not passed this frozen compatibility fixture. The tolerance is established on synthetic/reference-vector runs **before** GSE outcomes are inspected and then version-locked.

### Atomic/cross-filesystem persistence

Runtime output is written to staging and published with ONTSeq's cross-filesystem-safe result-commit protocol.

This requirement is explicit because Windows/WSL-mounted filesystems may not support the same atomic directory operations as ext4. Failed publication must not leave a directory that resume logic can mistake for a completed result.

## Class-annotation and confidence semantics

The 42 network outputs are raw model-unit scores and are always preserved.

The locked `marlin_v1.class_annotations.xlsx` defines deterministic groupings used by the upstream realtime plotting code:

```text
42 raw model-unit scores
    -> class_name_current groups
    -> methylation-family (`mcf`) groups
    -> lineage groups
```

ONTSeq derives each grouped score by summing the raw softmax scores assigned to that group in the locked annotation file. No second hand-written mapping is allowed.

### Primary v1 classification target

For v1, `top_class` is the `class_name_current` group with the largest grouped score.

The high-confidence rule is:

```text
max(class_name_current grouped score) >= 0.8 -> HIGH_CONFIDENCE
max(class_name_current grouped score) <  0.8 -> UNKNOWN
```

`UNKNOWN` therefore means **valid model inference without a high-confidence current-class assignment**. It is not a technical failure and does not erase the top candidate or its score.

Lineage and methylation-family grouped scores are retained as secondary summaries; they do not override the v1 primary class decision.

Any future change in which score level is primary or which threshold is used requires a new decision-policy version and separate validation. GSE280090 outcomes cannot be used to tune the v1 threshold post hoc.

## Result contract

Recommended normalized result:

```text
MarlinClassificationReport
├── schema_version
├── sample_id
├── status                       # existing ModuleRunStatus
├── classification_decision      # HIGH_CONFIDENCE | UNKNOWN | null
├── decision_policy
├── source_kind
├── input_schema
├── genome_build
├── artifact_lock_id
├── input_fingerprint
├── feature_summary
├── model_unit_scores[42]
├── current_class_scores
├── methylation_family_scores
├── lineage_scores
├── top_class
├── top_class_score
├── warnings
├── limitations
├── runtime_provenance
└── research_only = true
```

### Generic module status

- `COMPLETED` — valid model inference exists, regardless of high-confidence versus unknown decision;
- `NO_CALL` — inference is intentionally suppressed by a safety guard despite structurally valid input, e.g. zero observed model features;
- `NOT_RUN` — classifier not requested or intentionally unavailable before execution;
- `FAILED` — invalid input, build mismatch, artifact mismatch, runtime failure or malformed output.

Low confidence does **not** map to `NO_CALL`; it maps to `COMPLETED + UNKNOWN`.

### Zero-evidence guard

If zero locked model features are actually observed and the entire model vector would be zero, ONTSeq does not expose a prior-like network prediction as sample evidence.

It returns:

```text
status = NO_CALL
classification_decision = null
reason = NO_MODEL_FEATURES_OBSERVED
```

This is an ONTSeq safety boundary and receives an explicit regression test.

## Error handling

### Invalid input -> `FAILED`

Examples:

- wrong column count/delimiter;
- bad coordinates;
- beta outside `[0,1]` or non-finite;
- duplicate/empty probe ID;
- declared build mismatch.

No malformed rows are silently skipped.

### Artifact/runtime mismatch -> `FAILED`

Examples:

- changed model;
- changed feature/class/probe resource;
- incompatible runtime profile;
- incomplete score output.

### Valid sparse input -> `COMPLETED`

Published missing-value semantics remain intact. Sparsity is reported quantitatively. Result may be `HIGH_CONFIDENCE` or `UNKNOWN`.

## CLI/operator flow

```text
1. Install/prepare MARLIN runtime and artifacts explicitly/offline.
2. `ontseq marlin-lock ...` creates/verifies the artifact lock.
3. `ontseq marlin-features ...` inspects canonical feature construction.
4. `ontseq marlin-classify ...` classifies one sample.
5. `ontseq marlin-validate ...` runs a registered external processed-data validation manifest.
```

Classification always requires an explicit build and artifact lock. Build is not guessed from filename; model is not fetched implicitly.

## Validation strategy

### Phase 0 — pure preprocessing oracle

Synthetic tests independently define expected vectors for:

- row-order invariance;
- `0.5` boundary;
- `<0.5 -> -1`;
- `>=0.5 -> +1`;
- `NA -> 0`;
- absent feature -> 0;
- feature-order independence;
- duplicate rejection;
- malformed numeric rejection;
- gzip/plain equivalence;
- vector-digest determinism.

Expected vectors are hand/independently specified, not generated by production code.

### Phase 1 — fixed real-runtime compatibility oracle

A deterministic synthetic/reference feature vector is executed through the actual locked MARLIN runtime.

Freeze:

- feature-vector SHA-256;
- 42-score vector;
- score tolerance;
- top-class requirement;
- runtime/artifact lock identity.

Repeated runs must satisfy the frozen profile. This happens before external biological validation.

### Phase 2 — AL_001 external gate

Input: `GSM8587229_AL_001.txt.gz`.

Requirements:

- public accession and local SHA-256 recorded;
- GRCh37/hg19 declaration;
- deterministic parse/feature digest;
- deterministic 42-score output within frozen runtime tolerance;
- decision follows locked `class_name_current >= 0.8` policy;
- complete provenance in normalized JSON/report;
- no upstream Dorado/modkit claim.

Expected biological comparison must be registered from publication/source-data context before the ONTSeq result is inspected. Expected values may not be derived from the first ONTSeq execution.

### Phase 3 — GSE280090 cohort

Run every compatible restricted probe-level sample under the same lock/policy.

Predeclared metrics:

- attempted/completed/failed counts;
- high-confidence count;
- unknown count;
- top-class concordance where independent expected class exists;
- high-confidence concordance;
- no-call fraction;
- feature-coverage distribution;
- runtime distribution;
- deterministic-rerun failures;
- discordant-sample list.

No model, threshold or preprocessing change is allowed after outcome inspection within the same registered validation.

### Phase 4 — native ONTSeq bridge

Only after processed-data reproducibility:

```text
same specimen/equivalent source
published or locked precomputed probe input
                 versus
ONTSeq modBAM -> modkit 0.6.4 -> probe-level CpG extraction
```

Compare:

- observed probe IDs;
- pre-binarization methylation fractions;
- feature vector/digest;
- 42 raw scores;
- grouped class/family/lineage scores;
- top class/confidence decision.

This is the first evidence for native ONTSeq-to-MARLIN interoperability. It is still not clinical validation.

### Phase 5 — streaming/realtime, deferred

Streaming is a separate future design because early threshold crossing can reverse. It requires cumulative-evidence and stability policies; first-threshold-crossing alone is insufficient.

## Testing strategy

Implementation is strict TDD.

### Unit tests

Add dedicated tests for:

- contracts;
- probe-BED parser;
- gzip transport;
- feature mapping/binarization;
- artifact lock;
- 357,340-feature invariant;
- 42-softmax-score invariant;
- annotation grouping;
- `HIGH_CONFIDENCE`/`UNKNOWN` policy;
- zero-feature `NO_CALL` guard;
- deterministic hashing;
- serialization round trips.

### Property/metamorphic tests

At minimum:

- input row permutation preserves feature vector;
- plain/gzip semantic equivalence;
- result roundtrip preserves every score;
- any locked artifact-byte change invalidates preflight;
- build change invalidates run;
- repeated locked runs are numerically equivalent under the frozen runtime profile;
- rendering does not mutate normalized result.

Hypothesis may generate numeric/boundary cases where it gives an independent oracle.

### Integration tests

Actual R/Keras/TensorFlow runtime is exercised with the fixed synthetic/reference vector when installed.

CI has two layers:

1. always-blocking contract/preprocessing tests;
2. blocking MARLIN-runtime integration in the designated packaged-runtime workflow. Where the external runtime/model is intentionally absent, the job may skip only with an explicit reason and cannot be reported as validation success.

### Existing-suite protection

Before merge run all applicable repository checks:

- pytest;
- Ruff lint/format;
- mypy;
- schema/version consistency;
- repository-safety checks;
- relevant runtime/Snakemake dry runs;
- existing methylation suites;
- packaging/install smoke tests.

Existing biological output fixtures must not change as a side effect unless a separate independently justified defect fix is reviewed.

## Security/local execution

Third-party model/runtime content is treated as executable/untrusted input.

Requirements:

- no network required during inference;
- model/runtime installed only by explicit action;
- hashes verified before use;
- isolated subprocess working directory;
- output constrained to run workspace;
- minimized inherited environment;
- complete process-tree cancellation via existing hardened primitives;
- captured logs without exporting patient read identifiers;
- no BAM/POD5/FASTQ or processed patient payloads in GitHub artifacts.

## Packaging

Runtime/model packaging is separate from Python source packaging.

If redistribution terms permit, a portable bundle may include a verified runtime archive. Otherwise it includes an explicit installer/import workflow, source/version expectations and lock creation.

Operator-visible missing/error states remain distinct:

```text
runtime unavailable
model unavailable
artifact mismatch
input invalid
inference failed
inference completed unknown
inference completed high confidence
```

Optional MARLIN absence must not break unrelated ONTSeq profiles.

## Reporting/UI boundary

Scientific core is implemented and validated before UI integration.

First mergeable increment:

- contracts;
- parser;
- feature construction;
- artifact lock;
- runtime adapter;
- CLI;
- normalized JSON;
- validation runner.

UI/report code consumes only normalized MARLIN output.

Minimum UI summary later:

- MARLIN status;
- top class/score;
- `HIGH_CONFIDENCE` or `UNKNOWN`;
- observed model-CpG count/fraction;
- source kind;
- Research-Use-Only label;
- upstream-validation limitation.

All raw 42 scores remain available for reviewer detail.

## Versioning

This work belongs to the 0.9.0 development line but does not itself declare a 0.9.0 release.

- package/runtime/desktop/report version surfaces must remain internally consistent;
- MARLIN artifact-lock and result schemas have their own schema versions;
- changelog must state RUO/non-clinical status;
- release requires the broader 0.9.0 validation/audit gates, not merely a working classifier.

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

Focused edits are expected in the canonical CLI/dispatcher, changelog and version/schema guard inputs. UI/report files are a later reviewed increment.

Large model/data artifacts are never committed.

## Implementation sequence

1. contracts/enums with failing tests;
2. strict precomputed probe-BED parser;
3. canonical feature construction/digest;
4. artifact lock/preflight;
5. runtime adapter with fixed synthetic/reference vector;
6. annotation grouping and confidence policy;
7. CLI and normalized JSON;
8. validation manifest/runner;
9. AL_001 gate;
10. compatible GSE280090 cohort gate;
11. native modkit-to-MARLIN bridge;
12. UI/report consumption;
13. separate future streaming design.

AL_001 is not used as the debugging oracle for steps 1–8.

## Acceptance criteria

The first MARLIN classification increment is complete when:

- [ ] MARLIN v1 preprocessing exists in one canonical ONTSeq implementation;
- [ ] model, feature, class-annotation and hg19 probe artifacts are checksum-locked;
- [ ] v1 requires explicit GRCh37/hg19;
- [ ] inference is offline and rejects artifact mismatch;
- [ ] strict five-column probe-BED supports plain/gzip input;
- [ ] 357,340-feature invariant is enforced;
- [ ] row order cannot change the feature vector;
- [ ] `0.5`, missing/NA and duplicate semantics are regression-tested;
- [ ] zero observed model features produce explicit `NO_CALL` without model prediction;
- [ ] runtime returns exactly 42 finite `[0,1]` softmax scores summing to 1 within `1e-5`;
- [ ] all 42 raw scores are retained;
- [ ] class/family/lineage aggregation derives only from locked annotations;
- [ ] v1 primary decision is grouped `class_name_current >= 0.8 -> HIGH_CONFIDENCE`, otherwise `UNKNOWN`;
- [ ] threshold/policy is not tuned on GSE280090;
- [ ] synthetic real-runtime compatibility profile is frozen before external validation;
- [ ] `GSM8587229_AL_001.txt.gz` runs reproducibly with accession/build/hash/provenance;
- [ ] GSE280090 validation runs from an external-data manifest without committing payloads;
- [ ] validation separates technical completion, classification confidence and biological concordance;
- [ ] `PRECOMPUTED_METHYLATION` cannot be mistaken for native `MODKIT_DERIVED` validation;
- [ ] full existing ONTSeq regression/quality gates remain green;
- [ ] no clinical-reportable status is introduced.

## Risks and mitigations

### Runtime/version drift

Old Keras/TensorFlow HDF5 behavior may vary. Mitigation: full runtime lock plus frozen fixed-vector compatibility profile; preprocessing stays outside ML runtime.

### Build drift

MARLIN publishes hg19/hg38/T2T probe resources. Mitigation: v1 GRCh37-only lock; no implicit liftover.

### Overclaim from public processed data

Strong GSE concordance can be mistaken for raw-pipeline validation. Mitigation: mandatory source-kind provenance and explicit evidence-boundary wording.

### Sparse-input false certainty

Softmax always returns a distribution. Mitigation: preserve feature coverage, `UNKNOWN` below confidence threshold and `NO_CALL` when zero model features are observed.

### Numerical nondeterminism

Backend/threading can perturb low-order scores. Mitigation: pre-GSE frozen runtime compatibility profile, not post-hoc tolerance adjustment.

### Patient-derived data handling

GSE files are public processed human methylation data. Mitigation: no payloads in Git, accession/checksum manifests only, local paths excluded from public reports, existing ONTSeq governance retained.

## Scientific boundary after completion

Successful completion supports:

> ONTSeq reproducibly implements the published MARLIN v1 CpG-to-classification semantics against a pinned model/runtime and can evaluate public processed Nanopore methylation inputs with explicit provenance and confidence handling.

It does not support:

> ONTSeq clinically diagnoses acute leukemia from Nanopore sequencing.

It also does not yet support:

> ONTSeq's raw-signal/Dorado/modkit path is analytically equivalent to the MARLIN publication.

Those stronger claims require the native bridge, analytical performance characterization, intended-use cohorts, orthogonal truth, governance and human review.
