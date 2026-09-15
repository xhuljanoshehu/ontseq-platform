# MARLIN Dual-Runtime Compatibility Gate — Design

## Status

Approved architectural design for implementation on `feat/marlin-classification-v1`.

This design strengthens the existing frozen MARLIN runtime fixture. The current `marlin-freeze-runtime` path remains valid as a reproducibility baseline for one concrete runtime, but ONTSeq additionally needs an explicit cross-runtime comparison when the frozen baseline is intended to support claims about numerical compatibility with a separate reference environment.

## Goal

Establish a fail-closed, Research-Use-Only numerical compatibility gate that compares the same immutable MARLIN v1 artifact set and the same deterministic 357,340-value non-biological feature fixture across two independently identified execution runtimes:

1. a **reference runtime** representing the chosen MARLIN reference environment; and
2. a **candidate runtime** representing the ONTSeq runtime proposed for biological validation.

A PASS supports only the statement that the candidate runtime reproduced the reference runtime's locked MARLIN model behavior on the fixed non-biological fixture within predeclared engineering tolerances.

It does not establish analytical validity, clinical validity, GSE280090 concordance, native modkit equivalence, or raw-signal pipeline equivalence.

## Scientific baseline

The inspected upstream MARLIN repository documents the following tested software baseline:

- R 4.1.3;
- R Keras 2.13.0;
- R TensorFlow 2.13.0;
- compatible TensorFlow with Python 3.10.

The current ONTSeq Linux CPU engineering runtime resolves to:

- R 4.2.3;
- R Keras 2.13.0;
- R TensorFlow 2.13.0;
- Python 3.10.21;
- Python TensorFlow 2.13.0;
- CPU backend.

ONTSeq must not assume these environments are numerically equivalent merely because model loading and inference succeed.

## Core architectural correction

### Artifact-set identity and execution-runtime identity are separate trust domains

`MarlinArtifactLock` intentionally binds immutable classifier artifacts to one concrete execution runtime. That remains useful and is not removed. A dual-runtime experiment therefore uses **two execution locks**:

- a reference execution lock whose runtime fields describe the live reference environment;
- a candidate execution lock whose runtime fields describe the live candidate environment.

The two locks may have different `lock_id`, runtime versions, backend and creation timestamp. They may be compared only if their immutable classifier artifacts are exactly identical under a canonical artifact-set digest.

This avoids rewriting a lock to impersonate another runtime and preserves backward compatibility with all existing MARLIN v1 execution paths.

### Canonical artifact-set identity

Add a deterministic function and strict contract for immutable classifier identity:

```text
MarlinArtifactSetIdentity
schema_version = "0.1.0"
artifact_set_sha256
model_version
model_sha256
code_version_or_commit
code_manifest_sha256
feature_sha256
canonical_feature_list_sha256
class_annotation_sha256
probe_resource_sha256
genome_build
expected_feature_count = 357340
expected_model_unit_count = 42
preprocessing_contract_version
runtime_contract_version
research_only = true
```

`artifact_set_sha256` is calculated from a canonical byte representation of the fields above. It explicitly excludes:

- `lock_id`;
- R/Python/Keras/TensorFlow versions;
- execution backend;
- timestamps;
- local file paths;
- source URI formatting that does not change artifact bytes.

Both execution locks must yield exactly the same `artifact_set_sha256` before any candidate inference is allowed.

### Runtime identity contract

Add `MarlinRuntimeIdentity` with at least:

```text
schema_version = "0.1.0"
runtime_id
R_version
keras_version
tensorflow_version
python_version
execution_backend
runtime_contract_version = "marlin-runtime-v1"
created_at
research_only = true
```

The identity is populated only from the live runtime probe. Manually typed version strings are not accepted as evidence.

`MarlinRuntimeIdentity` must also be internally consistent with the corresponding `MarlinArtifactLock` runtime fields. A mismatch is a hard failure.

## Existing freeze behavior remains

`ontseq marlin-freeze-runtime` remains the canonical command to create a single-runtime frozen baseline. It continues to:

1. load one execution artifact lock;
2. probe the live runtime;
3. require the live runtime to match that lock's runtime identity;
4. generate the deterministic non-biological feature fixture;
5. execute the locked model once;
6. write the frozen fixture and 42-score runtime profile.

Therefore the first reference profile is frozen under the **reference execution lock**, not under the later candidate lock. No existing frozen profile is silently overwritten.

## New dual-runtime workflow

Introduce a separate command:

```text
ontseq marlin-compare-runtimes
```

The command compares a previously frozen reference profile with a live candidate runtime execution while proving that both execution locks represent the same immutable classifier artifact set.

Required inputs:

```text
--reference-artifact-lock
--candidate-artifact-lock
--reference-runtime-identity
--reference-profile
--runtime-fixture
--model
--inference-script
--candidate-runtime-probe-script
--candidate-rscript
--comparison-id
--output
```

The `--model` bytes must satisfy both execution locks' model SHA-256 because the locks are required to share the same artifact-set identity.

### Execution order

The command must execute in this exact logical order:

```text
load reference execution lock
-> load candidate execution lock
-> derive canonical artifact-set identity from each lock
-> require equal artifact-set SHA-256
-> load reference runtime identity
-> verify reference runtime identity matches reference execution lock
-> load reference profile
-> verify profile belongs to reference execution lock
-> verify fixed fixture SHA against profile
-> verify model SHA and immutable artifact identity
-> live-probe candidate runtime
-> build candidate runtime identity from probe
-> verify candidate runtime identity matches candidate execution lock
-> execute same frozen 357,340-value fixture through same locked model
-> validate exactly 42 finite softmax scores
-> compare all 42 candidate scores to reference scores
-> compare top model-unit identity
-> compare softmax-sum invariant
-> atomically write PASS/FAIL comparison report
```

Any artifact-set mismatch fails **before** candidate model execution.

A numerical mismatch is a failed compatibility comparison. It must not be converted to `UNKNOWN`, `NO_CALL`, or a biological classification state.

## Comparison contract

Add `MarlinDualRuntimeCompatibilityReport` with at least:

```text
schema_version = "0.1.0"
comparison_id
artifact_set_identity
reference_artifact_lock_id
candidate_artifact_lock_id
feature_vector_sha256
reference_profile_id
reference_runtime_identity
candidate_runtime_identity
absolute_score_tolerance
score_sum_tolerance
reference_scores[42]
candidate_scores[42]
absolute_differences[42]
max_absolute_difference
reference_top_model_unit_index
candidate_top_model_unit_index
all_scores_within_tolerance
top_model_unit_matches
softmax_invariants_pass
verdict = PASS | FAIL
created_at
research_only = true
```

The report stores no biological sample information.

### PASS semantics

PASS requires all of:

- reference profile belongs to the reference execution lock;
- reference and candidate execution locks yield the same canonical artifact-set SHA-256;
- candidate model bytes match the shared model SHA-256;
- fixture digest equals the reference profile fixture digest;
- reference runtime identity matches the reference execution lock;
- live candidate runtime identity matches the candidate execution lock;
- exactly 42 reference and 42 candidate scores exist;
- all values are finite probabilities;
- each score's absolute difference is `<=` the frozen absolute tolerance;
- candidate and reference top raw model-unit indices are identical;
- both score vectors satisfy the declared softmax-sum tolerance.

Any failed requirement yields overall `FAIL` or, for malformed/inconsistent evidence that prevents a valid comparison, a hard command failure. A report must never claim PASS after a failed precondition.

## Tolerance policy

The current engineering policy remains:

```text
absolute tolerance per raw score = 1e-7
softmax-sum tolerance            = 1e-5
```

These values are fixed before biological validation and may not be tuned after inspecting AL_001 or GSE280090 outcomes.

If the reference environment cannot reproduce the deterministic fixture stably enough to support the declared tolerance, that is an engineering finding. The tolerance must not be widened merely to force PASS without an explicit new policy/version and documented justification.

## Reference-runtime provenance

The dual-runtime report distinguishes:

- **reference runtime identity** — the live environment that generated the frozen 42-score oracle;
- **candidate runtime identity** — the live environment currently being qualified;
- **shared artifact-set identity** — the immutable model/code/features/annotations/probe resources common to both runs.

For the first compatibility study, the preferred reference runtime is the closest reproducible environment to upstream MARLIN's documented tested stack, centered on R 4.1.3 / Keras 2.13 / TensorFlow 2.13 / Python 3.10.

If exact reconstruction of that environment is impossible, ONTSeq must record the actually resolved reference runtime rather than falsely label it as the published environment. The report must state exact observed versions.

## Evidence hierarchy

The evidence ladder becomes:

```text
synthetic parser/feature contract tests
        <
real runtime import/backend smoke
        <
single-runtime frozen reproducibility baseline
        <
dual-runtime numerical compatibility
        <
AL_001 processed-CpG external validation
        <
GSE280090 cohort validation
        <
same-specimen native modkit bridge validation
        <
analytical intended-use validation
        <
clinical validation
```

No PASS may be promoted to a higher evidence level.

## Error handling

Fail closed on:

- reference/candidate artifact-set SHA mismatch;
- model SHA mismatch;
- reference profile belonging to a different reference execution lock;
- wrong fixture digest;
- missing or malformed runtime identity;
- reference runtime identity inconsistent with its lock;
- candidate runtime probe failure;
- candidate runtime identity inconsistent with its lock;
- wrong score count;
- non-finite score;
- score outside `[0,1]`;
- softmax sum outside tolerance;
- top-unit mismatch;
- any score drift greater than the fixed tolerance;
- partial output publication.

A comparison failure is evidence of runtime incompatibility under the tested conditions, not a leukemia classification outcome.

## Persistence and immutability

Comparison reports are local engineering evidence artifacts. They may be stored under a validation results directory but are not committed with model bytes, patient data, GSE payloads, credentials, or large runtime artifacts.

Output is written atomically. Existing comparison outputs are not silently overwritten unless a future explicit versioned overwrite policy is designed.

## Testing requirements

Implementation uses strict TDD.

Unit tests must cover:

- canonical artifact-set digest is identical for locks differing only in runtime identity/lock ID/timestamp;
- artifact-set digest changes when any immutable artifact SHA or build/contract field changes;
- runtime identity construction from a live probe;
- rejection of runtime identity inconsistent with its execution lock;
- exact preservation of reference and candidate runtime identities;
- artifact-set mismatch rejected before candidate model execution;
- fixture SHA mismatch rejected before candidate model execution;
- 42-score PASS within tolerance;
- one-score drift `>1e-7` -> FAIL;
- exact boundary drift `==1e-7` -> PASS;
- top-unit mismatch -> FAIL even if every numeric score individually satisfies tolerance;
- non-finite/invalid candidate score rejected;
- softmax invariant failure -> FAIL;
- deterministic serialization/round-trip;
- atomic output behavior;
- no patient/GSE payload dependency.

Integration tests must verify the CLI argument contract and a fake-runner end-to-end comparison using two execution locks with the same artifact-set identity but different runtime identities.

The real reference-vs-candidate model execution remains an external engineering qualification step because the trained MARLIN model is not stored in GitHub CI.

## CI and release boundaries

Repository CI must continue to pass:

- Python 3.11/3.12/3.13 Core CI;
- Ruff;
- mypy;
- repository safety checks;
- Desktop tests/build;
- Desktop bundle contract;
- MARLIN runtime smoke;
- CodeQL where configured.

The dual-runtime gate does not bump ONTSeq from 0.8.2 by itself.

PR #77 remains Draft. No merge to `main`, release, clinical claim, or publication occurs without explicit user approval.

## Acceptance criteria

The feature is complete when ONTSeq can produce a machine-readable comparison report proving, for two separately runtime-locked executions of one identical classifier artifact set and one fixed non-biological fixture, whether the live candidate runtime reproduces the frozen reference runtime's 42 MARLIN scores within the predeclared tolerance while preserving both runtime identities separately.

Only after that report is PASS should the project proceed to AL_001 and GSE280090 biological downstream validation under that candidate runtime.
