# MARLIN Dual-Runtime Compatibility Gate — Design

## Status

Approved architectural design for implementation on `feat/marlin-classification-v1`.

This design strengthens the existing frozen MARLIN runtime fixture. The current `marlin-freeze-runtime` path remains valid as a reproducibility baseline for one concrete runtime, but ONTSeq additionally needs an explicit cross-runtime comparison when the frozen baseline is intended to support claims about numerical compatibility with a separate reference environment.

## Goal

Establish a fail-closed, Research-Use-Only numerical compatibility gate that compares the same locked MARLIN v1 model and the same deterministic 357,340-value non-biological feature fixture across two independently identified runtimes:

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

### Artifact identity and runtime identity are separate trust domains

`MarlinArtifactLock` currently records both immutable classifier artifacts and one concrete runtime identity. That is useful for fail-closed execution of a specific locked environment, but it is insufficient as the sole provenance model for a dual-runtime experiment because both runtimes must consume the same artifact set while preserving their different runtime identities.

The dual-runtime gate therefore introduces an explicit runtime identity object rather than mutating the artifact identities between executions.

### New runtime identity contract

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

### Artifact set identity

The model-side identity remains governed by the existing artifact lock and must be identical for both reference and candidate executions:

- model SHA-256;
- upstream MARLIN code identity/manifest SHA-256;
- feature RData SHA-256;
- canonical ordered feature-list SHA-256;
- class-annotation SHA-256;
- hg19 probe-resource SHA-256;
- GRCh37 build;
- preprocessing contract;
- expected 357,340 inputs;
- expected 42 outputs.

The existing runtime fields in `MarlinArtifactLock` remain supported for backward compatibility and for the normal single-runtime execution boundary in MARLIN v1. The dual-runtime comparison must not rewrite those fields to impersonate another runtime.

## Existing freeze behavior remains

`ontseq marlin-freeze-runtime` remains the canonical command to create a single-runtime frozen baseline. It continues to:

1. load the artifact lock;
2. probe the live runtime;
3. require the live runtime to match the runtime identity recorded by that execution lock;
4. generate the deterministic non-biological feature fixture;
5. execute the locked model once;
6. write the frozen fixture and 42-score runtime profile.

No existing frozen profile is silently overwritten.

## New dual-runtime workflow

Introduce a separate command:

```text
ontseq marlin-compare-runtimes
```

The command compares a previously frozen reference profile with a candidate runtime execution.

Required inputs:

```text
--artifact-lock
--reference-runtime-identity
--reference-profile
--runtime-fixture
--model
--inference-script
--candidate-runtime-probe-script
--candidate-rscript
--output
```

Optional explicit candidate identity output may be supported, but the comparison report itself must contain the live-probed candidate identity.

### Execution order

The command must execute in this exact logical order:

```text
load artifact lock
-> load reference runtime identity
-> load reference profile
-> verify profile belongs to artifact lock
-> verify fixed fixture SHA against profile
-> verify model SHA and locked artifacts
-> live-probe candidate runtime
-> build candidate runtime identity from probe
-> execute same frozen 357,340-value fixture through same locked model
-> validate exactly 42 finite softmax scores
-> compare all 42 candidate scores to reference scores
-> compare top model-unit identity
-> compare softmax-sum invariant
-> atomically write PASS/FAIL comparison report
```

A mismatch is a failed compatibility comparison. It must not be converted to `UNKNOWN`, `NO_CALL`, or a biological classification state.

## Comparison contract

Add `MarlinDualRuntimeCompatibilityReport` with at least:

```text
schema_version = "0.1.0"
comparison_id
artifact_lock_id
model_sha256
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

- artifact lock ID matches the reference profile;
- candidate uses the same model bytes verified by artifact SHA-256;
- fixture digest equals the reference profile fixture digest;
- exactly 42 reference and 42 candidate scores exist;
- all values are finite probabilities;
- each score's absolute difference is `<=` the frozen absolute tolerance;
- candidate and reference top raw model-unit indices are identical;
- both score vectors satisfy the declared softmax-sum tolerance;
- the candidate runtime identity is live-probed and preserved in the report.

Any failed requirement yields overall `FAIL`.

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

- **reference runtime identity** — the environment that generated the frozen 42-score oracle;
- **candidate runtime identity** — the environment currently being qualified.

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

- artifact/model SHA mismatch;
- reference profile belonging to a different artifact lock;
- wrong fixture digest;
- missing or malformed runtime identity;
- runtime probe failure;
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

- runtime identity construction from a live probe;
- rejection of manually inconsistent runtime identity data;
- exact preservation of reference and candidate runtime identities;
- same artifact/model identity across both executions;
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

Integration tests must verify the CLI argument contract and a fake-runner end-to-end comparison.

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

The feature is complete when ONTSeq can produce a machine-readable comparison report proving, for one locked artifact set and one fixed non-biological fixture, whether a live candidate runtime reproduces the frozen reference runtime's 42 MARLIN scores within the predeclared tolerance while preserving both runtime identities separately.

Only after that report is PASS should the project proceed to AL_001 and GSE280090 biological downstream validation under that candidate runtime.
