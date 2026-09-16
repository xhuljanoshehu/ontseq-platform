# MARLIN frozen and dual-runtime compatibility gates

## Status

**Research Use Only. Engineering compatibility evidence only.**

These gates separate software/runtime reproducibility from biological performance. A PASS does not establish analytical validity, clinical validity, GSE280090 concordance, native modkit equivalence, or raw-signal equivalence.

## Fixed non-biological fixture

ONTSeq generates one deterministic 357,340-value input with generator version:

```text
sha256-index-mod3-v1
```

For zero-based feature index `i` from `0` through `357339`:

```text
digest  = SHA256(b"ONTSeq-MARLIN-runtime-fixture-v1\0" + i.to_bytes(8, "big"))
residue = digest[0] % 3
0 -> -1
1 ->  0
2 -> +1
```

Regression-locked identities:

```text
feature count:          357340
-1 values:              119915
 0 values:              118324
+1 values:              119101
non-zero values:        239016
feature-vector SHA256:  b7f7c11d52f777d5e9c45f5d61912d42c21fb2a6dc403a7303d8eb0bf5140cde
fixture-text bytes:     834595
fixture-text SHA256:    c41b18446193c5ca93e64e054762356203a020b9a53d9330b31e979d61b27d30
```

The fixture is explicitly non-biological and independent of model outcomes.

## Frozen engineering tolerances

```text
absolute tolerance per raw model score: 1e-7
softmax-sum tolerance:                  1e-5
```

These values are declared before biological validation and must not be widened after inspecting AL_001 or GSE280090 merely to improve concordance.

## Single-runtime freeze

`ontseq marlin-freeze-runtime` creates the frozen reference oracle for one execution lock. The command:

1. loads the `MarlinArtifactLock`;
2. live-probes the same `Rscript` environment used for inference;
3. requires exact agreement with the lock for R, R-Keras, TensorFlow, full Python patch version and CPU/GPU backend;
4. generates the deterministic fixture;
5. runs the locked model once and preserves all 42 softmax scores;
6. writes the fixture and `MarlinRuntimeCompatibilityProfile` atomically;
7. when requested, writes a `MarlinRuntimeIdentity` produced from **that same live probe and same freeze timestamp**.

The optional runtime-identity output is the provenance source required by the dual-runtime gate. It must not be recreated by manually typing version strings.

Example:

```bash
ontseq marlin-freeze-runtime \
  --artifact-lock results/marlin-validation/reference-artifact-lock.json \
  --model /controlled/marlin_v1.model.hdf5 \
  --inference-script scripts/marlin_infer_locked.R \
  --runtime-probe-script scripts/marlin_runtime_probe.R \
  --profile-id MARLIN_V1_REFERENCE_FROZEN \
  --rscript /reference-runtime/bin/Rscript \
  --output-fixture results/marlin-validation/marlin-runtime-fixture.txt \
  --output-profile results/marlin-validation/marlin-runtime-profile.json \
  --output-runtime-identity results/marlin-validation/reference-runtime-identity.json
```

All requested output paths must be distinct and must not already exist. Partial publication is cleaned up on failure.

## Why a dual-runtime gate is required

A single-runtime freeze is excellent for detecting later drift in that same locked environment, but it does not by itself prove numerical equivalence between two different environments.

For MARLIN v1, the preferred reference environment is the closest reproducible environment to the upstream documented tested stack, centered on:

```text
R                  4.1.3
R-Keras            2.13.0
R-TensorFlow       2.13.0
Python             3.10 compatible
```

The current ONTSeq engineering candidate resolves to:

```text
R                  4.2.3
R-Keras            2.13.0
R-TensorFlow       2.13.0
Python             3.10.21
Python TensorFlow  2.13.0
backend            CPU
```

ONTSeq does not assume these runtimes are numerically equivalent merely because both can load the model.

## Artifact-set identity versus execution-runtime identity

The comparison deliberately uses **two execution locks**:

- a reference `MarlinArtifactLock` created from the live reference runtime;
- a candidate `MarlinArtifactLock` created from the live candidate runtime.

They may differ in lock ID, runtime versions, backend, source-URI formatting and timestamp. Before candidate inference, ONTSeq derives a canonical immutable artifact-set identity from both locks and requires the SHA-256 to be identical.

The artifact-set digest binds:

- model version and model SHA-256;
- MARLIN code version/commit and code-manifest SHA-256;
- feature RData SHA-256;
- canonical ordered feature-list SHA-256;
- class-annotation SHA-256;
- hg19 probe-resource SHA-256;
- GRCh37 build;
- expected 357,340 input features and 42 model outputs;
- preprocessing and runtime contract versions.

It excludes runtime versions/backend, lock IDs, timestamps, local paths and source-URI formatting.

## Required operator sequence

```text
A. Reconstruct the reference MARLIN environment and live-probe it.
B. Create the reference execution lock with marlin-lock.
C. Freeze fixture + reference profile + reference runtime identity in the reference runtime.
D. Create the candidate execution lock from the candidate live probe using identical artifact bytes.
E. Run marlin-compare-runtimes with both locks and the frozen reference evidence.
F. Require report verdict PASS before AL_001 or GSE280090 is executed.
```

If the reconstructed reference environment resolves to versions different from the upstream documentation, record the **actual observed versions**. Never relabel a different runtime as R 4.1.3.

## Dual-runtime comparison command

```bash
ontseq marlin-compare-runtimes \
  --reference-artifact-lock results/marlin-validation/reference-artifact-lock.json \
  --candidate-artifact-lock results/marlin-validation/candidate-artifact-lock.json \
  --reference-runtime-identity results/marlin-validation/reference-runtime-identity.json \
  --reference-profile results/marlin-validation/marlin-runtime-profile.json \
  --runtime-fixture results/marlin-validation/marlin-runtime-fixture.txt \
  --model /controlled/marlin_v1.model.hdf5 \
  --inference-script scripts/marlin_infer_locked.R \
  --candidate-runtime-probe-script scripts/marlin_runtime_probe.R \
  --candidate-rscript /candidate-runtime/bin/Rscript \
  --comparison-id MARLIN_V1_REFERENCE_VS_ONTSEQ_CPU \
  --output results/marlin-validation/runtime-comparison.json
```

Logical gate order:

```text
load both execution locks
-> derive both canonical artifact-set identities
-> require equal artifact-set SHA-256
-> verify frozen reference runtime identity against reference lock
-> verify reference profile belongs to reference lock
-> load fixture and require frozen fixture digest
-> verify model SHA-256
-> live-probe candidate runtime
-> require candidate probe identity to match candidate lock
-> run the same fixture through the candidate runtime
-> require exactly 42 finite probability scores
-> compare all 42 scores, top model unit and softmax invariant
-> atomically write PASS/FAIL report
```

Artifact-set or fixture mismatch fails **before candidate model inference**.

A structurally valid numerical mismatch generates a machine-readable `FAIL` report. Malformed or inconsistent evidence is a hard command failure; it is never converted into `UNKNOWN` or `NO_CALL`.

## PASS semantics

Dual-runtime `PASS` requires all of:

- same canonical immutable artifact-set identity;
- same frozen feature-vector SHA-256;
- reference profile bound to the reference execution lock;
- reference runtime identity consistent with the reference lock;
- live candidate runtime identity consistent with the candidate lock;
- exactly 42 valid scores in both vectors;
- every absolute score difference `<= 1e-7`;
- identical top raw model-unit index;
- both softmax vectors within the declared `1e-5` sum tolerance.

The report preserves both exact runtime identities, both execution-lock IDs, the common artifact-set identity, all 42 reference/candidate scores and differences, the maximum absolute difference, and the final PASS/FAIL verdict.

## Evidence boundary

A dual-runtime PASS supports only:

> Under the tested conditions, the live candidate runtime reproduced the frozen reference runtime's numerical MARLIN v1 behavior on the fixed non-biological fixture within the predeclared engineering tolerances.

It does **not** demonstrate correct leukemia classification. After PASS, the evidence ladder continues with:

```text
AL_001 processed-CpG external gate
-> GSE280090 cohort validation
-> same-specimen native modkit bridge validation
-> analytical intended-use validation
-> clinical validation
```

## Data and repository policy

The trained model, generated runtime profile, runtime identities, comparison reports, GEO payloads and patient-derived data are controlled local evidence and are not committed to Git. Repository safety rejects force-tracked files under `results/marlin-validation/`.
