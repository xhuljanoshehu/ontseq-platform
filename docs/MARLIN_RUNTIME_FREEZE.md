# MARLIN frozen runtime compatibility fixture

## Status

**Research Use Only. Engineering compatibility gate only.**

This document describes how ONTSeq creates the fixed, non-biological MARLIN runtime fixture used
before any external biological validation. Creating or passing this fixture does not establish
analytical or clinical validity.

## Purpose

The freeze separates software/runtime reproducibility from biological performance. A fixed
357,340-value input is passed through the locked MARLIN v1 model once. ONTSeq stores the resulting
42 raw softmax scores in a `MarlinRuntimeCompatibilityProfile`. Before AL_001 or any GSE280090
sample is executed by the external-validation harness, the same fixture is re-run and every score
must remain within the predeclared tolerance.

This detects runtime/model drift without using patient-derived outcomes to choose the tolerance.

## Fixed generator contract

Generator version:

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

The generated vector is therefore independent of biological data, MARLIN class labels and model
predictions.

Regression-locked identities for this generator are:

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

The text representation is one integer (`-1`, `0`, or `1`) per line, LF terminated, with a final
newline.

## Frozen tolerances

The engineering tolerances are declared before biological validation:

```text
absolute tolerance per raw model score: 1e-7
softmax-sum tolerance:                  1e-5
```

They are part of the technical policy and must not be changed after inspecting AL_001 or GSE280090
results to improve concordance.

## Live runtime identity gate

`marlin-freeze-runtime` does not trust the runtime identity written in the artifact lock without
rechecking it. The command first executes `scripts/marlin_runtime_probe.R` using the same `Rscript`
that will execute model inference and requires exact agreement with the artifact lock for:

- R version;
- R-Keras version;
- Python TensorFlow version;
- full Python patch version;
- CPU/GPU backend.

A mismatch is a hard failure before the frozen model profile is written.

## Command

The trained model and artifact lock remain local controlled resources. They are not downloaded or
committed by this command.

```bash
ontseq marlin-freeze-runtime \
  --artifact-lock results/marlin-validation/marlin-artifact-lock.json \
  --model /path/to/marlin_v1.model.hdf5 \
  --inference-script scripts/marlin_infer_locked.R \
  --runtime-probe-script scripts/marlin_runtime_probe.R \
  --profile-id MARLIN_V1_CPU_FROZEN \
  --rscript /path/to/ontseq-marlin-runtime/bin/Rscript \
  --output-fixture results/marlin-validation/marlin-runtime-fixture.txt \
  --output-profile results/marlin-validation/marlin-runtime-profile.json
```

The two output paths must be different and must not already exist. ONTSeq refuses silent overwrite
of a previously frozen reference profile.

## Outputs

`marlin-runtime-fixture.txt` contains only the deterministic non-biological `-1/0/+1` vector.

`marlin-runtime-profile.json` binds:

- artifact/runtime lock ID;
- feature-vector SHA-256;
- all 42 reference softmax scores;
- top raw model-unit index;
- execution backend;
- absolute per-score tolerance;
- score-sum tolerance;
- timezone-aware creation timestamp;
- `research_only=true`.

The generated fixture/profile are local validation artifacts and are not committed to the
repository.

## What PASS means

A subsequent `ontseq marlin-validate` run must receive both `--runtime-profile` and
`--runtime-fixture`. It re-runs the locked fixture before the first biological sample.

A compatibility PASS supports only this statement:

> The candidate MARLIN runtime reproduced the pre-frozen numerical behavior of the locked model
> on the fixed non-biological ONTSeq fixture within the declared engineering tolerances.

It does **not** demonstrate correct leukemia classification, validate GSE280090, validate the
native modkit bridge, or establish clinical performance.
