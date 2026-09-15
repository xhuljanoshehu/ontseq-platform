# MARLIN classification in ONTSeq

## Status

**Research Use Only. No clinical release.**

This document describes the engineering integration of the published MARLIN v1 acute-leukemia
methylation classifier into ONTSeq. It does not claim analytical or clinical validity and it does
not claim that ONTSeq reproduces the complete raw-signal workflow used in the publication.

Current implementation scope:

- strict import of the published five-column probe-level methylation representation;
- explicit GRCh37/hg19-only v1 build contract;
- deterministic construction of the 357,340-value MARLIN v1 feature vector;
- checksum-locked model, feature, class-annotation and probe resources;
- locked R/Keras/TensorFlow inference boundary returning exactly 42 softmax scores;
- live runtime probing of R, R-Keras, Python TensorFlow, Python and CPU/GPU backend identity;
- deterministic non-biological runtime-fixture generation and one-time reference-profile freeze;
- grouped current-class, methylation-family and lineage summaries with explicit deterministic
  score/label tie ordering;
- published `0.8` high-confidence threshold with explicit `UNKNOWN` below threshold;
- frozen non-biological runtime-fixture verification before external cohort validation;
- manifest-driven external processed-CpG validation with deterministic repeat runs;
- a native modkit-to-MARLIN bridge contract reproducing MARLIN's combined 5mC+5hmC pileup
  semantics, while remaining unusable without a separately created same-specimen bridge-evidence
  lock.

The external GSE280090 / `GSM8587229_AL_001.txt.gz` validation gate has **not** been completed in
this repository state. Public/patient-derived validation payloads are not committed to Git.

The runtime-freeze implementation is present, but this repository does **not** contain a generated
compatibility profile from the real Zenodo MARLIN model. That profile remains a local controlled
validation artifact and must be created before biological validation begins.

## Scientific reference

Primary publication:

- Steinicke TL, Benfatto S et al. *Rapid epigenomic classification of acute leukemia.*
  Nature Genetics (2025). DOI: `10.1038/s41588-025-02321-z`.

Reference implementation inspected for this integration:

- repository: `https://github.com/hovestadt/MARLIN`
- pinned inspected commit: `442aa603415a54f62e7367794f9a31c6bc20fc2d`
- trained model: Zenodo record `15565404`, MARLIN model v1.0.0
- public training data: Zenodo record `15566584`
- processed nanopore validation series: GEO `GSE280090`

The published prediction semantics reproduced by ONTSeq are:

```text
ordered MARLIN v1 CpG feature present with beta >= 0.5  -> +1
ordered MARLIN v1 CpG feature present with beta <  0.5  -> -1
explicit NA or feature not observed                     ->  0
```

MARLIN v1 has exactly **357,340 input features** and **42 neural-network output units**.

## Evidence ladder

ONTSeq deliberately separates engineering evidence from analytical and clinical evidence:

```text
synthetic parser / contract / feature tests
        <
real R/Keras/TensorFlow import and backend smoke
        <
locked non-biological runtime compatibility
        <
GSE processed-CpG downstream reproduction
        <
same-specimen native modkit bridge comparison
        <
analytical intended-use validation
        <
clinical validation
```

A PASS at one level must never be relabelled as evidence for a higher level.

In particular:

- GSE280090 contains processed methylation data, not the original patient POD5/modBAM/full BAM
  material needed to validate ONTSeq's upstream raw-signal path.
- A successful GSE run can support only the statement that ONTSeq reproducibly executes the
  locked downstream CpG-to-MARLIN classification workflow.
- It does **not** validate Dorado basecalling, alignment, MM/ML generation, ONTSeq modkit
  extraction, CNV, SV, fusion calling, or clinical reportability.
- A successful runtime smoke proves that the selected R/Keras/TensorFlow stack imports and
  initializes. It does **not** prove numerical equivalence to the published MARLIN runtime or
  correct classification of biological samples.

## Architecture

The integration deliberately separates ONTSeq-owned preprocessing from MARLIN model inference.

```text
PRECOMPUTED_METHYLATION
five-column MARLIN probe table
          |
          v
strict ONTSeq parser
          |
          v
locked 357,340-feature construction
          |
          +-------------------------------+
                                          v
                               locked MARLIN runtime
                                          |
                                          v
                                  42 softmax scores
                                          |
                                          v
                         ONTSeq class/family/lineage
                         confidence + RUO report
```

A future native path is separate:

```text
MODKIT_DERIVED
modBAM
  |
  v
modkit 0.6.4 pileup
--modified-bases 5mC 5hmC
--combine-mods
  |
  v
combined cytosine bedMethyl (raw code C)
  |
  v
locked hg19 MARLIN probe resource
  |
  v
per-probe beta = sum(N_mod) / sum(N_valid)
  |
  v
same 357,340-feature constructor
```

Under this MARLIN-specific bridge, `N_mod` is the combined 5mC+5hmC modified-cytosine count,
`N_canonical` is the unmodified cytosine count and `N_other_mod` must be zero. Separate `m`/`h`
bedMethyl rows are refused by this bridge because they do not prove the published MARLIN
`--combine-mods` semantics.

This rule is **specific to the MARLIN bridge**. ONTSeq's general methylation lane continues to
preserve 5mC and 5hmC as separate modification codes and is not reinterpreted by this integration.

The native path is not enabled by the existence of code alone. It requires a
`MarlinBridgeLock` recording separate same-specimen comparison evidence. No valid bridge lock is
shipped as a default configuration.

## Input contract: published processed probe data

Adapter identity: `marlin_probe_bed_v1`.

The accepted logical schema is exactly five tab-separated fields:

```text
chromosome    start    end    methylation_fraction_or_NA    probe_id
```

Rules are fail-closed:

- v1 accepts only GRCh37/hg19;
- no inferred delimiter or build;
- no extra columns;
- malformed rows are not skipped;
- methylation values must be finite and within `[0,1]`, or literal `NA`;
- probe IDs must be non-empty and unique;
- input bytes are SHA-256 fingerprinted before parsing and rechecked afterwards;
- duplicate probe identifiers are rejected rather than resolved by row order.

## Feature construction

Feature order is controlled only by the locked canonical export of
`marlin_v1.features.RData`.

For each locked feature:

- observed beta `>= 0.5` -> `+1.0`;
- observed beta `< 0.5` -> `-1.0`;
- explicit `NA` -> `0.0`;
- absent feature -> `0.0`.

The feature summary records:

- expected feature count (`357340`);
- observed model-feature count;
- explicit-NA count;
- absent-feature count;
- non-model probe count;
- observed fraction;
- canonical feature-vector SHA-256;
- locked feature-artifact SHA-256.

The vector digest uses an explicit one-byte representation for `-1`, `0` and `+1`, avoiding
platform-dependent floating-point text formatting.

## Artifact lock

`MarlinArtifactLock` binds at least:

- trained model identity and SHA-256;
- upstream MARLIN code identity;
- upstream feature RData SHA-256;
- canonical exported ordered feature-list SHA-256;
- class-annotation SHA-256;
- hg19 probe-resource SHA-256;
- GRCh37 build;
- expected 357,340 features and 42 model units;
- preprocessing/runtime contract versions;
- live-probed R, R-Keras, Python TensorFlow and Python versions;
- live-probed execution backend;
- creation timestamp.

Inference performs no automatic artifact download. Hash mismatch is a hard failure.

`marlin-lock` no longer accepts manually typed runtime-version/backend values. It requires the
runtime probe script and the `Rscript` executable that will identify the actual candidate runtime.
This prevents a lock from claiming one software stack while executing another.

## Engineering runtime and live probe

ONTSeq deliberately keeps the control plane and the legacy MARLIN model runtime in separate Python
environments:

```text
ONTSeq control plane
Python >= 3.11
        |
        | argv-only Rscript execution
        v
MARLIN model runtime
R 4.2.x
Python 3.10.x
R-Keras 2.13.0
R-TensorFlow 2.13.0
Python tensorflow-cpu 2.13.0
```

The environment definition is `workflow/envs/marlin_runtime.yaml`. The CI smoke does not install
ONTSeq/Pydantic into the MARLIN Python 3.10 environment; this avoids contaminating the legacy
TensorFlow dependency set. ONTSeq invokes the MARLIN `Rscript` from its supported Python control
plane.

`scripts/marlin_runtime_probe.R` then:

1. locates the Python executable belonging to the same MARLIN environment as `Rscript`;
2. binds `reticulate` explicitly to that Python interpreter;
3. loads R-Keras and R-TensorFlow;
4. verifies that Keras is using the TensorFlow backend;
5. initializes Python TensorFlow with a real tensor operation;
6. queries CPU/GPU physical-device availability;
7. emits the observed R, R-Keras, Python TensorFlow and full Python patch versions.

The current Linux CPU engineering smoke resolves and verifies:

```text
R                 4.2.3
R-Keras           2.13.0
R-TensorFlow      2.13.0
Python            3.10.21
Python TensorFlow 2.13.0
backend            CPU
```

This is an **engineering runtime identity**, not a biological validation result. The inspected
upstream MARLIN repository documents R 4.1.3 as its tested R version. ONTSeq therefore does not
assume that R 4.2.3 is numerically equivalent. Numerical compatibility with the locked trained
model must be established independently by the frozen non-biological runtime fixture before any
biological validation is permitted.

The smoke workflow stores only a small non-biological runtime-identity JSON artifact. It does not
contain the trained MARLIN model, GSE280090 payloads or patient-derived data.

## Runtime compatibility gate

External biological validation requires more than a successful runtime import or a stored profile.
Before any validation sample is run, ONTSeq re-executes a fixed, non-biological 357,340-value
`-1/0/+1` feature fixture through the current locked model/runtime.

The canonical ONTSeq generator is `sha256-index-mod3-v1`. It deterministically maps each feature
index to `-1`, `0`, or `+1` under a fixed SHA-256 salt. The generated vector and text
representation have regression-locked SHA-256 identities. Full generator details and the locked
hashes are documented in `docs/MARLIN_RUNTIME_FREEZE.md`.

The engineering tolerances are fixed before biological validation:

```text
absolute tolerance per raw model score: 1e-7
softmax-sum tolerance:                  1e-5
```

`ontseq marlin-freeze-runtime` first probes the same live `Rscript` environment that will run
model inference and requires exact agreement with the artifact lock for R, R-Keras, Python
TensorFlow, full Python patch version and CPU/GPU backend. It then runs the fixed vector exactly
once through the locked model and writes the local fixture plus the 42-score compatibility
profile. Existing output paths are refused rather than overwritten.

The candidate runtime later must match the frozen `MarlinRuntimeCompatibilityProfile` for:

- artifact/runtime lock identity;
- execution backend;
- feature-vector digest;
- top model-unit identity;
- all 42 scores within the predeclared absolute tolerance;
- softmax-sum tolerance.

The tolerance is frozen before biological validation. It must not be tuned from AL_001 or
GSE280090 outcomes.

`ontseq marlin-validate` therefore requires an explicit `--runtime-fixture` path. A direct
`marlin-classify` invocation is not a substitute for this external-validation gate.

## Classification semantics

The 42 raw model-unit scores are retained. ONTSeq then uses the locked class annotation resource
to aggregate scores into:

- current methylation class;
- methylation family;
- lineage.

Grouped scores are ordered deterministically by `(-score, label)`. Therefore an exact score tie
is resolved by the lexicographically smaller locked label rather than by model-unit or workbook
row order. The normalized report validator enforces the same rule.

The primary classification decision uses the top grouped current-class score.

```text
top grouped current-class score >= 0.8 -> HIGH_CONFIDENCE
top grouped current-class score <  0.8 -> UNKNOWN
```

A valid low-confidence model result is therefore:

```text
ModuleRunStatus.COMPLETED + MarlinClassificationDecision.UNKNOWN
```

It is **not** `NO_CALL`.

`NO_CALL` is reserved for technically valid but non-interpretable input, currently including zero
observed MARLIN model features. Malformed input, lock mismatch, wrong runtime shape, non-finite
scores or partial output are failures, not biological negatives.

## Native modkit bridge

The bridge follows the published MARLIN probe aggregation rather than reusing ONTSeq's
chromosome/target-level methylation summaries.

The inspected upstream MARLIN workflow runs modkit over the hg19 probe resource with combined
modification semantics and subsequently reads genomic position, `N_valid` and `N_mod` to calculate
for each probe:

```text
beta = sum(N_mod) / sum(N_valid)
```

ONTSeq mirrors that probe-level rule with a dedicated MARLIN parser. It deliberately does **not**
reuse the generic ONTSeq 5mC-only parser, because the general methylation lane keeps 5mC and 5hmC
separate while MARLIN's native workflow combines them before computing beta.

Accepted MARLIN bridge rows therefore require:

```text
raw modification code = C
N_valid = N_mod + N_canonical + N_other_mod
N_other_mod = 0
```

A probe with no valid calls remains `NA` and becomes `0` only at MARLIN feature construction.

The bridge is gated by `MarlinBridgeLock`, which records:

- GRCh37/hg19 identity;
- adapter identity;
- modkit version `0.6.4`;
- fixed `5mC+5hmC-combine-mods-v1` pileup semantics;
- probe-resource and feature-resource fingerprints;
- same-specimen precomputed and native validation input fingerprints;
- precomputed and native feature-vector fingerprints;
- number of compared/concordant features and agreement fraction;
- evidence reference and timestamp.

No default bridge lock is shipped. Therefore the presence of the adapter does not by itself enable
or validate native modkit classification.

## CLI surface

Available engineering commands:

```text
ontseq marlin-lock
ontseq marlin-runtime-probe
ontseq marlin-freeze-runtime
ontseq marlin-features
ontseq marlin-classify
ontseq marlin-validate
```

The runtime probe can be executed independently:

```bash
ontseq marlin-runtime-probe \
  --probe-script scripts/marlin_runtime_probe.R \
  --rscript /path/to/ontseq-marlin-runtime/bin/Rscript \
  --output results/marlin-runtime-probe.json
```

`marlin-lock` requires the same live probe boundary:

```bash
ontseq marlin-lock \
  --model /path/to/marlin_v1.model.hdf5 \
  --feature-rdata /path/to/marlin_v1.features.RData \
  --feature-list /path/to/marlin_v1.features.txt \
  --class-annotations /path/to/marlin_v1.class_annotations.xlsx \
  --probe-bed /path/to/marlin_v1.probes_hg19.bed.gz \
  --lock-id MARLIN_V1_GRCH37_LOCAL \
  --code-version 442aa603415a54f62e7367794f9a31c6bc20fc2d \
  --code-manifest-sha256 <sha256> \
  --genome-build GRCh37 \
  --runtime-probe-script scripts/marlin_runtime_probe.R \
  --rscript /path/to/ontseq-marlin-runtime/bin/Rscript \
  --output results/marlin-artifact-lock.json
```

After the artifact lock has been created and before any biological validation, freeze the local
runtime compatibility reference:

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

The fixture/profile files are generated locally and are not repository assets.

### External validation

The validation runner requires a checksummed manifest, a locked artifact/runtime identity and the
frozen runtime fixture:

```bash
ontseq marlin-validate \
  --manifest results/marlin-validation/gse280090-v1/manifest.json \
  --artifact-lock results/marlin-validation/marlin-artifact-lock.json \
  --runtime-profile results/marlin-validation/marlin-runtime-profile.json \
  --runtime-fixture results/marlin-validation/marlin-runtime-fixture.txt \
  --model /path/to/marlin_v1.model.hdf5 \
  --feature-rdata /path/to/marlin_v1.features.RData \
  --feature-list /path/to/marlin_v1.features.txt \
  --class-annotations /path/to/marlin_v1.class_annotations.xlsx \
  --probe-bed /path/to/marlin_v1.probes_hg19.bed.gz \
  --inference-script /path/to/marlin_infer_locked.R \
  --output results/marlin-validation/gse280090-v1/report.json
```

Paths above are examples. Artifact identities are accepted only when they match their locks.

## First external validation target

The first registered target is intended to be:

```text
GEO series: GSE280090
accession:  GSM8587229
sample:     AL_001
file:       GSM8587229_AL_001.txt.gz
build:      hg19 / GRCh37
```

Before classification:

1. calculate and register the exact local file SHA-256;
2. verify the file matches the supported five-column processed-probe schema;
3. lock the expected publication-supported comparison target before inspecting ONTSeq output;
4. pass the non-biological runtime compatibility fixture;
5. register AL_001 as a one-sample validation manifest and execute it once through
   `ontseq marlin-validate`;
6. let the validation harness perform its required repeated classifications and require
   deterministic input/feature identity, decision and score agreement within the frozen runtime
   tolerance.

The repository currently contains no claim that this gate has passed.

## GSE280090 cohort validation

When the available restricted processed files are run, all samples must share the same:

- artifact lock;
- runtime compatibility profile;
- build;
- fixed `0.8` confidence policy;
- preprocessing contract.

The cohort report separates:

- `HIGH_CONFIDENCE`;
- `UNKNOWN`;
- `NO_CALL`;
- technical failures;
- concordant/discordant results where an independently locked comparison label exists;
- feature coverage;
- runtime;
- deterministic rerun failures.

Thresholds or preprocessing must not be altered after observing cohort outcomes to improve
concordance.

## Packaging and data boundary

The ONTSeq package may include the small ONTSeq-owned MARLIN adapter/configuration scripts. It must
not include:

- the trained MARLIN model unless distribution rights and packaging policy are separately
  reviewed and explicitly approved;
- generated runtime fixture/profile outputs from a local model freeze;
- GSE280090 patient-derived processed payloads;
- local validation outputs containing sample-level data;
- private/institutional methylation datasets.

Large/external artifacts remain local controlled resources and are referred to through checksums.

## Release boundary

The branch remains on the 0.8.2 release identity while implementation and validation gates are
open. A coordinated 0.9.0 engineering version bump is permitted only after the planned code,
packaging, full repository verification and review gates pass.

Even after a 0.9.0 engineering release, MARLIN remains Research Use Only until independent
analytical/intended-use and clinical validation requirements are separately satisfied.
