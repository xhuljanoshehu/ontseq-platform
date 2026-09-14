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
- grouped current-class, methylation-family and lineage summaries;
- published `0.8` high-confidence threshold with explicit `UNKNOWN` below threshold;
- frozen non-biological runtime-fixture verification before external cohort validation;
- manifest-driven external processed-CpG validation with deterministic repeat runs;
- a native modkit-to-MARLIN bridge contract that remains unusable without a separately created
  same-specimen bridge-evidence lock.

The external GSE280090 / `GSM8587229_AL_001.txt.gz` validation gate has **not** been completed in
this repository state. Public/patient-derived validation payloads are not committed to Git.

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
modBAM -> modkit 0.6.4 bedMethyl
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
- declared R, Keras and TensorFlow versions;
- execution backend;
- creation timestamp.

Inference performs no automatic artifact download. Hash mismatch is a hard failure.

## Runtime compatibility gate

External biological validation requires more than a stored profile. Before any validation
sample is run, ONTSeq re-executes a fixed, non-biological 357,340-value `-1/0/+1` feature fixture
through the current locked model/runtime.

The candidate runtime must match the frozen `MarlinRuntimeCompatibilityProfile` for:

- artifact/runtime lock identity;
- execution backend;
- feature-vector digest;
- top model-unit identity;
- all 42 scores within the predeclared absolute tolerance;
- softmax-sum tolerance.

The tolerance is frozen before biological validation. It must not be tuned from GSE280090
outcomes.

`ontseq marlin-validate` therefore requires an explicit `--runtime-fixture` path.

## Classification semantics

The 42 raw model-unit scores are retained. ONTSeq then uses the locked class annotation resource
to aggregate scores into:

- current methylation class;
- methylation family;
- lineage.

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

The upstream MARLIN workflow selects the probe resource positions, reads modkit `N_valid` and
`N_mod`, and calculates for each probe:

```text
beta = sum(N_mod) / sum(N_valid)
```

ONTSeq mirrors that probe-level rule using its pinned modkit 0.6.4 bedMethyl parser. A probe with
no valid calls remains `NA` and becomes `0` only at MARLIN feature construction.

The bridge is gated by `MarlinBridgeLock`, which records:

- GRCh37/hg19 identity;
- adapter identity;
- modkit version `0.6.4`;
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
ontseq marlin-features
ontseq marlin-classify
ontseq marlin-validate
```

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
5. run the sample twice;
6. require deterministic feature identity, decision and score agreement within the frozen runtime
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
