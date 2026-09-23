# MARLIN classification in ONTSeq

## Status

**Research Use Only. No clinical release.**

This document describes the engineering integration of the published MARLIN v1 acute-leukemia
methylation classifier into ONTSeq. It does not claim analytical or clinical validity and it does
not claim that ONTSeq reproduces the complete raw-signal workflow used in the publication.

Two independently versioned paths are available. The automatic desktop/pipeline path is
`marlin-native-research-v1`, schema `1.0.0`, with the mandatory validation status
`UNVALIDATED_RESEARCH`. The older standalone R commands retain their stricter GRCh37-only
contracts and validation gates. Installing the native adapter creates no `MarlinBridgeLock`
and does not satisfy or bypass those legacy command gates.

Legacy standalone command scope:

- strict import of the published five-column probe-level methylation representation;
- explicit GRCh37/hg19-only v1 build contract;
- deterministic construction of the 357,340-value MARLIN v1 feature vector;
- checksum-locked model, feature, class-annotation and probe resources;
- locked R/Keras/TensorFlow inference boundary returning exactly 42 softmax scores;
- live runtime probing of R, R-Keras, Python TensorFlow, Python and CPU/GPU backend identity;
- deterministic non-biological runtime-fixture generation and one-time reference-profile freeze;
- explicit dual-runtime numerical compatibility qualification using separate reference and
  candidate execution locks over one identical immutable MARLIN artifact set;
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

The runtime-freeze and dual-runtime comparison implementations are present, but this repository
does **not** contain a generated compatibility profile, dual-runtime PASS report, or candidate
self-baseline from the real Zenodo MARLIN model. Those remain local controlled validation
artifacts and must be established before biological validation begins.

## Automatic native research path

Selecting methylation also requests the MARLIN stage in the normal analysis. Desktop readiness
checks configuration, expected resource presence and a declared independently qualified runtime
profile without executing that runtime. Its message explicitly states that installed bytes and
execution are not yet verified. The analysis stage verifies all bytes before runtime preflight. A missing
installation becomes `NOT_RUN` with a reason; a parser, checksum, modkit or worker failure becomes
`FAILED`. MARLIN failure does not turn available CNV/SV results into failures. Deselecting
methylation excludes MARLIN, including leftover artifacts, from the current report.

The adapter accepts directly aligned GRCh37/hg19 and GRCh38/hg38 BAMs with verified modified-base
information. It chooses the original probe map bound to that build by a fixed SHA-256. No
realignment or liftover is performed. The native adapter checks sample, reference and intake
identity before executing a separate model-probe pileup. The selected adjacent BAI/CSI must
match the intake's required SHA-256 fingerprint; competing distinct BAM-adjacent indexes are
refused rather than left to modkit's index search order. Its exact bytes and resolved selection
are checked again immediately before and after pileup and before any final outcome, including
`NO_CALL`. The reference's adjacent `.fai` is likewise mandatory, fingerprinted and rechecked;
compressed references are refused because their additional `.gzi` dependency is not qualified.
Both index fingerprints are mandatory in concluded native reports and resume signatures:


```text
modkit 0.6.4 pileup --modified-bases 5mC 5hmC --combine-mods --cpg
                  --ref <matching FASTA> --include-bed <model probe positions>
                  --filter-threshold 0.8
```

The two modification probabilities are combined by modkit before count aggregation. For each
model probe, ONTSeq calculates `sum(N_mod) / sum(N_valid)` over its mapped positions and strands.
It preserves per-strand calls instead of requesting `--combine-strands`; duplicate site/strand
records and a mixture of combined and separate strand rows are rejected. Count identities,
finite fractions, geometry and exact feature order are checked. Measured beta zero encodes as
`-1`; an unobserved feature encodes as `0`. The model input has exactly 357,340 float32 values.

The same modkit guard used by the regional methylation lane remains active: stock 0.6.4 is
refused for independent cytosine MM groups. Only the exact locally qualified PR #709 executable
SHA-256 receives the exception. A version banner or an adjacent receipt alone cannot grant it.
A zero exit status with failed-record diagnostics is a failure. Regional 5mC measurements keep
their own modification and coverage policy; MARLIN's combined counts do not replace them.

When no model feature is observed, the result is `NO_CALL`, with decision `UNKNOWN`, no scores,
and no model execution. Otherwise the separate Python 3.10 / TensorFlow CPU 2.13.1 / Keras 2.13.1
worker checks the original model SHA-256, exact tensor hash and input/output shapes. It loads
with `compile=False` and evaluates with `training=False`. Linux seccomp denies network syscalls
before TensorFlow imports; execution preflight exercises that kernel rule only after full
runtime byte verification. GPU visibility is disabled,
model intra-op threads are capped at eight, and inter-op threads are fixed at one.

Exactly 42 finite model scores are grouped using the original XLSX annotations into current
classes, families and lineages. `model_score_threshold` is fixed at `0.8` and
`model_score_threshold_met` records whether the leading grouped class reaches it. This is a
raw model-score comparison, not specimen confidence. Native predictions always have decision
`UNKNOWN` and `assay_assessability: NOT_ESTABLISHED`, regardless of score or observed CpG count:
no independently validated ONTSeq specimen-assessability policy exists. Additional observed
CpGs alone cannot establish analytical validity, including coverage of all 357,340 model features.
The threshold result is absent for `NO_CALL`, `NOT_RUN` and `FAILED`. Raw scores, top class and
observed feature count/fraction remain available for completed research inference. Every outcome
retains `UNVALIDATED_RESEARCH`; scores are model rankings, not clinical diagnosis probabilities.

The coverage test grid (1, 10,720 and 357,340 observed model features) guards against promoting a
high score to specimen confidence. The intermediate count is a test point, not a clinical cutoff.
The [MARLIN publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC12513838/) reports benchmarking
under its study conditions; those results do not qualify an ONTSeq specimen-assessability policy.

### Installation and provenance

The resource kit supplies `<resource-root>/marlin/installation.json`; a separately qualified
installation can be selected with `--marlin-installation` for `run`, `analyze` or `serve`.
The currently qualified durable Linux layout uses `marlin-native-1.0.0/env` under the local
ONTSeq runtime root. TensorFlow is not added to
the ONTSeq control-plane environment. Local installation retains the model's CC BY 4.0
attribution, MARLIN's MIT license text and the runtime dependency notices. Archive extraction
checks member paths, types, links and size before restoring the runtime; relocation completes
before the final file inventory is generated. The installation contract records:

- original model, ordered feature list, original class workbook and both build-keyed probe maps,
  each with an exact checksum and source URI;
- selected Python executable, exact Python patch version, TensorFlow/Keras versions and the
  network-confinement contract;
- original runtime archive identity and a complete manifest of installed runtime files,
  regenerated after relocation to the final absolute prefix.

The model identity is `6210a674a0e7690b7b6c03184f5732a65037cf00ff4383e1892f26e90fdcb217`.
Runtime trust is independent of that installation JSON. The approved archive is fixed at
SHA-256 `b812927398645d3abaa3a56622bbf1b3279b70ff25067f42040d450083c822c6`. Two separately
verified installed-inventory identities are supported: `marlin-native-durable-linux-cpu-v1`
(manifest SHA-256 `1f7d1d996b06b851ab4116162bfff86ef9e86ec13ed1cf9e86de7ebfcb09250b`) and
`marlin-native-source-linux-cpu-v1`
(`b97f0b2468923415584b756306ed95ec5941d3b26ecb423095a2ac9135b7a461`). Those manifests bind
absolute runtime prefixes as well as bytes. Matching version banners and regenerating a
self-supplied manifest cannot qualify changed code or a new prefix. A new relocation requires
separate archive restoration, byte/inference qualification and reviewed registration of its
inventory digest; generic portability to arbitrary prefixes is not claimed.

The native report (`native-marlin.json`) binds run/sample/build, feature summary, tools,
parameters, installation signature, and exact BAM/index/reference/FAI/pileup/tensor/worker fingerprints.
Every signature request verifies actual installed bytes, including resume; added or missing
runtime files invalidate the inventory. The exact relative directory-symlink map from the approved
archive is also checked, preventing links from bypassing file enumeration. Source hashes identify the adapter, contracts and worker.
An existing nonempty stage directory is not overwritten. Failed outputs remain diagnostic
artifacts and contribute no current prediction. Archived reports remain unchanged.

HTML with and without JavaScript, JSON and XLSX consume the current typed MARLIN outcome.
An `UNKNOWN` leading class is shown as model ranking, not a confirmed classification. Missing
values remain missing. Run/sample/build and artifact hashes must agree before score import.

### Technical acceptance and remaining validation

The 2026-09-23 local synthetic acceptance generated a real modified-base BAM at four original
hg38 model probe positions, with forward/reverse calls and unequal 10/90 depths. The exact
qualified modkit executable produced the independently specified betas `0.49`, `0.50`, `0`, `1`.
The adapter's ordered tensor matched the independent oracle, and the original model returned
42 scores identical to a separate direct inference call. Its maximum difference from TensorFlow
`predict` was `3.3527612686157227e-08`, below the existing `1e-7` engineering tolerance. A separate
sample with tags but no model-probe coverage produced `NO_CALL` without inference.

These results establish synthetic data transfer and original-model execution on the inspected
runtime. They do not establish biological concordance, assay coverage suitability, sensitivity,
clinical validity, equivalence of historical R bindings, or a same-specimen validated bridge.
Independent intended-use validation and expert review remain required. The legacy biological
bridge and R runtime gates below retain their original meanings.

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

## Legacy validation evidence ladder

ONTSeq deliberately separates engineering evidence from analytical and clinical evidence:

```text
synthetic parser / contract / feature tests
        <
real R/Keras/TensorFlow import and backend smoke
        <
single-runtime frozen reproducibility baseline
        <
dual-runtime numerical compatibility
        <
candidate same-lock frozen reproducibility baseline
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
- A single-runtime frozen profile proves reproducibility against its own frozen numerical oracle;
  it does not by itself prove equivalence between two different runtime environments.
- A dual-runtime PASS proves only that the live candidate reproduced the frozen reference model
  output on the fixed non-biological fixture within the predeclared engineering tolerances.
- The candidate same-lock freeze after dual-runtime PASS is a reproducibility prerequisite for
  `marlin-validate`, not a new biological-validation claim.

## Legacy standalone architecture

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

The legacy validated native-bridge path is separate:

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

This legacy validated-bridge command is not enabled by the native research installation. It requires a
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

For dual-runtime qualification, reference and candidate environments receive separate execution
locks. ONTSeq derives a canonical `MarlinArtifactSetIdentity` from each lock and requires the same
artifact-set SHA-256 before candidate inference. The canonical identity includes immutable
model/code/feature/annotation/probe/build/contract fields but excludes runtime versions, backend,
lock IDs, timestamps, local paths and source-URI formatting.

## Legacy R engineering runtime and live probe

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
model must be established independently by the frozen non-biological runtime fixture and the
dual-runtime comparison before any biological validation is permitted.

The smoke workflow stores only a small non-biological runtime-identity JSON artifact. It does not
contain the trained MARLIN model, GSE280090 payloads or patient-derived data.

## Legacy runtime compatibility gates

External biological validation requires more than a successful runtime import or a stored profile.
ONTSeq therefore uses a reference reproducibility gate, a cross-runtime numerical gate, and a
candidate same-lock reproducibility gate before AL_001/GSE280090.

### Single-runtime reference baseline

ONTSeq executes a fixed, non-biological 357,340-value `-1/0/+1` feature fixture through the
reference locked model/runtime. The canonical generator is `sha256-index-mod3-v1`; its feature
vector and text representation have regression-locked SHA-256 identities documented in
`docs/MARLIN_RUNTIME_FREEZE.md`.

The engineering tolerances are fixed before biological validation:

```text
absolute tolerance per raw model score: 1e-7
softmax-sum tolerance:                  1e-5
```

`ontseq marlin-freeze-runtime` probes the same live `Rscript` environment used for inference,
requires exact agreement with the reference execution lock, runs the fixed vector once, and writes
the local fixture plus 42-score profile. With `--output-runtime-identity`, it also persists the
`MarlinRuntimeIdentity` created from that same live probe and freeze timestamp. Existing output
paths are refused rather than overwritten.

### Dual-runtime numerical compatibility

A separate `ontseq marlin-compare-runtimes` command qualifies the live candidate environment
against the frozen reference oracle. Before candidate inference it requires:

- reference and candidate execution locks with identical canonical artifact-set SHA-256;
- reference runtime identity consistent with the reference execution lock;
- reference profile bound to the reference execution lock;
- frozen fixture digest matching the reference profile;
- model bytes matching the shared artifact-set model SHA-256;
- a live candidate runtime probe consistent with the candidate execution lock.

The candidate then executes the same fixed fixture once. A machine-readable PASS requires all 42
candidate scores to differ from the frozen reference by `<= 1e-7`, the same top raw model-unit
index, and both softmax sums within `1e-5`. A structurally valid numerical incompatibility produces
`verdict=FAIL`; malformed/inconsistent evidence fails the command. Neither state is a biological
`UNKNOWN` or `NO_CALL`.

### Candidate same-lock validation baseline

A dual-runtime PASS does not change the existing `marlin-validate` identity rule. Biological
validation still requires a `MarlinRuntimeCompatibilityProfile` whose
`reference_runtime_lock_id` equals the **candidate** execution lock ID. Therefore, after the
cross-runtime report is PASS, run `marlin-freeze-runtime` once under the candidate execution lock
to create the candidate fixture/profile self-baseline. This second freeze uses the same
deterministic generator and fixed tolerances; it is not a biological qualification and must not be
used to retune any policy.

`marlin-validate` then uses the candidate execution lock, candidate profile and candidate fixture.
A reference profile supplied with the candidate lock remains a hard failure by design.

The tolerances are frozen before biological validation. They must not be tuned from AL_001 or
GSE280090 outcomes. A direct `marlin-classify` invocation is not a substitute for external
validation.

## Legacy standalone classification semantics

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

## Legacy validated modkit bridge

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
ontseq marlin-compare-runtimes
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

Freeze the **reference** runtime before biological validation and persist its live identity:

```bash
ontseq marlin-freeze-runtime \
  --artifact-lock results/marlin-validation/reference-artifact-lock.json \
  --model /path/to/marlin_v1.model.hdf5 \
  --inference-script scripts/marlin_infer_locked.R \
  --runtime-probe-script scripts/marlin_runtime_probe.R \
  --profile-id MARLIN_V1_REFERENCE_FROZEN \
  --rscript /reference-runtime/bin/Rscript \
  --output-fixture results/marlin-validation/reference-runtime-fixture.txt \
  --output-profile results/marlin-validation/reference-runtime-profile.json \
  --output-runtime-identity results/marlin-validation/reference-runtime-identity.json
```

Then create the candidate execution lock over the **same artifact bytes** and qualify its live
runtime against the frozen reference:

```bash
ontseq marlin-compare-runtimes \
  --reference-artifact-lock results/marlin-validation/reference-artifact-lock.json \
  --candidate-artifact-lock results/marlin-validation/candidate-artifact-lock.json \
  --reference-runtime-identity results/marlin-validation/reference-runtime-identity.json \
  --reference-profile results/marlin-validation/reference-runtime-profile.json \
  --runtime-fixture results/marlin-validation/reference-runtime-fixture.txt \
  --model /path/to/marlin_v1.model.hdf5 \
  --inference-script scripts/marlin_infer_locked.R \
  --candidate-runtime-probe-script scripts/marlin_runtime_probe.R \
  --candidate-rscript /candidate-runtime/bin/Rscript \
  --comparison-id MARLIN_V1_REFERENCE_VS_ONTSEQ_CPU \
  --output results/marlin-validation/runtime-comparison.json
```

After and only after `runtime-comparison.json` reports `PASS`, freeze the candidate runtime under
its own execution lock to create the same-lock validation baseline:

```bash
ontseq marlin-freeze-runtime \
  --artifact-lock results/marlin-validation/candidate-artifact-lock.json \
  --model /path/to/marlin_v1.model.hdf5 \
  --inference-script scripts/marlin_infer_locked.R \
  --runtime-probe-script scripts/marlin_runtime_probe.R \
  --profile-id MARLIN_V1_CANDIDATE_FROZEN \
  --rscript /candidate-runtime/bin/Rscript \
  --output-fixture results/marlin-validation/candidate-runtime-fixture.txt \
  --output-profile results/marlin-validation/candidate-runtime-profile.json
```

The generated fixture/profile/runtime-identity/comparison files are local engineering evidence and
are not repository assets. AL_001/GSE validation begins only after the controlled real-model
dual-runtime comparison is PASS and the candidate same-lock profile has been frozen.

### External validation

The validation runner requires a checksummed manifest, the candidate execution lock/profile and the
candidate frozen runtime fixture:

```bash
ontseq marlin-validate \
  --manifest results/marlin-validation/gse280090-v1/manifest.json \
  --artifact-lock results/marlin-validation/candidate-artifact-lock.json \
  --runtime-profile results/marlin-validation/candidate-runtime-profile.json \
  --runtime-fixture results/marlin-validation/candidate-runtime-fixture.txt \
  --model /path/to/marlin_v1.model.hdf5 \
  --feature-rdata /path/to/marlin_v1.features.RData \
  --feature-list /path/to/marlin_v1.features.txt \
  --class-annotations /path/to/marlin_v1.class_annotations.xlsx \
  --probe-bed /path/to/marlin_v1.probes_hg19.bed.gz \
  --inference-script /path/to/marlin_infer_locked.R \
  --output results/marlin-validation/gse280090-v1/report.json
```

Paths above are examples. Artifact identities are accepted only when they match their locks. The
controlled validation procedure retains the dual-runtime PASS report as preceding engineering
evidence while `marlin-validate` enforces the candidate's own same-lock profile.

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
4. require the dual-runtime comparison report to be `PASS`;
5. freeze and verify the candidate same-lock runtime fixture/profile;
6. register AL_001 as a one-sample validation manifest and execute it once through
   `ontseq marlin-validate`;
7. let the validation harness perform its required repeated classifications and require
   deterministic input/feature identity, decision and score agreement within the frozen candidate
   runtime tolerance.

The repository currently contains no claim that this gate has passed.

## GSE280090 cohort validation

When the available restricted processed files are run, all samples must share the same:

- candidate artifact lock;
- candidate runtime compatibility profile;
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
- generated runtime fixture/profile/runtime-identity/comparison outputs from local qualification;
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
