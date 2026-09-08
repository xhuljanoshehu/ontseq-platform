# In-silico dilution series and technical limit of detection

## Status

Research use only. This is a technical characterisation of pipeline behaviour under varying
tumour read fraction. It is **not** an analytical or clinical limit-of-detection study, and no
number it produces may be quoted as an assay sensitivity.

## Where this sits

`docs/BENCHMARKING.md` puts "coverage and tumor/blast-fraction dilution" on rung 3 of the
validation ladder, above CI fixtures and public truth sets and below intended-use validation on
characterised AML specimens. `docs/ROADMAP.md` has carried "benchmark CNV candidates across
coverage and tumor/blast-fraction dilution series" as an open Milestone 1 item. This module is
the machinery for that rung; running it on real characterised material is still the open work.

## The three steps

### 1. Plan — pure arithmetic, no BAM

`plan_dilution_series()` turns two read counts and a policy into the entire series as data:
per-level read budgets, subsample fractions, derived seeds and the exact `samtools view -s`
arguments. It touches no file and runs no tool, so the *design* of a titration can be reviewed,
diffed and unit tested before anyone spends a night of compute on it.

```bash
ontseq dilution-plan \
  --policy configs/benchmark/dilution_series.technical.yaml \
  --series-id AML_TITRATION_001 \
  --tumor-bam tumor.bam --normal-bam normal.bam \
  --tumor-sample-id TUMOR_001 --normal-sample-id NORMAL_001 \
  --genome-build GRCh38 \
  --output results/dilution/plan.json
```

Depth is held constant across levels: a titration is supposed to vary tumour fraction and
nothing else. With `total_read_target: null` the budget is the largest every level can fund from
the two sources, and the plan warns that two series built this way from different inputs are not
comparable. Pin `total_read_target` before comparing series.

Seeds are derived from `policy.seed` alone, and no two subsamples of the same BAM anywhere in the
series share one, so the whole design is reproducible from a single integer.

A source taken whole carries **no** subsample argument. samtools expresses "seed and fraction" as
one float, and there is no representation of *all reads* in that format — `.999999` is not all
reads — so taking a source whole means not subsampling it, and the plan records `null` rather
than an argument that would silently drop reads.

### 2. Mix — and verify the label

```bash
ontseq dilution-mix results/dilution/plan.json \
  --tumor-bam tumor.bam --normal-bam normal.bam \
  --output-dir results/dilution/levels \
  --output results/dilution/series.json
```

Subsampling is random, so a level asked for 5 % tumour reads does not contain exactly 5 %. Both
numbers are recorded per level, and a level that drifts past
`observed_fraction_tolerance` (default 2 percentage points) **fails the series**. A limit of
detection computed from mislabelled levels is worse than no limit of detection.

Each level is counted before the merge, because after it nothing in the mixture says which read
came from which source. Every mixed BAM is fingerprinted.

### 3. Evaluate — a limit that says what it is

Run the pipeline on each mixed level, benchmark each result against the truth set with
`ontseq benchmark`, and label each benchmark case with the level it came from:

```yaml
strata:
  dilution_series_id: AML_TITRATION_001
  tumor_fraction: 0.05
  replicate: 1
```

Those are the strata keys `examples/benchmarks/synthetic_cnv.yaml` already uses, so no second
labelling scheme is needed. `tumor_fraction` is required; `replicate` defaults to 1, which means
two reports for the same level and replicate are refused — a repeated level is not an extra
replicate, and treating it as one inflates the detection rate.

```bash
ontseq lod results/benchmarks/*.json \
  --policy configs/benchmark/lod.technical.yaml \
  --series-id AML_TITRATION_001 \
  --output results/dilution/lod.json
```

## How detection is decided

| Per replicate | Condition |
| --- | --- |
| `DETECTED` | `recall >= policy.minimum_recall` (default `1.0`: every truth event recovered) |
| `NOT_DETECTED` | recall defined and below the threshold |
| `NO_CALL` | recall undefined — the case carried no truth event, so it asked no question |

`NO_CALL` replicates are excluded from the detection rate rather than counted as failures.
Folding them in would depress the rate with evidence that does not exist.

A level passes when it has at least `minimum_replicates` evaluable replicates and its detection
rate reaches `minimum_detection_rate`. With `require_monotonic` (the default) the limit is the
lowest passing level that has no failing level *above* it, so one lucky low level cannot be
reported as the limit.

## Bracketing: the honest part

If the lowest fraction tested still detects, the series has not found where detection stops. It
has shown only that the stopping point is somewhere below the lowest level. `LodReport.bracketed`
is `true` only when a failing level was actually observed below the limit; otherwise the report
warns in plain words that the number is a bound from above, not a limit. Extend the series
downwards before quoting it.

When no level passes at all, `detection_limit_fraction` is `null` and the report says explicitly
that this did not establish that detection is impossible either.

## What an in-silico series cannot tell you

- It reproduces read-fraction effects and **nothing** about library preparation, input mass, or
  capture behaviour at low tumour content. A wet-lab dilution and a read mixture are different
  experiments.
- Levels share reads with one another and with the undiluted control, so replicates within a
  series are not statistically independent specimens and the detection rate is not a confidence
  statement about future samples.
- The estimate can be no finer than the fractions tested: the true limit lies between the lowest
  passing level and the highest failing one.
- One replicate per level exercises the machinery and characterises nothing. The technical
  default of `minimum_replicates: 1` exists so the lane runs on a demonstration series and must
  be raised before any result is quoted.

## Files

| Path | What it is |
| --- | --- |
| `configs/benchmark/dilution_series.technical.yaml` | Fractions, seed, replicates, drift tolerance |
| `configs/benchmark/lod.technical.yaml` | Detection criterion |
| `schemas/dilution-series-plan.schema.json` | The planned series |
| `schemas/dilution-series-report.schema.json` | What was actually mixed |
| `schemas/lod-report.schema.json` | Per-level outcomes and the limit statement |
