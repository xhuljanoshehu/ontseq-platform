# CNV truth, comparison and benchmarking

## Status

Research use only. The comparison harness is tested; no CNV method has been validated,
selected or promoted. Every number this subsystem produces is an engineering measurement
against a declared truth set, never a clinical performance claim.

## Why the existing benchmark layer was not enough for CNV

`src/ontseq_platform/benchmark.py` matches normalized events one-to-one under a
reciprocal-overlap threshold. That contract is right for structural variants, where a
breakend pair is the unit of truth. It misreports copy number for three reasons:

1. **Segmentation is not part of the biological claim.** One truth deletion emitted by a
   caller as three adjacent deleted segments scores 1 TP and 2 FP under one-to-one
   matching. The caller was right; the segmentation differed. Bin size and segmentation
   algorithm are implementation details, not findings.
2. **Event counting discards event size.** A whole-chromosome gain and a 200 kb
   duplication contribute equally to precision and recall.
3. **Unmatched is treated as wrong.** `false_positive = len(query) - TP` assumes every
   base was assessable. Centromeres, assembly gaps, regions below the coverage floor and
   regions outside an adaptive-sampling target design are not, and scoring them converts
   a known blind spot into a fabricated error rate.

The SV path is unchanged. `ontseq benchmark` and its fixtures still work exactly as
before; the CNV path is a separate subsystem under `src/ontseq_platform/cnv/`.

## Core idea: compare per-base states, not event lists

Copy number is a *segmentation of the genome into states*. The comparator therefore
builds an exact breakpoint partition of the reference from all truth boundaries, all call
boundaries and all mask boundaries, then assigns each elementary interval a truth state
and a called state and accumulates a base-pair-weighted confusion matrix.

Base-level agreement is invariant to how either side segmented. Event-level detection is
still reported, but derived from base-level concordance with a many-to-many rule, so
fragmentation cannot be penalised.

### The state vocabulary

`homozygous_loss`, `loss`, `neutral`, `gain`, `high_amplification`, `copy_neutral_loh`,
`no_call`.

Two rules matter:

- `no_call` is never concordant with anything, including itself. Those bases leave the
  evaluable genome rather than counting as agreement.
- `copy_neutral_loh` is only ever concordant with itself, even in directional mode.
  Collapsing it into `neutral` would let a dosage-only caller take credit for detecting
  an event it is structurally incapable of seeing.

### Concordance modes

- `directional` (default): loss-family agrees with loss-family, gain with gain. This is
  the clinically relevant question for karyotyping and it is robust at low coverage.
- `strict`: exact state equality. For method development where the exact state is the
  question.

## The evaluable genome

Every metric is computed over an explicitly constructed evaluable region, and everything
removed is accounted for by reason:

| Reason | Meaning |
| --- | --- |
| `outside_analysis_scope` | Outside the panel/ROI or declared evaluation scope |
| `assembly_gap`, `centromere`, `low_mappability`, `blacklist` | Structurally unusable |
| `below_coverage_floor` | Observed depth under the configured floor |
| `caller_no_call` | The method declined to call |
| `truth_not_informative` | The truth source could not assess this region |
| `contig_not_in_reference` | Contig absent from the reference lock |

This is what makes the three states the project cares about distinguishable:

- **negative** - observable, nothing found;
- **no-call** - not observable, no statement possible;
- **failure** - module did not run, carried by `ModuleRunStatus`.

Bases excluded by several tracks are attributed to the first track that removed them, so
per-reason counts sum exactly to the total removed.

### Comparing several methods fairly

`compare_methods()` removes the **union** of every method's no-call regions from the
shared mask before scoring any of them. Without this, a cautious method is rewarded for
its own blind spots, because declining to call shrinks only its own denominator.

## Open world versus closed world

Every truth set and every call set declares a `background_state`, which says what silence
means:

- `neutral` (**closed world**) - regions not listed are asserted to be unaltered. A SNP
  array within its probe map, a karyotype at band resolution, a genome-wide segmenter.
- `no_call` (**open world**) - regions not listed carry no assertion. FISH, a targeted
  panel, an alteration-only caller.

Getting this wrong is silent and severe. Treating an open-world truth as closed-world
manufactures a false positive for every genuine finding outside its scope; the reverse
hides real false positives. There is deliberately no default.

Closed-world truth sets must also declare `resolution_bp`: the smallest event the source
could detect. Below it the truth is *silent*, not negative, and calls smaller than the
truth's resolution are flagged in the report rather than counted as errors.

## Truth sources and what they can support

| Source | Breakpoint resolution | Background | CN-LOH | Notes |
| --- | --- | --- | --- | --- |
| Simulation | exact | closed | by construction | Harness behaviour only |
| ISCN karyotype | band width, often 5-20 Mb | closed, ~10 Mb | no | Breakpoint metric suppressed |
| FISH | probe locus | open | no | Probe-scoped only |
| SNP array | probe spacing | closed within probe map | yes | Sharp breakpoints |
| CGH array | probe spacing | closed within probe map | no | Dosage only |
| Short-read WGS/panel | base pair | closed / ROI | with allele data | Orthogonal pipeline |

### Breakpoint accuracy is suppressed when the truth cannot support it

Boundary uncertainty travels on each truth segment. When a truth event's uncertainty
exceeds `maximum_truth_boundary_uncertainty_bp`, the breakpoint metric is reported as
`null` with reason `truth_boundary_resolution_insufficient`, never as a number. Measuring
a caller's breakpoint error against `5q13` would be measuring the width of a Giemsa band.

A **call set's** declared uncertainty is deliberately ignored: a method must not be able
to excuse its own breakpoint error.

### ISCN karyotype conversion

`ontseq cnv-karyotype-truth` converts a karyotype into band-resolved intervals using a
versioned, checksummed cytoband table.

Supported: `+N`, `-N`, `del(N)(band)`, `del(N)(bandband)`, `dup(...)`, `i(N)(q10)`.
Recognised as balanced and contributing no segments: `t(...)`, `inv(...)`, `ins(...)`.

Everything else - derivative chromosomes, markers, `add`, anything containing `?` - is
recorded as an **unsupported construct with a reason** and the command exits non-zero. A
truth set that quietly lost half its findings is worse than no truth set: it turns real
events into apparent false positives.

The repository ships no cytoband table. Band definitions are build-specific annotation
that belongs with the reference bundle on approved storage; see `docs/DATA_SECURITY.md`.
`examples/references/synthetic.cytobands.txt` is a coarse invented fixture for tests.

## Metrics and why each one is reported

| Metric | Unit | Interval | Rationale |
| --- | --- | --- | --- |
| Base-level concordance | bp fraction | no | Segmentation-independent overall agreement |
| Per-state confusion matrix | bp | no | Shows *how* a method is wrong, not only that it is |
| Per-state bp recall / precision | fraction | no | Directional performance by state |
| Event detection rate | events | **yes** | Sensitivity on the clinically meaningful unit |
| Event confirmation rate | events | **yes** | Specificity proxy over assessable calls |
| Detection by size class | events | yes | Resolution is a strong function of event length |
| Detection by state | events | yes | Losses and gains do not behave alike |
| Copy-number MAE / RMSE | copies | no | Quantitative accuracy, bp-weighted |
| Breakpoint deltas | bp | no | Only where truth resolution supports it |
| Partition accounting | bp | no | How much genome each number applied to |

### Confidence intervals are event-level only

Wilson score intervals are reported on **event counts**, never on base-pair proportions.
Base pairs inside one segment are not independent observations; a binomial interval over
900 million bases would be absurdly narrow and would assert a precision that does not
exist. Base-level figures are descriptive point estimates and are labelled as such.

The Wilson interval is used rather than the normal approximation because benchmark
proportions routinely sit at 0 or 1 with small denominators, where the normal interval
collapses to zero width.

### Undefined stays undefined

A proportion with a zero denominator is `null`, never `0.0` or `1.0`. An unevaluated
stratum must not be able to look like a perfect or a failed one.

## Blast fraction, dilution series and limit of detection

A specimen with blast fraction `f` presents `CN_mix = f * CN_t + (1 - f) * CN_n`. Every
alteration compresses towards the baseline as `f` falls, which is the fundamental reason
low-blast specimens are hard.

`simulate_dilution_series()` generates reproducible series over blast fraction and
coverage. Counts are drawn from a negative binomial (gamma-Poisson), because sequencing
depth is consistently overdispersed relative to Poisson.

`estimate_limit_of_detection()` reports **two** estimates because they fail differently:

- **Empirical**: the lowest tested level whose *lower* confidence bound still meets the
  target. Conservative, cannot extrapolate, limited to tested levels.
- **Model-based**: from a logistic fit. Interpolates, but is withheld (`null`) when the
  fit did not converge, fewer than two levels were supplied, or the design is perfectly
  separated - because those designs cannot support a finite detection limit.

## The baseline caller is a control, not a product

`ontseq-baseline-readdepth` is a plain binned read-depth segmenter: median normalisation,
log2 ratios, recursive binary segmentation on a maximum standardised-difference statistic
with a robust noise scale from median absolute successive differences.

It exists for three reasons: it closes the loop so the harness is tested rather than
merely written; it is a null model any candidate must beat before its complexity is
justified; and it is small enough to debug when the harness reports something surprising.

It has no GC correction, no mappability correction, no allele-fraction information, no
ploidy search and no subclonal deconvolution. It assumes bins are comparable, which is
true for uniform low-coverage whole-genome data and **false for on-target capture or
adaptive-sampling enrichment**. It is never reportable.

## Assay-agnostic by construction: `data_basis`

Adaptive sampling produces two very different read populations in one run:

- **Rejected (off-target) reads** still occupy the flow cell and form a near-uniform
  low-coverage whole-genome background. This is the population that depth-based CNV
  methods such as ichorCNA, QDNAseq and Spectre actually assume, and it is the most
  promising basis for genome-wide CNV from an adaptive-sampling run.
- **On-target reads** are deeply but non-uniformly enriched. Their depth is dominated by
  enrichment efficiency rather than copy number, which violates the core assumption of
  every method above.

`CnvDataBasis` makes this an explicit, mandatory stratification key
(`adaptive_sampling_off_target`, `adaptive_sampling_on_target`,
`adaptive_sampling_combined`, `low_coverage_wgs`, `whole_genome`, `simulated`). Pooling
them into one benchmark stratum compares incomparable things.

**Unvalidated assumption:** that the off-target fraction of the local adaptive-sampling
assay is uniform enough for depth-based CNV at the achieved yield. Nothing in this
repository tests that yet; see "Next steps".

## Connection to target-coverage QC

PR #7 adds a Mosdepth target-coverage adapter. `coverage_floor_track()` consumes
per-interval depth and converts intervals below a floor into
`below_coverage_floor` exclusions, so a region the assay could not observe becomes a
no-call rather than a copy-number-neutral finding.

The floor is an engineering parameter for benchmarking. It is **not** a validated
adequacy threshold and must not be reused as a clinical no-call limit.

## Adding a method

Column mappings, not parsers. `ColumnMapping` describes a tool's segment table by header
**name**; `parse_segment_table()` raises when a required column is absent.

Positional parsing is the classic silent scientific error: a tool adds a column in a
patch release, every downstream number shifts by one field, and the pipeline keeps
producing plausible results. A loud failure is strictly preferable.

Shipped mappings: generic IGV `SEG`, and ichorCNA (`Corrected_Copy_Number` /
`Corrected_Call`, derived from the ichorCNA output documentation).

**Not shipped:** Spectre and QDNAseq/ACE mappings. Their exact column layouts were not
confirmed against upstream sources while this module was written, and shipping a guessed
mapping as a verified default is exactly the failure mode the design avoids. Supply a
mapping explicitly and record it in provenance.

These adapters only *parse*. Execution belongs behind the repository's existing adapter
boundary with version pinning and argument-vector invocation, as done for Sniffles2 and
Mosdepth.

## Running it

```bash
# Fully synthetic end-to-end benchmark: simulate, call, evaluate, aggregate.
ontseq cnv-demo-benchmark --output-dir results/cnv-demo

# Score one locked case.
ontseq cnv-evaluate examples/benchmarks/cnv_case_synthetic.yaml \
  --output results/cnv/evaluation.json

# Pool evaluations of one method.
ontseq cnv-aggregate results/cnv-demo/evaluations/*.json \
  --aggregate-id MY_RUN --output results/cnv/aggregate.json

# Convert a karyotype into band-resolved truth.
ontseq cnv-karyotype-truth \
  --karyotype "46,XY,del(5)(q13q33),-7,+8" \
  --cytobands examples/references/synthetic.cytobands.txt \
  --cytoband-resource-id SYNTHETIC_CYTOBANDS_V1 \
  --genome-build GRCh38 \
  --truth-id KARYO_001 --sample-id SYNTHETIC_AML_001 \
  --source-version synthetic-v1 \
  --output results/cnv/karyotype-truth.json
```

## Pre-registration requirement

Thresholds - `detection_overlap_fraction`, `minimum_assessable_fraction`,
`concordance_mode`, `copy_number_tolerance`, the coverage floor and the exclusion tracks -
must be locked **before** comparative results are inspected. Every one of them changes
reported performance, which is why the full option set is echoed into each report.

## What has and has not been shown

**Tested** (see `tests/test_cnv_*.py`):

- fragmented calls score identically to a single equivalent call;
- events in unobservable regions are `NOT_ASSESSABLE`, not false negatives;
- open-world truth does not manufacture false positives, closed-world truth does count
  them;
- breakpoint metrics are withheld under band-resolution truth and reported under exact
  truth, including when a caller overshoots;
- CN-LOH is not satisfied by a neutral call;
- undefined proportions stay `null`;
- the simulator is deterministic and its dilution series reproducible;
- the whole loop runs end to end with no external tool.

**Not shown:**

- that any method performs adequately on real ONT data;
- that the off-target fraction of the local adaptive-sampling assay supports depth-based
  CNV;
- that the synthetic noise model resembles real ONT bin-count behaviour;
- anything about GC bias, mappability, alignment error or subclonal structure;
- any clinical or analytical validity whatsoever.
