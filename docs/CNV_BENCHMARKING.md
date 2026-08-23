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

### The partition reconciles, and that is enforced

An accounting that does not add up cannot be audited, so both identities below are
asserted in the core, re-validated by the contract, and covered by a test:

```
reference_bases == mask_bases + excluded_bases
mask_bases      == evaluable_bases + truth_silent_bases + query_no_call_bases
```

`mask_bases` is what the observability mask allowed through; `evaluable_bases` is what
was actually scored, after removing bases where either side was silent. The two silence
counters are attributed exclusively (truth first) and are counted **only inside the
mask** — a caller declining outside the mask is already accounted for as an exclusion,
and counting it twice would make the table unreadable.

This is the table that turns "recall 0.8" into a statement with a denominator.

### Comparing several methods fairly

`compare_methods()` removes the **union** of every method's no-call regions from the
shared mask before scoring any of them. Without this, a cautious method is rewarded for
its own blind spots, because declining to call shrinks only its own denominator.

For choosing between two methods, `paired_detection_comparison()` goes further and pairs
outcomes **per truth event**, keeping only events assessable under both. Comparing two
independent detection rates ignores that the same events drive both numbers and therefore
overstates the uncertainty of the *difference*; McNemar's exact test on the discordant
pairs is the appropriate paired test. It returns `null` when there are no discordant
pairs, because perfect agreement is an absence of evidence, not evidence of equivalence —
and with the cohort sizes realistic here, a non-significant result is almost always a
power problem rather than a demonstrated tie.

```bash
ontseq cnv-compare-methods \
  --method-a results/spectre/*.json \
  --method-b results/ichorcna/*.json \
  --output results/cnv/spectre-vs-ichorcna.json
```

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

### Resolution: what a source cannot deny

Closed-world truth sets must also declare `resolution_bp`: the smallest event the source
could detect. Below it the truth is *silent*, not negative — and that is enforced in the
scoring, not merely mentioned in a warning. Bases under a call finer than the truth can
resolve leave the evaluable genome and are accounted as `truth_resolution_silent_bases`,
the fourth term of the partition identity:

```
mask_bases == evaluable_bases
            + truth_silent_bases
            + truth_resolution_silent_bases
            + query_no_call_bases
```

The affected calls are `NOT_ASSESSABLE`: neither confirmed nor false positives.

**The rule is asymmetric, and that is the substance of it.** Resolution limits what a
source can **deny**, never what it can **affirm**. A karyotype read at 10 Mb bands cannot
rule out a 200 kb duplication, so it cannot make that call wrong. Where the same karyotype
explicitly reports a deletion it has made a positive claim, and a small call agreeing with
it is confirmed on its merits. Applying the rule to affirmations as well would quietly
suppress true positives and depress sensitivity.

**This is the only exclusion in the design that flatters the caller.** Everything removed
here is something nobody can hold against it, so it is a named term rather than a silent
filter, and the report states the base count and the number of affected calls in words.
Specificity read from such an evaluation has to be read together with it: it describes what
the truth was able to see, not what is there. Leaving the false positives in would not
avoid the problem — it would move it to the other side of the ledger and call a method
wrong for seeing something the truth was never able to look for.

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
| Partition accounting | bp | no | How much genome each number applied to (reconciles exactly) |
| Paired method comparison | events | McNemar exact p | The only honest basis for choosing between two methods |

### Confidence intervals are event-level only

Wilson score intervals are reported on **event counts**, never on base-pair proportions.
Base pairs inside one segment are not independent observations; a binomial interval over
900 million bases would be absurdly narrow and would assert a precision that does not
exist. Base-level figures are descriptive point estimates and are labelled as such.

The Wilson interval is used rather than the normal approximation because benchmark
proportions routinely sit at 0 or 1 with small denominators, where the normal interval
collapses to zero width.

**Events are not independent either, and the report says so.** Several events routinely
come from one specimen, and they share its purity, library, coverage and artefacts. An
event-level interval therefore describes a population of independent events that does not
exist, and it is narrower than the data support. The aggregate reports the number of
specimens, the largest number of events contributed by one of them, and a flag —
`intervals_are_anticonservative` — that says this in one field. `discordant_specimens` does
the same for McNemar, which treats each discordant pair as an independent coin flip.

A specimen-weighted detection rate is reported alongside the event-level one. It counts a
specimen as a success only when every assessable event in it was detected: a deliberately
crude summary whose purpose is to make the gap between the two numbers visible. It is
**not** the cluster-robust endpoint an analytical validation would pre-specify — that is a
study-design decision (see ADR-020).

### A direction is not a finding

`favours` names a method only when McNemar's exact test is significant at a pre-specified
`alpha`. Which way the counts lean is reported separately as `observed_direction`.

The two are separated because they are different claims and only one needs a test. With
four discordant pairs the smallest attainable two-sided exact p-value is 0.125, so a 4-0
split looks decisive and could never have been significant at any conventional threshold.
`minimum_attainable_p_value` and `underpowered` report exactly that, which is what
separates *we compared them and found no difference* from *this comparison could not have
found one* — readings a bare non-significant p-value cannot tell apart, and the second of
which is a design fault rather than a result.

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

- **Rejected (off-target) reads** still occupy the flow cell and **may** form a
  near-uniform low-coverage whole-genome background. That is the population depth-based
  CNV methods such as ichorCNA, QDNAseq and Spectre assume, which makes these reads the
  most promising candidate basis for genome-wide CNV from an adaptive-sampling run — a
  candidate, not an established one. Whether the assumption holds locally depends on
  uniformity, GC behaviour, mappability and usable genome fraction at the achieved yield.
  **None of that has been measured on local GridION data**, and until it is, every result
  derived from this basis rests on an untested premise.
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

Shipped mappings: generic IGV `SEG`, ichorCNA (`Corrected_Copy_Number` /
`Corrected_Call`, derived from the ichorCNA output documentation), and `QDNASEQ_ACE_MAPPING`
for the table `scripts/run_qdnaseq_ace.R` writes.

The QDNAseq mapping is a narrower claim than the other two. It is not a mapping for QDNAseq
output in general: it describes the columns *this repository's* runner chooses, and it is
shipped for that reason — the layout is ours and its checksum is in every run's provenance,
which is exactly what could not be said when this module was first written. A differently
configured QDNAseq or ACE installation still needs its own mapping.

**Not shipped:** a Spectre mapping, for the original reason — its column layout was not
confirmed against the upstream source, and shipping a guessed mapping as a verified default
is exactly the failure mode the design avoids.

These adapters only *parse*. Execution belongs behind the repository's existing adapter
boundary with version pinning and argument-vector invocation, as done for Sniffles2 and
Mosdepth.

## The QDNAseq lane, measured through this subsystem

The runtime CNV lane is QDNAseq + ACE. It is wired *through* the benchmark architecture
rather than being promoted beside it: `call_set_from_qdnaseq_report()` turns one run's
`*.qdnaseq.json` into a `CnvCallSet`, which is the same contract every other candidate
method is scored under, over the same evaluable-genome mask.

```bash
ontseq cnv-callset-from-qdnaseq results/runs/RUN_001/S1/evidence/cnv/S1.qdnaseq.json \
  --call-set-id QDNASEQ_S1_500KBP \
  --data-basis adaptive_sampling_off_target \
  --reference-lock configs/reference/grch38.lock.json \
  --output results/cnv/S1.qdnaseq.callset.json
```

Four properties of that conversion are worth stating, because each one is a place where a
convenient default would have produced a wrong number quietly:

**`data_basis` has no default.** An adaptive-sampling run contains two read populations
whose depth behaviour is not comparable, and a run that pooled them is a third case. The
caller states which one the report came from; nothing guesses it from the manifest.

**Uncovered regions become declared no-calls.** QDNAseq drops bins it cannot correct, and
the runner keeps only chr1–22. Both limits are invisible in a segment table — the rows that
would say so are simply absent. Given contig lengths, every uncovered base is emitted as a
`no_call_region`, which keeps it out of the denominator instead of being scored as agreement
with whatever the truth set asserts there. Without a reference lock the conversion still
runs, but it says in a warning that it could not do this.

**The version names both packages.** `QDNAseq 1.38.0+ACE 1.20.0`, not `QDNAseq`. QDNAseq
decides the bins and the correction; ACE decides the purity/ploidy fit the absolute copy
numbers are expressed in. Either one moving changes the answer, so a result attributed to
one of them could not be reproduced from its own label.

**The quantitative column wins over the rounded one.** `absolute_copy_number` is used, not
`call`. Rounding is the step at which a shallow gain and a neutral region stop being
distinguishable, and the scorer should see that distinction rather than a band boundary
someone else chose.

What this does **not** do is promote anything. `CnvCallSet` fixes `reportable` to `False`
and no argument changes it. The bin size, the ACE penalty and the ploidy grid arrive from
the run's policy and are recorded as engineering parameters; a benchmark result does not
turn one of them into a validated threshold. Historical values from the laboratory's
previous pipeline — the 1000 kbp lane, ACE penalty 0.6, the 0.66 affected-band threshold —
are reference points for comparison, and must not become production defaults by being
mentioned here. Selecting a preferred CNV method still requires real cohort data, which
this repository does not have.

## Running it

```bash
# Fully synthetic end-to-end benchmark. Simulates depth, runs two configurations of the
# baseline caller on byte-identical data, evaluates and aggregates each, then compares
# them pairwise. Writes cnv-demo.aggregate.{default,conservative}.json and
# cnv-demo.comparison.json.
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

## Observed behaviour of the baseline in CI

`ontseq cnv-demo-benchmark` runs 63 evaluations per caller configuration (7 blast
fractions x 3 coverages x 3 replicates) against a synthetic del(5q) / -7q / +8 / del(20q)
profile, for two configurations of the baseline on byte-identical simulated data. The
figures below are the **default** configuration, from the CI run that introduced this
subsystem:

```
method: ontseq-baseline-readdepth 0.1.0
evaluations: 63
overall detection: 70/252 = 0.278 [0.226, 0.336]
  tumor_fraction= 0.05: detection=0.000 (0/36) base_concordance=0.595
  tumor_fraction=  0.1: detection=0.000 (0/36) base_concordance=0.594
  tumor_fraction= 0.15: detection=0.000 (0/36) base_concordance=0.595
  tumor_fraction= 0.25: detection=0.000 (0/36) base_concordance=0.595
  tumor_fraction=  0.4: detection=0.000 (0/36) base_concordance=0.594
  tumor_fraction=  0.6: detection=0.944 (34/36) base_concordance=0.990
  tumor_fraction=    1: detection=1.000 (36/36) base_concordance=0.999
  LoD95 by tumor_fraction: empirical=none model=withheld
```

Three things in this output are worth reading carefully, because they are the harness
behaving correctly rather than a caller performing well or badly.

**The cliff between 0.4 and 0.6 is a threshold artifact, not a detection limit.** It is
fully explained without reference to noise. The baseline assigns states with a fixed
neutral band of ploidy +/- 0.5, so a heterozygous loss is called only once the mixture
copy number drops below 1.5:

| Blast fraction | Mixture CN of a het loss | Below 1.5? |
| --- | --- | --- |
| 0.40 | 1.60 | no |
| 0.60 | 1.40 | yes |

The baseline is limited by its own rounding rule, not by counting statistics. A method
that estimates purity and ploidy instead of assuming them would not have this cliff. This
is precisely the kind of failure a null model is supposed to expose, and it is the reason
the demo deliberately does **not** tell the caller the tumor fraction: doing so would
invert the mixture by 1/f and make the result depend on one supplied number rather than
on the data.

**Base-level concordance of ~0.595 at low blast fraction is not 60% correctness.** The
truth profile alters roughly 36% of the simulated genome, so a caller that reports
"neutral everywhere" scores about 0.64 by default. Reading base concordance without the
partition accounting would badly overstate performance.

**Both limits of detection are withheld, and that is the correct answer.** The empirical
value is withheld because even 36/36 has a Wilson lower bound near 0.90, which does not
reach the 0.95 target; the model-based value is withheld because the near-separated
design gives the logistic fit no finite maximum. A harness that printed a number here
would be inventing a detection limit the data cannot support.

### The paired comparison shows why "consistent" is not "proven"

The same run compares the default configuration against a conservative one (split
threshold raised from 4.0 to 8.0) on byte-identical simulated data:

```
ontseq-baseline-readdepth vs ontseq-baseline-readdepth-conservative:
252 paired event(s), only-A=4 only-B=0, p=0.1250
```

The default configuration won **every** discordant pair — a perfectly one-sided result.
It is still not significant, and it could not have been: with 4 discordant pairs the
smallest attainable two-sided exact p-value is `2 x (1/2)^4 = 0.125`. At least **6**
discordant pairs are needed before `p < 0.05` is reachable at all.

This is the trap an unpaired comparison hides. Two aggregate detection rates that differ
slightly look like a difference; the paired view shows the entire difference rests on
four events, and that four events cannot settle the question no matter how one-sided they
are. A real method-selection study has to be powered before it is run, not interpreted
afterwards.

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
- the genome partition reconciles exactly, and no-call bases outside the mask are not
  double counted;
- the paired method comparison recovers exact McNemar p-values, excludes events not
  assessable under both methods, and refuses to call full agreement a tie;
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
