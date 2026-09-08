# Modified-base (methylation) lane

## Status

Research use only. The adapter is wired into the canonical run graph and covered by unit tests,
but it has **never been executed against the real modkit binary** in this repository or in CI.
`StageSpec.verification` records this as `unverified_adapter`, and a run whose methylation stage
completes is reported under `UNVERIFIED ADAPTERS COMPLETED`. Nothing in this lane is analytically
or clinically validated.

## Why the lane exists

The information was already in the data and was being thrown away. Dorado writes `MM`/`ML`
modified-base tags, `align.py` deliberately carries them through the FASTQ round trip unchanged,
and `preflight.py` warns when a basecalling policy would produce reads without them. Until now
nothing downstream read them.

This lane reads them. What it does not change is the limitation underneath: CI proves `MM`/`ML`
tags survive alignment, including on reverse-strand records, and proves nothing about whether a
caller interprets them correctly there.

## What it produces

`modkit pileup` over the aligned BAM, normalized into `MethylationReport`
(`schemas/methylation-report.schema.json`) at `evidence/methylation/<sample>.methylation.json`
inside the run envelope. The report holds one row per region and modification code:

- `sites_total` and `sites_at_minimum_coverage` — how much of the region the pileup saw, and
  how much of it was deep enough to use;
- `valid_call_count` / `modified_call_count` — the raw counts every fraction is derived from;
- `mean_modified_fraction` — call-weighted, so a deep site counts for more than a shallow one;
- `median_site_modified_fraction` — unweighted, so a single very deep site cannot dominate;
- `mean_valid_coverage`.

Rows aggregate either over canonical chromosomes (`region_source: chromosome`) or over the
locked target design (`region_source: target_bed`). Adaptive Sampling leaves the off-target
genome at a depth where a chromosome-wide fraction mixes measured targets with barely observed
background, so an enriched run should aggregate over the design.

## Three refusals

**An empty pileup is never reported as "unmethylated".** A BAM basecalled without a
modified-base model carries no `MM`/`ML` tags, and modkit answers that with an empty file — which
looks exactly like a sample with no methylation. The adapter probes the BAM first
(`samtools view -c -e '[MM]'`) and refuses to run when the count is zero, naming re-basecalling as
the remedy. A samtools too old to evaluate a filter expression returns *no answer* rather than
zero: `reads_with_modified_base_tags` stays `null`, the report carries a warning, and the run
continues. Not knowing is a distinct answer from knowing there are none.

**A region with no qualifying site reports `null`, not `0.0`.** A fraction of zero is a
measurement. The absence of a measurement is not one, and the report model enforces that
fractions exist exactly when a site met the coverage floor. A pileup where no site anywhere
reached the floor is `NO_CALL`, with the reason spelled out.

**The confidence threshold is pinned, not estimated.** modkit can pick a filter threshold from
the data in front of it. That makes the parameter a function of the sample, so two runs of the
same pipeline are no longer running the same pipeline. `MethylationPolicy.filter_threshold` is
explicit and versioned; the technical default of `0.8` is an engineering starting point.

A modification code present in the pileup but absent from the policy fails the run rather than
being dropped. A 5mC fraction computed from a model that also emitted 5hmC means something
different from one that did not, and silently discarding the `h` rows would change the meaning
without changing the number's name.

## Configuration

`configs/methylation/modkit.technical.yaml`, profile `modkit-cpg-technical-v1`:

| Field | Default | What it means |
| --- | --- | --- |
| `expected_version` | `0.4.1` | modkit version lock. **Re-pin deliberately** against the installed binary; pileup output is not comparable across modkit majors |
| `modification_codes` | `[m]` | Codes the report may contain |
| `ignored_codes` | `[h]` | Folded into the canonical count via `--ignore` |
| `cpg_only` / `combine_strands` | `true` / `true` | CpG restriction (needs the reference FASTA) and strand folding |
| `filter_threshold` | `0.8` | Per-call confidence floor |
| `minimum_valid_coverage` | `5` | Sites below this are counted but excluded from aggregates |
| `region_source` | `chromosome` | `chromosome` or `target_bed` |
| `verify_modified_base_tags` | `true` | The `MM`/`ML` probe; one extra pass over the BAM |

With `ignored_codes: [h]` the reported 5mC fraction is a *5mC-versus-everything-else* fraction,
not a share of total cytosine modification.

## Running it

Inside the canonical runner, when the manifest requests the module:

```yaml
analysis:
  profile: lcwgs
  modules: [qc, methylation, report]
```

```bash
ontseq preflight <manifest.yaml> --reference-lock <lock.json> --run-id RUN_001 \
  --reference-fasta <ref.fa> --methylation-policy configs/methylation/modkit.technical.yaml

ontseq run <manifest.yaml> --reference-lock <lock.json> --run-id RUN_001 \
  --reference-fasta <ref.fa> --methylation-policy configs/methylation/modkit.technical.yaml
```

Standalone, against an intake artifact that already passed the gate:

```bash
ontseq call-methylation <manifest.yaml> \
  --intake results/intake.json \
  --policy configs/methylation/modkit.technical.yaml \
  --reference-fasta <ref.fa> \
  --output-dir results/methylation \
  --output results/methylation.json
```

The stage is deselectable like any other component: `--without methylation`, or a
`configs/components/` selection naming `modkit` and its version.

## Limits

- Not validated. No threshold, region set or classifier here has analytical or clinical
  performance data behind it.
- Aggregated fractions depend on the basecalling model that produced the tags. Runs basecalled
  with different models are not comparable, and the report says so.
- Strand-folded CpG values combine both strands of one dinucleotide and are not per-strand
  measurements.
- `align.py` already records that modified-base tag interpretation on reverse-strand alignments
  has not been validated against a downstream caller here. That limitation now has a downstream
  caller and still has no validation.
- No read names, per-read modification probabilities or source BAM paths reach the report.
- Differential methylation, methylation-based classification and tissue-of-origin inference are
  all outside this adapter. It produces descriptive fractions and stops.
