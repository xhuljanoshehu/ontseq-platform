# Modified-base (methylation) lane

## Status

Research use only. The adapter is wired into the canonical run graph and targets modkit
0.6.x semantics: the pileup declares the reported modifications explicitly via
`--modified-bases`, which requires the locked reference FASTA. Continuous integration runs
the real pinned modkit 0.6.4 binary against synthetic MM/ML fixtures — a 75% modified
region, a measured zero, low-depth/empty targets and the missing-tag refusal
(`tests/test_modkit_real_tool.py`, enabled by the `local-real-tool-smoke` CI job). That is
tool interoperability, not analytical recovery or biological validation.
`StageSpec.verification` records this as `verified_with_real_tool`. Nothing in this lane is
analytically or clinically validated.

## Why the lane exists

The information was already in the data and was being thrown away. Dorado writes `MM`/`ML`
modified-base tags, `align.py` deliberately carries them through the FASTQ round trip
unchanged, and `preflight.py` warns when a basecalling policy would produce reads without
them. Until now nothing downstream read them.

This lane reads them. What it does not change is the limitation underneath: CI proves
`MM`/`ML` tags survive alignment, including on reverse-strand records, and proves nothing
about whether a caller interprets them correctly there.

## What it produces

`modkit pileup` over the aligned BAM, normalized into `MethylationReport`
(`schemas/methylation-report.schema.json`) at `evidence/methylation/<sample>.methylation.json`
inside the run envelope. The report holds one row per region and modification code:

- `sites_total` and `sites_at_minimum_coverage` — how much of the region the pileup saw,
  and how much of it was deep enough to use;
- `valid_call_count` / `modified_call_count` — the raw counts every fraction is derived
  from;
- `canonical_call_count` / `other_mod_call_count` — the two count classes that complete
  the valid denominator (`N_valid = N_mod + N_canonical + N_other_mod`, enforced before
  aggregation);
- `fail_call_count`, `nocall_call_count`, `delete_call_count`, `diff_call_count` — calls
  that did **not** enter the denominator, reported so that "not measured" can never be
  read as "not modified";
- `mean_modified_fraction` — call-weighted, so a deep site counts for more than a shallow
  one;
- `median_site_modified_fraction` — unweighted, so a single very deep site cannot
  dominate;
- `mean_valid_coverage`.

Rows aggregate either over canonical chromosomes (`region_source: chromosome`) or over the
locked target design (`region_source: target_bed`). Adaptive Sampling leaves the off-target
genome at a depth where a chromosome-wide fraction mixes measured targets with barely
observed background, so an enriched run should aggregate over the design.

## modkit 0.6.x semantics

Version 0.6.0 removed `--ignore` and changed the threshold algorithm. The lane therefore:

- declares the reported modifications explicitly via `--modified-bases`
  (`m` → `5mC`, `h` → `5hmC`, `a` → `6mA`); a code outside the policy is not tabulated
  at all, and a code in the pileup but not in the policy still fails the run;
- requires the locked reference FASTA for **every** methylation run, because
  `--modified-bases` requires it — not only for CpG-restricted runs;
- reports 5mC and 5hmC as separate rows sharing the valid-call denominator. No
  probability is folded from one class into another, and a model that emits only 5mC
  simply produces no 5hmC rows;
- drops the removed `--only-tabs` flag; the parser still tolerates the space-delimited
  trailing-count layout older builds emitted.

A re-pin across modkit majors is a semantic migration, not a version-string edit. The
0.4.1 adapter used `--ignore h`. At the probability-transformation stage, this
added half of the 5hmC probability to canonical C and half to 5mC. The effect on final
hard-call counts and reported fractions depended on competing probabilities and
thresholds; it was not a universal fixed increase in the final 5mC fraction. That option
no longer exists upstream, and the lane no longer offers it.

## Five refusals

**Known modkit 0.6.4 same-base MM-group risk fails closed.** The pinned
release can silently lose or misassign calls when one read carries independent modification
groups for the same canonical cytosine (for example separate `C+m` and `C+h` groups). ONTSeq
checks the complete BAM with the pinned samtools filter-expression engine before pileup and
refuses this representation. A single multi-code group is not rejected. This is an upstream
tool limitation, not a methylation threshold.

**A nominally successful pileup with failed records is not accepted.** modkit 0.6.4 can exit
zero while its log reports records that failed processing. ONTSeq inspects that diagnostic;
any non-zero `failed processing` count deletes the partial bedMethyl result and fails the
stage.

**An empty pileup is never reported as "unmethylated".** A BAM basecalled without a
modified-base model carries no `MM`/`ML` tags, and modkit answers that with an empty file —
which looks exactly like a sample with no methylation. The adapter probes the BAM first
(`samtools view -c -e '[MM]'`) and refuses to run when the count is zero, naming
re-basecalling as the remedy. A samtools too old to evaluate a filter expression returns
*no answer* rather than zero: `reads_with_modified_base_tags` stays `null`, the report
carries a warning, and the run continues. Not knowing is a distinct answer from knowing
there are none.

**A region with no qualifying site reports `null`, not `0.0`.** A fraction of zero is a
measurement. The absence of a measurement is not one, and the report model enforces that
fractions exist exactly when a site met the coverage floor. A pileup where no site anywhere
reached the floor is `NO_CALL`, with the reason spelled out. Failed and no-call counts are
kept separate from canonical counts, so a site whose calls all failed its threshold is
`NO_CALL`, never a measured zero.

**The confidence threshold is pinned, not estimated.** modkit can pick a filter threshold
from the data in front of it. That makes the parameter a function of the sample, so two
runs of the same pipeline are no longer running the same pipeline.
`MethylationPolicy.filter_threshold` is explicit and versioned; the technical default of
`0.8` is an engineering starting point.

## Configuration

`configs/methylation/modkit.technical.yaml`, profile `modkit-cpg-technical-v1`:

| Field | Default | What it means |
| --- | --- | --- |
| `expected_version` | `0.6.4` | modkit version lock. **Re-pin deliberately** against the installed binary; pileup output is not comparable across modkit majors |
| `modification_codes` | `[m, h]` | Codes declared via `--modified-bases`; one report row per code |
| `cpg_only` / `combine_strands` | `true` / `true` | CpG restriction and strand folding |
| `filter_threshold` | `0.8` | Per-call confidence floor |
| `minimum_valid_coverage` | `5` | Sites below this are counted but excluded from aggregates |
| `region_source` | `chromosome` | `chromosome` or `target_bed` |
| `verify_modified_base_tags` | `true` | The `MM`/`ML` probe; one extra pass over the BAM |

`configs/methylation/modkit-targets.technical.yaml` is the Adaptive-Sampling counterpart
with `region_source: target_bed`.

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
- The real-binary CI lane proves selected interoperability cases against pinned modkit
  0.6.4. It does not establish correctness for every valid MM/ML representation, and
  current upstream 0.6.4 limitations require the explicit same-base-group and failed-
  processing refusals above. It proves nothing about analytical recovery on biological
  data.
- The pileup tabulates only the modification codes the policy declares. Modifications
  outside the policy are not tabulated and cannot be discovered from this report.
- Aggregated fractions depend on the basecalling model that produced the tags. Runs
  basecalled with different models are not comparable, and the report says so.
- Strand-folded CpG values combine both strands of one dinucleotide and are not
  per-strand measurements.
- `align.py` already records that modified-base tag interpretation on reverse-strand
  alignments has not been validated against a downstream caller here. That limitation now
  has a downstream caller and still has no validation.
- No read names, per-read modification probabilities or source BAM paths reach the report.
- Differential methylation, methylation-based classification and tissue-of-origin
  inference are all outside this adapter. It produces descriptive fractions and stops.
