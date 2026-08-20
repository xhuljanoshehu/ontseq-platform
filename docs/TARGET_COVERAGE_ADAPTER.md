# Adaptive Sampling target coverage adapter

## Status

Research use only. This adapter provides descriptive technical QC and has not been analytically
or clinically validated.

## Purpose

The adapter answers a narrow question before CNV or fusion interpretation: **was each intended
analysis region actually covered, and to what descriptive depth?** It consumes an aligned BAM
that already passed the repository intake gate plus the assay's versioned target BED.

It deliberately does not infer a biological negative from absent or weak coverage.

## BED roles must remain separate

The manifest field `assay.target_bed` is interpreted by this adapter as the **unbuffered analysis
ROI BED**. A BED used by MinKNOW or another Adaptive Sampling controller may include selection
buffers or other operational transformations. Those operational regions must not silently replace
the analysis ROI because doing so changes denominators and can hide coverage gaps at the intended
clinical/research loci.

The software cannot prove whether a supplied BED was buffered. The file is fingerprinted and the
role `analysis_roi_unbuffered` is recorded so that provenance remains explicit.

## Tool boundary

The implementation is locked to Mosdepth `0.3.14` under the current technical policy. It executes
Mosdepth without a shell and requests:

- no per-base output;
- mean depth for each BED region;
- counts of bases at or above configured thresholds;
- an explicit mapping-quality threshold;
- an explicit excluded-SAM-flag mask.

The default technical thresholds are `1x`, `10x`, `20x`, and `30x`. They are descriptive bins,
not validated adequacy or reportability limits.

## Normalized output

For each exact BED interval the normalized report records:

- chromosome, zero-based half-open start and end;
- the fourth BED column as `region_id` when present;
- mean depth;
- number and fraction of interval bases at each configured threshold.

The summary records interval count, total interval bases, interval-weighted mean depth,
minimum/median/maximum regional mean depth, threshold fractions, and the number of overlapping
BED intervals.

Overlapping BED intervals are allowed but explicitly warned because interval-weighted summaries
then count overlapping genomic bases more than once. The adapter never merges intervals silently.

## Fail-closed conditions

The run is rejected when, among other conditions:

- the manifest is not `aligned_bam` plus `adaptive_sampling`;
- aligned-BAM intake previously failed;
- sample, reference or genome-build provenance disagree;
- the target BED is missing, empty, malformed, duplicated, or uses unsupported contigs;
- Mosdepth version differs from the policy lock;
- Mosdepth output is missing, duplicated, malformed, or does not match every exact BED interval;
- configured threshold columns differ from the observed output;
- a run would overwrite existing raw Mosdepth artifacts.

## Privacy boundary

The normalized JSON does not contain read names, inserted sequences, per-read evidence, or the
source BAM path. The target BED is represented by its version identifier, declared role and SHA256
fingerprint. Raw runtime files remain on the approved execution system and are not Git artifacts.

## What this adapter does not yet do

- infer gene semantics from arbitrary BED labels;
- calculate a validated on-target enrichment ratio;
- determine whether both partners of a candidate fusion are adequately observable;
- establish CNV or fusion reportability;
- define clinical no-call thresholds;
- replace local validation on the actual 111-gene design and AML samples.

These are intentionally separate validation steps so that an unvalidated technical threshold
cannot become a hidden clinical rule.

## Validation plan

Before promotion beyond `technical_defaults_only`, evaluate the locked local target design across:

1. replicate runs and sequencing batches;
2. coverage dilution and run-yield strata;
3. target-region size and genomic context;
4. known positive fusion partner pairs;
5. known negative and technical-failure samples;
6. blast/tumor-fraction strata;
7. comparison with orthogonal local assays and the existing analysis pipelines.

Acceptance criteria must be defined before looking at the validation results.
