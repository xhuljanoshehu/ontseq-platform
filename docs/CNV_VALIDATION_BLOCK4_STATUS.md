# CNV Validation Block 4 — QDNAseq+ACE Multi-Resolution Adapter Status

## Status

Block 4 implements the Research Use Only QDNAseq+ACE adapter that maps the locked
100/500/1000-kbp caller outputs into the Block-2 full-evidence and traceability contracts.

The adapter is additive. It does **not** alter QDNAseq/ACE execution, change the technical
primary resolution, tune caller thresholds, select a production caller, or establish analytical
or clinical validity.

## Native artifact retention

Every path declared by `QDNAseqCallReport.output_files` is retained as a
`CnvNativeArtifactReference` with:

- deterministic artifact ID;
- caller-native relative path;
- SHA-256 checksum;
- byte size;
- role;
- media type.

Paths must remain safely inside the supplied QDNAseq output directory. Missing files, unsafe
paths and duplicate declared outputs fail closed.

The run-summary full-evidence record references the complete retained artifact set, so sealing
the evidence manifest cannot silently orphan a caller-native file.

## Caller and provenance lock

The adapter verifies before mapping evidence that:

- caller ID is `qdnaseq_ace`;
- report specimen ID matches the registered validation specimen;
- genome build matches the registered specimen;
- report bin resolutions equal the locked QDNAseq policy resolutions;
- the selected primary fit resolution equals the locked primary resolution;
- QDNAseq and ACE versions are explicitly present;
- the caller-version lock equals the observed QDNAseq/ACE provenance.

The full caller policy snapshot and its canonical content SHA-256 remain attached to every
full-evidence record.

## Multi-resolution fit retention

For every registered resolution, including 100, 500 and 1000 kbp, Block 4 retains:

- the selected ACE fit;
- all retained ACE alternative/candidate fits;
- cellularity;
- ploidy;
- fit error;
- candidate count;
- segment count;
- references to the corresponding caller-native fit artifacts.

The current technical primary resolution remains 500 kbp only because that value is already
locked in the QDNAseq policy. The adapter does not infer or select a new preferred resolution.

Selected non-primary resolutions and retained ACE alternatives remain secondary evidence.
They are not deleted or overwritten.

## Caller-native bin evidence

Where a QDNAseq bin table is emitted, each caller-native row becomes a
`CALLER_BIN` full-evidence record.

The adapter preserves the exact zero-based half-open locus and available finite caller-native
quantitative values such as read count, GC, mappability, blacklist, residual, loess, corrected,
copy-number and segmented values when those columns are present.

The adapter does not invent missing numerical values.

## Caller-native segment evidence

Every retained QDNAseq+ACE segment row becomes a `CALLER_SEGMENT` full-evidence record at
its native 100/500/1000-kbp resolution.

Retained fields include:

- exact zero-based half-open interval;
- absolute copy number;
- segment bin count;
- ACE call value;
- finite `qnorm_log10` when numeric.

The known caller-native `-Inf` qnorm sentinel is left in the immutable native TSV rather than
converted into a non-finite structured numeric value.

Primary-resolution segments are marked as primary-analysis evidence; non-primary resolution
segments remain secondary evidence.

No new event class is inferred from a segment row by this adapter.

## Chromosome-level evidence

Per-resolution chromosome copy-number summaries are retained independently at 100, 500 and
1000 kbp.

The cross-resolution chromosome consensus already produced by the locked QDNAseq runtime is
also mapped into full evidence with:

- median copy number;
- rounded copy number;
- agreeing resolution count;
- contributing resolution count;
- minimum copy number;
- maximum copy number.

The consensus remains secondary evidence and does not replace any per-resolution values.

## Normalized-event traceability

The existing QDNAseq report's normalized primary CNV events remain additive derived records.

For each positive normalized event, Block 4 requires exactly one observed, locus-bearing
primary-resolution caller segment with the same interval and compatible copy number.

The resulting `CnvNormalizedEventRecord` retains:

- the primary 500-kbp resolution;
- the exact source full-evidence record ID;
- event type;
- interval;
- normalized copy number;
- caller/build/specimen/reference/input identity;
- adapter/normalization policy SHA-256 provenance.

If a normalized event cannot be mapped to exactly one primary caller segment, manifest
construction fails closed instead of guessing a source.

## Sealed evidence manifest

`build_qdnaseq_ace_validation_manifest` seals:

1. all retained caller-native artifacts;
2. all multi-resolution full-evidence records;
3. all additive normalized primary events.

The existing Block-2 manifest contract enforces:

- unique native artifact IDs;
- unique full-evidence analytical addresses;
- no orphan native artifacts;
- no normalized event with a missing full-evidence source;
- identity consistency from normalized event back to source evidence;
- primary normalized events sourced from primary full evidence;
- canonical content SHA-256 locking;
- `retain_all_evidence = true`.

## TDD evidence

Block 4 was implemented in RED → GREEN slices.

### Slice 1 — retained fits and artifacts

RED required:

- all 100/500/1000-kbp resolutions;
- selected ACE fit plus retained alternative fit per resolution;
- complete native artifact retention.

The initial RED failure was the intentionally missing QDNAseq validation adapter module.

### Slice 2 — caller-native row evidence

RED required:

- `CALLER_BIN` evidence at all three resolutions;
- `CALLER_SEGMENT` evidence at all three resolutions;
- per-resolution chromosome summaries;
- cross-resolution chromosome consensus.

The RED assertion showed the expected empty set of bin-resolution evidence before row mapping
was implemented.

### Slice 3 — sealed manifest and fail-closed traceability

RED required:

- a sealed evidence manifest;
- normalized primary-event source linkage;
- complete artifact reachability;
- failure when a declared native artifact is missing.

The RED failure was the intentionally missing
`build_qdnaseq_ace_validation_manifest` API.

## Explicit non-claims

Block 4 does not:

- execute QDNAseq+ACE on real patient specimens;
- establish analytical sensitivity, specificity, LoD or clinical validity;
- select QDNAseq+ACE as the production CNV caller;
- compare QDNAseq+ACE against ichorCNA or Spectre;
- tune registered thresholds or bin-size rules after observing performance;
- promote an ACE alternative fit to the selected fit;
- discard 100-kbp or 1000-kbp evidence;
- change the canonical production workflow;
- merge PR #79;
- change the package version.

Package version remains `0.8.2`.

Research Use Only. Human review remains required.
