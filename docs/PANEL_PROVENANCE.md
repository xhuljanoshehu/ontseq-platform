# GRCh38 Adaptive Sampling panel provenance

## Active bundle

The manifest at `configs/panels/AML_AS_111_GRCh38_v1/bundle.yaml` is the source definition
for the GRCh38 AML Adaptive Sampling panel. Loose panel files elsewhere under `configs/`
are legacy compatibility artifacts and are not registry-activated resources.

The bundle contains two immutable laboratory source artifacts:

| Artifact | Actual format | SHA256 |
| --- | --- | --- |
| `250611_fusion_panel_with_buffer.bed` | four-column, 1-based inclusive target table | `f454644f18d8728c03678f4c6e969da7067879367c894c274e8c44be9352ef7e` |
| `250611_fusion_panel_with_buffer.interval_list` | plain `chr:start-end` lines, not Picard/GATK IntervalList | `f9ebbfbaa555b05d42fdd9edfda93eb9556661a8c940dbbd64b938769c40b441` |

Both sources contain the same 111 numeric intervals in the same order. They contain target
coordinates and labels only; no sample, patient, run, or read data is committed.

## Coordinate contract

The `.bed` suffix does not establish BED semantics. Comparison with GRCh38 gene locations
and the exact 10 kb flanks establishes that both laboratory artifacts record 1-based
inclusive gene coordinates with a buffer. The source bytes are never used directly by an
analysis.

`panel_bundle.import_panel_sources` applies exactly one conversion:

```text
source [start, end] -> active [start - 1, end)
```

The resulting `derived/selection_panel.normalized.bed` is a standard 0-based half-open BED:

- 111 intervals;
- 17,028,488 total bases;
- one-base reduction of every source start;
- unchanged source end and interval order;
- build `GRCh38`;
- role `selection_panel_buffered`.

The original 1-based source and normalized derivative have separate checksums and resource
IDs. No liftover or GRCh37 fallback exists.

## Build and role evidence

The coordinates match GRCh38 rather than GRCh37, including loci such as `ABL1`, `ALK`, and
`RUNX1T1`. The interval boundaries are approximately 10 kb outside the named genes, so this
resource describes the regions selected by Adaptive Sampling. It does not describe an
unbuffered analysis ROI and cannot be used as one.

The Analysis ROI and panel transcript cache are generated only after the locked
`GRCh38_GENCODE50_MANE1.5_v1` SQLite annotation cache is installed:

- a target label must resolve to exactly one GENCODE gene on its declared source chromosome;
- its Analysis ROI is the GENCODE gene body, not an interval inferred by trimming 10 kb;
- every transcript is retained and ranked deterministically;
- MANE Select precedes MANE Plus Clinical, canonical/APPRIS principal,
  protein-coding/basic, and all remaining transcripts;
- CDS length, transcript length, and transcript ID break ties.

Pending generated resources keep the panel bundle incomplete and therefore unresolvable for
analysis. The registry activates it only after the derivative paths and checksums exist.

## Unresolved `IGH` source row

The source label `IGH` is paired with `chr5:143396959-143417420`. IGH is located on chr14q32
in GRCh38. The software does not reinterpret the chr5 coordinate as IGH or as a different
gene.

The source label remains byte-for-byte unchanged in source provenance. Every active
derivative uses `IGH_REVIEW_REQUIRED`; the transcript/ROI compiler omits it. This means:

- no invented IGH Analysis ROI;
- no negative observability statement for this unresolved target;
- no automatic correction to `NR3C1` or another chr5 label.

Resolution requires an authoritative corrected laboratory target definition.

## Other generated open issues

Compilation against the complete GENCODE 50 cache records every unresolved source label in the
activated panel manifest. In addition to `IGH_REVIEW_REQUIRED`, the current official release
leaves `GPR128` and `MKL1` unresolved:

- the source `GPR128` interval is on chr17 and does not match the current GENCODE `ADGRG7`
  (`GPR128`) locus on chr3;
- `MKL1` is a historical symbol for `MRTFA`, but ONTSeq does not silently substitute gene aliases
  without a curated, versioned alias resource.

`P2RY8` occurs on both chrX and chrY in the pseudoautosomal region. Its source row explicitly
declares chrX, so label plus source chromosome resolves it uniquely to the chrX GENCODE record.
The generated open-issue list is provenance, not a negative biological finding; unresolved rows
produce neither Analysis ROI nor observability claims.

## Selection coverage versus Analysis ROI

Selection coverage measures whether enrichment worked across the buffered sequencing
design. Analysis-ROI coverage measures observability over resolved GENCODE gene bodies.
They are separate resource roles, separate files, and separate report concepts. A missing or
unresolved Analysis ROI is not a biological negative and must remain visible as unresolved.

Historical values in `configs/qc/target_coverage_expectations.grch38.tsv` remain descriptive
technical observations only. They do not establish an adequacy, reportability, or no-call
threshold.

## Validation impact

The coordinate correction changes every active interval start by one base and changes the
total target span by 111 bases compared with treating the source as standard BED. Panel
coverage and breakpoint-in-target boundary results can therefore change at interval edges.

The new ROI/transcript compiler and CNV cytoband summarizer can change gene/transcript/band
annotations and review order. Tests lock the conversion, source checksums, unresolved target
behavior, transcript ranking, minus-strand exon/intron logic, CDS phase, the 66% cytoband
threshold, centromere separation, and whole-chromosome handling. These are deterministic
software contracts, not validated analytical performance claims.
