# GRCh38 and GRCh37/hg19 Adaptive Sampling panel provenance

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

## GRCh37/hg19 build-time projection authority

`configs/panels/AML_AS_111_GRCh37_v1` is a separate GRCh37 panel bundle. It does not relabel the
GRCh38 coordinates and it does not perform liftover during an analysis. Its checked-in source
and derivatives reproduce this controlled maintainer operation:

1. The 111-row GRCh38 laboratory design was normalized to 0-based half-open BED
   (`SHA256 9338e66a52310d4181e2eb121d2f3954c16bb7bf482a8c5d2ec3c7b1e5381ee1`).
2. The official UCSC `hg38ToHg19.over.chain.gz` snapshot was used with
   `liftOver -minMatch=0.99`; the chain had
   `SHA256 14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1`.
3. The Linux `liftOver` executable used for that controlled build had
   `SHA256 6529f91135d579f18d1d86681810cdd0150bf6d7531457baa79a6632706368c0`.
4. The raw result contains 110 intervals / 16,951,494 bases
   (`SHA256 909b9552e4ee50daa04a5e22934ef425833b44000e1086c7d92595974a53ec5d`).
   The unmapped output contains only the split `CT45A2` row
   (`SHA256 55038b8ce5b8aa4600e85e84e9cec23a108ad3da12bd444c8706ccc6d6ecf9e1`).
5. The official 227,698-byte UCSC `hg19ToHg38.over.chain.gz` snapshot
   (`MD5 35887f73fe5e2231656504d1f6430900`, independently matched to UCSC `md5sum.txt`;
   `SHA256 5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d`) was used for a
   reverse `liftOver -minMatch=0.99` check. Its raw mapped result has 109 rows / 3,274 bytes
   (`SHA256 a12f763f3ddbca94e28a29e432658755df553612cb454f2e9dc10593486c15bc`) and each returns to
   the exact original GRCh38 chromosome/start/end. The 44-byte reverse-unmapped output
   (`SHA256 c1b5bc5483cc78daf87f15d46174fe132444943d65ef95d1c5b864726a31d483`) contains only
   `#Split in new` plus `ACACA`. ACACA is retained for GRCh37 selection/ROI analysis but marked
   `ROUNDTRIP_REVIEW_REQUIRED`; it is not counted as reciprocally validated.

The bundle stores the normalized source, both checksum-pinned forward outputs, both checksum-pinned
reverse-audit outputs, the typed mapping lock and the final selection authority. Activation binds
the lock to the declared source bundle/resource and verifies every retained audit artifact. It
intentionally does **not** store or redistribute either chain or the executable. Rebuilding
requires independently obtaining the files from the
[UCSC hg38 liftOver directory](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) and
verifying the recorded hashes before executing the locked command. UCSC describes non-commercial
and commercial licensing separately on its [license](https://genome.ucsc.edu/license/) and
[store](https://genome.ucsc.edu/store.html) pages; institutional/commercial use must satisfy the
applicable terms. Those pages do not explicitly settle redistribution of derived coordinate BEDs,
so open or commercial distribution requires written UCSC clarification (or regeneration from an
appropriately licensed independent mapping authority) plus confirmation of rights in the original
laboratory panel.

### Independent NCBI audit

As an independent check, the NCBI Release 3.0 GRCh38.p14-to-GRCh37.p13 alignment snapshot was
inspected without making it the production build path. The official files were:

| NCBI artifact | SHA256 |
| --- | --- |
| `GCF_000001405.40-GCF_000001405.25.gff` | `8f487ba51a0a02de02187979cafd3f710094895737005611291ead3a399c52c2` |
| `GCF_000001405.40-GCF_000001405.25.report` | `2b229446b09d219a990c90ba398b4be9975cc4de9ab0898f79146732c179fbcb` |
| `GCF_000001405.40-GCF_000001405.25.asn` | `73c934b45d8a0e37af68259283d602edd803226bc4fcdc657bf1f99537c525f1` |

The [NCBI Release 3.0 directory](https://ftp.ncbi.nlm.nih.gov/pub/remap/Homo_sapiens/3.0/GCF_000001405.40_GRCh38.p14/GCF_000001405.25_GRCh37.p13/)
confirms `.40` as query and `.25` as target. Its alignments independently confirm that CT45A2
does not have a complete primary/Canonical-25 target (its full target is on a GRCh37 patch), so
ONTSeq must not invent a canonical interval. Other loci can require fragment merging or contain
additional patch alignments, reinforcing the need for explicit mapping policy and review rather
than a blanket 111/111 claim.

NCBI retired the Remap web service/API in 2023. The archived standalone clients were also tested:
the Linux binary lacks legacy shared libraries on the current WSL runtime, while the Windows
client cannot read the current v3 ASN alignment (`CSerialException::eOverflow`). Therefore these
clients are not a supported ONTSeq build/runtime dependency. NCBI static data remain useful as
validation evidence only; see the [NCBI retirement notice](https://ncbiinsights.ncbi.nlm.nih.gov/2023/05/04/remap-tool-to-retire/).

## GRCh37 target and ROI resolution

The projected selection intervals retain the original target labels and remain the coverage
authority. ROI/transcript compilation is independent and must resolve against the installed
native `GRCh37_GENCODE19_HG19_v2` annotation cache. A versioned target mapping fixes historical
GENCODE-19 symbols/stable IDs; it does not move an interval or override a chromosome conflict.

- `CT45A2_REVIEW_REQUIRED`: no single projected interval exists, so it is absent from selection
  and ROI.
- `IGH_REVIEW_REQUIRED`: the mapped selection row remains on chr5 while the IGH locus is on chr14;
  no ROI is fabricated.
- `GPR128_REVIEW_REQUIRED`: the mapped selection row remains on chr17 while the GENCODE-19 gene is
  on chr3; no ROI is fabricated.
- `ACACA_ROUNDTRIP_REVIEW_REQUIRED`: the forward GRCh37 interval and native ROI remain usable for
  engineering analysis, but the reverse chain splits and therefore does not prove exact reciprocal
  boundary conservation.

The resulting contract is therefore 111 intended source targets, 110 selection intervals and
108 native GENCODE-19 ROI/transcript targets; 109 selection intervals are reciprocal-exact. These
counts, the three mapping/ROI gaps and the independent ACACA roundtrip flag are validated before
bundle activation and recorded in the run-bound JSON, HTML and XLSX provenance. Neither absence
from the ROI nor low coverage may be reported as a negative biological finding.

## GRCh37 validation impact

The new bundle enables deterministic engineering execution for two GRCh37 Adaptive-Sampling
profiles and may change selection-coverage, observability and downstream review content compared
with an lcWGS run. It does not establish that enrichment performance is equivalent between
GRCh38 and GRCh37/hg19, or validate analytical sensitivity, specificity, LoD, reportability or
clinical use. Tests establish file identity, build/profile isolation, compiler behavior and
fail-closed dictionary handling only.
