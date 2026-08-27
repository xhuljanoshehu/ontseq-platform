# Structural-variant evidence and review layer

## Scope and safety boundary

This layer converts raw Sniffles2 and cuteSV calls into a compact, traceable review queue. It is
technical triage only. Caller agreement, a high score and an AML-pattern match are not analytical
or clinical validation. The implementation therefore keeps `reportable=false` for every event and
never turns two breakpoint-overlapping genes into a confirmed fusion.

The status fields are deliberately separate:

1. `detected` — at least one normalized caller record;
2. `technically_supported` — technical evidence such as multi-caller concordance;
3. `biologically_prioritized` — a versioned AML pattern matches the annotated genes;
4. `analytically_validated` — reserved for future locked benchmark acceptance;
5. `reportable` — reserved for the future controlled clinical release policy.

Only the first three may currently be produced. Tumor-only BAM data cannot establish whether an
event is somatic, germline or constitutional without an appropriate comparator and intended-use
validation.

## Productive flow

```text
aligned BAM
  -> Sniffles2 2.8.0 + cuteSV 2.1.3
  -> caller-specific normalization
  -> build-aware breakpoint clustering and caller consensus
  -> genes, cytobands and technical context
  -> Adaptive Sampling breakpoint observability
  -> local AML rearrangement-pattern lookup
  -> versioned, explainable technical score
  -> priority review queue + complete technical table
```

The run envelope retains both VCFs, both normalized caller reports and the consolidated report:

```text
evidence/sv/<sample>.sniffles.vcf
evidence/sv/<sample>.sniffles.json
evidence/sv/<sample>.cutesv.vcf
evidence/sv/<sample>.cutesv.json
evidence/sv/<sample>.consensus.json
```

No read names are copied into normalized or reviewer-facing artifacts.

## Locked technical policies

| Policy | Purpose | Current status |
| --- | --- | --- |
| `configs/sv/sniffles2.conservative.technical.yaml` | Sniffles2 execution and normalization | technical defaults only |
| `configs/sv/cutesv.conservative.technical.yaml` | cuteSV execution and normalization | technical defaults only |
| `configs/sv/sniffles2_cutesv.consensus.technical.yaml` | breakpoint matching and within-caller deduplication | technical defaults only |
| `configs/sv/evidence-priority.technical.yaml` | review score thresholds and weights | technical defaults only |

The consensus matcher distinguishes DEL/DUP/INV/INS from BND/TRA. It uses breakpoint distance,
reciprocal overlap, length similarity and orientation when both calls expose it. Translocation
breakends are canonicalized, so `chr2 -> chr21` and `chr21 -> chr2` can represent the same event.
The configured 500-bp distance is an engineering default, not a validated tolerance.

Within-caller representations may be clustered, but only distinct caller names create caller-
consensus evidence. This prevents three nearby Sniffles2 records from masquerading as three
independent confirmations. The B418-like regression fixtures are synthetic coordinate-only cases;
they do not claim a biological identity for the original specimen.

## Gene and cytoband resources

Reference resources are intentionally not committed. Download them into an institution-approved
local reference directory, retain the original file, then normalize and checksum-lock it with
`scripts/prepare_sv_annotations.py`.

Recommended primary sources:

| Build | Genes | Cytobands |
| --- | --- | --- |
| GRCh37/hg19 | [GENCODE release 19](https://www.gencodegenes.org/human/release_19.html), `gencode.v19.annotation.gtf.gz` | [UCSC hg19 `cytoBand.txt.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz) |
| GRCh38/hg38 | [GENCODE release 50](https://www.gencodegenes.org/human/release_50.html), `gencode.v50.annotation.gtf.gz` | [UCSC hg38 `cytoBand.txt.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz) |

Example for GRCh37:

```bash
python scripts/prepare_sv_annotations.py gencode.v19.annotation.gtf.gz genes.grch37.tsv \
  --resource-type genes --resource-id gencode-v19-genes \
  --source-name GENCODE --release v19 --genome-build GRCh37

python scripts/prepare_sv_annotations.py cytoBand.txt.gz cytobands.grch37.tsv \
  --resource-type cytobands --resource-id ucsc-hg19-cytobands \
  --source-name UCSC --release hg19 --genome-build GRCh37
```

The generated `*.lock.json` records the original and normalized SHA-256 values, build, release and
coordinate system. The pipeline refuses a checksum mismatch or a GRCh37/GRCh38 mismatch. GTF
coordinates are converted from 1-based inclusive to 0-based half-open; interval inputs must
already be zero-based half-open.

Both breakpoints retain their own `Locus.gene` and cytoband fields. Intergenic breakpoints keep a
nearest-gene distance in notes; the software does not relabel the breakpoint as genic.

## Repeat, blacklist and mappability context

Any checksum-locked four-column interval resource can be supplied as a repeated
`--sv-context-resource DATA LOCK` pair. Supported types are `repeatmasker`, `tandem_repeat`,
`segmental_duplication`, `blacklist`, `mappability`, `centromere` and `telomere`. Appropriate
publisher-native UCSC/ENCODE files must first be converted to the normalized four-column contract
`chrom, start, end, label`; the lock records that normalized artifact and the source hash.

Overlaps add side-specific flags such as `primary:segmental_duplication`. Calls are not deleted.
Each distinct context flag applies the penalty declared in the evidence policy, capped by that
policy's maximum. A strongly supported event therefore remains visible with its caveat.

## Adaptive Sampling observability

For Adaptive Sampling, the layer consumes the existing locked Mosdepth target-coverage report and
stores the mean depth at each breakpoint. It distinguishes:

- `OBSERVED_ADEQUATELY`
- `PARTIALLY_OBSERVED`
- `INSUFFICIENT_COVERAGE`
- `OUTSIDE_TARGET`
- `NOT_APPLICABLE`

The default 10x mean-depth floor is explicitly unvalidated and versioned in run parameters. A
buffered selection panel can never produce `OBSERVED_ADEQUATELY`; the unbuffered analysis ROI is
required for that label. Absence of a call outside or below observable target space is never
rendered as a negative biological result. Relative whole-genome background coverage is not yet a
validated gate and is not inferred.

## AML knowledge layer

`configs/knowledge/aml_rearrangements.v0.1.json` is a small, local, checksum-locked prioritization
resource. It covers exact recurrent pairs and open-partner patterns for RUNX1::RUNX1T1,
CBFB::MYH11, PML::RARA, KMT2A, NUP98, DEK::NUP214, RBM15::MRTFA, BCR::ABL1 and MECOM.

Its initial source frame is the
[WHO 5th edition myeloid classification](https://doi.org/10.1038/s41375-022-01613-1) and the
[2022 International Consensus Classification](https://doi.org/10.1182/blood.2022015850). This
curated list is a software prioritization aid, not a complete disease ontology and not a
substitute for current professional classification, transcript/orientation analysis or expert
review. A match is named `known_rearrangement_pattern`; neither `fusion_supported` nor
`fusion_validated` is asserted.

## Report and provenance

The report starts with raw-record, consolidated-event and confidence-tier counts. High and
moderate events appear in a filterable priority table with Gene A/B, cytobands, per-caller support,
caller consensus, breakpoint depth, observability, technical flags, AML relevance, confidence and
validation status. The complete technical table remains below it and is independently filterable.
JSON and XLSX retain source record IDs and full caller evidence.

Parameters, tool versions, reference/resource fingerprints and policy profile IDs participate in
the stage signature. Changing any of them invalidates content-addressed resume and reruns the SV
stage.

## Validation path

1. Synthetic fixtures verify parsing, canonicalization, clustering, annotation and rendering.
2. Version-pinned HG002/HG008 datasets assess transferable technical behavior using a format-aware
   comparator such as Truvari where appropriate.
3. Orthogonally characterized AML specimens provide intended-use truth from karyotype, FISH,
   PCR/RT-PCR, RNA fusion testing or another validated method.
4. Replicates and dilution/coverage series characterize precision, LoD and no-call behavior.
5. Acceptance criteria are predeclared per build, assay, SV type, coverage and tumor/blast fraction.

The existing normalized-event benchmark contract records TP, FP, FN, precision, recall, F1 and
breakpoint error. Separate cases carry coverage, VAF, assay and event-type strata. It is an
engineering framework; no acceptance cut-off is encoded.

## Known limitations and validation impact

- No local AML truth cohort or assay-specific acceptance criteria are included.
- Orientation is used only when caller VCF evidence exposes a supported strand field; transcript
  structure and in-frame fusion support are not inferred.
- Context resources are user-supplied and must be converted and locked locally.
- The AML list is deliberately small and requires controlled expert maintenance.
- A technical confidence tier may change the review order, never `reportable` or ISCN release.

This change alters biological prioritization and report order. It therefore requires the added
regression suite, schema/provenance updates, changelog entry and validation-impact review before a
validated release. Until the studies above pass, the layer remains research-only.
