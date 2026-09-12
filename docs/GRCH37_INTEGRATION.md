# Native GRCh37.p13 integration — 0.6.2 engineering candidate

This local executable integration extends the historical PR47 baseline. Every engineering
package records its exact source revision and source-tree inventory; the version label alone is
not provenance. It is not a released GitHub artifact, diagnostic device or clinical validation.

## Implemented contract

| Selection | Native resources | BAM dictionary |
| --- | --- | --- |
| `AML_LCWGS_GRCh37` | GENCODE 19 / GRCh37.p13, native GTF, build-specific cache and cytobands | exact full publisher FASTA, including scaffolds/patches/haplotypes |
| `AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25` | same pinned nuclear GRCh37 / GENCODE 19 resources | exactly chr1–22, chrX, chrY, chrM=16571 in UCSC hg19 order; no extras |
| `AML_AS_111_GRCh37` | same reference/knowledge plus dedicated `AML_AS_111_GRCh37_v1` selection and native GENCODE-19 ROI | exact full publisher FASTA, including scaffolds/patches/haplotypes |
| `AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25` | same build-matched panel/reference/knowledge resources | exactly chr1–22, chrX, chrY, chrM=16571 in UCSC hg19 order; no extras |
| Existing GRCh38 full profiles | GENCODE 50 / GRCh38.p14 + MANE 1.5 | existing full Primary-Assembly lock |
| Existing GRCh38 Canonical-25 profiles | same pinned GRCh38 resources | exactly chr1–22, X, Y, M |

Build choice propagates through Desktop, resource registry, BAM intake, CLI/service,
manifest, breakpoint/CNV annotation, result provenance and reports. Mismatches fail
before analysis; there is no retry against another genome and no automatic liftover.
Missing resources for an unselected build do not block a ready selected profile.

The GRCh37 family includes checksum-pinned native FASTA/GTF, generated FASTA index,
reference dictionary/lock, SQLite annotation cache, named nuclear cytobands, nuclear
repeat/simple-repeat/segmental-duplication tracks, blacklist and uniquely-mappable
context. Source URLs, checksums and license/terms caveats live alongside the
[recipe](../configs/reference_bundles/GRCh37_GENCODE19_HG19_v2/README.md).
The UCSC hg19 choice is an explicit dictionary contract, not an alias. Its 16,571-base
mitochondrial record is intentionally different from the 16,569-base record in the native
GENCODE 19 reference. The software derives the effective BAM lock only for this named profile;
it does not alter the BAM, substitute a FASTA, reuse mitochondrial annotation, or silently
fall back between profiles. Nuclear annotation/context remains pinned to the GRCh37 bundle.

The inherited annotation compiler indexes canonical chromosomes from the native GTF;
the complete 297-contig FASTA dictionary does not imply gene-context coverage on every
patch, scaffold or haplotype. No such coordinates are implicitly projected onto primary
chromosomes.

GRCh37 has no native MANE resource. Native GENCODE transcripts remain available;
absence of MANE is recorded, not replaced by relabeled GRCh38 coordinates. All four GRCh37
profiles remain separate validation lanes; sharing a reference bundle does not weaken their
dictionary or assay contracts.

## GRCh37 Adaptive-Sampling panel

The laboratory supplied only a GRCh38 selection design. Version 0.6.2 therefore preserves its
lineage through a controlled, checksum-pinned **build-time** projection rather than inventing a
different gene-body-plus-flank design. The source normalization contains 111 intervals
(`SHA256 9338e66a52310d4181e2eb121d2f3954c16bb7bf482a8c5d2ec3c7b1e5381ee1`). UCSC
`hg38ToHg19.over.chain.gz` with `liftOver -minMatch=0.99` maps 110 intervals
(`SHA256 909b9552e4ee50daa04a5e22934ef425833b44000e1086c7d92595974a53ec5d`); `CT45A2`
splits in the target assembly and is deliberately omitted rather than guessed.

A separately pinned reverse projection maps 109 of those 110 intervals back to their exact
original GRCh38 chromosome and boundaries. `ACACA` splits on the reverse projection and is
therefore `ROUNDTRIP_REVIEW_REQUIRED`; it remains in the forward GRCh37 selection and native ROI
but is not represented as reciprocally validated. The bundle status remains
`controlled_build_time_unvalidated`.

The mapped selection intervals remain the enrichment-coverage authority. A separate compiler
uses the installed native GENCODE-19 cache and a versioned target-to-stable-gene mapping to create
108 native analysis ROIs/transcript records. The mapped chr5 row labelled `IGH` and chr17 row
labelled `GPR128` contradict their native gene loci and remain `IGH_REVIEW_REQUIRED` and
`GPR128_REVIEW_REQUIRED`, without fabricated ROIs or negative observability statements.
Together with `CT45A2_REVIEW_REQUIRED` and the independent ACACA roundtrip flag, these gaps are
part of the manifest and reportable technical provenance.

The chain and `liftOver` binary are not shipped in the panel bundle, wheel or Desktop runtime.
No coordinate projection runs during BAM processing, no GRCh38 coordinate is queried against a
GRCh37 BAM, and neither profile falls back to the GRCh38 panel. Full command, source, checksum and
licensing provenance is recorded in [`PANEL_PROVENANCE.md`](PANEL_PROVENANCE.md).

`HEMATOLOGY_GRCh37_v1` contains audited coordinate-free gene/association knowledge;
it retains the original curated selection and attribution, including panel-derived
selection bias. It is not exhaustive hematological knowledge. Licensing notes do not
grant rights to redistribute restricted third-party material.

The GENCODE 19 symbol audit found direct matches for 43 of 46 symbols used by the
38 curated rules. Three rules use modern symbols absent under that spelling in this
native cache: ADAD2-AS1, CEP43 and MRTFA. Their native loci exist, but the matcher has
no controlled historical-symbol resolver. Consequently those rules cannot reliably
match by symbol yet; no alias, negative finding or complete knowledge coverage is
inferred. A versioned stable-ID/symbol curation step remains required.

## Real local UI

The new `/workspace` is served by `ontseq serve`, with embedded local assets and a
same-origin session token. It lists installed profiles and permitted BAM/BAI files,
starts the canonical pipeline, polls actual run state and displays real QC/events,
stages and provenance. HTML, XLSX and JSON are downloaded from fixed artifact paths.
Run, sample, build and profile binding is checked. Downloads require a valid result;
path traversal, cross-origin requests and tokenless API calls are rejected.

The prior standalone UX/report examples remain demonstrations. No synthetic finding,
signature, blast fraction, paired normal or executed state is introduced into live UI.
An empty event list is not interpreted as an exclusion of disease. Cooperative
cancellation remains unavailable until the backend can guarantee safe interrupted
state semantics.

Frontend source is in `web/`; run `npm ci && npm run build` there to regenerate the
single packaged `src/ontseq_platform/service/workspace.html`. The build has no runtime
CDN dependency. React/license notices are bundled with the UI.

## Verification and honest boundaries

Unit/contract checks cover both build families, native no-MANE cache lineage, both GRCh37
Adaptive-Sampling dictionary contracts, the 110-selection/108-ROI cardinalities, unresolved
target preservation, cross-build rejection, service authentication/exports and real-client contracts.
Desktop tests cover selected-family readiness and shared service/workspace behavior.
These tests do not replace executing native references and real bioinformatics tools.

`scripts/grch37_profile_smoke.py` resolves and verifies an installed native bundle,
creates a clearly synthetic BAM with the full publisher dictionary, then runs the
actual profile with installed tools. Its artificial depth pattern is chr7 loss and
chr8 gain. It checks QDNAseq/ACE 100/500/1000-kbp agreement, canonical result/report
artifacts, release checksums and content-addressed resume. It never relaxes scientific
thresholds. A sparse fixture may legitimately cause an SV NO_CALL.

QDNAseq stores native bins as one-based, closed intervals. The R export boundary
converts bin and collapsed-segment starts to zero-based, half-open intervals before
Python normalization; ends remain unchanged. Native RDS/ACE calculations are not
shifted. This correction applies to both reference builds. Early local engineering
outputs that retained one-based starts must be regenerated with the corrected script;
their successful execution alone is not final coordinate-validation evidence.

Generated smoke evidence belongs outside the repository. Record the exact output
path, version, tool versions and source-tree checksums when distributing a local
test package. `UNKNOWN`/`LOCAL_WORKTREE` provenance must not be replaced with the
historical base runtime commit. A new clinical/analytical validation cohort,
workstation compatibility tests and any needed licensing review remain separate gates. The
UCSC hg19 Canonical-25 profiles add input-compatibility paths; they do not themselves establish
analytical sensitivity, specificity or clinical validity. The new Adaptive-Sampling profiles
likewise establish deterministic software resolution only. They do not validate enrichment
equivalence across builds, target observability, LoD, sensitivity, specificity or reportability
on clinical specimens.

## Packaging

The engineering bundle uses the existing pinned relocatable tool archive plus a new
Core 0.6.2 wheel. Setup checks both SHA256 hashes, extracts into a new unique prefix,
runs `conda-unpack`, installs the wheel offline without dependency resolution and
requires `ontseq --version` to report 0.6.2 before saving settings. The old installation
is not overwritten. The wheel replaces old configuration and commit metadata; its
source-tree inventory is a content identity, not an invented git commit.
