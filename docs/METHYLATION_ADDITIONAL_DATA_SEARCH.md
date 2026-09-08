# Additional public ONT methylation data search

**Acquisition update:** After explicit user authorization the two selected GSE218679
objects were downloaded and structurally checked locally. Verified public-object hashes
are in `configs/methylation/acquired_public_sources.v1.json`; current status is in
`METHYLATION_INTAKE_STATUS.md`. The metadata-only catalogue and review below are the
preserved pre-acquisition snapshot; remaining assay provenance gaps are not filled by
the successful download.

Research use only. Metadata catalogue `additional_data_candidates.v1`, assessed
2026-09-05. Companion: `configs/methylation/additional_data_candidates.v1.json`.
No genomic or methylation measurement file was downloaded, including no BAM header
range, TSV body, FAST5/POD5, SRA object or genomic index. Public metadata APIs,
directory listings, articles and HTTP HEAD were used. No local genomic SHA-256
has been computed; every such catalogue field remains `null`.

## Immediate small intake proposal: GSE218679 native iPSC pair

The most practical next acquisition is **two full Nanopolish TSV.gz files totalling
208,521 bytes (0.199 MiB)**. This is a narrow FMR1-source feasibility pilot, separate
from GSE210260 and the GIAB sources. It is ready for an intake decision, not admitted
to recovery or donor-independent validation.

| Role | Public accession | Biological source | Exact compressed bytes |
| --- | --- | --- | ---: |
| A | [GSM6754749](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6754749) | Native genomic DNA from an unengineered FXS iPSC parent line | 87,638 |
| B | [GSM6754755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6754755) | Native genomic DNA from a normal-length FMR1 iPSC line | 120,883 |

Only these two objects are proposed:

- [GSM6754749_135_Nanopolish.tsv.gz](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6754nnn/GSM6754749/suppl/GSM6754749_135_Nanopolish.tsv.gz)
- [GSM6754755_barcode01_Nanopolish.tsv.gz](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6754nnn/GSM6754755/suppl/GSM6754755_barcode01_Nanopolish.tsv.gz)

The [GEO file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218679/suppl/filelist.txt)
and HTTP HEAD independently agree on those sizes. Both HEAD responses report
2022-11-23 19:59:00 UTC as Last-Modified. The series was public on 2023-12-21;
the two GSM metadata records were last updated 2025-12-19. The series metadata
was updated 2026-04-13. Public identifiers are GSE218679 / PRJNA904800 / SRP409538.
These timestamps describe public metadata/object revisions, not immutable checksums.

Both GSM records explicitly describe iPSC genomic DNA, MinION, the same extraction,
Cas9 FMR1 targeting, 5–12 kb BluePippin size selection and EXP-NBD104 native barcoding.
Their library descriptions contain no whole-genome amplification, PCR amplification,
bisulfite treatment or exogenous CpG methylation step. The Taq step is dA-tailing;
it is not a PCR thermocycling protocol. The original primary
[Cell study methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10794044/) identify
SQK-LSK109 ligation sequencing. This is assay-level evidence, not a per-file kit receipt.
Select the disease parent and normal line, not CRISPR-engineered cutback clones.

| Provenance field | Source A and B | Evidence limit |
| --- | --- | --- |
| Canonical caller | Guppy 6.2.1 | Full model/configuration absent in reviewed GSM metadata |
| Aligner | minimap2 2.22-r1101 | Same declared software |
| Reference | hg38 / GRCh38 | Exact FASTA release, contig dictionary and SHA-256 absent |
| Methylation caller | Nanopolish 0.13.2 | Exact model artifact/hash absent |
| Target modification | CpG 5mC | No separate 5hmC channel; STRique results are a different evidence layer |
| Instrument | MinION | Same GSM instrument; exact pore/flowcell chemistry absent |
| Library | Native Cas9 targeting; EXP-NBD104; study SQK-LSK109 | Do not infer R9.4.1 solely from STRique model names |
| Calibration/test independence | One file per distinct iPSC source | Read holdout possible after admission; no independent preparation/donor holdout |

The [ENA run metadata](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA904800&result=read_run&fields=run_accession,experiment_accession,read_count,instrument_model,experiment_title&format=json)
reports 8,603 archived reads for A (SRR22398697 / SRX18368393) and 1,201 for B
(SRR22398703 / SRX18368399). These are **not** usable Nanopolish read-group counts.
The TSV may contain only a subset. Do not fill read-depth budgets from these counts.

The deposited caller window is described as hg38 chrX:147902117–147960927. The
published analysis uses 19 CpGs in the 500 bp FMR1 promoter at
chrX:147911419–147911919. It is unknown whether each export contains the complete
caller window or only the selected promoter. A wider calling window does not establish
20 eligible markers. Preserve the repository's minimum marker requirement: a
`NO_CALL` would be an informative outcome. Nearby CpGs at a single locus remain
correlated even if their numerical count exceeds the floor.

The GSM descriptions mention reverse-strand and sequence-content filtering around
the repeat. The exact point at which these filters affected the deposited TSVs is
unclear. They also describe a likelihood cutoff of 0.1. That cutoff must not replace
the repository's policy, and the raw score units/header need verification.

### Reviewable acquisition and structural intake sequence

The proposed destination is
`C:/Users/sxhul/Documents/ChatGPT/ONTSeq/work/methylation_validation/gse218679/intake_v1`;
`work/` is already ignored by Git. No acquisition has been performed. Genomic-data
intake is subject to the owner/institutional approval boundary in AGENTS.md rule 7.

1. Authorize exactly the two complete compressed TSV files listed above. No whole-series
   archive, genomic signal, SRA download or unrelated sample is required.
2. Check the destination remains ignored and refuse overwrite. Record retrieval time,
   URL, expected byte size and complete compressed-file SHA-256 locally.
3. Perform a bounded gzip/header/schema check, record decoded size, count distinct read
   groups and coordinate/strand representations, and check cross-source read-group
   overlap. Do not calculate source methylation differences, marker rankings or recovery.
4. Resolve exact caller/model, reference and chemistry metadata; record anything still
   absent as an admission gap. A copied default is not provenance. Preserve the original
   TSVs and their parent hashes.
5. The candidate matrix `fmr1_targeted_holdout.v1.json` already fixes budgets 20/50/100,
   fractions, repetitions, seeds, marker policy and acceptance criteria before acquisition.
   After structural inventory and provenance admission, bind actual input hashes and the
   whole-read split policy in an immutable registration. Feasibility remains unverified
   until usable group counts are known. Preserve v1 and record NO_CALL for unmet gates;
   do not adjust it after outcomes. No biological outcome was inspected to choose files.
6. Run once against held-out read groups under the frozen rules. Report retained-read-group
   recovery and its conditional interval. This pair cannot establish tumor purity, a
   physical DNA fraction, genome-wide accuracy, independent-donor generalisation or a
   10–20% limit of detection.

## Other human candidates

### ONT RRMS 2022.07: COLO829 versus COLO829BL

The [official release](https://epi2me.nanoporetech.com/rrms2022.07/) provides per-read
modification BAMs for five flowcell replicates of each cell line. Both use R9.4.1 and
the same approximately 310 Mb adaptive-sampling target design. Public paths identify
Bonito/Remora processing; the exact versions and full model remain unspecified.
The hg38 target config is public, but the exact aligned reference must be verified.

The complete [S3 modBAM inventory](https://ont-open-data.s3.eu-west-1.amazonaws.com/?list-type=2&prefix=rrms_2022.07%2Fbasecalling%2Fbonito_remora%2F)
contains ten BAMs. Every full URL, size, Last-Modified and ETag is in the JSON catalogue.

| Intake | Objects | Bytes / GiB |
| --- | --- | ---: |
| Minimum two-source pair | COLO829_3 + COLO829BL_2 | 44,670,539,944 / 41.64 |
| Additional held-out flowcells | Also COLO829_4 + COLO829BL_4 | 95,893,968,014 total / 89.34 |
| All ten BAMs | Five per source | 257,380,503,998 / 239.70 |

The minimum pair is selected by smallest complete-file size before observing biology.
A four-file proposal can hold out different flowcells, but independent cultures,
extractions and donors have not been established. The archived AS rejection/short-read
policy must be documented; the coefficient is conditional on retained read groups.
Matched RRBS can assess CpG concordance, not the truth of a physical mixture.

The [later ONT poster](https://nanoporetech.com/api/assets/f/196663/735a4dd4c9/j2127_a1-poster_rrms-v3_digital.pdf)
also describes HCC1395 and 5hmC via a Megalodon/Remora workflow. Those details must not
be transferred to the 2022 Bonito files: HCC1395 BAMs were absent from the complete
reviewed prefix. No download paths were invented.

### K562 GSE173688 plus native H1ESc PRJNA876781

[K562 GSM5277007](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5277007&targ=self&form=text&view=full)
has a public [Nanopolish file listing](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5277nnn/GSM5277007/suppl/).
The Nanopolish TSV.gz is exactly 380,186,764 bytes (HTTP HEAD). Metadata specifies
GridION, R9.4.1 FLO-MIN106, SQK-RAD004, Guppy 4.2.2 HAC and hg38. The associated
[benchmark study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8524990/) reports Nanopolish 0.13.2.

The [Rockfish primary study](https://www.nature.com/articles/s41467-024-49847-0)
introduces native H1 embryonic stem-cell data in PRJNA876781. Six GridION native runs
are verified in ENA metadata and listed in the catalogue. A compatible ready per-read
export was not located. The study's Guppy 5.0.14/CHM13 processing differs from the
published K562 export; both must be processed consistently before pooling. Native H1
runs must not be confused with REPLI-g or M.SssI-treated controls. Multiple runs under
one BioSample do not establish independent preparations.

### nanoNOMe breast cell lines: GSE155791 / Zenodo 3969567

The [primary GSM metadata](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4712210)
documents Guppy 3.0.3 HAC, NGMLR 0.2.8, hg38 without alt contigs and Nanopolish 0.11.1.
[Zenodo's record metadata](https://zenodo.org/api/records/3969567) lists
`MCF10A_MCF7_MDAMB231_hetSVs_1kb_cpg_calls.txt.gz` at 87,419,865 bytes for three breast
cell lines; supplier MD5 is `a503979bdddef1057b8453c6bd739cfb` (not locally verified).

This is a custom per-read analysis table, not confirmed as the Nanopolish raw contract.
It uses exogenous GpC labeling and selected heterozygous-SV regions. Native CpG/GCG
disambiguation, preserved source/read identities, instrument comparability and selection
bias need a separate adapter/policy. It is not a simpler substitute for the iPSC intake.

### COLO829 / COLO829BL WGS release 2023.04

The [official Kit V14 release](https://epi2me.nanoporetech.com/colo-2023.05/) supplies a
matched PromethION 24/R10.4.1/DNeasy/GRCh38 pair with CRAMs and POD5. Two standard
tumor flowcells and one normal flowcell are described. Separate Monarch/ULK114 libraries
must be kept in different assay strata. The selected standard CRAMs are
54,275,227,484 and 66,233,833,667 bytes in the public S3 metadata.

This variant-focused release does not explicitly establish preserved MM/ML or a
methylation model. Exact caller/model/tag provenance must be resolved, or both samples
recalled consistently from approved signal data. It reuses the RRMS cell lines and is
therefore new technical data, not an independent source-class donor panel.

## Out-of-scope mouse alternatives and screened dead ends

Two real mouse alternatives are fully catalogued but are not the current human intake:

- [PRJNA1172918 / ONT methylation benchmark](https://registry.opendata.aws/ont_basemod_data/):
  public Dorado modBAMs for brain and ESC, R10.4.1/5 kHz/mm39. A matched v4r1 CpG pair
  occupies about 223 GiB; v5r1 about 186 GiB. Several model recalls of one preparation
  are not biological replicates. The [release documentation](https://github.com/SowpatiLab/ont-basemod-benchmark-data/blob/main/documentation.md)
  is V0.2 dated 2026-03-28; exact objects are dated December 2025.
- [GSE248540](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248540):
  murine Th1/Th2 cell preparations with three replicate labels each,
  MinION R9.4.1/LSK109, Dorado 0.4.2, mm10 and modkit 0.2.2.
  Public supplements are BED pileups; compatible per-read modBAM availability and
  biological origin of the replicate labels were not confirmed.

[GSE173687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173687&targ=self&form=text&view=full),
sometimes cited for HL60, was explicitly deleted by GEO on 2021-09-20. It is not an
available public pair. Tiny RNA poly(A) TSVs and binary tumor-classifier BED files found
during the targeted search were excluded because they do not represent native DNA
per-read CpG likelihoods.
