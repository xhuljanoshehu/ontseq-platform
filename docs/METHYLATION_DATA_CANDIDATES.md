# Public ONT methylation source candidates

Research use only. Catalogue version `public_data_candidates.v1`, assessed 2026-09-05.
This is a metadata review, not a completed dilution experiment or biological validation.
No BAM, CRAM, POD5, FAST5, FASTQ, index, read-name list or genomic byte range was acquired.
Only public articles, README documents, directory/S3 listings, checksum text and HTTP HEAD
metadata were read. The machine-readable companion is
`configs/methylation/public_data_candidates.v1.json`.

No candidate below currently satisfies independent biological calibration/test admission.
Every genomic object's local SHA-256 is `null` because it has not been acquired. Published
MD5 values are recorded as supplier metadata; they are neither locally verified nor SHA-256.
S3 multipart ETags are object metadata and must not be relabelled as file hashes.

## Preferred next intake: HG008 normal pancreas and normal duodenum

**Purpose:** a new, biologically distinct tissue-source feasibility experiment, independent
of GSE210260. Source A is HG008-N-P and source B is HG008-N-D. Both tissues derive from one
publicly consented reference donor. This supplies tissue contrast, but no independent donors,
purified cell types, tumor-fraction truth or physical dilution truth. Public dissemination is
described by [NIST Cancer GIAB](https://www.nist.gov/programs-projects/cancer-genome-bottle).

The [Northeastern ONT-std README](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern_ONT-std_20240422/README_NE_ONT-std.md)
documents matched PromethION R10.4.1/LSK114 sequencing, Dorado 0.5.3 and
`dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mC_5hmC@v1`. Both 5mC and 5hmC channels exist;
their probabilities must remain distinguishable. It specifies minimap2 2.26-r1175,
samtools 1.21 and GRCh38-GIABv3 alignments. The February 2025 alignment revision retains
modification tags; older alignments must not be substituted. Two DNA extractions per tissue
were pooled before sequencing. Three pancreatic and four duodenal flow cells therefore do
not supply independent extraction or donor replicates. The first pancreatic flow cell had
lower quality and a different loading amount. Raw POD5 archives are not on this FTP site.
These are supplier statements, not locally checked BAM-header observations.

### Concrete object choices and storage

Exact byte sizes were verified by HTTP HEAD on 2026-09-05. Full object URLs, modification
dates and supplier MD5 values are in the JSON catalogue. All listed objects are in the
[official NCBI data directory](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern_ONT-std_20240422/).

| Option | Source A bytes | Source B bytes | Total | Intended use |
| --- | ---: | ---: | ---: | --- |
| Individual unaligned flow-cell BAMs: pancreas `_2`, duodenum `_3` | 49,607,859,349 | 77,856,391,489 | 127,464,250,838 bytes; 118.71 GiB | Preferred initial feasibility intake; common local alignment needed |
| Full GRCh38-GIABv3 BAMs: pancreas 41x, duodenum 94x | 185,972,227,157 | 450,755,274,622 | 636,727,501,779 bytes; 593.00 GiB | All supplier reads; substantially more storage needed |
| GRCh38-GIABv3 BAMs filtered at 25 kb: pancreas 20x, duodenum 13x | 78,641,825,680 | 58,932,876,550 | 137,574,702,230 bytes; 128.13 GiB | Explicit length-restricted sensitivity experiment only |

The full pair exceeds the approximately 263 GiB currently free on C: reported during this
task. The smaller uBAM pair may fit as source files, but does not by itself establish enough
space for aligned copies, sorting scratch, per-read extracts, reference files and results.
Approve a storage budget before acquisition; avoid a full TSV export of every CpG.

Pancreas `_2` avoids the explicitly documented first-flow-cell loading issue; duodenum `_3`
is the smallest listed duodenal uBAM. This is a storage/QC based candidate choice, made before
observing estimator recovery. It is not a claim that these files recover proportions.
An initial seeded, whole-read-group subsample before common local alignment can limit the
derived-data volume. Record the input hash, sampling algorithm, read budget and seed before
sampling, and retain the all-length source denominator. This step requires implementation
and verification; it has not been run here.

Do not use the 25 kb files merely because they are smaller. Their estimand is conditional on
length selection, and a source-dependent read-length distribution changes which molecules
enter the experiment. Region-restricted remote extraction would likewise change inclusion;
it is not a free substitute for an all-length, genome-wide dilution experiment.

## Alternatives and admission gaps

| Candidate | Comparable aspects | Decision and unresolved admission |
| --- | --- | --- |
| ONT `giab_2023.05`, HG001 versus HG002 | Same release, PromethION R10.4.1, Kit14, 5 kHz; two physical flow cells per line; Dorado SUP and HAC analyses; 5mC/5hmC workflow | Useful independent technical control from distinct reference cell lines; marker contrast is unproven. Preserve physical flow-cell grouping across logical restarts. Exact caller/model, reference digest and modification tags still need verification. No independent donors within a source class. |
| HG008-T versus HG008-N-P or HG008-N-D | Matched cancer/normal reference material and R10.4.1 options | Do not mix published products as a compatible pair: ONT-std-1 normal uses Dorado 0.5.3/v4.3; ONT-std-2 tumor uses 0.3.4/v4.2 with different preparation/shearing. Recalling to one model does not erase preparation or cell-line/tissue confounding. |

The [official ONT release](https://epi2me.nanoporetech.com/giab-2023.05/)
describes its sequencing and logical run restarts. S3 ListObjectsV2 confirms per-run SUP
CRAMs. A small initial pair (one HG001 logical run plus one HG002 flow cell) is
90,410,403,222 bytes, but supplies read-split feasibility only. Four files covering both
HG001 physical flow cells total 147,144,169,971 bytes; two HG002 files total
181,157,379,736 bytes. Thus a complete physical-flow-cell comparison totals
328,301,549,707 bytes before derivatives. CRAM-to-BAM conversion needs a verified matching
reference and tag-preserving validation. The [ONT registry](https://registry.opendata.aws/ont-open-data/)
describes mutable releases; record exact objects and current applicable terms rather than
treating a bucket prefix as immutable.

The [HG008 primary data descriptor](https://www.nature.com/articles/s41597-025-05438-2)
documents the tumor/normal caller and preparation mismatches above and states that methylation
was not assessed in its data validation. Somatic variant benchmarks from this material are
not methylation-mixture truth. Later tumor passages, ultralong libraries and model versions
must remain separate strata. None of these alternatives is an AML validation cohort.

## Intake checklist before any genomic acquisition or benchmark execution

1. Record owner/institutional authorization and a local destination outside versioned data.
   The repository remains metadata/code/synthetic fixtures only; no external upload, cloud
   analysis, public deployment or sharing is part of this catalogue.
2. Choose feasibility versus biological-validation intent explicitly. To admit the latter,
   secure independently collected source-A and source-B samples for calibration and testing,
   with donor/extraction/library assignments known before marker selection. Reusing a donor,
   pooled extraction, flow cell or read pool is not an independent biological replicate.
3. Freeze exact file URLs, release/revision, acquisition list and space estimate. Recheck HEAD
   size/date and supplier checksum document. Record transport date and compute local MD5 plus
   SHA-256 over each complete acquired file; mismatches reject intake. Obtain and hash a
   reference FASTA matching the full dictionary, not merely the label `GRCh38`.
4. Inspect BAM/CRAM integrity and headers locally: `@RG`/`@PG`, caller and full model identifiers,
   flow cell/library identity, genome dictionary, sequence orientation and intact `MM`/`ML`.
   Freeze a local modkit version and extraction policy. Missing modification tags, incompatible
   builds/models, unresolved chemistry or missing per-read identity reject quantitative use.
5. Define CpG coordinates, primary-alignment filters, read/group identity, canonical/5mC/5hmC
   handling and missing-call semantics before extracting any markers. A 5mC+5hmC model is not
   a 5mC-only model; summed modifications must be a separate explicit target. BedMethyl
   aggregates alone cannot preserve disjoint read groups for dilution.
6. Lock the validation matrix, accuracy/coverage/no-call limits and seeds before test results.
   Lock markers using calibration only. Record every scheduled cell including failed intake,
   `NO_CALL` and failed recovery; do not tune policy on held-out outcomes and call it validation.
7. Keep operational completion separate from recovery and biological admission. A numeric
   estimate or software `PASS` on read-resampled feasibility cannot promote this dataset to
   donor-generalizable source fraction, DNA mass, cell fraction, tumor purity or a 10-20% LoD.

## Validation impact

This catalogue changes no estimator or reportability threshold. It nominates future evidence
and deliberately leaves biological admission unavailable. No bias, MAE, RMSE, interval coverage
or `NO_CALL` rate has been measured for these candidates. Subsequent numerical recovery can
only be reported with its actual independence level and locked experiment provenance.
