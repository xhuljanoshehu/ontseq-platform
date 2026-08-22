# Comparing against the historical Hannover pipeline

The predecessor of this repository is a Snakemake workflow at
`git.finf.uni-hannover.de/aal/ONTseq`. It runs QDNAseq + ACE, mosdepth, Sniffles2, Cramino,
SnpEff/SvAnna/ANNOVAR and an ISCN renderer inside Apptainer containers, and writes a set of
loose files per sample.

The two are **not** the same kind of software and one is not a port of the other. The old
workflow scans a directory of BAMs and produces outputs; this one builds one auditable run
envelope per sample with preflight, reference locking, checksums, resume and an explicit
`NOT_RUN` / `NO_CALL` distinction. What they can be compared on is results, and that
comparison is worth doing: same input, same reference, same panel, two implementations.

## Output mapping

| Legacy file | Contents | Counterpart here |
| --- | --- | --- |
| `CN.csv` | `Chromosome,Copies,Ploidy,CNA` | CNV events in `normalized/{sample}.result.json` |
| `cellularity.txt` | single ACE cellularity value | ACE fit in the CNV stage provenance |
| `dels_dups.csv` | `chromosome,name,event,frac_abr` | normalized copy-number events |
| `bins.bed`, `segments.bed` | UCSC `track` format, per-bin and per-segment | QDNAseq bin and segment artifacts |
| `copyplot.png`, `errorplot.png` | ACE plots | CNV stage plot artifacts |
| `cramino_stats.json` | `file_info`, `alignment_stats`, `read_stats`, `identity_stats`, `karyotype_stats`, `histograms` | `qc/cramino.json` |
| `mosdepth/*.regions.bed.gz` | per-region depth | `qc/target-coverage.json` |
| `sv.vcf`, `sv.snf` | Sniffles 2.4 output | `evidence/sv/{sample}.sniffles.vcf` and `.json` |
| `fusions.tsv` | 13 columns including HGVS | not implemented; fusion interpretation is a separate validation domain |
| `ISCN.txt` | rendered karyotype string | proposal-only renderer, explicitly unvalidated |
| `report.html` | ezCharts report | `reports/{sample}.report.html` |

## Differences that must be held constant before anything is concluded

* **Sniffles version.** The legacy VCFs were written by **Sniffles 2.4**; this repository
  pins **2.8.0**. Run `--components configs/components/legacy_sniffles_2.4.yaml` to hold the
  version constant, or accept that a difference in calls may be the version rather than the
  pipeline. See [`COMPONENT_SELECTION.md`](COMPONENT_SELECTION.md).
* **Aligner.** Some legacy outputs came from ngmlr, this repository uses minimap2.
* **Region restriction.** The legacy SV run was restricted to the buffered panel with
  `--regions`. A whole-genome run here will find breakpoints the legacy run could not.
* **Reference.** The legacy runs used `data/reference/hg38.fa`. Which distribution that is
  must be established and locked before any comparison is meaningful.

## What the legacy outputs are and are not

They are a **regression reference**: does this pipeline find the same breakpoints and the
same copy-number segments on the same input? A disagreement is a question, not a defect on
either side.

They are **not ground truth**. The fusion tables are unfiltered caller output — the highest
read-support entries are MUC6–MUC6 self-"fusions" and intergenic pseudogene artifacts. Using
them as a truth set would encode caller artifacts as expectations.

## MOLM13 as a positive control

The laboratory records include repeated runs of **MOLM13**, an AML cell line, alongside
patient material. It is the right first positive control for this pipeline:

* its karyotype is published, so an expectation exists that no patient data is needed to
  state;
* it carries a KMT2A–MLLT3 fusion and an FLT3 internal tandem duplication, and both genes
  are on the panel;
* it can be re-sequenced, so a failed comparison can be repeated;
* it raises no governance question, which means the work can start now.

What still has to be pre-specified before a MOLM13 run counts as a control rather than an
anecdote: which events are expected, at which coverage, with which caller versions, and what
result would count as a failure.

## Obtaining the legacy material

The code is at `git.finf.uni-hannover.de/aal/ONTseq` and is an access question rather than a
search question. Worth requesting, in order of usefulness:

1. `data/reference/fusion_panel_with_buffer.bed` — coordinates, no patient data, and it
   settles the open questions in [`PANEL_PROVENANCE.md`](PANEL_PROVENANCE.md);
2. the exact reference distribution, so it can be locked;
3. `config.yaml` and the rule definitions, for parameter comparison;
4. released outputs for samples this repository can also process.

Ownership and licensing of the thesis code must be clarified before any of it is copied
into this repository. Nothing in this document depends on that having happened.
