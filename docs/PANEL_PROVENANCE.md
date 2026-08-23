# Where the Adaptive Sampling panel came from

`configs/panels/aml_fusion_adaptive_sampling.grch38.buffered.bed` is the first target design
in this repository that is not synthetic. This document records exactly what it is, what
evidence supports that, and what has to be resolved before it may be treated as locked.

## What it is

111 uniquely labelled target intervals covering 17,028,377 bases on GRCh38. The labels
come from the laboratory workbooks; they are not evidence that all 111 labels have been
independently confirmed against an authoritative gene annotation. This is an AML fusion and
karyotyping panel. It is a **buffered selection panel**: each interval extends roughly
10 kb beyond the gene it names.

## Where it came from

Two laboratory workbooks, neither of which is in this repository and neither of which may
be: they carry sample identifiers.

* an Adaptive Sampling coverage and fusion workbook, which lists the targets together with
  per-run coverage;
* an ONT experiments workbook, which lists the same targets with mean coverage across
  further runs.

`scripts/build_adaptive_sampling_panel.py` reads both from paths the operator supplies and
writes only coordinates and aggregated numbers. The lock file records the SHA256 of each
source workbook, so a later reviewer can prove which revision produced the committed BED.

## Why we believe the build and the role

Four independent observations agree, which is the only reason this file exists at all
rather than a request to the laboratory:

1. **Both workbooks describe an identical set of 111 intervals.** They were produced by
   different people at different times.
2. **The chromosome set matches the legacy pipeline exactly.** A released Sniffles VCF from
   the historical Hannover pipeline carries contigs for chr1–chr19, chr21, chr22 and chrX —
   no chr20, no chrY. The panel contains no chr20 or chrY gene. That VCF was produced with
   `--regions data/reference/fusion_panel_with_buffer.bed`, so the interval set in these
   workbooks is the region file that run was restricted to.
3. **The contig lengths in that VCF are GRCh38**, not GRCh37 (chr1 = 248,956,422).
4. **Interval ends sit exactly 10,000 bp beyond the Ensembl GRCh38 gene end** for MYC,
   ABL1, KMT2A and PRDM16. Starts vary by a few hundred bases, consistent with a different
   transcript source for the gene start. Together with the legacy file name, this is what
   makes the design buffered rather than an analysis ROI.

Every interval was additionally checked against the GRCh38 primary assembly lengths, for
duplicates, for overlaps and for non-positive length. None was found.

## What is not established

* **That this file is byte-identical to the panel the sequencer selected on.** It is
  reconstructed from coverage tables, not copied from `fusion_panel_with_buffer.bed`.
  Obtaining that file remains the correct next step; this BED is what allows work to
  continue in the meantime.
* **The coordinate convention.** The source is a region-string column, so whether the start
  is 0-based half-open or 1-based inclusive is unverified. The difference is one base per
  interval — immaterial to a coverage mean, material to a locked contract.
* **Any threshold.** The panel says where to look. It says nothing about how much depth is
  enough.

## The open question that blocks promotion

One row is labelled **IGH** and placed at `chr5:143,396,959-143,417,420`. IGH is on chr14q32
in GRCh38; that chr5 interval lies inside the NR3C1 region. Both workbooks carry the same
value, so the error is upstream of this repository — it is not a transcription mistake made
here. The interval is kept in the BED, renamed `IGH_REVIEW_REQUIRED`, and named in the lock
file's `open_questions`. Nothing is silently dropped and nothing is silently accepted.

Until somebody in the laboratory says which of the two is wrong, the label or the interval,
this panel stays `status: derived_unconfirmed`.

## Coverage expectations

`configs/qc/target_coverage_expectations.grch38.tsv` gives the observed per-target mean
depth across nine historical runs, as minimum, median and maximum. Run labels are reduced
to a count so no sample can be identified from the file.

It is a sanity reference: a new run whose per-target depths fall far outside these ranges is
worth investigating before its results are believed. It is **not** an adequacy gate, a
reportability threshold or a no-call definition, and it must not become one without
pre-specified validation.
