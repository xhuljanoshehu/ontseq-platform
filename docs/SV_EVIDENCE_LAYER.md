# SV evidence prioritization

## Purpose

The SV evidence layer reduces the review burden from a large raw caller output to a compact,
deterministic review queue. It is a **technical triage layer**, not analytical or clinical
validation.

Every normalized SV remains available in the full result. The prioritizer never converts an
SV into a clinical assertion and never sets `reportable=true`.

## Current evidence used

The first implementation uses fields already normalized from the caller VCF:

- number of independent callers represented in the event evidence;
- supporting reads;
- variant allele fraction when available;
- caller precision flag;
- caller FILTER state;
- extra review priority for balanced rearrangements with meaningful read support.

These features generate `high`, `moderate`, or `low` **technical confidence**. High and
moderate events are surfaced in the HTML `SV review queue`; low events remain in the complete
event table.

The score is intentionally transparent and deterministic. Its current thresholds are technical
defaults only and must not be described as sensitivity/specificity validated cut-offs.

## Why reportability stays false

Confidence and reportability answer different questions:

- **technical confidence**: is this candidate internally well supported enough to prioritize?
- **reportability**: has the assay and decision rule been validated against independent truth
  data for the intended use?

Until an assay-specific benchmark establishes acceptance criteria, `reportable` stays false even
for technically high-confidence events.

## Validation path

The intended calibration path is:

1. synthetic positive/negative fixtures for parser and event-normalization correctness;
2. public structural-variant benchmark samples where applicable;
3. local AML specimens with independent cytogenetic or molecular truth, such as karyotype,
   FISH, PCR/RT-PCR, or another validated method;
4. reproducibility and dilution/coverage series to characterize limits of detection;
5. locked acceptance thresholds derived from the resulting sensitivity, precision, false-positive
   rate, breakpoint tolerance, and no-call behavior.

## Next implementation gates

The review queue is designed so stronger evidence can be added without changing the report
contract. The next useful extensions are:

1. run and normalize cuteSV in addition to Sniffles2;
2. merge equivalent breakpoints across callers using explicit distance/orientation rules;
3. attach cytoband and gene annotations to both breakpoints;
4. add repeat/blacklist and mappability context;
5. quantify Adaptive Sampling target observability rather than inferring it from whole-genome
   mean coverage;
6. prioritize recurrent AML rearrangements using a version-locked somatic knowledge resource;
7. calibrate the technical score on orthogonally characterized AML samples before any rule is
   allowed to influence `reportable`.

## Design principle

A clinician should review a small set of interpretable, evidence-rich candidates rather than
hundreds or thousands of raw SV calls. Automation should reduce review volume while retaining
traceability to every original candidate and every evidence field.
