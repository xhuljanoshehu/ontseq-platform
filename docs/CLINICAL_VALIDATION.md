# Analytical and clinical validation plan

## Current classification

The software is a research prototype. No output is validated for diagnosis, prognosis,
treatment selection or patient release.

## Validation units

Validation must be stratified; an aggregate accuracy number is insufficient.

| Dimension | Required strata |
| --- | --- |
| Assay | lcWGS, adaptive sampling, future RNA or methylation assays |
| Reference | GRCh37 and GRCh38 independently |
| CNV | whole chromosome, arm-level, focal deletion, focal duplication/amplification |
| SV | deletion, duplication, inversion, insertion, translocation/BND |
| Fusion | each clinically intended fusion class and breakpoint coverage pattern |
| Sample | fresh, archived, purity/cellularity range, coverage range, DNA quality range |

## Required studies

1. **Reference truth set:** orthogonally characterized positives and negatives.
2. **Accuracy:** sensitivity, specificity, precision, false-positive burden and no-call rate.
3. **Limit of detection:** cellularity/VAF, coverage and supporting-read thresholds per event
   class.
4. **Precision:** within-run, between-run, operator, instrument, flow cell and reagent lot.
5. **Robustness:** DNA input, quality, read N50, alignment, reference and target-BED changes.
6. **Interference:** repeats, segmental duplications, mapping ambiguity, clonal complexity and
   adaptive-sampling coverage gradients.
7. **Software regression:** locked fixtures for every validated event and every known failure.
8. **ISCN verification:** authorized reference examples, parser/renderer round trip, ambiguity
   handling and expert sign-off.

## Orthogonal comparators

- conventional karyotyping for chromosome- and arm-level events;
- FISH and/or chromosomal microarray for selected CNVs;
- RT-PCR, RNA sequencing or validated FISH for fusion confirmation;
- expert IGV review for breakpoint evidence, never as the sole truth standard.

## Release gates

A release candidate must provide:

- a signed validation report tied to the exact code release and reference bundle;
- locked containers/environments and reference checksums;
- a completed risk assessment and intended-use statement;
- acceptance criteria for every QC and reportability gate;
- two-person review for changes that can alter biological output;
- rollback, incident-response and post-release monitoring procedures.

`PASS` from a caller or QC module is not equivalent to clinical validity.
