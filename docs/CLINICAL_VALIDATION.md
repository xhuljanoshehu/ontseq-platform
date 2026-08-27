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
| Input path | aligned BAM, uBAM and POD5-to-report independently |
| Reportability | positive, negative, equivocal and no-call outcomes |

## Required studies

1. **Reference truth sets:** public HG002/HG008 material for transferable technical behavior,
   followed by orthogonally characterized AML positives and negatives for intended use.
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
9. **Coverage and purity:** predefined whole-genome, per-target and breakpoint coverage gates,
   plus tumor/blast-fraction dilution studies with explicit no-call thresholds.

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

Acceptance criteria, truth definitions, dataset versions and exclusions must be locked before
caller comparison. Public cancer benchmarks are useful engineering controls, but do not replace
validation in the intended AML specimens and workflow.

`PASS` from a caller or QC module is not equivalent to clinical validity.

## SV evidence-layer validation impact

The two-caller consensus, annotations, Adaptive Sampling observability, AML pattern lookup and
technical score alter which candidates a reviewer sees first. They do not alter reportability:
all automatically produced SV candidates remain `reportable=false` and no confirmed fusion is
emitted.

Promotion of any current technical default requires, at minimum:

- frozen Sniffles2, cuteSV, consensus and score policies;
- build-specific gene, cytoband and artifact-context resource checksums;
- per-caller and consensus results on locked HG002/HG008 technical truth sets;
- AML positive and negative specimens characterized by karyotype, FISH and/or PCR/RNA methods;
- breakpoint-error, false-positive burden and no-call results by SV class, coverage, VAF or
  tumor/blast fraction and assay mode;
- explicit verification of events outside or partly inside the Adaptive Sampling ROI;
- validation of reviewer workload and traceability from each cluster to every source caller record.

Until those studies and predefined acceptance criteria pass, `high` means high technical review
priority only. It must not be mapped to `analytically_validated`, `reportable`, a negative result or
a clinical action.
