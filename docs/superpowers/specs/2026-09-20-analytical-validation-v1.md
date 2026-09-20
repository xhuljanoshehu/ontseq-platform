# ONTSeq Analytical Validation Program v1

Date: 2026-09-20
Status: approved for staged implementation
Branch: `feat/analytical-validation-v1`
Base: `caca29feab87d6f58dc7d172e42a51a70f036fa0`
Scope: Research Use Only; Human Review Required

## 1. Purpose

Create the first prospectively registered ONTSeq CNV validation study framework spanning the four primary analytical axes:

- coverage;
- tumour/blast fraction;
- CNV size;
- caller.

The program must reuse the existing CNV validation contracts, registration lock, metric engine, evidence-retention layer and dilution machinery rather than creating parallel validation semantics.

This v1 repository implementation defines and locks the study design and generated synthetic truth declarations. It does not claim analytical or clinical validity.

## 2. Scientific boundary

Repository v1 uses synthetic truth/specimen declarations only.

Real analytical-validation evidence requires authorised, orthogonally characterised material in the appropriate execution environment. No patient BAM/FASTQ/POD5/VCF payloads or restricted truth resources are committed to Git.

Synthetic registration is a software/study-design qualification, not an assay-performance result.

## 3. Primary factor grid

### Coverage levels

Planned nominal coverage levels:

- 0.5x
- 1x
- 2x
- 5x
- 10x
- 20x

These are planned study levels, not universal clinical thresholds.

### Tumour-fraction levels

Planned fractions:

- 1.00
- 0.50
- 0.30
- 0.20
- 0.10
- 0.05

A pure-normal 0.00 control remains a negative control and is not treated as a positive dilution level.

### Subchromosomal CNV sizes

Planned event sizes:

- 100 kb
- 250 kb
- 500 kb
- 1 Mb
- 5 Mb
- 10 Mb

Whole-chromosome gain/loss remain distinct event classes, not another numeric size bin.

### Caller families

Primary CNV comparison families are registered only where the data regime is compatible:

- QDNAseq+ACE
- Spectre depth-only
- SAVANA CNA where the registered tumour-analysis inputs exist
- Wakhan phased CNA where the registered phased/breakpoint-assisted inputs exist
- ichorCNA only in its dedicated cfDNA/ULP-WGS study, not silently in marrow/lcWGS

The initial genome-wide lcWGS study does not force incompatible callers into the same denominator.

## 4. Event classes

v1 synthetic truth templates cover:

- subchromosomal deletion;
- subchromosomal duplication/gain;
- whole-chromosome loss;
- whole-chromosome gain.

Subchromosomal sizes are crossed with deletion and duplication only.

## 5. Repeats

Default software-study materialisation uses three deterministic replicates per positive design cell.

Replicate IDs and fingerprints are derived deterministically from the study seed and cell identity.

This exercises reproducibility contracts but does not make synthetic replicates statistically independent biological specimens.

## 6. Truth-series architecture

Add `src/ontseq_platform/cnv_validation_v1.py`.

It contains:

- typed factor/design contract;
- deterministic design-cell expansion;
- synthetic truth-event templates;
- synthetic cohort materialisation;
- validation-matrix builder;
- registration builder using existing `preregister_cnv_validation`;
- summary counts for planned cells/events.

No outcome data enter the builder.

## 7. Registration semantics

The v1 builder must produce existing canonical objects:

- `CnvValidationMatrix`
- `CnvValidationCohort`
- `CnvValidationRegistration`

Registration must remain content-addressed and `outcomes_unseen=True`.

Changing any of the following must change the registration lock:

- factor levels;
- caller lock;
- event template;
- matching thresholds;
- replicate count;
- reference identity;
- truth/mask identity;
- acceptance questions;
- code identity.

## 8. Acceptance questions

Repository defaults are prospective engineering candidates, not clinical cutoffs.

v1 registers questions for:

- sensitivity;
- precision;
- false-positive burden;
- no-call rate;
- technical-failure rate;
- copy-number error where quantitative truth is present.

Specificity is omitted from the default v1 matrix unless an explicit negative universe is materialised and registered.

No acceptance threshold may be tuned after looking at outcomes.

## 9. Dilution integration

Tumour-fraction execution must reuse `DilutionPolicy`, `plan_dilution_series` and `execute_dilution_series`.

The v1 design may project desired tumour fractions into a dilution policy, but it does not implement a second mixing engine.

Coverage downsampling is a separate controlled factor and must not be conflated with tumour-fraction mixing.

## 10. Evidence rules

All existing evidence-retention rules remain authoritative:

- caller-native evidence survives;
- normalised evidence is additive;
- exclusions do not delete records;
- QDNAseq 100/500/1000 kbp evidence remains retained;
- alternative fits remain retained;
- no-call/failure/not-assessable remain distinct;
- caller concordance is not truth.

## 11. TDD acceptance

At minimum:

1. factor levels are unique, ordered and valid;
2. the Cartesian positive study-cell count is deterministic;
3. subchromosomal event templates exactly match requested size;
4. whole-chromosome events never receive a numeric size template;
5. synthetic cohort specimens retain exact coverage/fraction continuous values;
6. each design cell/replicate receives a unique deterministic specimen ID and input fingerprint;
7. matrix cutpoints are prospectively derived from factor levels rather than outcomes;
8. incompatible cfDNA ichorCNA is not silently inserted into lcWGS caller locks;
9. registration lock changes if any factor or caller lock changes;
10. dilution policy generated from v1 tumour fractions preserves the existing dilution engine semantics;
11. no biological/clinical validity field is introduced.

## 12. Outputs

Repository-level outputs:

- versioned v1 design JSON schema;
- deterministic example design JSON;
- deterministic synthetic registration example;
- human-readable v1 status document.

Actual caller outputs/evidence remain generated only when a registered study is executed.

## 13. Non-claims

This program does not:

- establish an analytical LoD;
- establish sensitivity/specificity;
- establish clinical validity;
- select a winning caller;
- promote a production default;
- merge automatically;
- change package version merely because software-study tests pass.
