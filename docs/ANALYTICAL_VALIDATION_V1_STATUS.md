# Analytical Validation Program v1 — Status

Date: 2026-09-20  
Branch: `feat/analytical-validation-v1`  
PR: #84 (Draft)  
Base real-tool qualification head: `caca29feab87d6f58dc7d172e42a51a70f036fa0`  
Package version: `0.8.2`

Research Use Only. Human review required.

## What v1 now defines

The repository now has a prospective, deterministic CNV validation study design for the primary analytical axes:

**coverage × tumour fraction × CNV size/event class × caller**

The design is materialised by `src/ontseq_platform/cnv_validation_v1.py` and reuses the existing CNV validation contracts, registration lock, evidence layer, metric engine and dilution machinery.

No second truth, matching, LoD or evidence system was introduced.

## Registered factor levels

### Coverage

- 0.5x
- 1x
- 2x
- 5x
- 10x
- 20x

### Tumour fraction

- 100%
- 50%
- 30%
- 20%
- 10%
- 5%

A 0% normal-only level is projected as a dilution negative control rather than a positive CNV truth cell.

### Subchromosomal event sizes

- 100 kb
- 250 kb
- 500 kb
- 1 Mb
- 5 Mb
- 10 Mb

### Event classes

- deletion
- duplication/gain
- whole-chromosome loss
- whole-chromosome gain

Whole-chromosome events remain a distinct class and are not assigned an arbitrary numeric size bin.

### Replicates

- 3 deterministic synthetic repeats per positive design cell

These repeats exercise reproducibility and registration contracts. They are not statistically independent biological specimens.

## Cell count

Subchromosomal cells:

`6 coverage × 6 fractions × 2 classes × 6 sizes × 3 replicates = 1296`

Whole-chromosome cells:

`6 coverage × 6 fractions × 2 classes × 3 replicates = 216`

Total positive truth cells:

**1512**

## Initial caller scope

The first v1 matrix is intentionally restricted to genome-wide lcWGS-compatible callers:

- QDNAseq+ACE
- Spectre depth-only

QDNAseq retains all registered resolutions:

- 100 kbp
- 500 kbp
- 1000 kbp

With QDNAseq at three resolutions plus one Spectre lane, the default two-caller design represents:

**6048 planned caller-resolution lanes**

This number is a planned evaluation-lane count, not a performance result.

## Why SAVANA, Wakhan and ichorCNA are not forced into this matrix

They are not omitted.

They require different registered input regimes:

- SAVANA: tumour-aware SV/CNA inputs;
- Wakhan: phased tumour CNA, with the qualified breakpoint-assisted path requiring exact upstream breakpoint provenance;
- ichorCNA: cfDNA/ULP-WGS read-depth regime.

Forcing these tools into the same lcWGS denominator would create a methodologically invalid comparison.

They should receive separate prospectively registered validation matrices that reuse the same validation engine.

## Truth and study materialisation

The repository v1 materialises **synthetic study declarations**:

- deterministic specimen IDs;
- deterministic input fingerprints;
- exact continuous coverage values;
- exact continuous tumour-fraction values;
- one locked CNV truth event per positive cell;
- one locked truth-source identity;
- one locked assessability identity;
- repeat-group identity.

Synthetic declarations contain no patient or restricted genomic payload.

## Registration

The builder produces the existing canonical:

- `CnvValidationMatrix`
- `CnvValidationCohort`
- `CnvValidationRegistration`

The existing preregistration lock remains authoritative.

Changing factor levels, truth templates, caller locks, matching policy, code identity or specimen declarations changes the content lock.

`outcomes_unseen=True` remains mandatory.

## Matching and engineering acceptance candidates

v1 uses the existing common comparator and prospectively locks:

- reciprocal overlap: 0.5;
- maximum breakpoint distance: 500 bp;
- copy-number tolerance: 0.5.

The repository registers candidate engineering questions for:

- sensitivity;
- precision;
- false-positive burden;
- no-call rate;
- technical-failure rate;
- copy-number mean absolute error.

These are **prospective engineering candidate limits**, not clinical cutoffs. They must not be tuned after outcome review.

## Tumour-fraction execution

`dilution_policy_from_validation_v1()` projects the v1 tumour-fraction levels directly into the existing dilution engine.

Therefore:
- no second read-mixing implementation exists;
- nominal and observed tumour fraction remain distinct;
- the normal-only control remains explicit;
- the existing technical-LoD limitations remain in force.

## Public repository artifacts

- `schemas/cnv-validation-v1-design.schema.json`
- `examples/validation/analytical_validation_v1.design.json`
- `docs/superpowers/specs/2026-09-20-analytical-validation-v1.md`
- `docs/superpowers/plans/2026-09-20-analytical-validation-v1.md`

The full 1512-specimen registration is generated deterministically from the design rather than committed as a giant duplicated JSON artifact.

## What is still required for actual analytical validation

The software study design is now implemented, but **analytical validation is not yet complete**.

The next execution phase requires authorised, orthogonally characterised source material, for example:
- characterised tumour/normal source BAMs;
- orthogonal CNV truth such as validated array/WGS/FISH/karyotype evidence as appropriate;
- permitted assessability masks;
- measured/verified coverage;
- authorised truth/resource fingerprints.

The registered coverage and tumour-fraction series can then be materialised using controlled downsampling and the existing dilution engine.

Only those outcome data can estimate sensitivity, precision, LoD, failure rates and copy-number error.

## Explicit non-claims

v1 does not:
- claim analytical sensitivity or specificity;
- establish a clinical LoD;
- establish clinical validity;
- select a winning caller;
- use caller majority as truth;
- qualify incompatible callers by forcing them into the lcWGS matrix;
- add patient genomic payloads to Git;
- merge itself;
- bump package version.

Package version remains `0.8.2`.
