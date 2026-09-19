# CNV validation program design

Date: 2026-09-17
Status: design approved in chat; implementation not started
Branch: `feat/cnv-validation-program`
Base: `main` at `7f7072c2e8033d0cad54c7b913ed6fc1da5d6e90`
Scope: Research Use Only

## 1. Purpose

Build a prospective, reproducible, fail-closed validation layer for ONTSeq copy-number analysis that evaluates multiple CNV callers under identical truth, specimen, reference, masking and matching rules.

The validation layer must support:

- orthogonally characterised samples and explicitly versioned truth;
- GRCh37 and GRCh38 as separate reference strata;
- coverage strata;
- tumour/blast-fraction strata with source/method/timepoint metadata;
- event-class and event-size strata;
- QDNAseq+ACE 100/500/1000 kbp resolution comparison;
- caller comparison, initially including QDNAseq+ACE, ichorCNA and Spectre where their assumptions match the evaluated data;
- cellularity/purity and ploidy error where a caller emits such estimates;
- sensitivity/recall, precision/PPV, F1, false-positive burden, no-call rate and technical failure rate;
- specificity only when a prospectively defined assessable negative denominator exists;
- LoD characterisation through the existing dilution infrastructure;
- within-run and between-run reproducibility;
- complete traceability from aggregate metrics back to every retained caller result.

No caller is promoted to a production default by this design alone.

## 2. Non-negotiable evidence-retention rule

Prospective stratification must never delete, overwrite or hide the underlying analytical evidence.

The system will separate three layers:

1. **Full caller evidence** — all retained outputs that the caller or adapter exposes, including primary and alternative fits, per-bin/per-segment values when available, QC, parameters, versions, reference identity and execution provenance.
2. **Normalised comparison evidence** — immutable derived records that map caller-specific output into the common ONTSeq event/validation contract without replacing the original values.
3. **Aggregated validation views** — overall and stratified metrics generated from the normalised evidence.

Exclusion from a primary metric is a classification, not deletion. Records may be labelled for example:

- `USED_FOR_PRIMARY_ANALYSIS`
- `SECONDARY`
- `EXPLORATORY`
- `EXCLUDED_FROM_PRIMARY_METRIC`
- `NOT_ASSESSABLE`
- `NO_CALL`
- `FAILED`

Every aggregate value must be traceable to the exact records contributing to its numerator and denominator.

## 3. Existing building blocks to reuse

The design deliberately reuses rather than rewrites the current ONTSeq components:

- `src/ontseq_platform/benchmark.py` for deterministic event-level one-to-one matching and TP/FP/FN/precision/recall/F1;
- `BenchmarkCase` and existing benchmark schemas for event-level comparison inputs;
- the existing in-silico dilution and LoD framework for controlled tumour-fraction stress testing;
- `configs/cnv/qdnaseq_ace.technical.yaml` for the current 100/500/1000 kbp QDNAseq+ACE technical profile;
- existing caller-specific raw/result models and provenance;
- analysis-profile candidate declarations for QDNAseq+ACE, ichorCNA and Spectre.

The current benchmark comparator remains a low-level primitive. The new validation layer sits above it.

## 4. Architecture

Create an independent `cnv_validation` subsystem with four conceptual contracts:

### 4.1 Validation matrix

A prospectively defined study matrix declaring the planned dimensions and acceptance semantics before outcome review.

Required dimensions:

- caller;
- genome build;
- coverage stratum;
- tumour/blast-fraction stratum;
- event class;
- event-size stratum;
- bin size where applicable;
- replicate/repeat identifier;
- specimen/input-path category where relevant.

The matrix must distinguish:

- primary strata;
- secondary strata;
- exploratory strata;
- minimum evaluable denominator per stratum;
- rules for `NOT_EVALUABLE`.

The logical cross-product is not required to be fully populated. Missing combinations must remain visible as missing or not evaluable, not silently omitted.

### 4.2 Cohort and truth contract

Each validation specimen must carry explicit, non-inferred metadata:

- pseudonymised validation specimen identifier;
- authorised/public access basis;
- specimen/material type;
- reference build;
- input fingerprint(s);
- orthogonal truth method(s);
- truth resource/version/fingerprint;
- tumour/blast fraction;
- tumour/blast-fraction source or measurement method;
- measurement timepoint;
- technical coverage measure and definition;
- assessable genomic territory/mask;
- positive truth CNV events;
- prospective negative denominator definition where specificity will be estimated.

Unknown is a valid explicit state. No default tumour fraction, ploidy, coverage, truth method or assessability may be guessed.

Truth data must remain logically separate from caller output. Caller concordance is never truth.

### 4.3 Registration lock

Before outcome evaluation, create a registration object that locks by content hash:

- validation matrix;
- cohort/truth declaration;
- matching policy;
- masks/exclusions;
- caller versions and relevant adapter policy identities;
- reference identities;
- code identity/software version;
- prospective acceptance rules;
- declaration that test outcomes were not used to modify the registered rules.

Registration timestamps are operator declarations unless an external trusted timestamp mechanism is introduced later.

### 4.4 Evidence and report

The execution/evaluation layer produces:

- one full-evidence record set;
- one normalised comparison record set;
- one aggregate validation report.

The aggregate report must never be the sole surviving result.

## 5. Full evidence table

A machine-readable full-evidence table is mandatory.

At minimum, each row/record must be addressable by:

`specimen × caller × build × coverage × tumour_fraction × bin_size × event_or_fit × replicate`

Fields vary by record type, but the retained evidence should include where available:

- original caller event/segment identifier;
- coordinates and event type;
- estimated copy number;
- raw or corrected depth-derived values;
- log2 ratio or equivalent caller-specific signal;
- bin size;
- selected fit;
- alternative fits;
- cellularity/purity estimate;
- ploidy estimate;
- fit error/objective;
- caller QC status;
- caller parameters;
- caller and dependency versions;
- reference and resource fingerprints;
- normalisation/adaptation status;
- matching status against truth;
- reason for exclusion/non-assessability/no-call/failure;
- linkage to the aggregate strata to which the record contributes.

Caller-native artefacts that are too large for the tabular report may remain as separately fingerprinted evidence files, but the report must retain their identities and relationship to derived rows.

## 6. Prospective stratification

### 6.1 Primary biological/analytical strata

The first validation design should prioritise:

- tumour/blast fraction;
- coverage;
- event class;
- event size.

These are expected to materially affect sensitivity and false-negative risk.

### 6.2 Technical strata

Technical comparison strata include:

- caller;
- GRCh37 vs GRCh38;
- bin size/resolution;
- caller parameterisation where prospectively registered.

### 6.3 Event classes

Initial CNV classes:

- whole-chromosome gain;
- whole-chromosome loss;
- subchromosomal deletion;
- subchromosomal duplication/gain.

Additional classes require explicit extension rather than being silently folded into these groups.

### 6.4 Event-size bins

Exact cut-points must be prospectively locked in the validation matrix before outcome review. The implementation must support configurable ordered bins and a distinct whole-chromosome class.

The software must not optimise size cut-points after seeing validation outcomes.

### 6.5 Coverage and tumour-fraction bins

Likewise, exact coverage and tumour/blast-fraction bins are registered values, not hard-coded scientific truths.

A measured continuous value must remain preserved in the evidence record even when a categorical stratum is derived from it.

## 7. Metrics

### 7.1 Event-level metrics

Reuse the existing deterministic event comparator for:

- true positives;
- false positives;
- false negatives;
- sensitivity/recall;
- precision/PPV;
- F1;
- false-positive burden.

Matching thresholds must be registered before comparative outcome review.

### 7.2 Specificity

Specificity must not be inferred from event-level TP/FP/FN counts.

It may be calculated only if the validation design prospectively defines an assessable negative denominator, for example a locked genomic-region or bin-level truth universe with explicit masking and exclusion rules.

If a valid TN denominator is absent, specificity is `null`/`NOT_EVALUABLE` with an explanatory reason.

### 7.3 No-call and technical failure

No-call, technical failure and biological negative findings remain distinct.

Required rates:

- no-call rate;
- technical failure rate;
- not-assessable rate where relevant.

These denominators must be reported explicitly.

### 7.4 Copy-number error

Where truth supports quantitative copy number, report quantitative error separately from event detection.

Potential summary measures include absolute error and bias. Exact acceptance limits are prospectively configured rather than embedded as clinical defaults.

### 7.5 Cellularity/purity and ploidy

When caller outputs include cellularity/purity or ploidy estimates, retain:

- selected estimate;
- all retained alternative fits;
- fit/objective value;
- error against orthogonal/reference value if available;
- whether the value contributed to primary, secondary or exploratory analysis.

A caller lacking these outputs is not penalised by invented values; the corresponding metric is `NOT_APPLICABLE` or `NOT_EVALUABLE` as defined by the report contract.

### 7.6 LoD

Use the existing dilution-series framework rather than creating a second LoD implementation.

The validation report may reference/import registered LoD evidence but must preserve the current distinction between technical in-silico characterisation and wet-lab/intended-use LoD.

### 7.7 Reproducibility

Support explicit within-run and between-run repeat groupings.

Reproducibility metrics may include:

- event concordance;
- copy-number difference;
- selected-ploidy/cellularity difference;
- no-call/failure discordance.

Exact acceptance criteria are registered study parameters.

## 8. Caller comparison

Initial candidate set:

- QDNAseq+ACE;
- ichorCNA;
- Spectre.

A caller is compared only in data regimes compatible with its documented assumptions.

All evaluated callers must share the same:

- specimen set;
- truth lock;
- reference build;
- assessability mask;
- event matching policy;
- primary strata;
- exclusion rules;
- evaluation code version.

Caller-specific execution parameters remain visible and versioned.

No caller may receive more favourable matching thresholds solely because its native representation differs. Representation adaptation happens before the common comparator and is itself auditable.

The first validation implementation must not select a production default automatically.

## 9. QDNAseq+ACE multi-resolution handling

QDNAseq+ACE 100/500/1000 kbp results are all retained.

The current technical primary view may remain 500 kbp for engineering continuity, but validation reports must expose all three resolutions and must not discard 100 or 1000 kbp values.

Required comparisons include:

- per-resolution event performance;
- chromosome-level concordance across resolutions;
- resolution-specific subchromosomal sensitivity/FP burden;
- resolution-specific no-call/failure status;
- selected and alternative ACE fits by resolution where available.

Any future decision to promote one resolution for a specific event class must be evidence-based and separately approved.

## 10. Report decisions and states

The report must distinguish at minimum:

- `PASS`
- `FAIL`
- `NO_CALL`
- `NOT_EVALUABLE`

These decisions apply to registered acceptance questions, not to individual biological truth claims.

A missing stratum is not a pass.

A stratum below its registered minimum denominator is `NOT_EVALUABLE`.

A technical run failure is not converted to a biological negative.

## 11. Outputs

Planned machine-readable outputs:

- `cnv-validation-matrix.json`
- `cnv-validation-cohort.json`
- `cnv-validation-registration.json`
- `cnv-validation-evidence.json` or equivalent indexed evidence manifest
- `cnv-validation-full-evidence.csv`/parquet-compatible representation
- `cnv-validation-report.json`
- human-readable summary tables generated from the same report object

JSON schemas are exported and versioned.

The full evidence table is an output requirement, not an optional debug artefact.

## 12. Privacy and data boundary

No patient genomic payloads, BAM/POD5/FASTQ/VCF content, identifiers or institutionally restricted truth data are committed to GitHub.

Repository tests use synthetic fixtures only.

Real validation runs occur only on authorised data in the appropriate environment. Versioned repository objects may contain schemas, study templates, code and non-sensitive hashes/identifiers when permitted, but not restricted biological payloads.

## 13. Failure behaviour

The validation system is fail-closed.

Examples:

- build mismatch -> reject comparison;
- reference identity mismatch -> reject comparison;
- truth/mask fingerprint mismatch -> reject comparison;
- unregistered caller/policy -> do not silently evaluate as registered;
- missing tumour fraction where required by a primary stratum -> explicit `NOT_EVALUABLE`/eligibility failure, not imputation;
- malformed/non-finite quantitative values -> validation error;
- absent negative denominator -> specificity unavailable, not zero;
- absent caller-native output required by the registered question -> explicit missing/failed state.

## 14. Implementation decomposition

Implementation proceeds in small reviewable blocks.

### Block 1 — contracts and preregistration

Implement:

- CNV validation matrix models;
- cohort/truth models;
- prospective strata definitions;
- acceptance-question models;
- content hashing and preregistration lock;
- eligibility validation;
- schema export;
- synthetic unit tests.

No real biological data and no caller execution changes.

### Block 2 — full evidence and traceability contract

Implement:

- full-evidence record model;
- caller-native artefact references/fingerprints;
- normalised event linkage;
- contribution/exclusion status;
- aggregate-to-source traceability tests.

### Block 3 — aggregation metrics

Implement:

- overall and stratified TP/FP/FN metrics;
- denominator accounting;
- no-call/failure/not-assessable rates;
- specificity only with valid negative universe;
- quantitative CN/cellularity/ploidy error where applicable;
- explicit missing-stratum reporting.

### Block 4 — QDNAseq+ACE multi-resolution adapter

Map all 100/500/1000 kbp outputs and retained ACE alternatives into full evidence without discarding caller-native artefacts.

### Block 5 — independent caller adapters

Add/evaluate ichorCNA and/or Spectre adapters under the same registered evidence contract, one caller at a time, subject to runtime/data-regime suitability.

### Block 6 — real validation execution

Only after the contracts and tests are stable:

- run authorised/public validation cohorts;
- execute pre-registered analyses;
- generate locked reports;
- review failures without changing acceptance criteria post hoc.

## 15. Test strategy

Repository tests must use synthetic data and verify at least:

- matrix rejects duplicate/overlapping invalid strata definitions;
- continuous source values are preserved when categorical strata are derived;
- missing strata appear explicitly;
- registration hash changes when any matrix/truth/matching rule changes;
- registration rejects mismatched content hashes;
- unknown tumour fraction is never imputed;
- specificity is unavailable without a valid negative denominator;
- full-evidence records survive exclusion from primary metrics;
- 100/500/1000 kbp records remain independently addressable;
- alternative cellularity/ploidy fits remain retained;
- aggregate numerators/denominators can be traced back to contributing record IDs;
- no-call, failed, not-assessable and biological negative remain distinct;
- caller comparison uses identical registered matching policy;
- non-finite values are rejected;
- no repository fixture contains real patient data.

## 16. Acceptance criteria for the subsystem

The engineering subsystem is acceptable for real study execution when:

1. study design and truth are content-locked before outcome evaluation;
2. all primary and secondary strata are declared prospectively;
3. every measured continuous value remains accessible in full evidence;
4. every aggregate metric exposes its denominator and traceable contributing records;
5. specificity cannot be produced without a valid negative denominator;
6. all caller/native and normalised evidence is linked by immutable identities;
7. multi-resolution QDNAseq+ACE evidence is preserved in full;
8. no caller is automatically selected as production default;
9. synthetic contract/aggregation tests pass;
10. the subsystem makes no clinical-validity claim.

## 17. Explicit non-goals

This design does not:

- validate a clinical assay by itself;
- define universal clinical thresholds;
- choose QDNAseq+ACE, ichorCNA or Spectre as the winner;
- replace orthogonal truth;
- treat caller concordance as biological truth;
- infer specificity from absence of calls alone;
- make Adaptive Sampling off-target equivalent to lcWGS;
- merge CNV validation with MARLIN, ISCN or fusion validation;
- modify package version or release state.

## 18. Design decision summary

The central design choice is to preserve the current low-level benchmark engine and add a separate prospective validation layer above it.

The validation layer is intentionally evidence-complete: stratified reports are views over the data, never destructive filters. A reviewer must always be able to move from a headline metric to the exact specimen/caller/build/coverage/tumour-fraction/bin/event/replicate records and, from there, to the fingerprinted caller-native artefact.

This traceability requirement is part of the validation contract, not a later UI enhancement.
