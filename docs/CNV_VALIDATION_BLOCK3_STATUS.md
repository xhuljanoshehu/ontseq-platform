# CNV Validation Block 3 — Metric Aggregation Status

## Status

Block 3 implements the Research Use Only metric and acceptance layer above the locked CNV preregistration and Block-2 full-evidence/traceability contracts.

This layer does **not** establish analytical or clinical validity, does not select a production CNV caller, and does not alter caller execution.

## One-time event matching

Primary normalized CNV events are matched to registered truth **once per analytical lane** using the preregistered `BenchmarkThresholds` and the existing deterministic maximum-cardinality CNV comparator.

Stratified metrics consume that locked assignment. They do not re-match events inside coverage, tumour-fraction, event-class, event-size or bin-size strata.

This prevents a stratum boundary from changing which events are considered matched.

## Truth-side and query-side TP semantics

Matched events are represented as explicit truth/query pairs.

For event metrics:

- sensitivity uses truth-indexed TP and FN;
- precision uses query-indexed TP and FP;
- F1 is derived only when both truth-side and query-side denominators are evaluable;
- false-positive burden uses unmatched primary query events per observed analytical lane.

A matched pair that crosses an event-size boundary can therefore contribute to the truth-size stratum for sensitivity and to the query-size stratum for precision without being re-matched.

## Exact mathematical provenance

Every `CnvMetricResult` retains both broad provenance unions and explicit mathematical memberships:

- `numerator_membership`
- `denominator_membership`
- `excluded_membership`

Membership is typed and may contain:

- registered truth-event IDs;
- normalized-event IDs;
- full-evidence IDs;
- analytical lane IDs;
- negative-unit assessment IDs;
- reproducibility-pair IDs;
- explicit truth/query event pairs.

Thus a headline metric can be traced to the exact records that formed its numerator and denominator, while all other retained caller-native evidence remains available in Block 2.

Metric provenance cannot reference truth events, normalized events, full-evidence records, lanes, negative assessments or reproducibility pairs outside the locked report/evidence set.

## Prospective strata

Coverage, tumour/blast fraction and event-size labels are derived only from preregistered cutpoints.

Intervals are deterministic and left-closed/right-open (`[lower, upper)`), with explicit outer bands.

The full prospective event view is:

`caller × genome_build × data_basis × coverage × tumour_fraction × event_class × event_size × bin_size`

Registered strata with no evaluable denominator remain explicit `NOT_EVALUABLE`; they are never silently converted to zero performance.

Continuous coverage and tumour-fraction values remain preserved in Block-2 evidence and are not replaced by their band labels.

## Run-state metrics

`OBSERVED`, `NO_CALL`, `FAILED` and `NOT_ASSESSABLE` remain distinct states.

No-call, technical-failure and not-assessable rates use the number of actually executed lanes as their denominator. The theoretical prospective matrix cross-product is not used as an execution denominator.

A registered technical stratum with no executed lane remains `NOT_EVALUABLE`.

## Specificity

Specificity is **not** inferred from event-level FP counts.

It requires an explicit `CnvNegativeUnitAssessment` that is bound to:

- registration SHA-256;
- evidence-manifest SHA-256;
- analytical lane;
- specimen;
- registered negative-universe ID and resource SHA-256;
- the same assessability-mask SHA-256;
- the full registered assessable-unit count;
- explicit full-evidence IDs.

For the current implementation, every observed lane contributing to a specificity view requires a complete negative-unit assessment. Otherwise specificity remains `NOT_EVALUABLE`.

Specificity is:

`true_negative_units / assessed_units`

where:

`true_negative_units = assessed_units - false_positive_units`.

No event-class/event-size-specific specificity is claimed without a separately registered negative-unit classification.

## Copy-number, cellularity and ploidy error

Copy-number error uses only matched truth/query CNV pairs for which both sides have quantitative copy-number values.

The metric reports mean absolute error and signed bias with an explicit eligible-pair denominator.

Cellularity and ploidy error are evaluable only when:

- orthogonal quantitative truth was explicitly registered and linked to a registered truth source; and
- exactly one caller fit in the lane is marked as selected.

Alternative fits remain retained in full evidence and are not promoted into the primary quantitative metric.

## Reproducibility

Reproducibility is an engineering concordance metric, not biological truth.

Only lanes in the same registered repeat group and compatible caller/build/data-basis/bin context are compared.

A repeat group may not silently mix different biological specimens.

Pairwise event concordance uses the preregistered CNV matching thresholds. Two empty call sets are `NOT_EVALUABLE`, not artificial 100% concordance.

The metric denominator retains the exact evaluable reproducibility-pair IDs; non-evaluable pairs remain explicit excluded provenance.

## Acceptance questions

Acceptance questions are defined prospectively in `CnvValidationMatrix.acceptance`.

Question filters must match a canonical metric stratum exactly. Evaluation uses:

- `PASS`
- `FAIL`
- `NO_CALL`
- `NOT_EVALUABLE`
- `NOT_APPLICABLE` where appropriate.

A metric below the registered `minimum_evaluable_denominator` is `NOT_EVALUABLE`; no acceptance value is invented.

Lower and upper acceptance bounds are applied deterministically without changing the underlying metric.

## Sealed validation report

`CnvValidationMetricReport` binds:

- registration SHA-256;
- Block-2 evidence-manifest SHA-256;
- preregistered matching thresholds;
- one-time lane event assignments;
- all emitted metric views;
- acceptance results;
- negative-unit assessments;
- reproducibility pairs;
- exact metric numerator/denominator/exclusion provenance.

The report has a canonical SHA-256 content lock. Any change to assignments, metrics, memberships, acceptance results or supporting assessment/pair objects changes the report content and invalidates the old lock.

No runtime timestamp participates in the report lock.

The public report schema is versioned at:

`schemas/cnv-validation-report.schema.json`

## Explicit non-claims

Block 3 does not:

- run QDNAseq+ACE, ichorCNA or Spectre on real specimens;
- establish clinical sensitivity, specificity or LoD;
- claim analytical validity from synthetic tests;
- choose a winning/production caller;
- mutate preregistered matching thresholds or stratum cutpoints;
- discard alternative caller-native values;
- wire a CNV caller into the canonical production pipeline;
- merge PR #79 or change the package version.

Package version remains `0.8.2`.

Real intended-use validation remains dependent on prospectively locked, orthogonally characterized material and the predefined acceptance/truth/mask/coverage/purity design.
