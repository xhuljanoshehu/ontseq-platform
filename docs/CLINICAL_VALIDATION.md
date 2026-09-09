# Analytical and clinical validation plan

## Current classification

The software is a research prototype. No output is validated for diagnosis, prognosis,
treatment selection or patient release.

### Target interval identity correction, 2026-09-09

The interval-identity fix reused from the methylation handover branch (`b2d5256`)
changes biological output for target BEDs containing repeated region labels. Each
reported row now aggregates only sites within its own chromosome and half-open
interval. Shared sites still contribute independently to overlapping targets; empty
targets retain null fractions rather than inheriting another target's measurement.
Chromosome summaries and targets with unique labels retain their numerical behavior.

New normalization records `ontseq_region_assignment=interval-identity-v1` in the
existing tool-parameter provenance map, alongside the BED and bedMethyl fingerprints
and reported interval coordinates. The method also enters the stage resume signature,
invalidating previous method-unversioned results even in a local uncommitted checkout.
Historical reports are not relabeled on loading;
their missing method marker must not be interpreted as evidence of this correction.
The serialized schema, modkit command, confidence/coverage policies and mixture
estimators remain unchanged. Prospective runs require a fresh source/runtime lock;
existing affected target reports require regeneration before comparison or reuse.

A synthetic regression covers repeated labels on disjoint and overlapping intervals,
the same coordinates on different chromosomes, missing measurements and preservation
of input tool metadata. An integration regression verifies recomputation after the
method upgrade and normal resume thereafter. This checks aggregation correctness, not biological recovery
or clinical performance. Analytical promotion still requires independent review and
the applicable study gates below.

### Nanopolish import hardening impact, 2026-09-05

Adapter `nanopolish-call-table-v2` changes acceptance of oversized or malformed inputs:
the existing byte budget now separately bounds compressed and decoded data, including
blank lines. Header/physical-line/column limits apply before CSV allocation; read-name
control characters are rejected to make the existing selection serialization unambiguous
within its accepted domain. Ordinary read selection hashes and biological estimator
equations remain unchanged. Source provenance adds the input security profile and actual
decoded byte count; exported schemas are updated. The fixed study matrices are unchanged.
Fresh registrations must bind v2 and the new source/runtime digest; historical evidence
must not acquire v2 provenance merely by being reloaded.

Tiny synthetic regressions cover exact and exceeded byte bounds, gzip expansion, blank
lines, header/data limits, multiline headers and read-identifier controls, including the
previously ambiguous selection pools. Recovery, adapter and CLI suites remain passing.
An approved real-data structural intake was then run locally without marker selection
or biological recovery. Its detailed files and outcomes remain in ignored `work/`.
The intact-read pool check refuses a matrix whose largest budget cannot be supplied;
neither resampling nor post-outcome threshold changes are introduced. This is software
hardening and source intake, not biological or clinical validation.

## Validation units

### 0.7.1 methylation discovery impact

The Desktop preparation correction distinguishes inaccessible input drives/directories,
output paths, resources and runtime components before starting a probe. Structured local
prerequisite codes and a preparation-phase flag replace the misleading worker-failure
classification. Synthetic filesystem failures and presentation regressions are required.
This correction does not change the reader, probe provenance, biological calculations,
reference contracts or analysis admission requirements; an unstarted probe stays unknown.

The v2 discovery path uses a direct BAM reader, indexed quick sampling and an optional
sequential background scan. Discovery now distinguishes supported paired 5mC information,
incomplete/unsupported tags, reader failures, cancellation and a successful exhaustive
negative. This can change whether the interface offers the methylation module and therefore
requires synthetic positive/negative, long-read, late-positive, file-change, index,
timeout/cancellation and service-scope regression checks. The limits are engineering
resource bounds, not clinical thresholds or evidence of assay sensitivity.

The analytical methylation adapter, estimator equations, reportability policy, thresholds,
reference families and result schema remain unchanged. Probe v2 metadata records its method,
reason, completeness, elapsed time and inspected-record count. A stat fingerprint detects
ordinary file changes; it is not represented as a content checksum. Cached positive evidence
is freshly inspected before an analysis start. A prospective analytical study still needs
its own locked software identity and independent clinical validation.

### 0.7 integration and optional methylation impact

The follow-up Desktop setup correction changes diagnostic text transport and presentation.
It does not change reference selection, analytical policies or biological calculations.
Existing local GRCh38 resources may be imported beside GRCh37 only after the ordinary
authority and checksum checks; DNS failures never justify a different build or relaxed hashes.

The integration retains the newer supported-build reference contracts, corrected CNV
coordinate exports, annotation context and structured CNV-only ISCN proposal path. It
adds the separately developed methylation modules without changing their estimator
equations or versioned technical thresholds. Software identifiers and source/runtime
hashes change; prospective studies require a fresh lock before outcome inspection.

Adaptive Sampling explicitly selects `modkit-cpg-targets-technical-v1`, using the
profile's locked analysis target BED instead of chromosome-wide summaries. It preserves
the 0.8 call-confidence and 5x valid-coverage engineering defaults from the lcWGS policy;
off-target intervals are outside this report. This assay-specific selection and its
policy identity require their own prospective validation, including missing targets.

The BAM probe detects the presence of candidate MM/ML annotations, not the biological
correctness of those annotations or the modification model. Its bounded negative result
is unknown unless EOF and a successful tool exit establish a complete inspection. The
operator's explicit inclusion decision changes the requested module set. Missing tools
and incomplete evidence must remain visible and may not silently enable methylation.

Regression checks must cover both opt-in and opt-out, file changes after a decision,
malformed/truncated inputs, resource and artifact identity, and stale-result rejection.
Previously generated methylation files are excluded after module deselection or failure;
this intentionally prevents stale evidence from being presented as a current result.

The Windows process-lifetime gate now uses a read-only synchronization handle instead
of POSIX signal-zero semantics. This is an execution-integrity fix, with no change to
the biological algorithms. Real binary interoperability and biological recovery remain
separate from synthetic software and user-interface verification.

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

## Paired-source methylation validation impact

The standalone Nanopolish mixture path adds a new quantitative research output. It estimates the
coefficient of source A in a two-source linear methylation model and evaluates it against the
realized fraction of held-out read groups drawn from source A. It does not estimate tumour
purity, blast fraction, cell fraction or DNA-mass fraction. A `COMPLETED` level means only that
the locked marker, coverage, fit and interval-width gates passed. `NO_CALL` is not zero signal or
a negative biological result. At report level, `COMPLETED` requires the entire requested grid,
`PARTIAL` means that only some levels were quantifiable and `NO_CALL` means that none were. Bias,
MAE and RMSE are explicitly limited to completed levels and must be interpreted together with the
reported quantifiable fraction.
The separate technical recovery `PASS`/`FAIL`/`NOT_EVALUABLE` field compares those summaries and
empirical interval coverage with versioned engineering limits. It is a regression aid, not an
analytical-validation verdict or a claim that the limits are suitable for patient samples.

The included interval propagates finite methylation-call counts conditional on one source pair,
its selected exact markers and its disjoint calibration/mixture split. It does not cover donor
variation, correlated CpGs on a molecule, cell composition, wet-lab dilution, library effects,
upstream caller or platform changes, CNV, ploidy or allelic imbalance. Seeded replicates reuse the
same held-out pools and are not independent specimens. The standardized residual RMSE is a
technical fail-closed heuristic; because methylated and unmethylated channels are correlated, it
is not a calibrated chi-square test or biological acceptance criterion.

Promotion to a biological source-fraction assay requires, at minimum:

- a frozen intended quantity, upstream protocol, genome build, marker policy and software
  version, with checksummed reference data;
- independent donors and held-out samples rather than only resampling two source files;
- wet-lab mixtures with orthogonal source proportions across input mass and coverage;
- predefined bias, precision, interval-coverage and no-call acceptance criteria;
- within-run, between-run, operator, instrument, flow-cell, reagent-lot and upstream-caller
  reproducibility;
- interference studies for cell composition, DNA quality, CNV, ploidy, allelic imbalance and
  correlated methylation on individual molecules;
- external validation of any mapping from the source coefficient to tumour or blast fraction.

The current policy values, including the 5, 10 and 20 percent mixture levels, are technical test
points. They are not detection limits. Public methylation datasets may support engineering and
generalisation studies, but they remain outside Git and do not replace intended-use validation.

## Independent methylation recovery and new adapter impact (2026-09-05)

Biological sample, donor, source and cell-type identity fields now reject Unicode
control (`Cc`) and format (`Cf`) characters before the existing whitespace
normalization. A NUL-suffixed duplicate donor label previously passed the independent
donor equality gate in a synthetic registered experiment. The new registration
regression rejects it; identifiers retain their original case. This is an admission
correction, not verification of donor identity from biological measurements. Valid
canonical cohorts and numerical estimator outputs are unchanged. Previously locked
source/runtime versions require a new registration before further outcomes.

The adapter `0.1.1` follow-up tightens the supported input boundary: explicit
Dorado duplex offspring are rejected using `dx`, including when a recognizable
program header is absent, and modkit alignment extents must fit the locked FASTA.
Previously accepted out-of-scope input can now fail intake. Valid simplex m/h
probabilities, estimator equations and the frozen validation matrices are unchanged.
Synthetic regressions cover absent-header duplex tags, simplex-parent acceptance,
malformed tags and an invalid span whose reported CpG itself is valid. Registrations
must bind the new adapter version and source/runtime checksum. Actual binary BAM
and modkit interoperability and biological recovery remain unestablished.
The TSV export cannot independently establish `dx` status; simplex-parent reads
also do not guarantee independent physical molecules. Correlation uncertainty remains
outside the current conditional interval.

Two complementary study scopes now exist. The prospective two-source holdout experiment
tests recovery within the two declared sources using disjoint read groups; it does not
require four independent donors. Its PASS cannot be promoted to donor-independent or
population validation, and its report fixes `donor_independent_validation=NOT_EVALUATED`.
The separate four-sample lane adds non-overlapping calibration/test donors when such
material is available. This distinction preserves the requested staged validation rather
than making a later generalisation study a prerequisite for the immediate experiment.
See `docs/METHYLATION_HOLDOUT.md` for the explicit split and scope contracts.

The per-read modkit path is pinned to the source-reviewed `extract full` schema in 0.6.1.
Probability midpoint reconstruction and the exclusive alignment end affect normalization
at boundaries and have synthetic regression coverage. An actual binary modkit comparison
and biological source recovery are still required before interoperability or analytical
performance can be claimed.

The next validation lane separates calibration source A/B from test source A/B in the experiment
contract. Biological sample and donor declarations, input SHA-256 values, protocol metadata,
software content and the complete fraction/depth/seed matrix must be locked before execution.
A registration timestamp and operator declaration are an audit record, not a trusted timestamp
service or proof that an operator has never examined the test data. Seeded mixtures remain
computational replicates; they do not increase the number of independent donors.

The Dorado/modBAM and modkit per-read adapters alter the biological input representation and
therefore require revalidation. In particular, strand orientation, CpG reference matching,
alignment operations, modification-probability thresholds, ambiguous/missing calls and the
distinction between 5mC and 5hmC affect the estimated coefficient. Aggregate bedMethyl cannot
reconstruct intact read groups and is not a substitute for these per-read inputs. Nanopolish
multi-motif markers must not be silently pooled with single-CpG Dorado markers.

The acceptance matrix is a prospective research design, not an empirically validated set of
clinical thresholds. A numerical estimate can fail recovery; absent independent evidence must
remain `NO_CALL`. Report bias, MAE, RMSE and conditional interval coverage together with the
unconditional `NO_CALL` rate and fraction/depth strata. Do not tune marker filters or acceptance
limits on the held-out experiment and then describe that same experiment as prospective.

The GSE210260 outcomes provided in the project handoff remain negative controls for interpretation:
8 qualifying markers under the strict policy gave 0/18 quantifiable levels; the relaxed feasibility
analysis gave 10/10 numerical levels but MAE 0.256, RMSE 0.286 and 10% interval coverage, hence
`FAIL`. These are historical handoff results, not a newly reproduced independent experiment.

Revalidation must include synthetic adapter conformance tests and actual independently prepared
biological inputs, followed by additional donor panels and reference mixtures for the intended
use. Calibration, donor, marker-selection and within-molecule correlation uncertainty remain
unresolved. No source coefficient may be relabelled DNA mass, cell fraction, tumour purity or
subclone fraction without independent reference measurements and a successful additional
calibration study. SNV, CNV and ploidy remain separate evidence layers; 10–20% remains a hypothesis.
## ISCN proposal-path validation impact

Result schema `0.3.0` represents ISCN proposal execution explicitly as `NOT_REQUESTED`,
`NOT_ASSESSED`, `NO_RENDERABLE_CANDIDATE` or `PARTIAL_EVENT_LEVEL`. The regular engineering path
can produce only CNV event fragments from typed CNV output with matching, checksum-bound GRCh37 or
GRCh38 reference-lock, cytoband and annotation-cache provenance. Missing upstream evidence, failed
QC or resource disagreement blocks assessment. An assessment with no renderable fragment is a
technical `NO_CALL`; it is not a normal or negative karyotype.
The CNV artifact is also bound to the manifest sample/build. `+chr`/`-chr` requires an exact
zero-to-contig-end span, not only the wider CNV whole-chromosome fraction threshold.

This proposal path has not been analytically or clinically validated. It does not infer or assert
chromosome count, sex-chromosome complement, clonality, phase, derivative structure, balance,
normality or a complete karyotype. Automatic SV/BND/translocation/inversion-to-ISCN conversion is
disabled. Every proposal therefore remains `requires_expert_review=true` and
`clinical_release_allowed=false`.

Before expanding or clinically using this path, validation must additionally include:

- authorized, edition- and errata-controlled ISCN source material and an expert-approved expected
  result corpus;
- independent GRCh37 and GRCh38 positive, negative, ambiguous and boundary cases tied to exact
  reference, cytoband and annotation-cache checksums;
- parser/renderer round trips and regression tests for every supported construct;
- event-to-fragment traceability, duplicate handling, partial-band and whole-chromosome edge cases;
- explicit validation of sex-chromosome, clone/mosaic, phase, reciprocal/balanced and derivative
  semantics before any of them is enabled;
- user-interface evaluation showing that partial, unassessed and no-call states cannot be mistaken
  for a complete or normal karyotype;
- expert sign-off, change control and a release mechanism that cannot bypass required review.

The introduction of structured proposal states improves safety and auditability but does not
satisfy any of these clinical release gates or establish full ISCN 2024 conformance.
