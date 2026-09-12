# ONTSeq Independent Validation Audit — Design

## Status
Draft for user review before implementation.

## Goal
Add an independent, reproducible validation layer that measures whether ONTSeq’s existing regression suite can detect biologically meaningful faults and that establishes a clean path toward external somatic truth benchmarking without changing current production analysis behavior.

## Scope
The first implementation will add three audit lanes:

1. **Property/metamorphic audit** for invariants that should hold under representation-preserving input transformations.
2. **Mutation-testing audit** focused on high-risk biological logic, initially opt-in/non-blocking until a baseline mutation score is established.
3. **External somatic benchmark harness** designed around versioned public Cancer Genome in a Bottle-style truth bundles; code and schemas are implemented now, while large benchmark payloads remain external to Git.

The blind AML holdout is specified as the next validation phase but is not implemented in this change because it requires independently characterized AML material and a pre-registered study lock.

## Non-goals
- No change to CNV, SV, methylation, ISCN or reportability algorithms.
- No new clinical claims or promotion of any existing output to validated/reportable.
- No committing large BAM/CRAM/POD5/truth payloads to Git.
- No tuning thresholds against the external benchmark outcomes in the same run used to report performance.

## Architecture

### A. Property and metamorphic audit
Create a dedicated audit module and test suite that exercises deterministic invariants across existing normalized contracts. Candidate invariants include:

- input ordering does not alter normalized biological results;
- BED row ordering does not alter target-level results;
- reference-build mismatch fails closed;
- genomic intervals remain valid (`0 <= start < end <= contig_length` when contig length is known);
- equivalent event ordering produces the same canonical normalized result;
- report rendering does not mutate normalized result objects;
- repeated execution with the same locked source/runtime identity is deterministic;
- resume does not recompute completed content-addressed stages;
- schema round-trips preserve semantic result identity.

Tests will use Hypothesis where generated data is useful and explicit metamorphic fixtures where realistic file structure matters. Tests must not encode the implementation under test as their oracle.

### B. Mutation-testing audit
Add a reproducible mutation-testing configuration using `mutmut` (or an equivalent maintained Python mutation tool if compatibility blocks it). Initial mutation targets are limited to high-risk Python packages with deterministic unit coverage:

- reference/coordinate locking;
- CNV normalization/quantitation helpers;
- SV consensus/evidence normalization;
- ISCN proposal/normalization logic;
- reportability/release gates;
- methylation quantitative normalization where pure-Python tests exist.

The mutation job will be opt-in/manual and artifact-producing first. It becomes blocking only after a baseline is measured, surviving mutants are triaged, and a threshold is explicitly accepted.

### C. External somatic benchmark harness
Add a small manifest-driven benchmark runner that consumes external files by path and checksum. The harness will not download patient or truth data automatically in the normal CI lane. A benchmark manifest binds:

- dataset/truth-set identifier and version;
- genome build;
- input BAM/CRAM/VCF paths and SHA-256 hashes;
- truth VCF/BED paths and hashes;
- region inclusion/exclusion masks;
- caller/ONTSeq version, git commit, reference bundle and policy identities;
- predefined comparison parameters;
- output directory and run identifier.

The runner produces machine-readable JSON plus TSV summaries stratified by event class, including where available: TP, FP, FN, precision, recall, F1, false positives per sample, no-call fraction and breakpoint-distance summaries. It must distinguish technical benchmark evidence from intended-use AML validation.

A public cancer truth set such as NIST Cancer Genome in a Bottle HG008-T is the preferred first external target, but its exact accepted release and comparison rules must be locked in the benchmark manifest before any outcome is inspected.

### D. Future blind AML holdout
A follow-on prospective lane will require:

- independently characterized AML positives and negatives;
- orthogonal cytogenetic/FISH/PCR/RNA evidence by event class;
- pre-registered software/reference/policy lock;
- fixed acceptance criteria before outcome access;
- explicit stratification by coverage, blast/tumour fraction, event type and assay mode.

This phase is intentionally separated from the engineering audit so public technical benchmarks are not mislabeled as clinical validation.

## CI and artifacts

Add an independent workflow rather than overloading the already broad main CI:

- **property-audit**: runs on pull requests and main, fast enough to be blocking;
- **mutation-audit**: manual/scheduled, non-blocking initially, uploads mutation report and surviving-mutant list;
- **benchmark-contract**: validates manifest/schema and runs a tiny synthetic truth fixture in CI; full external benchmark execution remains manual/local or in a dedicated runner with externally supplied data.

Every audit run records commit, Python version, tool versions, configuration hashes and relevant reference/policy identities.

## Versioning

Current package version is `0.8.1`. This change introduces a meaningful new validation/audit subsystem but does not alter the biological result schema or analysis algorithms. Recommended version: **`0.9.0`** (minor-version increment).

The version bump must be applied consistently to all declared/versioned surfaces and checked by the existing version-consistency guard. Changelog wording must state that `0.9.0` is an engineering/audit expansion, not an analytically validated clinical release.

## Testing strategy

Implementation is test-driven:

1. Add failing tests for audit manifest validation and core invariants.
2. Implement the minimum audit contracts/runners to pass them.
3. Add Hypothesis property tests and curated metamorphic regression fixtures.
4. Add a deliberately tiny synthetic benchmark with known TP/FP/FN behavior and assert exact metrics.
5. Configure mutation testing and capture a baseline without setting a blocking threshold.
6. Run full existing unit/pytest, lint, format, mypy, schema/version checks and relevant Snakemake dry-runs.
7. Verify no production algorithm output changes on existing locked fixtures unless an independently identified defect requires a separately reviewed fix.

## Acceptance criteria

The change is complete when:

- property/metamorphic audit tests are deterministic and pass on supported Python versions;
- benchmark manifest rejects build/checksum/identity mismatches and produces exact metrics on the synthetic fixture;
- mutation testing can be reproduced and emits machine-readable or stable text artifacts listing surviving mutants;
- existing CI remains green;
- version `0.9.0` is consistent across package metadata/documentation/version surfaces;
- documentation clearly separates software audit evidence from analytical/clinical validation;
- no large external benchmark datasets are committed to Git.
