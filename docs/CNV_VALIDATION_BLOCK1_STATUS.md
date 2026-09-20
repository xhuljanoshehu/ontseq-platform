# CNV validation Block 1 status

Research Use Only.

Status: PASS.

Final verified Block 1 head: `afb721396f6f880290861dfd3abeaf81c794beef`.

Implemented on the isolated `feat/cnv-validation-program` branch:

- prospective CNV validation matrix and stratification contracts;
- explicit caller, genome-build and CNV data-basis locks;
- explicit coverage, tumour/blast-fraction and event-size strata;
- fail-closed requirement for registered cutpoints and QDNAseq+ACE bin sizes;
- orthogonal truth-source, assessability-mask and optional negative-universe contracts;
- explicit specimen coverage and tumour/blast-fraction metadata without imputation;
- repeat-group semantics for within-run and between-run reproducibility;
- fail-closed cohort eligibility checks;
- canonical SHA-256 preregistration of matrix and cohort before outcome review;
- specificity eligibility only when a prospectively registered negative denominator exists;
- versioned JSON schemas for matrix, cohort and registration;
- synthetic contract tests only; no real biological data are committed.

The matrix retains `retain_all_evidence=true` as a non-negotiable study rule. Block 2 extends that rule with full caller evidence, derived normalized events and aggregate-to-source traceability without reopening the verified Block 1 contracts.

Final same-head Block 1 verification on `afb721396f6f880290861dfd3abeaf81c794beef`:

- CI #732 — SUCCESS on Python 3.11, 3.12 and 3.13; 1716 passed, 6 skipped, 407 subtests; Ruff, mypy, schema freshness, repository safety, version agreement, wheel-resource checks, Snakemake DAGs and synthetic report generation passed;
- MARLIN runtime smoke #118 — SUCCESS;
- Desktop bundle contract #218 — SUCCESS;
- Desktop CI #492 — SUCCESS, including relocatable Linux runtime, stock-Ubuntu system smoke, WPF build, BAM-index/reference-identity tests, self-contained `win-x64` publish and first-run bundle verification.

No caller is selected as a production default. No biological or clinical performance claim is made. Package version remains `0.8.2`. Block 1 did not merge the draft PR.
