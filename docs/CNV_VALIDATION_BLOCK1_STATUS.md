# CNV validation Block 1 status

Research Use Only.

Current phase: implementation complete; final same-head verification pending.

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

The matrix retains `retain_all_evidence=true` as a non-negotiable study rule. Block 1 does not yet create full-evidence rows; full caller evidence, contribution/exclusion states and aggregate-to-source traceability are Block 2.

No caller is selected as a production default. No biological or clinical performance claim is made. The package version remains unchanged. The draft PR must remain unmerged until the final same-head CI/review evidence is recorded.
