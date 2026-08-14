# Inputs needed for the validated implementation

The repository foundation is executable without external literature. Moving from the
synthetic contract demo to a scientifically and clinically defensible pipeline requires
controlled source material and validation data.

## Highest-priority inputs

| Input | Why it is needed | Safe form for collaboration |
| --- | --- | --- |
| Approved thesis code or rule archive | Reproduce exact commands and parameters without guessing | Private repository; no sample data |
| Exact tool and reference versions | Lock reproducibility and resolve version-dependent fields | Text manifest with checksums |
| Institution-approved ISCN 2024 access | Build and test syntax beyond the deliberately limited subset | Authorized local test cases; do not copy the standard publicly |
| Adaptive-sampling target BEDs | Validate target coverage and fusion search space | Versioned, non-patient assay design |
| Orthogonally characterized truth cohort | Estimate sensitivity, specificity, LoD and no-call rate | Approved on-premises data access only |
| Expected HTML/Excel/SOP examples | Match reviewer workflow and required audit fields | De-identified or synthetic examples |
| Intended-use and release SOP | Define what the software may report and who may sign it | Controlled internal document |

## Literature workstream

For every production adapter, record the primary tool publication, official manual,
software version, parameter rationale, benchmark evidence and known limitations. Public
literature may be linked in Git; copyrighted standards and patient-derived material must
remain in approved controlled systems.

The first implementation priority is the **aligned-BAM MVP**. POD5 basecalling, optional
small-variant calling, methylation and app/LIMS integration are separate validation domains
and follow only after the core CNV/SV/fusion path is reproducible.
