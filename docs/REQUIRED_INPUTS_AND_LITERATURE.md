# Inputs needed for the validated implementation

The repository foundation is executable without external literature. Moving from the
synthetic contract demo to a scientifically and clinically defensible pipeline requires
controlled source material and validation data.

## Highest-priority inputs

| Input | Why it is needed | Safe form for collaboration |
| --- | --- | --- |
| Intended-use statement and release SOP | Define specimen types, reportable event classes, users and who may sign | Controlled internal document |
| Exact assay design | Lock library kit, flow cell, basecaller model, reference build and adaptive-sampling targets | Versioned, non-patient manifest with checksums |
| Institution-approved ISCN 2024 access | Build and test syntax beyond the deliberately limited subset | Authorized local test cases; do not copy the standard publicly |
| Orthogonally characterized AML truth cohort | Estimate accuracy, LoD and no-call rate by event class and tumor/blast fraction | Approved on-premises data access only |
| Public benchmark plan | Separate technical caller correctness from AML intended-use validation | Versioned HG002/HG008 accessions, truth versions and scoring rules |
| Expected HTML/Excel/SOP examples | Match reviewer workflow and required audit fields | De-identified or synthetic examples |
| Institutional compute profile | Size resources and verify runtime, isolation and data residency | Non-sensitive CPU/GPU/RAM/storage description |

## Optional historical material

Approved thesis code, parameter files and benchmark outputs can be used for regression comparison
after ownership and licensing are clarified. They are not required to select the architecture and
must not be treated as ground truth or copied into this repository without authorization.

## Literature workstream

For every production adapter, record the primary tool publication, official manual, software
version, parameter rationale, benchmark evidence, applicability and known limitations. The
living review and current decisions are in `docs/EVIDENCE_BASE.md`. Public literature may be
linked in Git; copyrighted standards and patient-derived material must remain in approved
controlled systems.

The first implementation priority is the **aligned-BAM MVP**. POD5 basecalling, optional
small-variant calling, methylation and app/LIMS integration are separate validation domains
and follow only after the core CNV/SV/fusion path is reproducible.
