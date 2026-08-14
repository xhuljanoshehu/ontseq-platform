# Repository guardrails for coding agents

This project processes potentially identifying human genomic data. Automated contributors
must preserve the following boundaries:

1. Never add patient data, genomic run files, credentials, direct identifiers or reports to
   Git. Use synthetic fixtures only.
2. Treat every ISCN string as an expert-reviewable proposal. Do not remove the research-use
   warning or introduce automatic clinical release.
3. Keep assay-specific QC and reportability rules versioned. A caller default is not a
   validated clinical threshold.
4. Any change that can alter biological output needs tests, provenance fields, a changelog
   entry and validation-impact review.
5. Prefer typed adapters and structured contracts over parsing presentation files.
6. Run `make safety`, `make lint`, and `make test` before proposing a change.
7. Do not add a public license, public deployment, cloud upload or external data transfer
   without explicit owner and institutional approval.

The source of truth for scope and limitations is `docs/ARCHITECTURE.md`,
`docs/CLINICAL_VALIDATION.md`, and `docs/THESIS_TRACEABILITY.md`.
