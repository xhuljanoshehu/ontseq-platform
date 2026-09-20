# Unreleased: optional GATK Mutect2 research adapter

Package version remains 0.8.2. Research Use Only. Human Review Required.

Added a separately executable, default-dry-run Mutect2/FilterMutectCalls adapter,
paired/tumor-only modes, optional locked resources and contamination commands,
fail-closed input/output handling, exact GATK 4.6.2.0 version gate, and native-allele
evidence with provenance. All candidates are non-reportable and NOT_VALIDATED.

Validation impact: this is a new optional evidence path, not a change to the existing
CNV/SV planner, thresholds, clinical report schema or Clair3 germline policy. It is
not registered in the canonical CNV/SV runtime and does not appear in the Desktop.
No existing biological output is silently reinterpreted. A new source/runtime and
assay-validation review is required before any analytical promotion.

Authoring verification: 58 new synthetic adapter tests passed locally on Python 3.13.
Full repository checks and real-tool qualification are separate gates; a fake runner
is not a verified GATK binary. No patient data, reference genomes, PoN datasets,
credentials or real analysis outputs are included. See GATK_MUTECT2_ADAPTER.md.
