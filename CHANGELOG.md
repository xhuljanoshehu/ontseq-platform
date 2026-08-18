# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

### Added

- Research-only DNA fusion evidence contracts with build-locked local gene annotation,
  privacy-safe BND adjacency descriptors and explicit breakpoint observability.
- Synthetic fusion software benchmark fixtures and executable VCF-to-review test paths that do
  not contain patient-derived genomic material.
- Exact duplicate/reciprocal breakpoint redundancy detection that preserves every source event
  and never auto-deduplicates candidate evidence.
- Adversarial fusion tests covering filtered and low-support records, malformed BND ALT fallback,
  reciprocal breakends and exact duplicate breakpoints.
- Typed fusion reviewer contract that keeps candidate evidence, `NO_CALL`, `NOT_RUN` and `FAILED`
  distinct and exposes review-required, research-only, non-reportable candidate summaries.
- Reviewer privacy tests that prohibit raw VCF ALT, inserted sequence, read names and source file
  paths from the serialized reviewer contract.
- Typed multi-caller SV concordance contract for independent caller comparison, with explicit
  exact-match, tolerance-based near-match, topology-conflict and unmatched evidence states.
- Synthetic Sniffles2/cuteSV concordance tests covering exact and near matches, out-of-tolerance
  calls, reciprocal breakpoint ordering, event-type conflicts and privacy-safe serialization.
- Conservative cuteSV v2.1.4 VCF normalization contract for DEL/DUP/INV/INS/BND/TRA evidence,
  including explicit filter/support accounting, version lock, `NO_CALL` semantics and a typed
  bridge into the multi-caller concordance layer.
- Synthetic cuteSV tests covering 1-based to 0-based coordinate normalization, BND mate parsing,
  support filtering, version mismatch, malformed breakends and exclusion of VCF IDs/read names.
- Shell-free cuteSV v2.1.4 execution adapter with explicit ONT software parameters, temporary
  work-directory cleanup, disabled read-ID reporting, suppressed inserted-sequence output and
  fail-closed FASTA-index/BAM-header contig-signature compatibility checking.
- Synthetic cuteSV execution tests covering CLI provenance, version mismatch, reference mismatch,
  partial-output cleanup, overwrite refusal and path-free normalized serialization.
- Direct Sniffles2-to-cuteSV concordance bridge that requires matching sample/build metadata and
  preserves caller-level `NO_CALL` as non-negative evidence status.

### Validation impact

The fusion branch now produces reviewer-facing structured projections of research-only DNA
rearrangement evidence, can compare normalized observations from independent SV callers, can
normalize a version-locked cuteSV VCF into the common candidate evidence model, and can execute
cuteSV locally against an aligned BAM plus an explicitly supplied reference FASTA. These
projections and executions can change how candidate evidence, observability, redundancy, module
status and caller concordance are presented to a reviewer, but they do not establish biological
truth, an expressed or functional fusion, or clinical reportability. `NO_CALL` explicitly remains
non-negative, known gene pairs remain annotation evidence only, genomic BND orientation does not
establish transcript direction, and duplicate/reciprocal records are preserved rather than
collapsed. Multi-caller exact or near agreement is labelled software evidence only; the breakpoint
tolerance is explicit configuration with `clinically_validated=false`, and conflicting caller
topology is preserved for review rather than coerced into agreement. The cuteSV adapter treats
BND/TRA records as DNA translocation evidence only, does not propagate RNAMES/read IDs or raw ALT
sequence, and locks software version 2.1.4 while its current thresholds remain
`technical_defaults_only` and `clinically_validated=false`. The local cuteSV runner uses explicit
ONT-oriented software parameters, never enables `--report_readid`, enables `--ignore_sequence`,
uses a temporary work directory, refuses output overwrite, and removes partial VCF output after a
non-zero caller exit. Reference compatibility is checked by matching the FASTA-index contig
name/length signature to the aligned-BAM header; this does not prove full reference sequence
identity and is not a substitute for a future sequence-level reference checksum lock.
Assay-specific analytical validation and authorized human review remain mandatory before any
clinical use, real-sample promotion, or connection to final HTML/XLSX/ISCN release logic.

## 0.3.0 - 2026-08-14

### Added

- Typed Sniffles2 policy and call-report contracts with generated JSON Schemas.
- Conservative Sniffles2 v2.8.0 adapter using shell-free execution, symbolic VCF and PASS-only
  normalization with filter-reason accounting.
- Defensive DEL/DUP/INV/INS/BND normalization with counted rejection reasons and `NO_CALL`.
- Runtime-generated synthetic long-read BAM fixture and real samtools/Cramino/Sniffles2 smoke test.
- CI job with a pinned Bioconda environment and selected raw-data-free reviewer artifacts.
- Candidate SV evidence in the common JSON, HTML and Excel result path.

### Changed

- Extended evidence provenance with supporting-read strands, coverage context and mean alignment
  NM.
- Connected the aligned-BAM Snakemake DAG to the candidate-only Sniffles2 adapter.
- Bumped the research software foundation to version 0.3.0.

### Validation impact

The aligned-BAM path can now change biological candidate output by executing Sniffles2. Every
normalized event is forced to `unclassified` and `reportable: false`; BND is not treated as a
fusion, and ISCN remains `NOT_RUN`. The included thresholds are engineering defaults only. This
change requires locked HG002/HG008 and intended-use AML validation before any promotion.

## 0.2.0 - 2026-08-14

### Changed

- Reframed the architecture as an independent, evidence-led implementation.
- Reclassified the thesis as historical context rather than a technical specification.
- Replaced unvalidated CNV/SV defaults with benchmark-gated candidate sets.
- Added mandatory coverage, tumor/blast-fraction and no-call validation axes.
- Bumped the research software foundation to version 0.2.0.
- HTML and Excel reports now expose explicit per-module execution status.

### Added

- Living literature and benchmark evidence base with explicit applicability and limitations.
- Fail-closed aligned-BAM/BAI/header/reference intake gate using samtools.
- Versioned reference locks generated from FASTA index files.
- Cramino JSON adapter with normalized, path-free descriptive QC.
- Deterministic normalized-event CNV/SV benchmark engine and synthetic fixtures.
- Snakemake entry points for the aligned-BAM MVP and synthetic benchmarks.

### Validation impact

The report schema and workbook gain module-status data. Aligned-BAM intake can now stop a run on
technical incompatibility. No scientific caller has been promoted to a production or clinical
default; CNV/SV/fusion output remains `NOT_RUN` in the MVP.
