# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

### Added

- Version-locked Mosdepth target-coverage adapter for Adaptive Sampling aligned-BAM runs.
- Strict target-BED and Mosdepth output normalization with exact interval reconciliation,
  per-region mean depth and configurable threshold fractions.
- Dedicated technical policy for target coverage and a synthetic real-tool CI path.
- Documentation that separates the unbuffered analysis ROI BED from operational Adaptive
  Sampling selection regions.
- End-to-end runtime wiring for the pre-existing `target_coverage` stage on Adaptive Sampling
  runs. The stage fingerprints both the aligned BAM and target BED for content-addressed resume,
  keeps raw Mosdepth files under non-exportable `work/`, and exports only the normalized
  observability JSON under `qc/`.
- Regression tests that keep lcWGS target coverage explicitly `NOT_RUN`, fail early on a Mosdepth
  version mismatch, and verify the normalized observability artifact does not expose source paths.

### Validation impact

The Adaptive Sampling path gains new descriptive coverage output. The default `1x`, `10x`, `20x`
and `30x` thresholds are engineering bins only and do not define assay adequacy, CNV/fusion
reportability, biological negativity or a clinical no-call. The target design and thresholds require
pre-specified validation on the locked local panel before any promotion.

The runtime now executes this descriptive observability stage automatically for Adaptive Sampling
runs and records failure explicitly if the target design, locked Mosdepth version, aligned BAM, or
normalization contract is unavailable. This improves traceability but does not alter CNV, SV,
fusion, ISCN, confidence or reportability logic. The normalized target-coverage sidecar is not yet a
field in `PipelineResult`; integration into the reviewer report remains a separate schema-migration
step so presentation cannot silently invent an observability conclusion.

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
