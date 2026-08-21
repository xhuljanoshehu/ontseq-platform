# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

### Added

- Canonical GRCh37/GRCh38 reference validation for the Desktop setup boundary.
- FAI-dictionary-addressed Desktop reference identity and full BAM/index/reference resume
  signatures.
- Actionable sequence-dictionary mismatch details retained in failed run artifacts and UI status.
- Explicit verification of the declared BAM index, without falling back to an unindexed BAM scan.
- Version-locked Mosdepth target-coverage adapter for Adaptive Sampling aligned-BAM runs.
- Strict target-BED and Mosdepth output normalization with exact interval reconciliation,
  per-region mean depth and configurable threshold fractions.
- Dedicated technical policy for target coverage and a synthetic real-tool CI path.
- Documentation that separates the unbuffered analysis ROI BED from operational Adaptive
  Sampling selection regions.

### Changed

- Bumped the execution core to 0.3.1 and embedded the exact Git commit in packed Desktop
  runtimes for run provenance.

### Validation impact

The aligned-BAM integrity gate remains strict. The Desktop now refuses partial, mixed-style or
wrong-build FAI dictionaries when they are configured as a full named GRCh37/GRCh38 reference,
and a changed BAM index or reference lock invalidates a prior intake resume. This is a dictionary
compatibility check, not proof that two FASTA files have identical bases. Region-extracted BAMs
still require the full sequence dictionary of the reference used for alignment; this change does
not treat sparse regional data as adequate lcWGS evidence.

The Adaptive Sampling path gains new descriptive coverage output. The default `1x`, `10x`, `20x`
and `30x` thresholds are engineering bins only and do not define assay adequacy, CNV/fusion
reportability, biological negativity or a clinical no-call. The target design and thresholds require
pre-specified validation on the locked local panel before any promotion.

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
