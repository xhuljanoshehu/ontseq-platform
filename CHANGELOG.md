# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

### Added

- `ontseq_platform.cnv` subsystem: a segmentation-independent CNV comparison core, an
  explicit observability mask, multi-source truth representation, deterministic
  simulation, a baseline read-depth caller and stratified aggregation.
- Base-level per-state confusion scoring over an exact genome breakpoint partition, so
  copy-number agreement no longer depends on how either side segmented.
- Closed vocabulary of exclusion reasons with per-reason base accounting, keeping
  biological negativity, no-call and technical failure distinguishable.
- `background_state` and `resolution_bp` on truth and call sets, encoding open-world
  versus closed-world semantics and each source's detection limit.
- Per-segment breakpoint uncertainty, with breakpoint metrics withheld when the truth
  source cannot support them.
- ISCN karyotype to copy-number conversion against a versioned cytoband resource, with
  unsupported constructs recorded and surfaced rather than dropped.
- Deterministic gamma-Poisson dilution and coverage simulation, plus empirical and
  model-based limit-of-detection estimation with explicit withholding.
- `ontseq-baseline-readdepth`, a non-reportable control caller.
- Header-driven segment-table adapters with declarative column mappings for the generic
  IGV `SEG` format and ichorCNA.
- Paired method comparison via McNemar's exact test on truth events assessable under
  both methods, with the p-value withheld when no discordant pair exists.
- A genome partition that reconciles exactly, enforced in the core and re-validated by
  the contract, so every metric carries an auditable denominator.
- CLI: `cnv-evaluate`, `cnv-aggregate`, `cnv-compare-methods`, `cnv-karyotype-truth`,
  `cnv-demo-benchmark`.
- JSON Schemas for the CNV truth set, call set, benchmark case, evaluation report,
  aggregate report and cytoband table.
- `docs/CNV_BENCHMARKING.md` and ADR-008 through ADR-012.

### Changed

- `scripts/export_schemas.py` also exports the CNV contracts.

### Validation impact

No CNV method is selected, promoted or validated by this change. Every figure the new
subsystem produces is an engineering measurement against a declared truth set on
synthetic data. Comparison thresholds, the coverage floor and the exclusion tracks are
engineering parameters that must be locked before comparative results are inspected;
none of them is a validated adequacy or reportability threshold. The existing SV
benchmark contract in `benchmark.py` is unchanged.

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
