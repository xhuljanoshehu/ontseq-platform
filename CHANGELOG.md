# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

The remainder of the withdrawn PR #8, refreshed onto current `main`: the CNV benchmark
subsystem and the ClinVar knowledge layer. The execution core from the same draft landed
separately through #17 through #21 and is not repeated here.

### Added

- `ontseq_platform.cnv` benchmark subsystem: a segmentation-independent CNV comparison
  core, an explicit observability mask, multi-source truth representation, deterministic
  simulation, a baseline read-depth caller and stratified aggregation. It never runs on a
  patient sample; it exists to measure a caller before one is trusted.
- Base-level per-state confusion scoring over an exact genome breakpoint partition, so
  copy-number agreement no longer depends on how either side segmented.
- Closed vocabulary of exclusion reasons with per-reason base accounting, keeping
  biological negativity, no-call and technical failure distinguishable.
- `background_state` and `resolution_bp` on truth and call sets, encoding open-world
  versus closed-world semantics and each source's detection limit.
- `truth_resolution_silent_bases`: a fourth term in the genome partition holding calls
  finer than the truth source can resolve, so they are neither confirmed nor counted as
  false positives. The count and the number of affected calls are stated in a warning,
  because this is the only exclusion in the design that flatters the caller.
- Per-segment breakpoint uncertainty, with breakpoint metrics withheld when the truth
  source cannot support them.
- ISCN karyotype to copy-number conversion against a versioned cytoband resource, with
  unsupported constructs recorded and surfaced rather than dropped.
- Deterministic gamma-Poisson dilution and coverage simulation, plus empirical and
  model-based limit-of-detection estimation with explicit withholding.
- `ontseq-baseline-readdepth`, a non-reportable control caller.
- Header-driven segment-table adapters with declarative column mappings for the generic
  IGV `SEG` format and ichorCNA.
- Paired method comparison via McNemar's exact test on truth events assessable under both
  methods, with the p-value withheld when no discordant pair exists.
- `observed_direction`, `underpowered`, `minimum_attainable_p_value` and `alpha` on the
  paired method comparison, separating the direction the counts lean from the direction a
  test supports.
- `SpecimenClustering` and a specimen-weighted detection rate on the aggregate report,
  plus `discordant_specimens` on the paired comparison, so event-level intervals are never
  shown without the clustering that qualifies them.
- `reports_biological_negative` on a call set, so a method that looked everywhere and
  found nothing is distinguishable from one that could not look.
- A genome partition that reconciles exactly, enforced in the core and re-validated by the
  contract, so every metric carries an auditable denominator.
- CLI: `cnv-evaluate`, `cnv-aggregate`, `cnv-compare-methods`, `cnv-karyotype-truth`,
  `cnv-demo-benchmark`.
- JSON Schemas for the CNV truth set, call set, benchmark case, evaluation report,
  aggregate report and cytoband table.
- `docs/CNV_BENCHMARKING.md` and ADR-008 through ADR-012, plus ADR-018 through ADR-020.
- `ontseq_platform.knowledge`: a ClinVar annotation layer that attaches records to
  copy-number and structural findings **without classifying them**. Each annotation
  carries the assertion verbatim together with the vocabulary it belongs to
  (`acmg_germline`), the record's own origin, the match type and reciprocal overlap,
  NCBI's review status and star rating, and the checksum of the exact weekly release it
  came from. Records whose origin does not match the assay's question are kept and marked
  as secondary findings rather than filtered, because filtering them would be a clinical
  decision made invisibly by code. Nothing in the package sets `reportable` or raises
  `confidence`: deciding that needs somatic criteria (ELN, ICC, a local gene list) this
  repository does not have.
- `ontseq annotate`: applies a locked ClinVar release to a result contract and writes the
  annotations back onto the events, leaving every judgement field untouched.
- `AnalysisSpec.intent` (`somatic` / `germline` / `both`), **without a default**. It
  decides how every knowledge-base assertion is read, and guessing would settle that
  silently. Left unset, scope alignment is reported as unknown rather than assumed.
- `KnowledgeResourceLock` and `EventAnnotation` contracts, plus a `GenomicEvent.annotations`
  list. `EventAnnotation` refuses an empty `caveats` list: without them a database
  classification reads as a finding about the sample, which is the one reading the whole
  design exists to prevent.
- Knowledge-base annotations are rendered in the HTML report under a heading stating that
  they are classifications of database records rather than findings about the sample, and
  in the Excel workbook as sheet `11_Annotations` under a banner in row 1 saying the same.
  The banner is *in* the sheet because a spreadsheet has no prose — a reviewer sorts by the
  assertion column and reads "Pathogenic" beside an event identifier, and nothing else on
  the grid would say what that means. Rows whose record origin does not match the assay's
  question are filled in the warning colour, because colour is read before column ten is.
- `db_records_matched` on the four event sheets. A reviewer who never opens the annotation
  sheet would otherwise not know there was anything there to open. It is a count, named for
  what it counts, so that it cannot be read as a classification of the finding itself.
- `examples/knowledge/synthetic.clinvar.variant_summary.txt`, an invented fixture in NCBI's
  layout containing no real record, with one row for each awkward case the loader must
  handle: another assembly, a variant with no matchable extent, and ClinVar's `-1`
  placeholder for unplaced records.
- ADR-021 and ADR-022.
- CI runs the synthetic CNV benchmark end to end, compares two baseline configurations
  pairwise through the CLI, and converts a synthetic ISCN karyotype into a truth set.

### Changed

- `scripts/export_schemas.py` also exports the CNV contracts.

### Fixed

- A truth source's declared resolution affected only a warning, not the score. Calls below
  it were counted as false positives while the warning beside them said they must not be —
  and the number, not the warning, is what aggregates. Resolution now removes those bases
  from the evaluable genome, and only where the truth asserts its background: resolution
  limits what a source can deny, never what it can affirm.
- `favours` named a winning method from the raw discordant counts, so a 4-0 split reported
  a winner at p=0.125 — a count at which no possible outcome could have been significant.
  A method is named only when the test is significant at a pre-specified alpha.
- A `COMPLETED` call set was required to contain segments, which forced a method that
  looked everywhere and found nothing to report `NO_CALL` — exactly the conflation between
  a biological negative and an inability to look that the vocabulary exists to prevent.
- The adaptive-sampling off-target read population was described as forming a near-uniform
  low-coverage background. That is a hypothesis about the local assay, not a measured
  property, and it is now written as one in the code and both documents.

### Validation impact

No CNV method is selected, promoted or validated by this change. Every figure the new
subsystem produces is an engineering measurement against a declared truth set on synthetic
data. Comparison thresholds, the coverage floor and the exclusion tracks are engineering
parameters that must be locked before comparative results are inspected; none of them is a
validated adequacy or reportability threshold. The existing SV benchmark contract in
`benchmark.py` is unchanged.

The runtime CNV lane that a patient run executes — QDNAseq + ACE, landed in #18 — is not
altered, selected or promoted here. Nothing in this change makes it reportable; that still
needs real cohort validation. What this change adds is the apparatus that could eventually
measure it.

The ClinVar layer attaches records; it classifies nothing. ClinVar speaks the ACMG germline
vocabulary and an AML workup asks a somatic question, so an assertion carried into such a
report is a statement about a database record and not about the specimen. `intent` has no
default for exactly that reason.

## 0.3.4 - 2026-08-22

### Added

- Per-run component selection: a `RunComponents` document names the provider and the exact
  tool version for each stage that runs an external program. The version is checked against
  the tool's own probe before the stage executes, and a mismatch fails that stage naming
  both versions. Deselected stages are recorded as `NOT_RUN` with the selection named.
- `configs/components/default.yaml` and `configs/components/legacy_sniffles_2.4.yaml`, plus
  a Sniffles 2.4.0 policy, so the structural-variant component of the historical pipeline
  can be reproduced for comparison without editing code.
- Target coverage is wired into the canonical runner as `StageId.TARGET_COVERAGE`. An
  adaptive-sampling run without a target-coverage policy now fails closed instead of
  recording the assay's central QC as absent; a non-adaptive run records that the stage does
  not apply, which is a scope statement rather than a coverage result.
- `AssaySpec.target_bed_role`, distinguishing an unbuffered analysis ROI from a buffered
  selection panel, and carried into the target-coverage report and its limitations.
- First non-synthetic target design: a 111-target, 17.03 Mb buffered GRCh38 panel under
  `configs/panels/`, with a lock file, a reproducible generator script and de-identified
  per-target coverage expectations. Provenance and open questions in
  `docs/PANEL_PROVENANCE.md`.
- CI proves per-target coverage end to end inside a run envelope, proves an lcWGS run
  reports the stage as out of scope, and proves a component pinned to an uninstalled version
  fails closed.
- `docs/COMPONENT_SELECTION.md` and `docs/LEGACY_COMPARISON.md`, the latter mapping the
  historical pipeline's outputs onto this one and proposing MOLM13 as the first positive
  control.
- Version-locked Mosdepth target-coverage adapter for Adaptive Sampling aligned-BAM runs.
- Strict target-BED and Mosdepth output normalization with exact interval reconciliation,
  per-region mean depth and configurable threshold fractions.
- Dedicated technical policy for target coverage and a synthetic real-tool CI path.
- Documentation that separates the unbuffered analysis ROI BED from operational Adaptive
  Sampling selection regions.

### Changed

- `ontseq` with no command lists both command groups. Execution commands were previously
  invisible to `--help`.

### Validation impact

Adaptive-sampling runs now produce per-target coverage as part of the canonical run rather
than only through a separate command. The default `1x`, `10x`, `20x` and `30x` thresholds
remain engineering bins only and do not define assay adequacy, CNV or fusion reportability,
biological negativity or a clinical no-call.

The coverage expectations table is a sanity reference derived from historical runs; it is
not an adequacy gate, not a reportability threshold and not a no-call definition. The panel
is `derived_unconfirmed`: it is reconstructed from laboratory coverage tables rather than
copied from the selection file used at sequencing time, and one target label contradicts its
coordinates.

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
