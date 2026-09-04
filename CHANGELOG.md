# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

### Added

- A modified-base (methylation) lane. `modkit pileup` is wired into the canonical run graph as
  the version-locked `methylation` stage, aggregating `MM`/`ML` calls into per-region 5mC
  fractions over either canonical chromosomes or the locked target design. It runs only when the
  manifest requests the module, is deselectable like any other component, and its normalized
  report is a validated artifact in the run envelope with its own module outcome, tool record and
  release-bundle checksum. See [`docs/METHYLATION_LANE.md`](docs/METHYLATION_LANE.md).
- A deterministic in-silico tumour dilution series. `plan_dilution_series` lays out the whole
  titration as reviewable data — per-level read budgets, derived seeds, exact subsample arguments
  — without touching a BAM; `execute_dilution_series` materialises it with version-locked
  samtools and verifies every level against the fraction it claims.
- A technical limit-of-detection evaluation over the benchmark reports of a series, reusing the
  `tumor_fraction`/`replicate` strata the benchmark cases already carry. It reports per-level
  detection rates, whether the limit is bracketed by an observed failing level, and refuses to
  report a low level that passed while a higher one failed.
- CLI: `call-methylation`, `dilution-plan`, `dilution-mix` and `lod`; `ontseq run` and
  `ontseq preflight` accept `--methylation-policy` and `--modkit`.
- Technical policies `configs/methylation/modkit.technical.yaml`,
  `configs/benchmark/dilution_series.technical.yaml` and `configs/benchmark/lod.technical.yaml`,
  plus seven exported schemas for the new contracts.

### Fixed

- **The SV consensus never reached any reviewer artifact.** The CNV lane arrives by
  registration and replaces the assemble stage; its copy read intake, QC and the raw Sniffles
  report only, so `assemble_aligned_bam_mvp` fell back to `sniffles_report.events`. Since the
  consensus layer landed in `6ce5d8f`, every `ontseq run` — the extension is registered for
  `run`, `serve` and `watch` — produced result JSON, HTML and XLSX built from raw Sniffles
  calls, silently dropping cuteSV entirely along with consensus merging, gene and cytoband
  annotation, repeat/blacklist/mappability context, Adaptive Sampling observability and AML
  rearrangement prioritization. The consensus artifact was still written to the envelope, so
  nothing looked wrong. Both implementations now read through one shared
  `load_assemble_inputs`, driven by a single `ASSEMBLE_SOURCE_ARTIFACTS` list.
- **Assemble could resume a stale result.** `stage_signature` hashes the artifacts of a
  stage's *declared* dependencies and assemble depends only on QC, so the SV, consensus and
  methylation reports it reads were outside its resume signature — the runner's copy
  fingerprinted none of them and the extension's copy missed the consensus. Both now
  fingerprint every artifact in `ASSEMBLE_SOURCE_ARTIFACTS`.
- The structural-variant stage now runs only when the manifest requests the `sv` module, the
  way CNV and methylation already do. It previously gated on *which policies were supplied*
  rather than on *what the run asked for*, so a manifest declaring `modules: [qc, cnv, report]`
  still drove an SV attempt — and a CNV-only run died with "cuteSV requires --reference-fasta"
  on a reference it had no reason to supply. A run that did not request the module now records
  `NOT_RUN` with the reason naming it a scope statement; a run that *did* request it with no
  caller policy still fails closed.

### Changed

- Preflight answers the methylation lane's preconditions before an envelope exists: policy
  present, reference FASTA present when the pileup is CpG-restricted, target BED present when
  aggregation is over the design, and a warning that MM/ML tag presence can only be established
  by reading the BAM.
- Preflight applies the same scope rule to the SV callers: a run that does not request the
  module is neither told about a missing sniffles or cuteSV nor held to their version locks,
  so the tool section keeps meaning something.
- Stage skip vocabulary is now consistent and distinct: `applicable: false` means the assay has
  nothing for the stage to measure, `requested: false` means the manifest did not ask for it.
  `docs/PIPELINE_EXECUTION.md` documents the three gates side by side.
- The unverified-adapter warning is no longer raised for a methylation stage the run never
  requested, so the line keeps meaning something.

### Validation impact

- Both lanes are new evidence surfaces and neither is validated. The modkit adapter has **never**
  been executed against the real binary here or in CI; it is declared `unverified_adapter` and a
  run completing that stage is reported as such. Fail-closed behaviour is deliberate and load
  bearing: a BAM without `MM` tags fails the stage rather than producing an empty pileup that
  would read as unmethylated DNA, a region with no site above the coverage floor reports `null`
  rather than `0.0`, and the modkit confidence threshold is pinned in policy rather than
  estimated from the sample.
- A detection limit from an in-silico series characterises software behaviour on one pair of
  BAMs. It reproduces read-fraction effects and nothing about library preparation, input mass or
  capture behaviour at low tumour content, its replicates are not independent specimens, and an
  unbracketed limit is reported as a bound rather than a limit. No number from this lane is an
  analytical or clinical sensitivity.
- No existing lane's output changes. Assembly gains an optional methylation module outcome;
  results without the lane are byte-identical apart from that absence.
- The consensus fix changes what every run with CNV registered reports: results now carry
  the consolidated, annotated, prioritized events instead of raw Sniffles calls. This is a
  correction, not a new capability — the evidence was already being computed and written to
  the envelope, only not read back. Reviewers of runs produced since `6ce5d8f` should know
  their reports understated the SV layer: single-caller, unannotated, without observability
  or AML relevance. Envelopes can be re-assembled by re-running with `--force`.
- The SV gating fix does change behaviour for one case, deliberately: a run whose manifest
  omits `sv` but whose configuration supplied caller policies previously produced structural
  variant evidence and now records `NOT_RUN`. That evidence was outside the declared analysis
  scope; a run asking for CNV was never asking for SV. The change is visible rather than
  silent — the module outcome, run report, HTML and XLSX all carry the reason — and a manifest
  that lists `sv` behaves exactly as before. Runs already in an envelope are unaffected: the
  stage signature change re-runs the stage rather than reinterpreting an existing artifact.

## 0.4.1 - 2026-08-27

### Added

- A productive cuteSV 2.1.3 adapter now runs beside Sniffles2 with version-locked parameters,
  atomic VCF finalization and caller-specific normalized evidence.
- A versioned consensus layer canonicalizes reversed BND/TRA representations and clusters
  compatible DEL/DUP/INV/INS/translocation records by type-aware breakpoint, overlap, length and
  available orientation rules. Source caller IDs and evidence remain traceable.
- Build- and checksum-locked interval annotation now supplies separate breakpoint genes,
  cytobands, nearest-gene distances and repeat/tandem-repeat/segmental-duplication/blacklist/
  mappability/centromere/telomere context flags. A preparation script records both original and
  normalized resource hashes.
- Adaptive Sampling SVs now carry explicit breakpoint depths and observability states. A small,
  versioned AML rearrangement-pattern resource prioritizes recurrent patterns without asserting a
  confirmed fusion.
- The technical evidence score is fully represented by a versioned policy. HTML and XLSX show a
  filterable high/moderate review queue, Gene A/B, caller support/consensus, coverage, artifact
  context, AML relevance and validation status while retaining the complete technical table.
- Synthetic regression coverage includes cuteSV-only normalization/execution, Sniffles-only
  behavior, two-caller consensus, canonical reversed BND orientation, within-caller clustering,
  separated nearby events, annotation build refusal, repeat context, Adaptive Sampling partial and
  insufficient coverage, AML patterns and the permanent non-reportable boundary.

### Fixed

- Corrected the productive cuteSV positional argument order to `BAM REF VCF WORKDIR` and exercised
  it through the real-tool CI smoke test.
- Resolved strict type-check findings in deterministic SV ordering and provenance resource locks.

### Validation impact

- This change can alter technical SV clustering, priority tiers and report order. All new matching,
  observability and scoring cut-offs are explicitly `technical_defaults_only`; none is a validated
  sensitivity, specificity, PPV, LoD or clinical threshold.
- Every automated event remains `reportable=false`. Multi-caller support and AML pattern matches do
  not set `analytically_validated`, assert a fusion or establish a somatic origin. Independent
  public technical benchmarks and orthogonally characterized AML specimens remain required before
  clinical use.

## 0.4.0 - 2026-08-27

### Changed

- Research-only/validation status no longer participates in QC verdict calculation.
  Configured measurable gates alone determine `PASS` or `FAIL`; absence of configured gates
  remains `WARN`.
- The research-use disclaimer remains visible as report text. It does not change QC status.
  Non-reportable findings and the separate expert-review workflow remain unchanged.

### Fixed

- Tool output is decoded as UTF-8 with replacement instead of through `text=True`, which
  used the platform locale. Under `LC_ALL=C` — the default in many containers, cron
  environments and freshly installed WSL distributions — that resolved to ASCII, and a
  single non-ASCII byte in a tool's own banner raised `UnicodeDecodeError` from inside the
  runner. A version probe could therefore fail a stage on one machine and pass on another
  from identical inputs. `run_to_file` already decoded this way; both paths now agree.
- The local service refused to escape its allowed roots for BAM paths but built the review
  envelope path by joining two raw URL segments onto the output directory, so
  `POST /api/review/../<sample>` reached an envelope outside `--output-dir` and could append
  a sign-off there. Both identifiers are now matched against the manifest's `run_id` /
  `sample_id` contract and the resolved path is confirmed to be inside the output directory.
- The service decided "one analysis at a time" by asking `Jobs.running()` and registering the
  job afterwards. Between those two steps it read the reference lock and resolved the BAM and
  its index, so simultaneous requests all passed the gate and started concurrent pipelines,
  each sized for the whole machine. The envelope lock does not catch this: two runs with
  different run ids take different envelopes. The slot is now claimed atomically, and reusing
  a run id already in the table is refused rather than dropping the earlier job's record.
- A dangling symlink or an unreadable entry in a browsed directory made `stat()` raise inside
  the request, which ended the response without a body and showed the operator a network
  error for a directory whose other files were usable. Such entries are now skipped, counted
  and reported to the page as `unreadable_entries`, so a BAM that is present but invisible is
  distinguishable from one that is absent.
- BAM name lookup now builds candidates from directory entries and resolves each candidate
  inside its declared allowed root before reading metadata. A BAM symlink that points outside
  the service boundary is omitted rather than exposing its target's size or presenting it as
  selectable input.
- `GET /` handed out the session token without the loopback `Host` check that every `/api/`
  route applies, so a page on an attacker-controlled name resolved to `127.0.0.1` was
  same-origin with the service and could read the token out of the response. The API check
  still refused the stolen token; the page route now applies the same check.

### Added

- `tests/test_cli_surface.py`. Both parsers, the dispatcher and the command overview had no
  test coverage at all. The tests tie the three descriptions of the command set together —
  the overview listing, the dispatch set and the parsers themselves — and check that a bad
  invocation leaves through `SystemExit` rather than a traceback.
- `scripts/check_version_consistency.py`, run by `make versions` and by CI, refuses a tree
  whose declared versions disagree. Five files state what a release is — package metadata,
  the package, the citation record, the desktop project and the README's status section —
  and nothing tied them together. They had drifted: the README described the desktop on
  `main` as the v0.2.1 engineering path while every other declaration said 0.3.4.

### Fixed (continued)

- Preflight was blind to everything adaptive sampling needs. `mosdepth` had no entry in
  `TOOLS_BY_STAGE`, so no run of any input kind checked for it, and nothing checked the
  target BED or the target-coverage policy. All three fail the target-coverage stage closed,
  and a `FAILED` target-coverage stage fails the whole run — but only after the envelope
  exists and the lock is taken, which is exactly what preflight was written to prevent.
  Preflight now probes Mosdepth with the adapter's own version parser, enforces its policy
  lock, and parses the declared target BED. A tool is escalated from warning to blocking when
  the stage needing it is one *this* run cannot do without, so a missing Mosdepth refuses an
  adaptive-sampling run and only warns on an lcWGS one, which correctly records targets as
  out of scope.
- `ontseq preflight` accepted `--target-coverage-policy` and `--components` and used
  neither. The policy was dropped on the way into `PreflightRequest`, and the component
  selection was resolved only for `run` and `serve` — so preflight could check the default
  policies while `ontseq run` would use the ones a selection names. Both are honoured now,
  for the alignment, basecall, SV and target-coverage policies alike.
- `ontseq` never listed `validate-reference` in its command overview, so a working command
  was reachable only by already knowing its name — the discoverability problem the overview
  exists to solve. It is listed now, and a test fails if any command is missing from the
  overview or from the dispatch set.

### Changed

- `target_coverage.mosdepth_version` is public, so preflight probes the binary with the same
  parser the stage uses rather than a second implementation of it.
- `docs/PIPELINE_EXECUTION.md` described target coverage as "not implemented" with no tool,
  and grouped it with CNV as "not wired in". It has been wired into the canonical runner and
  verified against real Mosdepth since 0.3.4; the stage table, the scope note and the
  preflight section now say what the code does.
- Core, Desktop, WSL runtime and Windows bundle identity advance together to 0.4.0. The
  consistency gate now covers the executable Desktop label, runtime installer, Desktop CI
  artifact version and current operator documentation in addition to the original five
  release declarations.

### Validation impact

The QC verdict update is a status-semantics change only. It does not alter QC measurements,
configured thresholds, CNV/SV evidence, reportability, ISCN output or clinical-release
safeguards.

The boundary fixes do not change an analytical threshold, caller parameter, normalization
rule or result contract, and no stage produces a different outcome from the same inputs on a
machine where the previous version ran to completion. The decoding fix removes a
locale-dependent failure mode, so a run that previously aborted mid-stage on a non-UTF-8
locale now proceeds; it does not change what that run reports. Nothing here makes any output
more validated than it was.

The 0.4.0 version advance changes release and runtime identity only; it does not itself add
analytical or clinical validation.

## 0.3.5 - 2026-08-25

### Added

- Pinned `QDNAseq.hg38` runtime resources for the GRCh38 CNV lane alongside the existing
  GRCh37/hg19 resources.
- A canonical real-tool GRCh38 CI gate that exercises 100/500/1000-kbp QDNAseq + ACE,
  verifies the expected synthetic chromosome 7 loss and chromosome 8 gain, checks the
  generated HTML/XLSX/JSON artifacts and proves content-addressed resume.

### Changed

- Desktop and Python Core version declarations were advanced together to 0.3.5.
- GRCh38 Desktop runs may request the integrated QDNAseq + ACE CNV stage now that the exact
  hg38 annotation resource is pinned and exercised in real-tool CI.

### Validation impact

The new GRCh38 lane has deterministic engineering verification only. It has not undergone
cohort-level analytical or clinical validation. The ACE cellularity/ploidy fit, CNV warning
boundaries and whole-chromosome classification fraction remain research-only engineering
parameters and must not be interpreted as validated clinical cut-offs.

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
