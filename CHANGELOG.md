# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

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
- `GET /` handed out the session token without the loopback `Host` check that every `/api/`
  route applies, so a page on an attacker-controlled name resolved to `127.0.0.1` was
  same-origin with the service and could read the token out of the response. The API check
  still refused the stolen token; the page route now applies the same check.

### Added

- `scripts/check_version_consistency.py`, run by `make versions` and by CI, refuses a tree
  whose declared versions disagree. Five files state what a release is — package metadata,
  the package, the citation record, the desktop project and the README's status section —
  and nothing tied them together. They had drifted: the README described the desktop on
  `main` as the v0.2.1 engineering path while every other declaration said 0.3.4.

### Changed

- README section 14 names the desktop version on `main` as 0.3.4, matching
  `desktop/ONTSeq.Desktop/ONTSeq.Desktop.csproj` and the desktop README.

### Validation impact

None. No analytical threshold, caller parameter, normalization rule or result contract
changed, and no stage produces a different outcome from the same inputs on a machine where
the previous version ran to completion. The decoding fix removes a locale-dependent failure
mode, so a run that previously aborted mid-stage on a non-UTF-8 locale now proceeds — it does
not change what that run reports. Nothing here makes any output more validated than it was.

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
