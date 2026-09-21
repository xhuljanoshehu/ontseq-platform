# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

- Restrict the legacy GATK research run directory to mode `0700` on POSIX and launch native child processes with umask `077`, preventing local disclosure of BAM-header/sample and variant artifacts. No caller parameters, thresholds, result semantics, or clinical-reportability rules change.

- Add opt-in GATK 4.6.2.0 research adapters, separate candidate/config/result schemas,
  tumor-only and paired modes, native artifact provenance and synthetic contract tests.
- Repair strict typing and exact-version checks; reject suffixed or extended GATK builds.
  Include both adapter entry points in focused CI without weakening full-project gates.
- Validation impact: adapters remain separately invoked, non-reportable and not
  real-GATK-qualified or analytically validated. No canonical runner, CNV/SV policy,
  reference bundle, package version or automatic clinical release changes.

- Integrate the owner-selected complete Befund design as `befund-interactive-v1`:
  offline React report, real run-bound evidence, original CNV plots, model selection,
  synchronized chromosome/dilution views, explicit hypothetical VAF anchoring,
  module states, methylation, ISCN and complete technical evidence. Package the
  reproducible report bundle in wheels and stamp the report resume signature.
- Add an authenticated current-layout export for existing runs; retain archived
  reports and analysis artifacts unchanged. Validate evidence identity, paths and
  recorded checksums. Preserve the static report as a no-JavaScript fallback.
- Validation impact: additional exploratory presentation calculations only; no new
  caller, threshold, classification or release behavior. No patient example from
  the supplied design is imported. Missing results never become measured zero.

- Local `core-or-sample-coverage-v1` fix: resolve both core and installed-runtime
  coverage filenames for SV observability and HTML/Excel/CLI reporting. Reject
  mismatched samples/builds and conflicting duplicates. Track coverage artifacts
  in adaptive-sampling SV/report resume signatures and fail before running SV callers
  when coverage is absent. Existing sample-named reports remain readable.
- Validation impact: existing per-target measurements can now reach downstream SV
  observability and reports. No depth threshold, caller setting, methylation policy
  or clinical release rule changes; candidate events remain research-only.

- Local presentation revision `befund-v8-local-1`: use the owner-selected report
  composition with a dark header, event/status cards and CNV evidence near the top.
  Preserve normalized results, release boundaries, complete evidence tables and offline use.
  Add responsive and print styles; show absent module records explicitly.
- Validation impact of this local layout revision: presentation only; no caller,
  classification, threshold, copy-number model or release-rule changes. The reference
  report's sample-specific interpretations and simulations are not imported.

- Add a pinned Spectre 0.2.1 `depth_only` CNV adapter that consumes only indexed
  Mosdepth coverage and the locked reference. The shell-free runtime rejects version,
  build, source and sample mismatches; never accepts Sniffles/SNFJ, SNV, population or
  cancer-mode inputs; and promotes output only after fail-closed normalization.
- Preserve every caller-native Spectre artifact with SHA-256 provenance, including empty
  native VCF output as an explicit technical `NO_CALL`. Pin the upstream wheel hash and
  source commit, and exercise both deterministic DEL/DUP calls and an empty result with
  the real executable on generated synthetic depth data.
- Validation impact: the Spectre lane is technically interoperability-qualified only.
  It is not analytically or clinically validated, does not establish a biological
  negative, remains Research Use Only and cannot make events reportable. No caller
  selection, voting, winner, threshold-sharing or automatic ISCN logic is introduced.
- Prevent CNV normalization from aborting when a segment covers all assessable bins but
  remains below the versioned chromosome-classification fraction. Such events remain
  duplications/deletions without whole-chromosome confirmation; thresholds are unchanged.
  Record `classified-event-span-v1` in tool provenance and the CNV resume signature.
- Preserve the interval end in methylation heatmap identity so nested targets sharing a
  label and start remain distinct columns and no longer abort HTML reporting.
- Synchronize the Live Workspace package and generated HTML with Core 0.8.2 so the
  reproducible web build no longer fails on an inherited 0.7.1 version declaration.

- Confirm a whole-chromosome ISCN span against the assessable QDNAseq bin extent instead of the
  raw contig ends. QDNAseq drops telomeric, blacklisted and residual-filtered bins, so the former
  exact zero-to-contig-end rule was unreachable through this lane: every real whole-chromosome
  candidate was classified as a gain or loss and then omitted as an unsupported construct. The
  basis is versioned in the CNV policy as `whole_chromosome_span_basis`, shipped as
  `assessable_bin_extent`, and is reported in the ISCN policy parameters. `exact_contig` keeps the
  previous behaviour, and a bin table without the exporter's `use` flag falls back to it.
- Record the assessed region on the event (`assessable_span_start`, `assessable_span_end`) so the
  model still verifies a confirmation instead of trusting the caller, and state in the event notes
  which basis confirmed or suppressed the span. Filtered regions are never asserted to be
  unchanged sequence.
- No caller threshold, classification fraction, reference bundle or reportability rule changes.
  Whole-chromosome proposals remain research-use-only and require cytogenetic expert review.

## 0.8.2 — 2026-09-12