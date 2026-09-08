# Two-source methylation holdout status — 2026-09-05

**Follow-up:** Approved local acquisition and bounded structural intake are complete.
See `METHYLATION_INTAKE_STATUS.md`. The snapshot below predates that intake and the
Nanopolish v2 security update; its code fingerprint is historical. Biological recovery
remains unestablished.

**Biological acceptance: NO_CALL / NOT_REGISTERED. No new biological recovery
experiment has run.** The existing uncommitted work and historical results were
preserved. No genomic input was downloaded or added to Git.

## Current implementation

The [two-source holdout](METHYLATION_HOLDOUT.md) now supports Nanopolish, SAM/BAM
MM/ML and modkit extract-full inputs through the same registered experiment runner.
Two genuine biological source records are split into calibration and test read
pools; those computational partitions are not additional biological specimens.
Source files, protocol metadata, whole-read split seeds, matrix, marker set and
software/runtime content are bound by SHA-256 and checked during execution/replay.
Registration records an operator declaration that test outcomes were unseen;
it is not an externally authenticated preregistration.

Reports contain JSON, per-level CSV and HTML with bias, MAE, RMSE, conditional
interval coverage and unconditional NO_CALL rate, including missing planned levels.
An optional stricter [four-source donor study](METHYLATION_VALIDATION.md) remains
available. The two-source report fixes donor-independent validation to
`NOT_EVALUATED`, including when synthetic recovery passes.

The [per-read modkit adapter](MODBAM_ADAPTER.md) is pinned to the source-reviewed
0.6.1 extract-full contract. It preserves the MM/ML probability bins despite TSV
decimal rounding and treats alignment ends as exclusive. The separate existing
pileup lane's version policy is unchanged.

## Concrete next experiment

The [public-source research and intake proposal](METHYLATION_ADDITIONAL_DATA_SEARCH.md)
identifies two native iPSC Nanopolish tables from GSE218679, separate from GSE210260.
The two complete compressed files total **208,521 bytes**. This is a targeted FMR1
pilot; exact pore chemistry, Guppy configuration and reference artifact provenance
remain unresolved. Local genomic SHA-256 fields are null until approved acquisition.

`configs/methylation/fmr1_targeted_holdout.v1.json` fixes nine fractions, budgets
20/50/100 and ten replicate seeds: 270 planned levels. Accuracy criteria and all
estimator gates, including the 20-marker minimum, are specified before acquisition.
The source may contain only 19 promoter CpGs, so NO_CALL remains a plausible outcome.
Structural intake must precede outcome inspection and must establish valid source
contracts; no unknown protocol detail may be replaced with a guessed default.

The ignored readiness artifact is
`results/methylation-holdout/20260905/fmr1.readiness.methylation-holdout.json`:

- Decision/status: `NO_CALL` / `NOT_REGISTERED`.
- Expected levels: 270, without biological observations.
- Matrix content SHA-256:
  `fe54b78dcef5f6fd12a70fb8f7735a25591412cc6add2d20a9999d7cba19731e`.
- Code/runtime content SHA-256:
  `0014e1f6c1dd6d89a7f1cfff22b3dceaff084fd52283f56f8d284fdeb444b1fa`.

These are prospective readiness fingerprints, not genomic file checksums or
evidence of recovery. Any subsequent code or policy change needs a fresh lock.

## Verification

`make safety`, `make versions`, `make lint` and `make test` were attempted; `make`
is unavailable in the inspected Windows and WSL runtimes. The underlying Makefile
commands were run directly with the existing Windows `.venv/Scripts` tools.

| Check | Result |
| --- | --- |
| Repository safety and version consistency | Passed |
| Ruff lint and format | Passed; 166 Python files |
| Mypy | Passed; 63 source files |
| Exported schema consistency | Passed |
| Methylation suites | 116 passed |
| Per-read MM/ML/modkit adapter | 47 passed; 1 skipped (pysam absent) |
| CLI surface | 15 passed |
| Repository safety regression tests | 6 passed |
| Full repository suite | 892 tests; 7 failures, 13 errors, 1 skip |

The full suite is **not green**. Its 20 failing/erroring tests are outside the
methylation suites and match the preceding run: Windows path and line-ending
expectations, symlink privileges, process/PID/lock behavior, permissions, service
manifest construction and XLSX cleanup. The untouched local test output is
`results/methylation-holdout/20260905/full-suite.log` (native subprocess encoding).
No tests were suppressed. Neither the real modkit binary nor binary BAM integration
was exercised in this environment; SAM and modkit TSV checks use synthetic fixtures.

## Remaining boundary

The exact two-file acquisition plan and ignored local destination are reviewable in
the research document. AGENTS.md rule 7 requires explicit owner and institutional
approval for external data transfer; genomic acquisition remains pending that input.
Approval does not resolve missing metadata or establish biological eligibility.

The estimand remains the fraction of retained whole read groups from source A.
Current intervals condition on calibration and marker choices. Calibration, donor,
marker and correlation uncertainty, physical dilution truth and conversion to DNA
mass/cell/tumour fractions remain unvalidated. SNV, CNV, ploidy and subclones are
separate evidence layers. No 10–20% detection-limit claim is supported.
