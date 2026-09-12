# Methylation validation status — 2026-09-05

**Update:** The two-source prospective mode is now separate from the stricter
donor-independent lane: see [`METHYLATION_HOLDOUT.md`](METHYLATION_HOLDOUT.md).
The four-donor gate below describes the first engineering stage, not a prerequisite
for the next two-source recovery experiment. This record and its readiness digest
are historical snapshots; subsequent implementation changes require a new lock.

**Biological acceptance: NO_CALL. No independent biological dilution experiment has run.**
The existing uncommitted work was retained. No genomic input was acquired, staged or committed.
The historical GSE210260 results were not re-run and remain negative technical evidence.

## Completed engineering work

- The public-source [catalogue](METHYLATION_DATA_CANDIDATES.md) records exact acquisition
  candidates, protocol differences, supplier checksums and storage requirements. Local genomic
  SHA-256 values are deliberately null until acquisition; no supplier MD5 or S3 ETag is
  represented as a verified SHA-256.
- The [Dorado adapter](MODBAM_ADAPTER.md) accepts per-read SAM MM/ML, optional BAM/pysam and
  modkit extract-full TSV. It checks coordinates, probability semantics, supported modification
  classes, reference integrity and available SAM/BAM header provenance. Nanopolish remains usable.
- The [prospective study](METHYLATION_VALIDATION.md) fixes 270 fraction/read-budget/seed cells
  with depth and fraction acceptance strata. Registration binds the actual source/runtime code,
  the four-source cohort and matrix. Missing evidence is NO_CALL; failed recovery is FAIL.
- Local CLI commands fingerprint inputs, register, execute and independently recompute reports.
  Reports include JSON, per-level CSV and HTML with stratified bias, MAE, RMSE, interval coverage
  and unconditional NO_CALL rates. Existing output files are protected from overwrite.
- Eleven new exported schemas, synthetic tests, safety checks and biological validation-impact
  documentation accompany the implementation.

The local readiness artifact is
`results/methylation-validation/20260905/readiness.methylation-validation.json`.
It states `NO_CALL` / `NOT_REGISTERED`; its 270 expected levels have no biological observations.
The matrix content digest is
`8476201c445144bb881f2ad7884cb2ef8985b1a9ffec2592f89994500afcb3f8`.
The checked source/runtime content digest is stored in that artifact. This is a prospective
template and readiness record, not an authenticated preregistration or completed validation.

## Verification on the available Windows runtime

The requested `make safety`, `make versions`, `make lint` and `make test` were attempted.
`make` is absent in both the Windows PATH and the inspected WSL runtime. Their underlying
commands were therefore executed directly with the existing `.venv/Scripts` tools:

| Check | Observed result |
| --- | --- |
| Repository safety | Passed |
| Version consistency | Passed |
| Ruff lint and format | Passed, 158 Python files formatted |
| Mypy | Passed, 60 source files |
| Exported schema consistency | Passed |
| Methylation suites, including independent runner | 85 tests passed |
| MM/ML/modkit adapter suite | 42 passed; one binary-BAM integration skipped because pysam is absent |
| CLI surface | 15 tests passed |
| Repository safety regression suite | 6 tests passed |
| Full repository suite | 856 tests; 7 failures, 13 errors, 1 skip |

The 20 full-suite failures/errors lie outside the new methylation tests: Windows path/line-ending
expectations, unavailable symlink privileges, process-lock/PID behavior, permissions, service
manifest construction and an XLSX file remaining open during temporary-directory cleanup.
The full suite is **not green**. No portability changes or test suppressions were introduced to
hide these failures. The complete local output is
`results/methylation-validation/20260905/full-suite.log`.

The synthetic four-source SAM experiment verifies that the new Dorado path reaches the same
locked estimator and report replay. It is not evidence of real-data accuracy. Neither a real
modkit executable comparison nor actual binary-BAM decoding ran in this environment.

## Remaining requirements

The reviewed public candidates provide further technical tissue/cell-line feasibility controls,
but no admitted independent calibration/test donor panel. A user/institution-provided local
cohort can fill that gap. Required roles are calibration A/B and held-out test A/B with coherent
biological source definitions, non-overlapping calibration/test donors, complete compatible
protocol metadata and actual file SHA-256 values. The conservative donor gate is the scope of
this new validation lane; technical repeated flowcells or seeded read splits do not satisfy it.

Any large acquisition needs a concrete local storage/derivative plan and the applicable owner
and institutional authorization under AGENTS.md. The default parser remains bounded in memory;
WGS subsets or a changed budget need a prospectively defined sampling domain and a new study
lock. Do not relax gates after seeing failed test outcomes and report the same data as validation.

Only actual eligible test data can yield the requested independent bias, MAE, RMSE, interval
coverage and NO_CALL-rate evidence. Calibration/donor/marker/correlation uncertainty and any
conversion to DNA mass, cell fraction or tumour fraction remain separate uncompleted studies.
No tumour purity, subclone fraction or 10–20% detection-limit claim is supported.
