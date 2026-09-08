# Public-source review and adapter follow-up — 2026-09-05

**Biological acceptance remains NO_CALL / NOT_REGISTERED.** No independent
biological recovery experiment ran. The existing uncommitted work, frozen
matrices and historical outcomes were retained. The user has no additional
released local inputs. No new genomic objects were acquired or uploaded.

The current [public-source review](METHYLATION_NEXT_SOURCE_REVIEW.md) identifies
test candidates and separates documented eligibility from unresolved metadata.
HG008 normal pancreas/duodenum remains a candidate for a scoped two-source tissue
experiment; no public-source nomination establishes recovery or tumour purity.
Acquisition and runtime preparation remain necessary before executing that study.

## Changes in this follow-up

- Modification adapter `0.1.1` rejects explicit Dorado duplex offspring using
  `dx`, including without a recognizable program header. It retains simplex
  parents as read groups and makes no independent-molecule claim.
- modkit extract-full rows must describe alignment spans contained within the
  locked FASTA contig. A valid individual CpG cannot override an invalid span.
- Biological sample, donor, source and cell-type identity fields reject invisible
  control/format characters; a synthetic NUL-suffixed duplicate-donor registration
  regression is now rejected. Existing space normalization and identity case persist.
- Documentation distinguishes numerical evidence replay from re-executing the
  seeded selections on the checksummed local inputs. Selection hashes alone do
  not authenticate original read choices.

All numerical estimator equations and versioned accuracy/marker thresholds are
unchanged. Validation impact and remaining limitations are recorded in
`CLINICAL_VALIDATION.md` and `MODBAM_ADAPTER.md`; the changelog records both fixes.
New registrations must bind the changed adapter and source/runtime versions.

## Current local evidence

Both acquired GSE218679 gzip sizes and SHA-256 values were rechecked against
`configs/methylation/acquired_public_sources.v1.json` and still match. This check
did not select markers or inspect recovery outcomes. The historical intact-read
intake recorded 68/77 test read groups; those counts were not recomputed here.
They remain insufficient for the frozen matrix's largest budget of 100. Exact
pore chemistry, Guppy configuration and the upstream reference artifact digest
remain unresolved; no substitute metadata or complete cohort was invented.

The new ignored readiness artifact is:

`results/methylation-validation/20260905-source-review-v1/source-review.readiness.methylation-holdout.json`

It retains 270 prospective levels with no experimental observations:

- Matrix content SHA-256:
  `fe54b78dcef5f6fd12a70fb8f7735a25591412cc6add2d20a9999d7cba19731e`.
- Source/runtime content SHA-256:
  `354ad890740666409a5d5a236312fbb27539aabd5a0da8e5b4093bb762d141fd`.
- Decision/status: `NO_CALL` / `NOT_REGISTERED`.

No new bias, MAE, RMSE, interval coverage or NO_CALL-rate measurement exists.
Prior GSE210260 failure evidence is preserved. Calibration, donor, marker and
correlation uncertainty remain future studies; no DNA-mass, cell, tumour,
subclone or 10–20% detection-limit claim follows from these software checks.

## Verification

`make safety`, `make versions`, `make lint` and `make test` were attempted;
`make` is absent from the current Windows PATH. All underlying Makefile commands
were then run using the existing `.venv/Scripts` runtime, without installations.

| Check | Result |
| --- | --- |
| Repository safety and version consistency | Passed |
| Ruff lint and format | Passed; 169 Python files |
| Mypy | Passed; 63 source files |
| Exported schema consistency | Passed; no schema structure change |
| Methylation suites | 130 passed |
| modBAM suite | 52 passed; 1 skipped because pysam is absent |
| Full repository suite | 911 tests; 7 failures, 13 errors, 1 skip |

The 20 full-suite failure/error identifiers exactly match the previous
`post-intake-safety-full-suite.log`; none was added or resolved. They remain
outside the methylation changes: Windows path/line-ending and permission
assumptions, symlinks, process/lock handling, service manifests and XLSX cleanup.
The complete suite is **not green**; no failing test was suppressed.

New logs and command exit codes are stored in the ignored directory
`results/methylation-validation/20260905-source-review-v1/`.
`pysam`, `samtools`, `modkit` and `dorado` are absent from the inspected current
Windows environment. SAM and TSV regressions use synthetic fixtures; binary BAM
and actual modkit executable conformance have not been demonstrated locally.
