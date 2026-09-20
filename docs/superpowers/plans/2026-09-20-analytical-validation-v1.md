# Analytical Validation Program v1 — Implementation Plan

> Branch: `feat/analytical-validation-v1`
> Base: `caca29feab87d6f58dc7d172e42a51a70f036fa0`
> Method: strict RED -> GREEN
> No patient data; no merge/version bump

## Task A1 — factor and design contracts

Files:
- `src/ontseq_platform/cnv_validation_v1.py`
- `tests/test_cnv_validation_v1.py`

RED requirements:
- explicit ordered coverage levels;
- explicit ordered tumour-fraction levels;
- explicit ordered CNV-size levels;
- replicate count and deterministic seed;
- invalid/duplicate levels fail closed.

## Task A2 — deterministic cell expansion

Positive study cells are the Cartesian product of:
- coverage;
- tumour fraction;
- event class;
- compatible event size;
- replicate.

Whole-chromosome classes do not cross numeric size levels.

Cell IDs must be stable and unique.

## Task A3 — synthetic truth/cohort materialisation

For every design cell:
- deterministic specimen/input identity;
- continuous coverage and tumour fraction retained;
- registered truth source;
- registered assessability mask;
- exact synthetic CNV truth event;
- no patient payload.

## Task A4 — v1 matrix builder

Build the existing `CnvValidationMatrix` using supplied caller locks.

Defaults:
- primary: coverage, tumour fraction, event class, event size;
- secondary: caller, genome build, data basis, bin size;
- QDNAseq 100/500/1000 kbp retained;
- no automatic caller ranking.

## Task A5 — registration builder

Use existing `preregister_cnv_validation`.

Tests:
- same design + same locks => same matrix/cohort content hashes;
- changed factor/caller lock => changed registration lock;
- outcomes remain unseen.

## Task A6 — dilution projection

Generate the existing `DilutionPolicy` from v1 tumour-fraction levels.

No new mixing implementation.

## Task A7 — schema/examples/status

Add:
- `schemas/cnv-validation-v1-design.schema.json`
- `examples/validation/analytical_validation_v1.design.json`
- `docs/ANALYTICAL_VALIDATION_V1_STATUS.md`

Example caller locks are explicitly synthetic and cannot be mistaken for real caller qualification evidence.

## Task A8 — full verification

Run:
- safety
- schema freshness
- version consistency
- Ruff
- mypy
- full tests
- local-real-tool-smoke
- existing caller real-tool workflows if touched transitively
- Desktop bundle/CI if packaging-relevant files change

No merge. No package version bump.
