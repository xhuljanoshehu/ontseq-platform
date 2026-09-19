# Multi-caller 5A: catalog and declared-input planning

Research Use Only. Human Review Required. Package version unchanged at 0.8.2.

## Delivered boundary

Eight callers have a common immutable catalog and typed, deterministic planning API.
This is not eight live native integrations. Existing QDNAseq+ACE, Sniffles2 and cuteSV
runtime code is unchanged. ONT-Spectre, SAVANA, Severus, Wakhan and ichorCNA remain
explicitly unimplemented native adapters in this increment. All plan lanes are NOT_RUN.
No CLI, desktop execution control, production provider table or report is changed.

The planner checks declared input, reference, assay, mask, coverage definition,
paired-normal identity/depth, study-policy and runtime-pin consistency. It distinguishes
NOT_SELECTED, BLOCKED, ADAPTER_PENDING and PLANNED. It does not inspect genomic files,
verify external binaries, authenticate a registration document or infer analytical validity.

## Observed local TDD evidence

Workspace: an isolated focused source workspace, not a complete repository checkout.
Command: PYTHONPATH=src python -m pytest tests/test_multicaller_plan.py -q

- Catalog RED: missing ontseq_platform.multicaller_catalog module before implementation.
- Catalog GREEN: 4 passed.
- Routing RED: missing planning API before implementation.
- Routing GREEN: 41 passed, 24 subtests passed.
- Scope-binding RED: new reference/mask/measurement-definition fields not yet supported.
- Scope-binding GREEN: 46 passed, 27 subtests passed.
- Schema RED: missing versioned multicaller-request schema snapshot.
- Schema GREEN and final focused rerun: 47 passed, 27 subtests passed.

Full repository CI, Ruff, mypy and runtime-smoke results must be read from the final
head's GitHub Actions runs and recorded in the PR, not inferred from these focused tests.
No native new-caller smoke or real biological cohort was run for this deliverable.

## Remaining integration

The design file specifies native runtime/import/evidence bridges for ONT-Spectre,
SAVANA, Severus/Wakhan and the separate ichorCNA cfDNA route, followed by reviewed
CLI/desktop/report integration. Each requires upstream format/version/resource/license
review, synthetic native fixtures, failure tests and real-tool smoke before activation.
Shared reads and borrowed breakpoints are not independent orthogonal truth.
