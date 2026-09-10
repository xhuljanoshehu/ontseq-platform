# ONTSeq audit repair implementation plan

**Goal:** repair the reproducible failures on f2776f13 without changing caller thresholds or claiming analytical validation.
**Architecture:** keep the evidence-first portable report, restore omitted typed coverage/SV/ISCN/resource views, and preserve pipeline failure provenance. Changes run on fix/audit-20260909 only.
**Tech stack:** Python/Pydantic, pytest/unittest, GitHub Actions, existing pinned bioinformatics toolchain and .NET desktop.
**Spec:** audit of 2026-09-09; AGENTS.md; docs/ARCHITECTURE.md and docs/CLINICAL_VALIDATION.md.

- [x] Obtain the exact source and run the complete pytest suite on GitHub (Python 3.11/3.12/3.13). Preserve all failures, not only the unittest subset.
- [x] Reproduce report, schema and version failures locally before implementation.
- [ ] Restore report signature with keyword-only target_coverage and selection_coverage; render their values without defaulting missing metrics to zero. Reject sidecars for another sample/build before writing.
- [ ] Restore SV review queue, full event review fields, ISCN fragment dispositions and resource/panel provenance while keeping the current overview/event cards.
- [ ] Repair existing report text/parameter export regressions. Use the existing tests as the contract; retain meaningful biological slash terms.
- [ ] Test unexpected exceptions at plan/execute/settle boundaries, persist FAILED and prevent a release; do not catch KeyboardInterrupt/SystemExit.
- [ ] Add adversarial target-coverage tests for nonfinite depth, nonmonotonic counts, counts/fraction disagreement and contradictory supplied summaries; retain genuinely missing summary values.
- [ ] Investigate the newly discovered methylation-service failures without relaxing result/artifact provenance checks.
- [ ] Synchronize current release identities to 0.8.1, regenerate stale schemas, fix cross-platform ctypes typing and document validation impact.
- [ ] Run all collected tests (not only unittest), typing/lint, real-tool integration, packed runtime and Windows compilation via PR CI. Save outcomes and failures separately.
- [ ] Keep analytical sensitivity/precision/specificity unclaimed until independent positive and negative reference data support them. Audit benchmark/LoD behavior without tuning thresholds on test outcomes.
