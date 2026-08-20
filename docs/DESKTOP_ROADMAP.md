# ONTSeq Desktop implementation gates

## Gate 0 — backend stability

- [x] Local loopback service starts under CI with bounded readiness polling.
- [x] Service uses PID-based cleanup in non-interactive CI.
- [x] Start/status/refusal/review browser path passes with real local tools.
- [x] Existing backend run #144 is green before the desktop branch is stacked on it.

## Gate 1 — Windows shell v0.1

- [x] Native Windows WPF project on .NET 10.
- [x] BAM picker and sample-ID suggestion.
- [x] GRCh37/GRCh38 selector.
- [x] lcWGS/Adaptive Sampling selector.
- [x] Automatic WSL backend launch; no terminal required for normal use.
- [x] Loopback API bootstrap and Start/Status integration.
- [x] Live stage display from atomically persisted `provenance/run.json`.
- [x] Preserve `COMPLETED`, `NO_CALL`, `FAILED`, `NOT_RUN` verbatim.
- [x] Open HTML, XLSX and run folder from the Windows UI.
- [x] Explicit RUO banner.
- [x] Windows CI build plus self-contained unsigned engineering artifact.
- [ ] Real Windows workstation smoke with the locally configured WSL distribution.
- [ ] Real mapped `P:` drive visibility check in the hospital environment.

## Gate 2 — safe cancellation

Do not enable the current disabled button until cancellation is a backend property rather than a UI illusion.

- [ ] `POST /api/runs/<id>/cancel`.
- [ ] job states `cancel_requested` and `cancelled`.
- [ ] command runner capable of terminating the process tree under orchestrator control.
- [ ] interrupted temporary outputs removed/quarantined; previously committed artifacts remain intact.
- [ ] downstream stages recorded as `NOT_RUN` with a cancellation reason.
- [ ] resume after cancellation proven.
- [ ] real-tool failure-injection CI.

## Gate 3 — installer

- [ ] decide signed MSIX/MSI packaging after enterprise-IT constraints are known.
- [ ] bootstrap/check WSL2 and a supported distribution.
- [ ] install/import the pinned offline backend image/environment.
- [ ] install locally approved reference bundles and write `desktop.settings.json`.
- [ ] self-test after installation.
- [ ] deterministic upgrade and rollback.
- [ ] code signing and SBOM attached to release artifacts.

## Gate 4 — review workspace

- [ ] in-app result summary from versioned JSON only.
- [ ] CNV chromosome view.
- [ ] SV/fusion evidence view with observability and caller concordance.
- [ ] explicit review/sign-off state.
- [ ] ISCN remains proposal-only until domain validation is completed.

The desktop layer must remain replaceable. Analytical truth lives behind contracts; UI code must never become a second implementation of the pipeline.
