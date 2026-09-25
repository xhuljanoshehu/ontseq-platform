# MARLIN Desktop Integration Implementation Plan

**Goal:** Automatically run the build-aware research MARLIN lane with selected methylation and deliver verified status/results through the existing Desktop workflow.
**Architecture:** Separate versioned native research adapter and isolated TensorFlow worker; optional pipeline stage whose failures do not block unrelated evidence. Typed current-run report shared by all presentations. Legacy R/v1 validation contracts retain their existing guarantees.
**Tech Stack:** Python/Pydantic, pinned modkit, separate TensorFlow CPU 2.13.1, WPF/C#, React offline report.
**Spec:** ../specs/2026-09-23-marlin-desktop-integration-design.md (user approved).
**Execution:** Parent integrates pipeline/UI/configuration. Independently scoped workers implement native adapter/tests and report consumers/tests; resource worker prepares local exact-byte resources and synthetic real-tool acceptance. Parent reviews changes and runs complete checks together, with independent final review.

## Global constraints

No patient data in Git or uploads. No automatic clinical release. Native execution is UNVALIDATED_RESEARCH; no invented same-specimen bridge lock. GRCh37 and GRCh38 use separately bound original probe maps. Preserve settings/results and detect active jobs before installation. No clinical threshold changes.

## Review focus

- Stale artifact after failed/reconfigured stage must never reappear in normalized result or export.
- Build/probe mismatch and modified resources must fail before inference and invalidate resume.
- Measured zero and no features must remain distinct; no features produces no model invocation.
- Stock modkit with independent MM groups must not bypass existing qualified-build protection.
- Missing installation must explain NOT_RUN while other analysis modules remain usable.

## Tasks

- [x] Adapter: new marlin_native_contracts.py, marlin_native.py, _marlin_native_worker.py. Tests first in test_marlin_native*.py for weighted 49/100 and 50/100 counts, strands, absent/zero, resources, 42-score invariants and process errors. API run_native_marlin returns NativeMarlinReport; check_native_marlin_readiness gives typed readiness; native_marlin_signature verifies identities for resume.
- [x] Pipeline: add optional StageId.MARLIN after intake, ordered after methylation and before assembly. Auto-request AnalysisModule.MARLIN when selecting methylation. Add installation path to RunConfiguration and profile/service/CLI configuration. New pipeline/marlin.py owns current-artifact loading, stage planning/execution and normalized outcome/provenance. Tests execute real pipeline bookkeeping with controlled external adapters; exercise missing config, success, stale/checksum mutation, failure and resume.
- [x] Presentation: render_html/render_workbook/interactive_payload accept marlin_report keyword. Shared helper validates run/sample/build/status and evidence provenance. Replace static placeholder. Verify escaping, below-threshold UNKNOWN, zero/no-data, all statuses and legacy missing records. Archived export validates recorded artifact checksum or reconstructs blocked status only from matching run record.
- [x] Operator flow: authenticated readiness endpoint reports selected-build prerequisites; Desktop/browser show MARLIN readiness alongside methylation and automatic inclusion. Per-service trusted installation config, never arbitrary browser-supplied executable/model paths. Use existing token/host/origin guards. Test endpoint protection, build errors, UI request serialization and missing-resource explanation.
- [x] Local acceptance: rehash existing original model/runtime/maps, construct exact installation manifest, run real synthetic BAM through qualified modkit and model; independently check weighted probe values/tensor identity and no-evidence refusal. No patient analysis. Preserve exact receipts and existing failed historical receipts.
- [ ] Final: required make safety versions lint test; report frontend and Windows build/tests; independent review and fixes; PR/exact-head CI/merge; build exact merged source and install with rollback and active-job check; actual installed Desktop/WSL acceptance. Document technical limits and validation impact in changelog/validation docs.

## Check commands

Use the existing developer Python environment and PYTHONPATH=src:tests plus the local development dependency directory. Focused pytest targets first; final `make safety versions lint test`. Build the embedded report using the existing web package scripts. Build and test Desktop with the locally bundled .NET SDK/offline NuGet configuration. Evidence goes outside Git under work/marlin-integration and outputs/Pruefnachweise/marlin-integration.

## Acceptance recorded during implementation

The actual synthetic BAM → qualified modkit → original-model stage passed with four observed
features and a maximum reference-predict score difference of 3.3527612686157227e-08. The actual
stage loader, HTML/XLSX/JSON payload and unchanged-signature resume were exercised; tampered
report bytes were rejected. A fresh run after durable runtime relocation reproduced the
exact tensor and direct-inference scores. The zero-feature case returned NO_CALL without
inference. These are synthetic technical checks, not analytical validation.

Independent review identified and prompted fixes for contradictory module/report states,
an explicit downloadable MARLIN JSON result, and configuration-only readiness wording plus
a configured modkit version check. Full asset/runtime byte checks still occur before
execution and resume; readiness does not claim those checks have already passed.
