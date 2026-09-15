# ONTSeq Clean-Room / One-Click Acceptance Design

## Status

**Research Use Only. Engineering acceptance only. No analytical or clinical release claim.**

This design defines a one-click workstation acceptance workflow whose purpose is to answer a practical question:

> Can a user take the complete ONTSeq engineering bundle onto a clean Windows workstation, provision the supported local runtime and one selected genome-build resource family, execute the installed synthetic system test, and receive an auditable readiness verdict without opening a terminal or manually troubleshooting Python, R, Conda, package paths, or reference files?

The workflow is an orchestration layer over existing ONTSeq installation, resource-management, self-test, and provenance contracts. It must not create a second bioinformatics pipeline or duplicate the existing `system-smoke` implementation.

## Scope

The first acceptance scope covers the currently packaged aligned-BAM ONTSeq platform:

- Windows Desktop application;
- WSL2 reachability;
- bundled relocatable Linux runtime installation;
- exact Core/Desktop version identity;
- runtime tools and packaged policies already enforced by Desktop preflight;
- one selected managed resource family: `GRCh37` or `GRCh38`;
- resource installation, repair, status and full validation;
- installed-runtime synthetic self-test;
- QC, BAM intake, CNV, SV, reporting, release checksums and content-addressed resume as exercised by the existing system smoke;
- persistent machine-readable acceptance evidence;
- operator-facing readiness summary and actionable remediation.

The first implementation does **not** silently claim full raw-signal or complete MARLIN readiness. The following remain separate gates:

- POD5/raw-signal -> Dorado -> alignment;
- independently installed modkit when methylation is requested;
- MARLIN model/resources and frozen numerical runtime-compatibility fixture;
- analytical validation on intended-use real/reference specimens;
- clinical validation.

Those capabilities may be displayed in the acceptance report as optional or separately gated modules, but absence of an optional module must never be mislabeled as a successful validation of that module.

## Design principles

1. **Reuse, do not fork.** Existing `WslServiceLauncher` methods remain the authority for WSL checks, runtime installation, resource installation/repair, resource status and `RunSelfTestAsync`.
2. **One operator action.** The setup window exposes one primary action: `Vollständig einrichten & testen`.
3. **Fail closed.** A required step that cannot be proven successful prevents a final PASS.
4. **No terminal dependency.** Expected setup conditions are handled through the Desktop UI. A stack trace, exit code or raw shell error may be preserved in diagnostic evidence but is not the primary operator message.
5. **No destructive replacement.** Existing installed runtimes and resource roots are not silently deleted. Runtime installation continues to use the existing separate-prefix contract; resource repair uses the existing transactional resource machinery.
6. **No build fallback.** A GRCh37 acceptance run validates only GRCh37 resources; a GRCh38 acceptance run validates only GRCh38 resources. No liftover or cross-build substitution is introduced.
7. **Evidence is tiered.** `technical PASS != analytical validity != clinical validity` remains explicit in UI and report.

## Architecture

A new Desktop orchestration component owns the acceptance state machine while delegating all real work to existing services.

```text
SetupWindow
   |
   | one click
   v
SystemAcceptanceRunner
   |
   +--> CheckWslAsync
   |
   +--> CheckBackendAsync
   |      |
   |      +--> InstallBundledRuntimeAsync when absent/mismatched
   |
   +--> CheckResourceFamiliesAsync
   |      |
   |      +--> InstallProfileResourcesAsync when absent
   |      +--> RepairProfileResourcesAsync when incomplete
   |
   +--> checksum-validate each managed bundle in the selected family
   |
   +--> RunSelfTestAsync
   |
   +--> verify self-test evidence exists and is PASS
   |
   +--> write acceptance.json + acceptance.txt
   v
SetupWindow readiness summary
```

`SystemAcceptanceRunner` must not know shell syntax and must not implement bioinformatics logic. It receives a small abstraction over the required Desktop service operations so the orchestration can be unit-tested without WSL, network downloads, or a real workstation.

## Acceptance request

The orchestration input is intentionally small:

```text
selected genome build: GRCh37 | GRCh38
resource root: normalized DesktopSettings.ResourceRootWsl
runtime bundle: AppContext.BaseDirectory/runtime/ontseq-linux-runtime.tar.gz
```

The first release accepts exactly one selected resource family per acceptance run. Installing both GRCh37 and GRCh38 in one click is not required for readiness because a real analysis uses one explicit profile/build contract at a time.

## Acceptance steps

### 1. Preflight WSL

Call the existing `CheckWslAsync`.

- PASS: configured WSL distribution starts and responds.
- FAIL: acceptance stops before any runtime/resource mutation.
- Operator message: WSL2 / the configured Ubuntu distribution must be installed and able to start.

The acceptance runner does not enable Windows optional features or reboot Windows automatically because those operations require elevated OS-level changes and can interrupt the workstation.

### 2. Verify or install bundled ONTSeq runtime

Call `CheckBackendAsync`.

If it already proves the exact release and required Desktop capabilities, record PASS and do not reinstall.

If absent, outdated or capability-incomplete:

1. verify the bundled runtime archive exists;
2. call `InstallBundledRuntimeAsync`;
3. call `CheckBackendAsync` again;
4. require PASS after installation.

A post-install mismatch is a hard failure. The acceptance report records whether installation was required.

### 3. Normalize and persist resource root

Use the same `DesktopSettings.NormalizeResourceRootWsl` contract already used by setup.

The normalized root is persisted before resource operations. Changing a root does not move old data and must not make an old installation appear valid in the new location.

### 4. Provision selected resource family

Call `CheckResourceFamiliesAsync` and inspect only the selected build.

- `Ready`: continue without download.
- `NotInstalled`: call `InstallProfileResourcesAsync` for the selected build.
- `Incomplete`: call `RepairProfileResourcesAsync` for the selected build.
- unavailable/unknown state: fail closed rather than guessing.

After install/repair, refresh status and require `CanAnalyze == true` for the selected build.

The UI must display that several GB and network access may be required before the operator starts the one-click process.

### 5. Full selected-family resource validation

A fast status result is not the final acceptance proof. The backend CLI already supports targeted full validation through:

```text
ontseq references validate <bundle_id> --resource-root <root>
```

The Desktop bridge is extended so `validate` may receive an explicit bundle ID. `SystemAcceptanceRunner` then obtains the selected family's exact managed IDs from `WslServiceLauncher.ManagedResourceBundleIds(genomeBuild)` and checksum-validates **each of those bundle IDs individually**.

All selected-family bundle validations must succeed. A bundle from the other genome-build family is outside this run's scope and must neither be installed implicitly nor cause failure merely because it is absent or independently incomplete.

After the targeted checksum validations, refresh `CheckResourceFamiliesAsync` and again require the selected family to report `CanAnalyze == true`. This proves both bundle integrity and profile-resolution readiness without turning the other build into a hidden prerequisite.

The report records:

- selected build;
- selected resource root;
- exact managed bundle IDs validated for that build;
- per-bundle validation outcome and detail;
- final selected-family status after validation.

No file hash is invented or duplicated in the Desktop layer; the backend resource manifests remain authoritative.

### 6. Installed-runtime system self-test

Call the existing `RunSelfTestAsync`.

The current self-test is already the canonical installed-runtime engineering test and exercises deterministic synthetic data through real runtime tools. The acceptance workflow reuses it rather than constructing another synthetic pipeline.

A PASS must leave the existing self-test evidence folder and the system-smoke report produced by the backend. Acceptance must retain the path to that evidence.

### 7. Evidence verification

The acceptance runner does not merely trust that `RunSelfTestAsync` returned a string. It verifies that the returned evidence location exists and that the expected machine-readable self-test report is present and declares the successful verdict expected by the current Desktop/Core contract.

If the backend self-test implementation changes its report schema in a later release, the corresponding Desktop acceptance parser must change with the same release rather than silently accepting an unknown schema.

### 8. Persist acceptance result

Write a new timestamped evidence directory under:

```text
%LOCALAPPDATA%\ONTSeq\acceptance\<UTC timestamp>\
```

At minimum it contains:

```text
acceptance.json
acceptance.txt
```

`acceptance.json` is the machine-readable authority. `acceptance.txt` is a concise operator summary generated from the same in-memory model.

The acceptance evidence references the self-test directory; it does not duplicate large BAM/report payloads.

## Result model

### Overall verdict

```text
PASS
FAIL
CANCELLED
```

There is no generic `WARNING` verdict that can hide a failed required step.

`PASS` means all required steps for the selected aligned-BAM Core scope succeeded on this workstation at this time.

It does **not** mean analytical or clinical validation.

### Step status

Each required step uses:

```text
PASS
FAIL
NOT_RUN
```

Cancellation is represented at the run level; downstream steps after cancellation are `NOT_RUN`.

Each step records:

- stable step ID;
- status;
- short operator summary;
- technical detail;
- whether an automatic action was taken;
- start and finish timestamps.

Stable first-version step IDs:

```text
wsl
runtime
resource_root
resource_family
resource_validation
system_self_test
self_test_evidence
```

## Module-readiness section

The report contains a separate module-readiness table so that a Core PASS cannot be misread as proof that optional/future modules are ready.

Initial module rows:

```text
aligned_bam_core
methylation_modkit
marlin
raw_signal_dorado
```

Initial states:

```text
READY
NOT_CONFIGURED
GATED
NOT_INCLUDED
```

For the current work package:

- `aligned_bam_core = READY` only after overall acceptance PASS;
- `methylation_modkit` reflects whether the separately configured supported modkit is available, but does not block Core PASS;
- `marlin = GATED` until model/resource lock plus frozen numerical runtime compatibility have passed;
- `raw_signal_dorado = NOT_INCLUDED` until the Desktop product exposes and validates that raw-signal workflow.

No optional module state may be converted into biological interpretation.

## UI

The existing `SetupWindow` remains the only setup surface.

Add one visually primary action below the current individual controls:

```text
Vollständig einrichten & testen
```

The selected `Profilressourcen` build controls which resource family the acceptance workflow provisions.

Before starting, the UI displays one confirmation explaining:

- WSL must already work;
- runtime installation may occur automatically;
- the selected resource family may download several GB;
- existing valid installations are reused;
- the synthetic self-test may take time;
- the result is engineering/RUO readiness, not clinical validation.

During execution the button is disabled and the existing busy/cancellation machinery is reused. A compact progress line shows the current stable step rather than raw command output.

On PASS:

```text
✓ ONTSeq Engineering Acceptance: PASS
Selected build: GRCh37 | GRCh38
Ready for authorized research aligned-BAM analysis: YES
Analytical/clinical validation: NO
```

On failure the UI shows exactly one primary remediation message associated with the failed step and the path to the saved acceptance evidence.

Examples:

```text
WSL2/Ubuntu is not available. Install or repair WSL, then rerun this acceptance test.
```

```text
GRCh38 resources could not be validated after automatic repair. No research analysis was enabled. The acceptance evidence is stored below %LOCALAPPDATA%\ONTSeq\acceptance.
```

The UI must not ask the operator to run pip, conda, R package installation, or ad-hoc shell commands for expected setup states.

## Automatic remediation boundary

The one-click workflow may automatically perform only operations already supported by current Desktop contracts:

- install bundled ONTSeq runtime;
- install selected managed resource family;
- repair selected managed resource family;
- save the normalized resource root;
- execute the existing system self-test.

It must **not** automatically:

- enable WSL Windows features;
- reboot Windows;
- delete existing runtimes or resources;
- modify mapped network-drive configuration;
- install an arbitrary external modkit binary;
- download/distribute MARLIN model artifacts before distribution/provenance policy is separately approved;
- modify patient/research inputs.

## Determinism and repeatability

Running acceptance twice without changing software or resources should not reinstall valid components.

The second run should:

- reuse the valid runtime;
- reuse validated resource bundles;
- rerun the deterministic system self-test;
- create a new timestamped acceptance record;
- retain the existing self-test's content-addressed resume assertions.

The acceptance result does not need to have byte-identical timestamps, but the technical identities and PASS/FAIL conclusions must remain consistent when the environment is unchanged.

## Error handling

Expected operational errors are mapped to stable acceptance steps and user-facing remediation text.

Unexpected exceptions are still fail-closed. Their exception text and launcher diagnostics may be stored in `technical_detail`, but the UI displays a concise message and evidence path.

No exception is swallowed to produce PASS.

If the process is cancelled:

- overall verdict is `CANCELLED`;
- completed steps retain their status;
- remaining steps are `NOT_RUN`;
- no `Ready for authorized research aligned-BAM analysis: YES` statement is displayed from that run.

## Test strategy

### Unit tests

The orchestration is tested against a fake acceptance service interface. Required cases:

1. clean system with working WSL, missing runtime, missing selected resources -> install runtime -> install resources -> validate selected bundles -> self-test -> PASS;
2. already-ready system -> no reinstall/re-download -> selected-bundle validation -> self-test -> PASS;
3. incomplete selected resources -> repair -> selected-bundle validation -> PASS;
4. WSL unavailable -> stop before mutation -> FAIL;
5. runtime install followed by version/capability mismatch -> FAIL;
6. resource install/repair followed by non-ready status -> FAIL;
7. any selected managed bundle failing checksum validation -> FAIL;
8. a broken or absent resource family from the other genome build does not block the selected family acceptance;
9. self-test failure -> FAIL;
10. self-test reports success but evidence/report is missing or invalid -> FAIL;
11. cancellation -> CANCELLED with downstream `NOT_RUN`;
12. optional modkit/MARLIN/raw-signal readiness never upgrades Core evidence beyond its declared scope.

### Desktop integration tests

Extend `desktop/ONTSeq.Desktop.Tests` to verify:

- result-model serialization;
- stable step IDs and verdict semantics;
- UI-facing summary generation;
- targeted `references validate <bundle_id>` argument construction;
- existing SetupWindow state rules remain compatible;
- acceptance evidence paths remain under `%LOCALAPPDATA%\ONTSeq\acceptance`.

### CI

Existing Desktop CI remains authoritative for packaging and current runtime smoke. Add acceptance-contract tests to the Windows Desktop test job.

The existing Linux runtime build continues to prove relocation and stock-Ubuntu execution. The new acceptance test does not replace these CI gates.

### Real clean-room workstation gate

Before declaring the one-click feature complete for operator use, perform one manual clean-room run on a fresh Windows 11 + WSL2/Ubuntu installation:

1. download the complete Desktop engineering artifact;
2. extract it;
3. start `ONTSeq.Desktop.exe`;
4. choose GRCh37 or GRCh38;
5. press only `Vollständig einrichten & testen`;
6. do not open PowerShell, CMD, WSL shell, Python, R, Conda or an editor;
7. require final PASS and inspect the generated acceptance evidence;
8. execute one authorized synthetic/non-patient aligned BAM through the normal analysis button and open HTML/XLSX/JSON from Desktop.

Any need for manual package installation, environment-variable editing, shell troubleshooting, resource-file browsing, or ad-hoc path correction during this gate is a product defect to fix before calling the workflow zero-troubleshooting.

## Acceptance criteria

The feature is engineering-complete when all of the following are true:

1. one Desktop action orchestrates the existing runtime/resource/self-test contracts;
2. a missing runtime is installed automatically from the checksummed bundled runtime;
3. the selected missing/incomplete resource family is installed/repaired automatically;
4. every managed bundle in the selected resource family is checksum-validated individually;
5. an absent or broken non-selected genome-build family does not become an implicit prerequisite;
6. the canonical installed-runtime self-test is mandatory;
7. self-test evidence is independently checked before final PASS;
8. machine-readable and human-readable acceptance evidence is persisted;
9. PASS has no failed or skipped required step;
10. optional methylation/MARLIN/raw-signal states are clearly separated from aligned-BAM Core readiness;
11. unit/Desktop tests cover happy paths, repair paths, failure paths, cross-build isolation and cancellation;
12. Core CI, Desktop CI and existing system-smoke gates remain green;
13. a fresh Windows workstation can complete the documented clean-room procedure without terminal troubleshooting.

## Release boundary

This feature is an engineering deployment/acceptance aid. It does not itself raise ONTSeq from RUO to a clinically validated system.

A workstation PASS supports only the statement that the packaged software, selected managed resources and deterministic synthetic Core workflow executed successfully under the tested environment.

Real/reference specimen performance remains a separate analytical validation program.