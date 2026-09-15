# Clean-Room Acceptance Implementation Plan Amendment

**Amends:** `docs/superpowers/plans/2026-09-15-clean-room-acceptance.md`

This amendment is normative for execution. It replaces only the specific Task 1, Task 6 and Task 7 details below. All other tasks, interfaces, constraints and commit boundaries in the base implementation plan remain unchanged.

## 1. Task 1 — validate the published expected profile set, not only indexed profiles

### Reason

`ResourceRegistry.profiles` contains only profiles whose manifest parsed and matched the active build. If one expected selected-build profile is missing or corrupted, iterating only `registry.profiles` could silently omit it from build-scoped validation. A workstation acceptance PASS must instead prove that every published profile expected for the selected build resolves with `verify_files=True`.

### Required imports

In `src/ontseq_platform/resource_commands.py`, extend the bootstrap imports to include:

```python
from .resource_bootstrap import (
    GRCH37_PROFILE_IDS,
    GRCH37_REFERENCE_BUNDLE_ID,
    PROFILE_IDS,
    REFERENCE_BUNDLE_ID,
    GRCh37ResourceBootstrapper,
    GRCh38ResourceBootstrapper,
)
```

### Replace `_profile_status` signature

Replace the base-plan helper contract with:

```python
def _profile_status(
    registry: ResourceRegistry,
    *,
    verify_checksums: bool,
    profile_ids: tuple[str, ...] | None = None,
) -> list[dict[str, object]]:
    """Resolve the requested profiles and preserve explicit failures.

    ``profile_ids=None`` retains legacy status/unscoped-validation behaviour by using the
    profiles currently indexed by the registry. Build-scoped acceptance supplies the complete
    published profile set so a missing/corrupt profile cannot disappear from validation.
    """

    selected = tuple(sorted(registry.profiles)) if profile_ids is None else profile_ids
    status: list[dict[str, object]] = []
    for profile_id in selected:
        try:
            context = registry.resolve_profile(profile_id, verify_files=verify_checksums)
        except (KeyError, OSError, ValueError) as error:
            status.append(
                {
                    "profile_id": profile_id,
                    "valid": False,
                    "error": str(error),
                }
            )
            continue
        status.append(
            {
                "profile_id": profile_id,
                "valid": True,
                "reference_bundle": context.reference_bundle_id,
                "panel_bundle": context.panel_bundle_id,
                "knowledge_bundle": context.knowledge_bundle_id,
            }
        )
    return status
```

### Replace `_validate_build_scope`

Use the published profile IDs and the published reference bundle ID for the selected build:

```python
def _validate_build_scope(
    args: argparse.Namespace,
    installer: ReferenceBundleInstaller,
    build: GenomeBuild,
) -> tuple[tuple[BundleValidationReport, ...], list[dict[str, object]]]:
    registry = ResourceRegistry(args.resource_root, active_build=build)
    if build is GenomeBuild.GRCH37:
        profile_ids = GRCH37_PROFILE_IDS
        reference_bundle_id = GRCH37_REFERENCE_BUNDLE_ID
    elif build is GenomeBuild.GRCH38:
        profile_ids = PROFILE_IDS
        reference_bundle_id = REFERENCE_BUNDLE_ID
    else:
        raise AssertionError(f"unsupported build: {build}")

    reports = (installer.validate(reference_bundle_id),)
    profile_status = _profile_status(
        registry,
        verify_checksums=True,
        profile_ids=profile_ids,
    )
    return reports, profile_status
```

This intentionally does not infer the expected family from filenames or from whatever survived registry indexing.

### Replace Task 1 tests

For GRCh37, configure the mock so every ID in `GRCH37_PROFILE_IDS` resolves and assert all calls:

```python
registry = MagicMock()
registry.profiles = {}
registry.diagnostics = ()
context = registry.resolve_profile.return_value
context.reference_bundle_id = GRCH37_REFERENCE_BUNDLE_ID
context.knowledge_bundle_id = GRCH37_KNOWLEDGE_BUNDLE_ID
context.panel_bundle_id = None

# execute build-scoped validation

assert registry.resolve_profile.call_args_list == [
    call(profile_id, verify_files=True) for profile_id in GRCH37_PROFILE_IDS
]
```

Add a fail-closed test in which one expected ID raises `KeyError` even though the other three resolve. The command must terminate with `SystemExit(2)` and JSON/profile status must contain an explicit invalid row for the missing profile.

Add the equivalent GRCh38 expected-set test using `PROFILE_IDS`.

Do not add a new pytest dependency/import solely for the mutual-exclusion test. Follow the existing unittest style:

```python
with self.assertRaisesRegex(ValueError, "mutually exclusive"):
    handle_references_command(args)
```

The existing unscoped status/validation test must continue to call `_profile_status(..., profile_ids=None)` implicitly and preserve its legacy all-build behaviour.

## 2. Task 6 — independently require every system-smoke check to PASS

The base plan already requires schema `0.1.0`, overall `verdict == "PASS"`, non-empty `checks`, and an `output_paths` object. Strengthen the independent evidence parser so an internally inconsistent report cannot pass merely because its top-level verdict says PASS.

After obtaining `checks`:

```csharp
var checks = root.GetProperty("checks");
if (checks.ValueKind != JsonValueKind.Array || checks.GetArrayLength() < 1)
    return (false, "System-Smoke enthält keine Prüfnachweise.");
foreach (var check in checks.EnumerateArray())
{
    if (check.ValueKind != JsonValueKind.Object ||
        !check.TryGetProperty("status", out var status) ||
        !string.Equals(status.GetString(), "PASS", StringComparison.Ordinal))
        return (false, "Mindestens ein System-Smoke-Prüfnachweis ist nicht PASS.");
}
```

Add a RED fixture with top-level `verdict="PASS"` but one check `status="FAIL"`; `VerifySelfTestEvidence` must refuse it.

This is evidence consistency checking only. Do not duplicate the backend's CNV/SV/report scientific assertions in Desktop.

## 3. Task 7 — preserve the final acceptance summary after the run

### Reason

The base-plan pseudocode calls `await RefreshAsync()` after writing the final acceptance summary. `RefreshAsync()` overwrites `DetailText` with the generic backend-readiness message, hiding the PASS/FAIL remediation and evidence directory that the operator needs.

### Required handler pattern

Do not call the full `RefreshAsync()` after `RunBusyAsync`.

Instead, refresh bundle status *inside* the successful action before writing the final acceptance message, then write the acceptance result last:

```csharp
private async void Acceptance_Click(object sender, RoutedEventArgs e)
{
    var build = SelectedResourceBuild;
    var confirm = MessageBox.Show(
        this,
        $"ONTSeq richtet die gebündelte Runtime und die ausgewählte {build}-Ressourcenfamilie ein und führt danach den vollständigen synthetischen Engineering-Selbsttest aus.\n\n" +
        "WSL2/Ubuntu muss bereits funktionieren. Fehlende Ressourcen können mehrere GB herunterladen. Vorhandene gültige Installationen werden wiederverwendet. Der Test kann einige Zeit dauern.\n\n" +
        "Das Ergebnis bestätigt nur Engineering/RUO-Bereitschaft, keine analytische oder klinische Validierung.\n\nFortfahren?",
        "ONTSeq Engineering-Acceptance",
        MessageBoxButton.YesNo,
        MessageBoxImage.Information);
    if (confirm != MessageBoxResult.Yes) return;

    await RunBusyAsync(async token =>
    {
        _settings.ResourceRootWsl = DesktopSettings.NormalizeResourceRootWsl(ResourceRootTextBox.Text);
        var service = new WslSystemAcceptanceService(_launcher, _settings);
        var evidenceRoot = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
            "ONTSeq", "acceptance");
        var runner = new SystemAcceptanceRunner(
            service,
            new SystemAcceptanceEvidenceWriter(evidenceRoot));
        var archive = Path.Combine(
            AppContext.BaseDirectory,
            "runtime",
            "ontseq-linux-runtime.tar.gz");
        var progress = new Progress<SystemAcceptanceProgress>(item =>
        {
            AcceptanceStatusText.Text = "● " + item.Message;
            DetailText.Text = item.Message;
        });

        var execution = await runner.RunAsync(
            new SystemAcceptanceRequest(build, _settings.ResourceRootWsl, archive),
            progress,
            token);

        if (execution.Report.Verdict == AcceptanceVerdict.PASS)
            await RefreshBundleStatusesAsync(token);

        AcceptanceStatusText.Text = execution.Report.Verdict switch
        {
            AcceptanceVerdict.PASS => "✓ PASS",
            AcceptanceVerdict.CANCELLED => "— CANCELLED",
            _ => "✕ FAIL"
        };
        DetailText.Text = execution.Report.ToOperatorSummary() +
            "Prüfnachweise: " + execution.EvidenceDirectory;
    }, "Engineering-Acceptance konnte nicht abgeschlossen werden.");
}
```

Do not append `await RefreshAsync();` afterward.

### Add a structural regression assertion

The Desktop contract test should assert the acceptance handler source does not contain a post-run `await RefreshAsync();` pattern. Prefer testing an extracted presentation helper if implementation naturally creates one; otherwise keep this as a focused source-contract test alongside the XAML structural checks.

## 4. Task 5 clarification — configured modkit remains GATED, never READY

A configured `ModkitExecutableWsl` path only proves that a path was configured. This Core acceptance workflow does not execute the methylation lane, so its module table must use:

```text
null/blank modkit path -> NOT_CONFIGURED
configured modkit path -> GATED
```

The existing methylation/run preflight remains the authority for whether that executable is present, pinned to the supported version and usable for a methylation-enabled analysis. Do not add a new modkit probe to the clean-room Core acceptance scope.

## 5. Updated final verification

In addition to the base plan's full test matrix, explicitly require these regression properties before calling the implementation complete:

```text
- every one of the four expected profiles for the selected build is checksum-resolved;
- a missing expected selected-build profile makes acceptance fail;
- a broken non-selected build remains irrelevant to the selected build acceptance;
- top-level system-smoke PASS with an internal FAIL check is rejected;
- the final UI summary still shows verdict + evidence path after the run;
- configured modkit is GATED, not READY, in Core acceptance evidence.
```
