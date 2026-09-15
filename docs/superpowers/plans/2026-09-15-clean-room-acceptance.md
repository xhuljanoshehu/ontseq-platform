# Clean-Room / One-Click Acceptance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a one-click Windows/WSL engineering acceptance workflow that can install or reuse the bundled ONTSeq runtime, provision exactly one selected GRCh37/GRCh38 resource family, fully validate that family, execute the existing installed-runtime system smoke, verify its evidence, and persist a machine-readable PASS/FAIL/CANCELLED readiness record without requiring terminal troubleshooting.

**Architecture:** `SetupWindow` remains the only operator setup surface. A new `SystemAcceptanceRunner` owns orchestration and depends on a small `ISystemAcceptanceService` abstraction; production delegates to the existing `WslServiceLauncher`, while tests use a fake. Bioinformatics execution remains entirely in the existing backend/runtime and `system-smoke`; the new code only orchestrates, validates evidence, persists acceptance records, and presents remediation.

**Tech Stack:** .NET 10, WPF, C# nullable reference types, existing executable Desktop contract tests, Python >=3.11, pytest, Pydantic/resource registry backend, WSL2, existing GitHub Desktop CI.

**Spec:** `docs/superpowers/specs/2026-09-15-clean-room-acceptance-design.md` plus `docs/superpowers/specs/2026-09-15-clean-room-acceptance-design-amendment.md`.

## Global Constraints

- Work only on `feat/marlin-classification-v1`; do not modify `main` directly.
- Research Use Only. Engineering acceptance is not analytical or clinical validation.
- Reuse the canonical `WslServiceLauncher` installation/resource/self-test paths; do not create a second bioinformatics pipeline.
- Required Core acceptance scope is aligned-BAM ONTSeq only.
- WSL2 must already be installed and able to start; the workflow must not enable Windows optional features or reboot Windows.
- A run validates exactly one selected genome-build family (`GRCh37` or `GRCh38`). The non-selected family is not a hidden prerequisite.
- No runtime liftover, build fallback, resource filename guessing, or cross-build substitution.
- Expected setup failures must be presented as concise remediation, not raw shell instructions.
- Existing valid runtime/resource installations are reused; do not silently delete them.
- `technical PASS != analytical validity != clinical validity` remains explicit.
- Optional methylation/modkit, MARLIN and raw-signal/Dorado readiness are reported separately and never upgrade aligned-BAM Core evidence.
- The existing system-smoke schema is `0.1.0`; final acceptance requires `system-smoke.report.json` with `verdict == "PASS"`.
- Keep Desktop/Core identity at `0.8.2` during this work; no release/version bump is part of this plan.

---

## File Structure

### New Desktop files

- `desktop/ONTSeq.Desktop/SystemAcceptanceModels.cs` — enums, stable step IDs, request/report/result models, summary formatting.
- `desktop/ONTSeq.Desktop/SystemAcceptanceRunner.cs` — orchestration state machine and evidence verification/persistence.
- `desktop/ONTSeq.Desktop/SystemAcceptanceService.cs` — `ISystemAcceptanceService` and production adapter around `WslServiceLauncher` + `DesktopSettings`.
- `desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs` — fake service and orchestration/evidence tests.

### Existing files modified

- `src/ontseq_platform/resource_commands.py` — build-scoped `references validate --genome-build`.
- `tests/test_grch37_resource_bootstrap.py` — GRCh37 scoped validation regression tests.
- `tests/test_resource_bootstrap.py` — GRCh38 scoped validation regression tests.
- `desktop/ONTSeq.Desktop/WslServiceLauncher.cs` — build-scoped validation bridge.
- `desktop/ONTSeq.Desktop/SetupWindow.xaml` — one-click acceptance control/status row.
- `desktop/ONTSeq.Desktop/SetupWindow.xaml.cs` — invoke runner, progress, result presentation, busy-state wiring.
- `desktop/ONTSeq.Desktop.Tests/Program.cs` — invoke the new contract-test suite and bridge assertions.
- `docs/DESKTOP_FIRST_RUN.md` — zero-troubleshooting operator path and evidence location.
- `desktop/README-FIRST-RUN.md` — bundle-local clean-room procedure.
- `.github/workflows/desktop-ci.yml` — assert the published first-run bundle contains the new one-click instructions; no new analytical lane.

---

## Task 1: Build-scoped full resource validation in the backend

**Files:**
- Modify: `src/ontseq_platform/resource_commands.py`
- Modify: `tests/test_grch37_resource_bootstrap.py`
- Modify: `tests/test_resource_bootstrap.py`

**Interfaces:**
- Existing: `ontseq references validate [bundle_id] --resource-root <root>`.
- Produces: `ontseq references validate --genome-build GRCh37|GRCh38 --resource-root <root> --json`.
- JSON payload for build-scoped validation:

```json
{
  "genome_build": "GRCh38",
  "references": [],
  "profile_status": []
}
```

- `bundle_id` and `--genome-build` are mutually exclusive.

- [ ] **Step 1: Write the GRCh37 RED test**

Add a test to `tests/test_grch37_resource_bootstrap.py` that constructs two mocked registries, selects `GRCh37`, and proves only the GRCh37 registry is instantiated/resolved with checksum verification:

```python
def test_build_scoped_validate_checks_only_grch37() -> None:
    registry = MagicMock()
    profile = MagicMock()
    profile.reference_bundle = GRCH37_REFERENCE_BUNDLE_ID
    registry.profiles = {"AML_LCWGS_GRCh37": profile}
    registry.diagnostics = ()
    context = registry.resolve_profile.return_value
    context.reference_bundle_id = GRCH37_REFERENCE_BUNDLE_ID
    context.knowledge_bundle_id = GRCH37_KNOWLEDGE_BUNDLE_ID
    context.panel_bundle_id = None
    installer = MagicMock()
    report = MagicMock()
    report.bundle_id = GRCH37_REFERENCE_BUNDLE_ID
    report.valid = True
    report.bundle_path = Path("resources/references") / GRCH37_REFERENCE_BUNDLE_ID
    report.manifest_valid = True
    report.errors = ()
    report.resources = ()
    installer.validate.return_value = report
    args = argparse.Namespace(
        command="references",
        references_command="validate",
        resource_root=Path("resources"),
        config_root=CONFIGS,
        as_json=True,
        bundle_id=None,
        genome_build="GRCh37",
    )
    with (
        patch("ontseq_platform.resource_commands.ReferenceBundleInstaller", return_value=installer),
        patch("ontseq_platform.resource_commands.ResourceRegistry", return_value=registry) as registry_type,
        redirect_stdout(io.StringIO()) as output,
    ):
        assert handle_references_command(args)
    registry_type.assert_called_once_with(Path("resources"), active_build=GenomeBuild.GRCH37)
    registry.resolve_profile.assert_called_once_with("AML_LCWGS_GRCh37", verify_files=True)
    installer.validate.assert_called_once_with(GRCH37_REFERENCE_BUNDLE_ID)
    payload = json.loads(output.getvalue())
    assert payload["genome_build"] == "GRCh37"
```

- [ ] **Step 2: Write matching GRCh38 and incompatibility RED tests**

In `tests/test_resource_bootstrap.py`, add equivalent GRCh38 coverage plus:

```python
def test_validate_rejects_bundle_id_and_genome_build_together() -> None:
    args = argparse.Namespace(
        command="references",
        references_command="validate",
        resource_root=Path("resources"),
        config_root=CONFIGS,
        as_json=True,
        bundle_id=REFERENCE_BUNDLE_ID,
        genome_build="GRCh38",
    )
    with pytest.raises(ValueError, match="mutually exclusive"):
        handle_references_command(args)
```

Also retain and re-run the existing test proving unscoped `validate` checks both builds.

- [ ] **Step 3: Run RED tests**

```bash
python -m pytest tests/test_grch37_resource_bootstrap.py tests/test_resource_bootstrap.py -q
```

Expected: new tests fail because `genome_build` is not implemented.

- [ ] **Step 4: Add parser argument**

In `add_references_parser()`:

```python
validate.add_argument("bundle_id", nargs="?")
validate.add_argument(
    "--genome-build",
    choices=[GenomeBuild.GRCH37.value, GenomeBuild.GRCH38.value],
)
validate.add_argument("--json", action="store_true", dest="as_json")
```

- [ ] **Step 5: Implement build-scoped validation**

Add a focused helper:

```python
def _validate_build_scope(
    args: argparse.Namespace,
    installer: ReferenceBundleInstaller,
    build: GenomeBuild,
) -> tuple[tuple[BundleValidationReport, ...], list[dict[str, object]]]:
    registry = ResourceRegistry(args.resource_root, active_build=build)
    if not registry.profiles:
        return (), []
    reference_ids = sorted({profile.reference_bundle for profile in registry.profiles.values()})
    reports = tuple(installer.validate(bundle_id) for bundle_id in reference_ids)
    profile_status = _profile_status(registry, verify_checksums=True)
    return reports, profile_status
```

In the `validate` branch:

```python
build_value = getattr(args, "genome_build", None)
if args.bundle_id is not None and build_value is not None:
    raise ValueError("bundle_id and --genome-build are mutually exclusive")
if build_value is not None:
    build = GenomeBuild(build_value)
    reports, validation_profile_status = _validate_build_scope(args, installer, build)
    payload = {
        "genome_build": build.value,
        "references": [_report_payload(report) for report in reports],
        "profile_status": validation_profile_status,
    }
    if args.as_json:
        print(json.dumps(payload, indent=2))
    else:
        _print_reports(reports, as_json=False)
        for item in validation_profile_status:
            marker = "OK" if item["valid"] is True else "INVALID"
            detail = "" if item["valid"] is True else f"\t{item['error']}"
            print(f"profile\t{marker}\t{item['profile_id']}{detail}")
    if (
        not reports
        or not validation_profile_status
        or any(not report.valid for report in reports)
        or any(item["valid"] is not True for item in validation_profile_status)
    ):
        raise SystemExit(2)
    return True
```

Leave the existing bundle-specific and unscoped branches unchanged after this new branch.

- [ ] **Step 6: Verify backend GREEN**

```bash
python -m pytest tests/test_grch37_resource_bootstrap.py tests/test_resource_bootstrap.py -q
python -m ruff check src/ontseq_platform/resource_commands.py tests/test_grch37_resource_bootstrap.py tests/test_resource_bootstrap.py
python -m mypy src/ontseq_platform/resource_commands.py
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/ontseq_platform/resource_commands.py tests/test_grch37_resource_bootstrap.py tests/test_resource_bootstrap.py
git commit -m "feat(resources): add build-scoped full validation"
```

---

## Task 2: Desktop bridge for build-scoped validation

**Files:**
- Modify: `desktop/ONTSeq.Desktop/WslServiceLauncher.cs`
- Modify: `desktop/ONTSeq.Desktop.Tests/Program.cs`

**Interfaces:**

```csharp
public static IReadOnlyList<string> ResourceFamilyValidationArguments(
    string genomeBuild,
    string resourceRootWsl);

public Task<(bool Ok, string Detail)> ValidateResourceFamilyAsync(
    DesktopSettings settings,
    string genomeBuild,
    CancellationToken cancellationToken);
```

- [ ] **Step 1: Add RED argument-construction tests**

In `Program.cs`:

```csharp
AssertSequenceEqual(
    ["references", "validate", "--genome-build", "GRCh38", "--resource-root",
     DesktopSettings.DefaultResourceRootWsl, "--json"],
    WslServiceLauncher.ResourceFamilyValidationArguments(
        "GRCh38", DesktopSettings.DefaultResourceRootWsl),
    "build-scoped validation uses structured argv and no shell interpolation");
AssertThrows<ArgumentOutOfRangeException>(() =>
    WslServiceLauncher.ResourceFamilyValidationArguments("hg38", "/tmp/resources"),
    "unknown build is rejected before process launch");
```

- [ ] **Step 2: Verify RED**

```powershell
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

Expected: compile failure because the method does not exist.

- [ ] **Step 3: Implement structured argv builder**

```csharp
public static IReadOnlyList<string> ResourceFamilyValidationArguments(
    string genomeBuild,
    string resourceRootWsl)
{
    _ = ManagedResourceBundleIds(genomeBuild);
    return new[]
    {
        "references", "validate",
        "--genome-build", genomeBuild,
        "--resource-root", DesktopSettings.NormalizeResourceRootWsl(resourceRootWsl),
        "--json"
    };
}
```

- [ ] **Step 4: Implement strict response validation**

Add private response records in `WslServiceLauncher.cs`:

```csharp
private sealed record ResourceValidationReference(bool Valid);
private sealed record ResourceValidationProfile(
    [property: JsonPropertyName("profile_id")] string ProfileId,
    bool Valid);
private sealed record ResourceFamilyValidationResponse(
    [property: JsonPropertyName("genome_build")] string GenomeBuild,
    ResourceValidationReference[] References,
    [property: JsonPropertyName("profile_status")] ResourceValidationProfile[] ProfileStatus);
```

Then:

```csharp
public async Task<(bool Ok, string Detail)> ValidateResourceFamilyAsync(
    DesktopSettings settings,
    string genomeBuild,
    CancellationToken cancellationToken)
{
    var result = await RunWslAsync(
        settings.WslDistribution,
        BackendInvocation(settings,
            ResourceFamilyValidationArguments(genomeBuild, settings.ResourceRootWsl).ToArray()),
        cancellationToken);
    var detail = string.IsNullOrWhiteSpace(result.StdOut) ? result.StdErr.Trim() : result.StdOut.Trim();
    if (result.ExitCode != 0) return (false, detail);
    try
    {
        var payload = JsonSerializer.Deserialize<ResourceFamilyValidationResponse>(
            result.StdOut, JsonDefaults.Options);
        if (payload is null ||
            !string.Equals(payload.GenomeBuild, genomeBuild, StringComparison.Ordinal) ||
            payload.References.Length == 0 || payload.ProfileStatus.Length == 0 ||
            payload.References.Any(item => !item.Valid) ||
            payload.ProfileStatus.Any(item => !item.Valid))
            return (false, "Build-spezifische Ressourcenvalidierung lieferte keinen vollständigen PASS.");
        return (true, detail);
    }
    catch (JsonException error)
    {
        return (false, "Ungültige JSON-Antwort der Ressourcenvalidierung: " + error.Message);
    }
}
```

- [ ] **Step 5: Add response-parser unit seam**

To avoid needing WSL in unit tests, extract:

```csharp
internal static (bool Ok, string Detail) InterpretResourceFamilyValidation(
    string requestedBuild,
    int exitCode,
    string stdout,
    string stderr)
```

Test nonzero exit, malformed JSON, wrong build, empty arrays, invalid reference and invalid profile all return `Ok == false`; a matching all-valid payload returns true.

- [ ] **Step 6: Verify GREEN**

```powershell
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

- [ ] **Step 7: Commit**

```bash
git add desktop/ONTSeq.Desktop/WslServiceLauncher.cs desktop/ONTSeq.Desktop.Tests/Program.cs
git commit -m "feat(desktop): bridge build-scoped resource validation"
```

---

## Task 3: Acceptance result model and persistent evidence

**Files:**
- Create: `desktop/ONTSeq.Desktop/SystemAcceptanceModels.cs`
- Create: `desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs`
- Modify: `desktop/ONTSeq.Desktop.Tests/Program.cs`

**Interfaces:**

```csharp
public enum AcceptanceVerdict { PASS, FAIL, CANCELLED }
public enum AcceptanceStepStatus { PASS, FAIL, NOT_RUN }
public enum ModuleReadinessState { READY, NOT_CONFIGURED, GATED, NOT_INCLUDED }

public sealed record SystemAcceptanceStep(
    string StepId,
    AcceptanceStepStatus Status,
    string Summary,
    string TechnicalDetail,
    bool AutomaticActionTaken,
    DateTimeOffset StartedAt,
    DateTimeOffset FinishedAt);

public sealed record ModuleReadiness(
    string ModuleId,
    ModuleReadinessState State,
    string Detail);

public sealed record SystemAcceptanceReport(
    string SchemaVersion,
    AcceptanceVerdict Verdict,
    string GenomeBuild,
    string ResourceRootWsl,
    DateTimeOffset StartedAt,
    DateTimeOffset FinishedAt,
    IReadOnlyList<SystemAcceptanceStep> Steps,
    IReadOnlyList<ModuleReadiness> Modules,
    string? SelfTestEvidenceDirectory,
    string ResearchUseOnlyNotice);
```

Stable step IDs are exactly:

```csharp
public static class SystemAcceptanceStepIds
{
    public const string Wsl = "wsl";
    public const string Runtime = "runtime";
    public const string ResourceRoot = "resource_root";
    public const string ResourceFamily = "resource_family";
    public const string ResourceValidation = "resource_validation";
    public const string SystemSelfTest = "system_self_test";
    public const string SelfTestEvidence = "self_test_evidence";
}
```

- [ ] **Step 1: Write RED serialization and invariant tests**

`SystemAcceptanceContractTests.RunAsync()` must assert:

```csharp
var report = SystemAcceptanceFixtures.Pass("GRCh38", "/tmp/resources");
var json = JsonSerializer.Serialize(report, JsonDefaults.OptionsIndented);
AssertContains(json, "\"verdict\": \"PASS\"", "verdict serializes as a string");
AssertEqual(7, report.Steps.Count, "all required steps are represented");
AssertEqual("aligned_bam_core", report.Modules[0].ModuleId,
    "core readiness has a stable module id");
```

Use `[JsonConverter(typeof(JsonStringEnumConverter))]` on each enum so evidence is human-readable and stable.

- [ ] **Step 2: Verify RED**

```powershell
dotnet build desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

Expected: missing types.

- [ ] **Step 3: Implement models and summary**

Add:

```csharp
public string ToOperatorSummary() =>
    $"ONTSeq Engineering Acceptance: {Verdict}\n" +
    $"Selected build: {GenomeBuild}\n" +
    $"Ready for authorized research aligned-BAM analysis: {(Verdict == AcceptanceVerdict.PASS ? "YES" : "NO")}\n" +
    "Analytical/clinical validation: NO\n";
```

Require constructor/factory validation that a PASS report has all seven required steps and all are PASS. A FAIL/CANCELLED report may have downstream NOT_RUN steps.

- [ ] **Step 4: Add evidence writer**

In the same file:

```csharp
public sealed class SystemAcceptanceEvidenceWriter
{
    private readonly string _root;
    public SystemAcceptanceEvidenceWriter(string root) => _root = Path.GetFullPath(root);

    public string Write(SystemAcceptanceReport report)
    {
        Directory.CreateDirectory(_root);
        var stamp = report.FinishedAt.UtcDateTime.ToString("yyyyMMdd_HHmmss_fff'Z'");
        var directory = Path.Combine(_root, stamp + "_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(directory);
        var jsonPath = Path.Combine(directory, "acceptance.json");
        var textPath = Path.Combine(directory, "acceptance.txt");
        WriteAtomic(jsonPath, JsonSerializer.Serialize(report, JsonDefaults.OptionsIndented) + Environment.NewLine);
        WriteAtomic(textPath, report.ToOperatorSummary());
        return directory;
    }
}
```

`WriteAtomic` writes to an adjacent `.tmp.<guid>` file and moves it into place. Production root is supplied as `%LOCALAPPDATA%\ONTSeq\acceptance`; tests use a temporary directory.

- [ ] **Step 5: Verify evidence writer**

Tests require both files, JSON round-trip, no output outside the injected root and no patient/sample payload in the summary.

```powershell
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

- [ ] **Step 6: Commit**

```bash
git add desktop/ONTSeq.Desktop/SystemAcceptanceModels.cs desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs desktop/ONTSeq.Desktop.Tests/Program.cs
git commit -m "feat(desktop): add acceptance evidence contract"
```

---

## Task 4: Acceptance service abstraction and production adapter

**Files:**
- Create: `desktop/ONTSeq.Desktop/SystemAcceptanceService.cs`
- Modify: `desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs`

**Interfaces:**

```csharp
public interface ISystemAcceptanceService
{
    Task<(bool Ok, string Detail)> CheckWslAsync(CancellationToken cancellationToken);
    Task<(bool Ok, string Detail)> CheckBackendAsync(CancellationToken cancellationToken);
    Task<string> InstallBundledRuntimeAsync(string runtimeArchiveWindows, CancellationToken cancellationToken);
    Task SaveResourceRootAsync(string resourceRootWsl, CancellationToken cancellationToken);
    Task<IReadOnlyDictionary<string, ResourceFamilyState>> CheckResourceFamiliesAsync(CancellationToken cancellationToken);
    Task<string> InstallResourceFamilyAsync(string genomeBuild, CancellationToken cancellationToken);
    Task<string> RepairResourceFamilyAsync(string genomeBuild, CancellationToken cancellationToken);
    Task<(bool Ok, string Detail)> ValidateResourceFamilyAsync(string genomeBuild, CancellationToken cancellationToken);
    Task<string> RunSelfTestAsync(CancellationToken cancellationToken);
    string? ConfiguredModkitExecutableWsl { get; }
}
```

Production adapter:

```csharp
public sealed class WslSystemAcceptanceService : ISystemAcceptanceService
{
    private readonly WslServiceLauncher _launcher;
    private readonly DesktopSettings _settings;
}
```

- [ ] **Step 1: Add RED adapter-contract tests**

The fake service used in tests must record operation names. Verify the interface makes no shell-string or bioinformatics argument part of the runner API.

- [ ] **Step 2: Implement adapter**

Mappings are exact:

```csharp
public Task<(bool Ok, string Detail)> CheckWslAsync(CancellationToken token) =>
    _launcher.CheckWslAsync(_settings, token);

public Task<(bool Ok, string Detail)> CheckBackendAsync(CancellationToken token) =>
    _launcher.CheckBackendAsync(_settings, token);

public async Task SaveResourceRootAsync(string root, CancellationToken token)
{
    token.ThrowIfCancellationRequested();
    _settings.ResourceRootWsl = DesktopSettings.NormalizeResourceRootWsl(root);
    _settings.SaveUserSettings();
    await Task.CompletedTask;
}
```

The remaining methods delegate to existing launcher methods including the new `ValidateResourceFamilyAsync`.

- [ ] **Step 3: Verify**

```powershell
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

- [ ] **Step 4: Commit**

```bash
git add desktop/ONTSeq.Desktop/SystemAcceptanceService.cs desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs
git commit -m "feat(desktop): add acceptance service boundary"
```

---

## Task 5: Orchestration runner — happy path, repair path and cross-build isolation

**Files:**
- Create: `desktop/ONTSeq.Desktop/SystemAcceptanceRunner.cs`
- Modify: `desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs`

**Interfaces:**

```csharp
public sealed record SystemAcceptanceRequest(
    string GenomeBuild,
    string ResourceRootWsl,
    string RuntimeArchiveWindows);

public sealed record SystemAcceptanceProgress(string StepId, string Message);

public sealed record SystemAcceptanceExecution(
    SystemAcceptanceReport Report,
    string EvidenceDirectory);

public sealed class SystemAcceptanceRunner
{
    public SystemAcceptanceRunner(
        ISystemAcceptanceService service,
        SystemAcceptanceEvidenceWriter evidenceWriter,
        TimeProvider? timeProvider = null);

    public Task<SystemAcceptanceExecution> RunAsync(
        SystemAcceptanceRequest request,
        IProgress<SystemAcceptanceProgress>? progress,
        CancellationToken cancellationToken);
}
```

- [ ] **Step 1: Write RED happy-path test**

Fake sequence:

```text
CheckWsl -> PASS
CheckBackend -> FAIL
InstallRuntime
CheckBackend -> PASS
SaveResourceRoot
CheckResourceFamilies -> GRCh38 NotInstalled
InstallResourceFamily GRCh38
CheckResourceFamilies -> GRCh38 Ready
ValidateResourceFamily GRCh38 -> PASS
CheckResourceFamilies -> GRCh38 Ready
RunSelfTest -> synthetic evidence folder with PASS report
```

Assert overall PASS and exact operation order. Assert runtime/resource `AutomaticActionTaken == true`.

- [ ] **Step 2: Write already-ready and repair RED tests**

Ready system must not call install/repair. Incomplete selected family must call repair exactly once. Include a deliberately broken non-selected family in the fake dictionary and prove GRCh38 acceptance still passes when GRCh37 is broken.

- [ ] **Step 3: Verify RED**

```powershell
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

- [ ] **Step 4: Implement sequential fail-closed runner**

The algorithm is exactly:

```text
1 wsl
2 runtime check; install only when needed; re-check
3 normalize + persist resource root
4 resource-family status; install/repair only selected family; re-check
5 build-scoped full validation; re-check selected fast status
6 run existing system self-test
7 independently validate self-test evidence
8 build module-readiness table and persist acceptance evidence
```

No step starts if the prior required step failed.

- [ ] **Step 5: Implement optional-module rows**

For this work package:

```csharp
new("aligned_bam_core",
    overallPass ? ModuleReadinessState.READY : ModuleReadinessState.GATED,
    overallPass ? "Selected aligned-BAM Core scope passed workstation acceptance."
                : "Aligned-BAM Core acceptance has not passed."),
new("methylation_modkit",
    string.IsNullOrWhiteSpace(_service.ConfiguredModkitExecutableWsl)
        ? ModuleReadinessState.NOT_CONFIGURED
        : ModuleReadinessState.GATED,
    string.IsNullOrWhiteSpace(_service.ConfiguredModkitExecutableWsl)
        ? "No separate modkit executable is configured."
        : "A separate modkit path is configured; methylation remains outside Core acceptance."),
new("marlin", ModuleReadinessState.GATED,
    "MARLIN requires locked model/resources and frozen numerical runtime compatibility."),
new("raw_signal_dorado", ModuleReadinessState.NOT_INCLUDED,
    "Raw-signal/Dorado is not part of this Desktop acceptance scope.")
```

Do not report methylation as READY from this Core acceptance run.

- [ ] **Step 6: Verify GREEN**

```powershell
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

- [ ] **Step 7: Commit**

```bash
git add desktop/ONTSeq.Desktop/SystemAcceptanceRunner.cs desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs
git commit -m "feat(desktop): orchestrate one-click engineering acceptance"
```

---

## Task 6: Failure, cancellation and self-test evidence gates

**Files:**
- Modify: `desktop/ONTSeq.Desktop/SystemAcceptanceRunner.cs`
- Modify: `desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs`

**Interfaces:**

```csharp
internal static (bool Ok, string Detail) VerifySelfTestEvidence(string directory);
```

- [ ] **Step 1: Write RED evidence tests**

Create four temporary folders:

1. no `system-smoke.report.json` -> fail;
2. malformed JSON -> fail;
3. `{"schema_version":"0.1.0","verdict":"FAIL"}` -> fail;
4. a minimal valid object containing `schema_version=0.1.0`, `verdict=PASS`, non-empty `checks`, `output_paths` and `limitations` -> pass.

Use `JsonDocument`; do not deserialize into a loose dynamic object.

- [ ] **Step 2: Write RED operational failure tests**

Required cases:

```text
WSL unavailable -> FAIL at wsl; no mutation calls
runtime install then backend still invalid -> FAIL at runtime
resource install/repair then selected family still non-ready -> FAIL at resource_family
build-scoped validation nonzero/malformed -> FAIL at resource_validation
self-test throws -> FAIL at system_self_test
self-test returns folder but evidence invalid -> FAIL at self_test_evidence
```

Every later required step must remain `NOT_RUN`.

- [ ] **Step 3: Write RED cancellation test**

Cancel during fake resource installation. Assert:

```csharp
AssertEqual(AcceptanceVerdict.CANCELLED, execution.Report.Verdict, "cancelled run stays distinct from FAIL");
AssertTrue(execution.Report.Steps.Skip(cancelledIndex + 1)
    .All(step => step.Status == AcceptanceStepStatus.NOT_RUN),
    "downstream steps are not fabricated after cancellation");
```

Acceptance evidence still has to be written for CANCELLED.

- [ ] **Step 4: Implement evidence parser**

Require:

```csharp
var reportPath = Path.Combine(directory, "system-smoke.report.json");
if (!File.Exists(reportPath)) return (false, "system-smoke.report.json fehlt.");
using var document = JsonDocument.Parse(File.ReadAllText(reportPath));
var root = document.RootElement;
if (root.GetProperty("schema_version").GetString() != "0.1.0") return (false, "Unbekanntes System-Smoke-Schema.");
if (root.GetProperty("verdict").GetString() != "PASS") return (false, "System-Smoke meldet keinen PASS.");
if (root.GetProperty("checks").GetArrayLength() < 1) return (false, "System-Smoke enthält keine Prüfnachweise.");
if (root.GetProperty("output_paths").ValueKind != JsonValueKind.Object) return (false, "System-Smoke output_paths fehlt.");
return (true, reportPath);
```

Catch `JsonException`, `KeyNotFoundException` and `InvalidOperationException` and convert them to fail-closed detail.

- [ ] **Step 5: Implement failure/cancellation finalization**

All exits use one finalization path that:

1. sets overall verdict;
2. fills untouched required steps with NOT_RUN;
3. generates module states without overclaiming;
4. persists `acceptance.json` and `acceptance.txt`;
5. returns the evidence directory.

Unexpected exceptions become FAIL at the currently active step, except `OperationCanceledException` with the supplied token, which becomes CANCELLED.

- [ ] **Step 6: Verify all contract tests**

```powershell
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

- [ ] **Step 7: Commit**

```bash
git add desktop/ONTSeq.Desktop/SystemAcceptanceRunner.cs desktop/ONTSeq.Desktop.Tests/SystemAcceptanceContractTests.cs
git commit -m "test(desktop): harden acceptance failures and evidence gate"
```

---

## Task 7: One-click SetupWindow integration

**Files:**
- Modify: `desktop/ONTSeq.Desktop/SetupWindow.xaml`
- Modify: `desktop/ONTSeq.Desktop/SetupWindow.xaml.cs`
- Modify: `desktop/ONTSeq.Desktop.Tests/Program.cs`

**Interfaces:**
- UI action name: `Vollständig einrichten & testen`.
- Handler: `Acceptance_Click`.
- Status control: `AcceptanceStatusText`.

- [ ] **Step 1: Add XAML contract RED checks**

In `Program.cs`, read the XAML source relative to the checkout and assert the button content/name/handler exist. This is a lightweight structural guard in addition to WPF compilation:

```csharp
var setupXaml = File.ReadAllText(Path.Combine(repoRoot,
    "desktop", "ONTSeq.Desktop", "SetupWindow.xaml"));
AssertEqual("True", setupXaml.Contains("Vollständig einrichten &amp; testen", StringComparison.Ordinal).ToString(),
    "setup exposes one-click acceptance action");
AssertEqual("True", setupXaml.Contains("Click=\"Acceptance_Click\"", StringComparison.Ordinal).ToString(),
    "acceptance action is wired to its handler");
```

- [ ] **Step 2: Add XAML row**

Increase the inner setup grid row count and add, below the existing self-test row:

```xml
<TextBlock Grid.Row="7" Grid.Column="0" Text="Engineering-Acceptance"
           FontWeight="SemiBold" Margin="0,12,0,8" />
<TextBlock Grid.Row="7" Grid.Column="1" x:Name="AcceptanceStatusText"
           Text="Noch nicht ausgeführt" Foreground="#5F6B7A"
           Margin="0,12,12,8" TextWrapping="Wrap" />
<Button Grid.Row="7" Grid.Column="2" x:Name="AcceptanceButton"
        Content="Vollständig einrichten &amp; testen"
        Click="Acceptance_Click" MinWidth="220" Margin="0,8,0,4" />
```

Move the existing `DetailText` border from row 7 to row 8 and add the ninth `RowDefinition`.

- [ ] **Step 3: Implement confirmation and runner invocation**

Handler starts with a confirmation containing all six points from the spec. On confirmation:

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
        var root = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
            "ONTSeq", "acceptance");
        var runner = new SystemAcceptanceRunner(service, new SystemAcceptanceEvidenceWriter(root));
        var archive = Path.Combine(AppContext.BaseDirectory, "runtime", "ontseq-linux-runtime.tar.gz");
        var progress = new Progress<SystemAcceptanceProgress>(item =>
        {
            AcceptanceStatusText.Text = "● " + item.Message;
            DetailText.Text = item.Message;
        });
        var execution = await runner.RunAsync(
            new SystemAcceptanceRequest(build, _settings.ResourceRootWsl, archive),
            progress,
            token);
        AcceptanceStatusText.Text = execution.Report.Verdict == AcceptanceVerdict.PASS
            ? "✓ PASS"
            : execution.Report.Verdict == AcceptanceVerdict.CANCELLED ? "— CANCELLED" : "✕ FAIL";
        DetailText.Text = execution.Report.ToOperatorSummary() +
            "Prüfnachweise: " + execution.EvidenceDirectory;
    }, "Engineering-Acceptance konnte nicht abgeschlossen werden.");

    await RefreshAsync();
}
```

- [ ] **Step 4: Wire busy state**

Add:

```csharp
AcceptanceButton.IsEnabled = !busy;
```

Do not enable an independent cancellation button; window close continues to cancel via `_cts`.

- [ ] **Step 5: Build/test WPF**

```powershell
dotnet build desktop/ONTSeq.Desktop/ONTSeq.Desktop.csproj -c Release
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add desktop/ONTSeq.Desktop/SetupWindow.xaml desktop/ONTSeq.Desktop/SetupWindow.xaml.cs desktop/ONTSeq.Desktop.Tests/Program.cs
git commit -m "feat(desktop): add one-click engineering acceptance UI"
```

---

## Task 8: Operator docs, bundle contract and complete verification

**Files:**
- Modify: `docs/DESKTOP_FIRST_RUN.md`
- Modify: `desktop/README-FIRST-RUN.md`
- Modify: `.github/workflows/desktop-ci.yml`
- Optional changelog entry in: `CHANGELOG.md`

**Interfaces:**
- Acceptance evidence root: `%LOCALAPPDATA%\ONTSeq\acceptance\<UTC timestamp>_<run id>\`.
- Main operator action: `Vollständig einrichten & testen`.

- [ ] **Step 1: Update both first-run documents**

The preferred path becomes:

```text
1. WSL2 + Ubuntu must already start.
2. Extract the complete engineering ZIP.
3. Start ONTSeq.Desktop.exe.
4. Open System einrichten.
5. Select GRCh37 or GRCh38.
6. Press Vollständig einrichten & testen.
7. Require Engineering Acceptance: PASS.
8. Only then use authorized research aligned-BAM data.
```

Explicitly state that no PowerShell/CMD/WSL/Python/R/Conda command is expected for normal setup and that any such requirement during a clean-room run is a product defect.

- [ ] **Step 2: Update local-files documentation**

Add:

```text
Engineering acceptance evidence:
%LOCALAPPDATA%\ONTSeq\acceptance\
```

State that it contains small JSON/text evidence and references, rather than copies, the self-test output.

- [ ] **Step 3: Strengthen Desktop bundle CI documentation check**

In `Verify complete first-run bundle`, after checking the first-run file exists:

```powershell
if (-not (Select-String -Path $firstRun -SimpleMatch "Vollständig einrichten & testen" -Quiet)) {
  throw "First-run documentation does not describe one-click engineering acceptance"
}
if (-not (Select-String -Path $firstRun -SimpleMatch "Engineering Acceptance" -Quiet)) {
  throw "First-run documentation does not state the acceptance verdict contract"
}
```

No additional runtime package is bundled by this task.

- [ ] **Step 4: Run complete local/static verification**

```bash
python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m mypy src
```

```powershell
dotnet build desktop/ONTSeq.Desktop/ONTSeq.Desktop.csproj -c Release
dotnet run --project desktop/ONTSeq.Desktop.Tests/ONTSeq.Desktop.Tests.csproj -c Release
```

Expected: all PASS.

- [ ] **Step 5: Commit docs/CI**

```bash
git add docs/DESKTOP_FIRST_RUN.md desktop/README-FIRST-RUN.md .github/workflows/desktop-ci.yml CHANGELOG.md
git commit -m "docs(desktop): document clean-room one-click acceptance"
```

If `CHANGELOG.md` was not changed, omit it from `git add` rather than creating an empty/noise edit.

- [ ] **Step 6: Require fresh GitHub gates on the final head**

Verify the final commit has:

```text
CI                       success
Desktop bundle contract  success
Desktop CI               success
MARLIN runtime smoke     success or unchanged/not-required by path filter
```

Within Desktop CI, require both `linux-runtime` and `windows-desktop` jobs to complete successfully, including stock-Ubuntu runtime smoke and the Windows Desktop test executable.

- [ ] **Step 7: Manual clean-room workstation qualification**

On a fresh Windows 11 VM or clean workstation containing only functional WSL2 + Ubuntu:

```text
A. Download the final complete Desktop engineering artifact.
B. Extract the entire artifact.
C. Start ONTSeq.Desktop.exe.
D. Open System einrichten.
E. Select GRCh38 for the first run.
F. Press only Vollständig einrichten & testen.
G. Do not open PowerShell, CMD, WSL shell, Python, R, Conda or an editor.
H. Require final PASS.
I. Open acceptance.json and acceptance.txt from the displayed evidence directory.
J. Confirm all seven required steps are PASS.
K. Confirm aligned_bam_core=READY, marlin=GATED, raw_signal_dorado=NOT_INCLUDED, and methylation_modkit is not promoted by Core acceptance.
L. Repeat acceptance without changing the machine and confirm the valid runtime/resources are reused rather than reinstalled.
M. Run one synthetic/non-patient aligned BAM through the normal Desktop analysis path and open HTML/XLSX/JSON.
```

Any manual dependency install, environment-variable edit, resource-file browsing, ad-hoc path patch, or shell troubleshooting is a **FAIL of the clean-room goal** and must become a new defect before the workflow is described as zero-troubleshooting.

- [ ] **Step 8: Final review boundary**

Keep PR #77 Draft until:

```text
- all automated tests/gates are green;
- the clean-room workstation run has PASS evidence;
- the PR diff has been independently reviewed for overclaiming and cross-build leakage;
- no patient/research payload, credentials, external model artifact or large validation data has entered Git.
```

Do not merge to `main` as part of this task without explicit user approval.

---

## Self-Review Results

### Spec coverage

- One-click orchestration: Tasks 5 and 7.
- Runtime reuse/install/re-check: Tasks 4 and 5.
- Build-scoped full resource validation including Knowledge/Panel: Tasks 1 and 2.
- Cross-build isolation: Tasks 1, 5 and 6.
- Existing canonical system smoke reuse: Tasks 4–6.
- Independent self-test evidence verification: Task 6.
- PASS/FAIL/CANCELLED + stable step IDs: Task 3.
- Persistent JSON/text evidence: Task 3.
- Optional module separation: Task 5.
- UI remediation/no terminal requirement: Task 7.
- Clean-room real workstation gate: Task 8.
- No analytical/clinical overclaim: Global constraints, Tasks 3, 5, 8.

### Placeholder scan

No `TBD`, `TODO`, “implement later”, unspecified test cases or undefined neighbor interfaces remain in this plan.

### Type consistency

The following names are canonical across tasks and must not be renamed independently during implementation:

```text
ISystemAcceptanceService
WslSystemAcceptanceService
SystemAcceptanceRunner
SystemAcceptanceRequest
SystemAcceptanceProgress
SystemAcceptanceExecution
SystemAcceptanceReport
SystemAcceptanceStep
SystemAcceptanceEvidenceWriter
AcceptanceVerdict
AcceptanceStepStatus
ModuleReadinessState
ValidateResourceFamilyAsync
ResourceFamilyValidationArguments
```
