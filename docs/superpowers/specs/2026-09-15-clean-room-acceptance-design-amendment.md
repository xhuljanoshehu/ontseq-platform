# Clean-Room Acceptance Design Amendment — Build-Scoped Full Validation

**Amends:** `docs/superpowers/specs/2026-09-15-clean-room-acceptance-design.md`

This amendment replaces only the resource-validation mechanism described in the approved clean-room acceptance design. All other scope, evidence, UI, automation and release boundaries remain unchanged.

## Reason for amendment

The existing backend command `ontseq references validate <bundle_id>` performs full checksum validation for a reference bundle, but the selected ONTSeq build family also depends on profile-resolved Knowledge and optional Panel resources. Those are checked through `ResourceRegistry.resolve_profile(..., verify_files=True)`, not by passing their IDs to `ReferenceBundleInstaller.validate()`.

Therefore, validating the three Desktop-managed IDs individually would not prove the complete selected family and would give the acceptance runner the wrong abstraction.

## Correct backend contract

Extend the existing command to support exactly one optional build selector:

```text
ontseq references validate --genome-build GRCh37 --resource-root <root>
ontseq references validate --genome-build GRCh38 --resource-root <root>
```

`bundle_id` and `--genome-build` are mutually exclusive. Existing forms remain backward compatible:

```text
ontseq references validate <reference_bundle_id> --resource-root <root>
ontseq references validate --resource-root <root>
```

Build-scoped validation must:

1. create `ResourceRegistry(resource_root, active_build=<selected build>)`;
2. identify the selected build's installed reference bundle through the registry/profile context rather than guessing by filename;
3. checksum-validate the installed reference bundle with the existing `ReferenceBundleInstaller.validate()` path;
4. call `_profile_status(registry, verify_checksums=True)` for the selected build only;
5. fail with exit code 2 if no selected-build profile is resolvable, the selected reference bundle is invalid, or any selected-build profile fails checksum-backed resolution;
6. not instantiate or validate the non-selected build's registry;
7. preserve the existing unscoped `validate` semantics for callers that intentionally request all installed builds.

The JSON form must return enough typed information for Desktop to distinguish the selected build, validated reference reports and checksum-backed profile status without parsing human text.

## Desktop bridge

Add a build-scoped method rather than teaching the acceptance runner CLI syntax:

```csharp
public Task<(bool Ok, string Detail)> ValidateResourceFamilyAsync(
    DesktopSettings settings,
    string genomeBuild,
    CancellationToken cancellationToken)
```

The `WslServiceLauncher` implementation invokes:

```text
ontseq references validate --genome-build <GRCh37|GRCh38> --resource-root <normalized root> --json
```

and returns `Ok = true` only on exit code 0 with a valid response for the requested build.

`SystemAcceptanceRunner` calls this method exactly once after install/repair and fast status readiness. It then refreshes `CheckResourceFamiliesAsync` and again requires `CanAnalyze == true` for the selected build.

## Updated acceptance evidence

The `resource_validation` step records:

- selected genome build;
- normalized resource root;
- build-scoped validation result;
- backend detail/JSON summary;
- final selected-family fast status after validation.

It does not claim that the other build was checked.

## Updated tests

The implementation plan must include regression tests proving:

- `--genome-build GRCh37` validates only GRCh37 registry/profile resources with `verify_files=True`;
- `--genome-build GRCh38` validates only GRCh38 registry/profile resources with `verify_files=True`;
- a deliberately broken non-selected build cannot fail the selected build command;
- `bundle_id` plus `--genome-build` is rejected;
- legacy unscoped `references validate` still checks all installed builds;
- Desktop constructs the build-scoped command without shell interpolation;
- acceptance fails closed if the build-scoped command returns nonzero or malformed JSON.

## Updated acceptance criterion

Replace the original criterion that every managed bundle ID is validated individually with:

> The selected genome-build family passes one full build-scoped backend validation that checksum-validates its installed reference bundle and resolves every selected-build profile with `verify_files=True`, thereby covering its pinned Reference, Knowledge and optional Panel resources without making the non-selected build an implicit prerequisite.
