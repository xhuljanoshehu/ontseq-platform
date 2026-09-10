# Isolated engineering test

`ONTSEQ_DESKTOP_SETTINGS` selects an absolute local `.json` settings file for that
Desktop process. Both reading and saving use only this file. A missing override file
starts with default settings; it does **not** read the existing user or machine-wide
configuration. Relative paths, UNC paths, control characters, and non-JSON extensions
are rejected. The value is a filesystem path, never a shell command.

Example in PowerShell, using a dedicated local test folder:

```powershell
$env:ONTSEQ_DESKTOP_SETTINGS = 'C:\ONTSeq-Test\desktop.settings.json'
Start-Process -FilePath 'C:\ONTSeq-Test\ONTSeq.Desktop.exe'
Remove-Item Env:ONTSEQ_DESKTOP_SETTINGS
```

The child Desktop inherits the override. Only that process's setup writes to the test
settings file. Closing the launching shell or removing the variable does not change an
already running Desktop. Without the variable, the original settings lookup is unchanged.

Fresh settings default to `~/.local/share/ontseq/resources-v0.8.1`, with the suffix derived
from the Desktop/Core version. An existing settings file keeps its selected resource path.
This naming convention does not update bundle contents or migrate an older installation.
If a local resource directory is deliberately renamed, finish the move and full checksum
validation before selecting its new absolute WSL path; retain historical provenance and
access to the old location as described in [the resource system](../docs/REFERENCE_SYSTEM.md#resource-root).

For a prepared local test, the JSON may specify `runtimeBinWsl`, `backendCommand`,
`resourceRootWsl`, `outputDirectoryWindows`, `defaultProfile`, and a free loopback `port`.
Point `resourceRootWsl` at the existing separately validated reference directory rather
than duplicating reference data. `AML_LCWGS_GRCh37` requires the exact full GENCODE 19
GRCh37.p13 assembly; `AML_AS_111_GRCh37` uses the same dictionary contract plus the pinned
GRCh37-specific Adaptive-Sampling panel. `AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25` and
`AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25` separately require exactly chr1–22/X/Y plus
`chrM=16571` in UCSC hg19 order. These contracts are not interchangeable, and neither AS
profile falls back to the GRCh38 panel. The GRCh37 panel contains 110/111 mapped selection
intervals, 109 reciprocal-exact intervals and 108 native GENCODE-19 ROIs; ACACA requires
roundtrip review, while CT45A2, IGH and GPR128 remain mapping/ROI gaps.

The runtime installer requires the unchanged base archive, matching versioned Core
wheel, and `runtime/SHA256SUMS`. It verifies both artifacts, creates a new
`~/.local/share/ontseq/runtime-v0.8.1-<unique-id>` directory, relocates the base runtime,
installs the wheel offline, and checks the exact Core version before saving the selected
settings file. Existing runtime directories are not removed or overwritten. Package
hashes are integrity records, not a publisher signature or clinical-validation claim.

Select the BAM before opening **Live-Workspace** to scope the service to its parent
directory. Without a selected BAM, the service is restricted to the configured output
directory. Opening the workspace does not start analysis. Keep Desktop open while using
the service; close only after running analyses finish. Changing service settings or its
input-root scope requires reopening Desktop and never silently restarts a live service.
