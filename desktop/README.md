# ONTSeq Desktop

Windows operator surface for the ONTSeq execution core. The desktop app does not reimplement bioinformatics; it launches the bundled Linux runtime through WSL2, talks to the loopback service, and opens the same persisted HTML/XLSX/result artifacts produced by the canonical backend.

> **Current engineering build: Desktop/Core v0.4.1. Research Use Only, unsigned, not clinically validated.**

## v0.4.1 user path

1. Extract the complete engineering bundle and start `ONTSeq.Desktop.exe`.
2. Open **System einrichten** on first use and install the bundled Linux runtime into WSL. The bundle contains Python, R, samtools, Cramino, Sniffles2, Mosdepth, minimap2, QDNAseq, QDNAseq.hg19, QDNAseq.hg38, ACE and the ONTSeq Core/resources.
3. Run **Selbsttest starten**. This is a full installed-runtime system smoke, not only a package import check.
4. Keep the user-writable WSL resource root at `~/.local/share/ontseq/resources` or configure another absolute WSL path. Use **Installieren** or **Reparieren** to provision and validate the complete resource family: `GRCh38_GENCODE50_MANE1.5_v1`, `HEMATOLOGY_v1`, `AML_AS_111_GRCh38_v1` and both profiles. Repair safely replaces divergent pinned resources through the backend transaction; manual deletion is not required.
5. Select `AML_LCWGS_GRCh38` (default) or `AML_AS_111_GRCh38`. Both resolve the pinned reference, annotation and knowledge bundles automatically; the AS profile additionally resolves `AML_AS_111_GRCh38_v1`. There is no GRCh37 fallback or cross-build resource selection in this work package.
6. Select one aligned `.bam` file. Both `sample.bam.bai` and `sample.bai` are accepted, with deterministic precedence for `sample.bam.bai` when both exist. The declared index is checked explicitly rather than rediscovered implicitly by samtools.
7. Sample ID is suggested from the filename and remains editable.
8. Select one of the two explicit GRCh38 analysis profiles. The expected build is shown immediately; the backend verifies the BAM sequence dictionary and the desktop replaces the display with the build persisted in run provenance when available.
9. Click **ANALYSE STARTEN**.
10. The application checks WSL, backend capabilities, input/index, resource root, required profile bundles and output location before starting the local service.
11. It starts the canonical `ontseq serve` backend, submits the run, polls progress, and enables the HTML, XLSX and result-folder buttons when outputs are present.

The UI keeps the backend's distinct stage outcomes (`COMPLETED`, `NO_CALL`, `FAILED`, `NOT_RUN`) instead of flattening them into done/pending. A technical `PASS` is not a biological negative result and is not clinical validation.

### Reference and resume integrity in v0.4.1

New profile runs resolve a manifest-driven `ReferenceBundle`, `PanelBundle` and `KnowledgeBundle` beneath `resourceRootWsl`. Bundle installation and repair remain backend operations; Desktop only invokes their typed CLI contracts and displays the result. The selected profile and BAM dictionary must both resolve to GRCh38.

The former explicit `referenceLocksWsl`, `adaptiveTargetBedWsl` and `adaptiveTargetBedVersion` settings remain deserializable for one compatibility release. They are labelled as legacy in setup and are never mixed into the new bundle-backed profiles.

Aligned-BAM intake signs the BAM, the explicitly declared BAI, manifest reference/build claims, the complete Reference-Lock identity, Core version and git commit into the stage resume signature. BAM and BAI are hashed before intake and again after validation; a change during the intake window fails the stage rather than permitting reuse of stale results. The packed runtime also contains the exact build git commit so Desktop-started runs can record provenance rather than `UNKNOWN`.

### Full installed-runtime self-test

`ontseq system-smoke` deliberately uses two independent deterministic synthetic fixtures so that one artificial signal does not distort the other analytical lane:

- a long-read fixture exercises BAM intake, samtools, Cramino, Sniffles2 and basic report rendering;
- a whole-genome-shaped GRCh37/hg19 read-depth fixture exercises the canonical run envelope with live QDNAseq + ACE. Baseline autosomes have two synthetic read starts per 100 kbp, chromosome 7 has half-depth and chromosome 8 has 1.5x depth.

A PASS requires the CNV path to recover chromosome 7 at approximately copy number 1 and chromosome 8 at approximately copy number 3, with agreement across the configured 100/500/1000-kbp QDNAseq resolutions. It also requires the CNV result to be present in normalized JSON, provenance, HTML and the `CNV Fits`, `CNV Consensus` and `CNV Segments` Excel sheets. Release SHA256 checksums are independently recalculated, the identical canonical run is executed a second time, the CNV stage must resume content-addressed, and release checksums must still verify afterward.

GRCh38 has a separate canonical real-tool CI gate using a whole-genome-shaped GRCh38 fixture with the same deterministic truth pattern. That gate must load the 100/500/1000-kbp hg38 annotations, recover chromosome 7 at CN~1 and chromosome 8 at CN~3, generate the normal report artifacts, verify release checksums and prove content-addressed resume.

The packaging workflow additionally relocates the packed Linux runtime and boots it inside a stock Ubuntu 24.04 container with no host Python or R. The packed runtime is checked for both build-specific QDNAseq annotation packages before the self-contained Windows x64 executable is bundled.

The desktop setup invokes the installed-runtime system smoke through WSL and stores `system-smoke.report.json`, the canonical synthetic run envelope and a text log under `%LOCALAPPDATA%\ONTSeq\self-test\<timestamp>`.

### Current CNV behavior

The bundled reproducible QDNAseq environment still contains the earlier hg19 package for
compatibility, but the two published Desktop profiles in this work package request only the
GRCh38 QDNAseq + ACE lane. They use the open-source `QDNAseq.hg38` 1.2.0 data package pinned to
upstream commit `cf7c07e39de0ac64a9c38cb030cba4626e2aae83`, with 100/500/1000-kbp views, a 500-kbp
primary view and multi-bin consensus. No hg19 resource is selected as a profile fallback.

The GRCh38 lane was enabled only after a canonical real-tool test recovered a deterministic synthetic chromosome 7 loss (CN~1) and chromosome 8 gain (CN~3), produced HTML/XLSX/JSON outputs, verified release checksums and resumed content-addressed on an unchanged rerun. This is engineering verification, not cohort-level analytical or clinical validation.

ACE cellularity/purity remains a fitted model parameter rather than a measured tumour/blast fraction. The configured warning boundaries and the whole-chromosome classification fraction are engineering defaults and are surfaced as limitations; they are not validated clinical cut-offs.

## Why .NET 10 + WPF

The desktop is a small Windows-native shell around the Linux bioinformatics backend. Filesystem integration, process control and a single self-contained Windows executable are first-class requirements; the analytical implementation stays independently testable in Python/R/Linux.

## Configuration

The app loads the first existing settings file from:

- `%LOCALAPPDATA%\ONTSeq\desktop.settings.json`
- `%PROGRAMDATA%\ONTSeq\desktop.settings.json`

See `desktop.settings.example.json` for the available keys. First-run setup writes the user-local settings automatically. A runtime installed by the v0.4.1 Desktop is stored under `~/.local/share/ontseq/runtime-v0.4.1` in the selected WSL distribution.

The profile-resource defaults are:

```json
{
  "resourceRootWsl": "~/.local/share/ontseq/resources",
  "defaultProfile": "AML_LCWGS_GRCh38"
}
```

### Network drives

A Windows path such as `P:\Lab\run.bam` is translated to `/mnt/p/Lab/run.bam`. This works only when that drive is visible inside the selected WSL distribution. The desktop preflight checks reachability before starting a run. If a mapped network drive is not available as `/mnt/<drive>`, mount it deliberately through WSL/DrvFs first; UNC paths are not guessed or silently remapped.

## Build

```powershell
dotnet build desktop/ONTSeq.Desktop/ONTSeq.Desktop.csproj -c Release
```

Self-contained single-file Windows publish:

```powershell
dotnet publish desktop/ONTSeq.Desktop/ONTSeq.Desktop.csproj -c Release -r win-x64 --self-contained true -p:PublishSingleFile=true
```

`.github/workflows/desktop-ci.yml` builds a relocatable Linux runtime from the pinned QDNAseq environment, proves the packed environment after relocation, runs the full system smoke in stock Ubuntu, then builds/tests the Windows WPF shell and creates the versioned engineering bundle. The final main-branch artifact should be preferred over a pre-merge PR artifact because its embedded runtime git commit then identifies the exact merged source revision.

## Deliberate v0.4.1 limitation: cancellation

The **Abbrechen** control remains disabled. A backend cancellation contract must first guarantee that an interrupted external tool cannot leave a stage looking complete. Closing a progress view is not treated as cancellation.

## Remaining workstation / validation gates

- first real Windows + WSL BAM run on the target workstation, including HTML/XLSX opening from the desktop;
- mapped-network-drive verification on the actual workstation if the BAM resides on `P:` or another mapped drive;
- service-level Windows smoke covering Desktop Start → loopback API → persisted outputs with a real `wsl.exe` host rather than Linux-container emulation;
- analytical validation of GRCh38 CNV on controlled real/reference samples for the active profiles;
- backend-owned cooperative cancellation with explicit interrupted-state semantics;
- analytical validation on controlled real/reference samples before any diagnostic use.

See also `README-FIRST-RUN.md` and `ONTSeq.Desktop/CHANGELOG.md`.
