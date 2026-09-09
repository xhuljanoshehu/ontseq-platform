# ONTSeq Desktop

Windows operator surface for the ONTSeq execution core. The desktop app does not reimplement bioinformatics; it launches the bundled Linux runtime through WSL2, talks to the loopback service, and opens the same persisted HTML/XLSX/result artifacts produced by the canonical backend.

> **Current engineering build: Desktop/Core v0.7.1. Research Use Only, unsigned, not clinically validated.**

## v0.7.1 user path

The local engineering candidate adds build-isolated GRCh37.p13/hg19 lcWGS and Adaptive
Sampling profiles to the four existing GRCh38 profiles, plus a live, same-origin report workspace. This is an executable software
integration, not the earlier standalone HTML demonstration. It is not a published or
clinically validated release.

1. Keep the complete package together. The runtime directory contains a pinned base tool
   archive, the Core v0.7.1 wheel, and their SHA256SUMS. Setup verifies both inputs, unpacks
   the tools into a new `runtime-v0.7.1-<unique-id>` prefix, then installs the new Core
   offline and verifies its version. It does not replace the existing runtime directory.
   Fresh settings use `~/.local/share/ontseq/resources-v0.7.1` for resources. The release
   suffix follows the Core version; previously saved resource paths remain selected.
2. Select the reference family in **System einrichten**, then install or validate that
   family. GRCh37.p13 and UCSC hg19 use `GRCh37_GENCODE19_HG19_v2`,
   `HEMATOLOGY_GRCh37_v1` and the
   strictly GRCh37-bound `AML_AS_111_GRCh37_v1` panel.
3. Choose `AML_LCWGS_GRCh37` or `AML_AS_111_GRCh37` only for BAMs aligned against the full GENCODE 19 publisher
   FASTA (chromosomes, scaffolds, patches and haplotypes). Its `exact_full` contract is not
   interchangeable with an arbitrary hg19 or primary-only BAM. Choose
   `AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25` or
   `AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25` only for the exact ordered UCSC hg19 dictionary
   `chr1`–`chr22`, `chrX`, `chrY`, `chrM=16571`, with no extra contigs. No reheadering,
   fallback or runtime liftover is performed. The lcWGS profiles remain panel-free; the two Adaptive
   Sampling profiles resolve only `AML_AS_111_GRCh37_v1` (110/111 selection intervals mapped,
   109 reciprocal-exact and 108 native GENCODE-19 ROIs; ACACA requires roundtrip review, while
   CT45A2, IGH and GPR128 remain mapping/ROI gaps).
4. The four GRCh38 choices keep their existing exact-full/Canonical-25 contracts. A
   missing GRCh37 installation must not block a ready GRCh38 profile, or vice versa.
5. Select a synthetic BAM and its BAI for the first test. Start through Desktop or use
   **Live-Workspace** for the connected browser UI. Both use the same pipeline and files.
6. Inspect the actual stage outcomes, QC, reference provenance and HTML/XLSX/JSON outputs.
   A NO_CALL or missing stage is not a negative result. Existing browser-started analyses
   are not silently stopped when returning to Desktop.
7. Where eligible CNV evidence exists, HTML/XLSX/JSON expose a partial ISCN event-fragment
   proposal with source-event and resource traceability. It is not a complete karyotype,
   requires cytogenetic expert review and cannot be released automatically.

Multiple engineering versions may run concurrently. When the configured loopback port is
already occupied, a new Desktop instance automatically starts its own backend on a free local
port and shows the effective port. It never connects to, takes over or stops the older service.
Every spawned backend must echo a fresh per-launch instance ID and the exact resource, output and
allowed-input roots before Desktop commits the connection. A simultaneous bind race is rejected
and retried on a fresh port; a failed candidate never becomes reusable Desktop state.
The profile selected in Desktop is handed to the Live Workspace and remains selected only when
that exact identifier is advertised by the connected backend.
Desktop intentionally keeps the eight supported profile contracts in a typed allowlist; backend
status is read dynamically to enable only complete local build families. An arbitrary profile ID
reported by another backend version is not inferred or exposed without a matching Desktop contract.

### Separately installed methylation tool

An independently installed local modkit executable can be selected with the optional
`modkitExecutableWsl` setting, for example `/opt/local-tools/modkit/bin/modkit`.
Use an absolute Linux file path inside the selected WSL distribution; spaces within the
path are preserved. Relative paths, shell expressions and Windows paths are refused.
The default `null` keeps the backend's existing `modkit` lookup through its runtime PATH.
Desktop passes the configured path as a separate `--modkit` process argument for both
profile-based and legacy service starts; it is not added to `system-smoke`.
Restart Desktop after changing this setting so that the running service retains its
original tool configuration. This setting neither installs modkit nor changes any runtime.
The backend checks the configured executable against the versioned methylation policy;
missing or mismatched tools keep methylation unavailable while genome-only analysis remains
available. The operator must still explicitly choose whether to include methylation.

GRCh37 supports lcWGS and Adaptive Sampling through four explicit profile identifiers. The two
dictionary contracts are intentionally separate, and the new panel is materialized and pinned for
GRCh37 at build time; it is never lifted or substituted during a run. The bundle exposes its
110/111 selection, 109 reciprocal-exact and 108-ROI cardinalities instead of presenting an
unresolved source target as covered. There is no GRCh37 MANE dataset. The new knowledge bundle is coordinate-free and
retains its original curated source-selection limits. Source-tree/checksum metadata in
the local package identifies its contents; an uncommitted build must not claim the base
archive's historical git commit. See [native integration](../docs/GRCH37_INTEGRATION.md).

## Historical v0.5.3 baseline (superseded operator instructions)

The remaining sections document the inherited implementation and earlier verification
claims. Use the v0.7.1 instructions above for new installation/build selection; historical
CI results are not evidence that this local candidate has passed a new release workflow.

## v0.5.3 user path

1. Extract the complete engineering bundle and start `ONTSeq.Desktop.exe`.
2. Open **System einrichten** on first use and install the bundled Linux runtime into WSL. The bundle contains Python, R, samtools, Cramino, Sniffles2, cuteSV 2.1.3, Mosdepth, minimap2, QDNAseq, QDNAseq.hg19, QDNAseq.hg38, ACE and the ONTSeq Core/resources.
3. Run **Selbsttest starten**. This is a full installed-runtime system smoke, not only a package import check.
4. Keep the user-writable WSL resource root at `~/.local/share/ontseq/resources` or configure another absolute WSL path. Use **Installieren** or **Reparieren** to provision and validate the complete resource family: `GRCh38_GENCODE50_MANE1.5_v1`, `HEMATOLOGY_v3`, `AML_AS_111_GRCh38_v1` and all four profiles. Repair safely replaces divergent pinned resources through the backend transaction; manual deletion is not required.
5. Select one of the four explicit profiles. `AML_LCWGS_GRCh38` (default) and `AML_AS_111_GRCh38` require the complete Primary-Assembly dictionary (`exact_full`). `AML_LCWGS_GRCh38_CANONICAL25` and `AML_AS_111_GRCh38_CANONICAL25` require exactly `chr1`-`chr22`, `chrX`, `chrY`, `chrM`. All four resolve the same pinned reference, annotation and knowledge bundles automatically; the AS profiles additionally resolve `AML_AS_111_GRCh38_v1`.
6. Select one aligned `.bam` file. Both `sample.bam.bai` and `sample.bai` are accepted, with deterministic precedence for `sample.bam.bai` when both exist. The declared index is checked explicitly rather than rediscovered implicitly by samtools.
7. Sample ID is suggested from the filename and remains editable.
8. Select the GRCh38 profile whose displayed dictionary contract matches the reference used for alignment. The backend verifies the complete BAM sequence dictionary and the desktop replaces the expected-build display with the build persisted in run provenance when available.
9. Click **ANALYSE STARTEN**.
10. The application checks WSL, backend capabilities, input/index, resource root, required profile bundles and output location before starting the local service.
11. It starts the canonical `ontseq serve` backend, submits the run, polls progress, and enables the HTML, XLSX and result-folder buttons when outputs are present.

v0.5.3 retains the working-directory-independent packaged policies from v0.5.1 and adds the two
explicit Canonical-25 choices. Canonical-25 is resolved from the already installed full GRCh38
`ReferenceLock`; it does not install or download a second multi-GiB reference. Neither contract
uses liftover or falls back to the other. Scientific thresholds and interpretation policies are
unchanged.

The UI keeps the backend's distinct stage outcomes (`COMPLETED`, `NO_CALL`, `FAILED`, `NOT_RUN`) instead of flattening them into done/pending. A technical `PASS` is not a biological negative result and is not clinical validation.

### Reference and resume integrity in v0.5.3

New profile runs resolve a manifest-driven `ReferenceBundle`, `PanelBundle` and `KnowledgeBundle` beneath `resourceRootWsl`. Bundle installation and repair remain backend operations; Desktop only invokes their typed CLI contracts and displays the result. The selected profile and BAM dictionary must both resolve to GRCh38 and must satisfy either the selected `exact_full` or exact Canonical-25 contract.

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
compatibility, but all four published Desktop profiles in this work package request only the
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

See `desktop.settings.example.json` for the available keys. First-run setup writes the user-local settings automatically. A runtime installed by the v0.5.3 Desktop is stored under `~/.local/share/ontseq/runtime-v0.5.3` in the selected WSL distribution.

The current profile-resource defaults for fresh settings are:

```json
{
  "resourceRootWsl": "~/.local/share/ontseq/resources-v0.7.1",
  "defaultProfile": "AML_LCWGS_GRCh38"
}
```

The directory suffix identifies the Desktop/Core release. Reference, panel and knowledge
bundle IDs and checksums retain their independent versions. Existing settings, including
custom names such as `resources-v0.6.2-hg19-as`, are preserved. Changing a setting does not
move or validate its contents; see the [resource-root transition](../docs/REFERENCE_SYSTEM.md#resource-root)
before selecting a renamed installation. Core/CLI retains its `/opt/ontseq` default.

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

## Deliberate v0.5.3 limitation: cancellation

The **Abbrechen** control remains disabled. A backend cancellation contract must first guarantee that an interrupted external tool cannot leave a stage looking complete. Closing a progress view is not treated as cancellation.

## Remaining workstation / validation gates

- first real Windows + WSL BAM run on the target workstation, including HTML/XLSX opening from the desktop;
- mapped-network-drive verification on the actual workstation if the BAM resides on `P:` or another mapped drive;
- service-level Windows smoke covering Desktop Start → loopback API → persisted outputs with a real `wsl.exe` host rather than Linux-container emulation;
- analytical validation of GRCh38 CNV on controlled real/reference samples for the active profiles;
- backend-owned cooperative cancellation with explicit interrupted-state semantics;
- analytical validation on controlled real/reference samples before any diagnostic use.

See also `README-FIRST-RUN.md` and `ONTSeq.Desktop/CHANGELOG.md`.
