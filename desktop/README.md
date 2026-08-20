# ONTSeq Desktop

Windows operator surface for the ONTSeq execution core. The desktop app does not reimplement bioinformatics; it launches the bundled Linux runtime through WSL2, talks to the loopback service, and opens the same persisted HTML/XLSX/result artifacts produced by the canonical backend.

## v0.2 user path

1. Start `ONTSeq.Desktop.exe` and open **System einrichten** on first use.
2. Install the bundled Linux runtime into WSL. The bundle contains Python, R, samtools, Cramino, Sniffles2, Mosdepth, minimap2, QDNAseq, QDNAseq.hg19, ACE and the ONTSeq package/resources.
3. Configure the Reference-Lock that matches the BAM (GRCh37 or GRCh38).
4. For Adaptive Sampling, select the actual analysis ROI BED once. The desktop copies it into its local resource store and records a SHA256-backed version string.
5. Select one aligned `.bam` file. Both `sample.bam.bai` and `sample.bai` are accepted, with deterministic precedence for `sample.bam.bai` when both exist.
6. Sample ID is suggested from the filename and remains editable.
7. Select GRCh38 or GRCh37 and low-coverage WGS or Adaptive Sampling.
8. Click **ANALYSE STARTEN**.
9. The application checks WSL, the input/index, configured reference, required runtime resources, output location and Adaptive Sampling BED when applicable.
10. It starts the canonical `ontseq serve` backend, submits the run, polls progress, and enables the HTML, XLSX and result-folder buttons when outputs are present.

The UI keeps the backend's distinct stage outcomes (`COMPLETED`, `NO_CALL`, `FAILED`, `NOT_RUN`) instead of flattening them into done/pending.

### Current CNV behavior

The bundled reproducible QDNAseq environment currently includes the Bioconductor `QDNAseq.hg19` annotation package. Therefore a GRCh37 desktop run requests the integrated QDNAseq + ACE CNV stage (100/500/1000 kbp, 500-kbp primary view and multi-bin consensus). GRCh38 runs continue through the other requested modules but do not request QDNAseq CNV until the hg38 annotation resource has its own pinned, real-tool-tested packaging path.

## Why .NET 10 + WPF

The desktop is a small Windows-native shell around the Linux bioinformatics backend. Filesystem integration, process control and a single self-contained Windows executable are first-class requirements; the analytical implementation stays independently testable in Python/R/Linux.

## Configuration

The app loads the first existing settings file from:

- `%LOCALAPPDATA%\ONTSeq\desktop.settings.json`
- `%PROGRAMDATA%\ONTSeq\desktop.settings.json`

See `desktop.settings.example.json` for the available keys. First-run setup writes the user-local settings automatically.

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

`desktop-ci.yml` additionally builds a relocatable Linux runtime from the pinned QDNAseq environment, proves it after relocation and inside stock Ubuntu without host Python/R, then embeds that archive beside the Windows executable.

## Deliberate v0.2 limitation: cancellation

The **Abbrechen** control remains disabled. A backend cancellation contract must first guarantee that an interrupted external tool cannot leave a stage looking complete. Closing a progress view is not treated as cancellation.

## Next desktop gates

- complete CI bootstrap proof on a clean Linux userspace and Windows WPF publish;
- HTTP/service contract smoke covering Start → Status → persisted outputs;
- first real Windows + WSL BAM run on the target workstation;
- robust mapped-network-drive setup for workstations where `/mnt/<drive>` is not already available;
- GRCh38 QDNAseq annotation packaging and real-tool test;
- backend-owned cooperative cancellation with explicit interrupted-state semantics.
