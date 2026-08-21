# ONTSeq Desktop

Windows operator surface for the ONTSeq execution core. The desktop app does not reimplement bioinformatics; it launches the bundled Linux runtime through WSL2, talks to the loopback service, and opens the same persisted HTML/XLSX/result artifacts produced by the canonical backend.

## v0.2.2 user path

1. Start `ONTSeq.Desktop.exe` and open **System einrichten** on first use.
2. Install the bundled Linux runtime into WSL. The bundle contains Python, R, samtools, Cramino, Sniffles2, Mosdepth, minimap2, QDNAseq, QDNAseq.hg19, ACE and the ONTSeq package/resources.
3. Run **Selbsttest starten**. This is now a full installed-runtime system smoke, not only a package import check.
4. Configure the Reference-Lock that matches the BAM (GRCh37 or GRCh38).
5. For Adaptive Sampling, select the actual analysis ROI BED once. The desktop copies it into its local resource store and records a SHA256-backed version string.
6. Select one aligned `.bam` file. Both `sample.bam.bai` and `sample.bai` are accepted, with deterministic precedence for `sample.bam.bai` when both exist.
7. Sample ID is suggested from the filename and remains editable.
8. Select GRCh38 or GRCh37 and low-coverage WGS or Adaptive Sampling.
9. Click **ANALYSE STARTEN**.
10. The application checks WSL, the input/index, configured reference, required runtime resources, output location and Adaptive Sampling BED when applicable.
11. It starts the canonical `ontseq serve` backend, submits the run, polls progress, and enables the HTML, XLSX and result-folder buttons when outputs are present.

The UI keeps the backend's distinct stage outcomes (`COMPLETED`, `NO_CALL`, `FAILED`, `NOT_RUN`) instead of flattening them into done/pending.

The Desktop accepts a file as named GRCh37/GRCh38 only when all canonical nuclear chromosomes
(`1`-`22`, `X`, `Y`) occur with the assembly's standard lengths in one consistent naming style.
Optional mitochondrial, decoy, ALT and random contigs remain permitted, but a partial,
mixed-style or wrong-build FAI dictionary cannot be presented as a complete named reference.
Setup shows build, contig count, total reference bases and the FAI-dictionary-addressed reference
ID for every valid lock. This establishes sequence-dictionary compatibility; an FAI does not
cryptographically prove the underlying FASTA bases are identical.

Region-extracted BAMs remain compatible when their header retains the complete dictionary of the
original alignment reference. A BAM with a full hs37d5 header must use the full hs37d5 lock; a
bundled chromosome-only analysis FASTA is not the reference that originally produced that BAM.

### Full installed-runtime self-test

`ontseq system-smoke` deliberately uses two independent deterministic synthetic fixtures so that one artificial signal does not distort the other analytical lane:

- a long-read fixture exercises BAM intake, samtools, Cramino, Sniffles2 and basic report rendering;
- a whole-genome-shaped GRCh37/hg19 read-depth fixture exercises the canonical run envelope with live QDNAseq + ACE. Baseline autosomes have two synthetic read starts per 100 kbp, chromosome 7 has half-depth and chromosome 8 has 1.5x depth.

A PASS requires the CNV path to recover chromosome 7 at approximately copy number 1 and chromosome 8 at approximately copy number 3, with agreement across the configured 100/500/1000-kbp QDNAseq resolutions. It also requires the CNV result to be present in normalized JSON, provenance, HTML and the `CNV Fits`, `CNV Consensus` and `CNV Segments` Excel sheets. Release SHA256 checksums are independently recalculated, the identical canonical run is executed a second time, the CNV stage must resume content-addressed, and release checksums must still verify afterward.

The desktop setup invokes exactly that installed-runtime command through WSL and stores `system-smoke.report.json`, the canonical synthetic run envelope and a text log under `%LOCALAPPDATA%\ONTSeq\self-test\<timestamp>`.

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

`desktop-ci.yml` builds a relocatable Linux runtime from the pinned QDNAseq environment, proves the packed environment after relocation, then boots it inside stock Ubuntu 24.04 without host Python or R and executes the full `system-smoke`. Only after that passes is the Windows WPF executable published and bundled with the exact tested Linux runtime.

## Deliberate v0.2.2 limitation: cancellation

The **Abbrechen** control remains disabled. A backend cancellation contract must first guarantee that an interrupted external tool cannot leave a stage looking complete. Closing a progress view is not treated as cancellation.

## Next desktop gates

- first real Windows + WSL BAM run on the target workstation, including HTML/XLSX opening from the desktop;
- mapped-network-drive verification on the actual workstation if the BAM resides on `P:` or another mapped drive;
- service-level Windows smoke covering Desktop Start → loopback API → persisted outputs with a real `wsl.exe` host rather than Linux-container emulation;
- GRCh38 QDNAseq annotation packaging and real-tool test;
- backend-owned cooperative cancellation with explicit interrupted-state semantics.
