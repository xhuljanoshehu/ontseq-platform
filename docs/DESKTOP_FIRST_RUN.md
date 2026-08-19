# ONTSeq Desktop v0.1.3 — first run (engineering / RUO)

ONTSeq Desktop is deliberately fail-closed. A BAM is not analysed until the local Linux runtime and an explicit reference lock are available.

## What the bundle contains

The Windows engineering ZIP contains:

- `ONTSeq.Desktop.exe` — native Windows operator UI.
- `runtime/ontseq-linux-runtime.tar.gz` — relocatable, pinned Linux analysis runtime containing the ONTSeq Python package plus the aligned-BAM toolchain used by this build.

The runtime is installed into the current WSL user's home directory. The operator does not need to edit `desktop.settings.json` by hand.

## Prerequisite

WSL2 and the configured Linux distribution (`Ubuntu` by default) must already be installed and able to start. WSL installation can require administrator rights and a Windows restart, so this engineering build checks it but does not silently modify Windows features.

## First-run workflow

1. Extract the **whole ZIP**. Do not copy only `ONTSeq.Desktop.exe`; the `runtime` folder is required.
2. Start `ONTSeq.Desktop.exe`.
3. Choose **System einrichten**.
4. Choose **System prüfen**.
5. If the backend is missing, choose **Runtime installieren**.
   When upgrading from v0.1.2, install the bundled runtime again so the desktop uses
   `runtime-v0.1.3` with the corrected Windows-to-backend request transport.
6. Configure the reference build needed for the BAM by selecting the **exact FASTA or FAI used for alignment**.
7. Run **Selbsttest starten**. This creates a synthetic alignment and runs the real bundled samtools/Cramino/Sniffles2 path before any research BAM is used.
8. Return to the main window, select BAM + `.bam.bai`, choose the matching build/profile, and start the analysis.

## Reference safety

`GRCh38` and `GRCh37` are assembly labels, not sufficient reference identity. Different FASTA distributions can differ in contig names, alternative/decoy contigs and sequence content. ONTSeq therefore does **not** download an arbitrary reference when a build is selected.

The setup wizard creates a versioned reference lock from the chosen FASTA index. The aligned-BAM intake subsequently compares the BAM header against that lock and fails closed on incompatible reference structure.

For a public benchmark such as HG002, use the exact reference distribution documented for the downloaded alignment. For local laboratory data, use the locally authorised reference bundle used to create the BAM.

## Adaptive Sampling

Adaptive Sampling additionally requires the approved analysis ROI BED and its version. The MinKNOW selection BED and the downstream analysis ROI are not assumed to be interchangeable. Until an authorised ROI is configured, an Adaptive Sampling run is stopped rather than silently falling back to WGS semantics.

## Local files written by the desktop app

Per-user configuration:

`%LOCALAPPDATA%\ONTSeq\desktop.settings.json`

Reference locks:

`%LOCALAPPDATA%\ONTSeq\references\`

Synthetic self-test results:

`%LOCALAPPDATA%\ONTSeq\self-test\`

Normal result envelopes default to:

`%USERPROFILE%\Documents\ONTSeq\results\`

## Scope

Research Use Only. This setup workflow proves deployability and technical execution. It does not constitute analytical or clinical validation, and `NO_CALL`, `FAILED` and `NOT_RUN` must not be interpreted as negative biological findings.
