# ONTSeq Desktop v0.1.1 — first run (engineering / RUO)

ONTSeq Desktop is deliberately fail-closed. A BAM is not analysed until the local Linux runtime and the pinned GRCh38 bundles required by its selected profile are available.

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
6. Confirm `resourceRootWsl` (default `~/.local/share/ontseq/resources`, below the WSL user's `$HOME`) and choose **Installieren** for `GRCh38_GENCODE50_MANE1.5_v1`. Interactive status checks manifests, pins, presence and declared sizes; `ontseq references validate` performs the explicit full SHA256 audit. Use **Reparieren** for damaged resources.
7. Run **Selbsttest starten** before any research BAM is used.
8. Return to the main window, select BAM + `.bam.bai`, choose `AML_LCWGS_GRCh38` or `AML_AS_111_GRCh38`, and start the analysis.

## Reference safety

This work package publishes only GRCh38 profiles. The registry activates manifest-bearing bundles,
the profile pins their exact IDs, and aligned-BAM intake compares the complete BAM dictionary with
the installed GRCh38 Reference-Lock before analysis. GRCh37, partial and mixed dictionaries fail
before the pipeline starts; there is no liftover or cross-build fallback.

Existing settings containing explicit GRCh37/GRCh38 Reference-Locks or Adaptive-Sampling BEDs
remain readable for one compatibility release. They are not combined with a new profile-resolved
resource context.

## Adaptive Sampling

`AML_AS_111_GRCh38` resolves the controlled selection panel and downstream analysis ROI as
separate artifacts from `AML_AS_111_GRCh38_v1`. The operator does not browse for either file.
Unresolved targets stay unresolved and are not turned into negative observability statements.

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
