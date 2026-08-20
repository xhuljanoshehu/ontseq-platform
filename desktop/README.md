# ONTSeq Desktop

Windows operator surface for the existing ONTSeq pipeline. The desktop app does **not** reimplement QC, CNV, SV, fusion, annotation or reporting logic. It starts the Linux backend through WSL2, talks only to the loopback service, and reads the same persisted `provenance/run.json` that underpins the release bundle.

## v0.1 user path

1. Select one aligned `.bam` file.
2. Sample ID is suggested from the filename and remains editable.
3. Select GRCh38 or GRCh37.
4. Select low-coverage WGS or Adaptive Sampling.
5. Click **ANALYSE STARTEN**.
6. The application checks WSL, the allowed input root, the configured reference lock and the Windows output directory.
7. It starts `ontseq serve` inside the configured WSL distribution, submits the run, and polls both the API and the atomically persisted run report.
8. When outputs exist, the user can open HTML, XLSX or the result folder directly.

The UI keeps the backend's distinct stage outcomes (`COMPLETED`, `NO_CALL`, `FAILED`, `NOT_RUN`) rather than flattening them into done/pending.

## Why .NET 10 + WPF

This branch intentionally uses a Windows-only native desktop shell. The deployment target is a Windows clinical/research workstation, local filesystem integration and process control are first-class requirements, and the bioinformatics already lives behind a local service boundary. WPF keeps the desktop layer small and avoids embedding a second browser runtime as the primary application architecture. The backend remains Python/Linux and therefore stays independently testable and deployable.

## Configuration

The installer, not the operator, should own the technical configuration. The app loads the first existing file from:

- `%LOCALAPPDATA%\ONTSeq\desktop.settings.json`
- `%PROGRAMDATA%\ONTSeq\desktop.settings.json`

See `desktop.settings.example.json`.

A production installer must replace the example paths with locally approved, checksummed reference and assay resources. Do **not** ship the literal example as a validated configuration.

### Network drives

The current backend translates a Windows drive such as `P:\Lab\run.bam` to `/mnt/p/Lab/run.bam`. That is only usable if the mapped drive is actually visible to WSL. The desktop preflight checks this before starting an analysis. UNC paths are intentionally refused until the share has a deliberate WSL mount.

## Build

```powershell
dotnet build desktop/ONTSeq.Desktop/ONTSeq.Desktop.csproj -c Release
```

Single-file self-contained publish:

```powershell
dotnet publish desktop/ONTSeq.Desktop/ONTSeq.Desktop.csproj -c Release -r win-x64 --self-contained true -p:PublishSingleFile=true
```

The resulting executable is an engineering artifact, not a clinical release.

## Deliberate v0.1 limitation: cancellation

The **Abbrechen** control is visible but disabled. Killing an arbitrary bioinformatics subprocess from the GUI without a backend cancellation contract could leave a stage half-written or misstate its outcome. The runner already uses atomic writes, but cancellation must be added at the orchestration/command-runner boundary and tested under failure injection before the control is enabled. The UI must never pretend that closing a progress view cancelled the analysis.

## Next desktop gates

- backend-owned cooperative cancellation with `CANCEL_REQUESTED/CANCELLED` job semantics while preserving stage `NOT_RUN` reasons;
- installer/bootstrap for WSL2 plus the pinned offline backend image/environment;
- signed Windows package and upgrade/rollback path;
- Windows CI smoke that launches a synthetic backend fixture and exercises Start → Status → Results;
- accessibility and keyboard-navigation review;
- in-app review of CNV/SV/fusion evidence after those result contracts are analytically locked.
