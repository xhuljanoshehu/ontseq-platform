# ONTSeq Desktop v0.1

ONTSeq Desktop is the clinician-/reviewer-facing Windows shell around the existing ONTSeq
research pipeline. The design goal is that routine users do not need to operate Python,
Snakemake, Conda or a Linux shell.

> **Research use only. Not clinically validated.** The desktop application does not convert
> candidate evidence into a clinical diagnosis or automatically release a report.

## User workflow

The default screen is intentionally small:

1. Select an **aligned BAM**.
2. The application proposes a pseudonymized sample ID from the filename.
3. Select **GRCh38/GRCh37**.
4. Select **Adaptive Sampling** or **low-coverage WGS**.
5. Press **ANALYSE STARTEN**.
6. Follow explicit module states.
7. Open the generated **HTML**, **Excel** or result folder.

The interface explicitly distinguishes:

- `PASS`
- `WARN`
- `FAIL`
- `NO_CALL`
- `NOT_RUN`

`NO_CALL` is not presented as a biological negative result. `NOT_RUN` is not presented as a
negative result. The current aligned-BAM MVP runs input/reference gates, Cramino QC,
Sniffles2 candidate SV evidence, result assembly and HTML/XLSX/JSON rendering. CNV, fusion
interpretation and ISCN remain visibly `NOT_RUN` until their benchmark/validation gates are
accepted and they are wired into the aligned-BAM production path.

## Architecture

```text
ONTSeqDesktop.exe
        |
        | local-only job/configuration layer
        v
Windows WSL2 boundary
        |
        v
ONTSeq CLI in Linux
        |
        +-- samtools input/reference gate
        +-- Cramino QC
        +-- Sniffles2 candidate SV evidence
        +-- typed ONTSeq result contract
        |
        v
HTML + XLSX + JSON + local run log
```

The Windows executable contains the graphical user interface. Bioinformatics tools remain in
the controlled WSL/Linux environment. This avoids attempting to repackage Linux-native
bioinformatics executables into the GUI binary and keeps the analytical environment
reproducible.

## First-time deployment

### 1. Prepare the WSL backend

The WSL distribution must contain the ONTSeq repository/environment and the aligned-BAM tool
stack. The repository's existing local smoke test remains the preferred technical verification
of samtools, Cramino, Sniffles2 and cuteSV wiring.

The desktop v0.1 runner currently invokes these CLI boundaries:

```text
ontseq inspect-bam
ontseq qc-cramino
ontseq call-sniffles
ontseq assemble-aligned-mvp
ontseq render
```

### 2. Configure Desktop settings

Open **Einstellungen** and configure:

- WSL distribution (empty means the Windows default distribution)
- ONTSeq project root in WSL, e.g. `~/ontseq-platform`
- local Windows result directory
- QC policy path relative to the ONTSeq project
- Sniffles2 policy path relative to the ONTSeq project
- one deployment-local reference profile per supported build

Each reference profile contains:

- genome build
- `reference_id`
- path to the locked reference JSON
- optional Adaptive-Sampling target BED path and explicit BED version

The application verifies that the selected reference lock's build and `reference_id` match the
desktop profile before starting analysis.

### 3. Run Systemcheck

The GUI **Systemcheck** verifies the output folder, WSL availability, required command-line
tools and configured local reference resources. A failed check should be corrected before
patient-derived research material is processed.

## Data boundary

The desktop configuration stores tool and reference paths only. Runtime BAM/BAI, VCF and
result files stay local/on-premises. The generated sample manifest always sets:

```yaml
privacy:
  pseudonymized: true
  contains_direct_identifiers: false
  cloud_upload_approved: false
```

Use pseudonymized sample IDs only. Do not place patient BAM/VCF files, clinical reports or
direct identifiers in GitHub.

## Build the Windows executable

From PowerShell on Windows:

```powershell
./scripts/build_windows_desktop.ps1
```

The build installs the optional `desktop` dependencies, executes repository tests/linting and
creates:

```text
dist/ONTSeqDesktop.exe
```

GitHub Actions also builds this executable on the desktop pull request and uploads it as a
workflow artifact. The executable is only the user-interface shell; a configured WSL backend
is still required for real analysis.

## Current v0.1 scope

Implemented:

- Windows-style GUI
- BAM selection
- automatic `.bam.bai` / `.bai` discovery
- pseudonymized sample ID handling
- GRCh37/GRCh38 selection
- lcWGS/Adaptive-Sampling selection
- local deployment configuration
- WSL2 backend abstraction
- explicit progress states
- local technical log
- start/cancel controls
- HTML/Excel/result-folder buttons
- reference-lock checks
- Windows `.exe` packaging workflow

Not yet claimed as complete:

- validated CNV calling in the desktop aligned-BAM path
- clinically validated fusion interpretation
- clinically conformant automated ISCN output
- user authentication/signature
- immutable clinical release bundle
- LIMS integration
- diagnostic/medical-device release

Those remain separate engineering and analytical validation gates rather than being hidden
behind a polished GUI.
