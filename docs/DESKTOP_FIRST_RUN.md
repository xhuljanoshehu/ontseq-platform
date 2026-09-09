# ONTSeq Desktop v0.7.1 — first run (engineering / RUO)

ONTSeq Desktop is deliberately fail-closed. A BAM is not analysed until the local Linux runtime
and the pinned GRCh37 or GRCh38 bundles required by its selected profile are available.

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
6. Confirm `resourceRootWsl` (fresh-settings default `~/.local/share/ontseq/resources-v0.7.1`, below the WSL user's
   `$HOME`) and choose **Installieren** for either `GRCh38_GENCODE50_MANE1.5_v1` or
   `GRCh37_GENCODE19_HG19_v2`, matching the profiles you intend to run. The GRCh37 bundle
   contains separate native full-build and exact UCSC hg19 Canonical-25 contracts. Interactive status checks
   manifests, pins, presence and declared sizes; `ontseq references validate` performs the
   explicit full SHA256 audit. Use **Reparieren** for damaged resources.
   Saved resource paths are preserved. The root suffix follows the software release;
   reference, panel and knowledge bundle versions remain independent. A changed setting
   does not move files; follow the [resource-root transition](REFERENCE_SYSTEM.md#resource-root)
   before selecting a deliberately renamed installation.
7. Run **Selbsttest starten** before any research BAM is used.
8. Return to the main window, select BAM + `.bam.bai`, choose the profile whose dictionary label
   matches the alignment reference, and start the analysis. Unsuffixed profiles require their
   complete pinned dictionary. GRCh38 `*_CANONICAL25` profiles require exactly `chr1`-`chr22`,
   `chrX`, `chrY`, `chrM`; both GRCh37 `*_UCSC_HG19_CANONICAL25` profiles require the same
   ordered names with UCSC-hg19 lengths and `chrM=16571`.

## Reference safety

This work package publishes four GRCh38 profiles and four GRCh37 profiles. The registry
activates manifest-bearing bundles, each profile pins exact build-matched IDs, and aligned-BAM
intake compares the complete BAM dictionary with the selected contract before analysis.
`AML_LCWGS_GRCh37` and `AML_AS_111_GRCh37` remain `exact_full`; their two
`*_UCSC_HG19_CANONICAL25` counterparts accept only the ordered UCSC-hg19 25-contig dictionary.
GRCh38 keeps its existing exact-full and Canonical-25 contracts. There is no run-time liftover,
cross-build or dictionary-contract fallback.

Existing settings containing explicit GRCh37/GRCh38 Reference-Locks or Adaptive-Sampling BEDs
remain readable for one compatibility release. They are not combined with a new profile-resolved
resource context.

## Adaptive Sampling

The two GRCh38 profiles resolve only `AML_AS_111_GRCh38_v1`; the two GRCh37 profiles resolve only
`AML_AS_111_GRCh37_v1`. The GRCh37 bundle is a checksum-pinned build-time projection of the
laboratory selection design: 110 of 111 source intervals map at the locked threshold, 109 are
reciprocal-exact and 108 targets obtain native GENCODE-19 analysis ROIs. `ACACA` requires
roundtrip review; `CT45A2`, `IGH` and `GPR128` remain explicit mapping/ROI gaps. The operator does
not browse for panel files, and unresolved targets are never
turned into negative observability statements. Coordinate projection is not executed during an
analysis and the projection tool/chain are not shipped as runtime dependencies.

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
The 0.7.1 report may contain a partial, CNV-derived ISCN event-fragment proposal. It requires
expert cytogenetic review, is not a complete karyotype and is never eligible for automatic
clinical release.
