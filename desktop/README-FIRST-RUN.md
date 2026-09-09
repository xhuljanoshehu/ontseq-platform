# ONTSeq Desktop v0.7.1 — first workstation run

This checklist is for the unsigned engineering/research build only. It is not analytical or clinical validation.

1. Confirm Windows WSL2 is enabled and the configured distribution opens normally.
2. Extract the complete `ontseq-desktop-v0.7.1-win-x64-setup-engineering` bundle; do not separate `ONTSeq.Desktop.exe` from the bundled `runtime` directory, Core wheel or SHA256SUMS.
3. Start `ONTSeq.Desktop.exe`, open **System einrichten**, and install the bundled ONTSeq Linux runtime into its new, separate prefix. The installed runtime must report Core v0.7.1 and expose the reference-validation, target-coverage, component-selection and complete SV-policy capabilities checked by the Desktop preflight.
   The preflight also verifies pinned cuteSV 2.1.3 and absolute packaged policy paths; this
   runtime correction does not change scientific thresholds.
4. Run **Selbsttest starten** before using any real sample. The self-test must finish successfully and produce its synthetic report bundle.
5. Fresh settings use `~/.local/share/ontseq/resources-v0.7.1`; previously saved paths remain selected. For a separate test, choose a new resource root or explicitly select an existing validated installation. Select and install either `GRCh38_GENCODE50_MANE1.5_v1` or `GRCh37_GENCODE19_HG19_v2` when that family is absent. The latter provides both the native full GRCh37.p13 contract and the exact UCSC hg19 Canonical-25 contract required for adaptive sampling. Fast bundle status must pass; use **Reparieren** for damaged resources. Run `ontseq references validate` for a full SHA256 audit. The root's software-version suffix does not change scientific bundle versions or move existing resources.
6. Choose the profile matching the BAM dictionary. `AML_LCWGS_GRCh38` and
   `AML_AS_111_GRCh38` remain `exact_full` and require the complete ordered Primary-Assembly
   lock. `AML_LCWGS_GRCh38_CANONICAL25` and `AML_AS_111_GRCh38_CANONICAL25` require exactly
   `chr1`-`chr22`, `chrX`, `chrY`, `chrM`. The AS profiles resolve their controlled panel and
   analysis ROI automatically. All four use the same installed bundles, so switching to
   Canonical-25 does not repeat the multi-GiB reference download. The new
   `AML_LCWGS_GRCh37` and `AML_AS_111_GRCh37` instead require the full native GENCODE 19
   GRCh37.p13 dictionary. `AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25` and
   `AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25` separately require exactly `chr1`–`chr22`,
   `chrX`, `chrY`, `chrM=16571` in UCSC hg19 order and no extras. These contracts are not
   aliases. The two GRCh37 Adaptive-Sampling profiles resolve only the build-time materialized,
   checksum-pinned `AML_AS_111_GRCh37_v1` panel; there is no reheadering, build fallback or
   runtime liftover. The UI states the actual contract: 110/111 selection intervals mapped,
   109 reciprocal-exact and 108 native GENCODE-19 ROIs; ACACA requires roundtrip review, while
   CT45A2, IGH and GPR128 remain mapping/ROI gaps.
7. Confirm the selected Windows BAM root is visible in WSL. For example, if input is on `P:`, `/mnt/p` must exist and point to the intended mapped drive. UNC paths are not guessed.
8. Select a **synthetic/non-patient** aligned BAM first together with its declared BAI (`sample.bam.bai` preferred; `sample.bai` also supported) and execute the complete Desktop path.
9. Open **Live-Workspace** to inspect the real run or use Desktop output buttons. Verify that HTML/XLSX/JSON and `provenance/run.json` are produced, build/profile binding agrees, and failed/`NOT_RUN`/`NO_CALL` stages are not interpreted as negative biological findings.
10. Inspect the ISCN section in the synthetic HTML/XLSX result. Any notation is a partial,
    unvalidated CNV-derived proposal requiring expert cytogenetic review; it is not a complete
    karyotype and cannot be released automatically.
11. Only after this workstation smoke succeeds should locally authorized research data be used.

If another ONTSeq engineering build already owns the configured loopback port, this Desktop
starts its backend on a free local port and displays that effective port. It does not adopt or
terminate the other process. The selected Desktop profile is passed to the Live Workspace and
is retained only when the connected service advertises the exact identifier in `/api/config`.
Desktop additionally verifies a fresh per-launch instance ID and the exact resource, output and
allowed-input roots before it accepts the service. A simultaneous bind race is retried on a new
port; a failed startup is not retained for the next attempt.

Current deliberate limitations: the build is unsigned and RUO. GRCh37.p13 and GRCh38 profiles
are selectable but scientifically distinct. Native GRCh37 has no MANE; its Adaptive-Sampling
profiles use a separate GRCh37-bound panel and never the GRCh38 panel. Legacy explicit path fields
remain readable for one compatibility release but are not mixed into profile runs. Cancellation
remains disabled until the backend can guarantee an explicit interruption state.
