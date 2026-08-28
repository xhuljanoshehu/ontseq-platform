# ONTSeq Desktop v0.5.1 — first workstation run

This checklist is for the unsigned engineering/research build only. It is not analytical or clinical validation.

1. Confirm Windows WSL2 is enabled and the configured distribution opens normally.
2. Extract the complete `ontseq-desktop-v0.5.1-win-x64-setup-engineering` bundle; do not separate `ONTSeq.Desktop.exe` from the bundled `runtime` directory.
3. Start `ONTSeq.Desktop.exe`, open **System einrichten**, and install the bundled ONTSeq Linux runtime. The installed runtime must report Core v0.5.1 and expose the reference-validation, target-coverage, component-selection and complete SV-policy capabilities checked by the Desktop preflight.
   The preflight also verifies pinned cuteSV 2.1.3 and absolute packaged policy paths; this
   runtime correction does not change scientific thresholds.
4. Run **Selbsttest starten** before using any real sample. The self-test must finish successfully and produce its synthetic report bundle.
5. Confirm the user-writable WSL resource root (default `~/.local/share/ontseq/resources`) and install `GRCh38_GENCODE50_MANE1.5_v1`. The fast bundle status (manifest, pins, presence and declared sizes) must pass; use **Reparieren** for damaged resources. Run `ontseq references validate` when an explicit full SHA256 audit is required.
6. Choose `AML_LCWGS_GRCh38` or `AML_AS_111_GRCh38`. The latter resolves its controlled panel and analysis ROI automatically. Do not add a GRCh37 resource or explicit BED to these profiles.
7. Confirm the selected Windows BAM root is visible in WSL. For example, if input is on `P:`, `/mnt/p` must exist and point to the intended mapped drive. UNC paths are not guessed.
8. Select a **synthetic/non-patient** aligned BAM first together with its declared BAI (`sample.bam.bai` preferred; `sample.bai` also supported) and execute the complete Desktop path.
9. Verify that HTML/XLSX/JSON and `provenance/run.json` are produced in the configured output root, that stage statuses agree, and that failed/`NOT_RUN`/`NO_CALL` stages are not interpreted as negative biological findings.
10. Only after this workstation smoke succeeds should locally authorized research data be used.

Current deliberate limitations: the build is unsigned and RUO. The active desktop profiles and
automatic resource resolution in this work package are GRCh38-only. Legacy explicit path fields
remain readable for one compatibility release but are not mixed into profile runs. Cancellation
remains disabled until the backend can guarantee an explicit interruption state.
