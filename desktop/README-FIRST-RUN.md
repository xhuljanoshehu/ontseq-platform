# ONTSeq Desktop v0.3.4 — first workstation run

This checklist is for the unsigned engineering/research build only. It is not analytical or clinical validation.

1. Confirm Windows WSL2 is enabled and the configured distribution opens normally.
2. Extract the complete `ontseq-desktop-v0.3.4-win-x64-setup-engineering` bundle; do not separate `ONTSeq.Desktop.exe` from the bundled `runtime` directory.
3. Start `ONTSeq.Desktop.exe`, open **System einrichten**, and install the bundled ONTSeq Linux runtime. The installed runtime must report Core v0.3.4 and expose the reference-validation, target-coverage and component-selection capabilities checked by the Desktop preflight.
4. Run **Selbsttest starten** before using any real sample. The self-test must finish successfully and produce its synthetic report bundle.
5. Configure the exact reference used to align the BAM. Select the corresponding full FASTA/FAI source; the Desktop creates a content-addressed Reference-Lock and requires complete canonical chromosomes 1–22, X and Y for the selected GRCh37/GRCh38 build.
6. If Adaptive Sampling is used, configure the controlled target BED and its version. Do not substitute a different or unversioned target design.
7. Confirm the selected Windows BAM root is visible in WSL. For example, if input is on `P:`, `/mnt/p` must exist and point to the intended mapped drive. UNC paths are not guessed.
8. Select a **synthetic/non-patient** aligned BAM first together with its declared BAI (`sample.bam.bai` preferred; `sample.bai` also supported) and execute the complete Desktop path.
9. Verify that HTML/XLSX/JSON and `provenance/run.json` are produced in the configured output root, that stage statuses agree, and that failed/`NOT_RUN`/`NO_CALL` stages are not interpreted as negative biological findings.
10. Only after this workstation smoke succeeds should locally authorized research data be used.

Current deliberate limitations: the build is unsigned and RUO; GRCh37 QDNAseq/ACE CNV is the bundled real-tool-tested CNV path, GRCh38 CNV remains fail-closed/not requested in the Desktop until a separately pinned and tested hg38 annotation runtime is available; cancellation remains disabled until the backend can guarantee an explicit interruption state.
