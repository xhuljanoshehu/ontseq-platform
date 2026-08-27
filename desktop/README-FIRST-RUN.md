# ONTSeq Desktop v0.4.1 — first workstation run

This checklist is for the unsigned engineering/research build only. It is not analytical or clinical validation.

1. Confirm Windows WSL2 is enabled and the configured distribution opens normally.
2. Extract the complete `ontseq-desktop-v0.4.1-win-x64-setup-engineering` bundle; do not separate `ONTSeq.Desktop.exe` from the bundled `runtime` directory.
3. Start `ONTSeq.Desktop.exe`, open **System einrichten**, and install the bundled ONTSeq Linux runtime. The installed runtime must report Core v0.4.1 and expose the reference-validation, target-coverage and component-selection capabilities checked by the Desktop preflight.
4. Run **Selbsttest starten** before using any real sample. The self-test must finish successfully and produce its synthetic report bundle.
5. Configure the exact reference used to align the BAM. Select the corresponding full FASTA/FAI source; the Desktop creates a content-addressed Reference-Lock and requires complete canonical chromosomes 1–22, X and Y for the selected GRCh37/GRCh38 build.
6. If Adaptive Sampling is used, configure the controlled target BED and its version. Do not substitute a different or unversioned target design.
7. Confirm the selected Windows BAM root is visible in WSL. For example, if input is on `P:`, `/mnt/p` must exist and point to the intended mapped drive. UNC paths are not guessed.
8. Select a **synthetic/non-patient** aligned BAM first together with its declared BAI (`sample.bam.bai` preferred; `sample.bai` also supported) and execute the complete Desktop path.
9. Verify that HTML/XLSX/JSON and `provenance/run.json` are produced in the configured output root, that stage statuses agree, and that failed/`NOT_RUN`/`NO_CALL` stages are not interpreted as negative biological findings.
10. Only after this workstation smoke succeeds should locally authorized research data be used.

Current deliberate limitations: the build is unsigned and RUO. QDNAseq/ACE CNV is bundled for both GRCh37 (`QDNAseq.hg19` 1.36.0) and GRCh38 (`QDNAseq.hg38` 1.2.0 pinned to upstream commit `cf7c07e39de0ac64a9c38cb030cba4626e2aae83`). Both assembly lanes are engineering-tested with deterministic synthetic truth but are not analytically or clinically validated on a real cohort. Cancellation remains disabled until the backend can guarantee an explicit interruption state.
