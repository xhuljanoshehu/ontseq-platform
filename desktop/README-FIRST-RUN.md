# ONTSeq Desktop — first workstation run

This checklist is for engineering/research setup only.

1. Confirm Windows WSL2 is enabled and the configured distribution opens normally.
2. Confirm `ontseq` is installed inside that distribution and `ontseq --help` works.
3. Copy `desktop.settings.example.json` to `%LOCALAPPDATA%\ONTSeq\desktop.settings.json` and replace every reference/assay placeholder with locally approved paths.
4. Confirm the selected Windows BAM root is visible in WSL. For example, if input is on `P:`, `/mnt/p` must exist and point to the intended mapped drive.
5. Start ONTSeq Desktop, select a **synthetic/non-patient** BAM first, and run the complete UI path.
6. Verify that HTML/XLSX/JSON are produced in the configured Windows output root and that stage statuses match `provenance/run.json`.
7. Only after this smoke succeeds should locally authorized research data be used.

Do not treat a successful first-run smoke as analytical or clinical validation.
