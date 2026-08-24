# ONTSeq Desktop changelog

## 0.3.4-engineering

- combine Core 0.3.4 with the working Desktop/reference-integrity hardening from v0.2.2
- report Desktop and assembly version 0.3.4 and install the bundled runtime under `runtime-v0.3.4`
- require canonical GRCh37/GRCh38 reference identity during setup and again before each Desktop run
- create content-addressed Reference-Locks from the selected FAI SHA-256 and publish new locks atomically
- preserve the previously active reference lock if a newly generated lock fails validation
- bind aligned-BAM intake resume identity to BAM, declared BAI, reference lock, pipeline version and git commit
- verify the explicitly declared BAI with `samtools idxstats -X` rather than allowing implicit index discovery
- pass Core 0.3.4 target-coverage policy and component selection through the Desktop service path
- embed the exact runtime git commit in the packed Linux backend
- package a relocatable Linux/R runtime and test it on stock Ubuntu without host Python or R
- build and test the self-contained Windows x64 WPF shell before creating the engineering bundle
- retain the RUO boundary; GRCh38 QDNAseq/ACE CNV and cooperative cancellation remain deliberately unavailable

## 0.1.3-engineering

- support both common BAM index names: `sample.bam.bai` and `sample.bai`
- prefer `sample.bam.bai` deterministically when both files exist
- validate the selected index again inside the local service boundary
- add Windows and backend regression tests for index selection

## 0.1.0-engineering

- initial .NET 10 WPF operator shell
- BAM selection and sample-ID suggestion
- GRCh37/GRCh38 and lcWGS/Adaptive Sampling profiles
- automatic WSL backend preflight and launch
- local API start/status integration
- live stage status from persisted run provenance
- HTML/XLSX/result-folder launch buttons
- RUO boundary and disabled cancellation pending backend support
