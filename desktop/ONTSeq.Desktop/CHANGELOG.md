# ONTSeq Desktop changelog

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
