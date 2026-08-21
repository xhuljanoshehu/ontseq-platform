# ONTSeq Desktop changelog

## 0.2.2-engineering

- validate a selected GRCh37/GRCh38 FAI against the complete canonical chromosome set before
  saving it as the active reference
- derive the dictionary ID and lock filename from the FAI SHA256 fingerprint, preventing a
  same-day dictionary replacement from reusing stale intake results (the FAI fingerprint does
  not prove underlying FASTA base identity)
- revalidate saved locks during setup and display concrete validation failures
- publish a newly validated lock atomically so a failed or changing FAI cannot overwrite the
  active lock
- detect a v0.2.1 backend as outdated and direct the operator to update the bundled Runtime
- preserve the detailed aligned-BAM intake failure and its `manifest/intake.json` diagnostic
  artifact in the Desktop job status
- keep strict BAM/reference dictionary matching; regional BAM fixtures still require the full
  reference dictionary used for their original alignment

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
