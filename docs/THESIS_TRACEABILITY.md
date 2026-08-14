# Traceability to Lea Evers' master's thesis

Source reviewed: Lea Evers, master's thesis on automated ONT karyotyping, dated
14 January 2026. This document records conceptual traceability only; thesis source code is
not present in this repository.

| Thesis requirement or result | Repository response | Status |
| --- | --- | --- |
| Low-coverage ONT WGS | Separate `lcwgs` assay profile | Contract ready; caller adapter pending |
| Single-sample analysis | One manifest and run envelope per sample | Implemented |
| Fully automated cellularity selection | ACE penalty represented and alternative fits retained | Profile ready; validation pending |
| Modular method comparison | Adapter boundary and modular Snakemake layout | Foundation implemented |
| Cramino quality metrics | Versioned QC contract | Model implemented; adapter pending |
| QDNAseq + ACE CNV | Primary lcWGS CNV profile with 100/500/1000 kb bins | Profile implemented; rules pending |
| Spectre comparison | Comparator marked benchmark-only | Profile implemented; rules pending |
| Sniffles2, cuteSV and NanoSV | Primary, secondary and research-comparator roles | Profile implemented; adapters pending |
| SnpEff/SvAnna annotation | Planned adapter and normalized fusion evidence | Contract ready |
| Adaptive-sampling fusion BED | Versioned BED required in manifest | Implemented |
| Cytoband merging | Build-aware service boundary | Planned |
| ISCN 2024 output | Explicitly unvalidated proposal with source-event traceability | Limited subset implemented |
| Self-contained HTML | Synthetic, self-contained HTML renderer | Implemented |
| Excel output | Ten-sheet review workbook | Implemented beyond thesis |
| Runtime/provenance | Tool, parameter, commit and reference records | Implemented |
| Approximately 30-35 min per sample in thesis | Benchmark contract; no performance claim in this repo yet | Pending reproduction |

## Evidence limitations carried forward

- The reported ACE penalty of `0.6` is a dataset-specific starting point, not a universal
  threshold.
- Smaller CNV bins produced more segmentation and require resolution-aware validation.
- Only five adaptive-sampling samples were available for SV comparison, without robust SV
  ground truth.
- Caller concordance measures agreement, not truth.
- The thesis started from pre-aligned BAM; Dorado/Minimap2 integration is a new validation
  domain.
