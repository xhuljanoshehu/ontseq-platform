# ONTSeq Platform

Safety-first foundation for an automated, single-sample Oxford Nanopore workflow that
produces structured JSON, a self-contained HTML report, an Excel workbook, and an
expert-reviewable ISCN proposal.

> **Research use only. Not clinically validated.** The repository deliberately separates
> automated evidence processing from clinical interpretation and release.

## Kurzfassung auf Deutsch

Dieses Repository ist die belastbare technische Grundlage für die spätere Software:

- genau eine Probe pro reproduzierbarem Auftrag;
- POD5, unaligned BAM oder aligned BAM als definierte Eingangspfade;
- modulare Snakemake-Architektur für QC, CNV, SV, Fusionen, Annotation und Reporting;
- ein gemeinsames, streng validiertes Ergebnisformat;
- automatisch erzeugte HTML-, Excel- und JSON-Ausgaben;
- ein nachvollziehbarer **ISCN-Vorschlag**, niemals eine automatische klinische Freigabe;
- keine Patienten- oder Genomdaten in GitHub.

The current milestone implements the contracts and reporting layer with synthetic data.
Bioinformatics tools are represented by versioned interfaces and analysis profiles; the
validated production rules will be integrated from approved pipeline source code in later
milestones.

## Why this repository exists

The architecture is an independent, evidence-led implementation based on peer-reviewed
literature, public benchmarks, official tool documentation and explicit intended-use
requirements. Lea Evers' 2026 master's thesis is retained only as useful local context. It does
not select algorithms, parameters or reportability thresholds, and no source code from the
thesis project is copied into this repository.

## What works now

```bash
python -m pip install -e .
ontseq demo --output-dir results/demo
```

This creates:

- `SYNTHETIC_AML_001.result.json`
- `SYNTHETIC_AML_001.report.html`
- `SYNTHETIC_AML_001.results.xlsx`

All demo values and coordinates are synthetic and must not be interpreted biologically.

Validate inputs or render an existing result contract:

```bash
ontseq validate-manifest examples/manifests/demo.yaml
ontseq validate-result results/demo/SYNTHETIC_AML_001.result.json
ontseq render results/demo/SYNTHETIC_AML_001.result.json --output-dir results/rerendered
```

Optional Snakemake demo after installing Snakemake:

```bash
snakemake --snakefile workflow/Snakefile --configfile workflow/config/demo.yaml --cores 1
```

## Repository map

| Path | Responsibility |
| --- | --- |
| `src/ontseq_platform/` | Typed contracts, ISCN proposal logic, HTML/Excel exporters, CLI |
| `workflow/` | Snakemake entry point, modular rules, environments, runtime profiles |
| `configs/` | Assay profiles and non-clinical default configuration |
| `schemas/` | Versioned JSON Schemas for manifests and results |
| `tests/` | Contract, safety, ISCN-subset, HTML and Excel regression tests |
| `docs/` | Architecture, security, validation, roadmap and thesis traceability |
| `.github/` | CI, dependency updates, issue forms and review templates |

## Planned production modules

1. Input integrity and reference-build validation
2. Dorado basecalling and Minimap2 alignment when starting from POD5/uBAM
3. Cramino QC plus coverage and adaptive-sampling target QC
4. Benchmark-gated CNV adapters for ichorCNA, QDNAseq + ACE and Spectre
5. Sniffles2 evidence plus independently evaluated somatic/consensus SV candidates
6. SnpEff/SvAnna annotation and fusion evidence normalization
7. Build-aware cytoband mapping and an authorized ISCN 2024 conformance test suite
8. Human review, signature, immutable release bundle and audit trail

Future optional modules are reserved for small variants, RNA fusions and modified-base/
methylation analysis. They will not be mixed into the karyotyping path without separate
assay validation.

## Data boundary

POD5, FASTQ, BAM/CRAM, patient VCFs, clinical reports, direct identifiers and reference
bundles are prohibited from Git. Runtime data belongs on an approved on-premises storage
system. See [Data security](docs/DATA_SECURITY.md).

## Documentation

- [Architecture](docs/ARCHITECTURE.md)
- [Evidence base and tool-selection record](docs/EVIDENCE_BASE.md)
- [Master-thesis traceability](docs/THESIS_TRACEABILITY.md)
- [Clinical validation plan](docs/CLINICAL_VALIDATION.md)
- [Roadmap](docs/ROADMAP.md)
- [Data security](docs/DATA_SECURITY.md)
- [Required inputs and literature](docs/REQUIRED_INPUTS_AND_LITERATURE.md)

## Project status and license

Version `0.1.0` is a testable repository foundation, not a finished diagnostic pipeline.
No open-source license has been assigned. Keep the repository private until intellectual
property, institutional governance and intended medical-device use have been reviewed.
