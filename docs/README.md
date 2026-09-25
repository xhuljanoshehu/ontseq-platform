# Documentation map

The sources of truth for scope and limitations are [`ARCHITECTURE.md`](ARCHITECTURE.md),
[`EVIDENCE_BASE.md`](EVIDENCE_BASE.md) and [`CLINICAL_VALIDATION.md`](CLINICAL_VALIDATION.md).
[`THESIS_TRACEABILITY.md`](THESIS_TRACEABILITY.md) is historical context only. Everything here is
research use only; no document describes a clinically validated procedure.

## Start here

| Document | What it answers |
| --- | --- |
| [`SCHNELLSTART.md`](SCHNELLSTART.md) | How to get a first synthetic run (German) |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | Design principles, run envelope, stage graph, adapter and ISCN boundaries |
| [`PACKAGE_STRUCTURE.md`](PACKAGE_STRUCTURE.md) | Where code lives, which dependency rules are enforced, where new code goes |
| [`PIPELINE_EXECUTION.md`](PIPELINE_EXECUTION.md) | Runner, resume, preflight, watchfolder, execution semantics |
| [`COMPONENT_SELECTION.md`](COMPONENT_SELECTION.md) | Choosing provider and version per stage |
| [`REFERENCE_SYSTEM.md`](REFERENCE_SYSTEM.md) | Reference bundles, profiles, dictionary contracts, offline install |
| [`DATA_SECURITY.md`](DATA_SECURITY.md) | What may never leave controlled infrastructure |
| [`ROADMAP.md`](ROADMAP.md) | Milestones and open research lanes |
| [`DECISIONS.md`](DECISIONS.md) | Recorded design decisions |

## Evidence and validation

| Document | What it answers |
| --- | --- |
| [`EVIDENCE_BASE.md`](EVIDENCE_BASE.md) | Literature and technical rationale for tool choices |
| [`CLINICAL_VALIDATION.md`](CLINICAL_VALIDATION.md) | Analytical/clinical validation plan and dated validation-impact notes |
| [`BENCHMARKING.md`](BENCHMARKING.md), [`DILUTION_SERIES.md`](DILUTION_SERIES.md), [`QUANTITATIVE_MODEL.md`](QUANTITATIVE_MODEL.md) | Synthetic benchmarks, in-silico dilution, detection limits |
| [`ANALYTICAL_VALIDATION_V1_STATUS.md`](ANALYTICAL_VALIDATION_V1_STATUS.md), `CNV_VALIDATION_BLOCK*_STATUS.md` | CNV validation programme status |
| [`REQUIRED_INPUTS_AND_LITERATURE.md`](REQUIRED_INPUTS_AND_LITERATURE.md) | Inputs a real validation still needs |
| [`LEGACY_COMPARISON.md`](LEGACY_COMPARISON.md), [`THESIS_TRACEABILITY.md`](THESIS_TRACEABILITY.md) | Relation to the historical pipeline |

## Analysis lanes

| Area | Documents |
| --- | --- |
| Intake and QC | [`ALIGNED_BAM_MVP.md`](ALIGNED_BAM_MVP.md), [`TARGET_COVERAGE_ADAPTER.md`](TARGET_COVERAGE_ADAPTER.md) |
| Adaptive Sampling panel | [`PANEL_PROVENANCE.md`](PANEL_PROVENANCE.md), [`PANEL_REACHABILITY.md`](PANEL_REACHABILITY.md), [`GRCH37_INTEGRATION.md`](GRCH37_INTEGRATION.md) |
| Structural variants | [`SNIFFLES2_ADAPTER.md`](SNIFFLES2_ADAPTER.md), [`SV_EVIDENCE_LAYER.md`](SV_EVIDENCE_LAYER.md), [`CUTESV_TROUBLESHOOTING.md`](CUTESV_TROUBLESHOOTING.md) |
| Multi-caller research | [`MULTICALLER_INTEGRATION_STATUS.md`](MULTICALLER_INTEGRATION_STATUS.md), [`MULTICALLER_5A_STATUS.md`](MULTICALLER_5A_STATUS.md), [`MULTICALLER_REALTOOL_STATUS.md`](MULTICALLER_REALTOOL_STATUS.md), [`GATK_MUTECT2_ADAPTER.md`](GATK_MUTECT2_ADAPTER.md) |
| Methylation lane | [`METHYLATION_LANE.md`](METHYLATION_LANE.md), [`MODBAM_ADAPTER.md`](MODBAM_ADAPTER.md), [`MODKIT_PR709_BUILD.md`](MODKIT_PR709_BUILD.md) |
| MARLIN classification | [`MARLIN_CLASSIFICATION.md`](MARLIN_CLASSIFICATION.md), [`MARLIN_RUNTIME_FREEZE.md`](MARLIN_RUNTIME_FREEZE.md) |
| Methylation research studies | [`METHYLATION_MIXTURE.md`](METHYLATION_MIXTURE.md), [`METHYLATION_HOLDOUT.md`](METHYLATION_HOLDOUT.md), [`METHYLATION_VALIDATION.md`](METHYLATION_VALIDATION.md) and their `*_STATUS.md` files |
| Methylation data sourcing | [`METHYLATION_DATA_CANDIDATES.md`](METHYLATION_DATA_CANDIDATES.md), [`METHYLATION_DATA_SAFETY_REVIEW.md`](METHYLATION_DATA_SAFETY_REVIEW.md), [`METHYLATION_ADDITIONAL_DATA_SEARCH.md`](METHYLATION_ADDITIONAL_DATA_SEARCH.md), [`METHYLATION_NEXT_SOURCE_REVIEW.md`](METHYLATION_NEXT_SOURCE_REVIEW.md), [`METHYLATION_SOURCE_REVIEW_STATUS.md`](METHYLATION_SOURCE_REVIEW_STATUS.md), [`METHYLATION_INTAKE_STATUS.md`](METHYLATION_INTAKE_STATUS.md) |

## Reports and interfaces

| Document | What it answers |
| --- | --- |
| [`REPORT_INTERFACE.md`](REPORT_INTERFACE.md), [`INTERACTIVE_REVIEW_WORKSPACE.md`](INTERACTIVE_REVIEW_WORKSPACE.md) | Reviewer report and live workspace |
| [`VISUALIZATION_DATA_MATRIX.md`](VISUALIZATION_DATA_MATRIX.md), [`VISUALIZATION_TECHNOLOGY_EVALUATION.md`](VISUALIZATION_TECHNOLOGY_EVALUATION.md) | What is plotted from which contract |
| [`DESKTOP_FIRST_RUN.md`](DESKTOP_FIRST_RUN.md), [`DESKTOP_API_CONTRACT.md`](DESKTOP_API_CONTRACT.md), [`DESKTOP_ROADMAP.md`](DESKTOP_ROADMAP.md) | Windows Desktop and its local service contract |

## Releases and repairs

[`INTEGRATION_0_7.md`](INTEGRATION_0_7.md), [`RELEASE_0_7_1.md`](RELEASE_0_7_1.md),
`merge-repairs/`, `changes/`, `gatk_adapter/` and `superpowers/` (design specs and plans of
individual work packages) record how specific changes were made; the changelog in the
repository root is the running record.
