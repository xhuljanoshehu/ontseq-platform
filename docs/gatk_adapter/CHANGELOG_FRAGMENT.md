# Proposed changelog fragment (not yet inserted into upstream CHANGELOG.md)

## Unreleased — experimental GATK adapter patch

Added an isolated, opt-in `ontseq_platform.gatk_adapter` API/CLI targeting the explicit
GATK 4.6.2.0 / Java 17 compatibility baseline. Includes paired/tumor-only command plans,
content-fingerprinted input handling, native validation steps, Mutect2/FilterMutectCalls
with statistics retention, optional contamination estimation, candidate-only extraction,
synthetic contract tests, isolated schemas and validation-impact documentation.

No changes to package version 0.8.2, canonical runner/registry, CNV/SV public schemas,
clinical reportability, consensus logic, license or deployment. No real-tool or analytical
qualification is claimed. Not committed or merged into the upstream repository.
