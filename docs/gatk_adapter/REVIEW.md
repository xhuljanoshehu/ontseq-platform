# PR #85 engineering review

Research Use Only. Human review required. No automatic clinical release.

## Reviewed scope

Review base: af9a1ba35f12fb9fbbeb3774789b259cd4794c2d.
Main integration base: 2c7aa2bf1a2d01f3f51aa75232c608407f8329b9.

Both explicitly invoked adapters are retained. `gatk_mutect2` is the legacy,
interval-bounded CLI; `gatk_adapter` is the separate versioned Pydantic API/CLI.
Their manifests, outputs and missing-resource policies are not interchangeable.
Neither is registered for automatic production analysis or clinical reporting.

## Findings and repairs

- Uploaded standalone files failed strict typing and repository lint/format gates.
  Add explicit literal, union, tuple, numeric and optional-input types; retain the
  existing Pydantic schemas and runtime opt-in defaults.
- Both version gates accepted some suffixed build identifiers as the pinned version.
  Add failing regressions before correction and reject extended/local/snapshot IDs
  before Mutect2 is called. Do not relax the fixed GATK 4.6.2.0 baseline.
- Focused CI omitted the newly added package. Include both implementations, tests
  and versioned schemas, with formatting enforced rather than merely displayed.
- Correct stale standalone-upload statements and add the root changelog entry.

## Verification required for this revision

The repair workflow runs the newly added regressions in RED before changes and
then the full project gates: `make safety`, `make versions`, `make lint`,
`make test`, schema freshness and wheel-resource verification. The final PR also
runs the repository's usual CI; only exact-head results are evidence for merging.
No ignore, skip, exclusion or reportability threshold is added by this repair.

## Unchanged qualification limits

Software contract tests use synthetic stand-ins. They do not execute GATK and
cannot establish native BAM/JAR interoperability, ONT sensitivity, specificity,
VAF bias, LoD, somatic origin or clinical validity. The versioned extractor
preserves native anchored alleles without left normalization or locus-by-locus
REF validation. Legacy and new outputs must not be joined by naive string
identity or treated as independent confirmation. All original native artifacts
remain the review authority. Real-tool and assay qualification remain separate.
