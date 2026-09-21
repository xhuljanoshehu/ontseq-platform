# Validation impact and integration boundary

Date: 2026-09-20. Adapter target schema: 0.1.0.
Research Use Only. No automatic clinical release.

## Implemented and exercised

- Local typed configuration, disabled-by-default opt-ins, API and standalone CLI.
- Separate paired/tumor-only command construction, no fallback between modes.
- Runtime argument planning against pinned GATK 4.6.2.0 source interfaces.
- Local hash/sidecar/reference dictionary checks and modeled native validation steps.
- Loss-conscious small-variant extraction, filter retention and explicit status/provenance.
- Synthetic adapter contracts, CLI subprocess tests and a harmless native Python transport test.

This is **not** a platform-wide qualification declaration. See delivery verification files
for the exact commands, environment, test count, exit codes and scope of each check.

## Original standalone-delivery boundary

This table describes the original standalone delivery. GitHub upload and software-gate status
are superseded by PR #85 checks on the exact reviewed revision. Native GATK and analytical
validation remain unqualified; see REVIEW.md.

### Historical status

| Area | Status | Consequence |
| --- | --- | --- |
| Native GATK 4.6.2.0 + Java 17 end-to-end execution | NOT_RUN | Real-tool interoperability unqualified. |
| Real BAM/FASTA/VCF/index fixtures | NOT_RUN | Synthetic contracts do not validate native parsing or indexing. |
| Full upstream `make safety` | NOT_RUN | Repository safety gate remains open. |
| Full upstream `make versions` | NOT_RUN | Repository version/configuration gate remains open. |
| Full upstream `make lint` / strict mypy | NOT_RUN | Style/type gates remain open; tools not present locally. |
| Full upstream `make test` | NOT_RUN | Existing upstream regression suite has not been exercised. |
| Canonical multi-caller catalog and runtime registry | NOT_IMPLEMENTED | Adapter must be invoked explicitly; no automatic selection. |
| Canonical small-variant schema and report integration | NOT_IMPLEMENTED | Isolated candidate JSON only; no existing report alteration. |
| Docker/Conda immutable environment lock | NOT_DELIVERED | JAR hash/version and JVM-major check are not a full environment lock. |
| Read/reference/resource biological provenance | OPERATOR_ATTESTATION | A checksum/label does not validate origin or specimen suitability. |
| ONT sensitivity/precision/VAF bias/LoD95 | NOT_ASSESSED | No performance or clinical claims. |
| Clinical validation/release | NOT_AUTHORIZED | Human review and institutional validation remain necessary. |
| GitHub branch/commit/PR/merge | NOT_CREATED | Delivery is local files plus an additive patch only. |

## Biological-output impact

This patch does not alter existing ONTSeq execution or reports. Explicit use of the new
CLI can create research variant candidates. Interpretation risk includes ONT-specific
errors, low-allele-fraction performance, matched-normal suitability, population/PoN
mismatch, ambiguous representation of indels, and default caller filtering assumptions.
The adapter cannot remove these risks by adding metadata or observing caller agreement.

In particular, a small-variant caller designed for germline detection must not be treated
as interchangeable with a somatic caller solely because both emit VCF. A future comparison
needs intended-use-matched callers, callable-region definitions, normalized equivalent
alleles, a locked truth set, false-positive/false-negative adjudication and independent
holdout specimens. No majority-vote truth construction is part of this patch.

## Next engineering gates, in order

1. Apply and review the additive patch in a full checkout; execute all upstream gates.
2. Qualify native CLI/output contracts on genuinely valid, generated synthetic resources.
   Verify tumor-only and paired modes, zero-call cases, disabled/missing resources, native
   filter behavior, contamination output, BAM validation, timeouts, and resource limits.
3. Design and review a canonical small-variant analytical domain, schema version and
   registry/runner bridge; retain real-tool-qualified=false until evidence supports a change.
4. Freeze chemistry/basecaller/alignment/runtime/resource policies and prospectively
   register an assay-appropriate truth/holdout validation program. Separate coverage,
   allele fraction, tumor fraction, copy-number state and genomic context.
5. Only independently established validation evidence may support changes to intended
   use or reportability. The current schema does not permit clinical release.
