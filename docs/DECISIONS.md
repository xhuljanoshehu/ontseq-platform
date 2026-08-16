# Architecture decision log

## ADR-001: Keep computation local by default

**Decision:** Raw and derived human genomic data remain on approved institutional compute.

**Reason:** Genomic files are sensitive and potentially identifying. Cloud orchestration is
metadata-only unless explicitly approved.

## ADR-002: Use Snakemake for the initial workflow

**Decision:** Use a file-oriented Snakemake DAG for the first validated implementation, with
every scientific tool isolated behind an adapter.

**Reason:** Snakemake supports explicit file dependencies, isolated software environments,
container execution, provenance and testable modular rules. The decision is independent of any
specific caller or prior implementation.

**Revisit when:** The operational target requires native cloud orchestration, a centrally
maintained Nextflow platform, or evidence shows that another engine materially improves
reproducibility or supportability. Scientific adapters and contracts must remain portable.

## ADR-003: Use one normalized event model

**Decision:** HTML, Excel, JSON and ISCN proposals are rendered from the same validated
event model.

**Reason:** Parallel report logic causes silent discrepancies. One source of truth makes each
displayed result traceable to caller evidence.

## ADR-004: Treat ISCN as a proposal until validated

**Decision:** No automated notation can be released without expert review.

**Reason:** Syntax alone does not establish semantic correctness, clonality, uncertainty or
clinical validity.

## ADR-005: Treat the thesis as context, not specification

**Decision:** Lea Evers' master's thesis remains a historical comparison source only. Public
primary literature, official documentation, controlled benchmarks and the intended-use
statement govern architecture and tool selection.

**Reason:** A single exploratory dataset cannot establish general performance, clinical
thresholds or a durable software architecture. Preserving traceability is useful, but inheriting
its tools or parameters without independent validation would create avoidable anchoring bias.

## ADR-006: Select callers through assay-specific benchmarks

**Decision:** No CNV, SV or fusion caller is a production default until it passes a predeclared
benchmark for the relevant assay mode, coverage, tumor/blast fraction, reference build and event
class. Candidate lists are not rankings.

**Reason:** Published results are not interchangeable across germline, tumor-normal, tumor-only,
low-coverage WGS and adaptive-sampling settings. Caller agreement is an evidence feature, not a
truth label.

**Initial candidates:** Compare ichorCNA, QDNAseq + ACE and Spectre for copy number; compare
Sniffles2 with at least one cancer-aware or independent long-read SV method. Evaluate SAVANA and
Severus only in data regimes matching their assumptions. See `docs/EVIDENCE_BASE.md`.

## ADR-007: Integrate Sniffles2 first as conservative candidate evidence

**Decision:** Pin Sniffles2 v2.8.0 in the executable environment and normalize its PASS records,
but force every event to `unclassified`, `reportable: false` under the technical policy.

**Reason:** Sniffles2 has strong primary evidence and a stable long-read VCF interface, making it
useful for exercising the complete adapter and reporting boundary. Published cross-platform
performance does not constitute AML, tumor-only, low-coverage or adaptive-sampling validation.

**Privacy and interpretation constraints:** Request symbolic alleles, never request read-name
output, count rejected records, map an empty accepted set to `NO_CALL`, and never infer a fusion or
ISCN assertion directly from a BND record.

## ADR-008: Score copy number per base, not per matched event

**Decision:** CNV comparison uses an exact breakpoint partition of the genome and a
base-pair-weighted per-state confusion matrix. Event-level detection is derived from
base-level concordance with a many-to-many rule. The one-to-one event matcher in
`benchmark.py` remains the SV contract and is unchanged.

**Reason:** Segment boundaries are an artifact of bin size and segmentation algorithm, not a
biological claim. One-to-one matching scores a correctly detected deletion emitted as three
adjacent segments as one true positive and two false positives. Event counting also discards
event size, so a whole-chromosome gain and a 200 kb duplication weigh equally. Base-level
scoring is invariant to both.

**Revisit when:** A representation-aware CNV comparator handling segmentation equivalence
directly becomes available, or allele-specific copy number requires a richer state space.

## ADR-009: Make the evaluable genome explicit and account for every excluded base

**Decision:** No CNV metric is computed without an explicitly constructed observability mask.
Excluded bases are attributed to a closed vocabulary of reasons whose counts sum exactly to
the total removed. Truth events in unobservable regions are `NOT_ASSESSABLE`, never false
negatives. When several methods are compared, the union of all their no-call regions is
removed from the shared mask before any of them is scored.

**Reason:** A sensitivity figure is meaningless without stating what fraction of the genome it
applied to. Counting unmatched calls as false positives assumes universal observability, which
converts known blind spots into fabricated error rates. Scoring each method on its own mask
rewards a cautious method for declining to call in hard regions.

## ADR-010: Truth sets declare their own background semantics and resolution

**Decision:** Every truth set and call set declares `background_state` (`neutral` = closed
world, `no_call` = open world) with no default, and closed-world truth additionally declares
`resolution_bp`. Breakpoint accuracy is withheld when a truth event's boundary uncertainty
exceeds the configured limit. A call set's own declared uncertainty is ignored.

**Reason:** FISH interrogates probes and says nothing elsewhere; an array asserts neutral
within its probe map; a karyotype asserts neutral only at band resolution. Encoding these
identically is silent and severe: an open-world truth treated as closed-world manufactures a
false positive for every genuine finding outside its scope. Reporting a breakpoint error
against `5q13` measures the width of a Giemsa band, not the caller. Letting a caller's own
uncertainty suppress the metric would let it excuse its own error.

## ADR-011: Ship a transparent baseline CNV caller as an experimental control

**Decision:** Implement `ontseq-baseline-readdepth`, a dependency-free binned read-depth
segmenter, explicitly marked non-reportable. It is not a candidate for production use.

**Reason:** A benchmark harness that has never scored a real call set is untested
infrastructure. The baseline closes the loop in CI without any external binary, reference
genome or genomic data, and supplies the null model a candidate method must beat before its
additional complexity is justified.

**Revisit when:** A version-pinned external caller runs in the executable environment. The
baseline then stays as a control rather than being removed.

## ADR-012: Treat the CNV data basis as a mandatory stratification key

**Decision:** Every call set declares `data_basis`. Adaptive-sampling off-target reads,
adaptive-sampling on-target reads, combined output and separate low-coverage WGS are never
pooled into one benchmark stratum.

**Reason:** Rejected reads in an adaptive-sampling run form a near-uniform low-coverage
whole-genome background, which is the population depth-based CNV methods assume. On-target
depth is dominated by enrichment efficiency rather than copy number and violates that
assumption. Which basis the local assay should use is an open empirical question, so the
architecture stays agnostic and makes the comparison possible instead of choosing
prematurely.
