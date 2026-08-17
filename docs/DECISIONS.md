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

## ADR-013: One command per run, into a self-describing run envelope

**Decision:** A run is invoked once per sample (`ontseq run`), writes everything into a single
directory keyed by run and sample, and records every artifact by its envelope-relative path
with a SHA-256. Writes are atomic. A watch-folder or queue, if it is ever needed, becomes a
caller of this command rather than a different execution path.

**Consequence:** The Snakemake workflow was reduced from five per-stage rules to one rule that
calls `ontseq run`. It had become a second execution path: flat files instead of a run
envelope, no per-artifact tool versions, no release bundle, and mtime-based resume — none of
it exercised by CI. Snakemake keeps cluster submission and multi-sample fan-out; it no longer
re-derives the stage graph.

**Reason:** The alternative — a long-lived service that discovers work — couples scheduling to
execution and makes a run irreproducible by hand. One command per run keeps the unit of work
identical whether a person, a cron job or a future watcher triggers it. Absolute paths are
excluded by construction rather than by review, because a reviewer artifact leaking the source
BAM location violates the data boundary in `docs/DATA_SECURITY.md`. Atomic writes exist so
that an interrupted run leaves either the previous artifact or none: a truncated artifact is
worse than a missing one, because resume would accept it.

## ADR-014: Resume on content, never on timestamps

**Decision:** A completed stage is reused only when its signature — upstream artifact
checksums, its own parameters, resolved tool versions and the fingerprints of external inputs
— is unchanged, *and* every artifact it claimed still verifies byte for byte. Anything else
re-runs. CI proves the property by running the same pipeline twice and failing if any stage
re-executes.

**Reason:** Timestamp-based resume silently accepts an artifact produced under different
parameters or a different tool version, which is how two incompatible results end up inside
one envelope with nothing recording that it happened. Resume is an optimisation, and an
optimisation that can corrupt a result is not one. Verifying the artifacts as well as the
signature closes the remaining gap where the inputs are unchanged but the output was edited,
truncated or deleted.

## ADR-015: Verification status is per adapter and is only claimed once CI runs the binary

**Decision:** Every stage declares a machine-readable `VerificationStatus`
(`verified_with_real_tool`, `verified_pure_python`, `unverified_adapter`, `not_implemented`).
A stage that completes on an unverified adapter is named in the run report and the release
bundle. The status is flipped in the same commit that adds the CI job executing the real tool
— never before.

**Reason:** "Implemented" and "known to work" are different claims, and the difference matters
most exactly where it is easiest to elide: an adapter that has never met its binary looks the
same in a diff as one that has. Dorado cannot be executed in this repository's environment, so
basecalling stays `unverified_adapter` and says so in its own output rather than in a comment
someone has to find. Alignment was `unverified_adapter` until the synthetic alignment-lane job
existed; the flip and the job landed together.

**Revisit when:** A GPU environment with a real Dorado model becomes available to CI, at which
point basecalling can be verified the same way alignment was.


## ADR-016: One run at a time per envelope, enforced by an exclusive lock

**Decision:** A run holds `.ontseq-run.lock` at its envelope root for its whole duration,
acquired with `O_CREAT | O_EXCL`. A lock whose holder is a dead process *on this host* is
reclaimed and the reclaim is recorded in the run report; a lock held from another host is
never reclaimed automatically. `ontseq run` exits 4 when the envelope is in use, distinct
from its failure exit code.

**Reason:** Nothing else in the design catches this. Atomic writes prevent a truncated
artifact; content-addressed resume prevents a stale artifact being accepted. Neither notices
a second process rewriting the same run report from a different set of stage records — the
loser vanishes from the history with nothing recording that it existed. That is invisible
rather than loud, which is what makes it worth preventing in code rather than in a warning.
The distinction only starts to matter once something other than a person triggers runs, so
it is a prerequisite for the watch folder rather than a nicety alongside it.

Cross-host reclaiming is refused because on shared storage a crashed remote run and a live
one look identical from here; guessing wrong produces exactly the concurrent-run state the
lock exists to prevent. PID reuse is likewise resolved towards refusing: a delayed run costs
time, a wrongly admitted one costs the envelope.

**Revisit when:** Runs need to be distributed across hosts sharing one envelope root, at
which point PID-and-host liveness is no longer sufficient and a lease with an explicit
heartbeat becomes necessary.
