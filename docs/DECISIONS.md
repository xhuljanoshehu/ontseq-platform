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

**Reason:** Rejected reads in an adaptive-sampling run **may** form a near-uniform
low-coverage whole-genome background — the population depth-based CNV methods assume —
but that is a hypothesis about the local assay, not an established property: it depends on
uniformity, GC behaviour, mappability and the usable genome fraction at the achieved yield,
none of which has been measured here. On-target depth is dominated by enrichment efficiency
rather than copy number and violates that assumption. Which basis the local assay should use is an open empirical question, so the
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

---

## ADR-017: Preflight is advisory, side-effect free, and never stricter than the run

**Decision:** `ontseq preflight` answers a run's preconditions before it starts. It creates
nothing — no envelope, no lock, no artifact, and not even the output directory; writability
is probed on the nearest existing ancestor. Five statuses are reported (`ok`, `FAILED`,
`WARNING`, `UNKNOWN`, `SKIPPED`) and only `FAILED` produces a non-zero exit. Version strings
are parsed by the adapters' own parsers, and a version lock is enforced only where a
*planned* stage would enforce it. Free disk space is reported as `UNKNOWN` unless the caller
supplies `--require-free-gb`.

**Reason:** Every gate preflight applies already exists inside the pipeline and already
fails closed. What was wrong was the timing: a POD5 run learns its Dorado model is missing
after the envelope exists and the lock is taken, an aligning run learns its reference index
is missing after intake. Moving those questions to the front converts hours of wasted
scheduling into two seconds.

That only helps if two things hold. **Preflight must never be stricter than the run**, or it
refuses work that would have succeeded — which is why the version parsing is shared rather
than duplicated, and why the alignment policy's samtools lock is not applied to an
aligned-BAM run that never aligns. And **preflight must leave nothing behind**, or a run
that was never started has an output tree afterwards suggesting it was.

Severity is derived, not curated: a missing tool is fatal exactly when a *required* stage
needs it, read from `StageSpec.required`. Sniffles serves only the optional SV stage, so its
absence warns that SV will record `NOT_RUN` rather than blocking a run that completes fine
without it. Curating a second list of "which tools really matter" would drift from the stage
graph the first time a stage changed.

Free space is the one check deliberately left unanswered. There is no measured relationship
in this repository between an input's size and the space a run consumes; it depends on the
lab's chemistry, depth and retention policy. A multiplier invented in the code would be
indistinguishable, in the output, from a validated figure. `UNKNOWN` with the raw number
attached is the honest form, and the check becomes real the moment somebody who has measured
it passes `--require-free-gb`.

**Revisit when:** Real run-size data exists for the intended assays, at which point the
space requirement can be derived from the declared input kind and depth rather than supplied
by the caller — and when a GPU-bearing host can be checked for basecalling capacity, which
is currently outside what this process can determine.

---

## ADR-018: Below its declared resolution a truth source is silent, and the partition says so

**Decision:** `evaluate()` takes the truth set's `resolution_bp`. Where the truth asserts
its *background* — its claim of absence — a called event shorter than that resolution
leaves the evaluable genome and is accounted as `truth_resolution_silent_bases`, a fourth
term in the partition identity. The affected calls are `NOT_ASSESSABLE`, neither confirmed
nor false positives, and the count is stated in a warning.

**Reason:** The documentation already said that a closed-world truth is *silent* below its
resolution rather than negative. The code did not implement it: it appended a warning
saying those calls "must not be read as false positives" while continuing to count them as
exactly that. A warning that contradicts the number beside it is worse than either alone,
because the number is what gets aggregated.

The rule is deliberately asymmetric, and the asymmetry is the substance of it. Resolution
limits what a source can **deny**, not what it can **affirm**. A karyotype read at 10 Mb
bands cannot rule out a 200 kb duplication, so it cannot make that call wrong. Where the
same karyotype explicitly reports a deletion it has made a positive claim, and a small call
agreeing with it is confirmed on its merits. Applying the rule to affirmations too would
quietly suppress true positives and depress sensitivity.

This exclusion is the only one in the design that *flatters* the caller: everything removed
here is something nobody can hold against it. That is why it is a named term in the
partition rather than a silent filter, why the warning states the base count and the number
of affected calls, and why specificity read from such a report has to be read together with
it. The alternative — leaving the false positives in — does not avoid the problem, it just
moves the dishonesty to the other side of the ledger and calls a caller wrong for seeing
something the truth was never able to look for.

**Revisit when:** A truth source can express per-region resolution rather than one figure
for the whole set — an array with variable probe density, or a karyotype where some
chromosomes were banded more finely than others. The single `resolution_bp` is then the
coarsest defensible summary of a structure the model cannot yet carry.

---

## ADR-019: A direction is only claimed when a pre-specified test supports it

**Decision:** `paired_detection_comparison()` takes an `alpha` and names a method in
`favours` only when McNemar's exact test is significant at it. The direction the discordant
counts happen to lean is reported separately as `observed_direction` and labelled
descriptive. `minimum_attainable_p_value` and `underpowered` state whether any observation
at that discordant count could have reached `alpha` at all.

**Reason:** The previous implementation set `favours` from the raw counts, so a 4-0 split
reported method A as the winner at p=0.125. That reads as a finding and is not one: with
four discordant pairs the smallest attainable two-sided exact p-value is 0.125, so no
possible outcome could have been significant. The study was decided before the data were
seen, and the report said the opposite.

Separating the two fields makes the distinction survive aggregation. "The counts lean
towards A" and "A is better" are different claims, and only the second needs a test; a
single field forces them into one word and the reader cannot recover which was meant.
`underpowered` is what separates *we compared them and found no difference* from *this
comparison could not have found one* — readings a bare non-significant p-value cannot tell
apart, and the second of which is a design fault rather than a result.

**Revisit when:** A validation study pre-registers a different inferential rule — a
non-inferiority margin, a one-sided test, a multiplicity correction across strata. `alpha`
is then one parameter of that rule rather than the whole of it.

---

## ADR-020: Event-level intervals are reported with their clustering, not without it

**Decision:** `aggregate()` reports how many specimens contributed the scored events, the
largest number contributed by one specimen, and a flag saying the intervals beside them are
anticonservative. A specimen-weighted detection rate is reported alongside the event-level
one. `paired_detection_comparison()` reports how many specimens the discordant pairs came
from and flags the p-value when there are fewer specimens than pairs.

**Reason:** Every interval in the aggregate is computed over events, and several events
routinely come from one specimen. Events within a specimen share its purity, its library,
its coverage and its artefacts, so they are not independent observations. Treating them as
independent narrows every interval: what is reported then describes a population of
independent events that does not exist. McNemar has the same problem one level down — it
treats each discordant pair as an independent coin flip.

The correct fix is a specimen-level endpoint or a cluster-robust test, and that is a
study-design decision this module has no standing to make: it depends on what the study is
trying to establish. What the module can do is refuse to present the problem as absent. The
specimen-weighted rate shipped here counts a specimen as a success only when every
assessable event in it was detected — a deliberately crude summary, labelled as one, whose
purpose is to make the gap between the two numbers visible rather than to be the endpoint.

**Revisit when:** An analytical validation is designed. At that point the endpoint, the
cluster-robust test (Durkalski's clustered McNemar or a GEE formulation) and the sample-size
calculation are pre-registered together, and these fields become inputs to that design
rather than caveats attached to an event-level number.

---

## ADR-021: A review is bound to content, appended, and honest about what it is not

**Decision:** `ontseq review record` appends one judgement — `accepted` or `rejected` — to
`review/review.log.jsonl` inside the run envelope. Every entry carries the SHA-256 of the
release bundle it was made against and the digest of the entry before it. Nothing is ever
overwritten. `ontseq review status` resolves the trail against the content currently on
disk into one of `pending`, `accepted`, `rejected`, `stale`, `broken` or `unreadable`, and
exits 0 / 2 / 6 on the same convention as `ontseq status`. `ontseq run` refuses, with its
own exit code 7, to write into an envelope whose latest review accepts its current content.

**Reason:** The pipeline could produce evidence and nobody could record that they had
examined it. `ReviewStatus` existed as a field on the ISCN proposal and nothing ever wrote
`REVIEWED` to it; the release bundle said `signature_status: "unsigned"` and left it there.
Without a trail, "this was reviewed" is a claim that cannot be checked afterwards, which is
exactly the claim a diagnostic workflow rests on.

**Bound to content, not to a directory.** A review that pointed at a path would keep
vouching for whatever later appeared there. Binding to the release bundle's digest means a
changed run makes the review `stale` — the judgement still stands for what it saw, and says
nothing about what is there now. This is the same rule the resume logic uses, for the same
reason: identity by content, never by name or timestamp.

**Refusing to re-run a reviewed envelope is the property that makes the rest mean
anything.** The lock stops two runs colliding now; content-addressed resume stops a stale
artifact being accepted. Neither notices that a human accepted this envelope yesterday and a
resumed run is about to rewrite what they accepted. That is deliberately *not* overridable
by a flag: a flag would be used, and the correct alternative costs nothing — a new run id,
after which the reviewed envelope keeps its review and the new run gets its own.

**Append-only, with a chain that shows it.** Each entry names the digest of its predecessor,
so removing, reordering or editing one breaks the chain there and at every later point. A
reviewer who accepted and then rejected leaves both facts behind; that is the whole purpose.

**And it says what it is not.** The reviewer identity is `asserted`: it comes from the
command line and nothing authenticated it. The chain is tamper-*evident*, not tamper-proof —
there is no key, so anyone who can rewrite the file can recompute the chain. It catches
accidental corruption and casual editing, and that is all it catches. Both facts are printed
on every human-readable report and carried as explicit `false` fields in the JSON, because a
record that quietly implies an authentication it never performed is worse than no record.

**Revisit when:** An identity provider exists and a signing key is available under a records
policy. `identity_source` then carries something other than `asserted`, the entry gains a
signature, and `signature_status` on the release bundle can stop saying `unsigned`. The log
format does not have to change for that; only the fields that are currently honest about
being empty.

---

## ADR-022: ClinVar annotates findings; it does not classify them

**Decision:** `ontseq_platform.knowledge` attaches ClinVar records to copy-number and
structural findings and stops there. Each annotation carries the source's assertion
*verbatim*, the vocabulary that assertion belongs to (`acmg_germline`), the record's own
origin, how the record matched (type and reciprocal overlap), NCBI's review status and star
rating, and the checksum of the exact release it came from. Nothing in the package sets
`reportable`, raises `confidence`, or translates one classification vocabulary into another.
Records whose origin does not match the assay's question are **kept and marked**, never
filtered out.

**Reason:** ClinVar was chosen as the knowledge base, and it is a good choice for what it
is — but what it is matters. ClinVar classifies under ACMG *germline* rules: *Pathogenic*
there means "this variant causes this inherited condition". The intended use of this
platform is AML, where the findings are *somatic*. A somatic driver is not "pathogenic" in
ClinVar's sense, and a somatic finding labelled *Pathogenic* in a report would read as
though it were.

So the assertion travels with its vocabulary attached and is never restated. This is the
same rule the CNV truth model already applies to `background_state`: a source's semantics
belong to the source, and losing them inverts the meaning of everything downstream.

**Origin is checked per record, not assumed per source.** ClinVar publishes `OriginSimple`
on every row, so alignment is computed rather than guessed. Where either side is silent —
including when the manifest does not declare whether the assay is looking for somatic or
germline variation — the alignment is `UNKNOWN` and says which side failed to declare. It is
never quietly treated as agreement.

**Mismatched records are kept.** A germline pathogenic deletion underlying a somatic finding
is a secondary finding a reviewer needs to see. Filtering it out would be a clinical decision
disguised as a technical filter, made by code, invisibly.

**A match is a measurement.** A 3 Mb ClinVar record inside a 90 Mb arm loss is a much weaker
statement than an exact coordinate match, and the two are distinguished by `MatchType` with
the reciprocal overlap beside it, rather than both appearing as "ClinVar: Pathogenic".

**The release is locked.** ClinVar republishes weekly. "ClinVar says Pathogenic" without
saying *which* ClinVar is not reproducible: the same BAM can produce different reports a
month apart with nothing recording why. The release checksum travels with every annotation,
as the reference lock and the cytoband resource already do.

**What this deliberately does not do:** decide reportability. That needs somatic criteria —
ELN, ICC, or a locally agreed gene list — which this repository does not have. Code that
promoted a ClinVar *Pathogenic* into a reportable finding would be inventing a clinical rule
nobody agreed to, in the least recoverable place: inside a report a physician is about to
sign.

**Where a reviewer meets this.** Annotations appear in the HTML report under a paragraph
stating how to read them, and in the Excel workbook on sheet `11_Annotations`. The workbook
needed more than the same sentence moved across: a spreadsheet has no prose, and a reviewer
who sorts by the assertion column sees *Pathogenic* one cell away from an event identifier.
So the reading rule occupies row 1 of the sheet itself — it survives the sheet being
exported, copied into another workbook or printed alone — and rows whose origin does not
match the assay's question are filled in the warning colour, because colour is read before
column ten is. The event sheets carry `db_records_matched`, a count and nothing more: a
reviewer who never opens sheet eleven would otherwise not know there was anything to open,
and a count cannot be misread as a classification of the finding.

**Revisit when:** A somatic knowledge source is added (OncoKB, CIViC, COSMIC) or ELN/ICC
criteria are encoded. `assertion_vocabulary` then carries more than one value and the
alignment check becomes genuinely useful in both directions, rather than mostly reporting
that a germline source was consulted for a somatic question.

## ADR-023: The runtime CNV lane is measured through the benchmark contract, not beside it

**Decision:** QDNAseq + ACE is the primary CNV candidate for lcWGS and Adaptive Sampling, and
it reaches the benchmark subsystem the same way any other candidate does — normalized into a
`CnvCallSet` by `call_set_from_qdnaseq_report()` and scored over the same evaluable-genome
mask. It is not given a private evaluation path, and the historical Lea adapter is not
promoted into the runtime in its place.

**Why not promote the historical adapter directly.** It has the strongest claim to
familiarity in the laboratory and the weakest claim to comparability: its numbers were
produced by a pipeline whose parameters, reference and tool versions are not recorded in a
form anything here can reproduce. Promoting it would import that gap into the runtime and
make the first "validated" lane the one whose provenance is least recoverable. Keeping it as
a comparator preserves what it is genuinely good for — telling us whether the new lane agrees
with what the laboratory has been reporting — without letting it decide anything.

**Why through the benchmark rather than beside it.** A candidate scored by its own harness is
scored on its own terms. The comparisons that eventually matter — QDNAseq against the
baseline control, against Spectre, against the historical values — are only meaningful if the
same mask, the same exclusion vocabulary and the same denominator apply to all of them. One
seam, so a method cannot look better by being measured differently.

**What is stated rather than inferred.** `data_basis` has no default. An adaptive-sampling run
holds an on-target population that is deeply and non-uniformly enriched and an off-target
population that may or may not behave like a shallow whole genome; pooling them is a third
case again. Placing a report in the wrong stratum because the code guessed would corrupt the
comparison silently, so the caller states it.

**What the mask absorbs.** QDNAseq drops bins it cannot correct and the runner segments only
chr1–22. Neither limit appears in a segment table — the rows that would say so are simply
absent — so every uncovered base is emitted as a declared no-call. Calls and declared
no-calls together account for each contig exactly, which is what makes the denominator
auditable; CI asserts that reconciliation on a real QDNAseq run rather than assuming it.

**What this does not decide.** Not reportability: `CnvCallSet` fixes `reportable` to `False`
and no argument changes it. Not a preferred runtime default beyond the current engineering
one — Spectre remains a legitimate benchmark comparator. And not the parameters: the bin
size, the ACE penalty of 0.6 and the ploidy grid arrive from the run's policy as engineering
values carried from the Lea evaluation. They are configured defaults, not measured ones, and
the policy file says so in as many words so that being the default cannot quietly become
being validated.

**Revisit when:** a real cohort with orthogonal characterization exists. Method selection is
that decision, and nothing before it — a good score on synthetic data included — is that
decision being made.
