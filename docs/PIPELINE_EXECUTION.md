# Pipeline Execution

Research use only. Nothing described here is validated for diagnostic use.

This document describes how one sample is taken from raw input to a checksummed release
bundle in a single command, and — just as important — what that command does *not* prove.

```
ontseq run <manifest.json> --reference-lock <lock.json> --reference-fasta <reference.fa> \
  --run-id <RUN>
```

The scope of automation ends where interpretation begins. The pipeline produces evidence
and provenance; a human decides what any of it means.

---

## 1. The stage graph

Stages are declared in `ontseq_platform/pipeline/stages.py`, which has no dependencies at
all — not even pydantic. The graph is data, so it can be inspected, tested and reasoned
about without executing anything.

| Stage | Tool | Applies to | Required | Verification |
| --- | --- | --- | --- | --- |
| `basecall` | Dorado | POD5 only | yes | **unverified adapter** |
| `align` | minimap2 + samtools | POD5, unaligned BAM | yes | verified with real tool |
| `intake` | samtools | all | yes | verified with real tool |
| `qc` | cramino | all | yes | verified with real tool |
| `target_coverage` | mosdepth | all | no | verified with real tool |
| `cnv` | — | all | no | not implemented |
| `sv` | Sniffles2 + cuteSV, consensus and annotations | all | no | adapters verified with synthetic contracts; real-tool CI |
| `methylation` | modkit | all | no | **unverified adapter** |
| `assemble` | — | all | yes | pure Python |
| `report` | — | all | yes | pure Python |
| `release` | — | all | yes | pure Python |

`EXECUTION_ORDER` is topologically validated at import time, so a dependency cycle or a
stage listed before its dependency is an import error rather than a runtime surprise.

### Applicability, not detection

Which stages run follows from the manifest's declared `input.kind`, never from what
happens to exist on disk. An aligned-BAM run does not "skip" basecalling; basecalling does
not apply to it. This distinction is load-bearing: a stage that does not apply must not
appear in the report as something that failed to happen.

### Three ways a stage can be out of scope

Applicability is not one question but three, and the run report keeps them apart because a
reader tracing an absent result needs to know which one they are looking at.

| Gate | Decided by | Records | Example |
| --- | --- | --- | --- |
| Input kind | `StageSpec.applicable_for` | absent from the plan entirely | `basecall` on an aligned-BAM run |
| Assay | the manifest's `assay.mode` | `applicable: false` | `target_coverage` on an lcWGS run |
| Requested analysis | the manifest's `analysis.modules` | `requested: false` | `sv` and `methylation` on a CNV-only run |

The third gate is the manifest acting as the run's scope contract. A stage that runs
anyway produces evidence nobody asked for and — the failure that motivated the gate — can
kill a run over a tool the operator had no reason to configure: a manifest declaring
`modules: [qc, cnv, report]` used to drive a structural-variant attempt regardless, so a
CNV-only run died on a missing cuteSV reference FASTA.

Skipping is never silent. Each of the three records a reason saying it is a scope
statement rather than a negative result, and a stage that *was* requested but cannot be
configured still fails closed: asking for SV evidence with no caller policy is an error,
not a skip.

### Bridging skipped stages

`intake` depends on `align`, and `align` depends on `basecall`. For an aligned-BAM run
neither applies, so `effective_dependencies(INTAKE, ALIGNED_BAM)` resolves transitively to
`()` rather than waiting forever on a stage that will never run.

### Failure propagates as NOT_RUN, never FAILED

If `qc` fails, `sv` is recorded as `NOT_RUN` with the reason naming `qc` as the cause.
Marking it `FAILED` would assert that structural-variant calling was attempted and did not
work — a claim nobody made. The distinction between the four outcomes is deliberate:

- `COMPLETED` — the stage ran and produced its artifacts.
- `NO_CALL` — the stage ran, looked, and declined to assert anything. **Not** a biological
  negative.
- `FAILED` — the stage ran and broke.
- `NOT_RUN` — the stage never started, because it does not apply, is not wired in, or a
  dependency did not deliver.

A run's verdict is `PASS` only when every *required* applicable stage reached `COMPLETED`
or `NO_CALL`. Optional stages that did not run are listed explicitly in the verdict text,
with the reminder that a stage which did not run is not a negative finding.

---

## 2. The run envelope

Every run gets a directory under `<output-dir>/<run-id>/<sample-id>/`. The file-level
layout is in [`ARCHITECTURE.md`](ARCHITECTURE.md); the directories are:

```
manifest/     sample manifest, reference lock, intake report
qc/           normalized QC
evidence/     cnv/, sv/, fusion/
alignment/    BAM and index          (never exportable)
normalized/   the validated result contract
reports/      HTML and Excel for a human reviewer
provenance/   run report, per-adapter provenance
release/      release.json + checksums.sha256
work/         scratch                (never exportable)
```

Three properties are enforced by `RunEnvelope` rather than by convention:

**Paths recorded are always envelope-relative.** `RunEnvelope.relative()` raises rather
than record a path outside the envelope. An absolute path would leak the source BAM
location and the local directory structure into reviewer artifacts, which
`docs/DATA_SECURITY.md` forbids. `RunEnvelope.path()` refuses traversal in the other
direction.

**Every write is atomic.** Content goes to a temporary name in the destination directory
and is `os.replace`d into position. An interrupted run leaves either the previous artifact
or none — never a truncated one that a later resume would accept as complete.

**Export is fail-closed twice over.** An artifact is exportable only if it is neither in an
intermediate directory (`work/`, `alignment/`) nor carries a raw genomic suffix (`.bam`,
`.vcf`, `.fastq`, `.gz`, …). The suffix list mirrors the one in
`scripts/check_repository_safety.py`: the rule that keeps raw data out of Git is the same
rule that keeps it out of a release bundle.

---

## 3. Resume

A stage is skipped on a repeat run only when **all three** conditions hold:

1. It **ran to a conclusion** last time — `COMPLETED` or `NO_CALL`. Both recorded a
   signature and both may carry artifacts, so both stand. A caller that looked and declined
   has concluded as surely as one that called something, and re-running Sniffles2 over a
   multi-hour BAM to reach the same `NO_CALL` buys nothing. `FAILED` and `NOT_RUN` always
   re-run.
2. Its **signature** is unchanged. The signature hashes the stage name, the checksums of
   every upstream artifact, the stage's own parameters, the resolved tool versions, and the
   fingerprints of any input from outside the envelope.
3. **Every artifact it claimed still verifies** — present, same size, same SHA-256.

Anything else re-runs. A timestamp comparison would silently accept an artifact produced
under different parameters or a different tool version; content addressing cannot. Resume
is an optimisation, and it must never become the reason two incompatible results share one
envelope. `--force` skips the check entirely.

CI proves this rather than asserting it: both lanes are run twice over the same inputs and
the second run's output is rejected if any stage reports anything other than `resumed`.

### Resuming a stage that changes what the next stage reads

Alignment does not only produce files; it changes what "the input" means. Everything after
it reads the aligned BAM, not the unaligned one, and the manifest is re-pointed accordingly.
A resumed stage produced its artifacts just as surely as one that just ran, so that
re-pointing has to happen on both paths.

Stages therefore declare an optional `settle` hook alongside `plan` and `execute`. The
runner calls it after a stage completes *and* after it resumes, passing the artifacts the
stage recorded — so adopting a multi-gigabyte BAM costs no second checksum pass. Doing this
inside `execute` instead is a quiet trap: the first run works, and the resumed run either
re-does everything downstream against the wrong file or fails at the intake gate.

---

## 4. What CI actually executes

Two lanes run against real binaries in the `local-real-tool-smoke` job, on synthetic data
that contains no genomic material of any kind.

**Aligned-BAM lane.** A synthetic BAM is built with samtools, then taken through
`intake → qc → sv → assemble → report → release` with real samtools, cramino, Sniffles2 and
cuteSV. CI requires both caller JSON artifacts and the consolidated consensus artifact.
The release bundle's `checksums.sha256` is then verified independently with `sha256sum -c`,
so the checksums are confirmed by a tool that shares no code with the one that wrote them.

**Unaligned-BAM lane.** `ontseq align-fixture` writes a deterministic synthetic reference
and an unaligned BAM whose reads are carved out of that reference, then the run starts one
stage earlier at `align`. Because the reads are genuine substrings of the reference,
alignment has a correct answer to find, and CI asserts on the result:

- all 24 primary reads map (a mapping-rate assertion, not just an exit code);
- at least the four reverse-complemented reads land on the reverse strand;
- the `@RG` header line and the per-read `RG:Z` tag both survive the FASTQ round trip;
- `MM`/`ML` modified-base tags and the `MD` tag are present after alignment;
- `MM` is still present on reverse-strand records specifically.

Structural-variant detection is deliberately *not* asserted in this lane. Some fixture
reads carry a 200 bp deletion so the aligner has a real gap to place, but whether either caller
calls it is the aligned-BAM lane's assertion. On the current fixture the SV stage can record
`NO_CALL`, which is a legitimate outcome and not a biological negative.

That job is what earns `align` its `verified_with_real_tool` status. Before it existed the
stage was marked `unverified_adapter`, and it would have gone back to that if the job were
removed.

### Read groups through alignment

Aligning through FASTQ is lossy by construction: `samtools fastq` discards the header and
minimap2 writes a fresh one. Read-group provenance is therefore carried in two halves — the
per-read `RG:Z` tag rides the FASTQ comment (`samtools fastq -T RG` plus minimap2 `-y`), and
the `@RG` header lines are re-attached with `samtools reheader` after sorting.

minimap2's own `-R` is deliberately unused. It accepts a single read-group line and stamps
it on every record, which would silently merge distinct read groups into one and contradict
the per-read tags. Re-attaching the source header preserves however many read groups the
input actually declared.

---

## 5. What is not verified

This is the part to read before trusting anything above.

**Basecalling has never been executed.** The Dorado adapter is written to the repository's
adapter boundary — explicit argument vector, no shell, fail-closed version and model gate —
but no Dorado binary exists in this repository's CI or development environment. It is marked
`unverified_adapter` in machine-readable form, `BasecallPolicy` refuses `status: validated`,
and any run that completes a basecalling stage carries an explicit warning in its run report
and release bundle. Treat POD5 runs as untested until someone executes one against a real
GPU and a real model.

**The methylation lane has never met modkit.** There is now a lane that reads the `MM`/`ML`
tags CI proves survive alignment (`docs/METHYLATION_LANE.md`), but no modkit binary exists in
this repository's CI or development environment, so the adapter is marked `unverified_adapter`
and a run completing that stage says so. Its bedMethyl parsing, region aggregation and refusals
are unit tested against synthetic pileups; its behaviour on real modkit output is an assumption.
What CI still does not prove is that a caller interprets modified-base tags on reverse-strand
records correctly — that needs a real run, not a lane.

**No stage output has clinical meaning.** Tool versions are pinned for reproducibility.
Thresholds are technical defaults. `qc` gates are `null` pending analytical validation. A
`PASS` verdict means the software executed as designed; it is not evidence about a sample.

**CNV is not wired in.** It is declared in the graph and records `NOT_RUN` with the reason
"No adapter is wired in for this stage." The CNV benchmarking subsystem
(`docs/CNV_BENCHMARKING.md`) exists to choose a caller on evidence before one is wired in.

**Target coverage is wired in, and only an adaptive-sampling run measures anything.** The
Mosdepth adapter runs in the canonical runner. For any other assay mode the stage records
`NOT_RUN` with the reason that per-target coverage does not apply — a scope statement, not a
coverage finding. An adaptive-sampling run without a target-coverage policy or a readable
target BED fails closed rather than producing a report that looks complete. The stage is
declared optional in the graph because an lcWGS run legitimately skips it; for an
adaptive-sampling run it is not optional, and a `FAILED` target-coverage stage fails the run.

**The release bundle is unsigned.** `signature_status` is the literal `"unsigned"`. It is a
checksum manifest, not a chain of custody.

---

## 6. One run at a time per envelope

A run holds an exclusive lock on its envelope for its whole duration: `.ontseq-run.lock` at
the envelope root, taken with `O_CREAT | O_EXCL` so two processes cannot both believe they
have it.

This is not defensive decoration. Atomic writes stop an artifact from being *truncated*, and
content-addressed resume stops a stale artifact from being *accepted* — neither notices a
second process rewriting the same run report from a different set of stage records. The
loser of that race simply vanishes from the history with nothing recording that it existed.
Under a person typing commands the race is theoretical; under a watch folder that
double-fires, it is not.

The lock records the host, PID, run and sample that took it, so a blocked run can say *who*
holds it rather than only *that* it is held:

```
ERROR: this run envelope is already in use: pid 8123 on seq-node-2,
run RUN_2026_014/AML_0031, acquired 2026-08-17T09:04:11+00:00.
If that process has stopped, remove …/.ontseq-run.lock
```

`ontseq run` exits **4** in that case, distinct from the failure exit code, so a scheduler
can tell "someone else already has this sample" apart from "this run failed" and move on
rather than retry.

Three rules govern reclaiming:

- **A crashed local run does not block forever.** If the holder is on this machine and that
  PID is gone, the lock is reclaimed — a power cut must not require manual cleanup before a
  run can be resumed. The reclaim is recorded as a warning **in the run report**, because
  stepping over another run's lock is not something that should only scroll past in a
  terminal.
- **A lock from another host is never reclaimed.** On shared storage, a crashed remote run
  and a running one are indistinguishable from here, and guessing wrong puts two live runs
  in one envelope. It fails closed and names the file to remove.
- **PID reuse fails closed.** A recycled PID makes the lock look held and the run refuses.
  Refusing a run that could have proceeded costs a delay; proceeding on one that should have
  refused costs the envelope.

---

## 7. The watch folder

```
ontseq watch /drop --manifest-template assay.manifest.yaml \
  --reference-lock GRCh38.lock.json --reference-fasta GRCh38.fa \
  --input-kind aligned_bam \
  --output-dir results/runs
```

Each sub-directory of the drop folder is one sample. The watcher decides two things per
directory, and both are stated rather than assumed.

**Is it finished being written?** A sequencer writes for hours, and analysing a run
mid-write yields a truncated result that looks complete. `--ready-marker` is a **glob**
matched anywhere beneath the sample directory; when set it is the *only* thing consulted,
because an explicit signal beats inferring completion from timestamps. Without one,
`--quiet-seconds` requires the directory to have been unmodified for that long — and the
reason recorded says in those words that quiescence is a heuristic. A file dated in the
future keeps a directory not-ready, so clock skew cannot make an in-progress run look long
finished.

A glob rather than a fixed name because that is what the instrument actually emits. MinKNOW
writes `final_summary_<flowcell>_<run>_<hash>.txt` when a GridION run completes: the name
carries per-run identifiers, and it lands in the run directory rather than the sample
directory above it. A literal top-level name would match neither, leaving the one
authoritative completion signal unusable and forcing every real run onto the heuristic.

### Running against a GridION

MinKNOW writes `<experiment>/<sample>/<run>/…`, so point the watcher at the **experiment**
directory — its sub-directories are then the samples, which is the layout above.

```
ontseq watch /data/AML_RUN_2026_08 \
  --manifest-template assay.manifest.yaml --reference-lock GRCh38.lock.json \
  --input-kind pod5 --pod5-subdir pod5_pass \
  --ready-marker 'final_summary_*.txt' --output-dir results/runs
```

`--pod5-subdir` exists because MinKNOW splits a run into `pod5_pass` and `pod5_fail` by
qscore, and **which of those enters the analysis is a scientific decision, not a filesystem
detail**: including failed reads changes the depth distribution that depth-based
copy-number methods assume. An undeclared split is refused with both directories named,
rather than resolved by picking one. A run written without qscore splitting is unambiguous
and needs no declaration.

The POD5 lane still depends on the unverified Dorado adapter — see section 5. The paths
above resolve correctly and are unit-tested; what happens after Dorado is invoked is not.

**Has it already been handled?** A ledger beside the *output* — never inside the drop
folder, which may be owned by the instrument or mounted read-only — records every attempt.
Completed samples are not repeated. A sample blocked by the run lock always retries, because
that says something about another process rather than about this sample. A **failed sample
is not retried automatically**: a deterministic failure does not become a success by being
repeated every minute, and the noise would bury the one sample somebody needs to look at.
`--retry-failed` re-attempts once the cause is understood.

### What the watcher will not infer

Sample identity, reference, genome build and assay mode are facts about a patient sample,
and a filename is not evidence of any of them. They come from `--manifest-template`, written
once per assay. Only two things are derived per sample:

- the **sample identifier**, from the directory name, and only when that name already
  satisfies the manifest contract. No cleaning, no truncation, no substitution — a repaired
  identifier ends up on a reviewer artifact under a name nobody chose, so an unusable name
  is rejected with an explanation instead.
- the **input path**, by looking for the declared `--input-kind` inside the directory.

Input kind is declared rather than sniffed, because sniffing means rules like "a BAM without
an index must be unaligned" — wrong the first time an index has not finished copying, and
wrong silently: the run would strip and re-align an already-aligned BAM and produce a
plausible result nobody asked for. Two BAMs where one is expected is likewise a refusal, not
a choice between them.

### Failure separation

A broken manifest template is one configuration mistake, not a failure per sample. Policies
and the template are resolved **once, before any directory is touched**; a problem there
exits 5 and attempts nothing. Exit codes: **2** a sample run failed, **4** the run envelope
was locked, **5** the configuration is unusable.

`SIGINT` and `SIGTERM` stop the watcher *after the current sample* rather than mid-run,
which costs one sample's runtime and avoids leaving a lock and a half-built envelope behind.

---

## 8. Asking what happened

```
ontseq status results/runs            # one line per envelope, plus the watch ledger
ontseq status results/runs --verbose  # every stage of every run
ontseq status results/runs --json     # for a monitoring check
```

A run leaves a complete record in `provenance/run.json` and a watcher leaves one in its
ledger. Both are machine-readable and neither is readable at a glance — which becomes a
problem the moment the pipeline runs unattended, because the person who needs to know
whether last night went well should not have to write a script to find out.

Two of the reported states cannot be seen from a directory listing at all:

| State | Meaning |
| --- | --- |
| `passed` / `failed` | a run finished and recorded its verdict |
| `running` | a lock is held by a process that still exists — or by an unreachable host |
| `interrupted` | a lock is held by a process that is **gone**: a run died here |
| `unfinished` | neither lock nor report; a run started and never reached a verdict |
| `unreadable` | a report exists but cannot be parsed |

`interrupted` is worth knowing before somebody deletes the directory to "start clean": the
next attempt reclaims the lock and resumes from where it stopped. For a lock taken on
another host, liveness is genuinely unknowable from here, and that is reported as unknown
(`null` in JSON) rather than collapsed into "probably fine".

Exit codes suit a monitoring check: **0** nothing wrong, **2** a run failed or its report is
unreadable, **6** a run was interrupted or never reached a verdict. A run merely *in
progress* is deliberately not a problem — a check that fires while the pipeline is doing its
job teaches people to ignore it.

A completed run that rests on an unverified adapter is called out in the summary rather than
folded into `passed`, so basecalling's status stays visible without opening the report.

The **review state** is shown alongside, so an operator scanning many envelopes sees run
health and sign-off in one pass rather than running a second command per envelope. It never
influences the exit code: `ontseq status` answers *did the runs work*, `ontseq review status`
answers *may this leave the system*. Folding the second into the first would make the check
fire on every fresh run, which is how a monitoring signal becomes noise.

---

## 9. Asking whether a run can succeed, before starting it

```
ontseq preflight <manifest> --reference-lock <lock> --run-id RUN_001
ontseq preflight ... --require-free-gb 400     # judge free space against a real figure
ontseq preflight ... --json                    # for a scheduler deciding whether to submit
```

Every gate this applies already exists inside the pipeline, and each one fails closed
correctly. The problem is *when* they fire. A POD5 run discovers a missing Dorado model
after the envelope exists and the lock is taken; an aligning run discovers a missing
reference index after intake; a run into a busy envelope discovers that only when it tries
to lock. Individually correct, collectively expensive — the feedback arrives after the run
has been queued, scheduled and partly executed.

Preflight asks the same questions up front, in a couple of seconds, **with no side effects
at all**: it creates no envelope, takes no lock, writes no artifact, and does not even
create the output directory — writability is probed on the nearest existing ancestor.

| Status | Meaning |
| --- | --- |
| `ok` | checked, and the precondition holds |
| `FAIL` | checked, and the run cannot succeed. The only status that blocks |
| `warn` | the run can proceed, and somebody should know anyway |
| `????` | genuinely not determinable from here |
| `--` | does not apply to this input kind |

What it checks: the declared input exists and has the shape its kind promises; the manifest
and the reference lock agree; the reference FASTA's `.fai` still hashes to the
`source_fai_sha256` the lock recorded; every binary the *planned* stages will invoke is
present, runnable and at its locked version; the Dorado model matches its lock and a
modified-base model was requested; for an adaptive-sampling run, that a target-coverage
policy was supplied and that the declared target BED parses into usable regions; the
envelope is free; the output location is writable.

It also reports two things about the run's *scope*, kept deliberately apart. A stage on an
`unverified_adapter` **will run**, on code nobody has executed against the real tool, and its
output is an assumption — that is `adapters.verification`. A `not_implemented` stage has no
adapter at all and will record `NOT_RUN` — that is `stages.not_implemented`, and it is not a
negative biological finding. Reporting the second as "an adapter that has never been
executed" would tell an operator that CNV rests on unexecuted code, when in fact no CNV
caller is wired in.

Three properties matter more than the list itself.

**Preflight must agree with the run.** A preflight that clears a run which then fails on the
very thing it checked converts a fast, honest failure into a slow, surprising one. So
version strings are parsed by the adapters' own parsers rather than re-implemented, and a
version lock is enforced only when a *planned* stage would enforce it — the alignment policy
locks a samtools version, but an aligned-BAM run never aligns and never applies that lock.

**A tool's absence is as fatal as its stage is required — for this run.** `sniffles` serves
only the optional SV stage, so a machine without it gets a warning saying SV will record
`NOT_RUN`, not a refusal. That is derived from `StageSpec.required`, not maintained by hand.

`required` is a property of the graph, though, and one stage needs more than the graph can
say: target coverage is optional because an lcWGS run records it as out of scope, but an
adaptive-sampling run neither skips it nor survives it failing. So preflight escalates a tool
to blocking when any stage needing it is one *this particular run* cannot do without, and a
missing Mosdepth is therefore a warning for lcWGS and a refusal for adaptive sampling. The
Mosdepth version lock is applied on the same condition, so an lcWGS run is never refused over
a tool it will not invoke.

**Not knowing is a distinct answer.** Free disk space is *reported*, not judged, unless the
caller states a requirement with `--require-free-gb`. There is no measured relationship in
this repository between an input's size and the space a run consumes — that depends on the
lab's chemistry, depth and retention policy, and nobody has measured it here. A multiplier
invented in the code would look like a validated figure and would not be one.

Exit codes: **0** nothing blocks, **2** at least one precondition makes the run impossible. A
warning or an unanswerable question never blocks, or the command would be unusable on
exactly the machines it exists to help.

CI runs preflight on both lanes before the runs it precedes, and asserts two negatives: a
stated space requirement that cannot be met, and the alignment lane pointed at the *other*
lane's reference lock — the mistake that otherwise produces confidently wrong coordinates.

### Producing the model checksum preflight compares against

```
ontseq model-lock /path/to/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
ontseq model-lock ... --list-files    # every file, in the order it enters the checksum
ontseq model-lock ... --json          # for a setup script
```

`basecall.model` may be a name or a path. A name resolves through Dorado's own cache and
cannot be checksummed, so `model_sha256` only means anything for a downloaded directory —
and until this command existed, nothing could produce that value except a preflight run
against a manifest that did not exist yet. The lock stayed `null`, and a check that is
never populated is a check that never fires.

The command prints `model` and `model_sha256` as the pair to record, because the checksum
belongs to *that directory*: a re-downloaded or moved model needs a fresh run, not the old
value copied across.

It also prints the file count and the total size, and this is the part that matters. Every
directory produces a valid-looking 64-character digest — including one holding three
zero-byte files from an interrupted download. Where the directory holds no files, empty
files, or broken symbolic links, the command exits **2** and *withholds* the checksum
rather than printing it beside a warning; a value on screen next to a caveat is a value
somebody pastes. Nothing beyond those defects is judged: this repository has never had a
real Dorado model, so it does not know what a correct one contains, and a file-count
threshold invented here would be a guess wearing the clothes of a check.

The digest is the one `basecall.model_signature` already computed and preflight already
compares against — one implementation, called by both, so the command that produces the
value and the check that enforces it cannot drift apart. CI asserts that equality directly.

---

## 10. Signing off a run

```
ontseq review record <envelope> --decision accepted --reviewer dr.mueller --note "…"
ontseq review status <envelope> --verbose
ontseq review status <envelope> --require-reviewers 2      # four-eyes release gate
ontseq review status <envelope> --json                     # for a release script
```

A run produces evidence. Someone has to look at it and say so, and that statement has to
survive as a record — otherwise "this was reviewed" is a claim nobody can check later.

| State | Meaning |
| --- | --- |
| `pending` | no review recorded |
| `accepted` | the latest review accepts **exactly the content on disk** |
| `rejected` | the latest review rejects it |
| `stale` | reviews exist, but the release bundle changed after the last of them |
| `broken` | the trail does not verify: an entry was removed, reordered or edited |
| `unreadable` | the log cannot be parsed |

Exit codes match `ontseq status`: **0** nothing in the way, **2** rejected or the trail does
not verify, **6** not reviewed yet or reviewed against different content.

### Bound to content, not to a directory

Every entry carries the SHA-256 of `release/release.json`, which already covers the run
report and every exportable artifact by checksum. A review that pointed at a *path* would
keep vouching for whatever later appeared there. Binding to content means a changed run
makes the review `stale`: the judgement still stands for what it saw and says nothing about
what is there now.

### A reviewed envelope cannot be re-run

`ontseq run` refuses with exit **7** when the envelope's latest review accepts its current
content. This is the property that makes the rest mean anything — the lock stops two runs
colliding *now*, content-addressed resume stops a stale artifact being *accepted*, and
neither notices that a human signed this off yesterday and a resumed run is about to rewrite
what they signed.

It is deliberately **not overridable by a flag**. A flag would be used, and the correct
alternative costs nothing: use a new `--run-id`. The reviewed envelope then keeps its review
and the new run gets its own. A rejected or stale review does not block — a rejection is
often precisely why somebody re-runs.

### Append-only, and the shape shows it

Each entry names the digest of the entry before it, so removing, reordering or editing one
breaks the chain there and at every later point. Nothing is overwritten: a reviewer who
accepted and then rejected leaves both facts behind, which is the entire purpose of an audit
trail. Appending to a log that no longer verifies is refused, because a record that looks
continuous and is not would be worse than one that admits it is broken.

### What this is **not**

Two limits are printed on every human-readable report and carried as explicit `false` fields
in the JSON:

- **`identity_is_authenticated: false`.** The reviewer name comes from the command line.
  Nothing authenticated it, so the record says `asserted` rather than implying otherwise.
- **`chain_is_tamper_proof: false`.** There is no key. Anyone who can rewrite the file can
  recompute the whole chain. It detects accidental corruption and casual editing; that is
  all it detects.

`signature_status` on the release bundle therefore still reads `unsigned`. A qualified
electronic signature needs an authorised key, an identity provider and a records policy,
none of which exist here — and a record that quietly implied one would be worse than no
record at all. See ADR-021.

CI signs off a real run, asserts the trail binds to the bundle digest, asserts a four-eyes
gate is not satisfied by one person, asserts that re-running the envelope exits 7, and
asserts that an edited trail is detected as `broken`.

---

## 11. Schedulers call the pipeline; they do not reimplement it

`workflow/aligned_bam.smk` is a single Snakemake rule that invokes `ontseq run`. It used to be
five rules calling the per-stage commands, which made Snakemake a second execution path —
flat files instead of an envelope, no per-artifact tool versions, no release bundle, and
mtime-based resume. Only the runner's path was ever proven in CI, so the second one was an
unproven way to produce results that looked the same.

Snakemake keeps what it is good at: cluster submission, resource declarations, fanning many
samples out through an executor plugin. Resume stays with the runner, because Snakemake
compares timestamps and the runner compares content, parameters and tool versions.

The same rule applies to anything added later. A watch folder, a LIMS trigger and a REST
endpoint are all callers of `ontseq run`. The moment one of them re-derives the stage order
itself, there are two behaviours again and only one of them is tested.

---

## 12. Extending the pipeline

To add a stage:

1. Add it to `StageId` and `STAGE_SPECS` with its dependencies, applicability, purpose and
   an honest `VerificationStatus`. Import will fail if the order is inconsistent.
2. Write plan and execute functions and register them in `IMPLEMENTATIONS`. The plan
   function declares everything the resume signature must cover; if a parameter is not in
   the plan, changing it will not re-run the stage.
3. Write the artifacts through the envelope, never with a bare `open()`. That is what makes
   them atomic, checksummed and correctly classified as exportable or withheld.
4. If the stage changes what a later stage reads, put that change in a `settle` hook rather
   than in `execute`, so it also happens when the stage is resumed.
5. Leave `VerificationStatus` at `unverified_adapter` until CI executes the real tool, then
   flip it in the same commit that adds the job.

Step 5 is the point of the field. It is only worth anything if it is not flipped early.
