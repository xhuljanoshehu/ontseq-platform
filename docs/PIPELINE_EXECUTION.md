# Pipeline Execution

Research use only. Nothing described here is validated for diagnostic use.

This document describes how one sample is taken from raw input to a checksummed release
bundle in a single command, and — just as important — what that command does *not* prove.

```
ontseq run <manifest.json> --reference-lock <lock.json> --run-id <RUN> [--reference-fasta …]
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
| `target_coverage` | — | all | no | not implemented |
| `cnv` | — | all | no | not implemented |
| `sv` | Sniffles2 | all | no | verified with real tool |
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
`intake → qc → sv → assemble → report → release` with real samtools, cramino and Sniffles2.
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
reads carry a 200 bp deletion so the aligner has a real gap to place, but whether Sniffles2
calls it is the aligned-BAM lane's assertion. On the current fixture the SV stage records
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

**Modified-base tags are carried, not interpreted.** CI proves `MM`/`ML` survive alignment,
including on reverse-strand records. It does not prove that a downstream methylation caller
reads them correctly, because there is no methylation lane yet to read them.

**No stage output has clinical meaning.** Tool versions are pinned for reproducibility.
Thresholds are technical defaults. `qc` gates are `null` pending analytical validation. A
`PASS` verdict means the software executed as designed; it is not evidence about a sample.

**CNV and target coverage are not wired in.** Both are declared in the graph and both record
`NOT_RUN` with the reason "No adapter is wired in for this stage." The CNV benchmarking
subsystem (`docs/CNV_BENCHMARKING.md`) exists to choose a caller on evidence before one is
wired in; target coverage is developed in the adaptive-sampling work stream and plugs into
the same seam.

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
  --reference-lock GRCh38.lock.json --input-kind aligned_bam \
  --output-dir results/runs
```

Each sub-directory of the drop folder is one sample. The watcher decides two things per
directory, and both are stated rather than assumed.

**Is it finished being written?** A sequencer writes for hours, and analysing a run
mid-write yields a truncated result that looks complete. `--ready-marker` names a file the
producer writes when it is done; when set it is the *only* thing consulted, because an
explicit signal beats inferring completion from timestamps. Without one, `--quiet-seconds`
requires the directory to have been unmodified for that long — and the reason recorded says
in those words that quiescence is a heuristic. A file dated in the future keeps a directory
not-ready, so clock skew cannot make an in-progress run look long finished.

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
modified-base model was requested; the envelope is free; the output location is writable.

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

**A tool's absence is as fatal as its stage is required.** `sniffles` serves only the
optional SV stage, so a machine without it gets a warning saying SV will record `NOT_RUN`,
not a refusal. That is derived from `StageSpec.required`, not maintained by hand.

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

---

## 10. Schedulers call the pipeline; they do not reimplement it

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

## 11. Extending the pipeline

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
