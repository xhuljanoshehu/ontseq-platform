# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

### Added

- `ontseq_platform.cnv` subsystem: a segmentation-independent CNV comparison core, an
  explicit observability mask, multi-source truth representation, deterministic
  simulation, a baseline read-depth caller and stratified aggregation.
- Base-level per-state confusion scoring over an exact genome breakpoint partition, so
  copy-number agreement no longer depends on how either side segmented.
- Closed vocabulary of exclusion reasons with per-reason base accounting, keeping
  biological negativity, no-call and technical failure distinguishable.
- `background_state` and `resolution_bp` on truth and call sets, encoding open-world
  versus closed-world semantics and each source's detection limit.
- Per-segment breakpoint uncertainty, with breakpoint metrics withheld when the truth
  source cannot support them.
- ISCN karyotype to copy-number conversion against a versioned cytoband resource, with
  unsupported constructs recorded and surfaced rather than dropped.
- Deterministic gamma-Poisson dilution and coverage simulation, plus empirical and
  model-based limit-of-detection estimation with explicit withholding.
- `ontseq-baseline-readdepth`, a non-reportable control caller.
- Header-driven segment-table adapters with declarative column mappings for the generic
  IGV `SEG` format and ichorCNA.
- Paired method comparison via McNemar's exact test on truth events assessable under
  both methods, with the p-value withheld when no discordant pair exists.
- A genome partition that reconciles exactly, enforced in the core and re-validated by
  the contract, so every metric carries an auditable denominator.
- CLI: `cnv-evaluate`, `cnv-aggregate`, `cnv-compare-methods`, `cnv-karyotype-truth`,
  `cnv-demo-benchmark`.
- JSON Schemas for the CNV truth set, call set, benchmark case, evaluation report,
  aggregate report and cytoband table.
- `docs/CNV_BENCHMARKING.md` and ADR-008 through ADR-012, plus ADR-018 through ADR-020.
- Knowledge-base annotations now reach the Excel workbook, not only the HTML report: sheet
  `11_Annotations` carries every matched record under a banner in row 1 stating that these
  are classifications of database records rather than findings about the sample. The banner
  is *in* the sheet because a spreadsheet has no prose — a reviewer sorts by the assertion
  column and reads "Pathogenic" beside an event identifier, and nothing else on the grid
  would say what that means. Rows whose record origin does not match the assay's question
  are filled in the warning colour, because colour is read before column ten is.
- `ontseq model-lock`, which produces the checksum `BasecallPolicy.model_sha256` locks a
  Dorado model to. The lock existed and preflight enforced it, but nothing could generate
  the value except a full preflight run against a manifest that did not exist yet — so in
  practice it stayed `null` and enforced nothing. The command also reports the file count,
  the total size and the listed defects, because every directory yields a valid-looking
  64-character digest, including one holding three zero-byte files from an interrupted
  download; it exits 2 and withholds the checksum rather than printing it beside a warning.
  The digest is the one `basecall.model_signature` already computed, now the single
  implementation both call, so the command and the check it feeds cannot drift apart.
- `windows/ONTSeq.bat`: one file on the desktop that starts the service in WSL and opens
  the browser. It checks whether the service is already running — a second double-click
  would otherwise fail on the busy port and look like a fault when nothing is wrong — and
  says in plain words what to do when WSL is missing, instead of showing an error nobody
  can act on. Four lines at the top are edited once per site and never again.
- Dragging a BAM onto the page selects it. A browser hands JavaScript only the file name,
  never the path, so the service searches the allowed roots for that name and reads the
  file where it already lies. Nothing is copied: at thirty gigabytes an upload would put
  the data on the same disk it started on, hours later. Where a name occurs more than once
  every match is offered — taking one silently would analyse a sample nobody chose — and a
  search that ran out of its directory budget says so rather than reporting "not found".
- The interface has two views, because it has two audiences. The operator picks a BAM and
  starts a run; the reviewing physician reads finished cases and signs them off. Both live
  in the same page behind a tab, and the second one lists **every** envelope — including the
  runs that failed and the sign-offs that went stale because the release changed underneath
  them. A list showing only what succeeded looks tidier and is less true.
- Sign-off is available in the interface after all. An earlier note said it did not belong
  there because a token-authenticated single-user page authenticates nobody — but neither
  does the command line, and the trail already records the name *as asserted*. The
  interface is no weaker than the path it mirrors, and it says so where the name is typed.
- `ontseq serve`, a local browser interface: one page, served by the pipeline itself, from
  which an operator picks a BAM, starts a run and watches it. It computes nothing — it calls
  the same `run_pipeline` the command line calls and reads progress from the same
  `provenance/run.json` that `ontseq status` reads.
- `service/guard.py` holds every decision where a wrong answer is dangerous — constant-time
  token comparison, the allowed-root boundary with symlinks resolved before the comparison,
  the `Host`/`Origin` checks against DNS rebinding and other open tabs, and Windows/WSL path
  translation. It carries no dependency beyond the standard library, deliberately: this is
  the part that must not be first executed on a CI runner, and 29 tests cover it locally.
- The page is served rather than opened from disk. A `file://` page has an opaque origin,
  cannot be told apart from any other local page, and would need the service to allow
  cross-origin requests at all. Serving it means the token can be handed over without
  anyone copying it, and a stray copy of the file has no access — the intended outcome.
- The file picker is server-side because a browser hands JavaScript only the file name,
  never the path. That also lets the listing say whether a `.bai` exists, instead of the
  run failing on it hours later.
- Progress carries all four stage outcomes. Collapsing `NOT_RUN` and `NO_CALL` into one
  empty circle would undo, in the one place everybody looks, the distinction the rest of
  the system exists to preserve — and CNV is `NOT_RUN` today.
- CI starts the service, obtains the token the way a browser does, and asserts the four
  refusals it must make: no token, wrong token, foreign `Origin`, forged `Host`, and a path
  outside the allowed roots. Then it runs an analysis through the HTTP interface end to end.
- `docs/SCHNELLSTART.md`, a step-by-step operating guide in German, and
  `examples/manifests/gridion_adaptive_sampling.example.yaml`, a template for the case this
  platform is actually aimed at: a GridION or PromethION run with adaptive sampling where
  the instrument has already basecalled, so the chain enters at `aligned_bam` and the
  untested POD5 lane is not involved. The guide is in German because it is an operating
  manual for the laboratory; the rest of `docs/` stays English as the technical record. It
  states up front which stages do *not* run, so nobody looks for a CNV result that the
  stage graph correctly reports as `NOT_RUN`.
- CI validates every shipped manifest template, so a template cannot go stale against the
  contract without the job failing.
- `containers/Dockerfile`, for an analysis machine with no network at all: one image holding
  the pipeline and the four pinned tools, built once where there is a network and carried
  across as a single file. It installs from the same environment file the CI smoke job uses,
  so image and CI cannot end up pinned differently, and `pip install --no-deps` keeps pip
  from resolving a second set of versions on top of the ones conda already pinned — which
  would leave the provenance record naming versions that were never executed. There is
  deliberately no separate Apptainer definition: `apptainer build` reads the same saved
  archive, and two recipes for one image are two things that drift.
- CI builds that image, asserts the four tool versions inside it, runs the pipeline in it,
  then saves, reloads and re-runs it — the way it actually travels. A build recipe that has
  never been built is a promise; this is the same rule ADR-015 applies to adapters.
- `db_records_matched` on the four event sheets. A reviewer who never opens the annotation
  sheet would otherwise not know there was anything there to open. It is a count, named for
  what it counts, so that it cannot be read as a classification of the finding itself.
- `truth_resolution_silent_bases`: a fourth term in the genome partition holding calls
  finer than the truth source can resolve, so they are neither confirmed nor counted as
  false positives. The count and the number of affected calls are stated in a warning,
  because this is the only exclusion in the design that flatters the caller.
- `observed_direction`, `underpowered`, `minimum_attainable_p_value` and `alpha` on the
  paired method comparison, separating the direction the counts lean from the direction a
  test supports.
- `SpecimenClustering` and a specimen-weighted detection rate on the aggregate report, plus
  `discordant_specimens` on the paired comparison, so event-level intervals are never shown
  without the clustering that qualifies them.
- `reports_biological_negative` on a call set, so a method that looked everywhere and found
  nothing is distinguishable from one that could not look.
- `ontseq_platform.pipeline` subsystem: a dependency-free stage graph, a self-describing
  run envelope with atomic writes and content-addressed resume, and a runner that records
  what ran, under which tool versions, producing which checksummed artifacts.
- `ontseq run`: one command takes a single sample from its declared input kind through to
  a checksummed release bundle, and resumes an interrupted run without re-executing work
  whose inputs are unchanged.
- Minimap2 alignment adapter, including read-group preservation across the FASTQ round
  trip (`samtools fastq -T RG` plus `samtools reheader`) and modified-base tag carry-over.
- Dorado basecalling adapter, shipped explicitly and machine-readably as an **unverified
  adapter**: it has never been executed against a real Dorado binary.
- `VerificationStatus` per stage, surfaced in the run report and the release bundle, so a
  stage that completed on an adapter nobody has executed says so in its own output.
- `ontseq align-fixture`: a deterministic synthetic reference plus an unaligned BAM whose
  reads are carved out of it, so the alignment lane can be exercised with real tools.
- `SubprocessRunner.run_to_file` for tools that emit binary output on stdout.
- JSON Schemas for the run report, release bundle, alignment policy and basecall policy.
- An exclusive per-envelope run lock, so two `ontseq run` invocations cannot work on the
  same run. A stale lock left by a crashed local run is reclaimed and the reclaim is
  recorded in the run report; a lock held from another host is never reclaimed. `ontseq
  run` exits 4 when the envelope is in use, distinct from its failure exit code.
- `ontseq watch`: processes every ready sample directory in a drop folder, once or
  continuously. Readiness is decided by an authoritative marker file or, failing that, a
  quiescence window that is documented as a heuristic. A ledger beside the output records
  every attempt, so completed samples are not repeated and failed ones are not retried
  until `--retry-failed` says the cause is understood. Sample identity and assay metadata
  come from a manifest template; only the identifier and the input path are derived, and an
  unusable directory name is rejected rather than repaired. Input kind is declared, not
  sniffed. SIGINT/SIGTERM stop after the current sample.
- `ontseq watch` handles a real GridION layout: `--ready-marker` is a glob matched anywhere
  beneath the sample directory, so MinKNOW's `final_summary_<flowcell>_<run>_<hash>.txt`
  completion signal is usable, and `--pod5-subdir` declares which of `pod5_pass` /
  `pod5_fail` to basecall. An undeclared split is refused with both directories named:
  including failed reads changes the depth distribution depth-based copy-number methods
  assume, so it is a declared decision rather than a guess.
- `ontseq status`: reports the state of every run envelope beneath an output directory,
  as text or JSON, plus the watch ledger. Distinguishes a run that is *running* from one
  that was *interrupted* — a lock whose holder is gone means a run died there and the next
  attempt will resume it — and reports liveness as unknown, not false, for a lock taken on
  another host. Exit codes suit a monitoring check: 0 nothing wrong, 2 a failed or
  unreadable run, 6 an interrupted or unfinished one; a run in progress is not an alert.
- `ontseq preflight`: checks a run's preconditions before it starts, with no side effects —
  no envelope, no lock, no artifact, not even the output directory. It verifies the declared
  input, that the reference FASTA's index still hashes to what the lock recorded, that every
  binary the *planned* stages will invoke is present and at its locked version, that the
  Dorado model matches its lock, and that the envelope is free. A tool serving only optional
  stages warns rather than blocks, derived from `StageSpec.required`. Version strings are
  parsed by the adapters' own parsers, and a version lock is applied only where a planned
  stage would apply it, so preflight cannot refuse a run that `ontseq run` would complete.
  Free disk space is reported rather than judged unless `--require-free-gb` states a
  requirement: no measured relationship between input size and space consumed exists here,
  and an invented multiplier would look like a validated figure. It separates a stage that
  will run on an unexecuted adapter from a stage with no adapter at all, because "code
  nobody has run" and "nothing is wired in, so this records NOT_RUN" are different claims.
- `ontseq review`: records who signed a run off, bound to the SHA-256 of the release bundle
  they saw rather than to a directory, so a changed run makes the review *stale* instead of
  silently vouching for something else. The trail is append-only and each entry names the
  digest of the one before it, so removing, reordering or editing an entry is detectable.
  A four-eyes gate is available via `--require-reviewers`, counting distinct reviewers of
  the current content. The record states what it is not: the identity is `asserted` and
  nothing authenticated it, and the chain is tamper-evident rather than tamper-proof because
  there is no key — both carried as explicit `false` fields in the JSON output.
- `ontseq run` refuses, with its own exit code 7, to write into an envelope whose latest
  review accepts its current content. Deliberately not overridable: a resumed run would
  rewrite what somebody signed off, and using a new run id costs nothing.
- `review/` in the run envelope.
- `ontseq status` shows each envelope's review state beside its run state, so one pass over
  an output directory answers both questions. It stays out of the exit code: an unreviewed
  run is not a fault, and a check that fires on every fresh run teaches people to ignore it.
- `ontseq_platform.knowledge`: a ClinVar annotation layer that attaches records to
  copy-number and structural findings **without classifying them**. Each annotation carries
  the assertion verbatim together with the vocabulary it belongs to (`acmg_germline`), the
  record's own origin, the match type and reciprocal overlap, NCBI's review status and star
  rating, and the checksum of the exact weekly release it came from. Records whose origin
  does not match the assay's question are kept and marked as secondary findings rather than
  filtered, because filtering them would be a clinical decision made invisibly by code.
  Nothing in the package sets `reportable` or raises `confidence`: deciding that needs
  somatic criteria (ELN, ICC, a local gene list) this repository does not have.
- `ontseq annotate`: applies a locked ClinVar release to a result contract and writes the
  annotations back onto the events, leaving every judgement field untouched. CI proves that
  by diffing `reportable` and `confidence` across the annotation step.
- `AnalysisSpec.intent` (`somatic` / `germline` / `both`), **without a default**. It decides
  how every knowledge-base assertion is read, and guessing would settle that silently. Left
  unset, scope alignment is reported as unknown rather than assumed.
- `KnowledgeResourceLock` and `EventAnnotation` contracts, plus a `GenomicEvent.annotations`
  list. `EventAnnotation` refuses an empty `caveats` list: without them a database
  classification reads as a finding about the sample, which is the one reading the whole
  design exists to prevent.
- Knowledge-base annotations are rendered in the HTML report under a heading stating that
  they are classifications of database records rather than findings about the sample.
- `examples/knowledge/synthetic.clinvar.variant_summary.txt`, an invented fixture in NCBI's
  layout containing no real record, with one row for each awkward case the loader must
  handle: another assembly, a variant with no matchable extent, and ClinVar's `-1`
  placeholder for unplaced records.
- `docs/PIPELINE_EXECUTION.md` and ADR-013 through ADR-017, plus ADR-021 and ADR-022.

### Changed

- `scripts/export_schemas.py` also exports the CNV and pipeline contracts.
- The per-adapter version parsers are public (`align.parse_version`, `qc.cramino_version`,
  `sniffles.sniffles_version`, `basecall.dorado_version`) so preflight reaches the same
  answer the run will, rather than re-implementing the parsing beside it.
- CI executes the whole pipeline end to end against real binaries in two lanes: from an
  aligned BAM, and from an unaligned BAM through real minimap2. It re-runs each pipeline
  to prove every stage resumes unchanged, and verifies the release bundle checksums with
  `sha256sum -c` rather than with the code that wrote them.
- `StageId.ALIGN` moved from `unverified_adapter` to `verified_with_real_tool`, in the
  same change that added the CI job executing minimap2.
- `workflow/envs/aligned_bam.yaml` pins `minimap2=2.28`, matching the version lock the
  alignment adapter enforces.

### Fixed

- Read-group preservation was only ever exercised with a *single* read group, so a bug
  collapsing every read onto the first group would have been invisible — the first group was
  the only group. A real ONT sample is routinely sequenced across several flowcells. The
  alignment fixture now carries two, and CI asserts both survive the FASTQ round trip *and*
  that the per-group read counts are unchanged, which a presence check alone would not catch.
  The assertion passed on its first run against real `minimap2` and `samtools`: the collapse
  this was written to detect does not occur. That is now a measured property of the current
  pinned versions rather than an assumption.
- The unit test asserting a *single* read group was left behind by that change and failed in
  CI. It also compared the `RG` tag by substring, and `SYNTHETIC_ALIGN_RG` is a prefix of
  `SYNTHETIC_ALIGN_RG2` — so the check would have passed on a record belonging to the other
  group. The tag is now read as a whole field, and which reads sit on which group is asserted
  rather than only how many there are.
- CI ran `Type check` and `Unit tests` only after every earlier step passed, so a lint error
  hid the type errors, which hid the test failures — one round trip per layer. Both steps now
  run whenever the install succeeded, and the job still reports failure. This matters more
  here than it would elsewhere: pydantic cannot be installed in the development sandbox, so
  21 of 36 test modules are first executed on the runner.

- A truth source's declared resolution affected only a warning, not the score. Calls below
  it were counted as false positives while the warning beside them said they must not be —
  and the number, not the warning, is what aggregates. Resolution now removes those bases
  from the evaluable genome, and only where the truth asserts its background: resolution
  limits what a source can deny, never what it can affirm.
- `favours` named a winning method from the raw discordant counts, so a 4-0 split reported
  a winner at p=0.125 — a count at which no possible outcome could have been significant.
  A method is named only when the test is significant at a pre-specified alpha.
- A `COMPLETED` call set was required to contain segments, which forced a method that
  looked everywhere and found nothing to report `NO_CALL` — exactly the conflation between
  a biological negative and an inability to look that the vocabulary exists to prevent.
- The adaptive-sampling off-target read population was described as forming a near-uniform
  low-coverage background. That is a hypothesis about the local assay, not a measured
  property, and it is now written as one in the code and both documents.

- The alignment adapter's limitations claimed read groups were inherited from the
  unaligned BAM, but `samtools fastq` drops them and minimap2 writes a fresh header. Read
  groups are now genuinely preserved, and CI asserts it.
- The alignment adapter required its caller to have created the output directory.
  `samtools sort` does not create it, so a run whose input declared no read groups would
  have failed on a fresh envelope. The adapter now creates it itself.
- A `NO_CALL` stage re-executed on every resume. It ran to a conclusion and recorded its
  artifacts like any completed stage, so it now resumes as one; re-running a caller to
  reach the same `NO_CALL` was pure cost.

### Validation impact

No CNV method is selected, promoted or validated by this change. Every figure the new
subsystem produces is an engineering measurement against a declared truth set on
synthetic data. Comparison thresholds, the coverage floor and the exclusion tracks are
engineering parameters that must be locked before comparative results are inspected;
none of them is a validated adequacy or reportability threshold. The existing SV
benchmark contract in `benchmark.py` is unchanged.

The pipeline subsystem changes how a run is executed and recorded, not what any result
means. A `PASS` verdict states that every required stage executed as designed on the
declared input; it is not evidence about a sample. Basecalling remains unexecuted and is
labelled as such in every artifact a POD5 run produces. CNV and target coverage are
declared in the stage graph but wired to nothing, and record `NOT_RUN` rather than a
negative finding. The release bundle is a checksum manifest, not a signed chain of
custody: `signature_status` is the literal `"unsigned"`.

## 0.3.0 - 2026-08-14

### Added

- Typed Sniffles2 policy and call-report contracts with generated JSON Schemas.
- Conservative Sniffles2 v2.8.0 adapter using shell-free execution, symbolic VCF and PASS-only
  normalization with filter-reason accounting.
- Defensive DEL/DUP/INV/INS/BND normalization with counted rejection reasons and `NO_CALL`.
- Runtime-generated synthetic long-read BAM fixture and real samtools/Cramino/Sniffles2 smoke test.
- CI job with a pinned Bioconda environment and selected raw-data-free reviewer artifacts.
- Candidate SV evidence in the common JSON, HTML and Excel result path.

### Changed

- Extended evidence provenance with supporting-read strands, coverage context and mean alignment
  NM.
- Connected the aligned-BAM Snakemake DAG to the candidate-only Sniffles2 adapter.
- Bumped the research software foundation to version 0.3.0.

### Validation impact

The aligned-BAM path can now change biological candidate output by executing Sniffles2. Every
normalized event is forced to `unclassified` and `reportable: false`; BND is not treated as a
fusion, and ISCN remains `NOT_RUN`. The included thresholds are engineering defaults only. This
change requires locked HG002/HG008 and intended-use AML validation before any promotion.

## 0.2.0 - 2026-08-14

### Changed

- Reframed the architecture as an independent, evidence-led implementation.
- Reclassified the thesis as historical context rather than a technical specification.
- Replaced unvalidated CNV/SV defaults with benchmark-gated candidate sets.
- Added mandatory coverage, tumor/blast-fraction and no-call validation axes.
- Bumped the research software foundation to version 0.2.0.
- HTML and Excel reports now expose explicit per-module execution status.

### Added

- Living literature and benchmark evidence base with explicit applicability and limitations.
- Fail-closed aligned-BAM/BAI/header/reference intake gate using samtools.
- Versioned reference locks generated from FASTA index files.
- Cramino JSON adapter with normalized, path-free descriptive QC.
- Deterministic normalized-event CNV/SV benchmark engine and synthetic fixtures.
- Snakemake entry points for the aligned-BAM MVP and synthetic benchmarks.

### Validation impact

The report schema and workbook gain module-status data. Aligned-BAM intake can now stop a run on
technical incompatibility. No scientific caller has been promoted to a production or clinical
default; CNV/SV/fusion output remains `NOT_RUN` in the MVP.
