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
- `docs/CNV_BENCHMARKING.md` and ADR-008 through ADR-012.
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
- `docs/PIPELINE_EXECUTION.md` and ADR-013 through ADR-016.

### Changed

- `scripts/export_schemas.py` also exports the CNV and pipeline contracts.
- CI executes the whole pipeline end to end against real binaries in two lanes: from an
  aligned BAM, and from an unaligned BAM through real minimap2. It re-runs each pipeline
  to prove every stage resumes unchanged, and verifies the release bundle checksums with
  `sha256sum -c` rather than with the code that wrote them.
- `StageId.ALIGN` moved from `unverified_adapter` to `verified_with_real_tool`, in the
  same change that added the CI job executing minimap2.
- `workflow/envs/aligned_bam.yaml` pins `minimap2=2.28`, matching the version lock the
  alignment adapter enforces.

### Fixed

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
