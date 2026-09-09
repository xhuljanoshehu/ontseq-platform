# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## 0.8.1 — 2026-09-09

### Fixed

- Restore optional, typed coverage arguments to HTML reporting and retain all coverage,
  SV review, fusion, ISCN disposition and reference/panel provenance sections alongside
  the new execution-state overview. Missing measurements remain unavailable, never zero.
- Reject cross-sample/build coverage exports and non-finite, contradictory or duplicated
  coverage evidence before writing either HTML or Excel. Valid numerical values and
  caller thresholds are unchanged.
- Record ordinary exceptions at adapter planning/execution/adoption boundaries as FAILED;
  interrupts still propagate. Remove previous release markers under the run lock before
  retrying, so failed attempts cannot leave a stale success bundle visible.
- Resolve Windows-only ctypes entry points inside their existing platform guards without
  masking typing failures. Align the methylation API's positive test fixture with the
  normalized tool provenance rather than weakening runtime identity checks.
- Synchronize Core/Desktop/current operator instructions at 0.8.1, regenerate the two
  stale readiness schemas and exercise pytest functions as well as unittest classes in CI.

- Replace overflowing binomial-coefficient arithmetic with a stable decreasing-tail
  recurrence and binary threshold search; reject invalid power/purity inputs. The
  binomial model and error/power cutoffs are unchanged. Fix maximum-matching recursion
  exhaustion on long ambiguous variant chains; compare all 512 three-by-three graphs
  against an independent exhaustive oracle.
- Add opt-in real modkit 0.4.1 integration in CI for 75%, measured zero, low-depth/empty
  targets and missing MM tags; this tests interoperability, not clinical accuracy.

### Validation impact

- Presentation and rejection behavior change; no calling threshold, assembly, basecalling
  model, biological classification or automatic clinical release is changed. No sensitivity
  or specificity improvement is asserted without an independent stratified benchmark.
- See `docs/CLINICAL_VALIDATION.md` for repair scope and accuracy acceptance requirements.

## Unreleased

## 0.8.0 - 2026-09-09 (local engineering candidate)

- Derive the default Desktop resource directory (`resources-v0.8.0`), runtime version and
  Desktop label from one release constant. Include current resource-path examples in
  `make versions` so future software releases cannot retain a stale default. Preserve
  explicitly configured installations and independently versioned reference, panel,
  knowledge and schema identities. This changes installation naming only, with no
  biological output or analytical validation impact.
- Reuse the interval-identity correction from the methylation handover branch
  (`b2d5256`): BED targets sharing a label retain their own chromosome, coordinates,
  call counts and missing-measurement status. Overlapping targets still independently
  include shared sites. Previously these rows could inherit pooled counts, including
  measurements from another chromosome or a measured value for an empty target.
- Record `ontseq_region_assignment=interval-identity-v1` in newly normalized
  methylation report tool parameters and stage resume signatures without rewriting
  source tool metadata. Old runs must recompute this stage even before the fix receives
  a new commit identity; unchanged corrected runs can still resume. Add a
  synthetic regression for repeated labels across disjoint, overlapping, empty and
  cross-chromosome intervals. See `docs/CLINICAL_VALIDATION.md` for validation impact.
## 0.8.0 - 2026-09-09 (local engineering candidate)

- Derive the default Desktop resource directory (`resources-v0.8.0`), runtime version and
  Desktop label from one release constant. Include current resource-path examples in
  `make versions` so future software releases cannot retain a stale default. Preserve
  explicitly configured installations and independently versioned reference, panel,
  knowledge and schema identities. This changes installation naming only, with no
  biological output or analytical validation impact.
- Reuse the interval-identity correction from the methylation handover branch
  (`b2d5256`): BED targets sharing a label retain their own chromosome, coordinates,
  call counts and missing-measurement status. Overlapping targets still independently
  include shared sites. Previously these rows could inherit pooled counts, including
  measurements from another chromosome or a measured value for an empty target.
- Record `ontseq_region_assignment=interval-identity-v1` in newly normalized
  methylation report tool parameters and stage resume signatures without rewriting
  source tool metadata. Old runs must recompute this stage even before the fix receives
  a new commit identity; unchanged corrected runs can still resume. Add a
  synthetic regression for repeated labels across disjoint, overlapping, empty and
  cross-chromosome intervals. See `docs/CLINICAL_VALIDATION.md` for validation impact.
- Advance Core, Desktop and the local workspace together to 0.8.0. Reuse the installed
  GRCh37/GRCh38 reference families without downloads or analytical policy changes.

### Fixed

- Sniffles2 and cuteSV now validate a bracketed BND ALT as exactly one of the four VCF forms
  before using its mate coordinate. Local sequence is limited to ASCII `ACGTN`, mate positions
  to positive ASCII digits, and whitespace, comma-separated alleles and trailing content are
  rejected rather than partially matched. Only explicit `<BND>` and `<TRA>` alleles may use the
  `CHR2`/`END` fallback. For a bracketed allele, `CHR2` must agree with the ALT whenever present;
  `END` is additionally compared when the complete `CHR2`/`END` mate pair is present. Safe
  `chr`-prefix/integer normalization is bounded and conversion errors become counted rejections.
  The ALT, record ID and internal bracket form remain absent from normalized evidence.

- `docs/EVIDENCE_BASE.md`: a second candidate survey (2026-08-31) covering small variants,
  methylation and clinical annotation, plus a packaging audit against the live Bioconda index.
  Seven new evidence-matrix rows, each with its applicability limit and the resulting decision.

  Three findings change the shortlist rather than extend it. **The largest gap is small
  variants, not another CNV caller**: ClairS-TO is tumour-only by design, which matches AML
  diagnostics having no matched normal, and it would make ten of the twenty-four criteria in
  `GUIDELINE_CRITERIA_DRAFT_v0` evaluable. **Methylation is evidence already being paid for
  and discarded**: the published nanopore leukaemia classifier consumes sparse genome-wide
  data, which is the shape of the off-target fraction this assay already produces, and a
  preprint characterises its failure mode at low blast fraction — the stratum that matters
  most here. **AnnotSV's ACMG/ClinGen ranking is a germline vocabulary**, joining ClinVar as a
  source whose content is useful and whose classification answers a different question than
  AML asks; ADR-022 applies unchanged.

  The survey also corrects a figure this document reasoned from: the EPI2ME row assumed "a
  roughly 3x lcWGS assay", but the measured local run produced 8.8x off-target. Several
  conclusions depend on which is right, so it is flagged for confirmation across runs rather
  than silently rewritten. And `Spectre`, referenced there as a CNV comparator, is **not** on
  Bioconda — it carries the same packaging cost as ClairS-TO rather than being the cheap option.

### Validation impact

- JSON and schema contracts are unchanged. Inputs previously accepted through a partial bracket
  match, an unsupported non-bracket ALT, or contradictory duplicate mate coordinates can now have
  fewer accepted events and explicit rejection counts. Valid four-form BND alleles and supported
  symbolic records retain their mate coordinates. No caller threshold, consensus-matching rule,
  confidence, fusion assertion or `reportable` value changes. This is defensive normalization,
  not biological validation; intended-use truth data and the open validation work remain required
  before promotion from draft review.

## 0.7.1 - 2026-09-09 (local engineering candidate)

- Correct Desktop preparation failures that were mislabeled as BAM worker failures.
  Identify inaccessible drive roots, input/output directories, resources and runtime
  components separately, with the failing path and a matching corrective action.
  A missing or disconnected drive no longer prompts reference-bundle installation.
- Replace SAM-text tag discovery with an isolated direct BAM reader, indexed quick
  sampling and an optional cancellable sequential scan. Keep memory/time bounds,
  input-change detection and early stopping on supported paired 5mC information.
- Add authenticated background scan start/status/cancel endpoints, progress counters,
  bounded display caches and typed reasons for incomplete, unsupported and failed scans.
  Before starting a run, refresh evidence from the selected input; a cached positive
  witness is read again rather than trusted solely from file metadata.
- Show the concrete reason, checked reads and elapsed time in Desktop and the local
  browser workspace. Changing the selection invalidates prior results and choices.
  Methylation remains an explicit opt-in; only a successful exhaustive scan can report
  absence, and unsupported or malformed tags remain unresolved.
- Advance Core, Desktop and the local workspace together to 0.7.1. Reuse the installed
  GRCh37/GRCh38 reference families without downloads or analytical policy changes.
- Validation impact: discovery and the operator's module-selection flow change;
  methylation estimators, technical thresholds and the 0.3.0 result schema do not.
  Probe provenance advances to v2. Synthetic boundary, concurrency, cancellation and
  interface tests are required; this is not biological or clinical validation.

## 0.7.0 - 2026-09-08 (local integration candidate)

### Unified platform

- Correct Windows/WSL text decoding so German diagnostics and Unicode paths survive
  redirected stdout/stderr. Decode Linux UTF-8 and native WSL startup diagnostics at
  the byte boundary; do not repair already misdecoded text by character substitution.
- Name the selected reference build on setup actions and retain failure details after
  the resource-status refresh. Keep the setup body scrollable at supported small window
  sizes so long resource statuses cannot hide actions. GRCh37 and GRCh38 remain
  independently selectable.
- Integrate the separate methylation development line with the GRCh37/GRCh38
  reference, panel, annotation, structured ISCN and reporting core from the local
  0.6.2 engineering candidate. The result contract remains independently versioned
  at 0.3.0; methylation and experiment schemas retain their own versions.
- Offer methylation as an explicit optional analysis of an aligned BAM. Before a
  run, the local workspace checks for modified-base information and asks whether
  to include it. Incomplete inspection is unknown, never evidence of absence.
- Present genome and methylation analysis in one guided local workspace while
  retaining the independent two-source mixture, holdout and validation commands.
- Advance Core, Desktop and the local workspace together to 0.7.0. Original source
  workspaces are preserved; this integration does not replace an installed runtime.
- Select the versioned target-BED methylation policy for Adaptive Sampling and retain
  chromosome summaries for lcWGS. Record policy, tool and input provenance independently.
- Re-probe immediately before each start; new selection, changed file metadata or tool
  availability clears the previous decision. Bind regional reports to current assembly
  and methylation-stage checksums, including `NO_CALL` reports and unavailable fractions.
- Support a separately installed modkit through the optional Desktop setting
  `modkitExecutableWsl`. Fix the empty workspace start state and use readable German
  setup text with technical details available on demand.
- Use a read-only Windows process handle for run locks, avoiding signal-zero side effects.

### Integration validation impact

- Optional methylation changes the requested module set and creates an additional
  derived report only after explicit selection. Regional fractions remain research
  measurements, and mixture coefficients are not tumour, cell or DNA-mass fractions.
- Existing reference/dictionary locks, CNV coordinate corrections, structured ISCN
  restrictions and clinical_release_allowed=false remain in force. Study software
  locks must be regenerated prospectively; historical evidence is not relabelled.
- A successful tag probe is not a check of caller/model correctness, CpG semantics
  or analytical performance. A local synthetic pipeline smoke now verifies the real
  modkit 0.4.1 path; broader interoperability, a real-binary CI gate and biological
  validation remain separate requirements.

The following methylation development changes are included in this integration:

### Security

- Validation sample, donor, source and cell-type identities reject Unicode control
  and format characters. A synthetic end-to-end regression showed that a NUL-suffixed
  donor label could previously evade the equality gate for a reused donor. Ordinary
  whitespace normalization and identity case remain unchanged. Unknown-provenance
  checks also normalize surrounding whitespace and case. Fresh source/runtime locks
  are required; declared identity still needs external evidence.
- Nanopolish adapter `nanopolish-call-table-v2` bounds decoded bytes as well as
  compressed input size, limits headers and physical lines before CSV parsing, and
  rejects control characters in read identifiers and read-selection digests. Valid
  input selection digests and estimator mathematics are unchanged. Source provenance
  now records decoded byte count and `bounded-tsv-gzip-v1`; old registrations need
  a fresh software lock. Validation impact is documented in `docs/CLINICAL_VALIDATION.md`.

### Added

- A separate prospective two-source read-holdout study for Nanopolish, SAM/modBAM and
  modkit inputs. It binds the intact-read split, seeds, cohort, matrix and source/runtime
  code before outcomes are examined, reports stratified recovery, and always states that
  donor-independent validation is `NOT_EVALUATED`. The stricter four-sample lane remains
  an additional study; it is not a prerequisite for scoped two-source recovery.
- A prospective independent methylation-validation lane with a versioned fraction/read-budget/
  seed matrix, source and software provenance, stratified recovery criteria and explicit missing-
  evidence `NO_CALL`. Its criteria are research proposals, not validated clinical thresholds.
- Separate Dorado/modBAM MM/ML and modkit per-read adapters, retaining the Nanopolish path.
  Biological-output impact and remaining uncertainty are reviewed in
  `docs/CLINICAL_VALIDATION.md`; only synthetic conformance fixtures belong in Git.
- A public methylation-data candidate catalogue distinguishing technical cell-line controls from
  eligible independent biological evidence and recording unresolved acquisition/provenance gaps.
- A modified-base (methylation) lane. `modkit pileup` is wired into the canonical run graph as
  the version-locked `methylation` stage, aggregating `MM`/`ML` calls into per-region 5mC
  fractions over either canonical chromosomes or the locked target design. It runs only when the
  manifest requests the module, is deselectable like any other component, and its normalized
  report is a validated artifact in the run envelope with its own module outcome, tool record and
  release-bundle checksum. See [`docs/METHYLATION_LANE.md`](docs/METHYLATION_LANE.md).
- A deterministic in-silico tumour dilution series. `plan_dilution_series` lays out the whole
  titration as reviewable data — per-level read budgets, derived seeds, exact subsample arguments
  — without touching a BAM; `execute_dilution_series` materialises it with version-locked
  samtools and verifies every level against the fraction it claims.
- A technical limit-of-detection evaluation over the benchmark reports of a series, reusing the
  `tumor_fraction`/`replicate` strata the benchmark cases already carry. It reports per-level
  detection rates, whether the limit is bracketed by an observed failing level, and refuses to
  report a low level that passed while a higher one failed.
- CLI: `call-methylation`, `dilution-plan`, `dilution-mix` and `lod`; `ontseq run` and
  `ontseq preflight` accept `--methylation-policy` and `--modkit`.
- Technical policies `configs/methylation/modkit.technical.yaml`,
  `configs/benchmark/dilution_series.technical.yaml` and `configs/benchmark/lod.technical.yaml`,
  plus seven exported schemas for the new contracts.
- A standalone paired-source Nanopolish methylation experiment with exact shared markers,
  disjoint calibration and held-out read pools, deterministic fixed-budget mixtures, constrained
  weighted least squares over separate methylated/unmethylated call-rate channels, conditional
  four-state conditional Dirichlet Monte Carlo intervals (methylated, unmethylated, ambiguous
  and missing) with the calibration marker lock held fixed,
  variance-standardized fit diagnostics, explicit per-level `NO_CALL`, aggregate
  `COMPLETED`/`PARTIAL`/`NO_CALL` semantics and a separate technical recovery assessment against
  known fractions. It emits
  schema-validated JSON, CSV and self-contained HTML; see
  [`docs/METHYLATION_MIXTURE.md`](docs/METHYLATION_MIXTURE.md).
- Technical policy `configs/methylation/paired_source_nanopolish.technical.yaml` for the paired-
  source fractions, seeds, marker filters, fit gate, conditional uncertainty and engineering-only
  recovery gates, plus an explicit Nanopolish source-metadata template.

### Fixed

- Per-read modification adapter version `0.1.1` rejects explicit Dorado duplex
  offspring (`dx:i:1`) even without a recognizable program header and rejects
  modkit alignment spans beyond the locked reference. Valid simplex probability
  classification is unchanged. New synthetic regressions cover both boundaries;
  previously registered adapter/code versions require a new prospective lock.
  See `docs/CLINICAL_VALIDATION.md` for validation impact and TSV limitations.
- The per-read modkit adapter now accepts only the source-reviewed 0.6.1 `extract full`
  contract. It rejects older/incompatible tables early, treats alignment ends as exclusive
  and reconstructs ML probability bins from rounded midpoint output so TSV and SAM call
  classification agree at policy boundaries. This pin applies to the new per-read path,
  independently of the existing pileup lane; actual binary interoperability remains untested.
- Replaced independent Dirichlet perturbation of both calibration endpoints with algorithm v2's
  fitted-mixture conditional interval. The former errors-in-variables construction displaced
  low-coverage endpoint intervals toward the centre even when the point fit was exact. Reloading
  a completed JSON report now recomputes its point estimate, fit and seeded interval from the
  stored marker counts, so coordinated result tampering is rejected. Reload validation also
  binds calibration-pool sizes and every level seed to the locked policy, including `NO_CALL`
  reports. Marker contrasts exactly on the configured delta-beta boundary now use a narrow
  numerical tolerance so binary floating-point representation cannot create a false `NO_CALL`.
- Repository safety now rejects Nanopolish call tables and sensitive methylation-mixture
  JSON/CSV/HTML even when renamed, and matching generated filenames are ignored by default.
- The structural-variant stage now runs only when the manifest requests the `sv` module, the
  way CNV and methylation already do. It previously gated on *which policies were supplied*
  rather than on *what the run asked for*, so a manifest declaring `modules: [qc, cnv, report]`
  still drove an SV attempt — and a CNV-only run died with "cuteSV requires --reference-fasta"
  on a reference it had no reason to supply. A run that did not request the module now records
  `NOT_RUN` with the reason naming it a scope statement; a run that *did* request it with no
  caller policy still fails closed.

### Changed

- Preflight answers the methylation lane's preconditions before an envelope exists: policy
  present, reference FASTA present when the pileup is CpG-restricted, target BED present when
  aggregation is over the design, and a warning that MM/ML tag presence can only be established
  by reading the BAM.
- Preflight applies the same scope rule to the SV callers: a run that does not request the
  module is neither told about a missing sniffles or cuteSV nor held to their version locks,
  so the tool section keeps meaning something.
- Stage skip vocabulary is now consistent and distinct: `applicable: false` means the assay has
  nothing for the stage to measure, `requested: false` means the manifest did not ask for it.
  `docs/PIPELINE_EXECUTION.md` documents the three gates side by side.
- The unverified-adapter warning is no longer raised for a methylation stage the run never
  requested, so the line keeps meaning something.
- The paired-source Nanopolish technical policy and documentation now match the implemented
  per-motif LLR threshold (`threshold * num_motifs`), M/U call-rate WLS and the dimensionless
  `maximum_standardized_model_fit_rmse` gate. They also state the required upstream provenance,
  sensitive-output handling and the parser's current policy-bounded in-memory scaling boundary.
  Quantifiability is now separated from recovery: bias, MAE, RMSE and empirical interval coverage
  produce `PASS`, `FAIL` or `NOT_EVALUABLE` under locked but unvalidated technical thresholds.
- The Nanopolish adapter now requires an explicit operator declaration that the inputs are
  biologically distinct, locks sequence-context hashes into exact marker identity, restricts
  known caller provenance to Nanopolish, checks a typed reference-build declaration and displays
  requested and realized fractions separately in HTML.

### Validation impact

- Both lanes are new evidence surfaces and neither is validated. The modkit adapter has **never**
  been executed against the real binary here or in CI; it is declared `unverified_adapter` and a
  run completing that stage is reported as such. Fail-closed behaviour is deliberate and load
  bearing: a BAM without `MM` tags fails the stage rather than producing an empty pileup that
  would read as unmethylated DNA, a region with no site above the coverage floor reports `null`
  rather than `0.0`, and the modkit confidence threshold is pinned in policy rather than
  estimated from the sample.
- A detection limit from an in-silico series characterises software behaviour on one pair of
  BAMs. It reproduces read-fraction effects and nothing about library preparation, input mass or
  capture behaviour at low tumour content, its replicates are not independent specimens, and an
  unbracketed limit is reported as a bound rather than a limit. No number from this lane is an
  analytical or clinical sensitivity.
- No existing lane's output changes. Assembly gains an optional methylation module outcome;
  results without the lane are byte-identical apart from that absence.
- The paired-source methylation coefficient is calibrated to the fraction of held-out read groups
  drawn from source A. It is not tumour purity, blast fraction, cell fraction or DNA-mass
  fraction. Its seeded replicates reuse the same two sources, and its conditional interval holds
  the calibration lock fixed; it does not account for calibration-source uncertainty, donor
  variation, wet-lab dilution, platform transfer, CNV, ploidy or correlated calls.
  Policy values are unvalidated engineering defaults; 5, 10 and 20 percent are test levels, not
  established detection limits. A SMALL public-data smoke can test only technical integration;
  it establishes neither purity nor LoD and does not exercise whole-genome memory scale. Outputs
  are sensitive derived genomic artifacts. The report embeds input fingerprints, software/Git
  identity and typed caller, reference, platform, flow-cell, library-kit and basecaller metadata;
  known cross-source conflicts fail closed, while omitted values remain explicitly unknown and
  produce reviewer-visible warnings. Byte and row caps bound accepted input scale, but the parser
  still holds accepted calls in RAM. Unknown provenance prevents cross-run or biological
  interpretation.
  A technically `COMPLETED` grid may still fail the separate recovery assessment; this prevents
  successful execution from being presented as accurate known-fraction recovery.
  The path changes no existing single-sample biological output.
- The SV gating fix does change behaviour for one case, deliberately: a run whose manifest
  omits `sv` but whose configuration supplied caller policies previously produced structural
  variant evidence and now records `NOT_RUN`. That evidence was outside the declared analysis
  scope; a run asking for CNV was never asking for SV. The change is visible rather than
  silent — the module outcome, run report, HTML and XLSX all carry the reason — and a manifest
  that lists `sv` behaves exactly as before. Runs already in an envelope are unaffected: the
  stage signature change re-runs the stage rather than reinterpreting an existing artifact.
## 0.6.2 - 2026-09-08 (local engineering candidate)

### Added

- Two build-isolated GRCh37 Adaptive-Sampling profiles are now available:
  `AML_AS_111_GRCh37` for the full GENCODE-19 publisher dictionary and
  `AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25` for the ordered UCSC-hg19 25-contig dictionary.
  Both resolve only the dedicated `AML_AS_111_GRCh37_v1` bundle.
- The immutable `GRCh37_GENCODE19_HG19_v2` reference family retains the native full
  GRCh37.p13 contract and adds a separately generated exact UCSC hg19 Canonical-25 analysis
  FASTA. The hg19 contract contains `chr1`-`chr22`, `chrX`, `chrY` and the original UCSC
  `chrM=16571`; its FASTA, FAI and source checksums are pinned independently.
- The GRCh37 panel preserves the exact laboratory GRCh38 selection-design lineage through a
  checksum-pinned, build-time coordinate projection. Of 111 intended source intervals, 110 map
  at the locked 0.99 threshold and 109 are reciprocal-exact; 108 obtain native GENCODE-19
  analysis ROIs. `ACACA` requires roundtrip review, while `CT45A2`, `IGH` and `GPR128` remain
  explicit mapping/ROI gaps that cannot produce negative observability claims.
- The panel authority now ships and validates the checksum-pinned forward and reciprocal-audit
  outputs, binds the mapping lock to the declared source bundle/resource, requires every native
  ROI to overlap its mapped selection, and records final-selection/ROI counts in JSON, HTML, XLSX
  and the local Live Workspace.

### Fixed

- Reference readiness and profile execution now resolve the exact contract-specific FASTA and
  FAI. A UCSC hg19 BAM can no longer be paired with the native 16,569-base mitochondrial
  reference, and no alias conversion, reheadering, liftover or silent fallback is attempted.
- Bundle publication stages files beside the final destination so Windows ACL inheritance stays
  valid when the resource root is shared with WSL. Explicit repair can replace an unreadable
  regular bundle tree while preserving fail-closed manifest-last publication and rollback.

### Validation impact

- The new profiles add an engineering execution path; they do not establish equivalence of
  enrichment performance between builds. The coordinate projection is performed only while
  constructing the immutable panel authority. Analyses never invoke liftover, reheader BAMs or
  fall back to the GRCh38 panel or another dictionary contract.
- The GRCh37 selection design is `110/111 mapped`, not a claim that all 111 source targets are
  represented or that every projected boundary is reciprocal-exact. Native GENCODE-19
  ROI/transcript analysis is limited to 108 targets; `ACACA` roundtrip review and the three
  mapping/ROI gaps remain report-visible technical review items.
- This remains unsigned Research Use Only engineering software. Adaptive Sampling on hg19 and
  the reference/panel changes require analytical validation and expert review before any
  clinical use.

### Release packaging

- Python Core, Windows Desktop, WSL runtime check, Live Workspace and operator documentation
  advance together to `0.6.2`; no functional change is shipped under the prior version number.
- The offline Windows test package installs Core 0.6.2 into a new isolated runtime prefix and
  uses a fresh resource root for the v2 GRCh37/hg19 family; existing 0.6.0/0.6.1 installations
  are not overwritten.
- UCSC source data and build-time mapping inputs retain their individual provenance and terms.
  Institutional or commercial distribution must independently satisfy the applicable source
  licences; absence from the runtime package is not a claim that licence review is complete.

## 0.6.1 - 2026-09-07 (local engineering candidate)

### Added

- Result schema `0.3.0` adds an explicit, structured ISCN proposal contract with the states
  `NOT_REQUESTED`, `NOT_ASSESSED`, `NO_RENDERABLE_CANDIDATE` and
  `PARTIAL_EVENT_LEVEL`. Every result event now carries a rendered, omitted,
  policy-excluded, not-requested or assessment-blocked disposition; legacy schemas remain as
  `LEGACY_UNSPECIFIED`.
- The regular CNV result path can emit traceable, event-level ISCN proposal fragments under
  `event-fragments-v0.2-unvalidated`. Eligible fragments originate only in the typed CNV report
  and retain their event identifiers, confidence, reportability and selection rationale.
- ISCN proposal assessment is bound to the selected GRCh37 or GRCh38 resource context, including
  the BAM dictionary contract, reference-bundle identity, reference-lock SHA256, cytoband
  release/SHA256 and annotation-cache SHA256. Missing, unreadable or inconsistent provenance
  blocks assessment instead of permitting a fallback.

### Fixed

- Prevent a near-whole-chromosome CNV classification (for example a segment meeting the 90%
  policy threshold) from becoming `+chr`/`-chr` notation unless its coordinates exactly span the
  locked reference contig. CNV artifacts are now also rejected before merge when sample ID or
  genome build disagrees with the manifest.
- Keep fragmentless workbook output aligned with `NOT_REQUESTED`, `NOT_ASSESSED` and
  `NO_RENDERABLE_CANDIDATE` instead of labelling every state `NO_CALL`; result validation now binds
  assessed ISCN states to CNV execution and enforces state-specific event dispositions.
- Isolate concurrent Desktop instances on separate loopback ports. If the configured port is
  already occupied, a new instance selects a free ephemeral port instead of bootstrapping
  against or terminating the existing ONTSeq service.
- Bind every Desktop launch to a random 128-bit instance ID and require `/api/config` to echo
  that ID together with the exact resource, output and allowed-input roots before committing
  the connection. A listener that wins a port race is rejected and retried on a fresh port;
  failed bootstrap candidates are never reused.
- Preserve the exact Desktop-selected profile when opening the authenticated Live Workspace.
  Both Desktop and browser now require that `GET /api/config` advertises the requested profile;
  there is no silent substitution with the service's first profile. Completed results without
  the exact run-bound `reference_context.profile_id` are also rejected.
- Keep the synthetic demonstration result's manifest and resolved resource context on the same
  profile identifier, so the authenticated workspace exercises the same fail-closed profile
  binding as real result envelopes.

### Validation impact

- This change can alter report content and the recorded ISCN module state, but it does not alter
  CNV or SV caller output, reportability or clinical-release eligibility. Result schema `0.3.0`
  is required for the structured proposal states; results using schema `0.1.0` or `0.2.0` retain
  legacy ISCN semantics and are not silently upgraded.
- The renderer is deliberately CNV-only. Automatic conversion of SV, BND, translocation or
  inversion evidence is disabled because reciprocity, balance and derivative structure are not
  established by the current evidence contract.
- A partial proposal is not a complete karyotype. The software does not infer chromosome count,
  sex-chromosome complement, clonality, phase, derivative structure, balance or normality, and
  an empty fragment set is never presented as a normal result.
- All generated fragments remain explicitly unvalidated, require expert cytogenetic review and
  have `clinical_release_allowed=false`. This change makes no claim of full ISCN 2024 conformance
  and does not constitute analytical or clinical validation.

### Release packaging

- Core, Windows Desktop, WSL installer identity, live workspace and operator documentation are
  advanced together to `0.6.1`; the result schema remains independently versioned at `0.3.0`.
- The unsigned Windows engineering package retains the pinned dual-build GRCh37/GRCh38 runtime
  contract and installs Core 0.6.1 into a separate runtime prefix without overwriting 0.6.0.
- The build-time UCSC chain and `liftOver` executable are provenance inputs, not redistributable
  runtime dependencies, and are not included in the wheel or Windows package. Institutional or
  commercial distribution must independently satisfy the applicable UCSC licensing terms.

## 0.6.0 - 2026-09-02 (local engineering candidate)

### Added

- Native GRCh37.p13 / GENCODE 19 reference family and an explicit lcWGS profile, isolated
  from the existing four GRCh38 profiles. Native GRCh37 annotations deliberately have no MANE.
- Explicit `AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25` profile for ordered BAM dictionaries
  containing exactly `chr1`-`chr22`, `chrX`, `chrY` and UCSC hg19 `chrM=16571`.
- Build-bound gene/transcript/cytoband cache consumption in SV and CNV annotation.
- A separately versioned coordinate-free GRCh37 hematology review bundle, preserving the
  original source scope and without transferring GRCh38 coordinates or target coverage.
- Local live report workspace with real profile/BAM selection, backend execution, actual
  stage/result/provenance presentation and authenticated original HTML/XLSX/JSON downloads.

### Safety and validation impact

- Correct QDNAseq CNV export coordinates for both builds: convert native one-based,
  closed bin/segment starts to zero-based, half-open starts, leaving the end unchanged.
  RDS/ACE calculations retain their native coordinates. Earlier engineering outputs
  with unconverted starts are superseded and must be regenerated, not relabeled.
- Correct positive Umap context handling for both builds: uniquely mappable interval overlaps
  remain traceable but no longer add an artifact-priority penalty. Missing Umap rows are not
  inverted into low-mappability evidence. Other artifact flags retain their technical penalties.
- This is a source/package integration candidate, not a clinical release or a claim that a
  version label proves execution. No caller threshold or reportability policy is promoted.
- GRCh37 dictionary selection is explicit and fail-closed. `AML_LCWGS_GRCh37` requires the
  complete GENCODE 19 publisher dictionary (including scaffolds/patches/haplotypes), while
  `AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25` requires exactly the UCSC hg19 25-contig order and
  `chrM=16571`. The latter intentionally differs from native GENCODE 19 `chrM=16569`; no BAM
  reheader, reference substitution, automatic fallback or liftover is performed.
- The hg19 contract changes only accepted input/reference binding. Nuclear GRCh37 annotation,
  caller policies and reportability boundaries are unchanged and require a separate analytical
  validation lane before any diagnostic use.
- Annotation can change review context; both builds require independent regression testing.
  Native GRCh37 Adaptive Sampling is not offered without a controlled build-matched design.
- The live workspace cannot fabricate data from demo fixtures, reassign a completed run's
  reference build or perform clinical sign-off. Source-mode service runs are marked
  LOCAL_WORKTREE instead of borrowing an older installed runtime's commit.
- Tests, source provenance and remaining package-validation gates are recorded in
  docs/GRCH37_INTEGRATION.md and the external engineering handoff artifacts.

## 0.5.3 - 2026-08-29

### Added

- `HEMATOLOGY_v3` adds 38 locked rearrangement review patterns and 74 source-attributed
  pathology associations. Exact public records are derived from a 2026-08-29 CIViC snapshot and
  filtered through Disease Ontology v2026-07-31 under `DOID:2531`; the curated
  `PICALM::MLLT10` literature record remains included.
- HTML and XLSX reports now expose a dedicated key-findings layer, fusion/rearrangement
  assessment, Adaptive-Sampling ROI coverage and separate technical SV review queues.
- cuteSV output provenance and tool identity are retained when the result is assembled.

### Changed

- Fusion module status now describes the assessment actually performed on consensus breakpoint
  evidence instead of remaining `NOT_RUN` whenever no separately validated fusion assertion exists.
- The report displays `BENCHMARK_REQUIRED` alongside review priority instead of presenting an
  unexplained `reportable=false` value as if it meant biologically irrelevant.
- HTML keeps the complete result in JSON/XLSX but limits the expanded technical appendix to a
  usable review subset. XLSX preserves all normalized events and adds focused review sheets.
- HTML and XLSX display associated pathology names and DOIDs beside matching rearrangement
  candidates, rather than requiring a reviewer to infer disease context from the gene pair.
- All four active GRCh38 profiles pin `HEMATOLOGY_v3`; versions 1 and 2 remain packaged for
  provenance and backward readability.

### Fixed

- The hematology knowledge builder now writes canonical LF JSON on every operating system, so
  the checked-in bundle size and SHA256 remain identical on Windows and Linux.
- Git now preserves every checksum-pinned knowledge, panel and reference-fixture byte stream on
  Windows, and repository safety rejects newly pinned resources that lack this protection.
- The repaired dependency lock now records the 0.5.3 root-package version used by runtime
  provenance and release checks.

### Validation impact

- The knowledge update can move a source-matched exact pair or order-reversed
  `MLLT10::PICALM` breakpoint candidate into the hematology review queue and attach relevant
  disease vocabulary. It does not establish transcript productivity, pathogenicity, diagnosis,
  prognosis, or analytical reportability. Existing SV/CNV caller thresholds are unchanged, and
  every surfaced candidate still requires expert and assay-specific validation.

## 0.5.2 - 2026-08-28

### Added

- Two explicit GRCh38 Canonical-25 profile variants now accept BAM dictionaries containing
  exactly `chr1`-`chr22`, `chrX`, `chrY` and `chrM`: `AML_LCWGS_GRCh38_CANONICAL25` and
  `AML_AS_111_GRCh38_CANONICAL25`.
- Desktop exposes the dictionary contract with each profile so that the operator can select the
  contract matching the reference used for alignment before starting a run.

### Changed

- The established `AML_LCWGS_GRCh38` and `AML_AS_111_GRCh38` profiles remain `exact_full`: their
  BAM dictionaries must still match the complete ordered Primary-Assembly `ReferenceLock`.
- Canonical-25 is a profile-level projection of the same pinned GRCh38 reference, annotation,
  panel and knowledge bundles. Selecting it performs no liftover, enables no fallback and does
  not download or install a second multi-GiB reference bundle.

### Validation impact

- This release adds an explicit input-dictionary compatibility contract. It changes no caller,
  CNV/SV threshold, panel coordinate, annotation release, evidence policy or reportability rule.
  A BAM matching neither exact contract still fails before pipeline execution; no additional
  analytical or clinical validation claim is introduced.

## 0.5.1 - 2026-08-28

### Fixed

- Packaged configuration defaults now resolve from the installed ONTSeq release rather than the
  process working directory. Desktop profile services additionally receive absolute paths for the
  cuteSV, Sniffles2/cuteSV consensus and SV evidence policies, preventing profile startup from
  failing when the application is launched outside a repository checkout.
- The Windows/WSL runtime preflight now verifies the complete policy and tool contract before a
  run. The packed environment includes pinned cuteSV 2.1.3, retains mosdepth, rejects duplicated
  `share/ontseq/configs/configs` trees and exercises the relocated service from an unrelated
  working directory in CI.

### Validation impact

- This patch repairs runtime discovery and packaging only. It changes no scientific threshold,
  caller policy, CNV/SV interpretation rule, profile resource identity or reportability boundary;
  no additional analytical or clinical validation claim is introduced.

## 0.5.0 - 2026-08-28

### Added

- The original 111-target GRCh38 Adaptive Sampling BED and companion region list are now
  byte-provenanced in `AML_AS_111_GRCh38_v1`; a deterministic importer creates a separate
  0-based half-open selection BED and keeps `IGH_REVIEW_REQUIRED` unresolved.
- A GRCh38 annotation-cache consumer compiles unbuffered gene-body ROIs and ranked panel
  transcripts without inferring targets from the selection buffers. Both SV breakpoints can retain
  gene, preferred transcript, exon/intron, CDS phase, cytoband, repeat, blacklist, and mappability
  context. Fusion evidence exposes Gene A/B, orientation, and an explicit `unknown` frame default.
- A manifest-pinned GRCh38 reference installer now stages downloads, verifies byte size and
  SHA256, builds the FASTA dictionary/reference lock and deterministic GENCODE/MANE/cytoband
  SQLite cache, then activates the bundle atomically. Status, checksum-limited repair and offline
  import never query a remote service. The installable GRCh38.p14/GENCODE 50/MANE 1.5 recipe pins
  exact sizes and SHA256 for all nine publisher artifacts; GENCODE transfers were additionally
  checked against the publisher MD5 index. A miniature bundle exercises the same path in CI and
  the multi-gigabyte installation smoke remains explicit opt-in.
- A CNV cytoband engine retains raw overlaps, applies the configurable 66% fraction-of-band rule,
  merges only adjacent same-direction bands on the same arm, and represents whole-chromosome calls
  separately. The existing AML rearrangement resource is manifest-pinned as `HEMATOLOGY_v1`.

### Validation impact

- Every active panel interval starts one base earlier than a literal BED interpretation of the
  laboratory source, increasing the locked span by 111 bases. Edge coverage and breakpoint target
  membership can change accordingly. Transcript and cytoband summaries can also change review
  ordering. Regression tests lock these transformations, but no analytical sensitivity,
  specificity, LoD, negative-observability, fusion-frame, or reportability claim is introduced.
- Reference/FAI/lock inconsistencies and coordinate-ambiguous panel manifests that older loaders
  could accept are now rejected before analysis. The CNV affected-band cutoff remains 0.66, but is
  now an explicit field of the versioned QDNAseq/ACE policy and recorded in stage/sidecar
  provenance instead of being an untracked engine default.

### Changed

- Interactive Desktop/service profile preflight now checks every pinned Reference, Panel and
  Knowledge resource for manifest validity, presence and declared size without re-hashing
  multi-gigabyte files. `ontseq references validate` remains the explicit full SHA256 audit.
- `/api/config` advertises only locally resolvable GRCh38 profiles, missing profiles fail as HTTP
  400, and Desktop uses the Core-derived `<sample>-<UTC timestamp>` run ID returned by the service.
- Explicit reference repair now transactionally restores the complete pinned GRCh38 profile
  family (reference, panel, knowledge and profile manifests) with rollback, so damaged dependent
  bundles no longer require manual deletion. Profile-backed service runs retain the configured
  QC/SV/coverage policies, minimum depth and component-version selection.

### Fixed

- `ontseq references repair GRCh38_GENCODE50_MANE1.5_v1` now repairs the complete pinned
  profile-resource family, including `HEMATOLOGY_v1`, `AML_AS_111_GRCh38_v1` and both profile
  manifests, with staged validation, path-atomic replacement and rollback instead of requiring
  operators to delete divergent resources manually. Repair and official-ID import also require
  the exact catalog Source-/Generator contract; changed sources or derivations require a new
  bundle ID/version.
- Native UCSC hg38 cytobands now ignore the unnamed chrM placeholder instead of rejecting the
  publisher table; named cytogenetic bands remain strictly validated.
- Pseudoautosomal panel symbols such as `P2RY8` are disambiguated by the explicitly declared
  source chromosome. Historical or coordinate-conflicting labels remain unresolved.
- Result assembly fingerprints the annotated SV consensus, so a changed breakpoint annotation or
  knowledge context cannot resume a stale `PipelineResult`. Annotated BND/translocation candidates
  with fusion evidence now also appear in the XLSX fusion worksheet.
- Full GRCh38 technical-context BEDs now use a path-backed, contig-lazy compact interval index with
  bisect/block-max point queries instead of retaining every row and rescanning it for each
  breakpoint.
- Legacy manifest runs without `--reference-fasta` no longer activate the optional cuteSV caller
  from its default policy path. Cramino histogram counts are written to an explicit temporary file
  so the primary stdout stream remains valid JSON before the numeric sidecar is normalized.
- SV and QC re-execution clears stale caller/consensus and histogram artifacts before new work;
  failed or evidence-free reruns therefore cannot expose an earlier run's sidecars. Official-ID
  repair/import also rejects changed source or generator contracts under an unchanged bundle
  identity, while custom bundle IDs remain importable.
- Wheels and containers now carry the immutable GRCh38 authority/configuration assets under the
  installation prefix, and Windows Python 3.11 rejects junction/reparse-point resource paths at
  the same mutation boundaries as Python 3.12.

## 0.4.1 - 2026-08-27

### Added

- A productive cuteSV 2.1.3 adapter now runs beside Sniffles2 with version-locked parameters,
  atomic VCF finalization and caller-specific normalized evidence.
- A versioned consensus layer canonicalizes reversed BND/TRA representations and clusters
  compatible DEL/DUP/INV/INS/translocation records by type-aware breakpoint, overlap, length and
  available orientation rules. Source caller IDs and evidence remain traceable.
- Build- and checksum-locked interval annotation now supplies separate breakpoint genes,
  cytobands, nearest-gene distances and repeat/tandem-repeat/segmental-duplication/blacklist/
  mappability/centromere/telomere context flags. A preparation script records both original and
  normalized resource hashes.
- Adaptive Sampling SVs now carry explicit breakpoint depths and observability states. A small,
  versioned AML rearrangement-pattern resource prioritizes recurrent patterns without asserting a
  confirmed fusion.
- The technical evidence score is fully represented by a versioned policy. HTML and XLSX show a
  filterable high/moderate review queue, Gene A/B, caller support/consensus, coverage, artifact
  context, AML relevance and validation status while retaining the complete technical table.
- Synthetic regression coverage includes cuteSV-only normalization/execution, Sniffles-only
  behavior, two-caller consensus, canonical reversed BND orientation, within-caller clustering,
  separated nearby events, annotation build refusal, repeat context, Adaptive Sampling partial and
  insufficient coverage, AML patterns and the permanent non-reportable boundary.

### Fixed

- Corrected the productive cuteSV positional argument order to `BAM REF VCF WORKDIR` and exercised
  it through the real-tool CI smoke test.
- Resolved strict type-check findings in deterministic SV ordering and provenance resource locks.

### Validation impact

- This change can alter technical SV clustering, priority tiers and report order. All new matching,
  observability and scoring cut-offs are explicitly `technical_defaults_only`; none is a validated
  sensitivity, specificity, PPV, LoD or clinical threshold.
- Every automated event remains `reportable=false`. Multi-caller support and AML pattern matches do
  not set `analytically_validated`, assert a fusion or establish a somatic origin. Independent
  public technical benchmarks and orthogonally characterized AML specimens remain required before
  clinical use.

## 0.4.0 - 2026-08-27

### Changed

- Research-only/validation status no longer participates in QC verdict calculation.
  Configured measurable gates alone determine `PASS` or `FAIL`; absence of configured gates
  remains `WARN`.
- The research-use disclaimer remains visible as report text. It does not change QC status.
  Non-reportable findings and the separate expert-review workflow remain unchanged.

### Fixed

- Tool output is decoded as UTF-8 with replacement instead of through `text=True`, which
  used the platform locale. Under `LC_ALL=C` — the default in many containers, cron
  environments and freshly installed WSL distributions — that resolved to ASCII, and a
  single non-ASCII byte in a tool's own banner raised `UnicodeDecodeError` from inside the
  runner. A version probe could therefore fail a stage on one machine and pass on another
  from identical inputs. `run_to_file` already decoded this way; both paths now agree.
- The local service refused to escape its allowed roots for BAM paths but built the review
  envelope path by joining two raw URL segments onto the output directory, so
  `POST /api/review/../<sample>` reached an envelope outside `--output-dir` and could append
  a sign-off there. Both identifiers are now matched against the manifest's `run_id` /
  `sample_id` contract and the resolved path is confirmed to be inside the output directory.
- The service decided "one analysis at a time" by asking `Jobs.running()` and registering the
  job afterwards. Between those two steps it read the reference lock and resolved the BAM and
  its index, so simultaneous requests all passed the gate and started concurrent pipelines,
  each sized for the whole machine. The envelope lock does not catch this: two runs with
  different run ids take different envelopes. The slot is now claimed atomically, and reusing
  a run id already in the table is refused rather than dropping the earlier job's record.
- A dangling symlink or an unreadable entry in a browsed directory made `stat()` raise inside
  the request, which ended the response without a body and showed the operator a network
  error for a directory whose other files were usable. Such entries are now skipped, counted
  and reported to the page as `unreadable_entries`, so a BAM that is present but invisible is
  distinguishable from one that is absent.
- BAM name lookup now builds candidates from directory entries and resolves each candidate
  inside its declared allowed root before reading metadata. A BAM symlink that points outside
  the service boundary is omitted rather than exposing its target's size or presenting it as
  selectable input.
- `GET /` handed out the session token without the loopback `Host` check that every `/api/`
  route applies, so a page on an attacker-controlled name resolved to `127.0.0.1` was
  same-origin with the service and could read the token out of the response. The API check
  still refused the stolen token; the page route now applies the same check.

### Added

- `tests/test_cli_surface.py`. Both parsers, the dispatcher and the command overview had no
  test coverage at all. The tests tie the three descriptions of the command set together —
  the overview listing, the dispatch set and the parsers themselves — and check that a bad
  invocation leaves through `SystemExit` rather than a traceback.
- `scripts/check_version_consistency.py`, run by `make versions` and by CI, refuses a tree
  whose declared versions disagree. Five files state what a release is — package metadata,
  the package, the citation record, the desktop project and the README's status section —
  and nothing tied them together. They had drifted: the README described the desktop on
  `main` as the v0.2.1 engineering path while every other declaration said 0.3.4.

### Fixed (continued)

- Preflight was blind to everything adaptive sampling needs. `mosdepth` had no entry in
  `TOOLS_BY_STAGE`, so no run of any input kind checked for it, and nothing checked the
  target BED or the target-coverage policy. All three fail the target-coverage stage closed,
  and a `FAILED` target-coverage stage fails the whole run — but only after the envelope
  exists and the lock is taken, which is exactly what preflight was written to prevent.
  Preflight now probes Mosdepth with the adapter's own version parser, enforces its policy
  lock, and parses the declared target BED. A tool is escalated from warning to blocking when
  the stage needing it is one *this* run cannot do without, so a missing Mosdepth refuses an
  adaptive-sampling run and only warns on an lcWGS one, which correctly records targets as
  out of scope.
- `ontseq preflight` accepted `--target-coverage-policy` and `--components` and used
  neither. The policy was dropped on the way into `PreflightRequest`, and the component
  selection was resolved only for `run` and `serve` — so preflight could check the default
  policies while `ontseq run` would use the ones a selection names. Both are honoured now,
  for the alignment, basecall, SV and target-coverage policies alike.
- `ontseq` never listed `validate-reference` in its command overview, so a working command
  was reachable only by already knowing its name — the discoverability problem the overview
  exists to solve. It is listed now, and a test fails if any command is missing from the
  overview or from the dispatch set.

### Changed

- `target_coverage.mosdepth_version` is public, so preflight probes the binary with the same
  parser the stage uses rather than a second implementation of it.
- `docs/PIPELINE_EXECUTION.md` described target coverage as "not implemented" with no tool,
  and grouped it with CNV as "not wired in". It has been wired into the canonical runner and
  verified against real Mosdepth since 0.3.4; the stage table, the scope note and the
  preflight section now say what the code does.
- Core, Desktop, WSL runtime and Windows bundle identity advance together to 0.4.0. The
  consistency gate now covers the executable Desktop label, runtime installer, Desktop CI
  artifact version and current operator documentation in addition to the original five
  release declarations.

### Validation impact

The QC verdict update is a status-semantics change only. It does not alter QC measurements,
configured thresholds, CNV/SV evidence, reportability, ISCN output or clinical-release
safeguards.

The boundary fixes do not change an analytical threshold, caller parameter, normalization
rule or result contract, and no stage produces a different outcome from the same inputs on a
machine where the previous version ran to completion. The decoding fix removes a
locale-dependent failure mode, so a run that previously aborted mid-stage on a non-UTF-8
locale now proceeds; it does not change what that run reports. Nothing here makes any output
more validated than it was.

The 0.4.0 version advance changes release and runtime identity only; it does not itself add
analytical or clinical validation.

## 0.3.5 - 2026-08-25

### Added

- Pinned `QDNAseq.hg38` runtime resources for the GRCh38 CNV lane alongside the existing
  GRCh37/hg19 resources.
- A canonical real-tool GRCh38 CI gate that exercises 100/500/1000-kbp QDNAseq + ACE,
  verifies the expected synthetic chromosome 7 loss and chromosome 8 gain, checks the
  generated HTML/XLSX/JSON artifacts and proves content-addressed resume.

### Changed

- Desktop and Python Core version declarations were advanced together to 0.3.5.
- GRCh38 Desktop runs may request the integrated QDNAseq + ACE CNV stage now that the exact
  hg38 annotation resource is pinned and exercised in real-tool CI.

### Validation impact

The new GRCh38 lane has deterministic engineering verification only. It has not undergone
cohort-level analytical or clinical validation. The ACE cellularity/ploidy fit, CNV warning
boundaries and whole-chromosome classification fraction remain research-only engineering
parameters and must not be interpreted as validated clinical cut-offs.

## 0.3.4 - 2026-08-22

### Added

- Per-run component selection: a `RunComponents` document names the provider and the exact
  tool version for each stage that runs an external program. The version is checked against
  the tool's own probe before the stage executes, and a mismatch fails that stage naming
  both versions. Deselected stages are recorded as `NOT_RUN` with the selection named.
- `configs/components/default.yaml` and `configs/components/legacy_sniffles_2.4.yaml`, plus
  a Sniffles 2.4.0 policy, so the structural-variant component of the historical pipeline
  can be reproduced for comparison without editing code.
- Target coverage is wired into the canonical runner as `StageId.TARGET_COVERAGE`. An
  adaptive-sampling run without a target-coverage policy now fails closed instead of
  recording the assay's central QC as absent; a non-adaptive run records that the stage does
  not apply, which is a scope statement rather than a coverage result.
- `AssaySpec.target_bed_role`, distinguishing an unbuffered analysis ROI from a buffered
  selection panel, and carried into the target-coverage report and its limitations.
- First non-synthetic target design: a 111-target, 17.03 Mb buffered GRCh38 panel under
  `configs/panels/`, with a lock file, a reproducible generator script and de-identified
  per-target coverage expectations. Provenance and open questions in
  `docs/PANEL_PROVENANCE.md`.
- CI proves per-target coverage end to end inside a run envelope, proves an lcWGS run
  reports the stage as out of scope, and proves a component pinned to an uninstalled version
  fails closed.
- `docs/COMPONENT_SELECTION.md` and `docs/LEGACY_COMPARISON.md`, the latter mapping the
  historical pipeline's outputs onto this one and proposing MOLM13 as the first positive
  control.
- Version-locked Mosdepth target-coverage adapter for Adaptive Sampling aligned-BAM runs.
- Strict target-BED and Mosdepth output normalization with exact interval reconciliation,
  per-region mean depth and configurable threshold fractions.
- Dedicated technical policy for target coverage and a synthetic real-tool CI path.
- Documentation that separates the unbuffered analysis ROI BED from operational Adaptive
  Sampling selection regions.
- Precision-first portable-report view model that keeps execution state, caller evidence,
  knowledge-resource annotations and interpretation boundaries separate.
- Synthetic regression tests for `NO_CALL`, `FAILED`, zero-valued evidence, HTML escaping,
  offline rendering and RUO/reportability boundaries.

### Changed

- Reorganized the portable HTML report around reviewer tasks: review blockers first, then module
  execution, QC, normalized event evidence, ISCN proposal, warnings and provenance.
- Expanded event evidence display to preserve support reads, local coverage, VAF, quality,
  strands, precision and filters without converting missing values into zero or vice versa.
- A `reportable: true` pipeline flag is now displayed with an explicit statement that the RUO
  report remains not clinically validated.

### Changed

- `ontseq` with no command lists both command groups. Execution commands were previously
  invisible to `--help`.

### Validation impact

Adaptive-sampling runs now produce per-target coverage as part of the canonical run rather
than only through a separate command. The default `1x`, `10x`, `20x` and `30x` thresholds
remain engineering bins only and do not define assay adequacy, CNV or fusion reportability,
biological negativity or a clinical no-call.

The coverage expectations table is a sanity reference derived from historical runs; it is
not an adequacy gate, not a reportability threshold and not a no-call definition. The panel
is `derived_unconfirmed`: it is reconstructed from laboratory coverage tables rather than
copied from the selection file used at sequencing time, and one target label contradicts its
coordinates.

The HTML changes are presentation-only: they do not alter callers, thresholds, normalized events,
module status, confidence, reportability or ISCN generation. They deliberately make `FAILED`,
`NO_CALL`, `NOT_RUN`, unavailable evidence and unvalidated interpretation boundaries harder to
misread. Clinical/analytical validation status is unchanged.

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
