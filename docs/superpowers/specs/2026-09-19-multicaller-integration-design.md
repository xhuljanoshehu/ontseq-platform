# ONTSeq Multi-Caller Integration Design

Date: 2026-09-19
Status: approved for staged implementation
Branch: `feat/multicaller-integration`
Base: `b6e0338b6f042239a1399e3438e6c2bab4e48b08`
Scope: Research Use Only; Human Review Required

## 1. Goal

Integrate the selected caller families into ONTSeq without turning caller agreement into biological truth and without replacing the existing single-provider production-style stages.

The selected families are:

- QDNAseq+ACE
- Spectre
- Sniffles2
- cuteSV
- SAVANA
- Severus
- Wakhan
- ichorCNA

"Integrated" has four separate meanings which must not be conflated:

1. known to the caller catalog;
2. eligible for a declared data regime and input set;
3. technically executable by a pinned runtime adapter;
4. analytically validated for a particular intended-use study.

A catalog entry is not a runtime qualification, and a runtime qualification is not analytical or clinical validation.

## 2. Existing system boundary

The existing `RunComponents` contract selects one provider per main pipeline stage. That remains intact. Multi-caller comparison is an additive subsystem with explicit parallel lanes rather than a rewrite of the canonical runner.

Existing QDNAseq+ACE, Sniffles2 and cuteSV behavior must remain stable. New comparison lanes may consume the same BAM but must have independent lane identity, policy identity, runtime identity, output directory, state and provenance.

## 3. Domain separation

The subsystem separates four analytical domains.

### Genome-wide copy number

- QDNAseq+ACE
- Spectre

Outputs may include depth-derived bins, segments, copy number, fit parameters and caller-native QC. QDNAseq multi-resolution evidence and ACE alternative fits remain retained.

### Structural variation

- Sniffles2
- cuteSV
- Severus
- SAVANA-SV

Breakpoint and SV evidence remains distinct from absolute copy number. A deletion-like SV does not imply an absolute CN value unless the caller actually reports one under a registered contract.

### Tumor CNA / allele-aware analysis

- SAVANA-CNA
- Wakhan

Purity, ploidy, allele-specific copy number and LOH are model-dependent outputs. Alternative solutions and upstream dependencies must be retained.

### cfDNA / ULP-WGS transfer lane

- ichorCNA

ichorCNA is a separate assay/data-regime lane. It is never silently treated as validated for marrow lcWGS or Adaptive Sampling.

## 4. Caller catalog

Every caller mode is identified by a stable provider ID and mode ID. Related modes that change analytical assumptions are distinct entries.

Initial entries:

- `qdnaseq_ace:multibin`
- `spectre:depth_only`
- `spectre:sniffles_supported` (future gated mode)
- `sniffles2:standard`
- `sniffles2:mosaic` (separately qualified)
- `cutesv:standard`
- `severus:paired`
- `severus:single_sample` (separately qualified; no automatic somatic claim)
- `savana:paired`
- `savana:tumor_only`
- `wakhan:phased_cna`
- `ichorcna:ulp_wgs`

ONT-Spectre, if added later, receives a distinct provider ID rather than being treated as another independent vote from upstream Spectre.

## 5. Eligibility inputs

Eligibility is explicit and fail-closed. Unknown or absent inputs are never imputed.

The planner can reason about:

- assay mode / data basis;
- genome build and locked reference identity;
- tumor BAM;
- matched normal BAM;
- generic single-sample BAM;
- coverage evidence;
- target design for Adaptive Sampling;
- SNP/BAF input;
- phased SNP or haplotype input;
- panel-of-normals input;
- directly supporting upstream caller evidence;
- requested analysis mode;
- registered runtime availability.

Examples:

- `savana:paired` requires tumor and matched normal inputs. Missing normal blocks the lane. It must not silently downgrade to `savana:tumor_only`.
- `wakhan:phased_cna` requires its registered allele/phasing inputs. If it is configured to use Severus breakpoints, the Severus lane is a declared parent dependency.
- `spectre:depth_only` must not consume Sniffles output. A supported Spectre mode is a separate lane and declares that dependency.
- `ichorcna:ulp_wgs` is eligible only for a registered low-pass/cfDNA-compatible study regime, not merely because a BAM exists.
- Adaptive Sampling does not become equivalent to unbiased genome-wide lcWGS.

## 6. Lane and dependency contract

A `MultiCallerPlan` contains independent `CallerLane` entries.

Each lane records at least:

- lane ID;
- provider and mode;
- analytical domain;
- specimen/run identity;
- genome build and data basis;
- required input roles and their fingerprints;
- caller version/runtime identity when qualified;
- policy identity;
- direct parent lane IDs;
- selection state;
- eligibility decision and reasons.

Dependencies form a directed acyclic graph. Cycles fail before any execution.

Changing a parent lane identity invalidates dependent-lane plan identity, but must not invalidate unrelated depth-only lanes.

## 7. States

Planning and execution states remain distinct.

Planning decisions:

- `ELIGIBLE`
- `INELIGIBLE`
- `RUNTIME_UNAVAILABLE`
- `NOT_REQUESTED`

Execution outcomes, once a runtime exists:

- `COMPLETED`
- `NO_CALL`
- `FAILED`
- `NOT_ASSESSABLE`
- `NOT_RUN`

An ineligible or unavailable lane never becomes a biological negative.

## 8. Evidence retention and identity

Caller-native files are fingerprinted and retained. Normalized records are additive.

A global multi-caller collection must not rely on record IDs that are unique only inside one single-sample manifest. Lane ID participates in collection-level identity so equal caller-native row names from different specimens or replicates cannot collide.

For all callers:

- raw/native artifacts survive;
- alternative fits/solutions survive;
- exclusions change contribution state, not retention;
- normalized events link to exact source evidence;
- aggregate metrics retain exact numerator/denominator/exclusion membership.

## 9. Dependency-aware comparison

Agreement is shown, not voted.

Examples:

- QDNAseq 100/500/1000 kbp are resolutions of one family.
- Spectre with Sniffles support is not an independent confirmation of Sniffles.
- SAVANA-CNA may depend on SAVANA SV breakpoints.
- Wakhan may consume Severus breakpoints; this relationship must remain visible.
- different callers sharing the same BAM still share measurement noise and sampling limitations.

No generic "N of M callers = true" score is introduced.

## 10. Validation integration

The existing CNV validation registration, full-evidence, traceability and metric contracts remain authoritative for CNV studies. New caller inclusion creates a new prospectively registered study version; an already locked study is not rewritten after outcome review.

SV comparisons use SV-specific breakpoint/type/orientation matching rather than reusing CNV overlap thresholds blindly.

Missing caller-specific outputs map to NOT_APPLICABLE/NOT_EVALUABLE as appropriate; values are not invented.

## 11. Staged implementation

### 5A — catalog, eligibility, lanes and dependency graph

Implement typed catalog definitions, input-role declarations, lane identities, eligibility evaluation, deterministic content lock and cycle detection. Existing canonical runner remains unchanged.

### 5B — Spectre

Add a depth-only Spectre adapter first. Native output and parameters remain retained. No hidden Sniffles input.

### 5C — existing Sniffles2/cuteSV integration

Expose existing adapters through the multi-caller plan without changing their current normalized results. Mosaic remains a separate Sniffles mode.

### 5D — allele/SNP/phasing input contract

Add typed, fingerprinted auxiliary inputs used by tumor CNA callers.

### 5E — Severus

Qualify paired mode first. Preserve native breakpoint graphs, clusters and haplotype assignments. Single-sample mode remains separately gated.

### 5F — SAVANA

Keep paired and tumor-only modes separate. Preserve both SV and CNA native outputs and alternative purity/ploidy solutions.

### 5G — Wakhan

Add phased CNA lane with explicit upstream SNP/phasing and optional registered breakpoint dependency. Prevent dependency cycles.

### 5H — ichorCNA

Add a distinct ULP-WGS/cfDNA-compatible lane. Any ONT transfer study requires its own preregistration and is not automatically promoted.

### 5I — comparison report and Desktop surface

Add a dependency-aware comparison view over retained evidence. It displays conflicts, dependencies, failures and missingness. It does not auto-promote a caller, auto-release ISCN or make a clinical-validity claim.

## 12. TDD acceptance requirements

At minimum:

1. SAVANA paired + missing normal => INELIGIBLE with explicit reason; no tumor-only substitution.
2. Unknown provider/mode => plan validation failure before execution.
3. Two specimens with identical native row names => globally unique collection identities.
4. Missing runtime => RUNTIME_UNAVAILABLE, not NO_CALL.
5. Parent-lane change => dependent lock changes; unrelated lane lock does not.
6. Dependency cycle => fail before execution.
7. Spectre depth-only plan has no Sniffles parent.
8. Adaptive Sampling cannot satisfy lcWGS eligibility merely from BAM presence.
9. QDNAseq 100/500/1000 and ACE alternatives remain retained under the existing Block-4 contract.
10. Existing pipeline tests remain green.
11. Each new executable runtime later requires a real-tool smoke test on synthetic data.
12. No real patient/genomic payload enters Git.

## 13. Release boundary

This program is Research Use Only.

It does not select a winning caller, define universal clinical thresholds, infer clinical validity from software tests, modify package version by itself, merge automatically, or permit automatic clinical release. Human review remains mandatory.
