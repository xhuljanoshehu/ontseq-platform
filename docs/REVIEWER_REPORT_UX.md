# Reviewer-first report information architecture

Status: research-only product-design contract; not a clinical reporting specification.

## Purpose

The ONTSeq reviewer report should help a physician, molecular pathologist, cytogeneticist, or bioinformatician answer four questions in order:

1. Was the sample technically analyzable?
2. What genomic evidence was observed?
3. How observable was each relevant locus, and where are the blind spots?
4. What remains unvalidated or requires expert review?

The interface must not collapse these questions into a single binary result.

## Design principles

### Evidence before interpretation

Display source evidence separately from interpretation. A Sniffles2 BND, a CNV segment, a gene overlap, and an ISCN proposal are distinct objects with distinct evidence levels.

### NO_CALL is visually distinct from negative

`NO_CALL`, `FAILED`, `NOT_RUN`, `COMPLETED with no event`, and a future validated negative state must never share the same badge, icon, or summary sentence.

### Observability is first-class

For targeted/adaptive assays, absence of a call is interpretable only in the context of locus observability. Both fusion breakpoints should therefore display an observability state and reason.

### Validation status is always visible

Research-only/unvalidated status must be visible at the page level and on any interpretation that could otherwise look clinical. Validation cannot be hidden in a footer.

### Progressive disclosure

The first screen should support rapid review. Raw technical detail should remain available one level deeper without overwhelming the clinical summary.

## Proposed report hierarchy

### 0. Persistent report banner

Always visible at top:

- sample pseudonym / run identifier
- assay mode
- genome build/reference lock
- analysis profile/version
- overall pipeline status
- validation state: `RESEARCH ONLY / UNVALIDATED`
- generated timestamp

If any critical module is `FAILED` or `NO_CALL`, the banner should surface that before genomic findings.

### 1. Review summary

Compact cards, ordered by clinical review flow rather than pipeline execution order:

- Technical analyzability
- CNV evidence
- Structural-variant/fusion evidence
- Observability / assay limitations
- ISCN proposal status
- Review status

Each card uses one of a small number of semantic states:

- Evidence observed
- No evidence observed in assessable scope
- Limited assessability
- NO_CALL
- Failed
- Not run
- Review required

No card should say simply `negative` until a module has a validated negative-result contract.

### 2. Quality and assay context

Show only reviewer-relevant metrics by default:

- aligned fraction
- mapped yield
- read-length summary
- mean/median genome-wide coverage where applicable
- target-coverage summary for Adaptive Sampling
- genome build/reference
- target BED identity/version and role

Detailed histograms and tool-specific metrics are expandable.

### 3. CNV review

Primary view:

- chromosome-level copy-number summary
- genome-wide copy-number plot
- cytoband/segment table
- alternative caller/consensus evidence when available
- cellularity/ploidy estimate with uncertainty/alternative solutions

Every CNV row should preserve:

- genomic interval
- cytoband mapping
- copy-number estimate
- caller/source evidence
- observability/coverage context
- validation/review status

Do not convert a caller output directly into a clinical karyotype assertion.

### 4. Fusion / structural-variant review

Primary fusion table columns:

- Candidate
- Breakpoint A
- Breakpoint B
- Gene overlap A
- Gene overlap B
- Support reads / caller evidence
- Genomic junction orientation
- Observability A
- Observability B
- Redundancy flag
- Evidence-source agreement
- Review state

Important interaction rules:

- `5′`/`3′` transcript direction is blank unless independently resolved.
- Genomic BND orientation is labelled explicitly as genomic adjacency, not transcript direction.
- Missing BND descriptor is shown as `orientation unavailable`, not inferred.
- Duplicate/reciprocal BND records are flagged as potentially redundant, never silently collapsed.
- Gene-intergenic and intergenic-intergenic candidates remain visible when scientifically useful.

A Circos-style overview may be offered as a secondary visualization, but the table remains the authoritative reviewer interface because it can expose uncertainty and provenance.

### 5. Observability map

Dedicated section answering: `Where could this assay reasonably have detected an event?`

For Adaptive Sampling this should distinguish:

- analysis ROI
- operational selection/buffered BED if known
- observed target coverage
- limited/not-assessable regions
- unknown state

Technical depth thresholds must not be presented as clinical adequacy thresholds unless analytically validated.

### 6. ISCN proposal

Display ISCN only as a proposal until the corresponding transformation has been validated.

Required labels:

- ISCN edition
- conformance profile
- source event IDs
- uncertainty markers
- review status

The UI should allow the reviewer to trace every ISCN element back to the evidence that generated it.

### 7. Evidence and provenance drawer

Expandable technical layer:

- tool name/version
- policy/profile identifier
- reference/build lock
- file fingerprints where appropriate
- caller parameters
- rejected-record counts/reasons
- module limitations
- schema version

No raw read names, raw inserted sequence, patient identifiers, or local filesystem paths should be exposed in reviewer artifacts.

## Fusion candidate detail view

A selected candidate should open a detail panel with four visually separated blocks:

1. `Observed evidence`
   - source caller(s)
   - support
   - exact breakpoints
   - BND form / genomic junction orientation

2. `Annotation`
   - overlapping genes/transcripts
   - distance to genes when flank annotation is used
   - annotation resource/version/checksum

3. `Assay observability`
   - status per breakpoint
   - coverage context
   - reason for limited/not-assessable state

4. `Interpretation boundary`
   - research-only
   - transcript orientation unresolved/resolved
   - expression evidence absent/present
   - orthogonal validation absent/present
   - reportability false until validated workflow permits otherwise

## Suggested visual semantics

Avoid red/green-only semantics. Use text + icon + shape/state label.

Recommended conceptual badges:

- `EVIDENCE` — observed source evidence
- `LIMITED` — technically limited observability
- `NO_CALL` — analysis could not produce an interpretable call
- `FAILED` — technical execution failure
- `REVIEW` — human assessment required
- `UNVALIDATED` — research-only method/result

The exact colors are implementation details and should be checked for accessibility/contrast.

## Product acceptance tests

The report design is not acceptable if any of the following are possible:

- a reviewer can confuse `NO_CALL` with `negative`
- a genomic BND orientation appears as transcript 5′/3′ direction
- a candidate looks clinically reportable despite `research_only=true`
- an off-target/poorly covered partner is hidden behind a simple absent-call statement
- duplicate or reciprocal BND records are counted as independent events without a warning
- an ISCN string cannot be traced to its source evidence
- raw genomic identifiers/sequence or local paths leak into a portable reviewer artifact

## Relationship to the legacy thesis report

The previous thesis report concept is useful as a functional baseline: QC cards, copy-number summary, ISCN text, ideogram/ideoplot, fusion table, and Circos visualization are all valuable reviewer components. ONTSeq extends that concept by making uncertainty, module state, observability, provenance, validation status, and human-review boundaries explicit first-class elements rather than implicit caveats.
