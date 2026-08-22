# ONTSeq interactive review workspace

**Status:** design contract / implementation proposal  
**Scope:** Research Use Only (RUO), reviewer-facing visualization and evidence navigation  
**Non-goal:** this document does not promote any caller, threshold, fusion, CNV, ISCN or other result to clinical reportability.

## 1. Purpose

The current portable `report.html` is deliberately conservative: it renders module status, a proposed ISCN string, aggregate QC metrics, normalized genomic events, warnings and provenance. That is a safe baseline, but it is not yet an evidence-review workspace.

The next visualization layer should let a physician, cytogeneticist, molecular pathologist or bioinformatician answer, in this order:

1. **Was the sample technically analyzable?**
2. **Which modules ran, failed, produced `NO_CALL`, or were not applicable?**
3. **What evidence was observed and where in the genome?**
4. **How observable were the relevant regions, especially under Adaptive Sampling?**
5. **What alternative CNV fits or caller-level evidence could materially change interpretation?**
6. **Which conclusions are still research-only, unresolved or expert-review dependent?**
7. **Can every displayed interpretation be traced back to normalized evidence and provenance?**

The interface MUST optimize for review and traceability, not for visual novelty.

## 2. Product split: portable report versus local evidence workspace

A single HTML artifact should not be forced to do everything.

### 2.1 Portable `report.html`

Purpose: immutable, archiveable, self-contained reviewer snapshot.

Properties:

- generated only from normalized/versioned result contracts;
- no remote network requests;
- no patient identifiers or raw read names;
- no dependency on a running service;
- essential values visible without hover;
- interactive filtering/selection may be included only when the interaction is fully local and does not alter scientific state;
- static/export fallback remains readable if JavaScript is unavailable;
- large raw BAM/VCF files are not embedded.

### 2.2 Local ONTSeq Review Workspace

Purpose: evidence exploration against the immutable run envelope.

Properties:

- opened from the Windows desktop shell after analysis;
- backed by the local ONTSeq service only;
- may load local BAM/BAI, VCF, BED, CNV segments and reference tracks through controlled loopback endpoints;
- supports synchronized genome navigation and read-level review;
- no cloud/CDN dependency in production;
- every selection is non-destructive; algorithmic output remains immutable;
- reviewer annotations, when implemented later, are a separate review layer and never overwrite source evidence.

This split preserves the existing self-contained-report requirement while allowing a much richer local workspace.

## 3. Information architecture

The main screen is an **evidence workspace**, not a grid of equal-weight cards.

```text
┌──────────────────────────────────────────────────────────────────────────────┐
│ RUO · Sample · Run · Assay · Build · Reference lock · Analysis version      │
│ Pipeline state · critical NO_CALL/FAILED states · generated timestamp        │
├──────────────┬───────────────────────────────────────────┬───────────────────┤
│ Review rail  │ Main evidence viewport                    │ Evidence inspector│
│              │                                           │                   │
│ Overview     │ synchronized genome / CNV / SV /          │ selected object   │
│ QC           │ observability / chromosome view           │ coordinates       │
│ Coverage     │                                           │ evidence          │
│ CNV          │                                           │ uncertainty       │
│ SV/Fusion    │                                           │ provenance        │
│ ISCN         │                                           │ limitations       │
│ Provenance   │                                           │ review boundary   │
├──────────────┴───────────────────────────────────────────┴───────────────────┤
│ Active reference · ROI role/version · result schema · toolchain provenance  │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 3.1 Persistent header

Always visible:

- `RESEARCH USE ONLY / NOT CLINICALLY VALIDATED`;
- sample pseudonym and run ID;
- assay mode (`lcwgs` or `adaptive_sampling`);
- genome build and exact reference identity;
- analysis profile/pipeline version;
- overall release state;
- critical module-state summary;
- generated timestamp.

The header MUST surface `FAILED` and `NO_CALL` before showing genomic findings.

### 3.2 Review rail

The rail is task-oriented rather than pipeline-execution-oriented:

1. Overview
2. Quality / assay context
3. Coverage & observability
4. CNV
5. SV / fusion evidence
6. ISCN proposal
7. Provenance / methods

Future modules such as methylation appear only if the result contract contains an applicable module. They are not shown as empty decorative navigation.

### 3.3 Main evidence viewport

One dominant view at a time. Supporting views remain synchronized but visually subordinate.

The user should not need to mentally reconcile different chromosome selections between sections. A shared genomic selection state is recommended:

```text
selectedEventId
selectedChromosome
selectedInterval
selectedBreakpointPair
selectedCnvFit
selectedCoverageRegion
activeFilters
```

Selection affects only presentation, never normalized results.

### 3.4 Evidence inspector

A selected CNV segment, SV, target region or ISCN element opens the same inspector framework:

- object type and status;
- genomic locus/loci;
- cytoband(s), if build-aware mapping is available;
- copy number / length / event type where applicable;
- source caller(s) and version(s);
- support reads, VAF, local coverage, quality, precision, strand/orientation fields when available;
- annotations with source release/checksum and scope caveats;
- reportability (`false` must remain explicit);
- limitations and warnings;
- links to local read-level evidence when available;
- source event IDs and provenance.

Missing values are rendered as **not available**, not inferred.

## 4. Non-negotiable scientific visual semantics

### 4.1 Module-state contract

These states MUST remain visually and textually distinct:

| State | Meaning | Allowed visual language |
| --- | --- | --- |
| `COMPLETED` + evidence | analysis completed; evidence observed | evidence label + event count |
| `COMPLETED` + no evidence, where supported | analysis completed; no event in assessable scope | descriptive wording only; never a generic green “negative” unless a validated negative contract exists |
| `NO_CALL` | no interpretable call could be produced | dedicated NO_CALL state and reason |
| `FAILED` | execution attempted and failed | dedicated technical failure state |
| `NOT_RUN` | module not run/not applicable | neutral not-run state |
| `REVIEW_REQUIRED` | human assessment required | separate review marker |
| `UNVALIDATED` | method/result remains research-only | visible at page and finding level when clinically confusable |

No color-only encoding. State is always represented by text plus icon/shape/pattern.

### 4.2 Evidence is not interpretation

Keep distinct visual objects for:

- raw/normalized caller evidence;
- genomic event object;
- gene overlap/annotation;
- fusion hypothesis;
- ISCN proposal element;
- reviewer interpretation.

A Sniffles2 BND must not visually become a confirmed fusion merely because two genes are nearby. A database assertion must not appear as a classification of the sample finding. `reportable=false` must not be visually hidden.

### 4.3 Observability before negative inference

For Adaptive Sampling, absence of a call can only be understood beside observability. The UI should therefore show:

- analysis ROI identity and role (`analysis_roi_unbuffered`);
- per-region mean depth;
- fractions at the configured technical bins (currently 1x/10x/20x/30x when present);
- limited/unknown observability state when available;
- both breakpoints separately for a paired event;
- a statement that technical depth bins are descriptive until validated adequacy thresholds exist.

### 4.4 Uncertainty is first-class

Do not collapse CNV model selection to one number when alternatives are available. Display:

- primary cellularity, ploidy and fit error;
- candidate count;
- alternative cellularity/ploidy/error solutions;
- bin-size dependence;
- chromosome-level agreement across bins;
- explicit primary-bin designation.

The selected model may be visually emphasized, but alternatives remain inspectable.

## 5. Visualization inventory

### 5.1 Overview: “Can I review this run?”

Primary artifacts:

- module-state strip;
- concise QC summary;
- assay/reference context;
- evidence counts by event class;
- critical warnings/limitations;
- compact chromosome alteration summary if CNV exists;
- compact SV breakpoint overview if paired SV evidence exists.

Avoid a KPI-card wall. The normal state should be readable as a linear review narrative.

### 5.2 QC

Current Cramino normalized metrics support:

- number of reads;
- number of alignments;
- aligned percentage;
- total yield (Gb);
- mean coverage;
- N50 and N75;
- median/mean read length;
- median/mean/modal identity when returned by Cramino.

Recommended representation:

- exact values in a compact table with in-cell bars/sparklines only when a meaningful comparison exists;
- one read-length distribution only when histogram data is actually retained in the run envelope;
- no fake gauge thresholds while the QC policy is `technical_defaults_only`;
- failed validated gates, if later present, drawn as threshold annotations on the relevant chart.

### 5.3 Target coverage / Adaptive Sampling

Current normalized `TargetCoverageReport` supports exact per-region:

- chromosome/start/end/region ID;
- mean depth;
- base counts and fractions at each configured technical threshold;
- BED version, role and checksum;
- weighted mean, min/median/max regional mean depth;
- aggregate interval fractions at each threshold;
- overlap count and limitations.

Recommended views:

1. **ROI coverage matrix** — rows = target regions; columns = mean depth and threshold fractions.
2. **Sorted dot plot** — regional mean depth distribution, with median and range directly labelled.
3. **Genome-position target track** — target intervals placed along chromosomes.
4. **Threshold profile** — stacked or aligned bars for fractions ≥1x/10x/20x/30x; labels must say “technical descriptive bins”.
5. **Coverage-region inspector** — exact interval, BED fingerprint and warning state.

Future off-target viability metrics such as CV/MAD, GC residuals, mappability dependence and autocorrelation must appear only when corresponding normalized contracts exist.

### 5.4 CNV

Current QDNAseq+ACE result data support:

- 100/500/1000 kbp fits when configured;
- explicit primary bin size;
- cellularity, ploidy, fit error, candidate count and segment count;
- alternative fit solutions;
- segment and chromosome files;
- per-chromosome multi-bin consensus;
- normalized gain/loss/deletion/duplication events with copy number;
- primary and fit/copy-number plots already emitted by the R lane.

Recommended views:

1. **Whole-genome copy-number track** — chromosome concatenation or faceted chromosomes, with segment levels and direct labels.
2. **Chromosome ideogram** — G-bands/cytobands with gain/loss overlays when build-aware cytobands are available.
3. **Fit explorer** — cellularity × fit-error profile with the selected fit and alternatives; ploidy encoded explicitly.
4. **Bin-size comparison** — 100/500/1000 kbp small multiples; never overlay three noisy profiles without clear differentiation.
5. **Chromosome consensus strip** — median CN, min–max range and agreeing/contributing bin counts.
6. **Segment table linked to genome** — selecting a row highlights the same genomic interval in every view.

Do not infer cytobands if no build-aware mapping was performed.

### 5.5 SV and fusion evidence

Normalized event/evidence fields can support:

- DEL/DUP/INV/INS/translocation/BND-like paired event presentation;
- one or two loci;
- event length;
- caller/version;
- support reads;
- local coverage;
- VAF;
- quality and filters;
- strand/orientation evidence;
- coverage context;
- alignment mismatch summary;
- positional/length standard deviation;
- precision flag;
- gene overlap and locked knowledge-resource annotations.

Recommended views:

1. **SV inspector table** as authoritative triage surface.
2. **Circular whole-genome breakpoint overview** for inter-chromosomal relationships and global structure.
3. **Breakpoint split view** for the selected paired event.
4. **Evidence inspector** preserving every caller independently.
5. **Local read-level genome browser** for BAM/VCF evidence; expert IGV/JBrowse review is evidence review, not truth.

Fusion-specific fields or transcript direction are displayed only if a future fusion contract actually resolves them. Genomic BND orientation must not be relabelled as transcript 5′/3′ direction.

### 5.6 ISCN proposal

ISCN is rendered as an evidence-linked proposal only.

Recommended interaction:

- each renderable proposal element has a stable source-event mapping;
- selecting a proposal element highlights the source CNV/SV event and locus;
- the inspector shows source event IDs, edition, conformance profile, uncertainty and review status;
- unsupported constructs remain visibly unsupported rather than approximated.

### 5.7 Provenance

Two complementary views:

1. **human-readable stage graph** — intake → QC → coverage → CNV/SV → assembly → reports/release, with stage state and runtime where available;
2. **exact provenance table** — tool/version/parameters, policy/profile IDs, reference/BED locks, schema versions and fingerprints.

The graph is explanatory. The table/JSON remains authoritative.

## 6. Recommended visualization stack

The stack should be modular; no single library needs to own every layer.

### 6.1 JBrowse 2 — preferred read-level/SV evidence browser

Use for the local Review Workspace, not as the portable report renderer.

Reasons:

- Apache-2.0 open source;
- BAM/CRAM, VCF, BED, BigWig and common genomic formats;
- linear, circular, breakpoint split and SV-inspector views;
- purpose-built structural-variant workflows;
- embeddable React components;
- local/controlled data can remain on-premises.

Candidate packages:

- `@jbrowse/react-linear-genome-view2` for focused evidence;
- `@jbrowse/react-circular-genome-view2` for paired SV overview;
- full `@jbrowse/react-app2` only if the product actually needs its complete workspace.

### 6.2 Gosling.js — preferred custom linked genomic summary layer

Use when ONTSeq needs a branded, declarative genome-scale visualization rather than a full genome browser.

Strengths:

- genomic linear/circular layouts;
- BAM/VCF/BED/JSON/CSV support;
- zoom, pan, brushing and linked views;
- semantic zoom and responsive visualization;
- concise declarative specs are easier to test than hand-built SVG for many track layouts.

Best candidates: CNV/coverage overview tracks and synchronized summary views.

### 6.3 Ideogram.js — preferred cytogenetic chromosome rendering candidate

Use for human chromosome/cytoband views where the chromosome itself is the domain surface.

Strengths:

- GRCh37/GRCh38 assembly selection;
- banded chromosomes;
- overlay/track/heatmap/histogram annotation layouts;
- brush/cursor interaction;
- compact genome-wide chromosome presentation.

Before production use, bundle the exact required band assets/version locally and verify build identity against ONTSeq reference locks.

### 6.4 Apache ECharts — preferred general statistical/QC chart engine

Use for non-genomic statistical charts.

Strengths:

- zoom, pan, brush, data filtering and export;
- Canvas or SVG renderer;
- accessible ARIA descriptions and decal patterns;
- good fit with the historical ezCharts ecosystem while allowing a modern custom visual system.

Use SVG for moderate-size export-oriented charts; Canvas for dense mark fields after profiling.

### 6.5 D3/SVG — custom glue, not default charting engine

Use only where ONTSeq-specific geometry is clearer than declarative libraries, for example:

- ISCN-element-to-chromosome trace links;
- bespoke breakpoint orientation diagrams;
- compact provenance/DAG annotations;
- specialized label placement.

React/WPF owns application state; D3 owns geometry, scales and SVG behavior inside a bounded component.

### 6.6 Legacy ezCharts

Keep as historical/reference compatibility, not as the sole design constraint. The thesis report benefited from ezCharts’ Nanopore-oriented report components and self-contained HTML generation. ONTSeq can preserve useful visual patterns without inheriting every component or layout decision.

## 7. Offline and security architecture

Production reviewer visualization MUST NOT depend on runtime CDN fetches.

Preferred packaging:

```text
ONTSeq runtime
  ├── normalized result JSON
  ├── local report assets / vendored JS bundles
  ├── local review-web bundle
  └── loopback-only evidence endpoints
          ├── reference sequence/index
          ├── BAM/BAI range access
          ├── VCF/index
          ├── BED/ROI
          └── CNV segment/coverage tracks

WPF desktop
  ├── analysis launcher/status
  ├── open portable report
  └── open local Review Workspace
```

If the Review Workspace is hosted in-process in WPF later, WebView2 is one possible shell, but adding it is a separate dependency/packaging decision. The current desktop project does not yet declare a WebView dependency, so it must not be assumed to exist.

Security requirements:

- loopback bind only unless an explicitly reviewed remote mode exists;
- unpredictable per-session token for local evidence URLs if a browser surface can access them;
- path allowlist rooted in the immutable run envelope;
- no arbitrary file browsing through the local HTTP layer;
- no raw read names copied into portable reports;
- no remote fonts, analytics, telemetry or third-party resources;
- Content Security Policy for the report/review bundle;
- version/checksum every vendored visualization dependency.

## 8. Responsive design contract

### Desktop / workstation

Primary intended reviewer surface. Target layout:

- persistent header: 56–72 px;
- review rail: approximately 180–220 px;
- main evidence viewport: largest flexible area;
- inspector: approximately 320–420 px, collapsible;
- tables can use sticky headers and column pinning where useful.

### Narrow/mobile portrait

Mobile is a review/triage fallback, not the primary read-level interpretation surface.

Order:

1. RUO + critical module state;
2. summary and warnings;
3. main selected visualization;
4. tabbed sections;
5. inspector as bottom sheet/drawer;
6. detailed provenance last.

No hover-only evidence. Dense genome views must provide tap/focus selection, larger hit regions, explicit zoom/reset controls and a text/table alternative.

## 9. Accessibility contract

Every meaningful visualization needs a non-visual path:

- direct labels for essential values;
- associated data table for charts where exact values matter;
- textual state summaries;
- keyboard-accessible selection/reset/filter controls;
- focus styles;
- color + shape/pattern/text redundancy;
- no red/green-only status system;
- `prefers-reduced-motion` support;
- accessible export/static representation.

The portable report should remain understandable as a screenshot/PDF even without tooltips.

## 10. Visual design direction

The target is a professional scientific instrument, not a generic SaaS dashboard.

Principles:

- white/neutral evidence canvas;
- restrained typography with strong numerical alignment;
- one primary interaction accent;
- event-class colors separate from UI interaction colors;
- low-saturation context and grid lines;
- no gradients/glow/3D unless the encoding itself requires them;
- direct labels over detached legends where possible;
- meaningful whitespace rather than nested card stacks;
- cytogenetic/genomic structures provide the domain-specific visual identity.

Suggested semantic roles (exact colors require contrast testing):

- neutral/reference;
- deletion/loss;
- duplication/gain;
- inversion;
- translocation/paired breakend;
- selection/focus;
- warning/limited;
- technical failure;
- unvalidated/review required.

The same hue must not simultaneously mean “selected” and “biological event class”.

## 11. Interaction contract

### Default

- no item selected;
- whole-genome or assay overview visible;
- critical warnings visible;
- essential values visible without hover.

### Preview

- mouse hover may preview a mark on desktop;
- hover never contains the only copy of essential information.

### Selected

- click/tap commits selection;
- selection persists while switching compatible views;
- inspector opens with exact data and provenance;
- linked views highlight the same event/locus.

### Expanded evidence

- explicit action opens read-level view or larger plot;
- close/back returns to the previous selected context.

### Reset

- one visible reset returns genomic zoom, filters and selection to the default review state.

### Share/export

Portable report state may be exported; local BAM-backed evidence URLs are not portable and must not be serialized into shareable external links.

## 12. QA gates

Visualization work is not complete until these tests pass.

### Scientific semantics

- `NO_CALL` cannot be mistaken for negative;
- `FAILED` cannot be mistaken for not run;
- technical 1x/10x/20x/30x coverage bins cannot read as validated adequacy thresholds;
- BND orientation cannot read as transcript 5′/3′ direction;
- `reportable=false` remains obvious;
- ISCN proposal always traces to source events;
- missing data is not fabricated;
- GRCh37/GRCh38 tracks use the matching locked build.

### Visual correctness

- exact axis units and coordinate conventions are shown;
- chromosome ordering is canonical 1–22, X, Y;
- copy-number baseline/ploidy is labelled, not implied;
- alternative CNV fits remain available;
- color-deficiency/grayscale review passes;
- static export keeps the same conclusion as the interactive view.

### Interaction

- keyboard selection/filter/reset;
- touch selection without hover;
- reset works from every state;
- zoom state is synchronized only where mathematically valid;
- tables and plots cross-highlight the same object ID;
- failed/missing track loading degrades visibly and does not blank the rest of the report.

### Performance

- do not instantiate a genome browser for every event row;
- one focal heavy genomic viewer at a time;
- use lazy initialization for read-level evidence;
- use Canvas only when data density justifies it;
- profile representative large BAM/VCF and dense CNV tracks locally before locking renderer choices.

## 13. Implementation sequence

### Phase V0 — design approval

- create large-screen and mobile visual concepts using only synthetic/fixture data;
- review with physician + bioinformatics/cytogenetics stakeholders;
- lock layout hierarchy, state semantics and visual tokens before product code.

### Phase V1 — portable report modernization

No biological logic changes.

- refactor `report.py` into view-model + renderer layers;
- preserve current tables and warnings;
- add module-state strip and evidence-linked event selection;
- add ECharts/D3 only through vendored pinned assets;
- render existing CNV artifacts and fit alternatives safely;
- add target-coverage plots from `TargetCoverageReport` when present;
- add print/PDF/static fallback;
- snapshot/DOM/accessibility tests.

### Phase V2 — linked genomic summary

- canonical genomic coordinate state;
- chromosome ideogram;
- linked CNV segments/coverage targets/SV marks;
- synchronized inspector;
- build-aware reference assets.

### Phase V3 — local read-level evidence

- add local evidence service endpoints with strict path scoping;
- integrate JBrowse 2 (preferred) or IGV.js after a focused technical spike;
- deep-link selected event → breakpoint/segment locus;
- load BAM/BAI/VCF only on demand;
- no changes to source evidence.

### Phase V4 — review/audit layer

- reviewer comments and disposition as separate immutable layer;
- explicit original algorithmic result versus reviewer interpretation;
- audit trail and release workflow only after governance design.

### Phase V5 — validation workspace

Separate from single-sample review:

- sensitivity/specificity/precision/no-call by predefined strata;
- LoD plots constrained to tested ranges;
- breakpoint-error distributions;
- truth resolution and not-assessable regions;
- specimen-aware uncertainty/intervals;
- caller comparisons and reproducibility views.

## 14. Acceptance criterion

The visualization program succeeds when a reviewer can move from a high-level run state to the exact evidence underlying any displayed CNV/SV/coverage/ISCN proposal without losing assay context, provenance or uncertainty — and without the interface ever turning an unvalidated technical observation into a clinical claim.