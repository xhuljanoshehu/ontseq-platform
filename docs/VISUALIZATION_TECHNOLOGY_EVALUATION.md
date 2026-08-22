# Visualization technology evaluation

**Status:** engineering research note  
**Decision scope:** visualization/rendering only; no scientific caller or reportability decisions.

## 1. Constraints specific to ONTSeq

A suitable visualization stack must satisfy more than “can draw a chart”.

Required:

- local/on-premises operation;
- GRCh37 and GRCh38 awareness;
- large genomic coordinate ranges;
- CNV, SV/BND, target intervals and later modified-base/methylation compatibility;
- read-level BAM/BAI evidence on demand;
- no cloud upload or runtime analytics;
- portable HTML snapshot plus richer local workspace;
- explicit `NO_CALL`/`FAILED`/`NOT_RUN`/unvalidated semantics;
- accessible non-hover paths;
- deterministic versioning and vendoring;
- static/export fallback;
- compatibility with a Windows WPF operator shell and local Linux/WSL analysis service.

The current WPF project has no browser/WebView package dependency in its project file. Therefore embedding a web review surface is a future engineering decision, not an existing capability.

## 2. Candidate matrix

| Technology | Best ONTSeq role | Genomic formats/views | Interaction | License | Primary risk | Recommendation |
| --- | --- | --- | --- | --- | --- | --- |
| **JBrowse 2** | local read-level/SV evidence workspace | BAM/CRAM, VCF, BED, BigWig; linear, circular, breakpoint split, SV inspector | mature navigation/filtering, linked views | Apache-2.0 | bundle/runtime integration and local file serving complexity | **preferred local evidence browser** |
| **IGV.js** | focused local BAM/VCF evidence | mature alignment, variant, copy-number, multi-locus views | mature genome-browser navigation | MIT | less native whole-workspace SV composition than JBrowse 2 | **mandatory comparator in technical spike** |
| **Gosling.js** | branded linked genomic summaries | BAM/VCF/BED/JSON/CSV; linear/circular genomic grammar | zoom, pan, brushing, linked views, semantic zoom | MIT | separate engine adds bundle/maintenance surface | **strong candidate for CNV/coverage summaries** |
| **Ideogram.js** | chromosome/cytoband domain surface | chromosome bands + annotations/heatmap/histogram/brush | chromosome-level interaction | CC0-1.0 repository license | exact band-asset provenance/build lock must be controlled | **strong candidate for karyotype-style view** |
| **Apache ECharts** | QC, coverage distributions, fit/validation plots | general tabular/statistical data | dataZoom, brush, filters, toolbox; SVG/Canvas | Apache-2.0 project | not genome-coordinate aware by itself | **preferred general chart engine** |
| **D3** | bespoke SVG glue | any normalized data after custom modeling | fully custom | ISC | high implementation/QA cost if overused | **use selectively** |
| **ezCharts** | historical report compatibility | Python wrappers around ECharts/Bokeh plus Nanopore report components | depends on backend | repository dependency review required before vendoring | older visual system / Python-side coupling | **reference, not sole future framework** |

## 3. JBrowse 2 assessment

Official JBrowse 2 documentation describes:

- linear genome view;
- circular view with BND/translocation chords;
- breakpoint split view;
- an SV inspector combining a filterable table and circular overview;
- BAM/CRAM, VCF, BED, BigWig and other genomic formats;
- embeddable React components;
- SVG export;
- a dedicated cancer/SV workflow and methylation tutorial.

### ONTSeq fit

Very strong for:

- selected SV → two breakpoint loci;
- BAM evidence at breakpoints;
- selected CNV → local alignment/coverage context;
- VCF event inspection;
- future modified-base/methylation evidence;
- whole-genome circular SV context.

### Packaging strategy

Do **not** use hosted public data services as part of a production run. Bundle exact NPM artifacts and configure references/tracks from local controlled resources.

A local evidence service can provide byte-range-capable endpoints inside the current run envelope. Only allowed run files should be addressable.

### Preferred integration level

Start with embedded focused components rather than the full application:

1. linear view for a selected locus;
2. circular view for SV overview;
3. full `react-app2` only if the smaller components cannot meet the review workflow.

This minimizes duplicated UI and keeps ONTSeq in control of status, provenance and interpretation boundaries.

## 4. IGV.js assessment

IGV.js is a mature embeddable genome viewer with BAM/VCF and copy-number examples, multi-locus navigation and a programmatic browser API.

### ONTSeq fit

Strong for:

- exact read-level review;
- rapid “open this event in evidence viewer” workflow;
- selected two-locus fusion/translocation inspection;
- ROI overlay.

### Why it remains a comparator rather than an automatic choice

ONTSeq’s target workspace needs whole-genome SV triage, paired-breakpoint workflows and potentially methylation in addition to ordinary alignment inspection. JBrowse 2 currently exposes more of those as integrated first-class view types. IGV.js may still win if its integration footprint and familiar reviewer behavior are materially better.

### Technical spike must compare

- offline bundle size;
- local BAM/BAI load latency;
- multi-gigabyte range access;
- two-breakpoint navigation;
- VCF/BND handling;
- reference/ROI configuration;
- keyboard/accessibility behavior;
- memory after repeated event selection;
- export/static evidence capture;
- packaging into WPF/local service.

## 5. Gosling.js assessment

Gosling is a declarative grammar designed specifically for genomics. Official documentation includes:

- CSV/GFF3/VCF/BED/JSON/BAM data;
- linear and circular layouts;
- zoom and pan;
- linked views;
- brushing;
- semantic zoom;
- responsive specifications.

### ONTSeq fit

Potentially ideal for a **custom scientific summary surface** where a full genome browser would provide too much chrome.

Candidate views:

- whole-genome CNV segment plot;
- target-coverage tracks;
- chromosome-linked overview/detail;
- event density / paired breakpoint overview;
- synchronized summary tracks above the evidence inspector.

### Boundary

Gosling is not a replacement for a read-level expert browser. It should visualize normalized/derived tracks and hand off to JBrowse/IGV for read evidence.

## 6. Ideogram.js assessment

Ideogram.js provides human banded chromosomes, assembly selection, annotation tracks, overlays, heatmaps, histograms and brush/cursor interactions.

### ONTSeq fit

Use the chromosome itself as the contextual substrate for:

- chromosome gain/loss;
- arm/focal CNV overlays;
- selected cytobands;
- location of structural events;
- target-panel distribution.

### Critical provenance requirement

The application must not silently use whatever band file a library happens to ship. The exact band geometry used for display and coordinate-to-band interpretation must be compatible with the ONTSeq locked reference/build. If an external library asset cannot be tied to a controlled build/version, ONTSeq should vendor/derive an approved cytoband dataset and render it itself.

## 7. Apache ECharts assessment

ECharts is well suited to the non-genome-specific analytical layer. Official documentation provides Canvas/SVG rendering, data zoom, brush interactions, export/toolbox functions and ARIA/decal accessibility features.

### ONTSeq fit

Recommended for:

- target-depth distributions;
- threshold-fraction matrices/bars;
- ACE fit alternatives;
- cellularity/ploidy exploration;
- bin-size comparisons;
- validation plots;
- runtime/memory distributions;
- compact tables with linked chart selection.

### Renderer rule

Prefer SVG when:

- the number of marks is moderate;
- crisp export/print is important;
- mobile memory is a concern.

Prefer Canvas only after profiling when dense marks or repeated redraw make SVG too expensive.

### Accessibility rule

Enable/configure ARIA intentionally and still provide an adjacent semantic table/text summary. Automated ARIA text is supportive, not sufficient for a complex genomics interpretation.

## 8. D3 assessment

D3 remains valuable for geometry that is genuinely ONTSeq-specific.

Use cases:

- genomic-adjacency orientation schematic;
- ISCN proposal-to-source-event trace connectors;
- custom chromosome/event glyphs if Ideogram cannot be safely build-locked;
- compact workflow/provenance graph;
- annotation connectors and label collision handling.

Avoid using D3 to rebuild standard histograms, bars and scatter plots already covered by ECharts.

## 9. Legacy thesis/ezCharts assessment

The thesis selected ezCharts because it provided Nanopore-oriented report components, ideoplot/ideogram support and self-contained HTML output. That remains useful design evidence.

ONTSeq should preserve the strengths:

- self-contained reviewer artifact;
- clear QC summary;
- ideogram/ideoplot idea;
- fusion/SV circular overview;
- readable tables.

But ONTSeq should extend the model with:

- module status semantics;
- observability;
- explicit uncertainty;
- source-event traceability;
- local read-level evidence;
- richer accessibility;
- synchronized selection.

The historical report is a baseline, not a technical specification.

## 10. Portable report rendering decision

The portable report should remain **small, local and deterministic**. Recommended initial architecture:

```text
PipelineResult + optional module artifacts
        ↓
ReviewerViewModel (pure Python/typed transformation)
        ↓
HTML renderer
   ├── semantic tables/text fallback
   ├── inline CSS
   ├── pinned/vendored visualization bundle
   └── serialized sanitized plot data
        ↓
self-contained report.html
```

Heavy genome-browser code is excluded from this artifact.

### Why a view-model layer is required

Current `report.py` reads `PipelineResult` directly and assembles HTML strings. A view model would:

- decide display order without altering scientific values;
- normalize units/labels;
- create explicit “not available” values;
- group evidence by event;
- expose plot-ready data separately from rendering;
- allow unit tests to validate semantics without parsing pixels.

## 11. Local Review Workspace decision

Recommended target architecture:

```text
WPF operator shell
       ↓ open review
local ONTSeq web UI bundle
       ↓
loopback review/evidence API
       ├── result contracts
       ├── CNV segments/fits
       ├── target coverage
       ├── VCF/index
       ├── BAM/BAI
       ├── controlled reference assets
       └── provenance
```

### Application framework

For the review web bundle, React + TypeScript is the leading candidate because JBrowse’s maintained embedded components are React packages and the workspace requires coordinated selection state. This decision does **not** require rewriting the existing WPF launcher.

WPF remains the desktop operator shell; the review web UI is a bounded evidence surface.

## 12. WPF embedding options

Evaluate separately:

### A. Open review workspace in system browser

Pros:

- lowest desktop dependency complexity;
- simple debugging;
- browser engine already maintained.

Cons:

- less integrated feel;
- must secure loopback session and lifecycle carefully.

### B. WebView2 inside WPF

Pros:

- unified application experience;
- direct navigation from run list to evidence.

Cons:

- new runtime/package dependency;
- packaging and environment checks;
- additional browser security surface;
- must test deployment across supported Windows configuration.

The first technical spike should support system-browser mode even if WebView2 becomes the final UX. This isolates scientific/evidence work from desktop-shell integration.

## 13. Dependency and supply-chain rules

Before any visualization library is promoted into the packaged application:

- pin exact version;
- record source repository and license;
- include required attribution/notice;
- vendor/build locally; no production CDN import;
- record content hashes in the distribution manifest/SBOM;
- run dependency/security review;
- test offline mode;
- verify no telemetry/network calls;
- validate that reference/annotation assets are version controlled independently from library code.

## 14. Recommended technical spikes

### Spike A — CNV/coverage summary

Build from synthetic fixtures only:

- ECharts fit explorer;
- Gosling whole-genome segment/coverage track;
- Ideogram cytoband overlay;
- linked selection by event ID and interval.

Pass criteria:

- GRCh37 and GRCh38 fixtures;
- zero/empty/no-call fixtures;
- responsive desktop/narrow layout;
- accessible table fallback;
- offline bundle.

### Spike B — read-level evidence

Load a synthetic/research-safe BAM/BAI and VCF into:

1. JBrowse 2;
2. IGV.js.

Compare the criteria in section 4 and make a recorded decision. Do not carry both dependencies indefinitely without a reason.

### Spike C — portable report

Render a complete synthetic `PipelineResult` with:

- module states;
- QC;
- target coverage;
- CNV fits/consensus;
- SV evidence;
- ISCN proposal;
- provenance.

Verify that disabling JavaScript leaves all essential values and warnings readable.

## 15. Current recommendation

**Portable report:** semantic HTML + vendored ECharts + limited custom SVG/D3; add chromosome rendering only after build-asset provenance is resolved.  
**Local evidence workspace:** React/TypeScript + JBrowse 2 as preferred genome evidence engine; IGV.js evaluated in a bounded comparator spike.  
**Genomic summary layer:** evaluate Gosling.js rather than hand-building every genomic track.  
**Cytogenetic surface:** evaluate Ideogram.js with a strict local/build-locked cytoband asset contract.  
**General rule:** keep renderers replaceable by feeding them typed ONTSeq view models rather than coupling scientific contracts to chart-library schemas.