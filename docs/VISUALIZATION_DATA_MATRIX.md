# ONTSeq visualization data matrix

**Status:** design/traceability artifact  
**Rule:** a visualization may only render a datum that is present in a versioned normalized contract or an explicitly identified generated artifact. Missing information is displayed as unavailable; it is never reconstructed from presentation files or inferred from visual context.

## 1. Data ownership hierarchy

1. **Normalized contracts / `result.json`** — authoritative structured reviewer data.
2. **Versioned module artifacts** — target coverage report, CNV call report, benchmark report, review record, etc.
3. **Raw evidence artifacts** — BAM/BAI, VCF/index, BED, generated segment tables and plots; opened only in the controlled local workspace.
4. **HTML/XLSX** — presentation/export surfaces, never parsed back into scientific state.

This implements the repository principle: typed adapters and structured contracts before presentation parsing.

## 2. Run and assay context

| Datum | Current source | Reviewer representation | Interaction | Do not imply |
| --- | --- | --- | --- | --- |
| sample ID | `SampleManifest.sample_id` | persistent header | copy/select only | direct patient identity |
| run ID | `SampleManifest.run_id` | persistent header / provenance | details | biological result |
| input kind | `InputSpec.kind` | provenance/context | details | identical validation across POD5/uBAM/aligned BAM |
| assay mode | `AssaySpec.mode` | prominent header | switches applicable sections only | lcWGS and Adaptive Sampling share QC assumptions |
| genome build | `AssaySpec.genome_build` | header + all genomic views | fixed per run | coordinates transferable between builds |
| reference ID | `AssaySpec.reference_id` + lock | header/provenance | details | reference identity from build name alone |
| target BED | manifest + coverage report | Adaptive Sampling context | open ROI view | selection BED equals analytical ROI unless explicitly locked as such |
| analysis profile | `AnalysisSpec.profile` | header/provenance | details | validated profile unless status says so |
| analysis intent | `AnalysisSpec.intent` | annotation context | details | somatic/germline scope if unset |

## 3. Module status

| Datum | Source | Visual | Required wording |
| --- | --- | --- | --- |
| module name | `ModuleOutcome.module` | module-state strip | exact module name |
| module status | `ModuleOutcome.status` | text + icon + pattern | COMPLETED / NO_CALL / FAILED / NOT_RUN |
| reason | `ModuleOutcome.reason` | inline critical reason / inspector | do not hide behind tooltip |
| QC verdict | `QCMetrics.verdict` | quality summary | distinguish descriptive WARN from validated PASS |
| release status | pipeline result | header | technical/release state only |

A module tile is never allowed to say “negative” solely because it has zero normalized events.

## 4. Cramino QC

Current normalized fields from `parse_cramino_json`:

| Field | Unit | Preferred visual | Secondary visual | Caveat |
| --- | --- | --- | --- | --- |
| `number_of_reads` | count | exact value | none | no quality judgement without policy |
| `number_of_alignments` | count | exact value | relation to reads if semantically valid | may exceed simple one-read/one-alignment assumptions |
| `aligned_percent` | % | exact value + restrained bar | policy threshold line only if configured | descriptive if no validated gate |
| `total_yield_gb` | Gb | exact value | comparison only with appropriate reference | no arbitrary green/red range |
| `mean_coverage_x` | × | exact value | genome-wide context | assay-specific meaning differs |
| `n50_bp` | bp | exact value | distribution marker if histogram exists | N50 is not median |
| `n75_bp` | bp | exact value | distribution marker if histogram exists | — |
| `median_length_bp` | bp | exact value | read-length distribution if retained | — |
| `mean_length_bp` | bp | exact value | read-length distribution if retained | — |
| `median_identity_percent` | % | exact value | identity distribution only if retained | may be estimated |
| `mean_identity_percent` | % | exact value | — | may be estimated |
| `modal_identity_percent` | % | exact value | — | may be estimated |
| `identity_estimated` | boolean-like text | explicit note | — | must remain visible with identity metrics |

### Current gap

The present canonical Cramino normalization stores aggregate metrics, not the full read-length histogram. Therefore the modernized report must not fabricate a read-length histogram from N50/mean/median alone. A histogram becomes available only if the pipeline intentionally retains and normalizes an appropriate distribution artifact.

## 5. Adaptive Sampling target coverage

Current `TargetCoverageReport` supports these displayable values.

### Per target region

| Field | Visual encoding | Interaction |
| --- | --- | --- |
| chromosome/start/end | genomic position | select/zoom |
| `region_id` | row/direct label | search/filter |
| `mean_depth` | dot/position + exact label | sort/select |
| `bases_at_threshold` | exact count in detail | details |
| `fraction_at_threshold` | aligned bars/heatmap cells | sort/filter |

### Summary

| Field | Preferred visual |
| --- | --- |
| region count | exact label |
| interval bases | exact label |
| interval-weighted mean depth | exact label + distribution context |
| minimum/median/maximum regional mean depth | dot/range summary |
| aggregate interval fraction at each threshold | threshold profile |
| overlapping interval count | warning/limitation callout |
| BED SHA-256/version/role | provenance drawer |
| tool version/policy | provenance drawer |

### Semantic guardrail

The currently configured 1×/10×/20×/30× bins are technical/descriptive. Visual labels should say, for example:

`Fraction of target bases ≥20× (descriptive technical bin)`

and not:

`20× coverage passed`.

## 6. CNV QDNAseq + ACE

### Fit-level data

`CnvFit` provides:

| Field | Visual |
| --- | --- |
| bin size (kbp) | small-multiple facet / selector |
| cellularity | fit plot axis + exact label |
| ploidy | fit point encoding + exact label |
| fit error | y-axis / sorted alternatives |
| candidate count | context label |
| segment count | context label |
| alternatives | explicit alternative fit table/points |
| `fit_plot` | existing artifact fallback/reference |
| `copy_number_plot` | existing artifact fallback/reference |

The new UI should prefer normalized numerical data where available. Existing generated plots remain useful as provenance-compatible fallback and side-by-side verification during migration.

### Chromosome consensus

`CnvChromosomeConsensus` provides:

| Field | Visual |
| --- | --- |
| chromosome | canonical 1–22, X, Y axis |
| median copy number | dot/segment position |
| rounded copy number | exact label only; distinguish from median |
| agreeing bins | numerator in agreement bar |
| contributing bins | denominator in agreement bar |
| min/max copy number | interval whisker |

Recommended single-chromosome glyph:

```text
chr5   1.82 ├────●────┤ 2.07     agreement 2/3
             min med max
```

Exact layout can differ; the key is to show both central estimate and cross-bin spread.

### Normalized CNV events

`GenomicEvent` adds:

- event ID;
- event class;
- exact primary interval;
- length;
- absolute copy number;
- caller evidence;
- confidence = currently often `unclassified`;
- `reportable=false`;
- notes containing bin size and ACE model context.

Visual requirements:

- event class color/pattern;
- copy number on a quantitative axis;
- event ID preserved for cross-linking;
- `unclassified` and `reportable=false` visible in inspector;
- no cytoband label unless a build-aware cytoband normalizer supplied it.

## 7. Structural variants

Current normalized `Evidence` supports richer review than the present HTML table exposes.

| Evidence field | Reviewer visual |
| --- | --- |
| caller + version | source chip/text |
| support reads | exact value + optional small bar when comparing callers |
| local coverage | exact value / context |
| VAF | exact fraction + position on 0–1 scale |
| quality | exact value with tool-specific label |
| filters | status text |
| supporting-read strands | explicit genomic strand notation |
| coverage context | mini-profile only when semantics/order are documented |
| mean alignment NM | exact technical evidence |
| position SD | breakpoint uncertainty interval/context |
| length SD | SV-length uncertainty context |
| `precise` | explicit precision state |

### Paired events

For translocation/fusion-type normalized events with `secondary` locus:

- circular chord = global orientation/context only;
- selected event opens two synchronized linear loci;
- both breakpoints get independent observability and annotation blocks when contracts provide them;
- the central connection diagram uses genomic adjacency, never inferred transcript direction.

### Single-locus events

DEL/DUP/INV/INS are generally clearer in linear/genome tracks and a table than in a Circos plot. The circular overview should not imply that every SV class is meaningfully represented as a chord.

## 8. Gene/knowledge annotations

`EventAnnotation` contains strong provenance/scope safeguards and the UI must preserve them.

| Field group | Display rule |
| --- | --- |
| source/release/SHA | always available in details |
| record ID/type | exact source reference |
| assertion | shown verbatim as source assertion |
| assertion vocabulary | adjacent to assertion |
| record origin | germline/somatic/unknown label |
| scope alignment | aligned/mismatched/unknown label |
| scope note | visible explanatory text |
| match type / reciprocal overlap | exact evidence of why record matched |
| review status/stars | source metadata only |
| genes/conditions | annotation details |
| caveats | mandatory visible section |

A source assertion such as germline “Pathogenic” must never be promoted into a somatic sample verdict by styling or wording.

## 9. ISCN proposal

Current proposal fields include notation, standard edition, conformance profile, review status and warnings.

### Portable report

- notation remains selectable text;
- edition/profile/review state immediately adjacent;
- warning if subset/unvalidated;
- source-event trace table when source mapping exists.

### Interactive workspace

Future requirement: structured proposal tokens rather than parsing the final notation string. Each token should carry:

```text
proposal_element_id
rendered_text
source_event_ids[]
source_loci[]
construct_type
supported_by_profile
uncertainty
review_status
```

Until such a contract exists, the UI must not pretend that arbitrary fragments of the ISCN string are safely clickable.

## 10. Provenance and integrity

Displayable from the current run/result contracts:

- pipeline version;
- tool names/versions/parameters;
- reference lock and FAI checksum;
- target BED version/checksum when applicable;
- normalized schema versions;
- module warnings/limitations;
- release/checksum state where available.

### Visualization

A pipeline DAG is explanatory and can show stage status, but exact provenance remains a table/JSON view. Never encode only by node color.

## 11. Validation workspace — separate data product

Single-sample review and method validation use different analytical units and must not be mixed into one dashboard.

Planned validation visualizations should be backed by benchmark contracts that preserve:

- truth/query events;
- event class;
- reciprocal overlap;
- breakpoint error;
- copy-number error;
- DETECTED/MISSED/NOT_ASSESSABLE/NO_CALL semantics;
- truth resolution;
- coverage and tumor/blast-fraction strata;
- specimen clustering/repeats;
- method/tool/version.

Recommended validation artifacts:

| Question | Visualization |
| --- | --- |
| sensitivity across event size | binned dot/interval plot |
| performance across coverage/purity | small-multiple heatmap or interval plot |
| breakpoint accuracy | error distribution / ECDF |
| copy-number agreement | truth vs estimate + difference plot |
| LoD | detection probability with observed levels and interval; never extrapolate past tested range |
| caller discordance | paired specimen/event comparison with explicit not-assessable category |
| reproducibility | paired replicate difference / concordance plot |
| failure/no-call burden | stacked state distribution by stratum |
| runtime/memory | distributions by tool and stage |

A single “accuracy” gauge is prohibited because the repository validation plan explicitly requires stratification.

## 12. Renderer ownership

| Layer | Preferred owner | Why |
| --- | --- | --- |
| application state/layout | web app / WPF shell | predictable product state |
| general QC/statistical charts | Apache ECharts | interaction + SVG/Canvas + accessibility |
| linked genomic summary | Gosling.js or custom D3 after spike | genomic coordinate semantics |
| chromosome cytobands | Ideogram.js or verified custom SVG | domain-specific chromosome geometry |
| read-level/SV evidence | JBrowse 2 preferred; IGV.js comparator | mature genomic evidence browsing |
| bespoke links/orientation diagrams | D3/SVG | exact custom geometry |
| tables | native semantic HTML/data grid | exact values, accessibility, sorting/filtering |

No chart library may become the source of scientific truth; renderers consume normalized data.

## 13. Migration test strategy

For each visualization added to the portable report:

1. fixture JSON has known normalized values;
2. renderer view-model is tested independently of HTML;
3. DOM/snapshot test checks labels, units, statuses and caveats;
4. no JavaScript test ensures fallback values remain present;
5. accessibility test checks table/text alternative and keyboard path;
6. static screenshot regression checks gross layout only, not scientific equality;
7. scientific equality is checked against the structured fixture, never against pixels;
8. GRCh37 and GRCh38 fixtures are separate;
9. empty/NO_CALL/FAILED/NOT_RUN fixtures are mandatory;
10. values such as VAF=0, CN=0, coverage=0 and empty event lists are tested explicitly so valid zeroes never disappear as falsey UI values.

## 14. First implementation backlog, ordered by value/risk

1. Introduce a report view-model layer so presentation stops reaching directly into `PipelineResult` everywhere.
2. Add a reviewer state strip with exact module semantics.
3. Add CNV fit and chromosome-consensus sections from the existing QDNAseq/ACE contract.
4. Add target-coverage matrix/distribution from `TargetCoverageReport` when present.
5. Expand SV evidence columns to expose available support/VAF/coverage/precision without inference.
6. Add evidence inspector/cross-selection in the portable HTML using local-only JS.
7. Add build-aware chromosome/segment visualization.
8. Perform JBrowse 2 versus IGV.js technical spike for local BAM/VCF review.
9. Add loopback-only evidence service and lazy read-level viewer.
10. Only after the above: reviewer annotation/audit UX.

This ordering deliberately extracts more value from data ONTSeq already produces before adding new biological inference.