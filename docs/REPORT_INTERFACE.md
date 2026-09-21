# Integrated Befund interface

The owner-supplied React design is integrated as `befund-interactive-v1`. It is the
default portable HTML presentation for newly rendered reports. Live Workspace also
offers **Interaktiver Befund herunterladen** for an existing run: this reformats its
recorded results in a temporary directory, without modifying the archived HTML,
normalized JSON, analysis outputs, checksums or review decisions. The original HTML,
Excel and JSON remain separately downloadable.

## Data boundary

The interface is generated from the typed `PipelineResult`, optional same-sample,
same-build `QDNAseqCallReport`, existing coverage and methylation reports, and
read-length histogram. The service checks available stage-output fingerprints before
reusing auxiliary files. Missing, conflicting, cross-sample, cross-build or corrupted
evidence is never silently substituted. PNG files must remain inside the CNV evidence
directory. Local paths are redacted from display data.

The uploaded design's embedded patient example, fixed interpretations, numerical
thresholds and prompts are not part of the source, wheel or report template. No
example records or screenshots may be committed. Tests use synthetic fixtures only.

## Presentation and interaction

The report covers key findings, module execution, original multi-bin CNV graphics,
purity/ploidy alternatives, exploratory VAF anchoring and dilution, chromosome tiles,
cytobands, all normalized events and caller evidence, QC, coverage, methylation,
ISCN, provenance, scope, limitations and review status. MARLIN, experimental dilution
series and clinical risk classification explicitly state when the result contract
contains no corresponding result; the interface does not run those analyses.

Selecting a stored solution updates the model summary, chromosome tiles and dilution
together. The original plots, normalized event calls and ISCN remain tied to the
recorded pipeline result. Changing bin size resets the selected model and dilution.

The explorer reconstructs a chromosome's relative signal from its consensus CN and
the **primary** stored fit:

`r = (f × CN + (1 − f) × 2) / (f × P + (1 − f) × 2)`

It inverts this expression under a user-selected stored alternative. These are
reconstructed consensus ratios, not raw segment measurements or a new fit. Invalid
fractions, non-finite numbers and impossible negative copy numbers yield no value.
X/Y are excluded. Zero is retained when actually present.

The VAF explorer requires explicit selection of a hypothetical neutral anchor and
labels VAF as hypothetical. It assumes a clonal heterozygous somatic variant in a
diploid neutral region (`f = 2 × VAF`). No fitted-error-to-detection-limit formula,
clinical threshold, automatic risk classification or release action is introduced.

## Offline packaging and fallback

`web/src/befund/` is compiled by `scripts/build_befund.mjs` with the same locked
React 19.2.0/esbuild 0.25.12 dependencies as Live Workspace. The generated
`src/ontseq_platform/report_bundle.html` is distributed in ordinary wheels. It
contains code, styles and third-party notices, never run data. No CDN, web font or
network request is needed. System font fallbacks preserve offline operation.

The complete escaped server-rendered report is retained. It remains visible without
JavaScript, can be opened from the interactive report, and reappears if React fails.
Print styling includes expanded evidence. The header records the presentation
revision separately from the analytical pipeline version. The report resume plan
also records the presentation revision, so an older report is not reused as current.

## Build and verification

```sh
cd web
npm ci
npm run build
npm run build:befund
npm run check:build
npm run check:befund
npm test
cd ..
make safety
make versions
make lint
make test
```

Browser verification must exercise a full CNV case, absent modules and values,
alternative selection, bin switch, VAF anchor, dilution, event search, evidence
expansion, offline loading, JavaScript-disabled fallback, desktop/mobile layouts and
print output. None of these technical tests establishes analytical or clinical validity.
