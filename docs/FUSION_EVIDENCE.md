# DNA fusion evidence architecture

**Status:** research-only design and implementation scaffold  
**Clinical reportability:** disabled  
**Primary lane:** structural-variant evidence from aligned Oxford Nanopore DNA

## Purpose

The platform must not equate a VCF breakend (`BND`) with a confirmed fusion gene. A DNA
breakpoint can support a fusion hypothesis, but gene overlap, transcriptional direction,
reading frame, expression, oncogenic relevance and assay observability are separate questions.

The intended evidence path is:

```text
Sniffles2 / independent SV caller
        |
        v
paired DNA breakpoints
        |
        +--> VCF breakend adjacency/orientation syntax
        +--> caller support / VAF / precision / local coverage
        +--> build-locked gene and transcript annotation
        +--> target / coverage observability for both partners
        |
        v
FusionCandidate (research_only=true, reportable=false)
        |
        +--> known-pair annotation (evidence feature, not truth)
        +--> transcript/orientation plausibility (later gate)
        +--> independent caller / RNA / FISH / RT-PCR evidence
        |
        v
expert review and assay-specific analytical validation
```

## Open-source components and why they are not copied wholesale

### NanoFG

NanoFG is directly relevant precedent for Nanopore **DNA** fusion detection. It combines
long-read structural-variant detection with fusion-gene interpretation and accepts BAM/VCF
inputs. The associated FUDGE/NanoFG publication demonstrated targeted long-read fusion
resolution and nucleotide-level breakpoint characterization.

Repository: <https://github.com/SdeBlank/NanoFG>  
Publication: <https://doi.org/10.1038/s41467-020-16641-7>

NanoFG is used here as an architectural and benchmark reference rather than vendored source
code. Reasons:

- the repository is GPL-3.0, while this repository currently has no public license and its
  institutional/IP disposition is unresolved;
- its pinned toolchain is substantially older than the current ONTSeq platform stack;
- our platform requires explicit `NO_CALL`, target observability, build/resource checksums,
  privacy-preserving reviewer artifacts and caller provenance;
- a clinical or near-clinical workflow needs current independent benchmarking rather than
  inheritance of historical defaults.

A later benchmark may execute NanoFG as an external comparator on synthetic/public or locally
approved validation material without incorporating its source code.

### Sniffles2

Sniffles2 v2.8.0 remains the first normalized DNA SV evidence source. Its `BND` records provide
paired genomic loci, but the current normalized event model previously discarded the VCF ALT
breakend form. `breakends.py` now preserves the four VCF adjacency forms as privacy-safe tokens
without retaining ALT sequence, read names or inserted sequence.

This syntax is retained as evidence only. It is not yet converted into 5-prime/3-prime gene
orientation until that mapping is covered by explicit tests against the VCF specification and
orthogonally characterized rearrangements.

VCF specification: <https://github.com/samtools/hts-specs>

### Ensembl VEP / SnpEff

Mature annotation engines should be preferred over reimplementing transcript consequence
annotation. Ensembl VEP and SnpEff are candidate external annotation backends for structural
variants. The platform-specific contract remains independent of either backend so resources can
be version-locked and benchmarked.

The initial implementation deliberately consumes a small BED6/BED7 gene interval contract to
exercise the fusion evidence model without silently downloading or embedding a reference bundle.
A production adapter must lock genome build, annotation release, exact source checksum and
transcript policy.

VEP structural-variant documentation: <https://www.ensembl.org/info/docs/tools/vep/vep_formats.html>

### JAFFAL / LongGF

JAFFAL and LongGF are long-read **RNA** fusion callers and are therefore not replacements for the
current DNA breakpoint lane. They are relevant later as an orthogonal RNA confirmation/research
lane. RNA evidence can answer questions that genomic DNA alone cannot, including whether a
candidate junction is expressed and which transcript configuration is present.

## Current contracts

`GeneAnnotationIndex`
: consumes a local build-declared BED6/BED7 file and records its SHA256 fingerprint.

`FusionBreakpoint`
: stores genomic position, overlapping/nearby genes, transcript IDs when supplied, observability
  and a deliberately unresolved breakend-orientation field.

`FusionGenePair`
: stores an unordered gene pair until orientation is genuinely resolved. The model rejects code
  that assigns 5-prime/3-prime labels while `orientation_resolved=false`.

`FusionCandidate`
: preserves source event ID and caller evidence and is permanently `research_only=true` and
  `reportable=false` in this contract version.

`BreakendDescriptor`
: preserves the four VCF BND ALT forms and mate coordinates while explicitly discarding inserted
  sequence and VCF record identifiers.

## Observability rule

A fusion-negative statement requires validated evidence that the relevant breakpoint space was
observable. In Adaptive Sampling, coverage of one selected gene does not prove coverage of an
unknown partner or of a breakpoint outside the selection design.

The fusion layer therefore accepts explicit per-region observability states:

- `observable`
- `limited`
- `not_assessable`
- `unknown`

PR #7 target-coverage output is the intended future connection point, but technical depth bins
must not become implicit clinical reportability thresholds.

## Known-pair handling

A known gene pair (for example a canonical hematologic fusion pair) may be annotated as an
evidence feature. It must not:

- turn a BND into a confirmed fusion automatically;
- imply the observed genomic orientation creates the canonical transcript;
- override inadequate partner coverage;
- infer expression or reading frame;
- set `reportable=true` before local validation and expert review.

The final controlled known-fusion resource should be versioned and source-attributed rather than
hard-coded in Python.

## Immediate next gates

1. Join `BreakendDescriptor` back into `FusionCandidate` and test all four VCF BND forms.
2. Add a version-pinned annotation backend (VEP or another benchmarked equivalent) behind the
   current contract.
3. Add synthetic canonical and non-canonical hematologic rearrangements with orientation edge
   cases; do not use patient-derived fixtures in Git.
4. Connect PR #7 coverage/target information to both breakpoints and explicitly represent partner
   non-observability.
5. Benchmark Sniffles2 against at least one independent DNA SV method under the actual assay lane.
6. Define an orthogonal truth cohort (FISH, RT-PCR, RNA sequencing and/or independently validated
   genomic breakpoints) before estimating sensitivity or negative predictive performance.
7. Only after validation, design expert-review evidence tiers and any clinical reportability
   policy.

## What this implementation does not prove

It does not prove that a candidate is expressed, in-frame, oncogenic, disease-defining or
clinically reportable. It does not validate sensitivity of the 111-gene Adaptive Sampling assay,
and it does not establish that absence of a BND is a negative fusion result.
