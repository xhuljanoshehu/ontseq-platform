# ISCN engine v0.2

## Status

Research Use Only. The generated ISCN string is an expert-reviewable proposal. Technical
validation must never be interpreted as authorization for automatic clinical release.

## Design principle

ONTSeq keeps the normalized genomic event as the source of truth and treats ISCN as one output
representation:

```text
BAM / caller outputs
        |
        v
normalized GenomicEvent
        |
        +--> coordinate-to-cytoband mapping (locked GRCh37/GRCh38 reference)
        |
        v
conservative ONTSeq ISCN renderer
        |
        +--> built-in subset validator
        +--> optional independent iscn-authenticator validator
        |
        v
expert-reviewable ISCN proposal
```

This avoids coupling biological evidence to a presentation string and makes it possible to improve
or replace the renderer without changing the underlying events.

## Implemented subset

The v0.2 renderer can emit only forms for which ONTSeq has an explicit structured event model and a
conservative mapping path:

- whole-chromosome gain (`+`)
- whole-chromosome loss (`-`)
- deletion (`del`)
- duplication (`dup`)
- inversion (`inv`)
- two-breakpoint translocation (`t`)

Insertions, fusions requiring a more specific cytogenomic representation, derivatives, marker
chromosomes, rings, dicentrics, mosaics/clones and other advanced constructs are not guessed. A
reportable event outside the implemented subset is omitted from the proposal and produces an
explicit warning.

## Cytoband reference layer

`src/ontseq_platform/iscn_reference.py` parses the five-column UCSC cytoband format into an immutable
`CytobandIndex`. Lookup uses zero-based, half-open genomic intervals. The index provides:

- coordinate -> cytoband lookup
- interval -> start/end cytoband lookup
- centromere bounds derived from `acen` records
- strict interval validation

The reference build is always explicit. GRCh37/hg19 and GRCh38/hg38 must never be mixed.

### Fetch locked source files

```bash
python scripts/fetch_iscn_reference_data.py \
  --output-dir references/iscn/cytobands
```

The script downloads the UCSC `cytoBand.txt.gz` tables for hg19 and hg38, decompresses them, checks
that a plausible number of human chromosome records was received, and writes a JSON manifest with
URL, SHA-256, record count and file size.

Official upstream locations used by the fetcher:

- `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz`
- `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz`

The downloaded files are reference inputs, not patient data. They should still be locked by checksum
for every validated software release.

## SQLite knowledge base

The knowledge base intentionally does **not** reproduce the ISCN 2024 publication. It stores only
open/reference coordinate data, provenance, implementation metadata, a limited set of reformulated
rules needed by the software, and synthetic/reference examples.

Build it with:

```bash
python scripts/build_iscn_kb.py \
  --hg19 references/iscn/cytobands/hg19_cytoBand.txt \
  --hg38 references/iscn/cytobands/hg38_cytoBand.txt \
  --output references/iscn/ontseq_iscn_kb_v0.2.sqlite
```

For schema development without downloading cytobands:

```bash
python scripts/build_iscn_kb.py \
  --bootstrap-only \
  --output /tmp/ontseq_iscn_kb_v0.2.bootstrap.sqlite
```

Main tables:

- `metadata`
- `sources`
- `genome_builds`
- `cytobands`
- `centromeres`
- `event_types`
- `rules`
- `examples`
- `external_validators`

## Validation strategy

### Layer 1: ONTSeq subset validator

`validate_subset()` verifies only the grammar that ONTSeq itself is allowed to emit. It is designed
to fail closed on unsupported fragments and formatting errors.

### Layer 2: independent validator

The optional dependency is pinned separately:

```bash
pip install -e '.[iscn]'
```

This installs `iscn-authenticator==0.2.1`, an MIT-licensed independent ISCN 2024 parser/rule engine.
If it is installed, `validate_iscn()` prefers it. If it is absent, the pipeline remains operable and
falls back to the narrower local validator.

The external package is a technical cross-check. It has not been analytically validated by this
project and must not be used as evidence that an ONTSeq result is clinically correct.

## Important safety behavior

The renderer deliberately refuses to invent information:

- no cytoband -> no structural ISCN fragment unless a coordinate reference is available
- simple deletion/duplication spanning p and q -> not auto-rendered
- unsupported event type -> omitted with warning
- validation failure -> explicit warning; never silently accepted
- ISCN proposal -> always retains expert-review requirement

## Validation work still required

Before promotion beyond research use, at minimum validate against locked, expert-reviewed cases for:

1. GRCh37 and GRCh38 coordinate-to-band mapping, including exact band boundaries.
2. Whole-chromosome gains/losses and chromosome-count consistency.
3. Large and focal deletions/duplications on both arms.
4. Inversions and translocations with independently confirmed breakpoints.
5. Complex karyotypes and event ordering.
6. Mosaic/subclonal cases and multiple cell lines (not implemented yet).
7. Derivative chromosomes and complex rearrangements (not implemented yet).
8. Uncertain breakpoints and no-call behavior.
9. Round-trip tests where supported: `GenomicEvent -> ISCN -> independent parser -> structured event`.
10. Cytogenetic expert sign-off on representative intended-use AML cases.

## Historical project context

The Lea Evers thesis describes the earlier pipeline's automatic ISCN generation and coordinate-aware
cytoband collapsing. It is useful historical context for traceability, but repository architecture,
tests, versioned reference locks and validation evidence remain the technical source of truth for the
current platform.
