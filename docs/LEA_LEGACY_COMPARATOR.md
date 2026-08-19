# Lea historical comparator

## Status

Research use only. This is a **compatibility and validation lane**, not a production caller.
It exists so results from the 2026 Lea Evers ONTseq workflow can be reproduced, normalized
and benchmarked against newer ONTSeq methods without copying the historical source code or
promoting its parameters to clinical defaults.

The historical source repository remains external. Raw BAM/VCF files, patient data and the
local historical truth tables must stay on approved storage and must not be committed to
Git.

## Why this lane exists

The historical workflow is scientifically valuable because it already contains an AML-oriented
QDNAseq + ACE path, cytoband-level CNV output, SV/fusion work and a readable HTML report.
The new platform has a stronger contract, observability and validation architecture but should
not throw away that prior work. The correct relationship is therefore:

```text
historical workflow output
        |
        v
strict compatibility adapter
        |
        v
canonical research-only CnvCallSet
        |
        +--> shared observability mask
        +--> shared truth set
        +--> same CNV benchmark engine as every newer caller
```

The adapter reproduces **output semantics**, not implementation. It does not import, vendor or
execute Lea's Python/R source.

## Frozen ACE compatibility profile

The first implemented profile corresponds to the supplied 2026 workflow snapshot:

| Property | Historical compatibility value | Meaning in ONTSeq |
| --- | --- | --- |
| reference | hg19 / GRCh37 | fixed; no lift-over |
| QDNAseq bin size | 1000 kbp | historical default profile, not a clinical default |
| ACE penalty | 0.6 | historical workflow parameter |
| assumed autosomal ploidy | 2 | compatibility metadata |
| affected cytoband threshold | `frac_abr >= 0.66` | historical reporting rule |
| whole-chromosome output | `CN.csv` | `Chromosome,Copies,Ploidy,CNA` |
| partial-band output | `dels_dups.csv` | `chromosome,name,event,frac_abr` |
| QDNAseq version | 1.38.0 in supplied `renv.lock` | declaration only |
| QDNAseq.hg19 version | 1.32.0 in supplied `renv.lock` | declaration only |
| ACE version | 1.20.0 in supplied `renv.lock` | declaration only |
| runtime R package versions | not proven by the container definition | explicit provenance warning |

The last line matters: the supplied R container installs Bioconductor packages without exact
package-version pins. A lock-file declaration therefore must not be presented as the proven
runtime version of a historical result.

## Fail-closed import rules

`src/ontseq_platform/cnv/lea_compat.py` implements the compatibility boundary. It deliberately
rejects ambiguous or internally inconsistent artifacts instead of repairing them.

The importer:

1. maps columns by header name, never position;
2. verifies `CNA = Copies - Ploidy` for every `CN.csv` row;
3. rejects duplicate/unsupported chromosomes and non-finite values;
4. requires the frozen GRCh37 profile and a matching build-locked cytoband table;
5. never performs coordinate lift-over;
6. converts whole-chromosome CNA rows to full reference-contig intervals;
7. converts partial calls only through exact cytobands from the locked resource;
8. does **not** invent an absolute copy number for `dels_dups.csv`, because that artifact
   does not contain one;
9. rejects a partial-band event on a chromosome already declared as a whole-chromosome CNA,
   because the historical ACE script suppresses that combination;
10. retains `reportable=false` and `research_only=true` through the shared `CnvCallSet`.

A structurally valid pair of outputs with no altered segment may be represented as a completed
historical-method negative. That statement is bounded by the method, its thresholds and the
supplied observability mask. It is never rendered as an assay-wide or clinical negative.

## Historical truth tables

`src/ontseq_platform/cnv/lea_truth_tables.py` parses the two local evaluation tables that were
supplied separately:

- `gt.tsv`: sample identifier plus cytogenetic ISCN string;
- `gt_full.csv`: cytogenetic karyotype, ONT karyotype, cellularity and the historical
  complex/monosomal/MRC/MRCA/MRA labels.

The parser preserves those strings and labels exactly. It does **not** silently correct uncertain
ISCN, recompute the classifications or turn old labels into production rules. Conversion to CNV
truth goes through the existing `cnv.truth` boundary, where unsupported constructs remain explicit.

The data files themselves are intentionally absent from the repository.

## Product / reviewer experience

The comparator belongs in a **Validation** workspace, not the routine operator path. A future UI
should make that distinction visible before any file is selected.

### Validation flow

```text
Validation > Historical comparator
    1. Select local historical result folder
    2. Select locked GRCh37 reference + cytoband resource
    3. Import CN.csv + dels_dups.csv
    4. Optional: attach local cytogenetic truth table
    5. Run compatibility gates
    6. Review normalized events and refused/unsupported content
    7. Run shared benchmark against another method
    8. Export comparison report
```

### Required visual states

The interface should never use a generic green check for every successful import. It should show
separate states:

- **READY FOR BENCHMARK** — artifacts are structurally compatible;
- **BUILD MISMATCH** — stop; no lift-over offered inside compatibility mode;
- **INCOMPLETE HISTORICAL ARTIFACTS** — stop; do not reinterpret as negative;
- **TRUTH CONVERSION INCOMPLETE** — benchmark blocked until unsupported cytogenetic constructs
  are adjudicated or an explicit partial-truth analysis is chosen;
- **RUNTIME VERSION UNVERIFIED** — import may continue for historical comparison, but the
  provenance limitation remains visible in every exported report.

### Comparison screen

The high-value screen is not a recreation of Lea's old report. It is a side-by-side method
comparison driven by the canonical contracts:

```text
Sample / stratum / truth source
-------------------------------------------------------------
Historical QDNAseq+ACE    | New candidate caller
version provenance        | version + container digest
coverage / observability  | coverage / observability
whole-chromosome events   | whole-chromosome events
partial events            | partial events
NO_CALL regions           | NO_CALL regions
-------------------------------------------------------------
shared evaluable genome
base-level concordance
truth-event detection
confirmation rate
copy-number error
unsupported truth constructs
```

Segmentation differences are shown for debugging but are not scored as biological false
positives. That remains the responsibility of the shared CNV benchmark engine.

## What is intentionally not implemented yet

### Spectre compatibility

The supplied historical code proves that the active workflow generated a Spectre CNV BED plus a
`spectre_CN.csv`, and a separate evaluation helper summarized the `avg` field by chromosome. The
exact representative raw output artifact has not yet been supplied. A production mapping is
therefore **not guessed**. When one or more real, non-sensitive example outputs are provided, the
format can be frozen and added behind the same adapter boundary.

### Cytoband-collapse / ISCN parity

Lea's historical cytoband-merging behavior is useful as a comparator, but it should not be wired
straight into final ISCN. The next step is an independently implemented, test-driven band-region
normalizer with explicit p/q-arm and centromere behavior. Its output will feed an expert-reviewed
ISCN proposal rather than an automatic clinical release.

### Clinical report parity

The old HTML report remains a UX reference. The new platform should reproduce useful information
architecture, not the old data flow. HTML/XLSX/desktop views must render from the canonical result
contract so that `NO_CALL`, `FAILED` and `NOT_RUN` cannot collapse into "nothing found".

## Acceptance criteria for the next increment

1. Import at least one locally authorized historical `CN.csv` + `dels_dups.csv` pair and compare
   the normalized events with the original report by expert review.
2. Add a representative raw Spectre output and freeze its actual schema before implementing the
   Spectre adapter.
3. Add the real checksummed hg19 cytoband resource to the on-prem reference bundle, not Git.
4. Run the historical 100-sample comparison table through the truth parser locally and quantify
   how many karyotypes are fully convertible versus explicitly unsupported.
5. Pre-register the evaluation mask, truth resolution and method-comparison thresholds before
   looking at caller rankings.
6. Only then compare QDNAseq+ACE, Spectre and newer candidate methods on the same sample/coverage/
   cellularity strata.
