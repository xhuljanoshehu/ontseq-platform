# Aligned-BAM MVP

## Intended use

This milestone is a research-only intake, descriptive-QC and candidate-SV path for exactly one
aligned ONT BAM. It establishes whether the file, index, header and locked reference are suitable
for downstream analysis, then optionally normalizes conservative Sniffles2 evidence. It does not
call CNVs, infer fusions, generate ISCN or establish a biological negative result.

## Inputs

- a validated single-sample manifest with `input.kind: aligned_bam`;
- coordinate-sorted BAM and matching BAI;
- a versioned reference lock generated from the exact FASTA index used for alignment;
- local `samtools`, `cramino` and Sniffles2 executables;
- the versioned QC policy in `configs/qc/defaults.yaml`.
- the technical Sniffles2 policy in `configs/sv/sniffles2.conservative.technical.yaml`.

Raw genomic data remain outside Git and on approved institutional compute.

## Reference lock

Create the lock beside the institutionally controlled reference bundle:

```bash
ontseq reference-lock \
  --fai /approved/references/GRCh38.fa.fai \
  --reference-id GRCh38-analysis-set-v1 \
  --genome-build GRCh38 \
  --output /approved/references/GRCh38-analysis-set-v1.lock.json
```

The lock records every contig and length plus the SHA-256 of the FAI. It does not copy the FASTA.
The manifest `reference_id` and genome build must match the lock.

## Fail-closed intake gate

```bash
ontseq inspect-bam sample.manifest.yaml \
  --reference-lock /approved/references/GRCh38-analysis-set-v1.lock.json \
  --output results/SAMPLE.intake.json
```

| Check | Behavior |
| --- | --- |
| BAM/BAI present | Hard failure when either file is missing |
| Manifest checksum | Hard failure when supplied SHA-256 differs |
| `samtools quickcheck` | Checks header and BAM EOF; hard failure on error |
| Header sort order | Requires `SO:coordinate` |
| Sequence dictionary | Compares names and lengths with the reference lock |
| Read groups | Missing/incomplete read groups produce a warning |
| `samtools idxstats` | Confirms the index is readable and compatible |

The report deliberately omits source paths and read-group sample values. Optional `--checksum`
adds BAM/BAI fingerprints for an immutable run envelope. According to the official
[`samtools quickcheck` documentation](https://www.htslib.org/doc/samtools-quickcheck.html), the
command does not read the middle of the file and therefore cannot exclude internal corruption.

## Descriptive Cramino QC

```bash
ontseq qc-cramino sample.manifest.yaml \
  --policy configs/qc/defaults.yaml \
  --output results/SAMPLE.qc.json
```

The adapter uses Cramino JSON output and normalizes read count, aligned percentage, yield, mean
coverage, N50/N75, length and identity metrics. Cramino's source filename, path and creation time
are not copied. The official [Cramino interface](https://github.com/wdecoster/cramino) is the
command contract.

Numeric gates remain `null` until analytical validation. Therefore a successful descriptive run
is `WARN`, not a clinically meaningful `PASS`. Configured but missing metrics fail visibly.

## Conservative Sniffles2 candidate adapter

```bash
ontseq call-sniffles sample.manifest.yaml \
  --intake results/SAMPLE.intake.json \
  --policy configs/sv/sniffles2.conservative.technical.yaml \
  --vcf results/SAMPLE.sniffles.vcf \
  --output results/SAMPLE.sniffles.json
```

The locked command contract uses Sniffles2 v2.8.0 with explicit `--minsupport`, `--minsvlen`,
`--mapq`, `--pass-only`, `--symbolic` and `--no-progress` settings. The adapter never enables
`--output-rnames`. Germline versus mosaic mode is explicit policy; the technical default exists
for engineering tests and is not an AML reportability threshold.

Accepted DEL, DUP, INV, INS and BND/TRA records are normalized into 0-based, half-open event
coordinates. Support, QUAL, VAF, strand orientation, coverage context and mean alignment NM are
retained when emitted. Raw IDs, read names, ALT sequences and source paths are not copied.
Malformed, filtered, low-support, disallowed or non-canonical records receive counted rejection
reasons. Zero accepted records becomes `NO_CALL`, never a biological negative.

The command interface follows the official
[Sniffles2 repository](https://github.com/fritzsedlazeck/Sniffles); VCF fields are normalized from
its [current VCF implementation](https://github.com/fritzsedlazeck/Sniffles/blob/master/src/sniffles/vcf.py).

## Assemble and render

```bash
ontseq assemble-aligned-mvp sample.manifest.yaml \
  --intake results/SAMPLE.intake.json \
  --qc results/SAMPLE.qc.json \
  --sniffles results/SAMPLE.sniffles.json \
  --git-commit COMMIT_SHA \
  --output results/SAMPLE.result.json

ontseq render results/SAMPLE.result.json --output-dir results/
```

The HTML and Excel reports include a module-status table. Candidate SVs are visible but remain
`unclassified` and `reportable: false`; CNV, fusion and ISCN remain `NOT_RUN` with a reason. An
empty event table must never be interpreted as a negative finding.

## Reproducible real-tool smoke test

```bash
PYTHONPATH=src python -m ontseq_platform local-smoke \
  --output-dir results/local-smoke
```

The command generates an identifier-free SAM fixture at runtime, converts/sorts/indexes it with
real samtools, runs Cramino, calls a known synthetic 200 bp deletion with Sniffles2 and renders the
same JSON/HTML/XLSX path used by the application. Intermediate SAM and unsorted BAM are removed;
all remaining genomic artifacts stay under ignored `results/`. CI uploads only selected synthetic
reviewer artifacts, never BAM or VCF.

## Snakemake entry point

After replacing the synthetic paths in `workflow/config/aligned_bam.example.yaml`:

```bash
snakemake --snakefile workflow/aligned_bam.smk --cores 4 --use-conda
```

The DAG stops scientific analysis when intake fails. Sniffles2 candidates can flow into the
report, but no event becomes reportable until its assay-specific benchmark gate passes. Without
`--use-conda`, the same rules use executables from the active environment.
