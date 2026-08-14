# Aligned-BAM MVP

## Intended use

This milestone is a research-only technical intake and descriptive-QC path for exactly one
aligned ONT BAM. It establishes whether the file, index, header and locked reference are suitable
for downstream analysis. It does not call CNVs, SVs or fusions and cannot establish a biological
negative result.

## Inputs

- a validated single-sample manifest with `input.kind: aligned_bam`;
- coordinate-sorted BAM and matching BAI;
- a versioned reference lock generated from the exact FASTA index used for alignment;
- local `samtools` and `cramino` executables;
- the versioned QC policy in `configs/qc/defaults.yaml`.

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

## Assemble and render

```bash
ontseq assemble-aligned-mvp sample.manifest.yaml \
  --intake results/SAMPLE.intake.json \
  --qc results/SAMPLE.qc.json \
  --git-commit COMMIT_SHA \
  --output results/SAMPLE.result.json

ontseq render results/SAMPLE.result.json --output-dir results/
```

The HTML and Excel reports include a module-status table. CNV, SV, fusion and ISCN are marked
`NOT_RUN` with a reason. An empty event table must never be interpreted as a negative finding.

## Snakemake entry point

After replacing the synthetic paths in `workflow/config/aligned_bam.example.yaml`:

```bash
snakemake --snakefile workflow/aligned_bam.smk --cores 4 --use-conda
```

The DAG stops before Cramino when intake fails. Scientific caller rules will be connected only
after their benchmark gates pass. Without `--use-conda`, the same rules use executables from the
active environment.
