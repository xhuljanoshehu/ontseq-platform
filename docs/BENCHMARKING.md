# CNV and SV benchmarking

## Purpose

The first benchmark layer tests representation, matching and metric calculation with fully
synthetic normalized events. It is deterministic and safe for CI. It does not estimate clinical
performance and does not select a production caller.

## Run the synthetic cases

```bash
ontseq benchmark examples/benchmarks/synthetic_cnv.yaml \
  --output results/benchmarks/SYNTHETIC_CNV_001.benchmark.json

ontseq benchmark examples/benchmarks/synthetic_sv.yaml \
  --output results/benchmarks/SYNTHETIC_SV_001.benchmark.json
```

Or run both through Snakemake:

```bash
snakemake --snakefile workflow/benchmark.smk --cores 1 --use-conda
```

## Matching contract

| Event class | Match conditions |
| --- | --- |
| CNV | Same event type and chromosome, minimum reciprocal overlap, optional copy-number tolerance |
| Single-locus SV | Same event type/chromosome and both breakpoints within the configured distance |
| Paired SV/fusion | Same event type and chromosome pair; direct or swapped breakend order accepted |

Matching uses a deterministic maximum-cardinality one-to-one assignment, so an attractive local
pair cannot reduce the achievable true-positive count. Reports contain TP, FP, FN, precision,
recall, F1, matched pairs and unmatched event IDs. Undefined metrics remain `null`; they are never
replaced with a misleading zero or one.

## What this layer does not do

The normalized-event comparator is not VCF-representation aware. Production SV benchmarking
will retain source VCFs and use a version-pinned tool such as
[Truvari](https://github.com/ACEnglish/truvari), whose `bench` command evaluates a truth VCF
against a comparison VCF. Thresholds must be locked before viewing comparative results.

## Validation ladder

1. **CI fixtures:** deterministic synthetic positive, discordant and negative cases.
2. **Public technical truth:** version-pinned GIAB HG002 and draft HG008 tumor/normal resources.
3. **Analytical stress tests:** coverage and tumor/blast-fraction dilution, event size/class,
   reference build and target observability.
4. **Intended-use validation:** orthogonally characterized AML specimens on the locked wet-lab
   and computational workflow.

Public benchmarks support engineering decisions but cannot replace AML intended-use validation.
Every comparison must record dataset version, caller/container, reference, target BED, parameters,
matching thresholds, exclusions, runtime, memory and no-call behavior.
