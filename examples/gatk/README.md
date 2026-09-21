# Local usage: experimental GATK adapter

These are Linux/WSL commands from the ONTSeq repository root with the reviewed adapter package present.
Use the project's existing approved Python environment with Pydantic 2.8–2.x.
The adapter itself requires no Python library beyond Pydantic and the standard library.
GATK and Java are separate local prerequisites, not installed by these commands.

All sample names below are illustrative pseudonyms. Replace the paths, names, reference
and assay labels with reviewed local values; do not publish a populated manifest to Git.
The tumor/normal sample names must equal their BAM read-group `SM` names.

## 1. Create a disabled, content-fingerprinted configuration

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
umask 077
python -m ontseq_platform.gatk_adapter lock \
  --run-id RESEARCH_RUN_001 \
  --assay-id REVIEWED_ONT_ASSAY_V1 \
  --reference-id REVIEWED_GRCH38_REFERENCE \
  --genome-build GRCh38 \
  --reference /secure/reference/GRCh38.fa \
  --tumor-bam /secure/research/TUMOR_001.bam \
  --tumor-sample TUMOR_001 \
  --germline /secure/resources/population-af.vcf.gz \
  --gatk-jar /opt/gatk/gatk-package-4.6.2.0-local.jar \
  --output /secure/research/gatk.config.json
```

For paired mode add both:

```text
--normal-bam /secure/research/NORMAL_001.bam --normal-sample NORMAL_001
```

Optional inputs, after their provenance and assay compatibility are reviewed:

```text
--pon /secure/resources/reviewed-ont-pon.vcf.gz
--contamination-sites /secure/resources/population-common-snps.vcf.gz
--intervals /secure/reference/research-targets.bed
```

`lock` computes actual input SHA-256 values; no invented hashes need to be edited into a
template. It records the supplied reference/assay association as an attestation, not a
provenance certification. It refuses to overwrite an existing configuration. It does not
run GATK and leaves both execution opt-ins false.

## 2. Inspect the command plan without executing it

```bash
python -m ontseq_platform.gatk_adapter plan \
  --config /secure/research/gatk.config.json \
  --output-dir /secure/research/gatk-run-001
```

Planning prints JSON to the terminal, reports `NOT_RUN`, does not create the output
directory and does not verify inputs or software. Plan output includes local file paths;
handle it as sensitive. These examples do not prove that any supplied BAM is valid.

## 3. Explicitly attempt a research run

```bash
python -m ontseq_platform.gatk_adapter run \
  --config /secure/research/gatk.config.json \
  --output-dir /secure/research/gatk-run-001 \
  --enable \
  --allow-experimental-ont
```

The output directory must not already exist. A rerun requires a new directory and should
use a new run ID. Omitting enablement returns `NOT_RUN`; CLI exit codes are 0=completed,
1=failed, 2=not run. Both opt-ins may alternatively be set explicitly in a reviewed
configuration. The executed effective configuration and its hash are retained.

The first native use should be a separately reviewed synthetic, valid BAM/FASTA
interoperability fixture, not a patient analysis. The supplied unit tests' text placeholders
are deliberately not suitable for that native qualification.

## 4. Inspect results without promoting them

Read `result.json`, `commands.json`, `config.lock.json`, step logs, native VCFs and statistics.
A completed result does not imply a validated assay; PASS does not authorize clinical
release. Zero candidates, failure and non-execution must not be summarized as negative.

Run only the delivered adapter's tests with:

```bash
PYTHONPATH=src python -m pytest -q tests/test_gatk_adapter.py
```

Before integrating into the complete project, separately run its required gates:

```bash
make safety
make versions
make lint
make test
```

The exact reviewed revision must pass all of those upstream gates before merge. The new schema
files are isolated documentation artifacts, not registrations in ONTSeq's canonical
schema exporter or multi-caller planner.
