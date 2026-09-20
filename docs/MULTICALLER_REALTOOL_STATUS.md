# ONTSeq Multi-Caller Real-Tool Qualification Status

Date: 2026-09-20  
Branch: `feat/multicaller-realtool-qualification`  
PR: #83 (Draft)  
Base integration head: `48614a1c41c5b7480f15dfe6246be53c72419499`  
Package version: `0.8.2`

Research Use Only. Human review remains required.

## Qualification meaning

**Real-tool-qualified** means the pinned external binary has executed successfully through the ONTSeq adapter on generated synthetic inputs in CI.

It does **not** mean:
- analytical sensitivity/specificity is established;
- limit of detection is established;
- biological truth is established;
- clinical validity is established;
- automatic clinical or ISCN release is authorized.

## Severus paired

Status: **real-tool-qualified**

Pinned runtime:
- Severus `1.7`;
- Bioconda package;
- paired tumor-normal mode only.

Verified path:
- generated synthetic tumor BAM;
- generated synthetic matched-normal BAM;
- indexed BAM inputs;
- standard alignment `NM` tags included;
- execution through `run_severus_paired`;
- version probe through the real binary;
- native somatic VCF and auxiliary caller outputs retained and fingerprinted;
- output-read-ID and alignment-export modes remain disabled;
- report remains Research Use Only;
- `COMPLETED` or `NO_CALL` are technical states only.

The first real-tool smoke correctly exposed that Severus 1.7 expects the standard BAM `NM` alignment tag. The synthetic fixture was corrected rather than weakening the adapter or external caller requirements.

Analytical validation: **not validated**.

## Wakhan phased CNA

Status: **real-tool-qualified for the pinned breakpoint-assisted tumor-only path**

Pinned runtime:
- Wakhan `0.4.4`;
- isolated Python 3.10 runtime for the external tool;
- ONTSeq test harness remains on supported Python 3.11;
- tumor-only phased-CNA mode;
- explicit Severus-style breakpoint input with exact parent lane ID and lane SHA-256;
- explicit checksum-locked centromere resource.

The real-tool qualification is deliberately narrower than the entire catalogued `wakhan:phased_cna` mode.

Qualified path:
- generated synthetic tumor BAM;
- generated phased tumor VCF;
- locked reference FASTA;
- locked centromere BED;
- explicit breakpoint VCF;
- exact breakpoint producer lane ID and SHA-256;
- Wakhan 0.4.4 binary invoked through `run_wakhan_phased_cna`;
- native solution ranking retained;
- caller-native integer CNA VCF retained;
- caller-native HP1 and HP2 integer copy-number BEDs retained separately;
- caller-native subclonal HP1/HP2 BEDs and VCF retained when emitted;
- no merge, overwrite or synthetic replacement of caller-native output.

Not separately real-tool-qualified:
- the change-point-detection-only path without breakpoint input;
- tumor-normal Wakhan mode;
- any alternative Wakhan version.

Interoperability findings resolved during qualification:
- the published/installable version contract was corrected from the stale `0.5.0` assumption to Wakhan `0.4.4`;
- Wakhan 0.4.4 requires its own Python 3.10 runtime, while ONTSeq itself remains on supported Python 3.11+;
- a single-contig synthetic reference triggers an upstream 0.4.4 header-parser edge case, so the generated fixture uses multiple `chr` contigs;
- the synthetic fixture was enlarged to provide sufficient coverage/phasing structure for CNA segmentation;
- the breakpoint-assisted lane uses an explicit checksum-locked centromere resource;
- Bioconda Wakhan 0.4.4 emits the legacy dual-haplotype layout rather than the later merged `integer_profile.bed` layout.

The ONTSeq adapter now accepts both:
- modern merged `integer_profile.bed`; or
- the native Wakhan 0.4.4 pair of HP1 + HP2 integer BEDs.

The adapter does not create a merged file or discard either haplotype-specific artifact.

Analytical validation: **not validated**.

## Verified implementation head

Implementation head: `ce652eab9c25ff4856796a6a7e57c7a15248083b`

Same-head verification:
- CI #994 — **SUCCESS**
  - Python 3.11 / 3.12 / 3.13;
  - repository safety — SUCCESS;
  - schema freshness — SUCCESS;
  - version agreement — SUCCESS;
  - wheel resource check — SUCCESS;
  - Ruff — SUCCESS;
  - mypy — SUCCESS;
  - **1895 passed, 11 skipped, 465 subtests**;
  - Snakemake DAG validation — SUCCESS;
  - synthetic report bundle — SUCCESS;
  - local-real-tool-smoke — SUCCESS;
- Tumor Severus real-tool #34 — **SUCCESS**;
- Tumor Wakhan real-tool #24 — **SUCCESS**;
- Desktop bundle contract #480 — **SUCCESS**;
- Desktop CI #739 — **SUCCESS**
  - relocatable Linux runtime — SUCCESS;
  - stock-Ubuntu full-system smoke without host Python/R — SUCCESS;
  - WPF build — SUCCESS;
  - BAM index/reference identity handling — SUCCESS;
  - self-contained win-x64 publish — SUCCESS;
  - complete first-run bundle verification — SUCCESS.

## Explicit non-claims

This qualification phase does not:
- rank callers;
- select a winning caller;
- use caller majority as truth;
- alter previously registered analytical thresholds;
- establish analytical or clinical validity;
- authorize clinical release;
- add patient/genomic payloads to Git;
- merge PR #80 or PR #83;
- bump the package version.

PR #83 remains Draft and unmerged. Package version remains `0.8.2`.
