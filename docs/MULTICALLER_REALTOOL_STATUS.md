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
- Severus `1.7`
- Bioconda package
- paired tumor-normal mode only

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

Status: **adapter-qualified; real-tool qualification pending**

The existing ONTSeq adapter remains fail-closed and retains:
- exact phased-input fingerprints;
- optional breakpoint parent lane ID and lane SHA-256;
- runtime-script fingerprint;
- caller-native outputs.

Version-contract issue:
- current ONTSeq adapter pins `0.5.0`;
- current Bioconda package advertises `0.4.4`;
- recent upstream usage reports a different public version line.

No silent repin is permitted. Wakhan real-tool qualification starts only after the executable/version contract is prospectively resolved.

Analytical validation: **not validated**.

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
