# ONTSeq Multi-Caller Real-Tool Qualification

Date: 2026-09-20
Status: staged implementation
Branch: `feat/multicaller-realtool-qualification`
Base: `48614a1c41c5b7480f15dfe6246be53c72419499`
Scope: Research Use Only; Human Review Required

## Goal

Qualify the remaining adapter-qualified multi-caller paths against real external binaries on generated synthetic inputs without changing analytical interpretation, thresholds, caller ranking or clinical release behavior.

This phase is an engineering/interoperability program only.

## Baseline

PR #80 remains the stable integration baseline and is not modified by this branch.

At baseline:
- Severus paired: adapter-qualified, not real-tool-qualified.
- Wakhan phased CNA: adapter-qualified, not real-tool-qualified.
- all other currently qualified real-tool lanes remain unchanged.

## Severus 1.7

Primary target: `severus:paired`.

The current ONTSeq policy lock is Severus 1.7 and matches the current Bioconda package version.

Qualification requirements:
- pinned reproducible environment;
- generated synthetic tumor and matched-normal BAMs only;
- no patient/genomic payload committed;
- real binary invoked through the existing ONTSeq adapter;
- missing-normal fail-closed behavior retained;
- native somatic VCF and native auxiliary outputs retained;
- real-tool smoke may legitimately return NO_CALL on synthetic data, but a successful caller execution must still produce a coherent ONTSeq report;
- no use of output-read-IDs or alignment export in the qualification lane unless explicitly required by the adapter contract.

Real-tool success is not analytical sensitivity/specificity or clinical validation.

## Wakhan

Wakhan requires a separate version-contract step before real-tool qualification.

Current observations:
- ONTSeq adapter currently pins `0.5.0`;
- current Bioconda package advertises `0.4.4`;
- current upstream issue reports use a different public version line.

Therefore:
- do not silently repin Wakhan;
- establish the actual executable/package version contract first;
- any adapter migration must be RED -> GREEN and preserve command semantics, provenance and optional Severus-parent locking;
- only after a reproducible pinned runtime exists may `real_tool_qualified` become true.

## TDD order

1. Severus real-tool environment and opt-in test contract.
2. Verify the current Severus adapter against the real binary.
3. Fix only concrete runtime incompatibilities found by that smoke.
4. Record Severus real-tool qualification only after same-head success.
5. Add Wakhan version-contract RED tests.
6. Resolve Wakhan runtime/version mismatch prospectively.
7. Add Wakhan real-tool smoke only if a reproducible runtime can be pinned.
8. Final same-head verification and status update.

## Non-claims

This phase does not:
- analytically validate any caller;
- select or rank callers;
- introduce majority voting;
- alter CNV/SV interpretation thresholds based on observed outcomes;
- authorize clinical release;
- merge PR #80;
- bump package version merely because interoperability tests pass.
