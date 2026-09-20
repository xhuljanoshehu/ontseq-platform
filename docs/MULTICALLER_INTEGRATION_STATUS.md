# ONTSeq Multi-Caller Integration Status

Date: 2026-09-20  
Branch: `feat/multicaller-integration`  
PR: #80 (Draft)  
Package version: `0.8.2`

Research Use Only. Human review remains required.

## Qualification ladder

The labels below are deliberately separate.

- **Catalogued** — the provider/mode has a typed ONTSeq catalog entry.
- **Planner-qualified** — the fail-closed planner understands its assay regime, required inputs
  and declared caller dependencies.
- **Adapter-qualified** — ONTSeq has a typed adapter/runtime or an explicit bridge to an existing
  ONTSeq runtime, with synthetic contract tests.
- **Real-tool-qualified** — the pinned technical path has executed against the real external tool
  on generated synthetic data in CI.
- **Analytically validated** — assay-specific analytical performance has been established against
  an appropriate prospectively registered truth set.

Real-tool qualification is an engineering/interoperability statement only. It does not imply
analytical sensitivity, specificity, limit of detection, biological truth, clinical validity or
clinical reportability.

## Current matrix

| Provider mode | Catalogued | Planner-qualified | Adapter-qualified | Real-tool-qualified | Analytically validated | Current boundary |
| --- | --- | --- | --- | --- | --- | --- |
| `qdnaseq_ace:multibin` | yes | yes | yes | yes | no | Existing QDNAseq+ACE runtime plus Block-4 full-evidence adapter; 100/500/1000-kbp evidence retained. |
| `spectre:depth_only` | yes | yes | yes | yes | no | Independent depth-only Spectre lane; no Sniffles parent is consumed. |
| `spectre:sniffles_supported` | yes | yes | no | no | no | Dependency is modelled prospectively; no executable ONTSeq adapter is qualified yet. |
| `sniffles2:standard` | yes | yes | yes | yes | no | Existing ONTSeq runtime; real Sniffles2 is exercised in `local-real-tool-smoke`. |
| `sniffles2:mosaic` | yes | yes | no | no | no | Kept separate from standard mode; no implicit promotion from the standard runtime. |
| `cutesv:standard` | yes | yes | yes | yes | no | Existing ONTSeq runtime; real cuteSV is exercised in `local-real-tool-smoke`. |
| `severus:paired` | yes | yes | yes | no | no | Paired tumor-normal adapter exists; missing normal fails closed; no real-binary qualification yet. |
| `severus:single_sample` | yes | yes | no | no | no | Catalogued separately; it must not inherit somatic paired claims. |
| `savana:paired` | yes | yes | yes | yes | no | Explicit paired tumor-normal SV+CNA mode; no tumor-only fallback. |
| `savana:tumor_only` | yes | yes | yes | yes | no | Explicit tumor-only SV+CNA mode; distinct plan/runtime identity from paired mode. |
| `wakhan:phased_cna` | yes | yes | yes | no | no | Phased CNA adapter exists; optional breakpoint parent is content-locked; policy remains real-tool unqualified. |
| `ichorcna:ulp_wgs` | yes | yes | yes | yes | no | Restricted to cfDNA/ULP-WGS; ONT/Adaptive-Sampling transfer is not implicitly qualified. |

No entry is analytically validated by this integration program.

## Shared planner and evidence guarantees

The multi-caller planner:

- keeps the existing canonical one-provider-per-stage runner unchanged;
- distinguishes `ELIGIBLE`, `INELIGIBLE`, `RUNTIME_UNAVAILABLE` and `NOT_REQUESTED`;
- never converts ineligibility or missing runtime into a biological negative;
- seals each lane with deterministic content identity;
- includes exact parent-lane hashes in dependent lane identity;
- rejects dependency cycles before execution;
- distinguishes required and optional caller parents;
- prevents paired modes from silently degrading to tumor-only modes.

Collection-level evidence addresses include the sealed lane identity, so identical native record names
from different specimens, runs or caller lanes cannot collide.

## Existing-caller bridges

The multi-caller layer reuses rather than rewrites the established ONTSeq implementations for:

- QDNAseq+ACE;
- Sniffles2 standard;
- cuteSV standard.

The bridge contract explicitly states that caller agreement is not biological truth and that native
output retention remains required.

## Spectre depth-only

The qualified Spectre lane is intentionally depth-only.

It:

- does not pass Sniffles evidence;
- fingerprints its inputs and native outputs;
- uses a pinned technical runtime;
- has a synthetic real-tool CI lane.

A future Sniffles-supported Spectre mode remains a separate dependency-bearing mode.

## Tumor auxiliary inputs

Tumor-aware caller inputs are bound to:

- analysis sample identity;
- source sample identity;
- genome build;
- reference ID;
- reference SHA-256;
- artifact SHA-256;
- exact producer lane ID and producer lane SHA-256 when an upstream caller produced the input.

This prevents silent reuse of SNP, BAF, phasing, breakpoint or matched-normal evidence from another
sample, build or upstream run.

## Severus

The current Severus implementation is paired tumor-normal only.

The adapter:

- requires the registered matched normal;
- checks BAM and index identity;
- pins runtime policy identity;
- retains native somatic VCF, breakpoint/cluster and graph evidence when emitted;
- does not enable optional read-ID/alignment outputs in the current research contract;
- never falls back to single-sample mode.

No real Severus binary qualification has been claimed yet.

## SAVANA

Paired and tumor-only SAVANA are separate policies and command paths.

The adapters retain independently:

- somatic SV evidence;
- CNA segments;
- selected purity/ploidy fit;
- ranked alternative fits;
- caller-native artifacts;
- sensitive read-support artifacts as non-exportable native evidence.

A valid CNA result can remain `COMPLETED` when the SV lane is `NO_CALL`; those states are not
collapsed into one biological interpretation.

The pinned SAVANA real-tool CI uses generated synthetic inputs only.

## Wakhan

The Wakhan phased-CNA adapter requires the registered phased input.

If breakpoint evidence is supplied, the adapter also requires the exact producing parent lane ID and
lane SHA-256. The multi-caller planner models the Severus parent as optional rather than pretending it
is independent supporting evidence.

The current Wakhan adapter is not marked real-tool-qualified. No analytical-validity claim follows
from its parser/runtime contract tests.

## ichorCNA

The ichorCNA lane is restricted to the registered cfDNA/ULP-WGS regime.

It binds and fingerprints:

- coverage WIG;
- GC WIG;
- mappability WIG;
- centromere resource;
- optional normal panel;
- pinned ichorCNA runtime/wrapper identity.

Structured output is additive to retained native output and includes available tumor-fraction,
ploidy, alternative solution and segment information.

Adaptive Sampling or arbitrary ONT BAM availability does not make this lane eligible. Any ONT
transfer study would require a new prospectively registered mode and validation program.

The current pinned ichorCNA 0.5.1 interoperability lane runs successfully on generated synthetic
cfDNA-like read-depth data.

## Dependency-aware comparison report

The comparison report is descriptive rather than adjudicative.

It preserves:

- planning decision and execution outcome independently;
- failed, unavailable, ineligible, no-call and not-assessable lanes;
- exact parent lane IDs and hashes;
- exact evidence collection addresses;
- native artifact references;
- conflicting numeric or textual observations without averaging them away.

The public report contract contains no winner, majority-vote or generic aggregate-value field.
Directly dependent caller evidence is labelled as dependent rather than independent confirmation.

## Public schemas

The public multi-caller contracts are exported through the normal ONTSeq schema generator:

- `schemas/multicaller-plan.schema.json`
- `schemas/multicaller-comparison-report.schema.json`

CI `scripts/export_schemas.py --check` is authoritative for schema freshness.

## Remaining work outside this integration branch

The following remain separate future validation work, not hidden TODOs inside the implemented
contracts:

- real-binary qualification for Severus;
- real-binary qualification for Wakhan;
- executable qualification for Spectre with Sniffles support;
- separate qualification for Sniffles2 mosaic;
- separate qualification for Severus single-sample;
- assay-specific analytical validation for every intended-use caller/mode;
- any clinical-validation or clinical-release program.

## Explicit non-claims

This branch does not:

- rank callers;
- select a winning caller;
- use majority voting as truth;
- infer biological truth from caller concordance;
- make a clinical diagnosis;
- authorize automatic ISCN release;
- establish analytical or clinical validity;
- add patient/genomic payloads to Git;
- merge PR #80;
- change the package version.

PR #80 remains Draft and unmerged. Package version remains `0.8.2`.
