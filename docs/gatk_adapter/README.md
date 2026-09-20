# Experimental GATK / Mutect2 adapter for ONTSeq

Adapter schema: 0.1.0. Prepared 2026-09-20.
Inspected upstream base: `xhuljanoshehu/ontseq-platform` at
`2c7aa2bf1a2d01f3f51aa75232c608407f8329b9` (package 0.8.2).

**Research Use Only. Human review required. No automatic clinical release.**

## Delivery boundary

This is an additive source patch implementing an opt-in Python API and standalone CLI.
It does **not** modify the canonical ONTSeq runner, its multi-caller registry, public
CNV/SV result schema, consensus engine, HTML/Excel report, desktop interface, or release path.
It has **not** been committed, pushed, merged, or tested in the full upstream repository.
A complete repository checkout was not available in this execution environment.

The Python adapter is implemented; native GATK interoperability and analytical validity
are **not established**. Contract tests use deliberately artificial text stand-ins for
BAMs, indices, resources and a JAR, plus a synthetic process executor. These stand-ins
are not valid sequencing files and must never be used to demonstrate native GATK success.
One separate transport test executes a harmless Python process, not GATK.

## Implemented flow

1. Require both explicit enablement and acknowledgement of experimental ONT use.
2. Create a fresh, private output directory. Never reuse an earlier run directory.
3. Check nonempty local inputs, SHA-256 identities, conventional index paths, reference
   dictionary/FAI agreement and declared reference/assay identity.
4. Probe the requested Java major and exact GATK version.
5. Run GATK `ValidateSamFile` and `GetSampleName` for each BAM. Compare native read-group
   sample names with the manifest. Do not silently substitute a different sample.
6. Run `Mutect2` in explicitly selected tumor-only or paired tumor-normal mode.
7. Retain the raw VCF, its index and the mandatory Mutect2 statistics sidecar.
8. Optionally run `GetPileupSummaries` and `CalculateContamination`, including the matched
   normal when present. Check the emitted contamination estimate before use.
9. Run `FilterMutectCalls`, with explicit statistics and optional contamination/segmentation.
10. Recheck all input hashes. Extract literal small-variant alleles from the filtered VCF
    into a separate research schema. Retain every FILTER state, including `.`.
11. Write structured status, commands, step exit codes, checksums and native artifacts.

`NOT_RUN` means disabled, not a biological negative. `FAILED` exposes no normalized
candidates. `COMPLETED` means this adapter's steps and parser completed, not that the
assay is validated or that a zero-candidate result is biologically negative.
`execution_evidence` records NONE, TEST_DOUBLE or NATIVE_EXECUTION; the latter means
native execution was attempted, not that interoperability qualification was passed.

## Compatibility target and runtime identity

The explicit compatibility target is **GATK 4.6.2.0 / Java 17**. This is not a statement
that 4.6.2.0 is the newest release, safest available version, or qualified ONT runtime.
A different GATK version is rejected rather than silently accepted. Re-targeting requires
review of arguments/output contracts and a real-tool qualification run.

The operator supplies a local `gatk-package-4.6.2.0-local.jar`. Its actual SHA-256 is
recorded and checked, and its runtime-reported version must match. Java's major version
and output log are recorded; this is **not** a fully immutable JVM/OS/container lock.
No GATK binary, Docker image, Conda environment, genomic resource, or network downloader
is bundled. No Docker digest or Conda solve is claimed. Obtain software and resources
from institutionally approved sources and independently verify provenance, terms and
security before supplying them. Merely computing a hash does not authenticate a source.

## Inputs and manifest meaning

Supported initial input is a local indexed BAM, not FASTQ, signal/POD5, SAM or CRAM.
Reference FASTA must be uncompressed with conventional `.fai` and `.dict` sidecars.
Population allele-frequency resources and optional ONT PoN must be indexed `.vcf.gz`
files with `.tbi` sidecars. BGZF/resource-content correctness is left to the native tools;
Python preflight does not certify it. No alternate unregistered BAM index is accepted.

A germline population resource is required by this adapter's policy, although this does
not assert that every possible Mutect2 invocation requires one. The optional PoN must
be declared ONT-specific and have the same assay ID. A typed label is an operator
attestation, not proof of matching chemistry, library preparation or basecaller.

Tumor and normal must have different sample names, BAM paths and hashes. A normal's
biological suitability cannot be established by these checks. Reference/assay labels
provided to `lock` are declarations: matching contig names and hashes do not prove the
historical alignment reference or the biological relationship between specimens.

Contamination estimation is optional but never silently assumed to be zero. It requires
an explicit population-sites resource. Insufficient eligible SNP sites may cause native
failure; a missing/non-finite estimate is not converted into a clean sample. With both
sites and target intervals, the pileup command uses their intersection.

## Candidate semantics

Each candidate retains run/sample identity, source VCF hash, original line identity,
allele index, literal REF/ALT, native site filter, unsplit `AS_FilterStatus`, and available
FORMAT/AF, AD and DP. AF is taken from the caller, not reconstructed from AD/DP and not
interpreted as tumor purity or probability of somatic origin. Missing values stay null.

Coordinates are zero-based, half-open **native reference spans including the VCF anchor**.
There is no left-normalization, reference-base validation of REF at every locus, or
cross-caller allele equivalence matching in this deliverable. Do not compare indels by
string equality or add this output to a vote-based consensus. Multi-ALT records retain
the complete native allele-filter string; it is not misrepresented as one allele's PASS.

Symbolic/spanning deletion alleles and unsupported records fail the entire extraction;
they are not silently dropped. Native VCF remains available for review. The parser has
candidate and per-line bounds but is intended for bounded research callsets; WGS-scale
memory/runtime performance has not been measured.

`clinically_reportable` and `analytically_validated` cannot be set true in the result
contract. Caller concordance is not independent/orthogonal evidence of biological truth.
Tumor-only candidates do not establish somatic versus germline origin.

## Deliberately absent behavior

No BQSR, duplicate marking/removal, short-read orientation-bias model, disabling of native
read filters, automatic sample relabeling, tumor-only fallback from paired mode, caller
voting, clinical classification, annotation download, cloud upload, or automatic release.
No report is silently populated with the experimental candidates.

## Use

See `examples/gatk/README.md` for local command examples. The package is discovered by
ONTSeq's existing `src` package layout, but this patch does not register an `ontseq gatk`
command. Its explicit CLI is `python -m ontseq_platform.gatk_adapter`.

Python API:

```python
from pathlib import Path
from ontseq_platform.gatk_adapter import GATKAdapter, Mutect2Config

config = Mutect2Config.model_validate_json(Path("gatk.config.json").read_text())
steps = GATKAdapter(config).plan(Path("/secure/research/run-new"))  # No execution.
# .run(...) remains NOT_RUN unless both manifest opt-ins are true.
```

All config/plan/native/result files can contain sample identifiers, local paths or genomic
information. Keep them outside Git and in approved local storage. Native processes use
shell-free argument vectors and restrictive output permissions on POSIX. Windows use is
through Linux/WSL; native Windows interoperability has not been tested. These controls do
not sandbox GATK or prove that an operator-supplied executable cannot access the network.

## Primary references used for interface review

- GATK source at the pinned tag:
  https://github.com/broadinstitute/gatk/tree/4.6.2.0
- Mutect2:
  https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2.java
- FilterMutectCalls (statistics, filtering and contamination interfaces):
  https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/filtering/FilterMutectCalls.java
- CalculateContamination:
  https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/contamination/CalculateContamination.java
- Upstream project guardrails:
  https://github.com/xhuljanoshehu/ontseq-platform/blob/2c7aa2bf1a2d01f3f51aa75232c608407f8329b9/AGENTS.md

Source inspection is not execution evidence. No third-party source code or genomic data
has been copied into this delivery. This patch adds no license or public deployment.
