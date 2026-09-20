# GATK Mutect2 experimental research adapter

Research Use Only. Human Review Required. No automatic clinical release.

## Scope and status

This is an opt-in, separately executable module in the ONTSeq Python package. It does
not change the canonical CNV/SV registry, the planning-only multicaller catalogue,
Clair3's germline policy, the Desktop selection, clinical reports, or the consensus
engine. Those integrations require separate contracts and validation. A separate
caller on the same reads is correlated evidence, not orthogonal confirmation.

The adapter accepts GATK **4.6.2.0** only. This is an intentionally pinned supported
interface, not a claim that it is the latest release. Its version string is checked;
an executable/container digest is not verified by this module. Samtools' version
is recorded. The optional environment recipe is version-pinned, not a resolved
Conda lock or a tested container image.

Both tumor-only and tumor-normal commands are implemented. Every output remains
`reportable=false`, `analytical_validation=NOT_VALIDATED`, and each allele has
`somatic_status=NOT_VALIDATED`. `technical_pass` preserves filter semantics; it does
not establish somatic origin, clinical significance, a limit of detection, or assay
sensitivity. Even a matched normal does not remove this boundary.

## Execution path

1. Verify local input hashes, reference dictionaries, declared reference identity,
   BAM read-group SM/sample identity, coordinate-sort header and nonempty BED scope.
   Run samtools quickcheck. These checks do not prove the biological identity of a
   specimen, fully scan alignments, or independently establish the BAM's original
   alignment reference. The latter remains an explicit checksum-bound declaration.
2. Run Mutect2, retaining the unfiltered VCF and its `.stats` file.
3. When a population-sites resource is supplied, run GetPileupSummaries for the
   tumor and, if present, normal, followed by CalculateContamination.
4. Run FilterMutectCalls with the mandatory Mutect2 stats and optional contamination
   table. Never substitute the unfiltered VCF if filtering fails.
5. Import literal SNV, indel and MNV alleles, verify their REF against the locked
   FASTA, recheck input hashes, and write the separate evidence JSON.

No automatic BQSR, duplicate marking, orientation-bias model, HaplotypeCaller, SV/CNV
calling, resource downloads, PoN construction, or clinical annotation is performed.
The optional contamination path has not been qualified for this ONT assay.

## Inputs and prerequisites

Use Python >=3.11 on a local Linux/WSL runtime with GATK 4.6.2.0 and samtools available.
Use only institutionally approved local files; never commit manifests, logs or results
from real specimens. Local output can contain specimen IDs and genomic variants.

Required: uncompressed reference FASTA, adjacent `.fai`, same-basename `.dict`, a
coordinate-sorted single-sample `.bam` and adjacent `.bam.bai`, and a nonempty `.bed`
with zero-based half-open intervals. Sample names must match BAM read-group `SM`
exactly. All paths are absolute and all files carry SHA-256 locks.

Optional normal: a different BAM, sample name and checksum; its declared reference
must equal the tumor's locked reference. In hematology, selection and verification
of a genuinely suitable constitutional normal is an upstream expert responsibility.

Optional VCF resources: bgzip-compressed `.vcf.gz` with adjacent `.tbi`. Each declares
the exact reference SHA-256. `germline` is a population allele-frequency resource;
`population_sites` is a separate common biallelic SNP/AF resource appropriate for
GetPileupSummaries. They are not interchangeable by name. `panel_of_normals` must
explicitly declare `technology="ONT"` and an `assay_id` identical to `profile_id`.
These declarations are checked, not independently certified. Resource VCF integrity,
AF semantics and dictionary interpretation also rely on GATK's native checks and
upstream resource qualification; dictionary validation is not disabled.

Absent resources are recorded as warnings, not substituted with imaginary defaults.
Without population sites, contamination is NOT_ESTIMATED, never assumed zero.
An empty object, false, empty string or array is not an acceptable optional resource:
use a populated object or null.

## Prepare a locked manifest locally

The example deliberately requires real approved paths rather than inventing checksums.
Adapt the paths, profile ID and exact SM value to your local assay. Confirm alignment
reference identity from the upstream run manifest before binding its checksum here.
Run this from a checkout containing the adapter, with `PYTHONPATH=src` or an installed
ONTSeq package.

```python
import json
from dataclasses import asdict
from pathlib import Path
from ontseq_platform.gatk_mutect2 import GATKConfig, LockedFile, Sample, VCFResource, digest


def locked(filename: str) -> LockedFile:
    path = Path(filename).expanduser().resolve(strict=True)
    return LockedFile(str(path), digest(path))


reference = locked('/approved/reference/GRCh38.fa')
config = GATKConfig(
    profile_id='ONT_AS_research_v1',
    reference=reference,
    reference_fai=locked('/approved/reference/GRCh38.fa.fai'),
    reference_dict=locked('/approved/reference/GRCh38.dict'),
    tumor=Sample(
        name='TUMOR_001',  # Must equal read-group SM, not merely the BAM filename.
        bam=locked('/approved/input/tumor.bam'),
        index=locked('/approved/input/tumor.bam.bai'),
        reference_sha256=reference.sha256,
    ),
    intervals=locked('/approved/reference/targets.bed'),
)
# For paired mode, set normal=Sample(...) with the verified normal's exact SM.
# For resources, use VCFResource(vcf=locked(...), index=locked(...),
# reference_sha256=reference.sha256). For a PoN also supply technology='ONT'
# and assay_id=config.profile_id. Missing resources remain explicit warnings.
Path('gatk.local.json').write_text(json.dumps(asdict(config), indent=2) + '\n')
```

Preview commands only; this is NOT_RUN, not successful preflight or execution:

```bash
PYTHONPATH=src python -m ontseq_platform.gatk_mutect2 \
  --config gatk.local.json --output-dir /approved/results/gatk-run-001
```

Execute explicitly; the output directory must not already exist:

```bash
PYTHONPATH=src python -m ontseq_platform.gatk_mutect2 \
  --config gatk.local.json --output-dir /approved/results/gatk-run-001 \
  --execute --allow-experimental-ont
```

A caller failure stops the run. Do not retry into the same directory or manually copy
an old VCF into it. Use a fresh run directory after resolving the failure.

## Evidence and failure semantics

Outputs include `config.lock.json`, a planning-only `command-plan.json`, numbered
command logs, native unfiltered/filtered VCFs and stats, optional contamination
artifacts, and `evidence.json`. The evidence records executed argv, versions, profile,
reference/input/config/artifact hashes and caveats. `failure.json` records a failed
attempt; no successful evidence is published from a failed run.

`COMPLETED` means at least one allele met native technical filter criteria. `NO_CALL`
means no technically passing allele; this is NOT a negative clinical finding. A
valid header-only VCF is NO_CALL; malformed VCF or missing output is FAILED.

All failed FILTER records remain available. `FILTER=.` is not PASS. At multiallelic
sites, allele-specific `AS_FilterStatus` is respected; absent allele-level filtering
is not enough to promote multiple ALT alleles. `caller_af` preserves FORMAT/AF,
whereas `ad_fraction` is ALT AD divided by sum(AD). FORMAT/DP is retained separately.
Missing values remain null and QUAL is allowed to be missing.

The evidence uses native one-based VCF coordinates. It is NOT left-normalized and
must not be joined to other callers without representation-aware normalization.
Symbolic ALT/gVCF records are outside this reader's scope and fail import rather
than disappearing silently. Indels carry an orthogonal-confirmation requirement;
all variants still require expert review.

## Verification and validation impact

At initial authoring, 58 focused tests passed locally under Python 3.13 using
synthetic fixtures and injected fake process outputs. They cover both modes, resource
routing, exact-version gating, provenance, malformed inputs, stale results, failed
processes, FILTER/AS semantics, quantities, REF checks and mid-run input changes.
These are adapter contract tests, NOT real GATK, real BAM, or biological tests.

The complete repository safety/version/lint/type/test commands must pass in CI before
merge. The local authoring workspace did not contain the full repository or installed
ruff/mypy. The dedicated GATK contract workflow complements, not replaces, normal CI.

Before runtime qualification: execute real pinned GATK/samtools on deliberately
synthetic indexed BAMs in both modes, including valid resources and failure paths.
Before analytical qualification: lock the assay/basecaller/chemistry, resources,
truth regions and normalization, then measure sensitivity/PPV/FP burden by depth,
VAF, variant class and context using held-out truth. No LoD95 or clinical performance
claim follows from the present tests. Existing CNV/SV results are unchanged.

## Primary interface references

- Mutect2: https://gatk.broadinstitute.org/hc/en-us/articles/35967636854939-Mutect2
- FilterMutectCalls: https://gatk.broadinstitute.org/hc/en-us/articles/4404607509403-FilterMutectCalls
- Pinned upstream release: https://github.com/broadinstitute/gatk/releases/tag/4.6.2.0
- ONT-focused somatic comparator family, not implemented by this adapter:
  https://github.com/HKU-BAL/ClairS

The upstream interfaces motivate the command plan; they do not validate this ONT assay.
