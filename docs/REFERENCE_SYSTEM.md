# Build-locked reference, panel and knowledge system

ONTSeq activates resources from manifests, not from filenames. Version 0.6.2 introduced four
GRCh38 profiles plus four build-isolated GRCh37 profiles (lcWGS and Adaptive Sampling, each with
exact-full and UCSC-hg19 Canonical-25 dictionary contracts). This document details the shared
reference/panel system; the GRCh37 implementation is documented in
[`GRCH37_INTEGRATION.md`](GRCH37_INTEGRATION.md). Neither family has a run-time liftover,
cross-build or dictionary-contract fallback.

## Resource root

The Core/CLI default root is `/opt/ontseq`; `ONTSEQ_RESOURCE_ROOT` or `--resource-root` may
select a different absolute root. Fresh Desktop settings use
`~/.local/share/ontseq/resources-v0.7.1`, derived from the Desktop/Core software version.
Previously saved resource paths remain selected, including unversioned and custom paths.
The release suffix names an installation directory; it does not update scientific contents.
Reference, panel, knowledge and profile versions remain separately pinned by their manifests.
One GRCh38 family layout is:

```text
/opt/ontseq/
  references/GRCh38_GENCODE50_MANE1.5_v1/bundle.yaml
  panels/AML_AS_111_GRCh38_v1/bundle.yaml
  knowledge/HEMATOLOGY_v3/bundle.yaml
  profiles/AML_LCWGS_GRCh38.yaml
  profiles/AML_AS_111_GRCh38.yaml
  profiles/AML_LCWGS_GRCh38_CANONICAL25.yaml
  profiles/AML_AS_111_GRCh38_CANONICAL25.yaml
```

The GRCh37 family mirrors this manifest layout with
`references/GRCh37_GENCODE19_HG19_v2`, `panels/AML_AS_111_GRCh37_v1`,
`knowledge/HEMATOLOGY_GRCh37_v1` and four `AML_*_GRCh37*.yaml` profiles.

Only a directory with a valid `bundle.yaml` is discoverable. Every active file has a size and
SHA256, every coordinate-bearing bundle declares one build, and a profile pins exact bundle IDs.
Loose files are ignored.

To deliberately rename an existing installation, first finish active analyses and close
services using that root. Record its bundle identities and verify its checksums. Move the
complete directory only to an explicitly selected, unused local destination, then run
`ontseq references validate --resource-root NEW_ABSOLUTE_ROOT` and check that the required
profiles resolve before saving the new Desktop path. Preserve the old location's usability
for historical provenance, for example with a compatibility link to the new location, and
do not rewrite historical result files. A new-name link to the old directory alone is not
a canonical rename: Core resolves symlinks before recording the resource root and paths.
There is no automatic migration, redownload or scientific bundle-version change merely
because the Desktop/Core version advances.

## Reference lifecycle

```bash
ontseq references list
ontseq references status --resource-root /opt/ontseq
ontseq references validate --resource-root /opt/ontseq
ontseq references install GRCh38_GENCODE50_MANE1.5_v1 --resource-root /opt/ontseq
ontseq references repair GRCh38_GENCODE50_MANE1.5_v1 --resource-root /opt/ontseq
ontseq references import /path/to/pinned-reference-tree --resource-root /opt/ontseq
```

Install downloads only in an explicit lifecycle command. Sources land in a private staging
directory, are checked before decompression, and the FASTA/FAI/full canonical dictionary is
validated before activation. GENCODE, MANE and cytobands compile into a checksum-pinned SQLite
cache. The panel Analysis ROI and transcript cache compile from that SQLite cache; all four
profiles are exposed only after reference, knowledge and panel validation succeeds.

The complete install is intentionally a multi-gigabyte one-time operation. Reserve at least
15 GiB of free disk for sources, staging, SQLite/VACUUM and the activated bundle, and 6 GiB of
available memory for GENCODE compilation. A full-release Windows smoke peaked near 4.3 GiB of
process memory; normal analyses do not recompile or rescan the GTF.

Both dictionary contracts use that same installed resource family. Canonical-25 is derived from
the pinned full `ReferenceLock` at run resolution; it is not a second FASTA, bundle installation
or coordinate system. Selecting a `*_CANONICAL25` profile therefore performs no additional
multi-GiB download and consumes no second reference installation.

`status` performs a fast presence/size check. `validate` computes full checksums. `repair`
retrieves or regenerates missing and invalid reference artifacts, then stages and repairs the
complete pinned Knowledge/Panel/profile family with rollback. It refuses a silent bundle-version
change or a changed Source-/Generator contract under the same version. `import` supports an
already activated bundle or a local pinned recipe/source tree and does not use the network. When
an imported ID is present in the packaged catalog, its Source-/Generator declarations must match
that authority exactly; unrelated custom bundle IDs remain importable.

The packaged `GRCh38_GENCODE50_MANE1.5_v1` recipe contains exact publisher byte sizes and
SHA256 values. GENCODE FASTA/GTF values were cross-checked against the publisher MD5 index as an
independent transfer check. The installer still refuses unresolved or malformed hashes; guessed
or MD5-substituted SHA256 locks are never accepted.

## Compiled annotation cache

The cache contains `metadata`, `genes`, `transcripts`, `exons`, `cds` and `cytobands` plus
chromosome/start/end indices. GTF/GFF coordinates are converted exactly once from 1-based
inclusive to 0-based half-open. Analysis code queries SQLite and does not rescan GENCODE.

Transcript preference is deterministic: MANE Select, MANE Plus Clinical,
canonical/APPRIS-principal, protein-coding/basic, then all remaining transcripts; CDS length,
transcript length and ID break ties.

Large RepeatMasker, simple-repeat, segmental-duplication, blacklist and mappability BEDs remain
separate checksum-pinned resources. At SV startup ONTSeq validates them in one streaming pass and
records only contig byte ranges and row counts. Breakpoint queries then lazy-load at most two
queried contigs per resource into compact coordinate arrays and use bisect/block-max pruning;
millions of genome-wide rows are not retained as Python interval objects or rescanned for every
breakpoint.

## Panel roles

`selection_panel_buffered` is the normalized 111-interval sequencing design. It remains distinct
from `analysis_roi_unbuffered`, which contains only uniquely resolved GENCODE gene bodies.
Adaptive Sampling measures both and writes separate coverage reports. `IGH_REVIEW_REQUIRED` has
no fabricated ROI and cannot produce a negative observability statement. See
[`PANEL_PROVENANCE.md`](PANEL_PROVENANCE.md).

The dedicated GRCh37 panel has the same two roles but a separate immutable authority. Its
selection file contains 110 of the 111 intended laboratory intervals after a checksum-pinned
build-time projection at `minMatch=0.99`; 109 are reciprocal-exact and its native GENCODE-19 ROI
contains 108 targets. `ACACA_ROUNDTRIP_REVIEW_REQUIRED`, `CT45A2_REVIEW_REQUIRED`,
`IGH_REVIEW_REQUIRED` and `GPR128_REVIEW_REQUIRED` remain explicit.
The projection chain and tool are provenance inputs only: they are neither installed nor invoked
by an analysis. No GRCh38 panel file is resolved by a GRCh37 profile.

## Profile analysis

```bash
ontseq analyze SAMPLE_GRCH38.bam --profile AML_LCWGS_GRCh38 \
  --resource-root /opt/ontseq
ontseq analyze SAMPLE_GRCH38.bam --profile AML_AS_111_GRCh38 \
  --resource-root /opt/ontseq
ontseq analyze SAMPLE_GRCH38_CANONICAL25.bam \
  --profile AML_LCWGS_GRCh38_CANONICAL25 --resource-root /opt/ontseq
ontseq analyze SAMPLE_GRCH38_CANONICAL25.bam \
  --profile AML_AS_111_GRCh38_CANONICAL25 --resource-root /opt/ontseq
ontseq analyze SAMPLE_GRCH37.bam --profile AML_AS_111_GRCh37 \
  --resource-root /opt/ontseq
ontseq analyze SAMPLE_HG19_CANONICAL25.bam \
  --profile AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25 --resource-root /opt/ontseq
```

The command locates `sample.bam.bai` before `sample.bai`, reads the full BAM sequence dictionary,
detects a complete consistently named GRCh38 assembly, and applies the selected profile's exact
dictionary contract:

- `AML_LCWGS_GRCh38` and `AML_AS_111_GRCh38` use `exact_full`. The BAM must match the complete,
  ordered Primary-Assembly `ReferenceLock`, including all 194 entries in that installed lock.
- `AML_LCWGS_GRCh38_CANONICAL25` and `AML_AS_111_GRCh38_CANONICAL25` require exactly
  `chr1`-`chr22`, `chrX`, `chrY`, `chrM`, with standard GRCh38 lengths and order. Missing,
  additional or reordered entries fail.

The contracts never fall back to one another. GRCh37, partial, mixed-style, reordered,
length-mismatched or otherwise non-matching dictionaries stop before the pipeline is created.
All four profiles pin the same `GRCh38_GENCODE50_MANE1.5_v1`, `HEMATOLOGY_v3` and, for Adaptive
Sampling, `AML_AS_111_GRCh38_v1` bundles. There is no coordinate conversion or liftover.

The two GRCh37 lcWGS/Adaptive-Sampling pairs follow the same dictionary separation:
unsuffixed profiles require the complete GENCODE-19 publisher dictionary, while both
`*_UCSC_HG19_CANONICAL25` profiles require exactly the ordered UCSC-hg19 25-contig dictionary.
The Adaptive-Sampling pair additionally pins `AML_AS_111_GRCh37_v1`. Panel/reference build and
dictionary compatibility is checked by the registry before the profile can resolve.

`HEMATOLOGY_v3` is a panel-scoped expert-review resource. It combines locally curated AML
patterns with exact, accepted-evidence CIViC fusion features that overlap the panel and have at
least one Disease Ontology association below `hematologic cancer` (`DOID:2531`). The checked-in
release contains 38 patterns and 74 pathology associations. Reports retain the pathology name,
DOID and source record, but a database match never changes `reportable` and never establishes a
productive transcript, diagnosis or prognosis. The deterministic maintainer tool
`scripts/build_hematology_knowledge.py` rebuilds a new version from explicitly supplied,
checksum-recorded source snapshots; analyses themselves remain offline.

The current result contract is `PipelineResult 0.3.0`. It records the resolved bundle context,
releases, checksums, large-table sidecars and the structured ISCN proposal state. Versions 0.1.0
and 0.2.0 remain readable and retain `LEGACY_UNSPECIFIED` ISCN semantics; they are not silently
upgraded. An assessed ISCN proposal additionally requires its reference lock, cytobands and
annotation-cache checksums to match this resolved context. lcWGS target coverage and
Adaptive-Sampling observability are `NOT_APPLICABLE`, not negative calls.

## Offline operation and updates

After activation, analyses read only the resource root and perform no network request. To update a
publisher release, create a new bundle ID/version and install it beside the old bundle; do not
edit an activated manifest. A profile update is a separate, reviewable change that pins the new
bundle. Existing result provenance therefore continues to resolve to the exact prior files.

GRCh37 is implemented as a separate bundle/profile/test series, including the complete
GRCh37.p13 dictionary and the UCSC hg19 Canonical-25 contract documented above. Its resources
must never be copied into GRCh38 directories or manifests; no profile performs an automatic
liftover or build fallback.
