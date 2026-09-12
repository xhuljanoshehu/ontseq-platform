# Native GRCh37.p13 / GENCODE19 reference lock

This is an installable source recipe, not a relabelled GRCh38 cache or a liftover result.
Every source was downloaded from its original provider on 2026-09-02, measured, SHA256-hashed
with Linux `sha256sum`, and independently rehashed with PowerShell `Get-FileHash`. Large source
files and generated databases are intentionally outside Git. Dates marked 2026-09-02 identify
the verified snapshot, not a claim that GENCODE19 was released on that date.

The original GENCODE release-19 `MD5SUMS` additionally agrees with both downloaded gzip streams:

| Publisher stream | Publisher and observed MD5 |
| --- | --- |
| `GRCh37.p13.genome.fa.gz` | `490904ef2561472d0801c704ffc85761` |
| `gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz` | `878f8d5c6d0e6603bcb2e4b61a6262a5` |

## Scope and non-equivalence

- [GENCODE release 19](https://www.gencodegenes.org/human/release_19.html) labels the FASTA and
  ALL-region GTF as GRCh37.p13, including scaffolds, patches and haplotypes. This is not the
  GRCh38 recipe's primary-only assembly, and not an arbitrary pre-existing hg19/B37 BAM dictionary.
- The profile `AML_LCWGS_GRCh37` requires the exact full publisher dictionary. Mismatches must
  fail visibly. No alias conversion, coordinate liftover or reference guessing is performed.
- This v1 authority contains no exact UCSC hg19 Canonical-25 analysis FASTA. It must not be used
  for a BAM whose dictionary has `chrM=16571`; doing so would give reference-dependent callers
  the incompatible native 16,569-base mitochondrial sequence. ONTSeq 0.6.2 profiles therefore
  pin the new immutable `GRCh37_GENCODE19_HG19_v2` bundle. No v1 bytes are silently changed.
- Native GENCODE19 does not supply native GRCh37 MANE. The cache declares that absence and ranks
  only the available native transcript evidence. No GRCh38 MANE coordinates are substituted.
- UCSC [hg19 database documentation](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/)
  notes that its mitochondrial sequence differs from the revised Cambridge sequence. Therefore
  UCSC RepeatMasker/simple-repeat/segmental-duplication and ENCODE blacklist derivatives are
  restricted to chr1-22/X/Y. The original files, including excluded rows, remain checksum-pinned
  for provenance. No mitochondrial, patch, scaffold or haplotype artifact context is inferred.
- The [Hoffman Umap source](https://bismap.hoffmanlab.org/) is a corrected zero-based hg19 k100
  single-read BED containing only chr1-22/X/Y. It marks uniquely mappable 100-mer regions, not
  nanopore analytical sensitivity. Missing zero-value rows cannot establish negative evidence.
- Adaptive Sampling and the UCSC hg19 Canonical-25 dictionary are published only through the
  v2 reference family. This legacy v1 recipe remains available for provenance and exact-full
  native GRCh37.p13 compatibility; it is not the 0.6.2 managed Desktop reference.
- `HEMATOLOGY_GRCh37_v1` contains gene-symbol-only review knowledge derived from the existing
  source-attributed HEMATOLOGY_v3 list. Its original panel-limited selection remains documented;
  matching depends on native annotation and does not create coverage or clinical validity.

## Installation and offline transfer

Run `ontseq references install GRCh37_GENCODE19_v1 --resource-root PATH` only for an explicitly
requested legacy/native-full installation. The managed 0.6.2 Desktop path installs
`GRCh37_GENCODE19_HG19_v2`, which also provides the exact hg19 analysis reference. Existing
GRCh38 resources are not replaced.

An activated bundle can be copied to another approved machine and imported with
`ontseq references import PATH --resource-root RESOURCE_ROOT`. The official ID remains bound
to the shipped authority recipe, not a caller-supplied configuration. Analyses after activation
are offline. Use `ontseq references validate --resource-root RESOURCE_ROOT` to verify installed
resources. Engineering success is not analytical or clinical validation.

## Source attribution and terms

These notes do not add a license to ONTSeq or authorize redistribution of the application.
Institutional release approval and third-party checks remain separate.

- GENCODE identifies its data as [open access](https://www.gencodegenes.org/pages/data_access.html).
  Retain GENCODE/GRC attribution, release identifiers and [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/),
  including attribution and any rights retained by original contributors.
- UCSC's hg19 database directory states that its files and tables are freely usable; retain the
  track/provider attribution, including RepeatMasker and segmental-duplication sources.
- The separately used UCSC chain and `liftOver` program are build-time provenance inputs and are
  not redistributed by ONTSeq. UCSC permits non-commercial use but requires separate commercial
  licenses; its published terms do not explicitly settle redistribution of derived coordinate
  BEDs. Institutional/open/commercial distribution must obtain written clarification as needed,
  and separately confirm rights in the original laboratory target design.
- ENCODE `ENCFF001TDO` remains attributed to the ENCODE exclusion-region resource. Consult the
  [ENCODE data-access policy](https://www.encodeproject.org/about/data-access/) before reuse or
  redistribution; this artifact-context layer is not a validated ONT clinical filter.
- Umap remains attributed to Karimzadeh et al. (2018), [doi:10.1093/nar/gky677](https://doi.org/10.1093/nar/gky677).
  A separate redistribution license for the specific downloaded BED was not established in this
  integration review; verify its applicable data terms before institutional redistribution.
  Software-package licenses are not assumed to be licenses for third-party datasets.
- Existing hematology source attribution and review-only restrictions are retained byte-for-byte
  for every knowledge record. No copyrighted classification manual or ISCN standard is bundled.
