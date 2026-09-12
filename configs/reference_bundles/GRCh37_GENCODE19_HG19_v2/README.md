# GRCh37.p13 / GENCODE19 plus exact UCSC hg19 Canonical-25

`GRCh37_GENCODE19_HG19_v2` is a new immutable reference authority. It retains the complete,
native GENCODE release-19 GRCh37.p13 assembly and annotation from the v1 bundle and adds a
separate analysis FASTA for BAMs aligned to the exact UCSC hg19 Canonical-25 dictionary.
Changing the reference bytes and generator contract under the v1 identity would be unsafe, so
v1 remains available as a legacy native-full recipe and all 0.6.2 GRCh37 profiles pin v2.

## Two explicit dictionaries

- `exact_full` uses `genome.fa`, the complete publisher GENCODE19 GRCh37.p13 assembly including
  patches, scaffolds and haplotypes.
- `grch37_ucsc_hg19_canonical_25` uses
  `analysis-reference/ucsc-hg19-canonical25.fa`, containing exactly chr1-22, chrX, chrY and
  chrM in that order. The nuclear sequences are streamed byte-for-byte from the pinned native
  GRCh37.p13 FASTA. `chrM` is taken from UCSC's original hg19 `chrM.fa.gz`: NC_001807,
  16,571 bases, publisher MD5 `dba65706c9bdd610f9aa3ebd950f2ac8` and pinned SHA256
  `d2368da512728003230b13a394b22b3197ffe368c6219cf62e5b6dfc4a5d4891`.
  The deterministic 3,095,697,110-byte derived FASTA is pinned as SHA256
  `1a814c7523f47e88d4fc1707a3bf9f18d51f2856a60477f56db13a833db3ca99`; its 25-row FAI is
  pinned as SHA256 `1fb875bbf510ef0c7bd5b9b92bd430afadccdd0b2bc2106d6ffce1d9f0f00339`.

UCSC documents that hg19 retained NC_001807 rather than replacing it with the later NC_012920
record. This is why the native GENCODE19 `chrM` (16,569 bases) cannot be used to analyze a UCSC
hg19 Canonical-25 BAM. ONTSeq validates the BAM against the selected profile's exact FAI and
passes the same FASTA to reference-dependent callers. It performs no BAM reheadering, runtime
liftover, alias conversion, implicit sequence substitution or cross-build fallback.

## Adaptive Sampling scope

Both GRCh37 Adaptive Sampling profiles use the separately checksum-pinned
`AML_AS_111_GRCh37_v1` target bundle. Its coordinate projection was performed once at build
time; no mapping tool or GRCh38 target coordinates are used during analysis. The profile and
report must retain the documented unresolved/review-required target gaps. Successful software
execution is not evidence of analytical or clinical validation.

## Installation and offline transfer

Run `ontseq references install GRCh37_GENCODE19_HG19_v2 --resource-root PATH` only after the
operator authorizes the public reference downloads. The source payload is approximately
1.03 GB and the expanded bundle contains both a full assembly and Canonical-25 FASTA; allow up
to 15 GB of working storage. Every source size and SHA256 is checked before generation, all
derived outputs are hashed into the activated manifest, and both FASTA/FAI pairs are validated.

For offline use, place all pinned source files at the recipe-relative `sources/` paths and use
the supported local-source/import workflow. Analysis after activation is local and offline.

## Attribution and licensing notes

- GENCODE release 19 is open-access; retain GENCODE/GRC attribution and applicable EMBL-EBI
  terms and original-contributor rights.
- The [UCSC hg19 chromosome directory](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/)
  identifies the assembly, documents the original mitochondrial record, publishes the MD5 list,
  and states that these chromosome files are freely available for public use. Retain UCSC/GRC
  attribution.
- UCSC database tables, ENCODE blacklist data and Hoffman Umap retain the same source-specific
  notices documented for v1. This reference bundle contains no `liftOver` executable or chain
  file. Any separately used coordinate-mapping software/data and laboratory panel design need
  their own institutional redistribution and commercial-license review.
- No copyrighted clinical classification manual or ISCN standard is bundled. This software
  remains RESEARCH USE ONLY and requires expert review.
