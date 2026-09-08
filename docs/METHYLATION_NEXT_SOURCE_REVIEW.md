# Public data that can support the next methylation test

Reviewed 2026-09-05. **Public test candidates exist; none has passed this
repository's biological recovery experiment.** This follow-up reviews public
documentation and existing source catalogues. It does not inspect new genomic
file contents, acquire data or claim an executable admitted cohort.

## Recommended next source pair: HG008 normal tissues

Use HG008-N-P (normal pancreas) and HG008-N-D (normal duodenum) as a candidate
two-source tissue experiment, separate from GSE210260. The matched ONT-std-1
protocol specifies PromethION R10.4.1, SQK-LSK114 and Dorado 0.5.3 with the
v4.3.0 5mC/5hmC model. The
[primary data descriptor](https://www.nature.com/articles/s41597-025-05438-2)
documents this protocol. Both modification channels must be preserved.

The tissues share one donor and contain mixed cell populations. Public genomic
dissemination is expressly supported by the
[NIST resource description](https://www.nist.gov/programs-projects/cancer-genome-bottle);
this is not evidence of anonymity or a local institutional authorization.
This pair can test scoped tissue-source recovery, without establishing purified
cell-type, donor-independent, physical DNA-mass or tumour-fraction recovery.

The existing [exact-object catalogue](../configs/methylation/public_data_candidates.v1.json)
records a 49,607,859,349-byte pancreas `_2` uBAM and a 77,856,391,489-byte duodenum
`_3` uBAM: **127,464,250,838 bytes total**, excluding reference, derivatives and
scratch space. These are previous HTTP HEAD measurements, not newly downloaded
or locally verified BAMs. Supplier MD5s are catalogued; local SHA-256 remains null.
The catalogue also identifies the February 2025 modification-preserving GRCh38
alignment products and the matching reference URL from the
[supplier README](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern_ONT-std_20240422/README_NE_ONT-std.md).

The uBAM pair needs a common local alignment before this repository's aligned
adapter. A bounded whole-read sample must have its algorithm, domain, budget and
seed fixed first. A region-limited extraction is a different, region-conditional
experiment; the current code does not implement an audited remote extraction.
Complete WGS TSV exports are unsuitable for the current bounded in-memory parser.

## Smaller full-file alternative: COLO829 / COLO829BL RRMS

The [official ONT RRMS release](https://epi2me.nanoporetech.com/rrms2022.07/)
provides per-read modification BAMs for a melanoma/B-lymphoblast pair, with five
technical replicates per source and a common adaptive-sampling design. It also
provides RRBS measurements. These do not establish physical dilution truth.

The [existing exact-object catalogue](../configs/methylation/additional_data_candidates.v1.json)
selects COLO829_3 plus COLO829BL_2 at **44,670,539,944 bytes**; these are previously
reviewed S3 metadata sizes. The Bonito/Remora versions, complete model, library
kit, aligned reference and source-specific filtering remain unresolved. The
current Dorado-only provenance contract cannot simply relabel this older caller.
This option needs its own compatibility review before execution. Technical
flowcells do not supply new donors.

## Already local small data: GSE218679

The two acquired files total only **208,521 compressed bytes** and have verified
SHA-256 values in the
[acquisition catalogue](../configs/methylation/acquired_public_sources.v1.json).
They remain useful for intake/parser checks. The
[recorded structural intake](METHYLATION_INTAKE_STATUS.md) found insufficient
holdout pools for the frozen matrix; exact pore, Guppy configuration and reference
artifact provenance are still missing. The existing targeted FMR1 `v1` matrix
remains unchanged and its admission remains NO_CALL.

GSE265944 was an additional metadata search lead, but a new compatible independent
pair was not established in this review. It is not an acquisition recommendation
or a substitute cohort; a new accession alone cannot rule out reused sources.

## Concrete remaining preparation

HG008 is preferred for the new Dorado path because its processing is best
documented among these candidates. Before acquiring it, fix the exact input list,
ignored local destination and a storage/sampling plan; retain the existing
owner/institutional transfer boundary in AGENTS.md. Check current sizes, revision
and supplier checksums, then compute actual complete-file SHA-256 after acquisition.
Locally verify reference and BAM integrity, headers, modification tags and input
scope before admitting sources. Metadata cannot certify file contents, malware
absence or biological validity. No new file-content inspection was performed here.

The inspected Windows environment lacks pysam, samtools and modkit; binary
interoperability is therefore still pending. Once runtime and sources are ready,
freeze the actual cohort, software, sampling domain and dilution matrix before
test outcomes. Existing count budgets are read groups, not genomic coverage.
Use held-out reads and retain every planned PASS/FAIL/NO_CALL outcome. Only a
successful scoped recovery study supports proceeding to independent DNA-mass
reference calibration; no current candidate establishes a 10–20% detection limit.

Current implementation and verification: [follow-up status](METHYLATION_SOURCE_REVIEW_STATUS.md).
