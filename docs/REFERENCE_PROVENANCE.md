# Reference sequence provenance

## Scope

This document defines a research-only, path-free reference provenance primitive for ONTSeq.
It addresses a specific safety gap in structural-variant analysis: matching only FASTA-index
contig names and lengths cannot prove that two references contain the same sequence bytes.
Two FASTA files can therefore have identical contig signatures while differing at the sequence
level.

## Implemented contract

`SequenceReferenceLock` records only:

- declared `reference_id`;
- declared genome build;
- SHA-256 of the FASTA bytes;
- SHA-256 of the FAI bytes; and
- SHA-256 of the ordered contig name/length signature derived from the FAI.

No reference path or sequence content is retained in the serialized lock.

`build_sequence_reference_lock()` creates the lock from a local FASTA/FAI pair.
`verify_sequence_reference_lock()` fails closed when the runtime FASTA, FAI, reference ID,
genome build, or contig signature differs from the precomputed lock.

The tests include a deliberate same-contig/same-length FASTA sequence substitution. The old
contig-signature-only check would not distinguish that condition; the sequence-level SHA-256
lock does.

## Scientific boundary

This is a provenance and reproducibility control only. A matching checksum does not validate a
reference build for clinical use, does not validate alignment or SV-calling performance, and does
not convert BND/TRA evidence into a confirmed fusion. `NO_CALL` remains distinct from a negative
result. All downstream ONTSeq biological outputs remain research-only and non-reportable until
assay-specific analytical validation and authorized human review are completed.

## Privacy boundary

The lock contains no patient data, read identifiers, sample identifiers, genomic reads, raw VCF
records, source paths, or inserted sequence. Reference FASTA/FAI files remain local and must not
be uploaded to GitHub, Monday.com, or external AI systems unless explicitly authorized.

## Next wiring step

The primitive is intentionally separate from the current cuteSV execution adapter. The next
bounded change is to require a precomputed `SequenceReferenceLock` at the cuteSV execution
boundary and verify it before any caller process is started. That wiring should be exercised in
the synthetic real-tool CI path before the previous contig-signature-only gate is retired.
