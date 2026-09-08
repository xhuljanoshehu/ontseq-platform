# Dorado MM/ML and modkit per-read input

Research use only. This adapter is a technical implementation, not a validated methylation
assay. `src/ontseq_platform/modbam.py` retains intact read groups for the existing source-A
call-rate estimator; it does not infer DNA mass or tumour fraction. The Nanopolish adapter
remains separate and unchanged in its biological interpretation.

## Supported representation

`parse_modbam_source` accepts aligned SAM or BAM with Dorado MM/ML tags. SAM is decoded using
the Python standard library; BAM additionally requires a local `pysam` runtime. CRAM, duplex
opposite-strand MM groups, paired-end alignments and unsupported cytosine modifications fail
closed. The directly supported modification models emit `C+m` or `C+mh`/equivalent separate
m/h groups. The caller and modification model are required declared provenance.

Adapter version `0.1.1` also checks the per-record Dorado `dx` tag independently of
program-header recognition. Explicit duplex offspring (`dx:i:1`) are rejected;
`dx:i:0` and `dx:i:-1` identify simplex reads and remain supported. Malformed or
repeated `dx` tags fail. A simplex parent can share a physical molecule with
another simplex parent, so retained read groups still do not establish molecule
independence. Missing `dx` is not verified simplex provenance. The
[official Dorado duplex metadata](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/#duplex-sequence-metadata)
describes the parent/offspring distinction. An extract-full TSV does not carry
this tag and requires upstream provenance establishing the supported input scope.

The adapter decodes deltas in original read orientation, interprets ML bytes as probability
intervals and applies CIGAR alignment to reference coordinates. Reverse-complemented SAM
records are mapped back to the original read orientation before MM decoding. The marker is
the actual modified reference base as a zero-based point, with strands kept separate:
the reference C on the plus strand and the reference G on the minus strand. It does not
collapse both strands into one CpG observation. The shared context hash is the hash of `CG`;
it is not the hash of a Nanopolish multi-motif sequence.

A local uncompressed FASTA is mandatory for SAM, BAM and modkit TSV. The parser checks its
SHA-256 against metadata and verifies CpG context against the actual reference. Reference
indexing reads the FASTA without creating or trusting an external `.fai`. Both sources must
use the same reference/build and adapter representation.

SAM/BAM sequence dictionaries (`@SQ` names and lengths) and mapped reference spans must
agree with the FASTA. Available Dorado `@PG` versions and structured `@RG` platform/model
fields are cross-checked. The summary records `header_verified_fields`; absent caller/model
fields remain declared provenance. An extract TSV cannot reconstruct the original header.

`configs/methylation/modbam.technical.json` locks the technical class-probability threshold
and resource limits. The default probability threshold of 0.8 is an engineering choice, not
a clinically validated cutoff. A confident 5mC call requires a 5mC probability lower bound
above the threshold; a confident canonical C requires sufficient residual probability after
**all** declared cytosine modifications. High 5hmC cannot be treated as canonical C.

Explicit `?`, implicit `.`, missing probabilities and omitted modification channels are
ambiguous, not confidently unmethylated. Malformed tags, stale MN lengths, impossible
probabilities, duplicate primary reads and duplicate read/site calls are rejected. Secondary,
supplementary, duplicate-flagged and QC-failed alignments are excluded and counted. Read names
remain local; the validation report contains counts and selection digests.

## modkit extraction

`parse_modkit_source` accepts the **modkit 0.6.1 extract full** per-read table, optionally
gzip compressed. Version 0.6.1 is the only accepted extraction release; both the metadata
model and the command builder enforce it. This is a source-reviewed compatibility pin,
not an executed-binary validation. The reviewed upstream commit is
[`481e3c9e7930f3f499eadf1ef441606f33e6881c`](https://github.com/nanoporetech/modkit/tree/481e3c9e7930f3f499eadf1ef441606f33e6881c).
The separate descriptive pileup policy's existing 0.4.1 pin does not apply to this adapter.
In particular, 0.4.1 extract-full tables omit `alignment_start` and `alignment_end` and are
rejected; the parser does not fabricate those coordinates. Other releases need their own
schema and executable conformance review before admission.

The parser requires a named header and the alignment/strand fields, original read positions,
modification codes and probabilities. Site-level m/h rows are combined before classification.
The modkit version, parent modBAM SHA-256 and the absence of probability transformation are
required provenance. BedMethyl/pileup files lack intact read identities and cannot be used.

The reviewed 0.6.1 writer emits 21 standard columns. Its optional `motifs` column is accepted
but is not required, including when `--cpg` is used; motif evidence is checked directly
against the locked FASTA. Named extra columns are tolerated with a consistent row width.
Headerless output (`--no-headers`) and call-only output (`extract calls`) are rejected.
`ref_position` and `forward_read_position` are zero-based; `alignment_end` is exclusive.
The alignment span must also fit the actual contig length in the locked FASTA;
an individually valid CpG cannot rescue a contradictory alignment extent.
For supported `C+` tags, `mod_strand` is `+`, and both reference-strand columns must agree
with SAM FLAG 0x10. Unmapped positions and excluded alignment flags remain exclusions.

Although the prose manual labels `mod_qual` as an integer, the reviewed implementation
outputs the floating midpoint `(ML + 0.5) / 256`. The parser recovers the byte bin with a
`1e-7` tolerance for printed float rounding, then uses the same probability interval as
the SAM/BAM adapter. Values outside that midpoint lattice fail, because they do not match
an untransformed 0.6.1 export. An `inferred=true` row must instead carry zero probability
and is retained as ambiguous. These rules preserve agreement at exact policy boundaries.

`build_modkit_extract_command` prepares a typed, shell-free command plan with a version-check
argument vector and a pinned expected version. It does not execute the command. Retain the
actual tool version, command, input/output fingerprints and success status in the controlled
upstream run record. Only the full table is appropriate; `--ignore`, probability
redistribution and `--pass-only` change its meaning and are not part of this adapter path.

A TSV cannot prove that the exporter retained every read/site. In particular, wholly omitted
reads cannot be reconstructed from `extract full`. The retained-read denominator is therefore
adapter-specific. The independent runner requires the same input format in all four roles;
it never mixes a SAM/BAM denominator with a TSV denominator. Filtered subsets require a
prospectively declared analysis domain and cannot establish a whole-specimen fraction.

## Integration and validation limits

The four-source runner accepts `input_format: modbam` for SAM/BAM and `input_format: modkit`
for per-read TSV in its local input bindings. Each entry supplies `metadata_path`,
`reference_path` and `adapter_policy_path`. The registered technical adapter name must match
the parsed `sam_mm_ml`, `modbam_mm_ml` or `modkit_extract_full_tsv` summary. The registered
parser version is `0.1.1` for SAM/TSV and `<pysam-version>/0.1.1` for BAM. The full adapter
policy has its own canonical JSON checksum. See
[`METHYLATION_VALIDATION.md`](METHYLATION_VALIDATION.md) for eligibility and execution.

Tests use synthetic reference sequences, SAM tags and modkit rows. Actual binary BAM
integration requires `pysam`; an unavailable runtime is explicitly reported rather than
silently substituting another parser. Public or institutional biological inputs remain
necessary to test interoperability, marker recovery and interval coverage. The parser is
bounded and in-memory for read calls; large WGS datasets require a separately reviewed,
versioned acquisition/subsetting plan within memory and storage limits.

The 2026-09-05 local conformance preflight found Python in WSL Ubuntu, but no `modkit` or
`samtools` in the Windows/WSL executable search paths or checked standard WSL tool paths.
No packages were installed. Consequently, a real extract execution has not been performed.
The executable gate must use a local, version-verified modkit 0.6.1 binary and freshly
generated synthetic SAM/BAM/FASTA with plus/minus alignments, explicit m/h probabilities,
implicit calls, CIGAR indels/clipping and excluded flags. Run the generated command without
a shell, then compare its parsed read/site/state map with the direct SAM adapter for the
explicitly represented sites. Include reference position zero and the last aligned base.
Unknown sites omitted by modkit must be assessed separately rather than added as canonical
calls or silently included in a cross-adapter denominator comparison.

Primary implementation references: [SAM optional tags specification](https://samtools.github.io/hts-specs/SAMtags.pdf),
[modkit extract documentation](https://nanoporetech.github.io/modkit/intro_extract.html),
[modkit 0.6.1 row writer](https://github.com/nanoporetech/modkit/blob/481e3c9e7930f3f499eadf1ef441606f33e6881c/modkit-core/src/read_ids_to_base_mod_probs.rs),
[modkit 0.6.1 probability conversion](https://github.com/nanoporetech/modkit/blob/481e3c9e7930f3f499eadf1ef441606f33e6881c/modkit-core/src/mod_bam.rs),
and [pysam API](https://pysam.readthedocs.io/en/latest/api.html).
