# Roadmap

## Milestone 0 - repository foundation

- [x] Typed sample manifest and result contracts
- [x] Synthetic end-to-end JSON/HTML/XLSX demo
- [x] Limited ISCN proposal renderer with explicit warnings
- [x] Unit tests, repository safety check and CI configuration
- [x] Assay profiles, architecture, security and validation plan

## Milestone 1 - aligned BAM MVP

- [ ] Freeze intended use, assay modes, truth definitions and acceptance criteria
- [x] Implement a normalized Cramino adapter from its official JSON interface
- [ ] Implement benchmarked CNV, SV and annotation adapters from official interfaces
- [x] Implement candidate-only Sniffles2 v2.8.0 adapter from its official CLI/VCF interfaces
- [x] BAM/BAI integrity, sort order, read-group and genome-build checks
- [x] Versioned reference lock generated from the exact FASTA index
- [ ] Versioned cytoband resource manager
- [ ] Normalize raw caller outputs into the event schema
- [x] Normalize Sniffles2 DEL/DUP/INV/INS/BND evidence into the event schema
- [x] Record per-module `COMPLETED`, `NOT_RUN`, `FAILED` and `NO_CALL` outcomes
- [ ] Interactive evidence-tiering review interface
- [x] Deterministic synthetic CNV/SV benchmark contract and CI fixtures
- [ ] Benchmark technical SV behavior on GIAB HG002 and draft HG008 tumor/normal resources
- [x] Represent CNV truth per source, with background semantics, resolution and breakpoint uncertainty
- [x] Score CNV per base over an explicit observability mask, independent of segmentation
- [x] Separate biological negativity, `NO_CALL` and technical failure in CNV accounting
- [x] Deterministic CNV dilution/coverage simulation and limit-of-detection estimation
- [x] Paired method comparison so CNV method selection rests on a paired statistic
- [x] Non-reportable baseline read-depth caller as an experimental control
- [x] ISCN karyotype to CNV truth conversion against a versioned cytoband resource
- [ ] Lock a real cytoband resource per genome build with checksum and provenance
- [ ] Add purity/ploidy estimation so detection is not limited by a fixed rounding band
- [ ] Add version-pinned Spectre, ichorCNA and QDNAseq/ACE execution adapters
- [ ] Establish whether adaptive-sampling off-target reads support depth-based CNV locally
- [ ] Benchmark CNV candidates across coverage and tumor/blast-fraction dilution series
- [ ] Evaluate the locked pipeline on an orthogonally characterized AML cohort
- [ ] Record no-call behavior and failure modes separately from negative results
- [x] Exercise samtools, Cramino and Sniffles2 in CI with a generated positive BAM fixture

## Milestone 2 - POD5 end-to-end

- [ ] Dorado model manifest, basecalling and resume support
- [ ] Demultiplexing and Minimap2 alignment
- [ ] GPU/CPU/HPC execution profiles
- [ ] Watchfolder or LIMS trigger with atomic handoff
- [ ] Failure recovery and operational dashboard

## Milestone 3 - ISCN and analytical validation

- [ ] Authorized ISCN 2024 test corpus and expert-reviewed rule inventory
- [ ] Coordinate/cytoband edge-case suite for GRCh37 and GRCh38
- [ ] Locked positive/negative reference cohort
- [ ] LoD, precision, reproducibility and interference studies
- [ ] Signed release and change-control workflow

## Evidence review gates

- [ ] Review the evidence registry quarterly and before changing a scientific dependency
- [ ] Record inclusion rationale, applicability and limitations for every production adapter
- [ ] Pre-register benchmark datasets, metrics and acceptance criteria before comparing callers
- [ ] Require validation-impact review for caller, model, reference or target-BED updates

## Milestone 4 - application layer

- [ ] Role-based API and authentication
- [ ] Single-sample upload/intake interface
- [ ] Interactive QC, evidence and ISCN review
- [ ] Electronic signature and immutable release bundle
- [ ] Institution-approved LIMS integration

## Optional research lanes

- [ ] Small variants/indels
- [ ] Modified-base and methylation analysis
- [ ] RNA fusion confirmation
- [ ] cfDNA tissue-of-origin research

Each lane gets a separate intended use, schema extension, QC contract and validation plan.
