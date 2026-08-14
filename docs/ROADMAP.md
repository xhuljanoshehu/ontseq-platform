# Roadmap

## Milestone 0 - repository foundation (current)

- [x] Typed sample manifest and result contracts
- [x] Synthetic end-to-end JSON/HTML/XLSX demo
- [x] Limited ISCN proposal renderer with explicit warnings
- [x] Unit tests, repository safety check and CI configuration
- [x] Assay profiles, architecture, security and validation plan

## Milestone 1 - aligned BAM MVP

- [ ] Import or reimplement approved Cramino, CNV, SV and annotation rules
- [ ] BAM/BAI integrity, sort order, read-group and genome-build checks
- [ ] Versioned cytoband and reference resource manager
- [ ] Normalize raw caller outputs into the event schema
- [ ] Evidence-tiering review interface
- [ ] Reproduce thesis benchmark cases on approved data

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
