# Security policy

## Supported versions

Only the latest tagged `0.x` release receives security updates during the research phase.

## Never report sensitive findings publicly

Do not open public issues containing sample identifiers, patient data, genomic data,
screenshots, local server paths, credentials, or institutional infrastructure details.
Use the institutionally approved private reporting channel instead.

## Repository data policy

The following are prohibited in Git history:

- POD5/FAST5, FASTQ, BAM/CRAM/SAM, patient VCF/BCF and methylation tracks;
- names, dates of birth, medical-record numbers or re-identification keys;
- credentials, API keys, private certificates and signed download URLs;
- production reference bundles or licensed annotation databases.

Synthetic fixtures must be visibly labelled `SYNTHETIC` and must not be derived from a
real patient through simple renaming.

Run `python scripts/check_repository_safety.py` before every commit.
