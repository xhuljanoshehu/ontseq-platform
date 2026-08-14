# Data security and deployment boundary

## Default deployment

The analysis engine and genomic data remain on approved institutional/on-premises compute.
GitHub stores code and synthetic tests only. Documentation systems store SOPs and validation
records, not patient results.

## Data classes

| Class | Examples | GitHub | External automation |
| --- | --- | --- | --- |
| Public | source code, public references, synthetic fixtures | Allowed | Allowed |
| Internal | architecture, non-sensitive validation plans | Private only | Approved services only |
| Sensitive genomic | POD5, BAM/CRAM, patient VCF, reports | Prohibited | Prohibited by default |
| Direct identifiers | names, DOB, MRN, re-identification keys | Prohibited | Prohibited |

## Controls

- pseudonymous sample identifiers in manifests;
- checksum and atomic-transfer gate before execution;
- least-privilege service accounts and read-only input mounting;
- encrypted transport and storage;
- immutable release bundle with provenance and checksums;
- access logging, retention policy, backup and tested restore;
- explicit cloud-upload approval field defaulting to `false`;
- secret scanning and dependency updates in CI.

Data de-identification is not achieved merely by renaming genomic files. Human genomic data
can remain intrinsically identifying and must be governed accordingly.
