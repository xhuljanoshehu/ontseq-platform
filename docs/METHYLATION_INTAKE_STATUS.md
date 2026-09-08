# Approved methylation intake and Nanopolish hardening — 2026-09-05

**Acquisition completed; biological readiness remains NO_CALL / NOT_REGISTERED.**
The user explicitly authorized downloading the two reviewed public NCBI files after
the safety review. Only those two objects were acquired, using HTTPS with certificate
verification, no redirects, exact byte limits and exclusive writes. Existing files,
study matrices and historical outputs were preserved.

## Provenance and local artifacts

`configs/methylation/acquired_public_sources.v1.json` catalogues the exact public URLs,
retrieval times, byte sizes and SHA-256 fingerprints. It contains no reads, measurement
tables or clinical reports. Supplier/reference checksums that remain unknown were not
invented. This acquisition supersedes the earlier metadata-only acquisition status,
without rewriting the historical candidate catalogue.

All genomic files and detailed intake results remain in the ignored directory
`work/methylation_validation/gse218679/intake_v1/`:

- `acquisition.json`: immutable two-object download provenance and actual hashes.
- `structural-intake.json`: bounded gzip/TSV inspection, counts and coordinate structure.
- `fixed-pool-inventory.json`: hardened parser and prospectively specified split results.
- `intake.readiness.methylation-holdout.json`: typed NO_CALL readiness and reasons.

Gzip integrity/EOF, complete TSV headers, finite numerical fields, expected coordinate,
strand and sequence encodings, identifier controls and cross-source read overlap were
checked. Raw read identifiers and sequences were not printed or placed in Git. The
bounded structural check did not calculate differential methylation, select markers
or run recovery. The hardened parser subsequently confirmed the retained read pools.

The unchanged matrix's largest read budget cannot be supplied by the actual held-out
pools. Exact pore chemistry, Guppy configuration and the upstream reference artifact
fingerprint also remain unresolved. Consequently no full cohort was invented or
registered and no recovery metrics are claimed. A subsequent feasibility design, if
appropriate, must have a new version and be fixed before biological outcome inspection;
the current matrix and its unsuccessful admission remain recorded.

The matrix fingerprint remains
`fe54b78dcef5f6fd12a70fb8f7735a25591412cc6add2d20a9999d7cba19731e`.
The current source/runtime fingerprint is
`97f72d0388028139906853dc64166fdc569a6ee7fdfb1ae969bedefa2a972a85`.

## Security change and verification

`nanopolish-call-table-v2` closes the two findings from
`METHYLATION_DATA_SAFETY_REVIEW.md`. Decoded byte, physical-line and header/column
limits apply before CSV allocation. Control characters are rejected in read identifiers
and read-selection digests, preserving hashes for valid existing identifiers. Source
summaries record the decoded byte count and `bounded-tsv-gzip-v1` security profile.
Biological estimator equations and frozen matrix values are unchanged. Validation
impact is reviewed in `CLINICAL_VALIDATION.md` and the changelog.

The requested make commands were attempted; `make` is absent. Their underlying
commands ran directly in the existing Windows virtual environment:

| Check | Result |
| --- | --- |
| Repository safety, versions, Ruff lint/format, mypy, schema consistency | Passed |
| Methylation suites | 124 passed |
| MM/ML/modkit suite | 47 passed, 1 skipped because pysam is absent |
| CLI surface and repository safety regression suites | 21 passed |
| Full repository | 900 tests; 7 failures, 13 errors, 1 skip |

The full suite remains unsuccessful for the same 20 tests outside methylation as
before: Windows paths/line endings, symlink privileges, permissions, process locks,
service manifests and XLSX cleanup. No test was suppressed. The complete new log is
`results/methylation-holdout/20260905/post-intake-safety-full-suite.log`.

The structural checks and importer regressions are not an antivirus certificate.
No separate antivirus scan was initiated. Read-only Defender settings showed cloud
reporting and automatic safe-sample submission enabled; a local-only manual scan was
therefore not assumed. No Defender settings were changed. Ordinary OS-level protection
or telemetry was not audited. See Microsoft's
[cloud protection and sample submission description](https://learn.microsoft.com/en-us/defender-endpoint/cloud-protection-microsoft-antivirus-sample-submission)
for the distinction.

Research-use restrictions persist. No tumour purity, physical DNA-mass fraction,
cell fraction, subclone abundance or detection-limit claim follows from this intake.
