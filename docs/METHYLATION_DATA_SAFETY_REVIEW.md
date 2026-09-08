# Safety review of the proposed GSE218679 files

**Follow-up:** The user subsequently authorized acquisition. The exact two files are
now stored locally with verified hashes and bounded structural checks. The two importer
findings below have been fixed in `nanopolish-call-table-v2` with synthetic regressions.
See `METHYLATION_INTAKE_STATUS.md` for current status. The assessment below records the
earlier, metadata-only review and must not be mistaken for the current acquisition state.

Assessed 2026-09-05. This is a metadata, transport and local-parser review, not a
malware scan, clinical report, biological validation or authorization to acquire data.
No genomic file body or byte range was fetched. No raw donor records were inspected.

## Conclusion

The exact objects have verified links from official public GEO records, and live
HTTPS HEAD requests succeed against the official NCBI host. This supports source
and transport credibility. It does not establish the content, checksum, malware
status, anonymity or scientific validity of either compressed file.

Direct import of uninspected gzip data through the current Nanopolish parser is
not sufficiently bounded. Two local security/integrity findings below were reproduced
using only tiny synthetic inputs. They are not evidence that these public files are
malicious. Hardening and bounded structural intake remain necessary before import.

## Exact objects and transport evidence

| GEO record | Exact basename | HEAD Content-Length |
| --- | --- | ---: |
| [GSM6754749](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6754749) | GSM6754749_135_Nanopolish.tsv.gz | 87,638 |
| [GSM6754755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6754755) | GSM6754755_barcode01_Nanopolish.tsv.gz | 120,883 |

Both records list the exact processed TSV supplements and public release on
2023-12-21. Exact HTTPS URLs are preserved in
`configs/methylation/additional_data_candidates.v1.json`.

Live HEAD checks at 2026-09-05 13:29:06 UTC used curl with HTTPS-only protocol,
normal certificate verification, no redirect following and no file body:

- HTTP 200 at the original `ftp.ncbi.nlm.nih.gov` URL for each object.
- TLS verification result 0 after successful connections; no insecure flag used.
- `Content-Type: application/x-gzip`; combined compressed size 208,521 bytes.
- `Last-Modified: Wed, 23 Nov 2022 19:59:00 GMT` on both objects.
- HSTS advertised; no cross-host redirect observed.

The initial sandboxed connections failed; subsequent authorized header-only network
checks succeeded. A zero TLS result from a failed connection was not treated as proof.
Names, MIME headers and size agreement do not substitute for byte inspection.
No authenticated publisher SHA-256 was found in the reviewed metadata/headers;
local SHA-256 remains uncomputed. A future local hash records acquired bytes and
does not independently prove authenticity or absence of malware.

## Access, privacy and scientific scope

[GEO's disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) permits access
and does not impose NCBI use/distribution restrictions, while explicitly declining
to guarantee third-party rights or independently verify biological quality. Public
availability must not be described as CC0 or an unrestricted legal clearance.

The [study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10794044/) and the two GSM records
describe human iPSC genomic DNA, including a patient-derived FXS line and a control.
Public cell-line identifiers do not establish anonymization. No source-specific
consent document, IRB number or anonymization statement was located in the reviewed
GSM records or full PMC main text. This is an evidence gap, not evidence of missing
consent. [GEO's human-subject rules](https://www.ncbi.nlm.nih.gov/geo/info/faq.html)
place appropriate consent/permission responsibilities on submitters.

No correction or retraction link was found in the reviewed
[PubMed record](https://pubmed.ncbi.nlm.nih.gov/38134876/) or
[publisher-deposited Crossref metadata](https://api.crossref.org/works/10.1016/j.cell.2023.11.019).
This limited status check cannot certify the study's results.

The data concern targeted FMR1 measurements. The described analysis covers 19 nearby
promoter CpGs; exact export scope, protocol gaps and eligibility remain unresolved.
They cannot independently validate genome-wide performance, tumour purity or a
detection limit. See `METHYLATION_ADDITIONAL_DATA_SEARCH.md` for the admission gaps.

## Local importer findings

Source reviewed: `src/ontseq_platform/methylation_mixture.py`.

1. **Missing decompressed-size and physical-line bounds.** `_open_nanopolish` uses
   gzip as a text stream. `parse_nanopolish_source` checks compressed disk size,
   then gives complete text records to `csv.DictReader`; its row limit counts only
   completed nonempty records. It does not cap total decoded bytes, physical line
   length or total memory/time. A bounded synthetic demonstration used 681 compressed
   bytes expanding to 10,532 bytes with a configured 1,024-byte input limit. The
   actual parser accepted 200 rows and 200 read groups. This demonstrates the limit's
   scope; no large decompression or resource-exhaustion test was run.
2. **Ambiguous read-selection serialization.** Read identifiers permit internal
   control characters, while `_selection_digest` uses newline/NUL delimiters without
   escaping or length prefixes. Two different two-name synthetic pools produced the
   same selection digest under the same context. This is a serialization ambiguity,
   not a cryptographic SHA-256 collision. Complete source-file hashes remain distinct.

No TSV-driven shell execution, eval, pickle deserialization, archive extraction or
embedded-path writing was found in this import path. Raw read names are not exported;
HTML strings are escaped and exported measurements are typed. These properties
reduce specific attack paths but do not resolve the findings above or certify every
dependency/runtime version as vulnerability-free.

No parser code was modified during this safety assessment. Before using untrusted
inputs, add bounded decoded streaming and line/header limits, reject unsupported
identifier control characters or use unambiguous serialization, and test those guards.
Then inspect approved files as data with size/time/memory bounds and immutable
fingerprints before running biological recovery. No upload to a third-party scanner
has occurred or been proposed.

## Persistent working rule

At the user's request, AGENTS.md rule 8 now requires source, access/use and technical
risk checks before future external data acquisition, with explicit distinctions
between metadata and content verification. Existing authorization rules are preserved.
