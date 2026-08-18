# Synthetic knowledge-base fixture

`synthetic.clinvar.variant_summary.txt` is an **invented** file in NCBI's ClinVar
`variant_summary.txt` layout. It contains no real ClinVar record, no real gene and no real
condition — the gene is literally `SYNTHETIC_GENE` and the condition `Synthetic condition`.

It exists so continuous integration can exercise the annotation path end to end without
downloading a 200 MB weekly release, and it deliberately contains one row of each awkward
kind the loader has to handle correctly:

| VariationID | What it is there to prove |
| --- | --- |
| 900001 | A germline region record that matches, so a somatic run gets a marked scope mismatch |
| 900002 | A somatic record with a weak review status, so the caveat about single-submitter evidence fires |
| 900003 | A row for the **other assembly**, which must be dropped and counted, never matched |
| 900004 | A single nucleotide variant, which has no matchable extent and must be counted and skipped |
| 900005 | The `-1` placeholder ClinVar writes for unplaced records, which must not become an interval |

**This is not ClinVar.** A real run needs the genuine release from NCBI, locked by checksum
with `--release` naming the publication date. See ADR-022.
