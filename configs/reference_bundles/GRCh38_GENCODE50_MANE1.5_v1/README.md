# GRCh38_GENCODE50_MANE1.5_v1 release lock

`bundle.recipe.yaml` is the reviewed, installable source lock. Every publisher byte stream has
an exact byte count and SHA256. The GENCODE FASTA and GTF hashes were additionally checked against
the MD5 values in the publisher's release-50 `MD5SUMS` file. No checksum was inferred from a URL,
filename or decompressed derivative.

The release-lock procedure was:

1. download every listed source into an approved staging directory;
2. record the exact byte count and SHA256 of each downloaded byte stream (`sha256sum` on Linux,
   `Get-FileHash -Algorithm SHA256` on PowerShell);
3. retain the publisher URL, release and source date beside those values;
4. normalize publisher-native repeat and duplication tables into the declared GRCh38 interval
   contract; the pinned Hoffman Umap k100 source is already a native hg38 BED.gz;
5. construct `bundle.recipe.yaml` using the `ReferenceBundle` 1.0.0 schema;
6. retain the generated source and derived hashes in the activated `bundle.yaml` and exercise the
   opt-in full installation smoke before changing the profile pin.

The production installer refuses a source without a lowercase 64-character SHA256, rejects a
size/hash mismatch before activation, and writes `references/<bundle-id>/bundle.yaml` only after
the derived FASTA index, reference lock and annotation cache pass validation. Once activated,
normal analyses do not perform internet requests.

For an offline transfer, install on a connected approved machine, copy the entire activated
bundle directory, then use `ReferenceBundleInstaller.import_bundle(PATH)`. Import recalculates
every checksum before it atomically activates the local copy.

The deterministic miniature installer contract is in
`tests/fixtures/reference_catalog/GRCh38_FIXTURE_v1/bundle.recipe.yaml`. It carries real SHA256 and
size values for every source fixture, and its tiny FASTA has a byte-exact matching FAI. Transaction
tests mock only the canonical-length gate; a separate test proves that the fixture cannot be
activated through the production path. This exercises installer/repair/cache behavior without
shipping a genome or adding a runtime validation bypass.
