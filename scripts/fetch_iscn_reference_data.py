from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path


MIN_ROWS = 800
USER_AGENT = "ontseq-platform/iscn-reference-fetcher"


@dataclass(frozen=True, slots=True)
class ReferenceFile:
    genome_build: str
    common_name: str
    url: str
    path: str
    sha256: str
    row_count: int
    size_bytes: int


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def count_supported_rows(text: str) -> int:
    supported = {str(value) for value in range(1, 23)} | {"X", "Y"}
    count = 0
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            fields = line.split()
        if len(fields) < 5:
            raise ValueError("cytoband table contains a malformed row")
        chromosome = fields[0].removeprefix("chr")
        if chromosome in supported:
            count += 1
    return count


def fetch_gzip_text(url: str, *, timeout: int = 60) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
        compressed = response.read()
    return gzip.decompress(compressed)


def fetch_reference(
    *,
    genome_build: str,
    common_name: str,
    url: str,
    output_dir: Path,
) -> ReferenceFile:
    payload = fetch_gzip_text(url)
    text = payload.decode("utf-8")
    row_count = count_supported_rows(text)
    if row_count < MIN_ROWS:
        raise ValueError(
            f"{common_name}: downloaded table contains only {row_count} supported rows"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{common_name}_cytoBand.txt"
    path.write_bytes(payload)
    return ReferenceFile(
        genome_build=genome_build,
        common_name=common_name,
        url=url,
        path=str(path),
        sha256=sha256_bytes(payload),
        row_count=row_count,
        size_bytes=len(payload),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch hg19/hg38 cytobands for the ONTSeq ISCN knowledge base"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--manifest",
        type=Path,
        help="JSON provenance manifest; defaults to <output-dir>/manifest.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    specifications = [
        (
            "GRCh37",
            "hg19",
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
        ),
        (
            "GRCh38",
            "hg38",
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
        ),
    ]
    records = [
        fetch_reference(
            genome_build=genome_build,
            common_name=common_name,
            url=url,
            output_dir=args.output_dir,
        )
        for genome_build, common_name, url in specifications
    ]

    manifest = args.manifest or args.output_dir / "manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(
        json.dumps(
            {
                "schema_version": "0.1.0",
                "purpose": "ONTSeq coordinate-to-cytoband reference provenance",
                "files": [asdict(record) for record in records],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(manifest)


if __name__ == "__main__":
    main()
