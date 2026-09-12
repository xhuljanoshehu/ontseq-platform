from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections.abc import Iterator
from pathlib import Path

_CANONICAL = re.compile(r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _lines(path: Path) -> Iterator[str]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        yield from handle


def _canonical_chromosome(value: str) -> str:
    return value if value.startswith("chr") else f"chr{value}"


def _gene_name(attributes: str) -> str | None:
    parsed: dict[str, str] = {}
    for item in attributes.rstrip(";").split(";"):
        key, _, raw_value = item.strip().partition(" ")
        parsed[key] = raw_value.strip().strip('"')
    return parsed.get("gene_name") or parsed.get("gene_id")


def normalize_genes(source: Path) -> list[tuple[str, int, int, str]]:
    records: list[tuple[str, int, int, str]] = []
    for line_number, raw_line in enumerate(_lines(source), start=1):
        if not raw_line or raw_line.startswith("#"):
            continue
        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"GTF line {line_number}: expected nine columns")
        if fields[2] != "gene" or _CANONICAL.fullmatch(fields[0]) is None:
            continue
        label = _gene_name(fields[8])
        if not label:
            raise ValueError(f"GTF line {line_number}: gene record has no gene name or ID")
        start, end = int(fields[3]) - 1, int(fields[4])
        if start < 0 or end <= start:
            raise ValueError(f"GTF line {line_number}: invalid coordinates")
        records.append((_canonical_chromosome(fields[0]), start, end, label))
    return sorted(set(records), key=lambda item: (item[0], item[1], item[2], item[3]))


def normalize_four_column(source: Path) -> list[tuple[str, int, int, str]]:
    records: list[tuple[str, int, int, str]] = []
    for line_number, raw_line in enumerate(_lines(source), start=1):
        line = raw_line.rstrip("\r\n")
        if not line or line.startswith(("#", "track ", "browser ")):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"interval line {line_number}: expected at least four columns")
        if _CANONICAL.fullmatch(fields[0]) is None:
            continue
        start, end = int(fields[1]), int(fields[2])
        if start < 0 or end <= start or not fields[3].strip():
            raise ValueError(f"interval line {line_number}: invalid coordinates or label")
        records.append((_canonical_chromosome(fields[0]), start, end, fields[3].strip()))
    return sorted(set(records), key=lambda item: (item[0], item[1], item[2], item[3]))


def _write_resource(
    records: list[tuple[str, int, int, str]],
    *,
    output: Path,
    source: Path,
    resource_id: str,
    resource_type: str,
    source_name: str,
    release: str,
    genome_build: str,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "".join(f"{chrom}\t{start}\t{end}\t{label}\n" for chrom, start, end, label in records),
        encoding="utf-8",
        newline="\n",
    )
    lock = {
        "schema_version": "0.1.0",
        "resource_id": resource_id,
        "resource_type": resource_type,
        "source": source_name,
        "release": release,
        "genome_build": genome_build,
        "sha256": _sha256(output),
        "source_sha256": _sha256(source),
        "coordinate_system": "zero_based_half_open",
        "columns": "chrom_start_end_label",
        "note": "Generated locally; source and normalized resources are checksum locked.",
    }
    output.with_suffix(output.suffix + ".lock.json").write_text(
        json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize build-specific SV annotation resources and emit checksum locks."
    )
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--resource-type",
        required=True,
        choices=[
            "genes",
            "cytobands",
            "repeatmasker",
            "tandem_repeat",
            "segmental_duplication",
            "blacklist",
            "mappability",
            "centromere",
            "telomere",
        ],
    )
    parser.add_argument("--resource-id", required=True)
    parser.add_argument("--source-name", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--genome-build", required=True, choices=["GRCh37", "GRCh38"])
    args = parser.parse_args()
    if not args.source.is_file():
        raise SystemExit(f"source resource not found: {args.source}")
    records = (
        normalize_genes(args.source)
        if args.resource_type == "genes"
        else normalize_four_column(args.source)
    )
    if not records:
        raise SystemExit("no canonical records survived normalization")
    _write_resource(
        records,
        output=args.output,
        source=args.source,
        resource_id=args.resource_id,
        resource_type=args.resource_type,
        source_name=args.source_name,
        release=args.release,
        genome_build=args.genome_build,
    )


if __name__ == "__main__":
    main()
