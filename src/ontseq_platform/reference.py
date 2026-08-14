from __future__ import annotations

import hashlib
from collections.abc import Iterable
from pathlib import Path

from .models import GenomeBuild, ReferenceContig, ReferenceLock


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def contig_signature(contigs: Iterable[tuple[str, int]]) -> str:
    digest = hashlib.sha256()
    for name, length in contigs:
        digest.update(f"{name}\t{length}\n".encode())
    return digest.hexdigest()


def reference_lock_from_fai(
    fai_path: Path,
    *,
    reference_id: str,
    genome_build: GenomeBuild,
    allow_extra_contigs: bool = False,
) -> ReferenceLock:
    contigs: list[ReferenceContig] = []
    with fai_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Invalid FASTA index line {line_number}: expected at least 2 fields"
                )
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid FASTA index length on line {line_number}: {fields[1]!r}"
                ) from exc
            contigs.append(ReferenceContig(name=fields[0], length=length))
    if not contigs:
        raise ValueError(f"FASTA index contains no contigs: {fai_path}")
    return ReferenceLock(
        reference_id=reference_id,
        genome_build=genome_build,
        contigs=contigs,
        allow_extra_contigs=allow_extra_contigs,
        source_fai_sha256=sha256_file(fai_path),
    )
