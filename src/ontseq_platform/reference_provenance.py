from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import Field

from .models import GenomeBuild, StrictModel
from .reference import contig_signature, reference_lock_from_fai, sha256_file


class SequenceReferenceLock(StrictModel):
    """Path-free sequence-level reference provenance for research-only analysis."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    reference_id: str = Field(min_length=1)
    genome_build: GenomeBuild
    source_fasta_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    source_fai_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    contig_signature_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


def build_sequence_reference_lock(
    fasta_path: Path,
    fai_path: Path,
    *,
    reference_id: str,
    genome_build: GenomeBuild,
) -> SequenceReferenceLock:
    """Fingerprint a FASTA and its FAI without retaining paths or sequence content."""

    if not fasta_path.is_file():
        raise ValueError("Reference FASTA is missing or unreadable")
    if not fai_path.is_file():
        raise ValueError("Reference FASTA index is missing or unreadable")

    fai_lock = reference_lock_from_fai(
        fai_path,
        reference_id=reference_id,
        genome_build=genome_build,
    )
    signature = contig_signature((contig.name, contig.length) for contig in fai_lock.contigs)
    return SequenceReferenceLock(
        reference_id=reference_id,
        genome_build=genome_build,
        source_fasta_sha256=sha256_file(fasta_path),
        source_fai_sha256=fai_lock.source_fai_sha256,
        contig_signature_sha256=signature,
    )


def verify_sequence_reference_lock(
    lock: SequenceReferenceLock,
    *,
    fasta_path: Path,
    fai_path: Path,
    reference_id: str,
    genome_build: GenomeBuild,
) -> None:
    """Fail closed when runtime reference bytes or declared identity differ from the lock."""

    if lock.reference_id != reference_id:
        raise ValueError("Reference ID does not match the sequence reference lock")
    if lock.genome_build != genome_build:
        raise ValueError("Genome build does not match the sequence reference lock")

    observed = build_sequence_reference_lock(
        fasta_path,
        fai_path,
        reference_id=reference_id,
        genome_build=genome_build,
    )
    if observed.source_fasta_sha256 != lock.source_fasta_sha256:
        raise ValueError("Reference FASTA SHA256 does not match the sequence reference lock")
    if observed.source_fai_sha256 != lock.source_fai_sha256:
        raise ValueError("Reference FAI SHA256 does not match the sequence reference lock")
    if observed.contig_signature_sha256 != lock.contig_signature_sha256:
        raise ValueError("Reference contig signature does not match the sequence reference lock")
