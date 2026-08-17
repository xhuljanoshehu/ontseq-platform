"""A synthetic unaligned-BAM fixture that exercises the alignment lane with real tools.

The aligned-BAM smoke test (``smoke.py``) proves intake, QC and SV calling against real
binaries, but it hands the pipeline a BAM that is already aligned. Everything upstream of
that — ``samtools fastq``, minimap2, the read-group re-attachment — was therefore never
executed by CI, which is why :data:`StageId.ALIGN` could not honestly be called verified.

This module closes that gap. It writes a small synthetic reference and an unaligned BAM
whose reads are carved out of that reference, so alignment has a correct answer to find:

* **The reference is deterministic, not random.** Bases come from a SHA-256 keystream, so
  the same fixture is byte-identical on every machine and every rerun. That matters for
  the resume test, which compares content hashes rather than timestamps.
* **The reads are real substrings of the reference.** Random sequence would map nowhere
  and a passing alignment step would prove only that minimap2 exits zero. Carved reads
  mean the test can assert on mapped-read counts.
* **Some reads carry a 200 bp deletion.** Not so the SV caller is asserted on — that is
  the aligned-BAM smoke test's job — but so the reads are not all perfect matches and the
  aligner has to place a real gap.
* **Some reads are reverse-complemented.** Reverse-strand alignment is where modified-base
  tag handling is easiest to get wrong, so the fixture forces that path to be taken.
* **Every read carries MM/ML and RG tags.** These are exactly the tags that alignment
  through FASTQ can silently drop; CI asserts that they survive.

Nothing here is biological. The reference is not GRCh38 and the reads are not sequencing
data; the fixture exists to prove that a tool chain executes and preserves information.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from .execution import CommandRunner, SubprocessRunner
from .io import write_json
from .models import (
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    SampleManifest,
)
from .reference import reference_lock_from_fai

SAMPLE_ID = "SYNTHETIC_ALIGN_001"
REFERENCE_ID = "SYNTHETIC_NOT_REAL_GRCH38_ALIGN_V1"
READ_GROUP_ID = "SYNTHETIC_ALIGN_RG"

#: Contig layout of the synthetic reference. chr1 carries every read; chr2 exists so the
#: sequence dictionary has more than one entry and the intake gate compares a real set.
CONTIG_LENGTHS: tuple[tuple[str, int], ...] = (("chr1", 200_000), ("chr2", 50_000))

#: Where the reads are carved from, and the deletion they span.
READ_START = 20_000
READ_ARM_BP = 5_000
DELETION_BP = 200
DELETION_READS = 12
REFERENCE_READS = 8
REVERSE_READS = 4

_BASES = "ACGT"


def deterministic_sequence(length: int, *, seed: str) -> str:
    """Generate reproducible pseudo-random DNA from a SHA-256 keystream.

    ``random`` would also be deterministic, but only for a given Python build; hashing an
    explicit counter is stable across interpreters, which is what a fixture that feeds a
    content-addressed resume check needs.
    """
    if length < 1:
        raise ValueError("length must be positive")
    out: list[str] = []
    counter = 0
    while len(out) < length:
        block = hashlib.sha256(f"{seed}:{counter}".encode()).digest()
        out.extend(_BASES[byte & 0b11] for byte in block)
        counter += 1
    return "".join(out[:length])


def synthetic_reference() -> dict[str, str]:
    """Return the synthetic contigs keyed by name."""
    return {
        name: deterministic_sequence(length, seed=f"ontseq-align-fixture:{name}")
        for name, length in CONTIG_LENGTHS
    }


def format_fasta(sequences: dict[str, str], *, line_length: int = 60) -> str:
    """Render contigs as FASTA with a fixed line length, as ``samtools faidx`` expects."""
    if line_length < 1:
        raise ValueError("line_length must be positive")
    lines: list[str] = []
    for name, sequence in sequences.items():
        lines.append(f">{name}")
        lines.extend(
            sequence[offset : offset + line_length]
            for offset in range(0, len(sequence), line_length)
        )
    return "\n".join(lines) + "\n"


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _modification_tags(sequence: str) -> tuple[str, str] | None:
    """Build an MM/ML pair addressing the first two cytosines of a read.

    ``C+m?`` marks 5mC on the forward strand with the implicit-probability convention
    ``?`` (unobserved positions are unknown rather than unmodified), which is the honest
    encoding for a fixture that is not claiming to know anything about the other bases.
    """
    positions = [index for index, base in enumerate(sequence) if base == "C"][:2]
    if len(positions) < 2:
        return None
    # MM deltas count skipped cytosines, not skipped bases.
    skipped = sum(1 for base in sequence[positions[0] + 1 : positions[1]] if base == "C")
    first = sum(1 for base in sequence[: positions[0]] if base == "C")
    return f"MM:Z:C+m?,{first},{skipped};", "ML:B:C,220,180"


def _unaligned_record(name: str, sequence: str) -> str:
    """One unmapped SAM record: flag 4, no reference, tags that must survive alignment."""
    fields = [
        name,
        "4",
        "*",
        "0",
        "0",
        "*",
        "*",
        "0",
        "0",
        sequence,
        "I" * len(sequence),
        f"RG:Z:{READ_GROUP_ID}",
    ]
    tags = _modification_tags(sequence)
    if tags is not None:
        fields.extend(tags)
    return "\t".join(fields)


def unaligned_sam_text(reference: dict[str, str] | None = None) -> str:
    """Return the unaligned SAM: no ``@SQ`` records, every read unmapped."""
    contigs = reference if reference is not None else synthetic_reference()
    chromosome = contigs["chr1"]
    lines = [
        "@HD\tVN:1.6\tSO:unknown",
        f"@RG\tID:{READ_GROUP_ID}\tSM:{SAMPLE_ID}\tPL:ONT",
        "@PG\tID:ontseq-align-fixture\tPN:ontseq-align-fixture\tVN:0.1",
    ]
    span = READ_ARM_BP * 2 + DELETION_BP

    for index in range(DELETION_READS):
        start = READ_START + (index % 3)
        left = chromosome[start : start + READ_ARM_BP]
        right = chromosome[start + READ_ARM_BP + DELETION_BP : start + span]
        lines.append(_unaligned_record(f"SYNTH_ALIGN_DEL_{index + 1:03d}", left + right))

    for index in range(REFERENCE_READS):
        start = READ_START + (index % 3)
        lines.append(
            _unaligned_record(f"SYNTH_ALIGN_REF_{index + 1:03d}", chromosome[start : start + span])
        )

    # Reverse-complemented reads map to the reverse strand, which is the path on which
    # modified-base tags are most easily mishandled.
    for index in range(REVERSE_READS):
        start = READ_START + (index % 3)
        lines.append(
            _unaligned_record(
                f"SYNTH_ALIGN_REV_{index + 1:03d}",
                reverse_complement(chromosome[start : start + span]),
            )
        )
    return "\n".join(lines) + "\n"


@dataclass(frozen=True)
class AlignmentFixture:
    """Paths written by :func:`build_alignment_fixture`."""

    reference_fasta: Path
    reference_fai: Path
    unaligned_bam: Path
    manifest: Path
    reference_lock: Path


def _run_checked(runner: CommandRunner, argv: list[str], *, label: str) -> None:
    result = runner.run(argv, timeout_seconds=600)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        tail = detail[-1] if detail else "no diagnostic output"
        raise ValueError(f"{label} failed with exit code {result.returncode}: {tail}")


def build_alignment_fixture(
    output_dir: Path,
    *,
    runner: CommandRunner | None = None,
    samtools: str = "samtools",
) -> AlignmentFixture:
    """Write the reference, the unaligned BAM, the manifest and the reference lock.

    Refuses to overwrite: a fixture silently rebuilt on top of a previous one would make
    a resume test meaningless, because the inputs it hashes would have changed underneath
    it without anyone noticing.
    """
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = output_dir / "synthetic.reference.fasta"
    fai_path = output_dir / "synthetic.reference.fasta.fai"
    sam_path = output_dir / "synthetic.unaligned.sam"
    bam_path = output_dir / "synthetic.unaligned.bam"
    manifest_path = output_dir / "synthetic.manifest.json"
    lock_path = output_dir / "synthetic.reference-lock.json"

    existing = [
        path.name
        for path in (fasta_path, fai_path, bam_path, manifest_path, lock_path)
        if path.exists()
    ]
    if existing:
        raise ValueError("Refusing to overwrite existing fixture files: " + ", ".join(existing))

    command_runner = runner or SubprocessRunner()
    reference = synthetic_reference()
    fasta_path.write_text(format_fasta(reference), encoding="utf-8")
    _run_checked(command_runner, [samtools, "faidx", str(fasta_path)], label="samtools faidx")

    sam_path.write_text(unaligned_sam_text(reference), encoding="utf-8")
    _run_checked(
        command_runner,
        [samtools, "view", "-b", "-o", str(bam_path), str(sam_path)],
        label="samtools BAM conversion",
    )
    sam_path.unlink()

    manifest = SampleManifest(
        sample_id=SAMPLE_ID,
        run_id="SYNTHETIC_ALIGN_RUN_001",
        input=InputSpec(kind=InputKind.UNALIGNED_BAM, path=str(bam_path)),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id=REFERENCE_ID,
        ),
        analysis=AnalysisSpec(
            profile="synthetic-alignment-lane-smoke",
            modules=[AnalysisModule.QC, AnalysisModule.SV, AnalysisModule.REPORT],
        ),
    )
    reference_lock = reference_lock_from_fai(
        fai_path, reference_id=REFERENCE_ID, genome_build=GenomeBuild.GRCH38
    )
    write_json(manifest, manifest_path)
    write_json(reference_lock, lock_path)

    return AlignmentFixture(
        reference_fasta=fasta_path,
        reference_fai=fai_path,
        unaligned_bam=bam_path,
        manifest=manifest_path,
        reference_lock=lock_path,
    )
