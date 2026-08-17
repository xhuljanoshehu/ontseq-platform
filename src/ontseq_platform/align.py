"""Minimap2 alignment adapter.

Takes an unaligned BAM and produces a coordinate-sorted, indexed BAM against the locked
reference, which is the input the rest of the pipeline already knows how to gate.

Parameter choices and why
-------------------------

``-x map-ont``
    The Oxford Nanopore preset. Using a preset built for another chemistry silently
    changes sensitivity, so it is pinned rather than left to minimap2's default.

``--MD``
    Emits the MD tag. Downstream long-read callers use it to reason about mismatches
    without re-reading the reference.

``-Y``
    Soft-clips supplementary alignments instead of hard-clipping them. Structural-variant
    callers need the clipped sequence to place a breakpoint; hard clipping discards it and
    quietly costs sensitivity at exactly the events this assay cares about.

``-y``
    Carries FASTQ comment tags through into the BAM. Combined with ``samtools fastq -T``
    this is what preserves modified-base tags (``MM``/``ML``) across alignment. Losing
    them would silently foreclose the methylation lane.

Read groups
-----------

Aligning through FASTQ is lossy by construction: ``samtools fastq`` drops the header, and
minimap2 writes a fresh one. Read-group provenance therefore has to be carried explicitly,
in two halves. The per-read ``RG:Z`` tag rides the FASTQ comment (``-T RG`` plus ``-y``),
and the ``@RG`` header lines are re-attached with ``samtools reheader`` after sorting.

minimap2's own ``-R`` option is deliberately not used: it accepts a single read-group line
and stamps it onto every record, which would silently merge distinct read groups into one
and contradict the per-read tags. Re-attaching the source header preserves however many
read groups the input actually declared.

Execution follows the repository's adapter boundary: no shell, an explicit argument
vector, a fail-closed version gate, and provenance for every tool involved.

Unlike the intermediate SAM, the sorted BAM and its index are treated as pipeline
artifacts. The SAM, FASTQ and pre-reheader BAM are scratch and are removed once the final
BAM exists.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .execution import CommandRunner, StreamingCommandRunner, SubprocessRunner
from .models import GenomeBuild, StrictModel, ToolRecord

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")

#: Modified-base tags preserved from the unaligned BAM through to the alignment.
MODIFIED_BASE_TAGS: tuple[str, ...] = ("MM", "ML")

#: Per-read tag that ties a record back to its ``@RG`` header line.
READ_GROUP_TAG = "RG"


class AlignmentPolicy(StrictModel):
    """Version-locked alignment configuration."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    expected_minimap2_version: str = Field(default="2.28", pattern=r"^\d+\.\d+(?:\.\d+)?$")
    expected_samtools_version: str = Field(default="1.24", pattern=r"^\d+\.\d+(?:\.\d+)?$")
    preset: Literal["map-ont", "lr:hq"] = "map-ont"
    emit_md_tag: bool = True
    soft_clip_supplementary: bool = True
    preserve_modified_base_tags: bool = True
    #: Minimum mapping quality retained downstream; recorded, not applied here.
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def soft_clipping_is_required_for_sv_work(self) -> AlignmentPolicy:
        if not self.soft_clip_supplementary:
            raise ValueError(
                "hard-clipped supplementary alignments discard the clipped sequence that "
                "structural-variant callers need; set soft_clip_supplementary"
            )
        return self


class AlignmentReport(StrictModel):
    """Normalized record of one alignment run."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    reference_id: str
    policy: AlignmentPolicy
    tools: list[ToolRecord] = Field(min_length=1)
    #: Envelope-relative paths; absolute source locations are never recorded.
    aligned_bam_relative_path: str = Field(min_length=1)
    aligned_bam_index_relative_path: str = Field(min_length=1)
    aligned_bam_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    modified_base_tags_requested: list[str] = Field(default_factory=list)
    #: ``@RG`` header lines carried over from the unaligned BAM. Zero means the input
    #: declared none, not that the adapter dropped them.
    read_group_count: int = Field(default=0, ge=0)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True


@dataclass(frozen=True)
class AlignmentInputs:
    """Everything the alignment stage needs from outside the envelope."""

    unaligned_bam: Path
    reference_fasta: Path


def _version(text: str, *, tool: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    raise ValueError(f"could not parse a version from {tool} output")


def probe_versions(
    runner: CommandRunner, *, minimap2: str = "minimap2", samtools: str = "samtools"
) -> dict[str, str]:
    """Probe both executables, failing closed when either cannot be identified."""
    versions: dict[str, str] = {}
    for name, argv in (
        ("minimap2", [minimap2, "--version"]),
        ("samtools", [samtools, "--version"]),
    ):
        result = runner.run(argv, timeout_seconds=60)
        if result.returncode != 0:
            raise ValueError(f"{name} version probe returned exit code {result.returncode}")
        versions[name] = _version(f"{result.stdout}\n{result.stderr}", tool=name)
    return versions


def _check_locked(observed: dict[str, str], policy: AlignmentPolicy) -> None:
    expected = {
        "minimap2": policy.expected_minimap2_version,
        "samtools": policy.expected_samtools_version,
    }
    for name, wanted in expected.items():
        found = observed.get(name)
        if found != wanted:
            raise ValueError(f"{name} version {found!r} does not match the policy lock {wanted!r}")


def read_group_lines(header_text: str) -> list[str]:
    """Extract the ``@RG`` lines of a SAM header, in the order they appear."""
    return [line for line in header_text.splitlines() if line.startswith("@RG\t")]


def header_with_read_groups(header_text: str, read_groups: Sequence[str]) -> str:
    """Insert ``@RG`` lines into a SAM header at their spec-mandated position.

    The SAM specification fixes the order ``@HD``, ``@SQ``, ``@RG``, ``@PG``, ``@CO``, so
    the lines go after the last ``@SQ`` record. Any ``@RG`` already present is replaced
    rather than appended to, so re-running this function is idempotent.
    """
    lines = [line for line in header_text.splitlines() if not line.startswith("@RG\t")]
    insert_at = 0
    for index, line in enumerate(lines):
        if line.startswith(("@HD\t", "@SQ\t")):
            insert_at = index + 1
    merged = lines[:insert_at] + list(read_groups) + lines[insert_at:]
    return "\n".join(merged) + "\n"


def build_minimap2_argv(
    *,
    minimap2: str,
    policy: AlignmentPolicy,
    reference_fasta: Path,
    reads_fastq: Path,
    output_sam: Path,
    threads: int,
) -> list[str]:
    """Build the alignment command as an explicit argument vector."""
    argv = [minimap2, "-a", "-x", policy.preset, "-t", str(threads)]
    if policy.emit_md_tag:
        argv.append("--MD")
    if policy.soft_clip_supplementary:
        argv.append("-Y")
    if policy.preserve_modified_base_tags:
        argv.append("-y")
    argv.extend(["-o", str(output_sam), str(reference_fasta), str(reads_fastq)])
    return argv


def run_alignment(
    inputs: AlignmentInputs,
    policy: AlignmentPolicy,
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    reference_id: str,
    scratch_dir: Path,
    output_bam: Path,
    runner: StreamingCommandRunner | None = None,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    threads: int = 4,
) -> AlignmentReport:
    """Align an unaligned BAM and produce a sorted, indexed BAM.

    The caller supplies the final BAM location, which lets the pipeline place it inside
    the run envelope while scratch files stay in a separate directory that can be cleaned
    without touching pipeline artifacts.
    """
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if not inputs.unaligned_bam.is_file():
        raise ValueError("unaligned BAM is missing or unreadable")
    if not inputs.reference_fasta.is_file():
        raise ValueError("reference FASTA is missing or unreadable")

    index_path = Path(f"{output_bam}.bai")
    for path in (output_bam, index_path):
        if path.exists():
            raise ValueError(f"refusing to overwrite an existing alignment artifact: {path.name}")

    command_runner = runner or SubprocessRunner()
    versions = probe_versions(command_runner, minimap2=minimap2, samtools=samtools)
    _check_locked(versions, policy)

    scratch_dir.mkdir(parents=True, exist_ok=True)
    fastq_path = scratch_dir / f"{sample_id}.reads.fastq"
    sam_path = scratch_dir / f"{sample_id}.aligned.sam"
    sorted_path = scratch_dir / f"{sample_id}.sorted.bam"
    header_path = scratch_dir / f"{sample_id}.header.sam"

    source_header = command_runner.run([samtools, "view", "-H", str(inputs.unaligned_bam)])
    if source_header.returncode != 0:
        raise ValueError(
            f"could not read the unaligned BAM header: exit code {source_header.returncode}"
        )
    read_groups = read_group_lines(source_header.stdout)

    carried_tags = list(MODIFIED_BASE_TAGS) if policy.preserve_modified_base_tags else []
    if read_groups:
        carried_tags.append(READ_GROUP_TAG)
    fastq_argv = [samtools, "fastq", "-@", str(threads)]
    if carried_tags:
        fastq_argv.extend(["-T", ",".join(carried_tags)])
    fastq_argv.extend(["-0", str(fastq_path), str(inputs.unaligned_bam)])
    _run_checked(command_runner, fastq_argv, label="samtools fastq")

    _run_checked(
        command_runner,
        build_minimap2_argv(
            minimap2=minimap2,
            policy=policy,
            reference_fasta=inputs.reference_fasta,
            reads_fastq=fastq_path,
            output_sam=sam_path,
            threads=threads,
        ),
        label="minimap2 alignment",
        timeout_seconds=14400,
    )
    # Sorting straight to the final path is safe only when nothing follows it; with read
    # groups to re-attach, the final path must not exist until the last step succeeds.
    sort_target = sorted_path if read_groups else output_bam
    _run_checked(
        command_runner,
        [samtools, "sort", "-@", str(threads), "-o", str(sort_target), str(sam_path)],
        label="samtools coordinate sort",
        timeout_seconds=14400,
    )

    if read_groups:
        aligned_header = command_runner.run([samtools, "view", "-H", str(sort_target)])
        if aligned_header.returncode != 0:
            raise ValueError(
                f"could not read the aligned BAM header: exit code {aligned_header.returncode}"
            )
        header_path.write_text(
            header_with_read_groups(aligned_header.stdout, read_groups), encoding="utf-8"
        )
        reheader = command_runner.run_to_file(
            [samtools, "reheader", str(header_path), str(sort_target)],
            output_bam,
            timeout_seconds=14400,
        )
        if reheader.returncode != 0:
            raise ValueError(
                f"samtools reheader failed with exit code {reheader.returncode}: "
                f"{(reheader.stderr or '').strip().splitlines()[-1:] or ['no diagnostic output']}"
            )

    _run_checked(
        command_runner,
        [samtools, "index", "-@", str(threads), str(output_bam), str(index_path)],
        label="samtools index",
    )

    for scratch in (fastq_path, sam_path, sorted_path, header_path):
        scratch.unlink(missing_ok=True)

    from .reference import sha256_file

    warnings: list[str] = []
    if not policy.preserve_modified_base_tags:
        warnings.append(
            "Modified-base tags were not requested; a later methylation lane cannot "
            "recover them from this alignment."
        )
    if not read_groups:
        warnings.append(
            "The unaligned BAM declared no @RG header line, so the alignment carries no "
            "read-group provenance; the intake gate reports this as a WARN."
        )

    return AlignmentReport(
        sample_id=sample_id,
        genome_build=genome_build,
        reference_id=reference_id,
        policy=policy,
        tools=[
            ToolRecord(
                name="minimap2",
                version=versions["minimap2"],
                parameters={
                    "preset": policy.preset,
                    "threads": threads,
                    "md_tag": policy.emit_md_tag,
                    "soft_clip_supplementary": policy.soft_clip_supplementary,
                    "carry_fastq_tags": policy.preserve_modified_base_tags,
                },
            ),
            ToolRecord(
                name="samtools",
                version=versions["samtools"],
                parameters={
                    "threads": threads,
                    "fastq_tags": carried_tags,
                    "read_groups_reattached": len(read_groups),
                },
            ),
        ],
        aligned_bam_relative_path=output_bam.name,
        aligned_bam_index_relative_path=index_path.name,
        aligned_bam_sha256=sha256_file(output_bam),
        modified_base_tags_requested=(
            list(MODIFIED_BASE_TAGS) if policy.preserve_modified_base_tags else []
        ),
        read_group_count=len(read_groups),
        warnings=warnings,
        limitations=[
            "Alignment parameters are technical defaults locked for reproducibility, not "
            "validated assay settings.",
            "Read groups are copied verbatim from the unaligned BAM; this adapter does not "
            "invent sample metadata and does not merge or rename read groups.",
            "Modified-base tags are carried through unchanged, as the SAM specification "
            "requires. Their interpretation on reverse-strand alignments has not been "
            "validated against a downstream methylation caller in this repository.",
            "No alignment quality gate is applied here; the aligned-BAM intake gate is "
            "the fail-closed boundary.",
        ],
    )


def _run_checked(
    runner: CommandRunner,
    argv: list[str],
    *,
    label: str,
    timeout_seconds: int = 3600,
) -> None:
    result = runner.run(argv, timeout_seconds=timeout_seconds)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        tail = detail[-1] if detail else "no diagnostic output"
        raise ValueError(f"{label} failed with exit code {result.returncode}: {tail}")
