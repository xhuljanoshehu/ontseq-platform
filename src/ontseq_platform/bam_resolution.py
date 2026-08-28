"""Fail-closed BAM/profile resolution for the manifest-free ``ontseq analyze`` path."""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from .bam_intake import ParsedBamHeader, parse_sam_header
from .execution import CommandRunner, SubprocessRunner, ToolExecutionError
from .models import GenomeBuild, ReferenceLock
from .reference import CanonicalReferenceSummary, validate_canonical_reference

_SAMPLE_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass(frozen=True)
class ResolvedBamInput:
    bam_path: Path
    index_path: Path
    sample_id: str
    genome_build: GenomeBuild
    naming_style: str
    contigs: tuple[tuple[str, int], ...]


def locate_bam_index(bam_path: Path) -> Path:
    """Use the documented BAI precedence and never guess beyond it."""

    bam = bam_path.expanduser().resolve()
    preferred = Path(f"{bam}.bai")
    alternate = bam.with_suffix(".bai")
    for candidate in (preferred, alternate):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "BAM index not found; expected "
        f"{preferred.name!r} or {alternate.name!r} beside {bam.name!r}"
    )


def sample_id_from_bam(bam_path: Path) -> str:
    """Create a stable manifest-safe sample ID from a validated BAM filename."""

    if bam_path.suffix.casefold() != ".bam":
        raise ValueError(f"aligned input must use the .bam extension: {bam_path.name}")
    candidate = _SAMPLE_SAFE.sub("_", bam_path.stem).strip("._-")
    if not candidate or not candidate[0].isalnum():
        candidate = f"sample_{candidate}".strip("_")
    if len(candidate) < 3:
        candidate = f"sample_{candidate or 'bam'}"
    return candidate[:64]


def default_run_id(sample_id: str, *, now: datetime | None = None) -> str:
    """Return ``<sample>-<UTC timestamp>`` while respecting the manifest length limit."""

    timestamp = (now or datetime.now(UTC)).astimezone(UTC).strftime("%Y%m%dT%H%M%SZ")
    suffix = f"-{timestamp}"
    return f"{sample_id[: 64 - len(suffix)]}{suffix}"


def detect_canonical_build(
    contigs: tuple[tuple[str, int], ...],
) -> tuple[GenomeBuild, CanonicalReferenceSummary]:
    """Detect one canonical build; partial, mixed or ambiguous dictionaries fail."""

    matches: list[tuple[GenomeBuild, CanonicalReferenceSummary]] = []
    errors: list[str] = []
    for build in GenomeBuild:
        try:
            summary = validate_canonical_reference(contigs, build)
        except ValueError as exc:
            errors.append(f"{build.value}: {exc}")
        else:
            matches.append((build, summary))
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise ValueError("BAM sequence dictionary is ambiguous between supported genome builds")
    raise ValueError(
        "BAM sequence dictionary is neither a complete canonical GRCh38 nor GRCh37 dictionary; "
        "partial, mixed-style and region-only dictionaries are rejected. " + " | ".join(errors)
    )


def validate_full_dictionary(
    parsed: ParsedBamHeader,
    reference_lock: ReferenceLock,
) -> None:
    """Require exact ordered agreement with the pinned ReferenceLock."""

    observed = parsed.contigs
    expected = tuple((item.name, item.length) for item in reference_lock.contigs)
    if observed == expected:
        return
    observed_map = dict(observed)
    expected_map = dict(expected)
    missing = [name for name in expected_map if name not in observed_map]
    extras = [name for name in observed_map if name not in expected_map]
    mismatched = [
        name
        for name, length in expected_map.items()
        if name in observed_map and observed_map[name] != length
    ]
    order_mismatch = not missing and not extras and not mismatched and observed != expected
    raise ValueError(
        "BAM sequence dictionary does not exactly match the profile's pinned ReferenceLock: "
        f"{len(missing)} missing, {len(extras)} extra, {len(mismatched)} length mismatches, "
        f"order_mismatch={order_mismatch}"
    )


def resolve_bam_header(
    *,
    bam_path: Path,
    header_text: str,
    reference_lock: ReferenceLock,
    required_build: GenomeBuild = GenomeBuild.GRCH38,
    sample_id: str | None = None,
) -> ResolvedBamInput:
    """Resolve a BAM from an already-read SAM header (pure and fixture-friendly)."""

    bam = bam_path.expanduser().resolve()
    if not bam.is_file():
        raise FileNotFoundError(f"BAM not found: {bam}")
    index = locate_bam_index(bam)
    parsed = parse_sam_header(header_text)
    detected_build, summary = detect_canonical_build(parsed.contigs)
    if detected_build != required_build:
        raise ValueError(
            f"BAM build is {detected_build.value}, but profile requires {required_build.value}; "
            "cross-build fallback is prohibited"
        )
    if reference_lock.genome_build != required_build:
        raise ValueError(
            f"ReferenceLock is {reference_lock.genome_build.value}, but profile requires "
            f"{required_build.value}"
        )
    validate_full_dictionary(parsed, reference_lock)
    resolved_sample = sample_id or sample_id_from_bam(bam)
    return ResolvedBamInput(
        bam_path=bam,
        index_path=index,
        sample_id=resolved_sample,
        genome_build=detected_build,
        naming_style=summary.naming_style,
        contigs=parsed.contigs,
    )


def resolve_bam_input(
    bam_path: Path,
    reference_lock: ReferenceLock,
    *,
    required_build: GenomeBuild = GenomeBuild.GRCH38,
    sample_id: str | None = None,
    samtools: str = "samtools",
    runner: CommandRunner | None = None,
) -> ResolvedBamInput:
    """Read the header with samtools and apply every pre-pipeline build gate."""

    bam = bam_path.expanduser().resolve()
    if not bam.is_file():
        raise FileNotFoundError(f"BAM not found: {bam}")
    locate_bam_index(bam)
    command_runner = runner or SubprocessRunner()
    try:
        result = command_runner.run([samtools, "view", "-H", str(bam)], timeout_seconds=120)
    except ToolExecutionError as exc:
        raise RuntimeError(f"samtools could not read the BAM header: {exc}") from exc
    if result.returncode != 0:
        detail = result.stderr.strip() or "no diagnostic output"
        raise RuntimeError(f"samtools view -H failed ({result.returncode}): {detail}")
    return resolve_bam_header(
        bam_path=bam,
        header_text=result.stdout,
        reference_lock=reference_lock,
        required_build=required_build,
        sample_id=sample_id,
    )


__all__ = [
    "ResolvedBamInput",
    "default_run_id",
    "detect_canonical_build",
    "locate_bam_index",
    "resolve_bam_header",
    "resolve_bam_input",
    "sample_id_from_bam",
    "validate_full_dictionary",
]
