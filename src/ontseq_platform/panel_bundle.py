"""Import and provenance helpers for the GRCh38 AML Adaptive Sampling panel.

The laboratory sources use 1-based inclusive coordinates even though one of them has a
``.bed`` suffix.  They are therefore never consumed directly by an analysis.  Import keeps
their bytes unchanged and creates a separate, explicitly 0-based half-open selection BED.

This module intentionally does not derive gene bodies from the 10 kb flanks.  Analysis ROIs
and transcript caches are compiled from the locked GENCODE/MANE SQLite cache by
``panel_compiler``.  The contradictory ``IGH`` label is retained in source provenance but is
called ``IGH_REVIEW_REQUIRED`` in every active derivative.
"""

from __future__ import annotations

import hashlib
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from .coordinates import normalize_contig, one_based_inclusive_to_interval

PANEL_BUNDLE_ID: Final = "AML_AS_111_GRCh38_v1"
SOURCE_BED_NAME: Final = "250611_fusion_panel_with_buffer.bed"
SOURCE_REGIONS_NAME: Final = "250611_fusion_panel_with_buffer.interval_list"
SOURCE_BED_SHA256: Final = "f454644f18d8728c03678f4c6e969da7067879367c894c274e8c44be9352ef7e"
SOURCE_REGIONS_SHA256: Final = "f9ebbfbaa555b05d42fdd9edfda93eb9556661a8c940dbbd64b938769c40b441"
SOURCE_INTERVAL_COUNT: Final = 111
NORMALIZED_INTERVAL_BASES: Final = 17_028_488
UNRESOLVED_SOURCE_LABEL: Final = "IGH"
UNRESOLVED_ACTIVE_LABEL: Final = "IGH_REVIEW_REQUIRED"

_REGION_PATTERN = re.compile(r"^(?P<chrom>[^:\s]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$")


class PanelBundleError(ValueError):
    """The supplied source cannot establish the locked panel provenance."""


@dataclass(frozen=True)
class SourceInterval:
    """One laboratory interval in its declared 1-based inclusive convention."""

    chromosome: str
    start: int
    end: int
    label: str

    @property
    def normalized_start(self) -> int:
        return one_based_inclusive_to_interval(self.chromosome, self.start, self.end).start

    @property
    def normalized_label(self) -> str:
        if self.label == UNRESOLVED_SOURCE_LABEL:
            return UNRESOLVED_ACTIVE_LABEL
        return self.label


@dataclass(frozen=True)
class PanelImportSummary:
    bundle_id: str
    source_bed_sha256: str
    source_regions_sha256: str
    normalized_bed_sha256: str
    interval_count: int
    interval_bases: int
    unresolved_targets: tuple[str, ...]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_source_bed(path: Path) -> tuple[SourceInterval, ...]:
    """Read a BED-looking source under the locked one-based-inclusive contract."""
    records: list[SourceInterval] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise PanelBundleError(f"panel source cannot be read: {path}") from exc
    for line_number, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 4:
            raise PanelBundleError(
                f"{path.name} line {line_number}: expected four tab-separated fields"
            )
        chromosome, raw_start, raw_end, label = fields
        try:
            start, end = int(raw_start), int(raw_end)
        except ValueError as exc:
            raise PanelBundleError(
                f"{path.name} line {line_number}: start/end must be integers"
            ) from exc
        if start < 1 or end < start or not label.strip():
            raise PanelBundleError(
                f"{path.name} line {line_number}: invalid 1-based inclusive interval"
            )
        records.append(SourceInterval(normalize_contig(chromosome), start, end, label))
    if not records:
        raise PanelBundleError(f"{path.name}: no panel intervals")
    return tuple(records)


def read_source_regions(path: Path) -> tuple[tuple[str, int, int], ...]:
    """Read the laboratory ``chr:start-end`` list (not a Picard IntervalList)."""
    records: list[tuple[str, int, int]] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise PanelBundleError(f"panel regions source cannot be read: {path}") from exc
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        match = _REGION_PATTERN.fullmatch(line.strip())
        if match is None:
            raise PanelBundleError(f"{path.name} line {line_number}: expected chromosome:start-end")
        start, end = int(match["start"]), int(match["end"])
        if start < 1 or end < start:
            raise PanelBundleError(
                f"{path.name} line {line_number}: invalid 1-based inclusive interval"
            )
        records.append((normalize_contig(match["chrom"]), start, end))
    if not records:
        raise PanelBundleError(f"{path.name}: no panel intervals")
    return tuple(records)


def validate_panel_sources(
    source_bed: Path,
    source_regions: Path,
    *,
    require_locked_sources: bool = True,
) -> tuple[SourceInterval, ...]:
    """Validate source identity and exact numeric equivalence in file order."""
    if require_locked_sources:
        observed_bed = sha256_file(source_bed)
        if observed_bed != SOURCE_BED_SHA256:
            raise PanelBundleError(
                f"{source_bed.name}: expected source SHA-256 {SOURCE_BED_SHA256}, "
                f"observed {observed_bed}"
            )
        observed_regions = sha256_file(source_regions)
        if observed_regions != SOURCE_REGIONS_SHA256:
            raise PanelBundleError(
                f"{source_regions.name}: expected source SHA-256 {SOURCE_REGIONS_SHA256}, "
                f"observed {observed_regions}"
            )
    bed_records = read_source_bed(source_bed)
    region_records = read_source_regions(source_regions)
    numeric_bed = tuple((record.chromosome, record.start, record.end) for record in bed_records)
    if numeric_bed != region_records:
        mismatch = next(
            (
                index
                for index, (bed_item, regions_item) in enumerate(
                    zip(numeric_bed, region_records, strict=False), start=1
                )
                if bed_item != regions_item
            ),
            min(len(numeric_bed), len(region_records)) + 1,
        )
        raise PanelBundleError(
            "panel BED and regions list do not contain identical numeric intervals "
            f"in the same order (first disagreement at record {mismatch})"
        )
    if require_locked_sources and len(bed_records) != SOURCE_INTERVAL_COUNT:
        raise PanelBundleError(
            f"locked panel must contain {SOURCE_INTERVAL_COUNT} intervals, "
            f"observed {len(bed_records)}"
        )
    return bed_records


def normalized_bed_text(records: tuple[SourceInterval, ...]) -> str:
    """Render active 0-based half-open intervals without altering source artifacts."""
    return "".join(
        f"{record.chromosome}\t{record.normalized_start}\t{record.end}\t{record.normalized_label}\n"
        for record in records
    )


def import_panel_sources(
    source_bed: Path,
    source_regions: Path,
    bundle_directory: Path,
    *,
    require_locked_sources: bool = True,
) -> PanelImportSummary:
    """Copy immutable sources and create the normalized selection derivative.

    The caller owns atomic activation of ``bundle_directory``.  This function writes into a
    staging directory and refuses to overwrite any existing source or derivative.
    """
    records = validate_panel_sources(
        source_bed, source_regions, require_locked_sources=require_locked_sources
    )
    source_directory = bundle_directory / "source"
    derived_directory = bundle_directory / "derived"
    outputs = (
        source_directory / SOURCE_BED_NAME,
        source_directory / SOURCE_REGIONS_NAME,
        derived_directory / "selection_panel.normalized.bed",
    )
    if any(path.exists() for path in outputs):
        raise PanelBundleError("panel import refuses to overwrite an existing bundle artifact")
    source_directory.mkdir(parents=True, exist_ok=True)
    derived_directory.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_bed, outputs[0])
    shutil.copyfile(source_regions, outputs[1])
    outputs[2].write_text(normalized_bed_text(records), encoding="utf-8", newline="\n")
    interval_bases = sum(record.end - record.normalized_start for record in records)
    if require_locked_sources and interval_bases != NORMALIZED_INTERVAL_BASES:
        raise PanelBundleError(
            f"normalized panel must span {NORMALIZED_INTERVAL_BASES} bp, observed {interval_bases}"
        )
    return PanelImportSummary(
        bundle_id=PANEL_BUNDLE_ID,
        source_bed_sha256=sha256_file(outputs[0]),
        source_regions_sha256=sha256_file(outputs[1]),
        normalized_bed_sha256=sha256_file(outputs[2]),
        interval_count=len(records),
        interval_bases=interval_bases,
        unresolved_targets=(UNRESOLVED_ACTIVE_LABEL,),
    )
