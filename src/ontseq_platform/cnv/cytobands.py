"""Versioned cytoband resource handling.

Cytogenetic truth is expressed in band coordinates, and band coordinates are
build-specific. Converting ``5q13`` into base pairs therefore requires an explicit,
versioned, checksummed resource, exactly like the reference lock does for contigs.

The repository does not ship a cytoband table. Band definitions are build-specific
annotation that belongs with the reference bundle on approved storage, and the data
boundary in ``docs/DATA_SECURITY.md`` keeps reference bundles out of Git. A small
synthetic table is used in tests instead, and production use loads a locked file whose
checksum is recorded in provenance.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from ..models import GenomeBuild, StrictModel
from .intervals import Interval, canonical_contig

BAND_PATTERN = re.compile(r"^(?P<arm>[pq])(?P<major>\d+)(?:\.(?P<minor>\d+))?$")


class Cytoband(StrictModel):
    """One banding interval."""

    contig: str = Field(pattern=r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    name: str = Field(min_length=1)
    stain: str = Field(min_length=1)

    @model_validator(mode="after")
    def end_after_start(self) -> Cytoband:
        if self.end <= self.start:
            raise ValueError("cytoband end must be greater than start")
        return self

    @property
    def arm(self) -> str:
        return self.name[0]


class CytobandTable(StrictModel):
    """A checksummed, build-specific set of cytobands."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    resource_id: str = Field(min_length=1)
    genome_build: GenomeBuild
    source: str = Field(min_length=1)
    source_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    bands: list[Cytoband] = Field(min_length=1)
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def bands_do_not_overlap(self) -> CytobandTable:
        by_contig: dict[str, list[Cytoband]] = {}
        for band in self.bands:
            by_contig.setdefault(canonical_contig(band.contig), []).append(band)
        for contig, items in by_contig.items():
            items.sort(key=lambda item: item.start)
            for previous, current in zip(items, items[1:], strict=False):
                if current.start < previous.end:
                    raise ValueError(
                        f"cytobands overlap on contig {contig}: "
                        f"{previous.name} and {current.name}"
                    )
        return self

    def contig_bands(self, contig: str) -> list[Cytoband]:
        target = canonical_contig(contig)
        return sorted(
            (band for band in self.bands if canonical_contig(band.contig) == target),
            key=lambda band: band.start,
        )

    def contig_length(self, contig: str) -> int:
        bands = self.contig_bands(contig)
        if not bands:
            raise KeyError(f"no cytobands for contig {contig}")
        return bands[-1].end

    def arm_interval(self, contig: str, arm: str) -> Interval:
        """Return the span of a whole chromosome arm."""
        if arm not in {"p", "q"}:
            raise ValueError("arm must be 'p' or 'q'")
        bands = [band for band in self.contig_bands(contig) if band.arm == arm]
        if not bands:
            raise KeyError(f"contig {contig} has no {arm} arm bands")
        return bands[0].start, bands[-1].end

    def band_interval(self, contig: str, band_name: str) -> Interval:
        """Return the span covered by a band designation.

        A designation may be less specific than the table, for example ``q13`` when the
        table stores ``q13.1``, ``q13.2`` and ``q13.3``. The union of all matching
        sub-bands is returned, which is the correct reading: the cytogeneticist named a
        region, not a point.
        """
        if BAND_PATTERN.fullmatch(band_name) is None:
            raise ValueError(f"unsupported band designation: {band_name!r}")
        matches = [
            band
            for band in self.contig_bands(contig)
            if band.name == band_name or band.name.startswith(f"{band_name}.")
        ]
        if not matches:
            raise KeyError(f"band {band_name} not found on contig {contig}")
        return min(band.start for band in matches), max(band.end for band in matches)

    def band_span(self, contig: str, first_band: str, second_band: str) -> Interval:
        """Return the span between two band designations, in either order."""
        left = self.band_interval(contig, first_band)
        right = self.band_interval(contig, second_band)
        return min(left[0], right[0]), max(left[1], right[1])

    def band_uncertainty(self, contig: str, band_name: str) -> int:
        """Return the width of a band, used as breakpoint uncertainty.

        A karyotype breakpoint written as ``q13`` is only localised to that band, so the
        band width is the honest uncertainty on the coordinate. Reporting a breakpoint
        error smaller than this against karyotype truth would be measuring the band, not
        the caller.
        """
        start, end = self.band_interval(contig, band_name)
        return end - start


def parse_ucsc_cytoband(
    lines: Sequence[str],
    *,
    genome_build: GenomeBuild,
    resource_id: str,
    source: str,
    source_sha256: str | None = None,
) -> CytobandTable:
    """Parse a UCSC ``cytoBand.txt`` table.

    The format is five tab-separated columns: chromosome, start, end, band name and
    Giemsa stain. Non-canonical contigs are skipped rather than rejected, because a real
    UCSC file contains alternate haplotypes and unplaced scaffolds that carry no
    cytogenetic meaning.
    """
    bands: list[Cytoband] = []
    skipped = 0
    for number, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise ValueError(f"cytoband line {number} has fewer than five columns")
        contig = fields[0]
        try:
            band = Cytoband(
                contig=contig,
                start=int(fields[1]),
                end=int(fields[2]),
                name=fields[3],
                stain=fields[4],
            )
        except ValueError:
            skipped += 1
            continue
        bands.append(band)
    if not bands:
        raise ValueError("cytoband table contains no canonical-contig bands")
    return CytobandTable(
        resource_id=resource_id,
        genome_build=genome_build,
        source=source,
        source_sha256=source_sha256,
        bands=bands,
    )


def load_cytoband_file(
    path: Path,
    *,
    genome_build: GenomeBuild,
    resource_id: str,
    source_sha256: str | None = None,
) -> CytobandTable:
    """Load a locked UCSC-format cytoband file from disk."""
    if not path.is_file():
        raise ValueError(f"cytoband resource is missing or unreadable: {path}")
    return parse_ucsc_cytoband(
        path.read_text(encoding="utf-8").splitlines(),
        genome_build=genome_build,
        resource_id=resource_id,
        source=path.name,
        source_sha256=source_sha256,
    )
