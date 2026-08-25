# ruff: noqa: I001
from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path

from .models import GenomeBuild


_VALID_CHROMOSOMES = {str(value) for value in range(1, 23)} | {"X", "Y"}


def normalize_chromosome(chromosome: str) -> str:
    value = chromosome.removeprefix("chr")
    if value not in _VALID_CHROMOSOMES:
        raise ValueError(f"unsupported chromosome: {chromosome}")
    return f"chr{value}"


@dataclass(frozen=True, slots=True)
class CytobandRecord:
    genome_build: GenomeBuild
    chromosome: str
    start: int
    end: int
    name: str
    stain: str

    @property
    def arm(self) -> str:
        if not self.name or self.name[0] not in {"p", "q"}:
            raise ValueError(f"invalid cytoband name: {self.name}")
        return self.name[0]


class CytobandIndex:
    """Immutable in-memory index for UCSC-style cytoband tables.

    Coordinates are interpreted as zero-based, half-open intervals, matching BED/UCSC.
    Only chromosomes 1-22, X and Y are loaded. The index intentionally contains no ISCN
    prose: it is a coordinate reference layer, not a copy of the ISCN publication.
    """

    def __init__(self, genome_build: GenomeBuild, records: list[CytobandRecord]) -> None:
        self.genome_build = genome_build
        grouped: dict[str, list[CytobandRecord]] = {}
        for record in records:
            if record.genome_build != genome_build:
                raise ValueError("cytoband record build does not match index build")
            grouped.setdefault(record.chromosome, []).append(record)

        self._records: dict[str, tuple[CytobandRecord, ...]] = {}
        self._starts: dict[str, tuple[int, ...]] = {}
        for chromosome, chromosome_records in grouped.items():
            ordered = tuple(sorted(chromosome_records, key=lambda item: item.start))
            self._validate_chromosome(chromosome, ordered)
            self._records[chromosome] = ordered
            self._starts[chromosome] = tuple(item.start for item in ordered)

    @staticmethod
    def _validate_chromosome(chromosome: str, records: tuple[CytobandRecord, ...]) -> None:
        previous_end: int | None = None
        for record in records:
            if record.start < 0 or record.end <= record.start:
                raise ValueError(f"invalid cytoband interval: {record}")
            if record.chromosome != chromosome:
                raise ValueError("cytoband grouping mismatch")
            if previous_end is not None and record.start < previous_end:
                raise ValueError(f"overlapping cytobands on {chromosome}")
            previous_end = record.end

    @classmethod
    def from_tsv_text(
        cls,
        text: str,
        genome_build: GenomeBuild,
    ) -> CytobandIndex:
        records: list[CytobandRecord] = []
        for line_number, raw_line in enumerate(text.splitlines(), start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"cytoband line {line_number} has fewer than five columns")
            chromosome, start_text, end_text, name, stain = fields[:5]
            try:
                chromosome = normalize_chromosome(chromosome)
            except ValueError:
                continue
            records.append(
                CytobandRecord(
                    genome_build=genome_build,
                    chromosome=chromosome,
                    start=int(start_text),
                    end=int(end_text),
                    name=name,
                    stain=stain,
                )
            )
        if not records:
            raise ValueError("cytoband table contains no supported human chromosomes")
        return cls(genome_build=genome_build, records=records)

    @classmethod
    def from_tsv_path(cls, path: str | Path, genome_build: GenomeBuild) -> CytobandIndex:
        return cls.from_tsv_text(Path(path).read_text(encoding="utf-8"), genome_build)

    def band_at(self, chromosome: str, position: int) -> CytobandRecord | None:
        if position < 0:
            raise ValueError("position must be non-negative")
        chrom = normalize_chromosome(chromosome)
        records = self._records.get(chrom)
        starts = self._starts.get(chrom)
        if not records or not starts:
            return None
        index = bisect_right(starts, position) - 1
        if index < 0:
            return None
        record = records[index]
        if record.start <= position < record.end:
            return record
        return None

    def bands_for_interval(
        self,
        chromosome: str,
        start: int,
        end: int,
    ) -> tuple[CytobandRecord, CytobandRecord] | None:
        if end <= start:
            raise ValueError("interval end must be greater than start")
        first = self.band_at(chromosome, start)
        last = self.band_at(chromosome, end - 1)
        if first is None or last is None:
            return None
        return first, last

    def centromere_bounds(self, chromosome: str) -> tuple[int, int] | None:
        chrom = normalize_chromosome(chromosome)
        records = self._records.get(chrom, ())
        acen = [record for record in records if record.stain == "acen"]
        if not acen:
            return None
        return min(record.start for record in acen), max(record.end for record in acen)

    def record_count(self) -> int:
        return sum(len(records) for records in self._records.values())
