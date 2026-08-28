"""Typed, build-locked breakpoint annotation over the compiled GRCh38 cache.

Coordinates are 0-based half-open throughout.  A breakpoint is represented as a single
0-based position.  Transcript annotations are retained independently for each breakpoint;
the helper does not invent a gene fusion, orientation, or coding frame.
"""

from __future__ import annotations

import gzip
import sqlite3
from array import array
from bisect import bisect_right
from collections import OrderedDict
from collections.abc import Iterable, Mapping, Sequence
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from threading import RLock
from typing import IO, TYPE_CHECKING, Literal, Protocol, TypeAlias, cast, runtime_checkable

from .coordinates import normalize_contig
from .models import (
    BreakpointAnnotation as BreakpointAnnotationModel,
)
from .models import (
    BreakpointContextAnnotation,
    BreakpointTranscriptAnnotation,
    EventType,
    FusionAnnotation,
    FusionPartnerAnnotation,
    GenomicEvent,
)
from .panel_compiler import CachedTranscript

if TYPE_CHECKING:
    from .fusion_evidence import FusionPartnerEvidence


class BreakpointAnnotationError(ValueError):
    """The breakpoint or annotation cache violates the GRCh38 annotation contract."""


@dataclass(frozen=True)
class Breakpoint:
    chromosome: str
    position: int

    def __post_init__(self) -> None:
        if self.position < 0:
            raise BreakpointAnnotationError("breakpoint position must be non-negative")


@dataclass(frozen=True)
class ContextInterval:
    chromosome: str
    start: int
    end: int
    label: str
    value: float | None = None

    def __post_init__(self) -> None:
        if self.start < 0 or self.end <= self.start:
            raise BreakpointAnnotationError("context interval must be 0-based half-open")


@runtime_checkable
class ContextIntervalQuery(Protocol):
    """Query contract accepted alongside the historical interval-sequence API."""

    def overlaps(self, chromosome: str, position: int) -> tuple[ContextInterval, ...]:
        """Return every half-open interval containing one genomic position."""


ContextIntervalSource: TypeAlias = Sequence[ContextInterval] | ContextIntervalQuery

_INDEX_BLOCK_SIZE = 128


@dataclass(frozen=True)
class _ContigIntervalIndex:
    """Compact point-query index for one normalized chromosome.

    Coordinates and label identifiers live in packed arrays rather than one Python object per
    BED row. ``bisect_right`` narrows by start; a small block-max tree then prunes expired
    intervals without losing nested or overlapping records.
    """

    chromosome: str
    starts: array[int]
    ends: array[int]
    label_ids: array[int]
    labels: tuple[str, ...]
    values: array[float] | None
    value_present: bytes | None
    block_leaf_count: int
    block_max_tree: array[int]

    @classmethod
    def build(
        cls,
        chromosome: str,
        records: Iterable[tuple[int, int, str, float | None]],
    ) -> _ContigIntervalIndex:
        starts = array("Q")
        ends = array("Q")
        label_ids = array("I")
        labels: list[str] = []
        label_lookup: dict[str, int] = {}
        values: array[float] | None = None
        value_present: bytearray | None = None
        previous_start: int | None = None
        already_sorted = True

        for start, end, label, value in records:
            if start < 0 or end <= start:
                raise BreakpointAnnotationError("context interval must be 0-based half-open")
            if previous_start is not None and start < previous_start:
                already_sorted = False
            previous_start = start
            starts.append(start)
            ends.append(end)
            label_id = label_lookup.get(label)
            if label_id is None:
                label_id = len(labels)
                label_lookup[label] = label_id
                labels.append(label)
            label_ids.append(label_id)
            if value is not None and values is None:
                values = array("d", [0.0]) * (len(starts) - 1)
                value_present = bytearray(len(starts) - 1)
            if values is not None:
                values.append(float(value) if value is not None else 0.0)
                assert value_present is not None
                value_present.append(value is not None)

        if not already_sorted:
            order = sorted(range(len(starts)), key=starts.__getitem__)
            starts = array("Q", (starts[index] for index in order))
            ends = array("Q", (ends[index] for index in order))
            label_ids = array("I", (label_ids[index] for index in order))
            if values is not None:
                values = array("d", (values[index] for index in order))
                assert value_present is not None
                value_present = bytearray(value_present[index] for index in order)

        block_count = (len(starts) + _INDEX_BLOCK_SIZE - 1) // _INDEX_BLOCK_SIZE
        block_leaf_count = 1 << max(0, (block_count - 1).bit_length())
        block_max_tree = array("Q", [0]) * (2 * block_leaf_count)
        for block in range(block_count):
            first = block * _INDEX_BLOCK_SIZE
            last = min(first + _INDEX_BLOCK_SIZE, len(ends))
            block_max_tree[block_leaf_count + block] = max(ends[first:last])
        for node in range(block_leaf_count - 1, 0, -1):
            block_max_tree[node] = max(block_max_tree[node * 2], block_max_tree[node * 2 + 1])

        return cls(
            chromosome=chromosome,
            starts=starts,
            ends=ends,
            label_ids=label_ids,
            labels=tuple(labels),
            values=values,
            value_present=bytes(value_present) if value_present is not None else None,
            block_leaf_count=block_leaf_count,
            block_max_tree=block_max_tree,
        )

    def _interval(self, index: int) -> ContextInterval:
        value = None
        if self.values is not None:
            assert self.value_present is not None
            if self.value_present[index]:
                value = self.values[index]
        return ContextInterval(
            chromosome=self.chromosome,
            start=self.starts[index],
            end=self.ends[index],
            label=self.labels[self.label_ids[index]],
            value=value,
        )

    def overlaps(self, position: int) -> tuple[ContextInterval, ...]:
        upper = bisect_right(self.starts, position)
        if upper == 0:
            return ()

        hits: list[ContextInterval] = []

        def scan_block(block: int, stop: int | None = None) -> None:
            first = block * _INDEX_BLOCK_SIZE
            last = min(first + _INDEX_BLOCK_SIZE, len(self.starts), stop or len(self.starts))
            for index in range(first, last):
                if self.ends[index] > position:
                    hits.append(self._interval(index))

        complete_blocks = upper // _INDEX_BLOCK_SIZE

        def visit(node: int, left: int, right: int) -> None:
            if left >= complete_blocks or self.block_max_tree[node] <= position:
                return
            if right - left == 1:
                scan_block(left)
                return
            middle = (left + right) // 2
            visit(node * 2, left, middle)
            visit(node * 2 + 1, middle, right)

        if complete_blocks:
            visit(1, 0, self.block_leaf_count)
        if upper % _INDEX_BLOCK_SIZE:
            scan_block(complete_blocks, upper)
        return tuple(hits)


class ContextIntervalIndex:
    """Reusable compact index for an existing sequence of context intervals."""

    def __init__(self, intervals: Iterable[ContextInterval]) -> None:
        by_contig: dict[str, list[tuple[int, int, str, float | None]]] = {}
        for interval in intervals:
            chromosome = normalize_contig(interval.chromosome)
            by_contig.setdefault(chromosome, []).append(
                (interval.start, interval.end, interval.label, interval.value)
            )
        self._contigs = {
            chromosome: _ContigIntervalIndex.build(chromosome, records)
            for chromosome, records in by_contig.items()
        }

    @property
    def contig_interval_counts(self) -> Mapping[str, int]:
        return {chromosome: len(index.starts) for chromosome, index in self._contigs.items()}

    def overlaps(self, chromosome: str, position: int) -> tuple[ContextInterval, ...]:
        index = self._contigs.get(normalize_contig(chromosome))
        return index.overlaps(position) if index is not None else ()


@dataclass(frozen=True)
class _BedByteRange:
    start: int
    end: int


def _open_binary_bed(path: Path) -> IO[bytes]:
    if path.suffix.casefold() == ".gz":
        return cast(IO[bytes], gzip.open(path, mode="rb"))
    return path.open(mode="rb")


def _parse_bed_record(
    raw: bytes,
    *,
    resource_type: str,
    line_number: int,
) -> tuple[str, int, int, str] | None:
    line = raw.decode("utf-8").rstrip("\r\n")
    if not line or line.startswith(("#", "track", "browser")):
        return None
    fields = line.split("\t")
    if len(fields) < 3:
        raise ValueError(f"{resource_type} line {line_number}: expected BED chrom/start/end")
    try:
        start, end = int(fields[1]), int(fields[2])
    except ValueError as exc:
        raise ValueError(
            f"{resource_type} line {line_number}: non-integer BED coordinates"
        ) from exc
    if start < 0 or end <= start:
        raise ValueError(f"{resource_type} line {line_number}: invalid half-open BED interval")
    label = fields[3].strip() if len(fields) > 3 and fields[3].strip() else resource_type
    return normalize_contig(fields[0]), start, end, label


class PathBackedContextIntervalIndex:
    """Validated BED index that materializes at most a few queried contigs.

    The initial pass retains only byte ranges and row counts. Uncompressed byte offsets and gzip
    uncompressed offsets are both seekable through their respective Python streams.
    """

    def __init__(
        self,
        path: Path,
        *,
        resource_type: str,
        max_cached_contigs: int = 2,
    ) -> None:
        if max_cached_contigs < 1:
            raise ValueError("max_cached_contigs must be positive")
        self.path = path.resolve()
        self.resource_type = resource_type
        self._max_cached_contigs = max_cached_contigs
        self._ranges, self._counts, self._signature = self._scan()
        self._cache: OrderedDict[str, _ContigIntervalIndex] = OrderedDict()
        self._cache_lock = RLock()

    @staticmethod
    def _stat_signature(path: Path) -> tuple[int, int, int]:
        stat = path.stat()
        return stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns

    def _scan(
        self,
    ) -> tuple[
        Mapping[str, tuple[_BedByteRange, ...]],
        Mapping[str, int],
        tuple[int, int, int],
    ]:
        if not self.path.is_file():
            raise ValueError(f"annotation resource is missing: {self.path}")
        before = self._stat_signature(self.path)
        ranges: dict[str, list[_BedByteRange]] = {}
        counts: dict[str, int] = {}
        with _open_binary_bed(self.path) as handle:
            line_number = 0
            while True:
                offset = handle.tell()
                raw = handle.readline()
                if not raw:
                    break
                line_number += 1
                record = _parse_bed_record(
                    raw,
                    resource_type=self.resource_type,
                    line_number=line_number,
                )
                if record is None:
                    continue
                chromosome = record[0]
                end_offset = handle.tell()
                contig_ranges = ranges.setdefault(chromosome, [])
                if contig_ranges and contig_ranges[-1].end == offset:
                    previous = contig_ranges[-1]
                    contig_ranges[-1] = _BedByteRange(previous.start, end_offset)
                else:
                    contig_ranges.append(_BedByteRange(offset, end_offset))
                counts[chromosome] = counts.get(chromosome, 0) + 1
        after = self._stat_signature(self.path)
        if before != after:
            raise ValueError(f"{self.resource_type} BED changed while it was being indexed")
        return (
            {chromosome: tuple(items) for chromosome, items in ranges.items()},
            counts,
            after,
        )

    @property
    def contig_interval_counts(self) -> Mapping[str, int]:
        return dict(self._counts)

    @property
    def cached_contigs(self) -> tuple[str, ...]:
        with self._cache_lock:
            return tuple(self._cache)

    def _load_contig(self, chromosome: str) -> _ContigIntervalIndex:
        def records() -> Iterable[tuple[int, int, str, float | None]]:
            with _open_binary_bed(self.path) as handle:
                for byte_range in self._ranges[chromosome]:
                    handle.seek(byte_range.start)
                    while handle.tell() < byte_range.end:
                        raw = handle.readline()
                        if not raw:
                            raise ValueError(
                                f"{self.resource_type} BED ended while loading {chromosome}"
                            )
                        record = _parse_bed_record(
                            raw,
                            resource_type=self.resource_type,
                            line_number=0,
                        )
                        if record is None or record[0] != chromosome:
                            raise ValueError(
                                f"{self.resource_type} BED index no longer matches its source"
                            )
                        _, start, end, label = record
                        yield start, end, label, None

        return _ContigIntervalIndex.build(chromosome, records())

    def overlaps(self, chromosome: str, position: int) -> tuple[ContextInterval, ...]:
        normalized = normalize_contig(chromosome)
        if normalized not in self._ranges:
            return ()
        with self._cache_lock:
            cached = self._cache.get(normalized)
            if cached is not None:
                self._cache.move_to_end(normalized)
                return cached.overlaps(position)
            if self._stat_signature(self.path) != self._signature:
                raise ValueError(f"{self.resource_type} BED changed after it was indexed")
            loaded = self._load_contig(normalized)
            self._cache[normalized] = loaded
            self._cache.move_to_end(normalized)
            while len(self._cache) > self._max_cached_contigs:
                self._cache.popitem(last=False)
            return loaded.overlaps(position)


@dataclass(frozen=True)
class ContextHit:
    resource_type: str
    label: str
    value: float | None


@dataclass(frozen=True)
class TranscriptBreakpointHit:
    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str | None
    strand: str
    preferred: bool
    rank_tier: int
    region: Literal["exon", "intron", "transcript"]
    exon_number: int | None
    intron_number: int | None
    cds_phase: int | None


@dataclass(frozen=True)
class AnnotatedBreakpoint:
    breakpoint: Breakpoint
    cytoband: str | None
    transcripts: tuple[TranscriptBreakpointHit, ...]
    contexts: tuple[ContextHit, ...]

    @property
    def genes(self) -> tuple[str, ...]:
        return tuple(sorted({hit.gene_name for hit in self.transcripts}))


@dataclass(frozen=True)
class AnnotatedBreakpointPair:
    primary: AnnotatedBreakpoint
    secondary: AnnotatedBreakpoint | None


def _require_grch38_cache(connection: sqlite3.Connection) -> None:
    try:
        row = connection.execute("SELECT value FROM metadata WHERE key = 'genome_build'").fetchone()
    except sqlite3.Error as exc:
        raise BreakpointAnnotationError("annotation cache metadata is unreadable") from exc
    observed = str(row[0]) if row is not None else "unspecified"
    if observed != "GRCh38":
        raise BreakpointAnnotationError(
            f"breakpoint annotation requires GRCh38 cache, observed {observed!r}"
        )


def _cached_transcript(row: sqlite3.Row) -> CachedTranscript:
    return CachedTranscript(
        gene_id=str(row[0]),
        transcript_id=str(row[1]),
        transcript_name=str(row[2]) if row[2] is not None else None,
        chromosome=str(row[3]),
        start=int(row[4]),
        end=int(row[5]),
        strand=str(row[6]),
        transcript_type=str(row[7]) if row[7] is not None else None,
        mane_status=str(row[8]) if row[8] is not None else None,
        mane_refseq_id=str(row[9]) if row[9] is not None else None,
        appris=str(row[10]) if row[10] is not None else None,
        is_canonical=bool(row[11]),
        is_basic=bool(row[12]),
        cds_length=int(row[13] or 0),
        transcript_length=int(row[14] or 0),
    )


def _all_gene_transcripts(
    connection: sqlite3.Connection, gene_id: str
) -> tuple[CachedTranscript, ...]:
    rows = connection.execute(
        """
        SELECT gene_id, transcript_id, transcript_name, chrom, start, end, strand,
               transcript_type, mane_status, mane_refseq_id, appris, is_canonical,
               is_basic, cds_length, transcript_length
        FROM transcripts WHERE gene_id = ?
        """,
        (gene_id,),
    ).fetchall()
    return tuple(sorted((_cached_transcript(row) for row in rows), key=lambda item: item.sort_key))


def _exon_or_intron(
    connection: sqlite3.Connection,
    transcript_id: str,
    position: int,
) -> tuple[Literal["exon", "intron", "transcript"], int | None, int | None]:
    rows = connection.execute(
        """
        SELECT exon_number, start, end
        FROM exons
        WHERE transcript_id = ?
        ORDER BY exon_number, start, end
        """,
        (transcript_id,),
    ).fetchall()
    exons = tuple((int(number), int(start), int(end)) for number, start, end in rows)
    for number, start, end in exons:
        if start <= position < end:
            return "exon", number, None
    by_number = sorted(exons)
    for (left_number, left_start, left_end), (
        right_number,
        right_start,
        right_end,
    ) in zip(by_number, by_number[1:], strict=False):
        gap_start = min(left_end, right_end)
        gap_end = max(left_start, right_start)
        if gap_start <= position < gap_end:
            return "intron", None, min(left_number, right_number)
    return "transcript", None, None


def _cds_phase(
    connection: sqlite3.Connection,
    transcript_id: str,
    position: int,
) -> int | None:
    rows = connection.execute(
        """
        SELECT phase FROM cds
        WHERE transcript_id = ? AND start <= ? AND ? < end
        ORDER BY start, end
        """,
        (transcript_id, position, position),
    ).fetchall()
    phases = {int(row[0]) for row in rows if row[0] is not None}
    if len(phases) == 1:
        return next(iter(phases))
    return None


def _transcript_hits(
    connection: sqlite3.Connection,
    chromosome: str,
    position: int,
) -> tuple[TranscriptBreakpointHit, ...]:
    try:
        rows = connection.execute(
            """
            SELECT t.gene_id, t.transcript_id, t.transcript_name, t.chrom, t.start, t.end,
                   t.strand, t.transcript_type, t.mane_status, t.mane_refseq_id, t.appris,
                   t.is_canonical, t.is_basic, t.cds_length, t.transcript_length, g.gene_name
            FROM transcripts AS t
            JOIN genes AS g ON g.gene_id = t.gene_id
            WHERE t.chrom = ? AND t.start <= ? AND ? < t.end
            ORDER BY g.gene_name, t.transcript_id
            """,
            (chromosome, position, position),
        ).fetchall()
    except sqlite3.Error as exc:
        raise BreakpointAnnotationError(
            "annotation cache transcript tables are unreadable"
        ) from exc
    preferred_by_gene: dict[str, str] = {}
    for row in rows:
        gene_id = str(row[0])
        if gene_id not in preferred_by_gene:
            all_transcripts = _all_gene_transcripts(connection, gene_id)
            if all_transcripts:
                preferred_by_gene[gene_id] = all_transcripts[0].transcript_id
    hits: list[TranscriptBreakpointHit] = []
    for row in rows:
        transcript = _cached_transcript(row)
        region, exon_number, intron_number = _exon_or_intron(
            connection, transcript.transcript_id, position
        )
        hits.append(
            TranscriptBreakpointHit(
                gene_id=transcript.gene_id,
                gene_name=str(row[15]),
                transcript_id=transcript.transcript_id,
                transcript_name=transcript.transcript_name,
                strand=transcript.strand,
                preferred=preferred_by_gene.get(transcript.gene_id) == transcript.transcript_id,
                rank_tier=transcript.rank_tier,
                region=region,
                exon_number=exon_number,
                intron_number=intron_number,
                cds_phase=_cds_phase(connection, transcript.transcript_id, position),
            )
        )
    return tuple(
        sorted(
            hits,
            key=lambda hit: (
                hit.gene_name,
                not hit.preferred,
                hit.rank_tier,
                hit.transcript_id,
            ),
        )
    )


def _cytoband(connection: sqlite3.Connection, chromosome: str, position: int) -> str | None:
    try:
        rows = connection.execute(
            """
            SELECT name FROM cytobands
            WHERE chrom = ? AND start <= ? AND ? < end
            ORDER BY start, end, name
            """,
            (chromosome, position, position),
        ).fetchall()
    except sqlite3.Error as exc:
        raise BreakpointAnnotationError("annotation cache cytobands table is unreadable") from exc
    names = tuple(dict.fromkeys(str(row[0]) for row in rows))
    return names[0] if len(names) == 1 else None


def _context_hits(
    breakpoint: Breakpoint,
    context_resources: Mapping[str, ContextIntervalQuery],
) -> tuple[ContextHit, ...]:
    chromosome = normalize_contig(breakpoint.chromosome)
    hits = [
        ContextHit(resource_type, interval.label, interval.value)
        for resource_type, index in context_resources.items()
        for interval in index.overlaps(chromosome, breakpoint.position)
    ]
    return tuple(sorted(hits, key=lambda hit: (hit.resource_type, hit.label, hit.value or 0.0)))


def _index_context_resources(
    context_resources: Mapping[str, ContextIntervalSource] | None,
) -> dict[str, ContextIntervalQuery]:
    return {
        resource_type: (
            intervals
            if isinstance(intervals, ContextIntervalQuery)
            else ContextIntervalIndex(intervals)
        )
        for resource_type, intervals in (context_resources or {}).items()
    }


def annotate_breakpoint_pair(
    annotation_cache: Path,
    primary: Breakpoint,
    secondary: Breakpoint | None = None,
    *,
    context_resources: Mapping[str, ContextIntervalSource] | None = None,
) -> AnnotatedBreakpointPair:
    """Annotate both SV breakpoints without collapsing their independent evidence."""
    contexts = _index_context_resources(context_resources)
    try:
        with closing(sqlite3.connect(annotation_cache)) as connection:
            _require_grch38_cache(connection)

            def annotate(breakpoint: Breakpoint) -> AnnotatedBreakpoint:
                chromosome = normalize_contig(breakpoint.chromosome)
                normalized = Breakpoint(chromosome, breakpoint.position)
                return AnnotatedBreakpoint(
                    breakpoint=normalized,
                    cytoband=_cytoband(connection, chromosome, breakpoint.position),
                    transcripts=_transcript_hits(connection, chromosome, breakpoint.position),
                    contexts=_context_hits(normalized, contexts),
                )

            return AnnotatedBreakpointPair(
                primary=annotate(primary),
                secondary=annotate(secondary) if secondary is not None else None,
            )
    except sqlite3.Error as exc:
        raise BreakpointAnnotationError(
            f"annotation cache cannot be read: {annotation_cache}"
        ) from exc


def _breakpoint_model(
    label: Literal["primary", "secondary"],
    annotation: AnnotatedBreakpoint,
) -> BreakpointAnnotationModel:
    return BreakpointAnnotationModel(
        label=label,
        chromosome=annotation.breakpoint.chromosome,
        position=annotation.breakpoint.position,
        cytoband=annotation.cytoband,
        transcripts=[
            BreakpointTranscriptAnnotation(
                gene_id=hit.gene_id,
                gene_name=hit.gene_name,
                transcript_id=hit.transcript_id,
                transcript_name=hit.transcript_name,
                strand=hit.strand,
                preferred=hit.preferred,
                rank_tier=hit.rank_tier,
                region=hit.region,
                exon_number=hit.exon_number,
                intron_number=hit.intron_number,
                cds_phase=hit.cds_phase,
            )
            for hit in annotation.transcripts
        ],
        contexts=[
            BreakpointContextAnnotation(
                resource_type=hit.resource_type,
                label=hit.label,
                value=hit.value,
            )
            for hit in annotation.contexts
        ],
    )


def _fusion_model(pair: AnnotatedBreakpointPair, orientation: str | None) -> FusionAnnotation:
    from .fusion_evidence import fusion_evidence_from_breakpoints

    evidence = fusion_evidence_from_breakpoints(pair, orientation=orientation)

    def partner(value: FusionPartnerEvidence) -> FusionPartnerAnnotation:
        # ``FusionPartnerEvidence`` is deliberately a dataclass rather than a Pydantic
        # output contract. Keeping this conversion here prevents that implementation type
        # from leaking into PipelineResult.
        return FusionPartnerAnnotation(
            gene=value.gene,
            preferred_transcript=value.preferred_transcript,
            region=value.region,
            exon_number=value.exon_number,
            intron_number=value.intron_number,
            strand=value.strand,
        )

    return FusionAnnotation(
        gene_a=partner(evidence.gene_a),
        gene_b=partner(evidence.gene_b),
        orientation=evidence.orientation,
        frame_status=evidence.frame_status.value,
    )


def annotate_events_from_cache(
    events: Sequence[GenomicEvent],
    annotation_cache: Path,
    *,
    context_resources: Mapping[str, ContextIntervalSource] | None = None,
) -> list[GenomicEvent]:
    """Attach independent, transcript-aware annotations to every available breakpoint.

    BND/translocation records already carry an explicit secondary locus. For an
    intrachromosomal deletion, duplication or inversion, the second breakpoint is the final
    covered base of the primary half-open span. Insertions have only one reference-side
    coordinate and therefore retain one breakpoint annotation.
    """

    result: list[GenomicEvent] = []
    indexed_contexts = _index_context_resources(context_resources)
    paired_span_types = {EventType.DELETION, EventType.DUPLICATION, EventType.INVERSION}
    for event in events:
        secondary: Breakpoint | None = None
        if event.secondary is not None:
            secondary = Breakpoint(event.secondary.chromosome, event.secondary.start)
        elif event.event_type in paired_span_types and event.primary.end > event.primary.start:
            secondary = Breakpoint(event.primary.chromosome, event.primary.end - 1)
        pair = annotate_breakpoint_pair(
            annotation_cache,
            Breakpoint(event.primary.chromosome, event.primary.start),
            secondary,
            context_resources=indexed_contexts,
        )
        breakpoint_models = [_breakpoint_model("primary", pair.primary)]
        if pair.secondary is not None:
            breakpoint_models.append(_breakpoint_model("secondary", pair.secondary))

        genes = sorted(
            set(event.genes)
            .union(pair.primary.genes)
            .union(pair.secondary.genes if pair.secondary is not None else ())
        )
        primary = event.primary.model_copy(
            update={
                "gene": pair.primary.genes[0] if len(pair.primary.genes) == 1 else None,
                "cytoband_start": pair.primary.cytoband,
                "cytoband_end": (
                    pair.secondary.cytoband
                    if pair.secondary is not None
                    and pair.secondary.breakpoint.chromosome == pair.primary.breakpoint.chromosome
                    else pair.primary.cytoband
                ),
            }
        )
        secondary_locus = event.secondary
        if event.secondary is not None and pair.secondary is not None:
            secondary_locus = event.secondary.model_copy(
                update={
                    "gene": (pair.secondary.genes[0] if len(pair.secondary.genes) == 1 else None),
                    "cytoband_start": pair.secondary.cytoband,
                    "cytoband_end": pair.secondary.cytoband,
                }
            )
        context_flags = [
            f"{annotation.label}:{hit.resource_type}"
            for annotation in breakpoint_models
            for hit in annotation.contexts
        ]
        orientation = next(
            (
                evidence.supporting_read_strands
                for evidence in event.evidence
                if evidence.supporting_read_strands is not None
            ),
            None,
        )
        fusion = (
            _fusion_model(pair, orientation)
            if event.event_type in {EventType.TRANSLOCATION, EventType.FUSION}
            and pair.secondary is not None
            else None
        )
        result.append(
            event.model_copy(
                update={
                    "primary": primary,
                    "secondary": secondary_locus,
                    "genes": genes,
                    "technical_flags": sorted(set([*event.technical_flags, *context_flags])),
                    "breakpoint_annotations": breakpoint_models,
                    "fusion_evidence": fusion,
                    "reportable": False,
                }
            )
        )
    return result
