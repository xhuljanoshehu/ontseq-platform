"""Deterministic local SQLite cache for GRCh38 transcript and cytoband annotation.

The cache is compiled once while a reference bundle is installed.  Analysis code opens the
finished database read-only; it never scans the full GENCODE or MANE sources and never needs an
internet connection.  Every genomic interval in this module is zero-based and half-open.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import sqlite3
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from urllib.parse import unquote

from .coordinates import normalize_contig, to_zero_based_half_open
from .models import CoordinateSystem
from .reference import sha256_file

ANNOTATION_CACHE_SCHEMA_VERSION = "1.0.0"

_CANONICAL_CONTIG = re.compile(r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$")
_REFSEQ_TRANSCRIPT = re.compile(r"\b(?:N[MR]|X[MR])_\d+(?:\.\d+)?\b")
_PRINCIPAL_APPRIS = re.compile(r"^appris_principal(?:_\d+)?$", re.IGNORECASE)


@dataclass(frozen=True)
class AnnotationCacheSummary:
    path: Path
    sha256: str
    genes: int
    transcripts: int
    exons: int
    cds: int
    cytobands: int
    mane_matched_transcripts: int


@dataclass(frozen=True)
class RankedTranscript:
    transcript_id: str
    gene_id: str
    gene_name: str
    transcript_name: str | None
    chromosome: str
    start: int
    end: int
    strand: str
    transcript_type: str | None
    tags: tuple[str, ...]
    mane_status: str | None
    mane_refseq_id: str | None
    appris: str | None
    is_canonical: bool
    is_basic: bool
    cds_length: int
    transcript_length: int


@dataclass(frozen=True)
class _Gene:
    gene_id: str
    gene_name: str
    chromosome: str
    start: int
    end: int
    strand: str
    gene_type: str | None


@dataclass(frozen=True)
class _Transcript:
    transcript_id: str
    gene_id: str
    transcript_name: str | None
    chromosome: str
    start: int
    end: int
    strand: str
    transcript_type: str | None
    tags: tuple[str, ...]
    mane_status: str | None
    mane_refseq_id: str | None
    appris: str | None
    is_canonical: bool
    is_basic: bool
    cds_length: int = 0
    transcript_length: int = 0


@dataclass(frozen=True)
class _Exon:
    transcript_id: str
    exon_id: str | None
    exon_number: int | None
    chromosome: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class _Cds:
    transcript_id: str
    exon_id: str | None
    chromosome: str
    start: int
    end: int
    strand: str
    phase: int | None


@dataclass(frozen=True)
class _Cytoband:
    chromosome: str
    start: int
    end: int
    name: str
    gie_stain: str


@dataclass(frozen=True)
class _ManeRecord:
    transcript_id: str
    status: str
    refseq_id: str | None
    chromosome: str
    start: int
    end: int
    strand: str


def _text_lines(path: Path) -> Iterator[str]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        yield from handle


def _normalize_chromosome(raw: str) -> str | None:
    value = normalize_contig(raw)
    if _CANONICAL_CONTIG.fullmatch(value) is None:
        return None
    return value


def _one_based_inclusive(start: str, end: str, *, source: str) -> tuple[int, int]:
    try:
        source_start = int(start)
        source_end = int(end)
    except ValueError as exc:
        raise ValueError(f"{source}: coordinates must be integers") from exc
    try:
        return to_zero_based_half_open(
            source_start, source_end, CoordinateSystem.ONE_BASED_INCLUSIVE
        )
    except ValueError as exc:
        raise ValueError(f"{source}: invalid one-based inclusive interval {start}-{end}") from exc


def _parse_gtf_attributes(raw: str) -> dict[str, list[str]]:
    attributes: dict[str, list[str]] = defaultdict(list)
    for item in raw.rstrip().rstrip(";").split(";"):
        stripped = item.strip()
        if not stripped:
            continue
        key, separator, value = stripped.partition(" ")
        if not separator:
            key, separator, value = stripped.partition("=")
        if not separator:
            continue
        attributes[key].append(value.strip().strip('"'))
    return dict(attributes)


def _parse_gff3_attributes(raw: str) -> dict[str, list[str]]:
    attributes: dict[str, list[str]] = defaultdict(list)
    for item in raw.rstrip().rstrip(";").split(";"):
        key, separator, value = item.partition("=")
        if not separator:
            continue
        for part in value.split(","):
            attributes[unquote(key.strip())].append(unquote(part.strip()))
    return dict(attributes)


def _first(attributes: Mapping[str, Sequence[str]], *keys: str) -> str | None:
    for key in keys:
        values = attributes.get(key)
        if values and values[0]:
            return values[0]
    return None


def _tags(attributes: Mapping[str, Sequence[str]]) -> tuple[str, ...]:
    return tuple(sorted(set((*attributes.get("tag", ()), *attributes.get("Tag", ())))))


def _mane_status(values: Iterable[str]) -> str | None:
    normalized = {re.sub(r"[^a-z]+", "_", value.casefold()).strip("_") for value in values}
    if any("mane_select" in value for value in normalized):
        return "MANE Select"
    if any("mane_plus_clinical" in value for value in normalized):
        return "MANE Plus Clinical"
    return None


def _stable_tag_data(attributes: Mapping[str, Sequence[str]]) -> tuple[tuple[str, ...], str | None]:
    tags = _tags(attributes)
    appris = next((tag for tag in tags if tag.casefold().startswith("appris_")), None)
    return tags, appris


def _parse_gencode(
    path: Path,
) -> tuple[list[_Gene], list[_Transcript], list[_Exon], list[_Cds]]:
    genes: dict[str, _Gene] = {}
    transcripts: dict[str, _Transcript] = {}
    exons: list[_Exon] = []
    cds: list[_Cds] = []

    for line_number, raw_line in enumerate(_text_lines(path), start=1):
        if not raw_line.strip() or raw_line.startswith("#"):
            continue
        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"GENCODE line {line_number}: expected nine GTF columns")
        chromosome = _normalize_chromosome(fields[0])
        if chromosome is None or fields[2] not in {"gene", "transcript", "exon", "CDS"}:
            continue
        start, end = _one_based_inclusive(
            fields[3], fields[4], source=f"GENCODE line {line_number}"
        )
        strand = fields[6]
        if strand not in {"+", "-"}:
            raise ValueError(f"GENCODE line {line_number}: unsupported strand {strand!r}")
        attributes = _parse_gtf_attributes(fields[8])
        gene_id = _first(attributes, "gene_id")
        if not gene_id:
            raise ValueError(f"GENCODE line {line_number}: missing gene_id")

        if fields[2] == "gene":
            gene_name = _first(attributes, "gene_name") or gene_id
            genes[gene_id] = _Gene(
                gene_id=gene_id,
                gene_name=gene_name,
                chromosome=chromosome,
                start=start,
                end=end,
                strand=strand,
                gene_type=_first(attributes, "gene_type", "gene_biotype"),
            )
            continue

        transcript_id = _first(attributes, "transcript_id")
        if not transcript_id:
            raise ValueError(f"GENCODE line {line_number}: {fields[2]} record has no transcript_id")
        if fields[2] == "transcript":
            tags, appris = _stable_tag_data(attributes)
            transcripts[transcript_id] = _Transcript(
                transcript_id=transcript_id,
                gene_id=gene_id,
                transcript_name=_first(attributes, "transcript_name"),
                chromosome=chromosome,
                start=start,
                end=end,
                strand=strand,
                transcript_type=_first(attributes, "transcript_type", "transcript_biotype"),
                tags=tags,
                mane_status=_mane_status(tags),
                mane_refseq_id=None,
                appris=appris,
                is_canonical=any(
                    tag.casefold() in {"ensembl_canonical", "gencode_primary"} for tag in tags
                ),
                is_basic=any(tag.casefold() == "basic" for tag in tags),
            )
        elif fields[2] == "exon":
            raw_exon_number = _first(attributes, "exon_number")
            try:
                exon_number = int(raw_exon_number) if raw_exon_number else None
            except ValueError as exc:
                raise ValueError(
                    f"GENCODE line {line_number}: invalid exon_number {raw_exon_number!r}"
                ) from exc
            exons.append(
                _Exon(
                    transcript_id=transcript_id,
                    exon_id=_first(attributes, "exon_id"),
                    exon_number=exon_number,
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    strand=strand,
                )
            )
        else:
            raw_phase = fields[7]
            if raw_phase not in {".", "0", "1", "2"}:
                raise ValueError(f"GENCODE line {line_number}: invalid CDS phase {raw_phase!r}")
            cds.append(
                _Cds(
                    transcript_id=transcript_id,
                    exon_id=_first(attributes, "exon_id"),
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    strand=strand,
                    phase=None if raw_phase == "." else int(raw_phase),
                )
            )

    if not genes:
        raise ValueError(f"GENCODE source contains no canonical gene records: {path}")
    if not transcripts:
        raise ValueError(f"GENCODE source contains no canonical transcript records: {path}")

    unknown_transcript_ids = sorted(
        {item.transcript_id for item in exons if item.transcript_id not in transcripts}
        | {item.transcript_id for item in cds if item.transcript_id not in transcripts}
    )
    if unknown_transcript_ids:
        raise ValueError(
            "GENCODE exon/CDS records refer to absent transcript records: "
            + ", ".join(unknown_transcript_ids[:5])
        )
    unknown_gene_ids = sorted(
        {item.gene_id for item in transcripts.values() if item.gene_id not in genes}
    )
    if unknown_gene_ids:
        raise ValueError(
            "GENCODE transcript records refer to absent gene records: "
            + ", ".join(unknown_gene_ids[:5])
        )

    exon_lengths: dict[str, int] = defaultdict(int)
    for exon in set(exons):
        exon_lengths[exon.transcript_id] += exon.end - exon.start
    cds_lengths: dict[str, int] = defaultdict(int)
    for coding in set(cds):
        cds_lengths[coding.transcript_id] += coding.end - coding.start
    transcripts = {
        transcript_id: replace(
            item,
            transcript_length=exon_lengths[transcript_id] or item.end - item.start,
            cds_length=cds_lengths[transcript_id],
        )
        for transcript_id, item in transcripts.items()
    }
    return (
        sorted(genes.values(), key=lambda item: item.gene_id),
        sorted(transcripts.values(), key=lambda item: item.transcript_id),
        sorted(
            set(exons),
            key=lambda item: (
                item.transcript_id,
                item.exon_number if item.exon_number is not None else 2**31,
                item.start,
                item.end,
                item.exon_id or "",
            ),
        ),
        sorted(set(cds), key=lambda item: (item.transcript_id, item.start, item.end)),
    )


def _transcript_identifier(attributes: Mapping[str, Sequence[str]]) -> str | None:
    identifier = _first(attributes, "transcript_id", "transcript", "ID")
    if identifier:
        return identifier.removeprefix("transcript:")
    parent = _first(attributes, "Parent")
    return parent.removeprefix("transcript:") if parent else None


def _parse_mane(path: Path) -> list[_ManeRecord]:
    records: dict[str, _ManeRecord] = {}
    for line_number, raw_line in enumerate(_text_lines(path), start=1):
        if not raw_line.strip() or raw_line.startswith("#"):
            continue
        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"MANE line {line_number}: expected nine GFF3 columns")
        if fields[2].casefold() not in {"mrna", "transcript"}:
            continue
        attributes = _parse_gff3_attributes(fields[8])
        chromosome = _normalize_chromosome(fields[0])
        if chromosome is None:
            continue
        start, end = _one_based_inclusive(fields[3], fields[4], source=f"MANE line {line_number}")
        strand = fields[6]
        if strand not in {"+", "-"}:
            raise ValueError(f"MANE line {line_number}: unsupported strand {strand!r}")
        transcript_id = _transcript_identifier(attributes)
        if not transcript_id:
            raise ValueError(f"MANE line {line_number}: transcript has no identifier")
        all_values = [value for values in attributes.values() for value in values]
        status = _mane_status(all_values)
        if status is None:
            continue
        refseq_id = next(
            (
                match.group(0)
                for value in all_values
                if (match := _REFSEQ_TRANSCRIPT.search(value)) is not None
            ),
            None,
        )
        previous = records.get(transcript_id)
        if previous and previous.status == "MANE Select":
            continue
        records[transcript_id] = _ManeRecord(
            transcript_id,
            status,
            refseq_id,
            chromosome,
            start,
            end,
            strand,
        )
    if not records:
        raise ValueError(f"MANE source contains no Select or Plus Clinical transcripts: {path}")
    return sorted(records.values(), key=lambda item: item.transcript_id)


def _apply_mane(
    transcripts: list[_Transcript], mane_records: list[_ManeRecord]
) -> tuple[list[_Transcript], int]:
    full_ids = {item.transcript_id: item.transcript_id for item in transcripts}
    base_ids: dict[str, list[str]] = defaultdict(list)
    for item in transcripts:
        base_ids[item.transcript_id.partition(".")[0]].append(item.transcript_id)

    matched: dict[str, _ManeRecord] = {}
    for mane in mane_records:
        target = full_ids.get(mane.transcript_id)
        if target is None:
            candidates = base_ids.get(mane.transcript_id.partition(".")[0], [])
            if len(candidates) == 1:
                target = candidates[0]
        if target is not None:
            previous = matched.get(target)
            if previous is None or mane.status == "MANE Select":
                matched[target] = mane

    integrated = []
    for item in transcripts:
        matched_mane = matched.get(item.transcript_id)
        if matched_mane is None:
            integrated.append(item)
            continue
        observed_locus = (
            matched_mane.chromosome,
            matched_mane.start,
            matched_mane.end,
            matched_mane.strand,
        )
        expected_locus = (item.chromosome, item.start, item.end, item.strand)
        if observed_locus != expected_locus:
            raise ValueError(
                f"MANE transcript {matched_mane.transcript_id!r} locus does not match "
                f"GENCODE transcript {item.transcript_id!r}"
            )
        tags = tuple(sorted(set((*item.tags, matched_mane.status.replace(" ", "_")))))
        integrated.append(
            replace(
                item,
                tags=tags,
                mane_status=matched_mane.status,
                mane_refseq_id=matched_mane.refseq_id,
            )
        )
    return integrated, len(matched)


def _parse_cytobands(path: Path) -> list[_Cytoband]:
    records: set[_Cytoband] = set()
    for line_number, raw_line in enumerate(_text_lines(path), start=1):
        line = raw_line.rstrip("\r\n")
        if not line or line.startswith(("#", "track ", "browser ")):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise ValueError(f"cytoband line {line_number}: expected five columns")
        chromosome = _normalize_chromosome(fields[0])
        if chromosome is None:
            continue
        try:
            start, end = to_zero_based_half_open(
                int(fields[1]),
                int(fields[2]),
                CoordinateSystem.ZERO_BASED_HALF_OPEN,
            )
        except ValueError as exc:
            raise ValueError(f"cytoband line {line_number}: invalid BED coordinates") from exc
        name = fields[3].strip()
        gie_stain = fields[4].strip()
        if start < 0 or end <= start or not gie_stain:
            raise ValueError(f"cytoband line {line_number}: invalid interval or label")
        # UCSC's native hg38 table includes a chrM placeholder with an empty band name.
        # It is not a cytogenetic band and therefore must not become an annotatable record.
        if not name:
            continue
        records.add(_Cytoband(chromosome, start, end, name, gie_stain))
    if not records:
        raise ValueError(f"cytoband source contains no canonical records: {path}")
    return sorted(records, key=lambda item: (item.chromosome, item.start, item.end, item.name))


_SCHEMA = """
PRAGMA foreign_keys = ON;
CREATE TABLE metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
) WITHOUT ROWID;
CREATE TABLE genes (
    gene_id TEXT PRIMARY KEY,
    gene_name TEXT NOT NULL,
    chrom TEXT NOT NULL,
    start INTEGER NOT NULL CHECK (start >= 0),
    end INTEGER NOT NULL CHECK (end > start),
    strand TEXT NOT NULL CHECK (strand IN ('+', '-')),
    gene_type TEXT
) WITHOUT ROWID;
CREATE TABLE transcripts (
    transcript_id TEXT PRIMARY KEY,
    gene_id TEXT NOT NULL REFERENCES genes(gene_id),
    transcript_name TEXT,
    chrom TEXT NOT NULL,
    start INTEGER NOT NULL CHECK (start >= 0),
    end INTEGER NOT NULL CHECK (end > start),
    strand TEXT NOT NULL CHECK (strand IN ('+', '-')),
    transcript_type TEXT,
    tags_json TEXT NOT NULL,
    mane_status TEXT CHECK (mane_status IN ('MANE Select', 'MANE Plus Clinical')),
    mane_refseq_id TEXT,
    appris TEXT,
    is_canonical INTEGER NOT NULL CHECK (is_canonical IN (0, 1)),
    is_basic INTEGER NOT NULL CHECK (is_basic IN (0, 1)),
    cds_length INTEGER NOT NULL CHECK (cds_length >= 0),
    transcript_length INTEGER NOT NULL CHECK (transcript_length > 0)
) WITHOUT ROWID;
CREATE TABLE exons (
    row_id INTEGER PRIMARY KEY,
    transcript_id TEXT NOT NULL REFERENCES transcripts(transcript_id),
    exon_id TEXT,
    exon_number INTEGER CHECK (exon_number IS NULL OR exon_number > 0),
    chrom TEXT NOT NULL,
    start INTEGER NOT NULL CHECK (start >= 0),
    end INTEGER NOT NULL CHECK (end > start),
    strand TEXT NOT NULL CHECK (strand IN ('+', '-')),
    UNIQUE (transcript_id, start, end, exon_id)
);
CREATE TABLE cds (
    row_id INTEGER PRIMARY KEY,
    transcript_id TEXT NOT NULL REFERENCES transcripts(transcript_id),
    exon_id TEXT,
    chrom TEXT NOT NULL,
    start INTEGER NOT NULL CHECK (start >= 0),
    end INTEGER NOT NULL CHECK (end > start),
    strand TEXT NOT NULL CHECK (strand IN ('+', '-')),
    phase INTEGER CHECK (phase IS NULL OR phase IN (0, 1, 2)),
    UNIQUE (transcript_id, start, end)
);
CREATE TABLE cytobands (
    chrom TEXT NOT NULL,
    start INTEGER NOT NULL CHECK (start >= 0),
    end INTEGER NOT NULL CHECK (end > start),
    name TEXT NOT NULL,
    gie_stain TEXT NOT NULL,
    PRIMARY KEY (chrom, start, end, name)
) WITHOUT ROWID;
CREATE INDEX idx_genes_name ON genes(gene_name);
CREATE INDEX idx_genes_locus ON genes(chrom, start, end);
CREATE INDEX idx_transcripts_gene ON transcripts(gene_id);
CREATE INDEX idx_transcripts_locus ON transcripts(chrom, start, end);
CREATE INDEX idx_exons_transcript ON exons(transcript_id, exon_number, start, end);
CREATE INDEX idx_exons_locus ON exons(chrom, start, end);
CREATE INDEX idx_cds_transcript ON cds(transcript_id, start, end);
CREATE INDEX idx_cds_locus ON cds(chrom, start, end);
CREATE INDEX idx_cytobands_locus ON cytobands(chrom, start, end);
"""


def _insert_cache(
    connection: sqlite3.Connection,
    *,
    genes: list[_Gene],
    transcripts: list[_Transcript],
    exons: list[_Exon],
    cds: list[_Cds],
    cytobands: list[_Cytoband],
    metadata: Mapping[str, str],
) -> None:
    connection.executescript(_SCHEMA)
    complete_metadata = {
        **metadata,
        "schema_version": ANNOTATION_CACHE_SCHEMA_VERSION,
        "coordinate_system": "zero_based_half_open",
        "gencode_coordinate_conversion": (
            "one_based_inclusive_to_zero_based_half_open:start_minus_one"
        ),
        "mane_coordinate_conversion": (
            "one_based_inclusive_to_zero_based_half_open:start_minus_one"
        ),
        "cytoband_coordinate_conversion": "zero_based_half_open:unchanged",
        "contig_normalization": "1-22,X,Y=>chr-prefixed;M,MT,chrM=>chrM",
    }
    connection.executemany(
        "INSERT INTO metadata(key, value) VALUES (?, ?)", sorted(complete_metadata.items())
    )
    connection.executemany(
        "INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?)",
        (
            (
                item.gene_id,
                item.gene_name,
                item.chromosome,
                item.start,
                item.end,
                item.strand,
                item.gene_type,
            )
            for item in genes
        ),
    )
    connection.executemany(
        "INSERT INTO transcripts VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (
            (
                item.transcript_id,
                item.gene_id,
                item.transcript_name,
                item.chromosome,
                item.start,
                item.end,
                item.strand,
                item.transcript_type,
                json.dumps(item.tags, ensure_ascii=False, separators=(",", ":")),
                item.mane_status,
                item.mane_refseq_id,
                item.appris,
                int(item.is_canonical),
                int(item.is_basic),
                item.cds_length,
                item.transcript_length,
            )
            for item in transcripts
        ),
    )
    connection.executemany(
        """INSERT INTO exons(
               transcript_id, exon_id, exon_number, chrom, start, end, strand
           ) VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (
            (
                item.transcript_id,
                item.exon_id,
                item.exon_number,
                item.chromosome,
                item.start,
                item.end,
                item.strand,
            )
            for item in exons
        ),
    )
    connection.executemany(
        """INSERT INTO cds(
               transcript_id, exon_id, chrom, start, end, strand, phase
           ) VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (
            (
                item.transcript_id,
                item.exon_id,
                item.chromosome,
                item.start,
                item.end,
                item.strand,
                item.phase,
            )
            for item in cds
        ),
    )
    connection.executemany(
        "INSERT INTO cytobands VALUES (?, ?, ?, ?, ?)",
        ((item.chromosome, item.start, item.end, item.name, item.gie_stain) for item in cytobands),
    )


def compile_annotation_cache(
    gencode_gtf: Path,
    mane_gff3: Path,
    cytobands: Path,
    output: Path,
    *,
    metadata: Mapping[str, str],
) -> AnnotationCacheSummary:
    """Compile a complete cache to a sibling temporary file and atomically publish it."""

    for label, path in (
        ("GENCODE", gencode_gtf),
        ("MANE", mane_gff3),
        ("cytoband", cytobands),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} source not found: {path}")
    if metadata.get("genome_build") != "GRCh38":
        raise ValueError("the active annotation compiler accepts only genome_build=GRCh38")

    genes, transcripts, exons, cds = _parse_gencode(gencode_gtf)
    mane_records = _parse_mane(mane_gff3)
    transcripts, mane_matched = _apply_mane(transcripts, mane_records)
    bands = _parse_cytobands(cytobands)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent, delete=False
    ) as temporary_handle:
        temporary_path = Path(temporary_handle.name)
    try:
        connection = sqlite3.connect(temporary_path)
        try:
            connection.execute("PRAGMA page_size = 4096")
            connection.execute("PRAGMA journal_mode = OFF")
            connection.execute("PRAGMA synchronous = OFF")
            connection.execute("PRAGMA temp_store = MEMORY")
            with connection:
                _insert_cache(
                    connection,
                    genes=genes,
                    transcripts=transcripts,
                    exons=exons,
                    cds=cds,
                    cytobands=bands,
                    metadata=metadata,
                )
            violations = connection.execute("PRAGMA foreign_key_check").fetchall()
            integrity = connection.execute("PRAGMA integrity_check").fetchone()
            if violations or integrity != ("ok",):
                raise ValueError("compiled annotation cache failed SQLite integrity validation")
            connection.execute("VACUUM")
        finally:
            connection.close()
        os.replace(temporary_path, output)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    return AnnotationCacheSummary(
        path=output,
        sha256=sha256_file(output),
        genes=len(genes),
        transcripts=len(transcripts),
        exons=len(exons),
        cds=len(cds),
        cytobands=len(bands),
        mane_matched_transcripts=mane_matched,
    )


def validate_annotation_cache(path: Path) -> AnnotationCacheSummary:
    """Validate schema, integrity and row counts without modifying the cache."""

    if not path.is_file():
        raise FileNotFoundError(f"annotation cache not found: {path}")
    connection = sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)
    try:
        integrity = connection.execute("PRAGMA integrity_check").fetchone()
        if integrity != ("ok",):
            raise ValueError(f"annotation cache integrity check failed: {integrity!r}")
        version = connection.execute(
            "SELECT value FROM metadata WHERE key = 'schema_version'"
        ).fetchone()
        coordinates = connection.execute(
            "SELECT value FROM metadata WHERE key = 'coordinate_system'"
        ).fetchone()
        build = connection.execute(
            "SELECT value FROM metadata WHERE key = 'genome_build'"
        ).fetchone()
        if version != (ANNOTATION_CACHE_SCHEMA_VERSION,):
            raise ValueError(f"unsupported annotation cache schema: {version!r}")
        if coordinates != ("zero_based_half_open",):
            raise ValueError(f"unsupported annotation coordinate system: {coordinates!r}")
        if build != ("GRCh38",):
            raise ValueError(f"active reference cache must be GRCh38, found: {build!r}")
        counts = {
            table: int(connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
            for table in ("genes", "transcripts", "exons", "cds", "cytobands")
        }
        mane_count = int(
            connection.execute(
                "SELECT COUNT(*) FROM transcripts WHERE mane_status IS NOT NULL"
            ).fetchone()[0]
        )
    except sqlite3.DatabaseError as exc:
        raise ValueError(f"invalid annotation cache {path}: {exc}") from exc
    finally:
        connection.close()
    if counts["genes"] == 0 or counts["transcripts"] == 0 or counts["cytobands"] == 0:
        raise ValueError("annotation cache is incomplete")
    return AnnotationCacheSummary(
        path=path,
        sha256=sha256_file(path),
        genes=counts["genes"],
        transcripts=counts["transcripts"],
        exons=counts["exons"],
        cds=counts["cds"],
        cytobands=counts["cytobands"],
        mane_matched_transcripts=mane_count,
    )


class AnnotationCache:
    """Small read-only query facade over an already validated compiled cache."""

    def __init__(self, path: Path) -> None:
        self.path = path
        validate_annotation_cache(path)

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(f"file:{self.path.as_posix()}?mode=ro", uri=True)
        connection.row_factory = sqlite3.Row
        return connection

    def ranked_transcripts(self, gene_label: str) -> list[RankedTranscript]:
        if not gene_label.strip():
            raise ValueError("gene label may not be empty")
        query = """
            SELECT t.*, g.gene_name
              FROM transcripts AS t
              JOIN genes AS g ON g.gene_id = t.gene_id
             WHERE upper(g.gene_name) = upper(?) OR g.gene_id = ?
             ORDER BY
               CASE
                 WHEN t.mane_status = 'MANE Select' THEN 0
                 WHEN t.mane_status = 'MANE Plus Clinical' THEN 1
                 WHEN t.is_canonical = 1 OR lower(coalesce(t.appris, '')) LIKE 'appris_principal%'
                   THEN 2
                 WHEN t.transcript_type = 'protein_coding' OR t.is_basic = 1 THEN 3
                 ELSE 4
               END,
               t.cds_length DESC,
               t.transcript_length DESC,
               t.transcript_id ASC
        """
        connection = self._connect()
        try:
            rows = connection.execute(query, (gene_label, gene_label)).fetchall()
        finally:
            connection.close()
        return [
            RankedTranscript(
                transcript_id=str(row["transcript_id"]),
                gene_id=str(row["gene_id"]),
                gene_name=str(row["gene_name"]),
                transcript_name=row["transcript_name"],
                chromosome=str(row["chrom"]),
                start=int(row["start"]),
                end=int(row["end"]),
                strand=str(row["strand"]),
                transcript_type=row["transcript_type"],
                tags=tuple(json.loads(str(row["tags_json"]))),
                mane_status=row["mane_status"],
                mane_refseq_id=row["mane_refseq_id"],
                appris=row["appris"],
                is_canonical=bool(row["is_canonical"]),
                is_basic=bool(row["is_basic"]),
                cds_length=int(row["cds_length"]),
                transcript_length=int(row["transcript_length"]),
            )
            for row in rows
        ]
