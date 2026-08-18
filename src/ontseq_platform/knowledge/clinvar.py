"""Reading a locked ClinVar release into records this pipeline can match against.

ClinVar is distributed as ``variant_summary.txt``, a tab-delimited export NCBI republishes
weekly. Three properties of that file drive everything here.

**It changes every week, so it must be locked.** An annotation carrying "ClinVar says
Pathogenic" without saying *which* ClinVar is not reproducible: the same command run a month
apart can produce different reports from the same BAM with nothing recording why. A release
is therefore identified by its own checksum, exactly as the reference genome and the cytoband
table are, and the identifier travels with every annotation derived from it.

**It carries both assemblies in one file.** Every row names its ``Assembly``, and rows for
the build you are not using are silently coordinate-incompatible rather than absent. Rows
whose assembly does not match are dropped and counted, because a quiet mismatch produces an
empty annotation set that looks exactly like a sample nobody knows anything about.

**Most of it is not what this pipeline produces.** ONTSeq emits copy-number and structural
findings; the bulk of ClinVar is single-nucleotide and small indel records with no useful
extent to match a segment against. Those rows are counted and skipped, so the report can say
how much of the release was applicable rather than implying the whole of it was consulted.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

from .scope import Interval, Origin, canonical_contig, review_stars

#: Column names in NCBI's ``variant_summary.txt``. Named rather than indexed: NCBI adds
#: columns between releases, and positional parsing would shift every field by one without
#: raising anything.
REQUIRED_COLUMNS: tuple[str, ...] = (
    "Type",
    "GeneSymbol",
    "ClinicalSignificance",
    "OriginSimple",
    "Assembly",
    "Chromosome",
    "Start",
    "Stop",
    "ReviewStatus",
    "PhenotypeList",
    "VariationID",
)

#: ClinVar ``Type`` values that describe an extent this pipeline can match a segment against.
#: Everything else — single nucleotide variants, small indels, microsatellites — has no
#: comparable extent, and matching a copy-number segment to one would assert a
#: correspondence that does not exist.
REGION_TYPES: frozenset[str] = frozenset(
    {
        "copy number gain",
        "copy number loss",
        "deletion",
        "duplication",
        "insertion",
        "inversion",
        "translocation",
        "complex",
        "fusion",
        "tandem duplication",
    }
)

#: ClinVar's assembly strings, as they appear in the file.
ASSEMBLIES: frozenset[str] = frozenset({"GRCh37", "GRCh38"})


class ClinVarError(ValueError):
    """Raised when a ClinVar release cannot be read as one."""


@dataclass(frozen=True)
class ClinVarRecord:
    """One ClinVar row, reduced to what can be matched and what must be shown."""

    variation_id: str
    record_type: str
    interval: Interval
    #: ClinVar's classification, verbatim and in ClinVar's own vocabulary.
    assertion: str
    #: ClinVar's own statement of what the record concerns.
    origin: Origin
    review_status: str
    #: NCBI's star rating for ``review_status``; ``None`` when the string is unrecognised.
    stars: int | None
    genes: tuple[str, ...]
    conditions: tuple[str, ...]

    @property
    def describe(self) -> str:
        stars = "?" if self.stars is None else str(self.stars)
        return (
            f"ClinVar {self.variation_id} ({self.record_type}): {self.assertion} "
            f"[{self.origin.value}, {stars}★]"
        )


@dataclass
class LoadSummary:
    """What the release contained and what was usable, so a reader can judge coverage."""

    rows_read: int = 0
    kept: int = 0
    wrong_assembly: int = 0
    not_a_region_type: int = 0
    unusable_coordinates: int = 0
    unknown_review_status: set[str] = field(default_factory=set)

    def describe(self) -> str:
        return (
            f"{self.rows_read} row(s) read, {self.kept} usable; "
            f"{self.wrong_assembly} for another assembly, "
            f"{self.not_a_region_type} without a matchable extent, "
            f"{self.unusable_coordinates} with unusable coordinates"
        )


def _origin(value: str) -> Origin:
    """Map ClinVar's ``OriginSimple`` onto the origin vocabulary.

    ClinVar uses ``germline``, ``somatic``, ``germline/somatic``, ``not provided`` and a few
    variants. Anything carrying both, or neither, resolves to ``UNKNOWN``: a record that
    covers both cannot settle which one this sample's finding is.
    """
    text = value.strip().lower()
    if text in {"germline", "inherited"}:
        return Origin.GERMLINE
    if text == "somatic":
        return Origin.SOMATIC
    return Origin.UNKNOWN


def release_sha256(path: Path, *, chunk_size: int = 1 << 20) -> str:
    """Checksum the release file, so an annotation can name the exact ClinVar it used."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def parse_rows(rows: Iterable[dict[str, str]], *, assembly: str) -> Iterator[ClinVarRecord]:
    """Yield the usable records from already-parsed rows. Separated for testability."""
    summary = LoadSummary()
    yield from _parse(rows, assembly=assembly, summary=summary)


def _parse(
    rows: Iterable[dict[str, str]], *, assembly: str, summary: LoadSummary
) -> Iterator[ClinVarRecord]:
    for row in rows:
        summary.rows_read += 1
        if row.get("Assembly", "").strip() != assembly:
            summary.wrong_assembly += 1
            continue
        record_type = row.get("Type", "").strip().lower()
        if record_type not in REGION_TYPES:
            summary.not_a_region_type += 1
            continue
        try:
            start = int(row["Start"])
            stop = int(row["Stop"])
            interval = Interval(
                contig=canonical_contig(row["Chromosome"]), start=start - 1, end=stop
            )
        except (KeyError, ValueError):
            # ClinVar uses -1 for records whose placement is unknown, and some rows have a
            # stop at or before the start. Neither can be matched, and both are counted
            # rather than dropped silently.
            summary.unusable_coordinates += 1
            continue

        review_status = row.get("ReviewStatus", "").strip()
        stars = review_stars(review_status)
        if stars is None and review_status:
            summary.unknown_review_status.add(review_status)

        summary.kept += 1
        yield ClinVarRecord(
            variation_id=row.get("VariationID", "").strip(),
            record_type=record_type,
            interval=interval,
            assertion=row.get("ClinicalSignificance", "").strip(),
            origin=_origin(row.get("OriginSimple", "")),
            review_status=review_status,
            stars=stars,
            genes=_split(row.get("GeneSymbol", "")),
            conditions=_split(row.get("PhenotypeList", "")),
        )


def _split(value: str) -> tuple[str, ...]:
    """ClinVar packs multiple values into one field with several separators."""
    items: list[str] = []
    for part in value.replace("|", ";").split(";"):
        text = part.strip()
        if text and text.lower() not in {"not provided", "not specified", "-"}:
            items.append(text)
    return tuple(dict.fromkeys(items))


def load(path: Path, *, assembly: str) -> tuple[list[ClinVarRecord], LoadSummary]:
    """Read a ClinVar release, keeping the records this pipeline can match against.

    Returns the records and a summary of what was skipped and why. The summary is not
    optional bookkeeping: without it a reader cannot tell an assay with nothing known about
    it from a release that was read with the wrong assembly.
    """
    if assembly not in ASSEMBLIES:
        raise ClinVarError(f"unknown assembly {assembly!r}; ClinVar publishes {sorted(ASSEMBLIES)}")
    if not path.is_file():
        raise ClinVarError(f"ClinVar release does not exist: {path}")

    summary = LoadSummary()
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fieldnames = [name.lstrip("#").strip() for name in (reader.fieldnames or [])]
        missing = [name for name in REQUIRED_COLUMNS if name not in fieldnames]
        if missing:
            raise ClinVarError(
                f"{path.name} is missing column(s) {', '.join(missing)}. Expected NCBI's "
                "variant_summary.txt; a different export will not carry the fields this "
                "annotation depends on."
            )
        reader.fieldnames = fieldnames
        records = list(_parse(reader, assembly=assembly, summary=summary))
    return records, summary
