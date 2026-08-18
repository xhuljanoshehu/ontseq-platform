from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from enum import StrEnum
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .models import StrictModel

_BND_ALT = re.compile(
    r"^(?P<prefix>[^\[\]]*)(?P<bracket>[\[\]])"
    r"(?P<chromosome>[^:\[\]]+):(?P<position>\d+)"
    r"(?P=bracket)(?P<suffix>[^\[\]]*)$"
)


class BreakendAltForm(StrEnum):
    """Privacy-safe representation of the four VCF breakend ALT forms.

    The names describe syntax only. They deliberately do not claim transcript direction,
    gene strand, or 5-prime/3-prime fusion orientation.
    """

    LOCAL_THEN_OPEN = "local_then_open"
    LOCAL_THEN_CLOSE = "local_then_close"
    OPEN_THEN_LOCAL = "open_then_local"
    CLOSE_THEN_LOCAL = "close_then_local"


class BreakendDescriptor(StrictModel):
    source_event_id: str = Field(pattern=r"^SNIFFLES2-\d{6}$")
    primary_chromosome: str
    primary_position_0based: int = Field(ge=0)
    mate_chromosome: str
    mate_position_0based: int = Field(ge=0)
    alt_form: BreakendAltForm
    inserted_sequence_retained: Literal[False] = False

    @model_validator(mode="after")
    def positions_are_single_breakpoints(self) -> BreakendDescriptor:
        if not self.primary_chromosome or not self.mate_chromosome:
            raise ValueError("breakend chromosomes must be non-empty")
        return self


def _open_vcf(path: Path) -> Iterator[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield from handle


def parse_breakend_alt(alternate: str) -> tuple[str, int, BreakendAltForm]:
    """Parse a VCF BND ALT allele without retaining local/inserted sequence.

    VCF encodes four adjacency forms using whether the remote locus appears before or after
    the local sequence and whether `[` or `]` brackets are used. This function preserves
    that four-state syntax and the mate locus only. It intentionally does not convert the
    syntax into biological 5-prime/3-prime direction.
    """

    match = _BND_ALT.fullmatch(alternate)
    if match is None:
        raise ValueError("unsupported or malformed VCF breakend ALT allele")
    prefix = match.group("prefix")
    suffix = match.group("suffix")
    bracket = match.group("bracket")
    if bool(prefix) == bool(suffix):
        raise ValueError("VCF breakend ALT must place local sequence on exactly one side")
    if prefix:
        alt_form = (
            BreakendAltForm.LOCAL_THEN_OPEN if bracket == "[" else BreakendAltForm.LOCAL_THEN_CLOSE
        )
    else:
        alt_form = (
            BreakendAltForm.OPEN_THEN_LOCAL if bracket == "[" else BreakendAltForm.CLOSE_THEN_LOCAL
        )
    return match.group("chromosome"), int(match.group("position")) - 1, alt_form


def breakend_descriptors_from_sniffles_vcf(path: Path) -> dict[str, BreakendDescriptor]:
    """Extract BND orientation descriptors keyed by normalized Sniffles event id.

    The normalized Sniffles adapter numbers non-header VCF records in encounter order.
    We mirror that numbering so the descriptor can be joined back to `GenomicEvent.event_id`.
    No VCF ID, read name, REF/ALT sequence, inserted sequence, or source path is retained.
    """

    if not path.is_file():
        raise ValueError("Sniffles VCF is missing or unreadable")
    descriptors: dict[str, BreakendDescriptor] = {}
    record_number = 0
    for raw_line in _open_vcf(path):
        line = raw_line.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        record_number += 1
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chromosome, raw_position, _record_id, _reference, alternate = fields[:5]
        info = fields[7]
        if "SVTYPE=BND" not in info and "SVTYPE=TRA" not in info:
            continue
        try:
            primary_position = int(raw_position) - 1
            mate_chromosome, mate_position, alt_form = parse_breakend_alt(alternate)
        except ValueError:
            continue
        event_id = f"SNIFFLES2-{record_number:06d}"
        descriptors[event_id] = BreakendDescriptor(
            source_event_id=event_id,
            primary_chromosome=chromosome,
            primary_position_0based=primary_position,
            mate_chromosome=mate_chromosome,
            mate_position_0based=mate_position,
            alt_form=alt_form,
        )
    return descriptors
