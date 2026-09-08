"""Conservative read-level Dorado MM/ML and modkit extract-full adapters.

Supported scope is aligned simplex DNA, C+ m or C+ m/h tags, and reference CpGs.
Both reference strands remain separate; coordinates are the actual modified base
(0-based, inclusive point), so a minus-strand CpG is at the reference G. CIGAR,
reverse-complemented SAM SEQ, and original-orientation MM deltas are respected.
Opposite-strand MM groups, unknown cytosine modifications, paired reads, stale MN,
and duplicate primary alignments are rejected. Secondary/supplementary, QC-failed
and duplicate-flagged records are excluded. No implicit/omitted probability is
ever treated as a confident canonical call. High 5hmC is ambiguous for 5mC vs C.

SAM and modkit TSV need only the standard library. BAM additionally needs pysam;
CRAM is deliberately unsupported. A local uncompressed FASTA is mandatory for
both paths and is fingerprinted and indexed in memory without writing sidecars.
modkit extract *full* must preserve every called cytosine modification, with no
--ignore, probability redistribution or call-only filtering. bedMethyl cannot
preserve read groups and is not an accepted input.

Implementation references (consulted 2026-09-05):
https://samtools.github.io/hts-specs/SAMtags.pdf
https://nanoporetech.github.io/modkit/intro_extract.html
https://pysam.readthedocs.io/en/latest/api.html
https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib
import math
import re
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Any, BinaryIO, Literal, TextIO

from pydantic import ConfigDict, Field, model_validator

from .methylation_mixture import (
    MethylationMixturePolicy,
    NanopolishCallState,
    _canonical_chromosome,
    _MarkerKey,
    _NanopolishCall,
)
from .models import FileFingerprint, GenomeBuild, StrictModel
from .reference import sha256_file

_PROVENANCE = r"^[A-Za-z0-9][A-Za-z0-9 ._:+()\-]{0,199}$"
_MODEL_PROVENANCE = r"^[A-Za-z0-9][A-Za-z0-9 ._:@+(),\-]{0,799}$"
MODBAM_ADAPTER_VERSION = "0.1.1"
# Reviewed release source, not a claim of executed-binary validation:
# nanoporetech/modkit 481e3c9e7930f3f499eadf1ef441606f33e6881c (v0.6.1).
# v0.4.1 extract-full lacks the mandatory alignment_start/alignment_end columns.
SUPPORTED_MODKIT_EXTRACT_VERSIONS = ("0.6.1",)
_IDENTIFIER = r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$"
_CG_HASH = hashlib.sha256(b"CG").hexdigest()
_EXCLUDE_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
_MM_GROUP = re.compile(r"([ACGTUN])([+-])([a-z]+|[0-9]+)([.?]?)(,\d+(?:,\d+)*)?;")
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")


class ModbamAdapterPolicy(StrictModel):
    """Technical defaults, never validated assay thresholds."""

    model_config = ConfigDict(extra="forbid", frozen=True)
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(default="dorado-explicit-cpg-v1", pattern=_IDENTIFIER)
    status: Literal["technical_defaults_only"] = "technical_defaults_only"
    minimum_class_probability: float = Field(default=0.8, gt=0.5, le=1)
    implicit_call_handling: Literal["ambiguous"] = "ambiguous"
    other_modification_handling: Literal["ambiguous"] = "ambiguous"
    maximum_bases_per_read: int = Field(default=2_000_000, ge=1)
    maximum_reference_bytes: int = Field(default=5_000_000_000, ge=1)


class ModbamSourceMetadata(StrictModel):
    """Declared upstream provenance; a declaration is not biological validation."""

    model_config = ConfigDict(extra="forbid", frozen=True)
    schema_version: Literal["0.1.0"] = "0.1.0"
    caller_name: Literal["dorado"] = "dorado"
    caller_version: str = Field(pattern=_PROVENANCE)
    reference_id: str = Field(pattern=_PROVENANCE)
    reference_genome_build: GenomeBuild
    reference_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    sequencing_platform: str = Field(pattern=_PROVENANCE)
    flow_cell_product_code: str = Field(pattern=_PROVENANCE)
    library_kit: str = Field(pattern=_PROVENANCE)
    basecaller_name: Literal["dorado"] = "dorado"
    basecaller_version: str = Field(pattern=_PROVENANCE)
    basecaller_model: str = Field(pattern=_MODEL_PROVENANCE)
    modification_model: str = Field(pattern=_MODEL_PROVENANCE)
    modification_context: Literal["CpG_5mC_vs_unmodified_C"] = "CpG_5mC_vs_unmodified_C"
    cytosine_modifications: Literal["5mC", "5mC_5hmC"]
    modkit_version: Literal["0.6.1"] | None = None
    source_modbam_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    probability_transform: Literal["none"] = "none"

    def unknown_fields(self) -> tuple[str, ...]:
        return tuple(
            name
            for name, value in self.model_dump().items()
            if value == "unknown" and name != "schema_version"
        )

    @property
    def modification_codes(self) -> frozenset[str]:
        return frozenset({"m", "h"} if self.cytosine_modifications == "5mC_5hmC" else {"m"})


MODBAM_COMPATIBILITY_FIELDS = (
    "caller_name",
    "caller_version",
    "reference_id",
    "reference_genome_build",
    "reference_sha256",
    "sequencing_platform",
    "flow_cell_product_code",
    "library_kit",
    "basecaller_name",
    "basecaller_version",
    "basecaller_model",
    "modification_model",
    "modification_context",
    "cytosine_modifications",
    "modkit_version",
    "probability_transform",
)


class ModbamSourceSummary(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    source_id: str = Field(pattern=_IDENTIFIER)
    genome_build: GenomeBuild
    input_format: Literal["modbam_mm_ml", "sam_mm_ml", "modkit_extract_full_tsv"]
    coordinate_system: Literal["reference_modified_base_0based_point_stranded"] = (
        "reference_modified_base_0based_point_stranded"
    )
    chromosome_normalization: Literal["canonical_chr_prefix"] = "canonical_chr_prefix"
    upstream: ModbamSourceMetadata
    adapter_policy: ModbamAdapterPolicy
    fingerprint: FileFingerprint
    reference_fingerprint: FileFingerprint
    parser_name: str
    parser_version: str
    rows_total: int = Field(ge=0)
    canonical_call_rows: int = Field(ge=0)
    informative_call_rows: int = Field(ge=0)
    ambiguous_call_rows: int = Field(ge=0)
    read_group_count: int = Field(ge=1)
    excluded_counts: dict[str, int] = Field(default_factory=dict)
    header_verified_fields: list[str] = Field(default_factory=list)
    provenance_verification: Literal["declared_with_available_header_cross_checks"] = (
        "declared_with_available_header_cross_checks"
    )

    @model_validator(mode="after")
    def consistent(self) -> ModbamSourceSummary:
        if self.canonical_call_rows != self.informative_call_rows + self.ambiguous_call_rows:
            raise ValueError("Adapter call counts are inconsistent")
        if self.upstream.reference_genome_build != self.genome_build:
            raise ValueError("Declared reference genome build contradicts the source")
        if self.reference_fingerprint.sha256 != self.upstream.reference_sha256:
            raise ValueError("Adapter reference fingerprint contradicts declared provenance")
        if self.fingerprint.sha256 is None or any(v < 0 for v in self.excluded_counts.values()):
            raise ValueError("Adapter provenance or exclusion counters are invalid")
        return self


@dataclass(frozen=True)
class ParsedModbamSource:
    summary: ModbamSourceSummary
    calls_by_read: Mapping[str, tuple[_NanopolishCall, ...]]


@dataclass(frozen=True)
class _FastaEntry:
    offset: int
    length: int
    line_bases: int
    line_bytes: int


class _Reference:
    """Read-only index for regular-line, uncompressed DNA FASTA; no sidecar trust."""

    def __init__(self, handle: BinaryIO) -> None:
        self.handle = handle
        self.entries: dict[str, _FastaEntry] = {}
        name: str | None = None
        offset = length = line_bases = line_bytes = 0
        short_line = False
        normalized_names: set[str] = set()
        while raw := handle.readline(1_000_001):
            if len(raw) > 1_000_000:
                raise ValueError("Reference FASTA line exceeds the indexing safety limit")
            if raw.startswith(b">"):
                if name is not None:
                    if length == 0:
                        raise ValueError("Reference contains an empty sequence")
                    self.entries[name] = _FastaEntry(offset, length, line_bases, line_bytes)
                try:
                    name = raw[1:].split()[0].decode("ascii")
                except (IndexError, UnicodeError) as exc:
                    raise ValueError("Reference has an invalid FASTA header") from exc
                if name in self.entries:
                    raise ValueError("Reference contains duplicate sequence names")
                normalized = _canonical_chromosome(name)
                if normalized is not None:
                    if normalized in normalized_names:
                        raise ValueError("Reference has ambiguous chromosome aliases")
                    normalized_names.add(normalized)
                offset, length, line_bases, line_bytes = handle.tell(), 0, 0, 0
                short_line = False
                continue
            bases = raw.rstrip(b"\r\n")
            if (
                name is None
                or not bases
                or re.fullmatch(b"[ACGTRYKMSWBDHVNacgtrykmswbdhvn]+", bases) is None
            ):
                raise ValueError("Reference must be an uncompressed DNA FASTA")
            if line_bases == 0:
                line_bases, line_bytes = len(bases), len(raw)
            elif short_line or len(bases) > line_bases:
                raise ValueError("Reference FASTA must use regular sequence line lengths")
            elif len(bases) < line_bases or len(raw) != line_bytes:
                short_line = True
            length += len(bases)
        if name is None or length == 0:
            raise ValueError("Reference contains no sequence")
        self.entries[name] = _FastaEntry(offset, length, line_bases, line_bytes)

    def base(self, chromosome: str, position: int) -> str:
        entry = self.entries.get(chromosome)
        if entry is None:
            raise ValueError("Aligned chromosome is absent from the locked reference")
        if position < 0 or position >= entry.length:
            return "N"
        offset = entry.offset + (position // entry.line_bases) * entry.line_bytes
        self.handle.seek(offset + position % entry.line_bases)
        return self.handle.read(1).decode("ascii").upper()

    def is_cpg(self, chromosome: str, position: int, strand: Literal["+", "-"]) -> bool:
        start = position if strand == "+" else position - 1
        return self.base(chromosome, start) + self.base(chromosome, start + 1) == "CG"


def _fingerprint(path: Path, maximum: int) -> FileFingerprint:
    if not path.is_file():
        raise ValueError("Adapter input or reference is missing or unreadable")
    size = path.stat().st_size
    if size > maximum:
        raise ValueError("Adapter input exceeds its locked byte safety limit")
    return FileFingerprint(size_bytes=size, sha256=sha256_file(path))


def _unchanged(path: Path, fingerprint: FileFingerprint) -> None:
    if path.stat().st_size != fingerprint.size_bytes or sha256_file(path) != fingerprint.sha256:
        raise ValueError("Adapter input or reference changed while it was being parsed")


@contextmanager
def _open_text(path: Path) -> Iterator[TextIO]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            yield handle


def _bounded_lines(handle: TextIO, maximum_bytes: int) -> Iterator[str]:
    total = 0
    while line := handle.readline(maximum_bytes + 1):
        total += len(line.encode("utf-8"))
        if total > maximum_bytes:
            raise ValueError("Adapter decompressed input exceeds its locked byte safety limit")
        yield line


def _integer(value: str, label: str, minimum: int = 0) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"Adapter has an invalid {label}") from exc
    if parsed < minimum:
        raise ValueError(f"Adapter has an invalid {label}")
    return parsed


@dataclass
class _Collector:
    policy: MethylationMixturePolicy
    calls: dict[str, list[_NanopolishCall]] = field(default_factory=dict)
    seen: set[tuple[str, str, str, int]] = field(default_factory=set)
    excluded: dict[str, int] = field(default_factory=dict)
    rows_total: int = 0
    call_count: int = 0
    ambiguous_count: int = 0
    header_verified_fields: set[str] = field(default_factory=set)

    def row(self) -> None:
        self.rows_total += 1
        if self.rows_total > self.policy.maximum_rows_per_source:
            raise ValueError("Adapter input exceeds its locked row safety limit")

    def exclude(self, reason: str) -> None:
        self.excluded[reason] = self.excluded.get(reason, 0) + 1

    def add(
        self,
        read_id: str,
        chromosome: str,
        position: int,
        strand: Literal["+", "-"],
        state: NanopolishCallState,
    ) -> None:
        key = (read_id, chromosome, strand, position)
        if key in self.seen:
            raise ValueError("Adapter duplicates a CpG within one read group")
        self.seen.add(key)
        self.call_count += 1
        if self.call_count > self.policy.maximum_rows_per_source:
            raise ValueError("Adapter calls exceed the locked in-memory row safety limit")
        self.ambiguous_count += state == NanopolishCallState.AMBIGUOUS
        marker = _MarkerKey(chromosome, strand, position, position, 1, _CG_HASH)
        self.calls.setdefault(read_id, []).append(_NanopolishCall(marker, state))


def _classify(
    probabilities: Mapping[str, tuple[float, float] | None],
    expected_codes: frozenset[str],
    threshold: float,
) -> NanopolishCallState:
    if set(probabilities) != expected_codes or any(p is None for p in probabilities.values()):
        return NanopolishCallState.AMBIGUOUS
    explicit = {code: pair for code, pair in probabilities.items() if pair is not None}
    if sum(pair[0] for pair in explicit.values()) > 1 + 1e-7:
        raise ValueError("Modification probabilities exceed total probability one")
    if explicit["m"][0] >= threshold:
        return NanopolishCallState.METHYLATED
    # Lower bound on canonical probability includes *every* called modification.
    canonical_lower = max(0.0, 1.0 - sum(pair[1] for pair in explicit.values()))
    if canonical_lower >= threshold:
        return NanopolishCallState.UNMETHYLATED
    return NanopolishCallState.AMBIGUOUS


@dataclass(frozen=True)
class ModbamAlignment:
    """Minimal typed SAM-orientation record used by the common MM/ML decoder."""

    read_id: str
    chromosome: str
    flag: int
    reference_start: int
    sequence: str
    cigar: str
    mm: str | None
    ml: tuple[int, ...] | None
    mn: int | None = None
    dx: int | None = None


def _alignment_map(record: ModbamAlignment, maximum_bases: int) -> list[int | None]:
    sequence = record.sequence
    if not sequence or sequence == "*" or len(sequence) > maximum_bases:
        raise ValueError("Adapter sequence is missing or exceeds the per-read safety limit")
    if re.fullmatch("[ACGTRYKMSWBDHVN]+", sequence) is None:
        raise ValueError("Adapter requires a DNA query sequence")
    if record.reference_start < 0:
        raise ValueError("Mapped alignment has an invalid reference position")
    operations = [(int(m[1]), m[2]) for m in _CIGAR.finditer(record.cigar)]
    if "".join(f"{length}{op}" for length, op in operations) != record.cigar:
        raise ValueError("Adapter has an invalid CIGAR")
    if not operations or any(length < 1 for length, _ in operations):
        raise ValueError("Adapter has an invalid CIGAR")
    if sum(length for length, op in operations if op in {"M", "I", "S", "=", "X"}) != len(sequence):
        raise ValueError("CIGAR query length disagrees with sequence")
    if record.mn is not None and record.mn != len(sequence):
        raise ValueError("MN tag disagrees with sequence length; MM/ML may be stale")
    if any(op == "H" for _, op in operations) and record.mn is None:
        raise ValueError("Hard-clipped MM/ML records require a current MN tag")
    # Reject internal clipping, allowing the standard terminal H/S combinations.
    unclip = list(operations)
    for side in (0, -1):
        if unclip and unclip[side][1] == "H":
            unclip.pop(side)
        if unclip and unclip[side][1] == "S":
            unclip.pop(side)
    if any(op in {"S", "H"} for _, op in unclip):
        raise ValueError("Adapter does not support internal CIGAR clipping")
    mapping: list[int | None] = []
    reference_position = record.reference_start
    for length, op in operations:
        if op in {"M", "=", "X"}:
            mapping.extend(range(reference_position, reference_position + length))
            reference_position += length
        elif op in {"D", "N"}:
            reference_position += length
        elif op in {"I", "S"}:
            mapping.extend([None] * length)
        elif op == "P":
            raise ValueError("Padded CIGAR is outside the adapter scope")
        if len(mapping) > len(sequence):
            raise ValueError("CIGAR query length disagrees with sequence")
    if len(mapping) != len(sequence):
        raise ValueError("CIGAR query length disagrees with sequence")
    return mapping


def decode_mm_ml(
    record: ModbamAlignment,
    *,
    expected_codes: frozenset[str],
    maximum_bases: int = 2_000_000,
) -> tuple[str, list[int | None], Mapping[int, Mapping[str, tuple[float, float] | None]]]:
    """Decode explicit probability intervals in original read orientation.

    MM '.'/'?' and absent skip flags all leave omitted probabilities unknown to
    this adapter; a caller's unspecified elision threshold cannot pass our gate.
    Other fundamental-base groups are validated/consumed but not returned.
    """
    if expected_codes not in {frozenset({"m"}), frozenset({"m", "h"})}:
        raise ValueError("Adapter requires the declared 5mC or 5mC/5hmC model")
    if record.dx is not None:
        if type(record.dx) is not int or record.dx not in {-1, 0, 1}:
            raise ValueError("Dorado dx must be one of -1, 0 or 1")
        if record.dx == 1:
            raise ValueError("Duplex reads are outside the simplex ONT adapter scope")
        # dx=-1 is a simplex parent, not its duplex offspring. Its presence
        # does not establish molecule independence for the retained read groups.
    mapping = _alignment_map(record, maximum_bases)
    reverse = bool(record.flag & 0x10)
    forward = record.sequence.translate(_COMPLEMENT)[::-1] if reverse else record.sequence
    if reverse:
        mapping.reverse()
    if record.mm is None or record.ml is None:
        raise ValueError("Aligned modBAM record lacks required MM/ML tags")
    if any(type(value) is not int or value < 0 or value > 255 for value in record.ml):
        raise ValueError("ML must contain unsigned byte probabilities")
    groups = list(_MM_GROUP.finditer(record.mm))
    if not groups or "".join(match[0] for match in groups) != record.mm:
        raise ValueError("MM tag does not follow the supported SAM grammar")
    events: dict[int, dict[str, tuple[float, float] | None]] = {}
    declared: set[tuple[str, str, str]] = set()
    cursor = 0
    for group in groups:
        base, strand, raw_codes, _skip_flag, deltas = group.groups()
        codes = tuple(raw_codes) if raw_codes.isalpha() else (raw_codes,)
        if strand != "+":
            raise ValueError("Opposite-strand MM groups are outside the Dorado adapter scope")
        if base in {"N", "U"}:
            raise ValueError("Wildcard/RNA MM groups are outside the DNA Dorado adapter scope")
        if base == "C" and set(codes) - expected_codes:
            raise ValueError("MM cytosine modifications contradict the declared model")
        for code in codes:
            signature = (base, strand, code)
            if signature in declared:
                raise ValueError("MM tag repeats a modification group")
            declared.add(signature)
        positions = [index for index, value in enumerate(forward) if value == base or base == "N"]
        index = -1
        for delta in deltas[1:].split(",") if deltas else ():
            index += int(delta) + 1
            if index >= len(positions):
                raise ValueError("MM delta extends beyond the query sequence")
            position = positions[index]
            for code in codes:
                if cursor >= len(record.ml):
                    raise ValueError("ML count does not match MM events")
                probability = record.ml[cursor]
                cursor += 1
                if base == "C":
                    events.setdefault(position, {})[code] = (
                        probability / 256.0,
                        (probability + 1) / 256.0,
                    )
    if cursor != len(record.ml):
        raise ValueError("ML count does not match MM events")
    if {code for base, _, code in declared if base == "C"} != expected_codes:
        raise ValueError("MM cytosine modification groups do not match the declared model")
    return forward, mapping, events


def _parse_sam_record(line: str) -> ModbamAlignment:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 11:
        raise ValueError("SAM alignment has fewer than eleven fields")
    flag = _integer(fields[1], "SAM flag")
    if flag > 65535:
        raise ValueError("SAM flag exceeds the unsigned 16-bit range")
    tags: dict[str, tuple[str, str]] = {}
    for value in fields[11:]:
        parts = value.split(":", 2)
        if len(parts) != 3:
            raise ValueError("SAM alignment has an invalid auxiliary tag")
        tag, dtype, raw = parts
        if tag in {"MM", "ML", "MN", "dx"}:
            if tag in tags:
                raise ValueError("SAM alignment repeats MM, ML, MN or dx")
            tags[tag] = (dtype, raw)
        elif tag in {"Mm", "Ml"}:
            raise ValueError("Legacy Mm/Ml tags require upstream normalization")
    mm: str | None = None
    ml: tuple[int, ...] | None = None
    mn: int | None = None
    dx: int | None = None
    if "MM" in tags:
        if tags["MM"][0] != "Z":
            raise ValueError("MM must have SAM type Z")
        mm = tags["MM"][1]
    if "ML" in tags:
        dtype, raw = tags["ML"]
        if dtype != "B" or (raw != "C" and not raw.startswith("C,")):
            raise ValueError("ML must have SAM unsigned byte-array type B:C")
        ml = tuple(_integer(value, "ML byte") for value in raw[2:].split(",")) if raw != "C" else ()
    if "MN" in tags:
        if tags["MN"][0] != "i":
            raise ValueError("MN must have SAM integer type i")
        mn = _integer(tags["MN"][1], "MN")
    if "dx" in tags:
        if tags["dx"][0] != "i":
            raise ValueError("Dorado dx must have SAM integer type i")
        dx = _integer(tags["dx"][1], "Dorado dx", minimum=-1)
    return ModbamAlignment(
        fields[0],
        fields[2],
        flag,
        _integer(fields[3], "SAM position") - 1,
        fields[9].upper(),
        fields[5],
        mm,
        ml,
        mn,
        dx,
    )


def _validate_alignment_header(
    lines: Sequence[str],
    reference: _Reference,
    upstream: ModbamSourceMetadata,
) -> tuple[set[str], set[str]]:
    """Check dictionaries and documented structured Dorado metadata when present.

    Raw headers can contain identifiers and command-line paths. Never put their
    contents in error messages or the returned field-name-only verification list.
    Header absence is not evidence that caller/model declarations were verified.
    """
    sequences: set[str] = set()
    read_groups: set[str] = set()
    verified: set[str] = set()
    for line in lines:
        fields = line.rstrip("\r\n").split("\t")
        if fields[0] not in {"@SQ", "@PG", "@RG"}:
            continue
        tags: dict[str, str] = {}
        for text in fields[1:]:
            parts = text.split(":", 1)
            if len(parts) != 2 or parts[0] in tags:
                raise ValueError("Alignment header contains malformed or duplicate fields")
            tags[parts[0]] = parts[1]
        if fields[0] == "@SQ":
            name = tags.get("SN", "")
            length = _integer(tags.get("LN", ""), "@SQ length", minimum=1)
            if name in sequences:
                raise ValueError("Alignment @SQ dictionary repeats a sequence")
            sequences.add(name)
            entry = reference.entries.get(name)
            if entry is None or entry.length != length:
                raise ValueError("Alignment @SQ dictionary disagrees with the locked reference")
        elif fields[0] == "@PG":
            role = tags.get("ID", "").split(".", 1)[0].lower()
            caller = tags.get("PN", "").lower()
            if role == "duplex":
                raise ValueError("Duplex basecalling is outside the simplex ONT adapter scope")
            if role in {"basecaller", "dorado"}:
                if caller and caller != "dorado":
                    raise ValueError("Alignment @PG caller contradicts declared Dorado provenance")
                version = tags.get("VN")
                if caller == "dorado" and version is not None:
                    if version != upstream.caller_version or version != upstream.basecaller_version:
                        raise ValueError(
                            "Alignment @PG Dorado version contradicts declared provenance"
                        )
                    verified.update({"caller_version", "basecaller_version"})
        else:
            name = tags.get("ID", "")
            if not name or name in read_groups:
                raise ValueError("Alignment @RG header is missing or repeats its ID")
            read_groups.add(name)
            if tags.get("PL", "ONT") != "ONT":
                raise ValueError("Alignment @RG platform is not ONT")
            description: dict[str, str] = {}
            for token in tags.get("DS", "").split():
                if "=" not in token:
                    continue
                key, value = token.split("=", 1)
                if key in description:
                    raise ValueError("Alignment @RG DS repeats a provenance field")
                description[key] = value
            base_model = description.get("basecall_model")
            if base_model is not None:
                if base_model != upstream.basecaller_model:
                    raise ValueError("Alignment @RG basecall model contradicts declared provenance")
                verified.add("basecaller_model")
            mod_models = description.get("modbase_models")
            if mod_models is not None:
                if set(mod_models.split(",")) != set(upstream.modification_model.split(",")):
                    raise ValueError(
                        "Alignment @RG modification model contradicts declared provenance"
                    )
                verified.add("modification_model")
    if not sequences:
        raise ValueError("Aligned SAM/BAM requires an @SQ reference dictionary")
    verified.add("reference_sequence_dictionary")
    return verified, sequences


def _sam_records(
    path: Path,
    maximum_bytes: int,
    reference: _Reference,
    upstream: ModbamSourceMetadata,
    collector: _Collector,
) -> Iterator[ModbamAlignment]:
    with _open_text(path) as handle:
        headers: list[str] = []
        checked_header = False
        sequences: set[str] = set()
        for line in _bounded_lines(handle, maximum_bytes):
            if line.startswith("@"):
                if checked_header:
                    raise ValueError("SAM header appears after an alignment record")
                headers.append(line)
                continue
            if not checked_header:
                verified, sequences = _validate_alignment_header(headers, reference, upstream)
                collector.header_verified_fields.update(verified)
                checked_header = True
            record = _parse_sam_record(line)
            if not record.flag & _EXCLUDE_FLAGS and record.chromosome not in sequences:
                raise ValueError("SAM alignment chromosome is absent from its @SQ dictionary")
            yield record
        if not checked_header:
            verified, _ = _validate_alignment_header(headers, reference, upstream)
            collector.header_verified_fields.update(verified)


def _bam_records(
    path: Path,
    pysam: Any,
    reference: _Reference,
    upstream: ModbamSourceMetadata,
    collector: _Collector,
) -> Iterator[ModbamAlignment]:
    with pysam.AlignmentFile(str(path), "rb", check_sq=True) as bam:
        verified, _ = _validate_alignment_header(str(bam.header).splitlines(), reference, upstream)
        collector.header_verified_fields.update(verified)
        for record in bam.fetch(until_eof=True):
            # Reuse strict SAM tag type validation. This is serialization of the
            # typed alignment, not parsing a human-facing presentation/report.
            yield _parse_sam_record(record.to_string())


def _validate_source(
    path: Path,
    reference_path: Path,
    source_id: str,
    genome_build: GenomeBuild,
    policy: MethylationMixturePolicy,
    upstream: ModbamSourceMetadata,
    adapter_policy: ModbamAdapterPolicy,
) -> tuple[FileFingerprint, FileFingerprint]:
    if re.fullmatch(_IDENTIFIER, source_id) is None:
        raise ValueError("source_id must be a safe pseudonymous identifier")
    if upstream.reference_genome_build != genome_build:
        raise ValueError("Declared genome build contradicts the requested source build")
    fingerprint = _fingerprint(path, policy.maximum_input_bytes_per_source)
    reference = _fingerprint(reference_path, adapter_policy.maximum_reference_bytes)
    if reference.sha256 != upstream.reference_sha256:
        raise ValueError("Reference SHA-256 differs from the locked upstream metadata")
    return fingerprint, reference


def _finish(
    collector: _Collector,
    *,
    path: Path,
    reference_path: Path,
    source_id: str,
    genome_build: GenomeBuild,
    upstream: ModbamSourceMetadata,
    adapter_policy: ModbamAdapterPolicy,
    fingerprint: FileFingerprint,
    reference: FileFingerprint,
    input_format: Literal["modbam_mm_ml", "sam_mm_ml", "modkit_extract_full_tsv"],
    parser_name: str,
    parser_version: str,
) -> ParsedModbamSource:
    _unchanged(path, fingerprint)
    _unchanged(reference_path, reference)
    if not collector.calls:
        raise ValueError("Adapter source contains no retained reference-CpG read groups")
    summary = ModbamSourceSummary(
        source_id=source_id,
        genome_build=genome_build,
        input_format=input_format,
        upstream=upstream,
        adapter_policy=adapter_policy,
        fingerprint=fingerprint,
        reference_fingerprint=reference,
        parser_name=parser_name,
        parser_version=parser_version,
        rows_total=collector.rows_total,
        canonical_call_rows=collector.call_count,
        informative_call_rows=collector.call_count - collector.ambiguous_count,
        ambiguous_call_rows=collector.ambiguous_count,
        read_group_count=len(collector.calls),
        excluded_counts=collector.excluded,
        header_verified_fields=sorted(collector.header_verified_fields),
    )
    return ParsedModbamSource(
        summary, {name: tuple(calls) for name, calls in collector.calls.items()}
    )


def parse_modbam_source(
    path: Path,
    *,
    reference_path: Path,
    source_id: str,
    genome_build: GenomeBuild,
    policy: MethylationMixturePolicy,
    upstream: ModbamSourceMetadata,
    adapter_policy: ModbamAdapterPolicy | None = None,
) -> ParsedModbamSource:
    """Read local SAM[.gz] or BAM; never load an alignment from a remote URI."""
    adapter_policy = adapter_policy or ModbamAdapterPolicy()
    fingerprint, ref_fingerprint = _validate_source(
        path,
        reference_path,
        source_id,
        genome_build,
        policy,
        upstream,
        adapter_policy,
    )
    parser_name, parser_version = "ontseq_mm_ml", MODBAM_ADAPTER_VERSION
    input_format: Literal["sam_mm_ml", "modbam_mm_ml"] = "sam_mm_ml"
    collector = _Collector(policy)
    records_factory: Callable[[_Reference], Iterator[ModbamAlignment]]
    if path.suffix.lower() == ".bam":
        try:
            pysam = importlib.import_module("pysam")
        except ImportError as exc:
            raise ValueError(
                "BAM MM/ML requires optional pysam; use local SAM or modkit TSV"
            ) from exc
        records_factory = partial(_bam_records, path, pysam, upstream=upstream, collector=collector)
        parser_name = "pysam/ontseq_mm_ml"
        parser_version = f"{pysam.__version__}/{MODBAM_ADAPTER_VERSION}"
        input_format = "modbam_mm_ml"
    elif path.name.lower().endswith((".sam", ".sam.gz")):
        records_factory = partial(
            _sam_records,
            path,
            policy.maximum_input_bytes_per_source,
            upstream=upstream,
            collector=collector,
        )
    else:
        raise ValueError("MM/ML adapter accepts only SAM, SAM.gz or BAM; CRAM is unsupported")
    primary_reads: set[str] = set()
    with reference_path.open("rb") as handle:
        reference = _Reference(handle)
        for record in records_factory(reference):
            collector.row()
            if record.flag & _EXCLUDE_FLAGS:
                collector.exclude("alignment_flags")
                continue
            if record.flag & 1:
                raise ValueError("Paired alignments are outside the simplex ONT adapter scope")
            if not record.read_id or record.read_id == "*":
                raise ValueError("Adapter requires a read identifier")
            if record.read_id in primary_reads:
                raise ValueError("Adapter found duplicate primary alignments for one read group")
            primary_reads.add(record.read_id)
            chromosome = _canonical_chromosome(record.chromosome)
            if chromosome is None:
                collector.exclude("noncanonical_chromosome")
                continue
            forward, mapping, events = decode_mm_ml(
                record,
                expected_codes=upstream.modification_codes,
                maximum_bases=adapter_policy.maximum_bases_per_read,
            )
            reference_span = sum(
                int(match[1])
                for match in _CIGAR.finditer(record.cigar)
                if match[2] in {"M", "D", "N", "=", "X"}
            )
            if reference_span < 1 or (
                record.reference_start + reference_span
                > reference.entries[record.chromosome].length
            ):
                raise ValueError("Alignment CIGAR extends beyond the locked reference")
            strand: Literal["+", "-"] = "-" if record.flag & 0x10 else "+"
            for query_position, base in enumerate(forward):
                if base != "C":
                    continue
                reference_position = mapping[query_position]
                if reference_position is None:
                    collector.exclude("unaligned_cytosine")
                    continue
                if not reference.is_cpg(record.chromosome, reference_position, strand):
                    collector.exclude("non_cpg_or_reference_mismatch")
                    continue
                state = _classify(
                    events.get(query_position, {}),
                    upstream.modification_codes,
                    adapter_policy.minimum_class_probability,
                )
                collector.add(record.read_id, chromosome, reference_position, strand, state)
    return _finish(
        collector,
        path=path,
        reference_path=reference_path,
        source_id=source_id,
        genome_build=genome_build,
        upstream=upstream,
        adapter_policy=adapter_policy,
        fingerprint=fingerprint,
        reference=ref_fingerprint,
        input_format=input_format,
        parser_name=parser_name,
        parser_version=parser_version,
    )


_MODKIT_REQUIRED = frozenset(
    {
        "read_id",
        "forward_read_position",
        "ref_position",
        "chrom",
        "mod_strand",
        "ref_strand",
        "ref_mod_strand",
        "read_length",
        "mod_qual",
        "mod_code",
        "canonical_base",
        "modified_primary_base",
        "inferred",
        "flag",
        "alignment_start",
        "alignment_end",
    }
)


@dataclass
class _ModkitSite:
    read_id: str
    chromosome: str
    position: int
    strand: Literal["+", "-"]
    query_position: int
    probabilities: dict[str, tuple[float, float] | None] = field(default_factory=dict)


def _reviewed_modkit_version(version: str) -> Literal["0.6.1"]:
    if version not in SUPPORTED_MODKIT_EXTRACT_VERSIONS:
        raise ValueError(
            "Read-level modkit extraction supports only source-reviewed version 0.6.1; "
            "other releases require a separate schema/conformance review"
        )
    return "0.6.1"


def _modkit_probability_interval(probability: float) -> tuple[float, float]:
    """Recover the ML bin from v0.6.1's rounded float midpoint representation."""
    byte = round(probability * 256 - 0.5)
    if not 0 <= byte <= 255 or abs(probability - (byte + 0.5) / 256) > 1e-7:
        raise ValueError(
            "Explicit modkit probability is not an untransformed 0.6.1 ML-bin midpoint"
        )
    return byte / 256, (byte + 1) / 256


def parse_modkit_source(
    path: Path,
    *,
    reference_path: Path,
    source_id: str,
    genome_build: GenomeBuild,
    policy: MethylationMixturePolicy,
    upstream: ModbamSourceMetadata,
    adapter_policy: ModbamAdapterPolicy | None = None,
) -> ParsedModbamSource:
    """Read modkit extract-full probabilities, joining m/h by read and position.

    All modification rows must be retained. Explicit probabilities are interpreted
    as rounded ML-bin midpoints from the source-reviewed modkit 0.6.1 schema and
    reconstructed as the same intervals as SAM/BAM. The TSV itself cannot prove
    omitted-class absence, upstream flags or
    completeness; its source BAM SHA and untransformed provenance are mandatory.
    """
    adapter_policy = adapter_policy or ModbamAdapterPolicy()
    fingerprint, ref_fingerprint = _validate_source(
        path,
        reference_path,
        source_id,
        genome_build,
        policy,
        upstream,
        adapter_policy,
    )
    if upstream.modkit_version is None or upstream.source_modbam_sha256 is None:
        raise ValueError("modkit TSV requires a modkit version and source modBAM SHA-256")
    _reviewed_modkit_version(upstream.modkit_version)
    collector = _Collector(policy)
    sites: dict[tuple[str, str, int, str], _ModkitSite] = {}
    read_alignments: dict[str, tuple[str, str, int, int, int, int]] = {}
    query_sites: dict[tuple[str, int], tuple[str, int, str]] = {}
    with reference_path.open("rb") as ref_handle, _open_text(path) as handle:
        reference = _Reference(ref_handle)
        reader = csv.DictReader(
            _bounded_lines(handle, policy.maximum_input_bytes_per_source),
            delimiter="\t",
        )
        names = reader.fieldnames
        if names is None or len(names) != len(set(names)) or not set(names) >= _MODKIT_REQUIRED:
            raise ValueError(
                "Input must be modkit extract full with named probability/read columns"
            )
        for row in reader:
            collector.row()
            if None in row or any(value is None for value in row.values()):
                raise ValueError("modkit row width contradicts its header")
            flag = _integer(row["flag"], "modkit flag")
            if flag > 65535:
                raise ValueError("modkit flag exceeds unsigned 16-bit range")
            if flag & _EXCLUDE_FLAGS:
                collector.exclude("alignment_flags")
                continue
            if flag & 1:
                raise ValueError("Paired alignments are outside the simplex ONT adapter scope")
            position = _integer(row["ref_position"], "reference position", minimum=-1)
            if position == -1:
                collector.exclude("unaligned_base")
                continue
            chromosome = _canonical_chromosome(row["chrom"])
            if chromosome is None:
                collector.exclude("noncanonical_chromosome")
                continue
            if row["mod_strand"] != "+":
                raise ValueError("Opposite-strand modifications are outside the Dorado scope")
            strand: Literal["+", "-"] = "-" if flag & 0x10 else "+"
            if row["ref_strand"] != strand or row["ref_mod_strand"] != strand:
                raise ValueError("modkit strand columns contradict the alignment flag")
            read_id = row["read_id"].strip()
            if not read_id or read_id == "*":
                raise ValueError("modkit requires a read identifier")
            read_length = _integer(row["read_length"], "read length", minimum=1)
            query_position = _integer(row["forward_read_position"], "query position")
            if read_length > adapter_policy.maximum_bases_per_read or query_position >= read_length:
                raise ValueError(
                    "modkit query position or read length is outside the allowed range"
                )
            alignment_start = _integer(row["alignment_start"], "alignment start")
            alignment_end = _integer(row["alignment_end"], "alignment end")
            # v0.6.1 obtains this from rust-htslib reference_end(): first base
            # after the alignment, unlike the loose wording in the prose docs.
            if not alignment_start <= position < alignment_end:
                raise ValueError("modkit reference position is outside its alignment")
            reference_entry = reference.entries.get(row["chrom"])
            if reference_entry is None:
                raise ValueError("Aligned chromosome is absent from the locked reference")
            if alignment_end > reference_entry.length:
                raise ValueError("modkit alignment extends beyond the locked reference")
            alignment = (row["chrom"], strand, read_length, alignment_start, alignment_end, flag)
            if read_alignments.setdefault(read_id, alignment) != alignment:
                raise ValueError("modkit contains inconsistent primary alignments for one read")
            if row["canonical_base"] not in {"A", "C", "G", "T"}:
                raise ValueError("Wildcard/RNA modification rows are outside the DNA adapter scope")
            if row["canonical_base"] != "C":
                collector.exclude("other_fundamental_base")
                continue
            if row["modified_primary_base"] != "C":
                raise ValueError("modkit modified primary base contradicts cytosine")
            code = row["mod_code"]
            if code not in upstream.modification_codes:
                raise ValueError("modkit cytosine modifications contradict the declared model")
            if row["inferred"] not in {"true", "false"}:
                raise ValueError("modkit inferred must be true or false")
            try:
                probability = float(row["mod_qual"])
            except ValueError as exc:
                raise ValueError("modkit mod_qual must be a probability in [0, 1]") from exc
            if not math.isfinite(probability) or not 0 <= probability <= 1:
                raise ValueError("modkit mod_qual must be a probability in [0, 1]")
            if row["inferred"] == "true" and probability != 0:
                raise ValueError("Implicit modkit probability must be zero")
            if not reference.is_cpg(row["chrom"], position, strand):
                collector.exclude("non_cpg_or_reference_mismatch")
                continue
            query_key = (read_id, query_position)
            ref_key = (chromosome, position, strand)
            if query_sites.setdefault(query_key, ref_key) != ref_key:
                raise ValueError("modkit maps one query base to multiple reference sites")
            key = (read_id, chromosome, position, strand)
            site = sites.setdefault(
                key, _ModkitSite(read_id, chromosome, position, strand, query_position)
            )
            if code in site.probabilities or site.query_position != query_position:
                raise ValueError("modkit duplicates a modification within one read/site")
            site.probabilities[code] = (
                None if row["inferred"] == "true" else _modkit_probability_interval(probability)
            )
        for site in sites.values():
            state = _classify(
                site.probabilities,
                upstream.modification_codes,
                adapter_policy.minimum_class_probability,
            )
            collector.add(site.read_id, site.chromosome, site.position, site.strand, state)
    return _finish(
        collector,
        path=path,
        reference_path=reference_path,
        source_id=source_id,
        genome_build=genome_build,
        upstream=upstream,
        adapter_policy=adapter_policy,
        fingerprint=fingerprint,
        reference=ref_fingerprint,
        input_format="modkit_extract_full_tsv",
        parser_name="ontseq_modkit_extract_full",
        parser_version=MODBAM_ADAPTER_VERSION,
    )


class ModkitExtractCommand(StrictModel):
    """Local plan only; callers must check version before execution.

    argv contains local paths and must not be included in public reports. The
    builder neither invokes modkit nor writes an output; use shell=False if a
    separate authorized runner executes this plan.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)
    schema_version: Literal["0.1.0"] = "0.1.0"
    tool_name: Literal["modkit"] = "modkit"
    expected_version: Literal["0.6.1"]
    version_check_argv: tuple[str, ...]
    argv: tuple[str, ...]
    probability_transform: Literal["none"] = "none"
    sensitive_output: Literal[True] = True


def build_modkit_extract_command(
    *,
    modbam_path: Path,
    reference_path: Path,
    output_path: Path,
    expected_version: str,
    threads: int = 1,
    executable: str = "modkit",
) -> ModkitExtractCommand:
    """Build a version-bound, shell-free extract-full CpG command; do not execute.

    No force, ignored modification class, sample, probability threshold or
    probability redistribution options are generated. Existing output must be
    protected by modkit's own no-overwrite default as well as this preflight.
    """
    if not executable or any(character in executable for character in "\r\n\x00"):
        raise ValueError("modkit executable is invalid")
    reviewed_version = _reviewed_modkit_version(expected_version)
    if type(threads) is not int or not 1 <= threads <= 256:
        raise ValueError("modkit threads must be an integer between 1 and 256")
    if not modbam_path.is_file() or not reference_path.is_file():
        raise ValueError("modkit command requires existing local BAM and reference inputs")
    if modbam_path.suffix.lower() != ".bam":
        raise ValueError("modkit extract command requires a BAM input")
    resolved_input = modbam_path.resolve()
    resolved_reference = reference_path.resolve()
    resolved_output = output_path.resolve()
    if resolved_output in {resolved_input, resolved_reference} or output_path.exists():
        raise ValueError("modkit command refuses to overwrite an existing input or output")
    return ModkitExtractCommand(
        expected_version=reviewed_version,
        version_check_argv=(executable, "--version"),
        argv=(
            executable,
            "extract",
            "full",
            str(resolved_input),
            str(resolved_output),
            "--ref",
            str(resolved_reference),
            "--mapped-only",
            "--cpg",
            "--threads",
            str(threads),
        ),
    )
