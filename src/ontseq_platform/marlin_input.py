from __future__ import annotations

import gzip
import math
from pathlib import Path
from typing import Literal, TextIO

from pydantic import Field, model_validator

from .marlin_contracts import MarlinProbeObservation, MarlinSourceKind
from .models import FileFingerprint, GenomeBuild, StrictModel
from .reference import sha256_file


class MarlinPrecomputedInput(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    adapter_id: Literal["marlin_probe_bed_v1"] = "marlin_probe_bed_v1"
    source_kind: Literal[MarlinSourceKind.PRECOMPUTED_METHYLATION] = (
        MarlinSourceKind.PRECOMPUTED_METHYLATION
    )
    genome_build: GenomeBuild
    input_fingerprint: FileFingerprint
    observations: list[MarlinProbeObservation] = Field(min_length=1)

    @model_validator(mode="after")
    def validate_v1_contract(self) -> MarlinPrecomputedInput:
        if self.genome_build is not GenomeBuild.GRCH37:
            raise ValueError("MARLIN probe-BED v1 currently requires GRCh37/hg19")
        probe_ids = [item.probe_id for item in self.observations]
        if len(probe_ids) != len(set(probe_ids)):
            raise ValueError("Duplicate MARLIN probe identifiers are not allowed")
        return self


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def _parse_int(value: str, *, line_number: int, field_name: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"Line {line_number}: {field_name} must be an integer") from exc
    return parsed


def _parse_fraction(value: str, *, line_number: int) -> float | None:
    if value == "NA":
        return None
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(
            f"Line {line_number}: methylation fraction must be a number or literal NA"
        ) from exc
    if not math.isfinite(parsed):
        raise ValueError(f"Line {line_number}: methylation fraction must be finite")
    if not 0 <= parsed <= 1:
        raise ValueError(f"Line {line_number}: methylation fraction must be between 0 and 1")
    return parsed


def parse_marlin_probe_bed(path: Path, *, genome_build: GenomeBuild) -> MarlinPrecomputedInput:
    """Parse the published MARLIN five-column probe representation without coercion."""
    path = Path(path)
    if genome_build is not GenomeBuild.GRCH37:
        raise ValueError("MARLIN probe-BED v1 currently requires GRCh37/hg19")
    if not path.is_file():
        raise ValueError(f"MARLIN probe input does not exist or is not a file: {path}")

    initial_sha256 = sha256_file(path)
    initial_size = path.stat().st_size
    observations: list[MarlinProbeObservation] = []
    seen_probe_ids: set[str] = set()

    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                raise ValueError(f"Line {line_number}: blank rows are not allowed")
            fields = line.split("\t")
            if len(fields) != 5:
                raise ValueError(
                    f"Line {line_number}: expected exactly 5 tab-separated fields"
                )
            if any(field != field.strip() for field in fields):
                raise ValueError(f"Line {line_number}: surrounding field whitespace is not allowed")

            chromosome, start_text, end_text, fraction_text, probe_id = fields
            if not probe_id:
                raise ValueError(f"Line {line_number}: MARLIN probe identifier is empty")
            if probe_id in seen_probe_ids:
                raise ValueError(f"Line {line_number}: Duplicate MARLIN probe {probe_id!r}")

            start = _parse_int(start_text, line_number=line_number, field_name="start")
            end = _parse_int(end_text, line_number=line_number, field_name="end")
            if start < 0:
                raise ValueError(f"Line {line_number}: start must be non-negative")
            if end <= start:
                raise ValueError(f"Line {line_number}: end must be greater than start")
            fraction = _parse_fraction(fraction_text, line_number=line_number)

            try:
                observation = MarlinProbeObservation(
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    methylation_fraction=fraction,
                    probe_id=probe_id,
                )
            except ValueError as exc:
                raise ValueError(
                    f"Line {line_number}: invalid chromosome or probe row: {exc}"
                ) from exc
            observations.append(observation)
            seen_probe_ids.add(probe_id)

    if not observations:
        raise ValueError("MARLIN probe input contains no MARLIN probe rows")

    final_sha256 = sha256_file(path)
    final_size = path.stat().st_size
    if final_sha256 != initial_sha256 or final_size != initial_size:
        raise ValueError("MARLIN probe input changed while it was being parsed")

    return MarlinPrecomputedInput(
        genome_build=genome_build,
        input_fingerprint=FileFingerprint(size_bytes=initial_size, sha256=initial_sha256),
        observations=observations,
    )
