from __future__ import annotations

import gzip
import math
from collections import defaultdict
from pathlib import Path
from typing import Literal, TextIO

from pydantic import Field, model_validator

from .marlin_contracts import MarlinBridgeLock, MarlinProbeObservation, MarlinSourceKind
from .methylation import ModificationCode, parse_bedmethyl
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


class MarlinModkitProbeInput(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    adapter_id: Literal["marlin_modkit_probe_v1"] = "marlin_modkit_probe_v1"
    source_kind: Literal[MarlinSourceKind.MODKIT_DERIVED] = MarlinSourceKind.MODKIT_DERIVED
    genome_build: GenomeBuild
    input_fingerprint: FileFingerprint
    probe_resource_fingerprint: FileFingerprint
    bridge_lock_id: str = Field(min_length=3)
    modkit_version: Literal["0.6.4"] = "0.6.4"
    observations: list[MarlinProbeObservation] = Field(min_length=1)

    @model_validator(mode="after")
    def validate_native_bridge_contract(self) -> MarlinModkitProbeInput:
        if self.genome_build is not GenomeBuild.GRCH37:
            raise ValueError("MARLIN native modkit bridge currently requires GRCh37/hg19")
        probe_ids = [item.probe_id for item in self.observations]
        if len(probe_ids) != len(set(probe_ids)):
            raise ValueError("MARLIN native bridge observations require unique probe identifiers")
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
                raise ValueError(f"Line {line_number}: expected exactly 5 tab-separated fields")
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


def _load_marlin_probe_positions(
    probe_resource: Path,
) -> tuple[dict[str, tuple[tuple[str, int], ...]], FileFingerprint]:
    """Load the upstream MARLIN probe resource using its published V1/V2/V4 semantics."""
    path = Path(probe_resource)
    if not path.is_file():
        raise ValueError("MARLIN probe resource is missing")
    initial_sha256 = sha256_file(path)
    initial_size = path.stat().st_size
    positions: dict[str, list[tuple[str, int]]] = defaultdict(list)
    seen_rows: set[tuple[str, int, str]] = set()

    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"MARLIN probe-resource line {line_number} requires at least 4 columns"
                )
            chromosome, start_text, _ignored_end, probe_id = fields[:4]
            if any(value != value.strip() for value in (chromosome, start_text, probe_id)):
                raise ValueError(
                    f"MARLIN probe-resource line {line_number} contains surrounding whitespace"
                )
            if not probe_id:
                raise ValueError(f"MARLIN probe-resource line {line_number} has empty probe ID")
            start = _parse_int(
                start_text,
                line_number=line_number,
                field_name="probe-resource start",
            )
            if start < 0:
                raise ValueError(
                    f"MARLIN probe-resource line {line_number} start must be non-negative"
                )
            try:
                MarlinProbeObservation(
                    chromosome=chromosome,
                    start=start,
                    end=start + 1,
                    methylation_fraction=None,
                    probe_id=probe_id,
                )
            except ValueError as exc:
                raise ValueError(
                    f"MARLIN probe-resource line {line_number} has invalid chromosome/probe"
                ) from exc
            row_key = (chromosome, start, probe_id)
            if row_key in seen_rows:
                raise ValueError("MARLIN probe resource contains a duplicate probe position")
            seen_rows.add(row_key)
            positions[probe_id].append((chromosome, start))

    if not positions:
        raise ValueError("MARLIN probe resource contains no probe positions")
    for probe_id, probe_positions in positions.items():
        if len({chromosome for chromosome, _start in probe_positions}) != 1:
            raise ValueError(f"MARLIN probe {probe_id!r} maps to multiple chromosomes")

    final_sha256 = sha256_file(path)
    final_size = path.stat().st_size
    if final_sha256 != initial_sha256 or final_size != initial_size:
        raise ValueError("MARLIN probe resource changed while it was being parsed")

    frozen = {probe_id: tuple(values) for probe_id, values in positions.items()}
    return frozen, FileFingerprint(size_bytes=initial_size, sha256=initial_sha256)


def parse_marlin_modkit_probe_input(
    path: Path,
    *,
    probe_resource: Path,
    bridge_lock: MarlinBridgeLock,
    modkit_version: str,
) -> MarlinModkitProbeInput:
    """Map locked modkit bedMethyl calls to MARLIN probes using the published aggregation."""
    if bridge_lock.genome_build is not GenomeBuild.GRCH37:
        raise ValueError("MARLIN native bridge currently requires GRCh37/hg19")
    if modkit_version != bridge_lock.modkit_version:
        raise ValueError("MARLIN native bridge modkit version differs from bridge lock")

    probe_resource_path = Path(probe_resource)
    if not probe_resource_path.is_file():
        raise ValueError("MARLIN probe resource is missing")
    if sha256_file(probe_resource_path) != bridge_lock.probe_resource_sha256:
        raise ValueError("MARLIN native bridge probe-resource SHA-256 differs from bridge lock")
    probe_positions, probe_fingerprint = _load_marlin_probe_positions(probe_resource_path)
    if probe_fingerprint.sha256 != bridge_lock.probe_resource_sha256:
        raise ValueError("MARLIN native bridge probe-resource SHA-256 changed during parsing")

    calls = Path(path)
    if not calls.is_file():
        raise ValueError("MARLIN native bridge bedMethyl input is missing")
    initial_calls_sha256 = sha256_file(calls)
    initial_calls_size = calls.stat().st_size
    sites, _skipped_non_canonical = parse_bedmethyl(
        calls,
        allowed_codes=[ModificationCode.FIVE_MC],
    )
    final_calls_sha256 = sha256_file(calls)
    final_calls_size = calls.stat().st_size
    if final_calls_sha256 != initial_calls_sha256 or final_calls_size != initial_calls_size:
        raise ValueError("MARLIN native bridge bedMethyl input changed while it was parsed")

    site_by_position = {(site.chromosome, site.start): site for site in sites}
    observations: list[MarlinProbeObservation] = []
    for probe_id, positions in probe_positions.items():
        chromosome = positions[0][0]
        starts = [start for _chromosome, start in positions]
        matched = [
            site_by_position[position] for position in positions if position in site_by_position
        ]
        total_valid = sum(site.valid_coverage for site in matched)
        total_modified = sum(site.modified_calls for site in matched)
        beta = total_modified / total_valid if total_valid else None
        observations.append(
            MarlinProbeObservation(
                chromosome=chromosome,
                start=min(starts),
                end=max(starts) + 1,
                methylation_fraction=beta,
                probe_id=probe_id,
            )
        )

    observations.sort(key=lambda item: (item.chromosome, item.start, item.probe_id))
    return MarlinModkitProbeInput(
        genome_build=GenomeBuild.GRCH37,
        input_fingerprint=FileFingerprint(
            size_bytes=initial_calls_size,
            sha256=initial_calls_sha256,
        ),
        probe_resource_fingerprint=probe_fingerprint,
        bridge_lock_id=bridge_lock.bridge_id,
        modkit_version=bridge_lock.modkit_version,
        observations=observations,
    )
