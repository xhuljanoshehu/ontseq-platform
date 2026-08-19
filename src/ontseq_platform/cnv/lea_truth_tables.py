"""Parsers for local historical evaluation tables accompanying Lea's ONTseq work.

The tables are validation inputs, not runtime rules.  This module preserves their raw
karyotype strings and historical classification labels without silently reconciling them or
turning them into clinical assertions.  The actual files remain local and must not be added
to Git.
"""

from __future__ import annotations

import csv
import io
import math
from collections.abc import Sequence
from dataclasses import dataclass


class LeaTruthTableError(ValueError):
    """Raised when a historical evaluation table is malformed or ambiguous."""


@dataclass(frozen=True)
class LeaGroundTruthRow:
    """One local cytogenetic truth record from the historical ``gt.tsv`` table."""

    sample_id: str
    karyotype: str


@dataclass(frozen=True)
class LeaEvaluationRow:
    """One local row from the historical ``gt_full.csv`` comparison table.

    Classification flags are preserved exactly as historical labels.  They are not
    recomputed and must not be used as hidden reportability rules by the production
    pipeline.
    """

    sample_id: str
    karyotype_cytogenetics: str
    karyotype_ont: str
    cellularity: float
    complex_karyotype: bool
    monosomal_karyotype: bool
    mrc: bool
    mrca: bool
    mra: bool


def _rows(
    lines: Sequence[str],
    *,
    required: set[str],
    delimiter: str,
    artifact_name: str,
) -> list[dict[str, str]]:
    text = "\n".join(line.rstrip("\n") for line in lines)
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    if reader.fieldnames is None:
        raise LeaTruthTableError(f"{artifact_name} contains no header")
    fieldnames = [name.lstrip("\ufeff").strip() for name in reader.fieldnames]
    if len(fieldnames) != len(set(fieldnames)):
        raise LeaTruthTableError(f"{artifact_name} contains duplicate column names")
    missing = sorted(required - set(fieldnames))
    if missing:
        raise LeaTruthTableError(
            f"{artifact_name} is missing required column(s): {', '.join(missing)}; "
            f"observed header: {fieldnames}"
        )

    result: list[dict[str, str]] = []
    for line_number, raw in enumerate(reader, start=2):
        if None in raw:
            raise LeaTruthTableError(
                f"{artifact_name} line {line_number} contains more fields than the header"
            )
        row = {
            key.lstrip("\ufeff").strip(): (value or "").strip()
            for key, value in raw.items()
        }
        if not any(row.values()):
            continue
        result.append(row)
    return result


def parse_lea_gt_tsv(lines: Sequence[str]) -> list[LeaGroundTruthRow]:
    """Parse ``gt.tsv`` while deliberately leaving ISCN interpretation to ``cnv.truth``."""
    rows = _rows(
        lines,
        required={"sample", "iscn"},
        delimiter="\t",
        artifact_name="gt.tsv",
    )
    parsed: list[LeaGroundTruthRow] = []
    seen: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        sample_id = row["sample"]
        karyotype = row["iscn"]
        if not sample_id or not karyotype:
            raise LeaTruthTableError(
                f"gt.tsv line {line_number} requires non-empty sample and iscn fields"
            )
        if sample_id in seen:
            raise LeaTruthTableError(f"gt.tsv contains duplicate sample {sample_id!r}")
        seen.add(sample_id)
        parsed.append(LeaGroundTruthRow(sample_id=sample_id, karyotype=karyotype))
    return parsed


def _finite_fraction(raw: str, *, line_number: int) -> float:
    try:
        value = float(raw)
    except ValueError as error:
        raise LeaTruthTableError(
            f"gt_full.csv line {line_number} has non-numeric cellularity {raw!r}"
        ) from error
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise LeaTruthTableError(
            f"gt_full.csv line {line_number} has cellularity outside [0, 1]: {raw!r}"
        )
    return value


def _binary_flag(raw: str, *, line_number: int, field: str) -> bool:
    if raw == "1":
        return True
    if raw == "0":
        return False
    raise LeaTruthTableError(
        f"gt_full.csv line {line_number} requires binary 0/1 for {field}, observed {raw!r}"
    )


def parse_lea_gt_full_csv(lines: Sequence[str]) -> list[LeaEvaluationRow]:
    """Parse ``gt_full.csv`` and preserve every historical label without reinterpretation."""
    required = {
        "sample_name",
        "karyotype_cg",
        "karyotype_ont",
        "cellularity",
        "complex",
        "monosomal",
        "mrc",
        "mrca",
        "mra",
    }
    rows = _rows(
        lines,
        required=required,
        delimiter=",",
        artifact_name="gt_full.csv",
    )
    parsed: list[LeaEvaluationRow] = []
    seen: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        sample_id = row["sample_name"]
        karyotype_cg = row["karyotype_cg"]
        karyotype_ont = row["karyotype_ont"]
        if not sample_id or not karyotype_cg or not karyotype_ont:
            raise LeaTruthTableError(
                f"gt_full.csv line {line_number} requires non-empty sample_name, "
                "karyotype_cg and karyotype_ont"
            )
        if sample_id in seen:
            raise LeaTruthTableError(
                f"gt_full.csv contains duplicate sample_name {sample_id!r}"
            )
        seen.add(sample_id)
        parsed.append(
            LeaEvaluationRow(
                sample_id=sample_id,
                karyotype_cytogenetics=karyotype_cg,
                karyotype_ont=karyotype_ont,
                cellularity=_finite_fraction(row["cellularity"], line_number=line_number),
                complex_karyotype=_binary_flag(
                    row["complex"], line_number=line_number, field="complex"
                ),
                monosomal_karyotype=_binary_flag(
                    row["monosomal"], line_number=line_number, field="monosomal"
                ),
                mrc=_binary_flag(row["mrc"], line_number=line_number, field="mrc"),
                mrca=_binary_flag(row["mrca"], line_number=line_number, field="mrca"),
                mra=_binary_flag(row["mra"], line_number=line_number, field="mra"),
            )
        )
    return parsed
