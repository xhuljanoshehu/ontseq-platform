"""Aggregate-only audit of historical cytogenetic truth convertibility.

The audit answers a narrow validation question: how much of a historical ``gt.tsv`` table
can the current CNV truth converter represent without guessing?  It deliberately emits no
sample identifiers and no per-sample karyotypes.  Unsupported constructs remain counts with
reason categories, so the output can be committed or attached to engineering review without
copying the local truth table itself.
"""

from __future__ import annotations

import argparse
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from ..models import GenomeBuild, StrictModel
from .cytobands import CytobandTable, load_cytoband_file
from .lea_io import sha256_file
from .lea_truth_tables import LeaGroundTruthRow, parse_lea_gt_tsv
from .truth import convert_karyotype


class LeaTruthAuditSummary(StrictModel):
    """Privacy-minimized aggregate result of historical truth conversion."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    source: Literal["lea-historical-gt.tsv"] = "lea-historical-gt.tsv"
    genome_build: Literal[GenomeBuild.GRCH37] = GenomeBuild.GRCH37
    total_rows: int = Field(ge=0)
    fully_convertible_rows: int = Field(ge=0)
    incomplete_rows: int = Field(ge=0)
    conversion_fraction: float | None = Field(default=None, ge=0, le=1)
    unsupported_construct_count: int = Field(ge=0)
    unsupported_reason_counts: dict[str, int] = Field(default_factory=dict)
    balanced_construct_count: int = Field(ge=0)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    contains_sample_identifiers: Literal[False] = False
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def counts_reconcile(self) -> LeaTruthAuditSummary:
        if self.fully_convertible_rows + self.incomplete_rows != self.total_rows:
            raise ValueError("truth-audit row counts do not reconcile")
        if sum(self.unsupported_reason_counts.values()) != self.unsupported_construct_count:
            raise ValueError("truth-audit unsupported reason counts do not reconcile")
        return self


def _reason_category(reason: str) -> str:
    lowered = reason.lower()
    if "uncertainty marker" in lowered:
        return "uncertainty_marker"
    if "sex-chromosome complement" in lowered:
        return "sex_chromosome_complement"
    if "cytoband" in lowered or "band " in lowered or "contig " in lowered:
        return "reference_or_band_mapping"
    if "construct is not supported" in lowered:
        return "unsupported_construct"
    return "other"


def audit_lea_ground_truth(
    rows: Sequence[LeaGroundTruthRow],
    cytobands: CytobandTable,
) -> LeaTruthAuditSummary:
    """Audit all rows without exporting their sample IDs or karyotype strings."""
    if cytobands.genome_build != GenomeBuild.GRCH37:
        raise ValueError("Lea historical ground-truth audit requires a GRCh37 cytoband table")

    complete = 0
    incomplete = 0
    unsupported_count = 0
    balanced_count = 0
    reason_counts: Counter[str] = Counter()
    for row in rows:
        conversion = convert_karyotype(row.karyotype, cytobands)
        balanced_count += len(conversion.balanced_constructs)
        if conversion.unsupported:
            incomplete += 1
            unsupported_count += len(conversion.unsupported)
            reason_counts.update(_reason_category(item.reason) for item in conversion.unsupported)
        else:
            complete += 1

    total = len(rows)
    return LeaTruthAuditSummary(
        total_rows=total,
        fully_convertible_rows=complete,
        incomplete_rows=incomplete,
        conversion_fraction=(complete / total if total else None),
        unsupported_construct_count=unsupported_count,
        unsupported_reason_counts=dict(sorted(reason_counts.items())),
        balanced_construct_count=balanced_count,
        warnings=(["The local truth table contained no rows."] if total == 0 else []),
        limitations=[
            "Convertibility is a software capability measurement, not agreement with ONT and "
            "not a clinical performance metric.",
            "A fully convertible karyotype can still be limited by cytogenetic resolution, "
            "clone flattening and constructs that assert balanced rather than dosage change.",
            "The aggregate intentionally contains no sample identifiers or karyotype strings.",
        ],
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m ontseq_platform.cnv.lea_audit",
        description=(
            "Audit how many rows of a local Lea historical gt.tsv are fully representable "
            "by the current CNV truth converter. Output is aggregate-only."
        ),
    )
    parser.add_argument("--gt-tsv", type=Path, required=True)
    parser.add_argument("--cytobands", type=Path, required=True)
    parser.add_argument("--cytoband-sha256", required=True)
    parser.add_argument("--cytoband-resource-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    args = _parser().parse_args()
    observed_sha256 = sha256_file(args.cytobands)
    if observed_sha256 != args.cytoband_sha256.lower():
        raise SystemExit("ERROR: cytoband SHA-256 does not match the deployment lock")
    table = load_cytoband_file(
        args.cytobands,
        genome_build=GenomeBuild.GRCH37,
        resource_id=args.cytoband_resource_id,
        source_sha256=observed_sha256,
    )
    rows = parse_lea_gt_tsv(args.gt_tsv.read_text(encoding="utf-8-sig").splitlines())
    summary = audit_lea_ground_truth(rows, table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(summary.model_dump_json(indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
