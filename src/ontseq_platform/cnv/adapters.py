"""Normalization of third-party CNV output into the shared call-set contract.

Design decision: column names, not column positions
---------------------------------------------------

Every parser here maps columns by *header name* through a declarative
:class:`ColumnMapping`, and raises when a required column is absent. Positional parsing
of a tool's tabular output is the classic source of silent scientific error: a caller
adds a column in a patch release, every downstream number shifts by one field, and the
pipeline keeps producing plausible-looking results. A loud failure on an unrecognised
header is strictly preferable.

Adding a method therefore usually means adding a :class:`ColumnMapping`, not writing a
parser. The mappings below cover the generic IGV ``SEG`` interchange format and ichorCNA.
Mappings for Spectre and QDNAseq/ACE are intentionally **not** shipped as verified
defaults, because their exact column layouts were not confirmed against the upstream
sources while this module was written; see ``docs/CNV_BENCHMARKING.md``. Supply a mapping
explicitly and record it in provenance.

None of these adapters executes a tool. Execution belongs behind the repository's
existing adapter boundary with version pinning and argument-vector invocation, as done
for Sniffles2 and Mosdepth.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ..models import GenomeBuild, ModuleRunStatus, ToolRecord
from .intervals import canonical_contig
from .models import (
    CnvCallSet,
    CnvDataBasis,
    CnvSegment,
    GenomicRegion,
)
from .states import CopyNumberState, state_from_copy_number

CANONICAL_CONTIGS = {str(index) for index in range(1, 23)} | {"X", "Y"}


@dataclass(frozen=True)
class ColumnMapping:
    """Declarative description of one tool's segment table."""

    contig: str
    start: str
    end: str
    #: At least one quantitative column must be present.
    copy_number: str | None = None
    log2_ratio: str | None = None
    copy_ratio: str | None = None
    categorical_call: str | None = None
    supporting_bins: str | None = None
    #: Some tools emit one-based inclusive starts; BED-style tools emit zero-based.
    one_based_start: bool = False
    #: Categorical call values mapped onto the state vocabulary.
    call_vocabulary: Mapping[str, CopyNumberState] | None = None
    comment_prefixes: tuple[str, ...] = ("#",)
    delimiter: str = "\t"

    def __post_init__(self) -> None:
        if not any((self.copy_number, self.log2_ratio, self.copy_ratio, self.categorical_call)):
            raise ValueError(
                "a column mapping must define at least one of copy_number, log2_ratio, "
                "copy_ratio or categorical_call"
            )


#: Generic IGV SEG interchange format. Widely produced by array and NGS segmenters.
SEG_MAPPING = ColumnMapping(
    contig="chrom",
    start="loc.start",
    end="loc.end",
    log2_ratio="seg.mean",
    supporting_bins="num.mark",
    one_based_start=True,
)

#: ichorCNA ``.cna.seg`` / ``.seg`` columns.
#:
#: Derived from the ichorCNA output documentation. ``Corrected_Copy_Number`` is preferred
#: over ``copy.number`` because it is the ploidy-corrected value, and ``logR_Copy_Number``
#: is deliberately not used as the primary quantity since it is a derived log-ratio
#: estimate rather than the model's integer call.
ICHORCNA_MAPPING = ColumnMapping(
    contig="chrom",
    start="start",
    end="end",
    copy_number="Corrected_Copy_Number",
    log2_ratio="seg.median.logR",
    categorical_call="Corrected_Call",
    one_based_start=True,
    call_vocabulary={
        "HOMD": CopyNumberState.HOMOZYGOUS_LOSS,
        "HETD": CopyNumberState.LOSS,
        "NEUT": CopyNumberState.NEUTRAL,
        "GAIN": CopyNumberState.GAIN,
        "AMP": CopyNumberState.HIGH_AMPLIFICATION,
        "HLAMP": CopyNumberState.HIGH_AMPLIFICATION,
    },
)


class SegmentParseError(ValueError):
    """Raised when a segment table cannot be parsed without ambiguity."""


def _resolve_header(header: Sequence[str], mapping: ColumnMapping) -> dict[str, int]:
    """Map required logical fields onto physical column indices."""
    lookup = {name.strip(): index for index, name in enumerate(header)}
    resolved: dict[str, int] = {}
    required = {"contig": mapping.contig, "start": mapping.start, "end": mapping.end}
    optional = {
        "copy_number": mapping.copy_number,
        "log2_ratio": mapping.log2_ratio,
        "copy_ratio": mapping.copy_ratio,
        "categorical_call": mapping.categorical_call,
        "supporting_bins": mapping.supporting_bins,
    }
    missing = [name for name, column in required.items() if column not in lookup]
    if missing:
        raise SegmentParseError(
            "segment table is missing required column(s) "
            + ", ".join(f"{name}={required[name]!r}" for name in missing)
            + f"; observed header: {list(header)}"
        )
    for name, column in required.items():
        resolved[name] = lookup[column]
    for name, column in optional.items():
        if column is not None and column in lookup:
            resolved[name] = lookup[column]
    if not any(
        key in resolved for key in ("copy_number", "log2_ratio", "copy_ratio", "categorical_call")
    ):
        raise SegmentParseError(
            "segment table contains no quantitative or categorical column from the "
            f"mapping; observed header: {list(header)}"
        )
    return resolved


def _optional_float(fields: Sequence[str], index: int | None) -> float | None:
    if index is None or index >= len(fields):
        return None
    raw = fields[index].strip()
    if raw in {"", ".", "NA", "NaN", "nan", "null"}:
        return None
    try:
        value = float(raw)
    except ValueError:
        return None
    return None if math.isnan(value) or math.isinf(value) else value


def parse_segment_table(
    lines: Sequence[str],
    mapping: ColumnMapping,
    *,
    baseline_ploidy: float = 2.0,
    skip_non_canonical_contigs: bool = True,
) -> tuple[list[CnvSegment], list[str]]:
    """Parse a segment table into CNV segments, returning parse warnings alongside.

    Copy number is taken directly when a mapping supplies it. Otherwise it is derived
    from a log2 ratio as ``ploidy * 2**log2``, which assumes the tool's ratio is relative
    to the sample's own baseline; that assumption is recorded as a warning so it can never
    be mistaken for a measured absolute copy number.
    """
    warnings: list[str] = []
    header: list[str] | None = None
    segments: list[CnvSegment] = []
    skipped_contigs: set[str] = set()
    derived_from_ratio = 0

    for number, raw in enumerate(lines, start=1):
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith(mapping.comment_prefixes) and header is not None:
            continue
        fields = [item.strip() for item in line.lstrip("#").split(mapping.delimiter)]
        if header is None:
            header = fields
            indices = _resolve_header(header, mapping)
            continue
        if len(fields) < 3:
            raise SegmentParseError(f"segment table line {number} has fewer than three columns")

        contig_raw = fields[indices["contig"]]
        contig = canonical_contig(contig_raw)
        if contig not in CANONICAL_CONTIGS:
            if skip_non_canonical_contigs:
                skipped_contigs.add(contig_raw)
                continue
            raise SegmentParseError(f"unsupported contig {contig_raw!r} on line {number}")
        try:
            start = int(float(fields[indices["start"]]))
            end = int(float(fields[indices["end"]]))
        except ValueError as error:
            raise SegmentParseError(
                f"segment table line {number} has non-numeric coordinates"
            ) from error
        if mapping.one_based_start:
            start -= 1
        if end <= start:
            raise SegmentParseError(
                f"segment table line {number} has an empty or inverted interval"
            )

        copy_number = _optional_float(fields, indices.get("copy_number"))
        log2_ratio = _optional_float(fields, indices.get("log2_ratio"))
        copy_ratio = _optional_float(fields, indices.get("copy_ratio"))
        bins_value = _optional_float(fields, indices.get("supporting_bins"))

        if copy_number is None and copy_ratio is not None:
            copy_number = baseline_ploidy * copy_ratio
            derived_from_ratio += 1
        if copy_number is None and log2_ratio is not None:
            copy_number = baseline_ploidy * (2.0**log2_ratio)
            derived_from_ratio += 1

        state: CopyNumberState | None = None
        call_index = indices.get("categorical_call")
        if call_index is not None and mapping.call_vocabulary and call_index < len(fields):
            state = mapping.call_vocabulary.get(fields[call_index].strip().upper())
        if state is None:
            if copy_number is None:
                raise SegmentParseError(
                    f"segment table line {number} carries neither an interpretable call "
                    "nor a quantitative value"
                )
            state = state_from_copy_number(copy_number, baseline_ploidy=baseline_ploidy)

        segments.append(
            CnvSegment(
                contig=contig,
                start=start,
                end=end,
                state=state,
                copy_number=None if copy_number is None or copy_number < 0 else copy_number,
                copy_ratio=copy_ratio,
                log2_ratio=log2_ratio,
                supporting_bins=int(bins_value) if bins_value is not None else None,
            )
        )

    if header is None:
        raise SegmentParseError("segment table contained no header row")
    if skipped_contigs:
        warnings.append(
            "Skipped non-canonical contigs: " + ", ".join(sorted(skipped_contigs))
        )
    if derived_from_ratio:
        warnings.append(
            f"{derived_from_ratio} segment(s) had no absolute copy number; it was derived "
            "from the reported ratio assuming the ratio is relative to this sample's own "
            "baseline. Treat the value as relative, not as measured absolute copy number."
        )
    segments.sort(key=lambda item: (item.contig, item.start, item.end))
    _reject_overlaps(segments)
    return segments, warnings


def _reject_overlaps(segments: Sequence[CnvSegment]) -> None:
    by_contig: dict[str, list[CnvSegment]] = {}
    for segment in segments:
        by_contig.setdefault(segment.contig, []).append(segment)
    for contig, items in by_contig.items():
        for previous, current in zip(items, items[1:], strict=False):
            if current.start < previous.end:
                raise SegmentParseError(
                    f"overlapping segments on contig {contig}: "
                    f"[{previous.start}, {previous.end}) and [{current.start}, {current.end}). "
                    "A copy-number call set must be a partition, not a list of candidates."
                )


def call_set_from_segment_table(
    lines: Sequence[str],
    mapping: ColumnMapping,
    *,
    call_set_id: str,
    sample_id: str,
    genome_build: GenomeBuild,
    method: str,
    method_version: str,
    data_basis: CnvDataBasis,
    closed_world: bool,
    baseline_ploidy: float = 2.0,
    bin_size_bp: int | None = None,
    estimated_tumor_fraction: float | None = None,
    estimated_ploidy: float | None = None,
    mean_coverage_x: float | None = None,
    no_call_regions: Sequence[GenomicRegion] = (),
    tool: ToolRecord | None = None,
    extra_limitations: Sequence[str] = (),
) -> CnvCallSet:
    """Parse a tool's segment table into a normalized, non-reportable call set."""
    segments, warnings = parse_segment_table(
        lines, mapping, baseline_ploidy=baseline_ploidy
    )
    status = ModuleRunStatus.COMPLETED if segments else ModuleRunStatus.NO_CALL
    return CnvCallSet(
        call_set_id=call_set_id,
        sample_id=sample_id,
        genome_build=genome_build,
        method=method,
        method_version=method_version,
        data_basis=data_basis,
        background_state=(
            CopyNumberState.NEUTRAL if closed_world else CopyNumberState.NO_CALL
        ),
        status=status,
        segments=segments,
        no_call_regions=list(no_call_regions),
        bin_size_bp=bin_size_bp,
        estimated_tumor_fraction=estimated_tumor_fraction,
        estimated_ploidy=estimated_ploidy,
        mean_coverage_x=mean_coverage_x,
        tool=tool,
        warnings=warnings,
        limitations=[
            "Normalized from a third-party segment table. No claim is made about the "
            "upstream tool's correctness, parameters or suitability for this assay.",
            "The call set is research-only and not reportable.",
            *extra_limitations,
        ],
    )
