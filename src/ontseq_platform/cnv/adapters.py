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
parser. The mappings below cover the generic IGV ``SEG`` interchange format, ichorCNA and
the QDNAseq/ACE lane this repository runs.

The QDNAseq/ACE mapping is a narrower claim than the other two. It is not a mapping for
"QDNAseq output" in general: it describes the table written by ``scripts/run_qdnaseq_ace.R``
in this repository, whose columns are chosen there rather than by the upstream package. It
is shipped because that layout is ours and is checksummed into every run's provenance —
which was exactly what could not be said when this module was first written. Output from a
differently configured QDNAseq or ACE installation still needs its own mapping, supplied
explicitly and recorded in provenance.

A mapping for Spectre is still deliberately absent for the original reason: its column
layout was not confirmed against the upstream source. See ``docs/CNV_BENCHMARKING.md``.

None of these adapters executes a tool. Execution belongs behind the repository's
existing adapter boundary with version pinning and argument-vector invocation, as done
for Sniffles2 and Mosdepth.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from ..models import GenomeBuild, ModuleRunStatus, ToolRecord
from .intervals import canonical_contig, subtract
from .models import (
    CnvCallSet,
    CnvDataBasis,
    CnvSegment,
    GenomicRegion,
)
from .qdnaseq import QDNAseqCallReport
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


#: The segment table written by ``scripts/run_qdnaseq_ace.R`` in this repository.
#:
#: ``absolute_copy_number`` is ACE's purity- and ploidy-adjusted estimate, taken as the
#: primary quantity rather than ``call``: ``call`` is a rounded band and rounding is the
#: step at which a shallow gain and a neutral region become indistinguishable. The
#: quantitative column keeps that distinction available to the scorer.
#:
#: Coordinates are one-based inclusive, because they come from QDNAseq bin annotations.
QDNASEQ_ACE_MAPPING = ColumnMapping(
    contig="chromosome",
    start="start",
    end="end",
    copy_number="absolute_copy_number",
    supporting_bins="bin_count",
    one_based_start=True,
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
    for name, optional_column in optional.items():
        if optional_column is not None and optional_column in lookup:
            resolved[name] = lookup[optional_column]
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
        warnings.append("Skipped non-canonical contigs: " + ", ".join(sorted(skipped_contigs)))
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
    segments, warnings = parse_segment_table(lines, mapping, baseline_ploidy=baseline_ploidy)
    # An empty segment table is genuinely ambiguous from the outside: the format records no
    # difference between a tool that ran and found nothing and one that produced nothing.
    # A caller that knows its own run can assert `reports_biological_negative`; an adapter
    # reading its output afterwards cannot, so NO_CALL is the only honest reading here.
    status = ModuleRunStatus.COMPLETED if segments else ModuleRunStatus.NO_CALL
    return CnvCallSet(
        call_set_id=call_set_id,
        sample_id=sample_id,
        genome_build=genome_build,
        method=method,
        method_version=method_version,
        data_basis=data_basis,
        background_state=(CopyNumberState.NEUTRAL if closed_world else CopyNumberState.NO_CALL),
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


def qdnaseq_method_version(tools: Sequence[ToolRecord]) -> str:
    """Compose a method version from the tool records the runtime lane already recorded.

    Both packages are named because either one moving changes the answer: QDNAseq decides
    the bins and the correction, ACE decides the purity/ploidy fit the absolute copy
    numbers are expressed in. A benchmark result attributed to "QDNAseq" alone could not be
    reproduced from that label.
    """
    versions = {record.name: record.version for record in tools}
    missing = [name for name in ("QDNAseq", "ACE") if name not in versions]
    if missing:
        raise SegmentParseError(
            "cannot attribute a QDNAseq/ACE call set without both tool versions; "
            f"missing {', '.join(missing)}"
        )
    return f"QDNAseq {versions['QDNAseq']}+ACE {versions['ACE']}"


def _uncovered_regions(
    segments: Sequence[CnvSegment],
    contig_lengths: Mapping[str, int],
) -> list[GenomicRegion]:
    """Every canonical base the segment table does not speak about, as explicit no-calls.

    QDNAseq drops bins it cannot correct, and the runner keeps only the autosomes. Both are
    real limits on where the method can answer at all, and both are invisible in a segment
    table — the rows that would have said so are simply absent. Turning the gaps into
    declared no-call regions is what keeps them out of the denominator instead of being
    scored as agreement with whatever the truth set happens to assert there.
    """
    covered: dict[str, list[tuple[int, int]]] = {}
    for segment in segments:
        covered.setdefault(canonical_contig(segment.contig), []).append(
            (segment.start, segment.end)
        )
    regions: list[GenomicRegion] = []
    for contig, length in sorted(contig_lengths.items()):
        canonical = canonical_contig(contig)
        if canonical not in CANONICAL_CONTIGS or length <= 0:
            continue
        for start, end in subtract([(0, length)], covered.get(canonical, [])):
            regions.append(
                GenomicRegion(
                    contig=canonical,
                    start=start,
                    end=end,
                    label="not covered by the QDNAseq segmentation",
                )
            )
    return regions


def call_set_from_qdnaseq_report(
    report: QDNAseqCallReport,
    *,
    call_set_id: str,
    data_basis: CnvDataBasis,
    output_dir: Path,
    contig_lengths: Mapping[str, int] | None = None,
    bin_size_kbp: int | None = None,
    mean_coverage_x: float | None = None,
    extra_limitations: Sequence[str] = (),
) -> CnvCallSet:
    """Normalize one QDNAseq/ACE run into the benchmark's call-set contract.

    This is the seam the CNV direction asked for: the runtime lane is measured *through*
    the existing benchmark architecture rather than being promoted alongside it. Nothing
    here selects QDNAseq, marks it preferred or makes it reportable — ``CnvCallSet`` fixes
    ``reportable`` to ``False`` and no argument can change that. What the function produces
    is a scoreable object with its provenance attached.

    ``data_basis`` has no default on purpose. An adaptive-sampling run yields two read
    populations whose depth behaviour is not comparable, and a run that pooled them is a
    third thing again; guessing which one a report came from would silently place it in the
    wrong benchmark stratum. The caller knows, so the caller states it.
    """
    fit = report.primary_fit
    if bin_size_kbp is not None:
        matching = [item for item in report.fits if item.bin_size_kbp == bin_size_kbp]
        if not matching:
            available = ", ".join(str(item.bin_size_kbp) for item in report.fits)
            raise SegmentParseError(
                f"no QDNAseq fit at {bin_size_kbp} kbp in this report; available: {available}"
            )
        fit = matching[0]

    segment_path = Path(output_dir) / fit.segment_file
    if not segment_path.is_file():
        raise SegmentParseError(f"QDNAseq segment table not found: {segment_path}")
    lines = segment_path.read_text(encoding="utf-8").splitlines()
    segments, warnings = parse_segment_table(lines, QDNASEQ_ACE_MAPPING)

    no_call_regions = _uncovered_regions(segments, contig_lengths) if contig_lengths else []
    if contig_lengths is None:
        warnings.append(
            "No contig lengths were supplied, so the regions this segmentation does not "
            "cover — dropped bins and the sex chromosomes the runner excludes — could not "
            "be declared as no-calls. Score against a mask that excludes them explicitly."
        )

    # A full-partition segmenter that produced nothing did not look; it is not a negative.
    status = ModuleRunStatus.COMPLETED if segments else ModuleRunStatus.NO_CALL

    return CnvCallSet(
        call_set_id=call_set_id,
        sample_id=report.sample_id,
        genome_build=report.genome_build,
        method="QDNAseq+ACE",
        method_version=qdnaseq_method_version(report.tools),
        data_basis=data_basis,
        # QDNAseq emits a segment for every bin it kept, neutral ones included, so silence
        # inside a covered region asserts neutrality rather than absence of information.
        # Everything it could not cover is named above instead of being left to inference.
        background_state=CopyNumberState.NEUTRAL,
        status=status,
        segments=segments,
        no_call_regions=no_call_regions,
        bin_size_bp=fit.bin_size_kbp * 1000,
        estimated_tumor_fraction=fit.cellularity,
        estimated_ploidy=fit.ploidy,
        mean_coverage_x=mean_coverage_x,
        tool=ToolRecord(
            name="QDNAseq+ACE",
            version=qdnaseq_method_version(report.tools),
            parameters={
                "bin_size_kbp": fit.bin_size_kbp,
                "fit_error": fit.fit_error,
                "candidate_count": fit.candidate_count,
                "segment_file": fit.segment_file,
                "primary_fit_bin_size_kbp": report.primary_fit.bin_size_kbp,
                "available_bin_sizes_kbp": [item.bin_size_kbp for item in report.fits],
            },
        ),
        warnings=[*warnings, *report.warnings],
        limitations=[
            "Normalized from this repository's QDNAseq/ACE runner. The call set is "
            "research-only and not reportable.",
            "Copy numbers are ACE's purity- and ploidy-adjusted estimates. They depend on "
            "the fit ACE selected, and a different penalty or ploidy grid would move them.",
            "The bin size, ACE penalty and ploidy grid used are engineering parameters "
            "carried from the run's policy. None of them is a validated threshold, and a "
            "benchmark result does not turn one into one.",
            "Only one fit is scored here. A report holding several resolutions is several "
            "candidate methods, not one method measured several ways.",
            *report.limitations,
            *extra_limitations,
        ],
    )
