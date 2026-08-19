"""Historical compatibility adapters for Lea Evers' ONTseq CNV outputs.

The goal is reproducibility and comparison, not code reuse. This module does not execute
Lea's workflow and does not copy its implementation. It normalizes the documented output
contracts of the historical QDNAseq + ACE path into ONTSeq's research-only CNV contract so
that the old workflow can be benchmarked on the same evaluable genome as newer methods.

The historical profile represented here is the supplied 2026 snapshot whose active
Snakemake configuration used GRCh37/hg19, QDNAseq bin size 1000 kbp, ACE penalty 0.6 and a
66% affected-band threshold. Those values are compatibility metadata, not validated
clinical thresholds.
"""

from __future__ import annotations

import csv
import io
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ..models import GenomeBuild, ModuleRunStatus, ToolRecord
from .cytobands import CytobandTable
from .intervals import canonical_contig
from .models import CnvCallSet, CnvDataBasis, CnvSegment, GenomicRegion
from .states import CopyNumberState

CANONICAL_CONTIGS = {str(index) for index in range(1, 23)} | {"X", "Y"}
AUTOSOMES = {str(index) for index in range(1, 23)}


class LeaCompatibilityError(ValueError):
    """Raised when historical output cannot be normalized without guessing."""


@dataclass(frozen=True)
class LeaAceCompatibilityProfile:
    """Frozen metadata for one reproducibility profile of Lea's ACE path.

    ``declared_*_version`` records versions found in the supplied ``renv.lock``. The
    supplied ``r-ace.def`` installed Bioconductor packages without version pins, so these
    values are deliberately labelled declarations rather than verified runtime versions.
    """

    profile_id: str = "lea-evers-ontseq-2026-ace-hg19"
    genome_build: GenomeBuild = GenomeBuild.GRCH37
    qdnaseq_bin_size_kbp: int = 1000
    ace_penalty: float = 0.6
    assumed_autosomal_ploidy: float = 2.0
    affected_band_fraction: float = 0.66
    declared_qdnaseq_version: str = "1.38.0"
    declared_qdnaseq_hg19_version: str = "1.32.0"
    declared_ace_version: str = "1.20.0"
    runtime_versions_verified: bool = False

    def validate(self) -> None:
        if self.genome_build != GenomeBuild.GRCH37:
            raise LeaCompatibilityError(
                "the frozen Lea ACE compatibility profile is GRCh37/hg19 only"
            )
        if self.qdnaseq_bin_size_kbp <= 0:
            raise LeaCompatibilityError("QDNAseq bin size must be positive")
        if self.ace_penalty < 0:
            raise LeaCompatibilityError("ACE penalty must not be negative")
        if self.assumed_autosomal_ploidy <= 0:
            raise LeaCompatibilityError("assumed autosomal ploidy must be positive")
        if not 0.0 <= self.affected_band_fraction <= 1.0:
            raise LeaCompatibilityError("affected-band fraction must be between 0 and 1")


LEA_ACE_2026_HG19 = LeaAceCompatibilityProfile()


@dataclass(frozen=True)
class LeaChromosomeCopy:
    """One row of the historical ``CN.csv`` output."""

    contig: str
    copies: float
    ploidy: float
    cna: float


@dataclass(frozen=True)
class LeaBandCall:
    """One row of the historical ``dels_dups.csv`` output."""

    contig: str
    band: str
    event: str
    affected_fraction: float


def _reader(
    lines: Sequence[str],
    *,
    required: set[str],
    delimiter: str,
    artifact_name: str,
) -> list[dict[str, str]]:
    """Read a small delimited table with strict header semantics."""
    text = "\n".join(line.rstrip("\n") for line in lines)
    stream = io.StringIO(text)
    reader = csv.DictReader(stream, delimiter=delimiter)
    if reader.fieldnames is None:
        raise LeaCompatibilityError(f"{artifact_name} contains no header")

    fieldnames = [name.lstrip("\ufeff").strip() for name in reader.fieldnames]
    if len(fieldnames) != len(set(fieldnames)):
        raise LeaCompatibilityError(f"{artifact_name} contains duplicate column names")
    missing = sorted(required - set(fieldnames))
    if missing:
        raise LeaCompatibilityError(
            f"{artifact_name} is missing required column(s): {', '.join(missing)}; "
            f"observed header: {fieldnames}"
        )

    rows: list[dict[str, str]] = []
    for line_number, raw_row in enumerate(reader, start=2):
        if None in raw_row:
            raise LeaCompatibilityError(
                f"{artifact_name} line {line_number} contains more fields than the header"
            )
        normalized = {
            key.lstrip("\ufeff").strip(): (value or "").strip()
            for key, value in raw_row.items()
        }
        if not any(normalized.values()):
            continue
        rows.append(normalized)
    return rows


def _finite_float(raw: str, *, artifact_name: str, line_number: int, field: str) -> float:
    try:
        value = float(raw)
    except ValueError as error:
        raise LeaCompatibilityError(
            f"{artifact_name} line {line_number} has a non-numeric {field}: {raw!r}"
        ) from error
    if not math.isfinite(value):
        raise LeaCompatibilityError(
            f"{artifact_name} line {line_number} has a non-finite {field}: {raw!r}"
        )
    return value


def _canonical_human_contig(raw: str, *, artifact_name: str, line_number: int) -> str:
    contig = canonical_contig(raw)
    if contig not in CANONICAL_CONTIGS:
        raise LeaCompatibilityError(
            f"{artifact_name} line {line_number} contains unsupported contig {raw!r}"
        )
    return contig


def parse_lea_cn_csv(lines: Sequence[str]) -> list[LeaChromosomeCopy]:
    """Parse and cross-check the historical ACE ``CN.csv`` output.

    The active R script writes ``Chromosome,Copies,Ploidy,CNA`` and defines
    ``CNA = Copies - Ploidy``. The adapter verifies that identity rather than trusting a
    plausible-looking CNA column.
    """
    rows = _reader(
        lines,
        required={"Chromosome", "Copies", "Ploidy", "CNA"},
        delimiter=",",
        artifact_name="CN.csv",
    )
    if not rows:
        raise LeaCompatibilityError(
            "CN.csv contains no chromosome rows; this is not interpretable as a negative call"
        )

    parsed: list[LeaChromosomeCopy] = []
    seen: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        contig = _canonical_human_contig(
            row["Chromosome"], artifact_name="CN.csv", line_number=line_number
        )
        if contig in seen:
            raise LeaCompatibilityError(f"CN.csv contains duplicate chromosome {contig}")
        seen.add(contig)
        copies = _finite_float(
            row["Copies"], artifact_name="CN.csv", line_number=line_number, field="Copies"
        )
        ploidy = _finite_float(
            row["Ploidy"], artifact_name="CN.csv", line_number=line_number, field="Ploidy"
        )
        cna = _finite_float(
            row["CNA"], artifact_name="CN.csv", line_number=line_number, field="CNA"
        )
        if copies < 0 or ploidy < 0:
            raise LeaCompatibilityError(
                f"CN.csv line {line_number} contains negative copy number or ploidy"
            )
        if not math.isclose(cna, copies - ploidy, abs_tol=1e-6):
            raise LeaCompatibilityError(
                f"CN.csv line {line_number} violates CNA = Copies - Ploidy"
            )
        parsed.append(LeaChromosomeCopy(contig, copies, ploidy, cna))
    return parsed


def parse_lea_dels_dups_csv(
    lines: Sequence[str],
    *,
    minimum_affected_fraction: float = 0.66,
) -> list[LeaBandCall]:
    """Parse the historical ``dels_dups.csv`` band-call output."""
    rows = _reader(
        lines,
        required={"chromosome", "name", "event", "frac_abr"},
        delimiter=",",
        artifact_name="dels_dups.csv",
    )
    parsed: list[LeaBandCall] = []
    seen_band: dict[tuple[str, str], str] = {}
    for line_number, row in enumerate(rows, start=2):
        contig = _canonical_human_contig(
            row["chromosome"], artifact_name="dels_dups.csv", line_number=line_number
        )
        band = row["name"]
        if not band:
            raise LeaCompatibilityError(
                f"dels_dups.csv line {line_number} contains an empty cytoband"
            )
        event = row["event"].lower()
        if event not in {"del", "dup"}:
            raise LeaCompatibilityError(
                f"dels_dups.csv line {line_number} has unsupported event {row['event']!r}"
            )
        fraction = _finite_float(
            row["frac_abr"],
            artifact_name="dels_dups.csv",
            line_number=line_number,
            field="frac_abr",
        )
        if not 0.0 <= fraction <= 1.0:
            raise LeaCompatibilityError(
                f"dels_dups.csv line {line_number} has frac_abr outside [0, 1]"
            )
        if fraction + 1e-12 < minimum_affected_fraction:
            raise LeaCompatibilityError(
                f"dels_dups.csv line {line_number} has frac_abr={fraction}, below the "
                f"declared historical threshold {minimum_affected_fraction}"
            )
        key = (contig, band)
        previous = seen_band.get(key)
        if previous is not None:
            if previous != event:
                raise LeaCompatibilityError(
                    f"dels_dups.csv assigns both deletion and duplication to {contig}{band}"
                )
            raise LeaCompatibilityError(
                f"dels_dups.csv contains duplicate band call {contig}{band} {event}"
            )
        seen_band[key] = event
        parsed.append(LeaBandCall(contig, band, event, fraction))
    return parsed


def _normalized_contig_lengths(contig_lengths: Mapping[str, int]) -> dict[str, int]:
    normalized: dict[str, int] = {}
    for raw, length in contig_lengths.items():
        contig = canonical_contig(raw)
        if contig not in CANONICAL_CONTIGS:
            continue
        if length <= 0:
            raise LeaCompatibilityError(f"reference contig {raw!r} has non-positive length")
        if contig in normalized and normalized[contig] != length:
            raise LeaCompatibilityError(
                f"reference supplies conflicting lengths for canonical contig {contig}"
            )
        normalized[contig] = int(length)
    if not AUTOSOMES.issubset(normalized):
        missing = sorted(AUTOSOMES - set(normalized), key=int)
        raise LeaCompatibilityError(
            "reference lengths are missing autosome(s): " + ", ".join(missing)
        )
    return normalized


def _contig_sort_key(contig: str) -> int:
    canonical = canonical_contig(contig)
    if canonical.isdigit():
        return int(canonical)
    return {"X": 23, "Y": 24}.get(canonical, 99)


def _directional_state(*, cna: float) -> CopyNumberState:
    if cna < 0:
        return CopyNumberState.LOSS
    if cna > 0:
        return CopyNumberState.GAIN
    return CopyNumberState.NEUTRAL


def _reject_segment_overlaps(segments: Sequence[CnvSegment]) -> None:
    by_contig: dict[str, list[CnvSegment]] = {}
    for segment in segments:
        by_contig.setdefault(canonical_contig(segment.contig), []).append(segment)
    for contig, items in by_contig.items():
        ordered = sorted(items, key=lambda item: (item.start, item.end))
        for previous, current in zip(ordered, ordered[1:], strict=False):
            if current.start < previous.end:
                raise LeaCompatibilityError(
                    f"normalized Lea calls overlap on chromosome {contig}: "
                    f"[{previous.start}, {previous.end}) and [{current.start}, {current.end})"
                )


def lea_ace_call_set_from_outputs(
    *,
    cn_lines: Sequence[str],
    dels_dups_lines: Sequence[str],
    cytobands: CytobandTable,
    contig_lengths: Mapping[str, int],
    call_set_id: str,
    sample_id: str,
    genome_build: GenomeBuild,
    data_basis: CnvDataBasis,
    profile: LeaAceCompatibilityProfile = LEA_ACE_2026_HG19,
    no_call_regions: Sequence[GenomicRegion] = (),
) -> CnvCallSet:
    """Normalize Lea's paired ACE outputs into a research-only :class:`CnvCallSet`.

    No lift-over is attempted. The historical profile, the locked cytoband table and the
    caller-supplied reference lengths must all agree on the genome build. Whole-chromosome
    CNA rows are converted to whole-contig intervals. Partial calls are mapped from the
    exact cytobands emitted in ``dels_dups.csv``. The latter file does not contain an
    absolute copy number, so partial calls retain only the historical gain/loss direction.
    """
    profile.validate()
    if genome_build != profile.genome_build:
        raise LeaCompatibilityError(
            f"requested build {genome_build.value} does not match historical profile "
            f"{profile.genome_build.value}; compatibility mode never performs lift-over"
        )
    if cytobands.genome_build != genome_build:
        raise LeaCompatibilityError(
            "cytoband table genome build does not match the historical call-set build"
        )
    lengths = _normalized_contig_lengths(contig_lengths)
    chromosome_rows = parse_lea_cn_csv(cn_lines)
    band_rows = parse_lea_dels_dups_csv(
        dels_dups_lines,
        minimum_affected_fraction=profile.affected_band_fraction,
    )

    chromosome_by_contig = {row.contig: row for row in chromosome_rows}
    segments: list[CnvSegment] = []
    warnings: list[str] = []

    missing_autosomes = sorted(AUTOSOMES - set(chromosome_by_contig), key=int)
    if missing_autosomes:
        raise LeaCompatibilityError(
            "historical CN.csv lacks autosome row(s): " + ", ".join(missing_autosomes)
            + "; the frozen whole-genome compatibility profile is incomplete"
        )
    missing_sex = sorted({"X", "Y"} - set(chromosome_by_contig))
    if missing_sex:
        warnings.append(
            "Historical CN.csv lacks sex-chromosome row(s): " + ", ".join(missing_sex) + "."
        )

    abnormal_whole_contigs: set[str] = set()
    for row in chromosome_rows:
        state = _directional_state(cna=row.cna)
        if state == CopyNumberState.NEUTRAL:
            continue
        length = lengths.get(row.contig)
        if length is None:
            raise LeaCompatibilityError(
                f"reference lengths do not contain chromosome {row.contig} from CN.csv"
            )
        abnormal_whole_contigs.add(row.contig)
        segments.append(
            CnvSegment(
                contig=row.contig,
                start=0,
                end=length,
                state=state,
                copy_number=row.copies,
                notes=[
                    "Historical Lea QDNAseq+ACE whole-chromosome call; direction derives "
                    "from CNA = Copies - Ploidy.",
                    f"Historical baseline ploidy for this chromosome: {row.ploidy:g}.",
                ],
            )
        )

    for row in band_rows:
        if row.contig not in chromosome_by_contig:
            raise LeaCompatibilityError(
                f"dels_dups.csv contains {row.contig}{row.band}, but CN.csv has no "
                f"chromosome {row.contig} row"
            )
        if row.contig in abnormal_whole_contigs:
            raise LeaCompatibilityError(
                f"dels_dups.csv contains a partial call on chromosome {row.contig}, while "
                "CN.csv already declares a whole-chromosome CNA; the historical ACE script "
                "suppresses this combination, so the artifacts are inconsistent"
            )
        try:
            start, end = cytobands.band_interval(row.contig, row.band)
        except (KeyError, ValueError) as error:
            raise LeaCompatibilityError(
                f"cannot map historical cytoband {row.contig}{row.band}: {error}"
            ) from error
        state = CopyNumberState.LOSS if row.event == "del" else CopyNumberState.GAIN
        segments.append(
            CnvSegment(
                contig=row.contig,
                start=start,
                end=end,
                state=state,
                copy_number=None,
                notes=[
                    "Historical Lea QDNAseq+ACE partial-band call.",
                    f"Affected fraction of cytoband: {row.affected_fraction:.6g}.",
                    "Absolute copy number is unavailable in dels_dups.csv and is not "
                    "invented by the compatibility adapter.",
                ],
            )
        )

    segments.sort(key=lambda item: (_contig_sort_key(item.contig), item.start, item.end))
    _reject_segment_overlaps(segments)

    if not profile.runtime_versions_verified:
        warnings.append(
            "The supplied renv.lock declares QDNAseq/ACE package versions, but the supplied "
            "r-ace container installs Bioconductor packages without exact version pins. "
            "Runtime package versions therefore remain unverified provenance."
        )

    reports_negative = not segments
    return CnvCallSet(
        call_set_id=call_set_id,
        sample_id=sample_id,
        genome_build=genome_build,
        method="lea-qdnaseq-ace-historical-comparator",
        method_version=profile.profile_id,
        data_basis=data_basis,
        background_state=CopyNumberState.NEUTRAL,
        status=ModuleRunStatus.COMPLETED,
        segments=segments,
        no_call_regions=list(no_call_regions),
        reports_biological_negative=reports_negative,
        bin_size_bp=profile.qdnaseq_bin_size_kbp * 1000,
        tool=ToolRecord(
            name="Lea Evers historical QDNAseq+ACE profile",
            version=profile.profile_id,
            parameters={
                "genome_build": genome_build.value,
                "qdnaseq_bin_size_kbp": profile.qdnaseq_bin_size_kbp,
                "ace_penalty": profile.ace_penalty,
                "assumed_autosomal_ploidy": profile.assumed_autosomal_ploidy,
                "affected_band_fraction": profile.affected_band_fraction,
                "declared_qdnaseq_version": profile.declared_qdnaseq_version,
                "declared_qdnaseq_hg19_version": profile.declared_qdnaseq_hg19_version,
                "declared_ace_version": profile.declared_ace_version,
                "runtime_versions_verified": profile.runtime_versions_verified,
            },
        ),
        warnings=warnings,
        limitations=[
            "Historical comparator only. The adapter reproduces output semantics; no "
            "source code from the thesis project is copied or executed here.",
            "The frozen profile is GRCh37/hg19 only. No coordinate lift-over is performed.",
            "Partial-band calls are thresholded historical outputs (frac_abr >= 0.66) and "
            "carry gain/loss direction but no absolute copy number.",
            "A COMPLETED empty call set means the two historical output artifacts were "
            "structurally valid and contained no alteration under that historical method; "
            "it is not an assay-wide or clinical negative result.",
            "Observability and coverage exclusions must be supplied separately as no-call "
            "regions before performance metrics are interpreted.",
            "Research-only and non-reportable until assay-specific analytical validation.",
        ],
    )
