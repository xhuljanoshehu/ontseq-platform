"""Multi-source CNV truth construction.

The project has four realistic truth sources, and they differ in ways that matter more
than their formats:

============================  ==================  ==================  ==================
Source                        Breakpoint          Background          Detects CN-LOH
============================  ==================  ==================  ==================
Cytogenetics (ISCN)           band width (Mb)     closed, ~5-10 Mb    no
FISH                          probe locus         open, probes only   no
SNP array                     probe spacing (kb)  closed, probe map    yes
CGH array                     probe spacing (kb)  closed, probe map    no
Short-read WGS/panel          base pair           closed / ROI        with allele data
Simulation                    exact               closed              by construction
============================  ==================  ==================  ==================

Encoding these differences is the whole point of this module. A karyotype-derived truth
set that claims base-pair breakpoints will make every caller look inaccurate; a FISH
truth set that claims a genome-wide neutral background will make every caller look
non-specific. Both failures are silent unless the truth carries its own metadata.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field

from ..models import GenomeBuild
from .cytobands import CytobandTable
from .intervals import canonical_contig
from .models import (
    CnvSegment,
    CnvTruthSet,
    CnvTruthSource,
    GenomicRegion,
)
from .states import CopyNumberState, state_from_copy_number

#: ISCN constructs that are balanced and therefore assert no copy-number change.
BALANCED_CONSTRUCTS = ("t(", "inv(", "ins(")

_CLONE = re.compile(r"^(?P<body>.*?)(?:\[(?P<cells>\d+)\])?$")
_WHOLE_CHROMOSOME = re.compile(r"^(?P<sign>[+-])(?P<contig>\d{1,2}|X|Y)$")
_STRUCTURAL = re.compile(
    r"^(?P<kind>del|dup)\((?P<contig>\d{1,2}|X|Y)\)"
    r"\((?P<bands>[pq0-9.]+)\)$"
)
_ISOCHROMOSOME = re.compile(r"^i\((?P<contig>\d{1,2}|X|Y)\)\((?P<arm>[pq])10\)$")
_BAND_TOKEN = re.compile(r"[pq]\d+(?:\.\d+)?")


@dataclass(frozen=True)
class UnsupportedConstruct:
    """An ISCN token the converter deliberately refused to interpret."""

    token: str
    reason: str


@dataclass
class KaryotypeConversion:
    """Result of converting an ISCN karyotype into copy-number intervals."""

    segments: list[CnvSegment] = field(default_factory=list)
    unsupported: list[UnsupportedConstruct] = field(default_factory=list)
    balanced_constructs: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    clone_cell_counts: list[int] = field(default_factory=list)

    @property
    def is_complete(self) -> bool:
        """Whether every token was either converted or recognised as balanced."""
        return not self.unsupported


def _copy_numbers(state: CopyNumberState, baseline_ploidy: float) -> float:
    if state == CopyNumberState.LOSS:
        return max(0.0, baseline_ploidy - 1.0)
    if state == CopyNumberState.GAIN:
        return baseline_ploidy + 1.0
    return baseline_ploidy


def _split_clones(karyotype: str) -> list[tuple[str, int | None]]:
    """Split a karyotype into clones with their observed cell counts."""
    clones: list[tuple[str, int | None]] = []
    for raw in karyotype.split("/"):
        text = raw.strip()
        if not text:
            continue
        match = _CLONE.match(text)
        if match is None:  # pragma: no cover - the pattern always matches
            clones.append((text, None))
            continue
        cells = match.group("cells")
        clones.append((match.group("body").strip(), int(cells) if cells else None))
    return clones


def _tokenize_clone(clone: str) -> list[str]:
    """Split one clone into ISCN tokens, respecting parentheses.

    A naive split on commas would break ``del(5)(q13q33)`` apart at nothing, but ISCN
    also contains constructs such as ``t(9;22)(q34;q11)`` where commas can appear inside
    a bracketed group in derivative notation. Tracking depth keeps those intact.
    """
    tokens: list[str] = []
    current: list[str] = []
    depth = 0
    for character in clone:
        if character == "(":
            depth += 1
        elif character == ")":
            depth = max(0, depth - 1)
        if character == "," and depth == 0:
            tokens.append("".join(current).strip())
            current = []
            continue
        current.append(character)
    if current:
        tokens.append("".join(current).strip())
    return [token for token in tokens if token]


def convert_karyotype(
    karyotype: str,
    cytobands: CytobandTable,
    *,
    baseline_ploidy: float = 2.0,
) -> KaryotypeConversion:
    """Convert an ISCN karyotype string into copy-number segments with band uncertainty.

    The converter supports whole-chromosome gains and losses, interstitial and terminal
    ``del``/``dup``, and isochromosomes. Everything else is recorded in
    :attr:`KaryotypeConversion.unsupported` with a reason and is never silently dropped,
    because a truth set that quietly lost half its findings is worse than no truth set:
    it converts real events into apparent false positives.

    Balanced rearrangements are recognised and recorded separately. They assert no
    copy-number change, so they contribute no segments, but their presence in the
    karyotype is retained for provenance.

    This is a research converter. It does not implement ISCN 2024 and is not a substitute
    for expert cytogenetic review; see ``docs/ARCHITECTURE.md`` on the ISCN boundary.
    """
    conversion = KaryotypeConversion()
    clones = _split_clones(karyotype)
    if not clones:
        conversion.warnings.append("The karyotype string contained no clone.")
        return conversion
    if len(clones) > 1:
        conversion.warnings.append(
            f"The karyotype describes {len(clones)} clones. Segments from all clones are "
            "merged into one truth profile, so subclonal structure is lost. Clone-aware "
            "truth requires a per-clone truth set."
        )

    seen: set[tuple[str, int, int, CopyNumberState]] = set()
    for clone, cells in clones:
        if cells is not None:
            conversion.clone_cell_counts.append(cells)
        tokens = _tokenize_clone(clone)
        for position, token in enumerate(tokens):
            # The first two fields are the modal number and the sex complement.
            if position == 0 and re.fullmatch(r"\d{2,3}(?:[-~]\d{2,3})?", token):
                continue
            if position == 1 and re.fullmatch(r"[XY]{0,4}", token):
                if token not in {"XX", "XY"}:
                    _record_unsupported(
                        conversion,
                        token,
                        "sex-chromosome complement differs from XX/XY; sex chromosome "
                        "dosage requires the constitutional complement, which the "
                        "karyotype alone does not supply",
                    )
                continue
            if token in {"", "?"}:
                continue
            if token.startswith(BALANCED_CONSTRUCTS):
                conversion.balanced_constructs.append(token)
                continue
            if "?" in token:
                _record_unsupported(
                    conversion,
                    token,
                    "contains an uncertainty marker, so no interval is defensible",
                )
                continue

            for segment in _convert_token(token, cytobands, baseline_ploidy, conversion):
                key = (segment.contig, segment.start, segment.end, segment.state)
                if key in seen:
                    # The same alteration in two clones is one genomic claim.
                    continue
                seen.add(key)
                conversion.segments.append(segment)

    conversion.segments.sort(key=lambda item: (item.contig, item.start, item.end))
    return conversion


def _record_unsupported(conversion: KaryotypeConversion, token: str, reason: str) -> None:
    """Record a refused construct once per distinct token and reason."""
    if any(item.token == token and item.reason == reason for item in conversion.unsupported):
        return
    conversion.unsupported.append(UnsupportedConstruct(token=token, reason=reason))


def _convert_token(
    token: str,
    cytobands: CytobandTable,
    baseline_ploidy: float,
    conversion: KaryotypeConversion,
) -> list[CnvSegment]:
    """Convert one ISCN token into zero or more copy-number segments.

    Returns a list because a single token can assert two changes: an isochromosome
    simultaneously loses one arm and gains the other.
    """
    whole = _WHOLE_CHROMOSOME.fullmatch(token)
    if whole:
        contig = canonical_contig(whole.group("contig"))
        try:
            length = cytobands.contig_length(contig)
        except KeyError:
            _record_unsupported(
                conversion, token, f"contig {contig} is absent from the cytoband table"
            )
            return []
        state = CopyNumberState.GAIN if whole.group("sign") == "+" else CopyNumberState.LOSS
        return [
            CnvSegment(
                contig=contig,
                start=0,
                end=length,
                state=state,
                copy_number=_copy_numbers(state, baseline_ploidy),
                # A whole-chromosome change has no intra-chromosomal breakpoint.
                start_uncertainty_bp=0,
                end_uncertainty_bp=0,
                notes=[f"Derived from ISCN token {token!r}."],
            )
        ]

    isochromosome = _ISOCHROMOSOME.fullmatch(token)
    if isochromosome:
        contig = canonical_contig(isochromosome.group("contig"))
        retained = isochromosome.group("arm")
        lost = "p" if retained == "q" else "q"
        try:
            gained_span = cytobands.arm_interval(contig, retained)
            lost_span = cytobands.arm_interval(contig, lost)
        except KeyError as error:
            _record_unsupported(conversion, token, str(error))
            return []
        return [
            CnvSegment(
                contig=contig,
                start=lost_span[0],
                end=lost_span[1],
                state=CopyNumberState.LOSS,
                copy_number=_copy_numbers(CopyNumberState.LOSS, baseline_ploidy),
                notes=[f"Lost arm implied by ISCN token {token!r}."],
            ),
            CnvSegment(
                contig=contig,
                start=gained_span[0],
                end=gained_span[1],
                state=CopyNumberState.GAIN,
                copy_number=_copy_numbers(CopyNumberState.GAIN, baseline_ploidy),
                notes=[f"Gained arm implied by ISCN token {token!r}."],
            ),
        ]

    structural = _STRUCTURAL.fullmatch(token)
    if structural:
        contig = canonical_contig(structural.group("contig"))
        bands = _BAND_TOKEN.findall(structural.group("bands"))
        if not bands:
            _record_unsupported(conversion, token, "no interpretable band designation")
            return []
        state = CopyNumberState.LOSS if structural.group("kind") == "del" else CopyNumberState.GAIN
        try:
            if len(bands) == 1:
                # One breakpoint denotes a terminal event running to the arm end.
                band_start, band_end = cytobands.band_interval(contig, bands[0])
                arm_start, arm_end = cytobands.arm_interval(contig, bands[0][0])
                start, end = (band_start, arm_end) if bands[0][0] == "q" else (arm_start, band_end)
                uncertainty = cytobands.band_uncertainty(contig, bands[0])
                start_uncertainty = uncertainty if bands[0][0] == "q" else 0
                end_uncertainty = 0 if bands[0][0] == "q" else uncertainty
            else:
                start, end = cytobands.band_span(contig, bands[0], bands[1])
                start_uncertainty = cytobands.band_uncertainty(contig, bands[0])
                end_uncertainty = cytobands.band_uncertainty(contig, bands[1])
        except (KeyError, ValueError) as error:
            _record_unsupported(conversion, token, str(error))
            return []
        if end <= start:
            _record_unsupported(conversion, token, "band designations produced an empty interval")
            return []
        return [
            CnvSegment(
                contig=contig,
                start=start,
                end=end,
                state=state,
                copy_number=_copy_numbers(state, baseline_ploidy),
                start_uncertainty_bp=start_uncertainty,
                end_uncertainty_bp=end_uncertainty,
                notes=[
                    f"Derived from ISCN token {token!r}. Breakpoints are band-resolved, "
                    "not base-pair resolved."
                ],
            )
        ]

    _record_unsupported(
        conversion,
        token,
        "construct is not supported by the copy-number converter; derivative and marker "
        "chromosomes need expert interpretation before their dosage effect is defined",
    )
    return []


def truth_from_karyotype(
    *,
    truth_id: str,
    sample_id: str,
    karyotype: str,
    cytobands: CytobandTable,
    source_version: str,
    resolution_bp: int = 10_000_000,
    tumor_fraction: float | None = None,
    baseline_ploidy: float = 2.0,
) -> tuple[CnvTruthSet, KaryotypeConversion]:
    """Build a truth set from an ISCN karyotype.

    ``resolution_bp`` defaults to 10 Mb, the conventional approximate resolution of
    metaphase banding at 400-550 band level. Events smaller than this are invisible to
    karyotyping, so the truth set is silent about them and a caller must not be charged
    with a false positive for reporting one.

    The truth set is marked closed-world: a karyotype does assert that the rest of the
    genome is unremarkable, but only at its own resolution. That is why the resolution
    must be recorded alongside the background state - the pair is what makes the claim
    interpretable.
    """
    conversion = convert_karyotype(karyotype, cytobands, baseline_ploidy=baseline_ploidy)
    limitations = [
        "Breakpoints are resolved to cytogenetic bands, not base pairs. Breakpoint "
        "accuracy against this truth set is not meaningful below the band width.",
        f"Events smaller than the declared resolution of {resolution_bp} bp are invisible "
        "to karyotyping; this truth set is silent, not negative, below that size.",
        "Copy-neutral loss of heterozygosity is undetectable by banding and is absent "
        "from this truth set by construction.",
        "Clone structure is flattened into a single profile.",
    ]
    if conversion.unsupported:
        limitations.append(
            "The karyotype contained "
            f"{len(conversion.unsupported)} construct(s) the converter refused to "
            "interpret; they are absent from the segments and must be reconciled "
            "manually before this truth set is used for scoring."
        )
    if conversion.balanced_constructs:
        limitations.append(
            "Balanced rearrangements present in the karyotype "
            f"({', '.join(conversion.balanced_constructs)}) assert no copy-number change "
            "and contribute no segments."
        )

    truth = CnvTruthSet(
        truth_id=truth_id,
        sample_id=sample_id,
        genome_build=cytobands.genome_build,
        source=CnvTruthSource.ISCN_KARYOTYPE,
        source_version=source_version,
        background_state=CopyNumberState.NEUTRAL,
        resolution_bp=resolution_bp,
        segments=conversion.segments,
        tumor_fraction=tumor_fraction,
        baseline_ploidy=baseline_ploidy,
        limitations=limitations,
    )
    return truth, conversion


def truth_from_segments(
    *,
    truth_id: str,
    sample_id: str,
    genome_build: GenomeBuild,
    source: CnvTruthSource,
    source_version: str,
    segments: Sequence[CnvSegment],
    resolution_bp: int,
    informative_regions: Sequence[GenomicRegion] = (),
    uninformative_regions: Sequence[GenomicRegion] = (),
    tumor_fraction: float | None = None,
    baseline_ploidy: float = 2.0,
    closed_world: bool = True,
    extra_limitations: Sequence[str] = (),
) -> CnvTruthSet:
    """Build a truth set from an already-segmented source such as an array or WGS run."""
    limitations = [
        f"Truth derived from {source.value} at version {source_version}.",
        f"Events below the declared resolution of {resolution_bp} bp are not asserted.",
        *extra_limitations,
    ]
    if source in {CnvTruthSource.CGH_ARRAY, CnvTruthSource.SHORT_READ_WGS} and any(
        segment.state == CopyNumberState.COPY_NEUTRAL_LOH for segment in segments
    ):
        limitations.append(
            "Copy-neutral LOH segments are present but the declared source does not "
            "normally resolve allelic imbalance; confirm the provenance of those calls."
        )
    return CnvTruthSet(
        truth_id=truth_id,
        sample_id=sample_id,
        genome_build=genome_build,
        source=source,
        source_version=source_version,
        background_state=(CopyNumberState.NEUTRAL if closed_world else CopyNumberState.NO_CALL),
        resolution_bp=resolution_bp,
        segments=list(segments),
        informative_regions=list(informative_regions),
        uninformative_regions=list(uninformative_regions),
        tumor_fraction=tumor_fraction,
        baseline_ploidy=baseline_ploidy,
        limitations=limitations,
    )


def truth_from_fish(
    *,
    truth_id: str,
    sample_id: str,
    genome_build: GenomeBuild,
    source_version: str,
    probes: Sequence[tuple[str, int, int, str, float]],
    baseline_ploidy: float = 2.0,
    tumor_fraction: float | None = None,
) -> CnvTruthSet:
    """Build an open-world truth set from FISH probe results.

    ``probes`` are ``(contig, start, end, label, copy_number)`` tuples.

    FISH interrogates named loci and says nothing whatsoever about the rest of the
    genome, so the background is ``NO_CALL`` and the probe footprints become the
    informative regions. Any other encoding would turn every genuine finding outside the
    probe set into a false positive.
    """
    segments = [
        CnvSegment(
            contig=canonical_contig(contig),
            start=start,
            end=end,
            state=state_from_copy_number(copy_number, baseline_ploidy=baseline_ploidy),
            copy_number=copy_number,
            notes=[f"FISH probe {label}."],
        )
        for contig, start, end, label, copy_number in probes
    ]
    regions = [
        GenomicRegion(contig=canonical_contig(contig), start=start, end=end, label=label)
        for contig, start, end, label, _ in probes
    ]
    return CnvTruthSet(
        truth_id=truth_id,
        sample_id=sample_id,
        genome_build=genome_build,
        source=CnvTruthSource.FISH,
        source_version=source_version,
        background_state=CopyNumberState.NO_CALL,
        resolution_bp=0,
        segments=segments,
        informative_regions=regions,
        tumor_fraction=tumor_fraction,
        baseline_ploidy=baseline_ploidy,
        limitations=[
            "FISH asserts copy number only at the interrogated probe loci. The rest of "
            "the genome is not negative, it is unexamined.",
            "Probe footprints are an approximation of the interrogated locus and do not "
            "define event breakpoints.",
        ],
    )
