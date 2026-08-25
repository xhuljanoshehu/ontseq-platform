from __future__ import annotations

import re
from dataclasses import dataclass

from .iscn_reference import CytobandIndex
from .iscn_validation import ISCNValidationStatus, validate_iscn
from .models import EventType, GenomicEvent, ISCNProposal, Locus


def _chromosome_number(chromosome: str) -> int:
    value = chromosome.removeprefix("chr")
    return {"X": 23, "Y": 24}.get(value, int(value) if value.isdigit() else 99)


def _event_order(event: GenomicEvent) -> tuple[int, int, int, str]:
    numerical = event.event_type in {EventType.CHROMOSOME_GAIN, EventType.CHROMOSOME_LOSS}
    return (
        _chromosome_number(event.primary.chromosome),
        0 if numerical else 1,
        event.primary.start,
        event.event_type.value,
    )


@dataclass(frozen=True, slots=True)
class _ResolvedBands:
    start: str
    end: str
    derived_from_coordinates: bool


def _resolve_locus_bands(locus: Locus, cytobands: CytobandIndex | None) -> _ResolvedBands | None:
    if locus.cytoband_start:
        return _ResolvedBands(
            start=locus.cytoband_start,
            end=locus.cytoband_end or locus.cytoband_start,
            derived_from_coordinates=False,
        )
    if cytobands is None:
        return None
    mapped = cytobands.bands_for_interval(locus.chromosome, locus.start, locus.end)
    if mapped is None:
        return None
    first, last = mapped
    return _ResolvedBands(
        start=first.name,
        end=last.name,
        derived_from_coordinates=True,
    )


def _band_range(
    event: GenomicEvent,
    cytobands: CytobandIndex | None,
) -> str | None:
    bands = _resolve_locus_bands(event.primary, cytobands)
    if bands is None:
        return None

    start = bands.start
    end = bands.end
    if event.event_type in {EventType.DELETION, EventType.DUPLICATION}:
        if start[0] != end[0]:
            return None
        # Genomic coordinates run telomere -> centromere on the p arm. For coordinate-derived
        # deletion/duplication intervals, render the cytobands in the conventional
        # centromere-to-telomere order used by the legacy cytoband-merging implementation.
        if bands.derived_from_coordinates and start[0] == "p" and start != end:
            start, end = end, start

    return start if start == end else f"{start}{end}"


def _breakpoint_band(locus: Locus, cytobands: CytobandIndex | None) -> str | None:
    if locus.cytoband_start:
        return locus.cytoband_start
    if cytobands is None:
        return None
    band = cytobands.band_at(locus.chromosome, locus.start)
    return band.name if band is not None else None


def _fragment(event: GenomicEvent, cytobands: CytobandIndex | None) -> str | None:
    chrom = event.primary.chromosome.removeprefix("chr")
    if event.event_type == EventType.CHROMOSOME_GAIN:
        return f"+{chrom}"
    if event.event_type == EventType.CHROMOSOME_LOSS:
        return f"-{chrom}"

    band = _band_range(event, cytobands)
    if event.event_type == EventType.DELETION and band:
        return f"del({chrom})({band})"
    if event.event_type == EventType.DUPLICATION and band:
        return f"dup({chrom})({band})"
    if event.event_type == EventType.INVERSION and band:
        return f"inv({chrom})({band})"
    if event.event_type == EventType.TRANSLOCATION and event.secondary:
        other = event.secondary.chromosome.removeprefix("chr")
        band_a = _breakpoint_band(event.primary, cytobands)
        band_b = _breakpoint_band(event.secondary, cytobands)
        if band_a and band_b:
            return f"t({chrom};{other})({band_a};{band_b})"
    return None


def _derived_chromosome_count(events: list[GenomicEvent]) -> int:
    count = 46
    for event in events:
        if not event.reportable:
            continue
        if event.event_type == EventType.CHROMOSOME_GAIN:
            count += 1
        elif event.event_type == EventType.CHROMOSOME_LOSS:
            count -= 1
    return count


def build_iscn_proposal(
    events: list[GenomicEvent],
    chromosome_count: int | None = None,
    sex_chromosomes: str = "XX",
    *,
    cytobands: CytobandIndex | None = None,
    validate: bool = True,
    prefer_external_validator: bool = True,
) -> ISCNProposal:
    """Render the auditable ONTSeq ISCN subset from normalized genomic events.

    The renderer keeps the genomic event model as the source of truth and treats ISCN as an
    output representation. When a cytoband index is supplied, genomic coordinates are mapped
    to bands at render time. Only chromosome gains/losses, del, dup, inv and reciprocal-style
    translocation notation are emitted. Unsupported events are never guessed; they are
    omitted with an explicit warning.

    The returned object remains a proposal requiring expert cytogenetic review. Validation is
    a technical gate, not authorization for automatic clinical release.
    """

    if not re.fullmatch(r"(?:XX|XY|X|XXY|XYY)", sex_chromosomes):
        raise ValueError("Unsupported sex chromosome complement for subset renderer")
    if chromosome_count is not None and not 1 <= chromosome_count <= 999:
        raise ValueError("chromosome_count must be between 1 and 999")

    ordered = sorted(events, key=_event_order)
    fragments: list[str] = []
    event_ids: list[str] = []
    warnings: list[str] = [
        "Automatically generated ISCN proposal; expert cytogenetic review is mandatory."
    ]

    for event in ordered:
        if not event.reportable:
            continue
        fragment = _fragment(event, cytobands)
        if fragment:
            fragments.append(fragment)
            event_ids.append(event.event_id)
        else:
            warnings.append(
                f"Event {event.event_id} is outside the implemented ISCN subset, lacks "
                "mappable cytobands, or cannot be represented safely and was omitted."
            )

    resolved_count = chromosome_count if chromosome_count is not None else _derived_chromosome_count(events)
    base = f"{resolved_count},{sex_chromosomes}"
    notation = base if not fragments else f"{base},{','.join(fragments)}"

    if validate:
        validation = validate_iscn(notation, prefer_external=prefer_external_validator)
        warnings.append(f"ISCN validation {validation.status.value} via {validation.engine}.")
        warnings.extend(validation.messages)
        if validation.status == ISCNValidationStatus.FAIL:
            warnings.append("Validation failed; the proposal must not be released as valid ISCN.")
    else:
        warnings.append("ISCN validation NOT_RUN by caller request.")

    return ISCNProposal(notation=notation, source_event_ids=event_ids, warnings=warnings)
