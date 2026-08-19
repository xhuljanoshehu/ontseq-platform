"""Build-aware normalization of adjacent cytoband events into reviewable regions.

This module is an independent implementation of a useful behavior found in the historical
ONTseq workflow: adjacent affected cytobands can be presented as one region, and a complete
set of sub-bands can sometimes be represented by their immediate parent designation.  The
output here is deliberately **not ISCN**.  It is a structured intermediate representation
with source-band traceability that can later be reviewed by the ISCN proposal layer.

Coordinates and adjacency always come from a locked :class:`CytobandTable`.  Band-name
arithmetic is used only for one-level parent-label compaction after coordinate coverage has
been proved against that reference.  Restricting one normalization pass to the immediate
parent avoids silently over-collapsing deep labels such as ``q31.21`` straight to ``q31``.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from .cytobands import Cytoband, CytobandTable
from .intervals import canonical_contig
from .states import CopyNumberState

_BAND_RE = re.compile(r"^(?P<arm>[pq])(?P<integer>\d+)(?:\.(?P<decimal>\d+))?$")
_ALLOWED_STATES = {
    CopyNumberState.LOSS,
    CopyNumberState.GAIN,
    CopyNumberState.HOMOZYGOUS_LOSS,
    CopyNumberState.HIGH_AMPLIFICATION,
}


class CytobandRegionError(ValueError):
    """Raised when band events cannot be normalized without ambiguity."""


@dataclass(frozen=True)
class CytobandEvent:
    """One source cytoband-level dosage event.

    ``source_id`` is an optional non-identifying event identifier used only for traceability.
    It must not contain patient identifiers or filesystem paths.
    """

    contig: str
    band: str
    state: CopyNumberState
    source_id: str | None = None


@dataclass(frozen=True)
class CytobandRegion:
    """One contiguous, same-arm, same-state cytoband region.

    ``start``/``end`` are zero-based half-open genomic coordinates.  ``display_first_band``
    and ``display_last_band`` are compacted cytogenetic labels for presentation.  They are
    not asserted to be an ISCN expression.
    """

    contig: str
    arm: str
    state: CopyNumberState
    start: int
    end: int
    display_first_band: str
    display_last_band: str
    source_bands: tuple[str, ...]
    source_ids: tuple[str, ...]

    @property
    def is_single_band_label(self) -> bool:
        return self.display_first_band == self.display_last_band


@dataclass(frozen=True)
class _ResolvedEvent:
    contig: str
    arm: str
    state: CopyNumberState
    start: int
    end: int
    source_band: str
    source_id: str | None


def _parse_band_name(name: str) -> tuple[str, str, str | None]:
    match = _BAND_RE.fullmatch(name)
    if match is None:
        raise CytobandRegionError(f"unsupported cytoband designation: {name!r}")
    return match.group("arm"), match.group("integer"), match.group("decimal")


def _parent_band(name: str) -> str | None:
    """Return the immediate useful parent label, stopping at the integer band.

    Examples: ``q31.21 -> q31.2`` and ``p12.3 -> p12``.  ``q31`` has no parent here:
    collapsing it to ``q3`` would change a band into a region-level shorthand that is not
    represented by the UCSC cytoband resource.
    """
    arm, integer, decimal = _parse_band_name(name)
    if decimal is None:
        return None
    if len(decimal) == 1:
        return f"{arm}{integer}"
    return f"{arm}{integer}.{decimal[:-1]}"


def _leaf_bands_under(table: CytobandTable, contig: str, designation: str) -> list[Cytoband]:
    """Return reference leaf bands represented by ``designation``.

    UCSC tables normally contain terminal band labels only.  A parent such as ``q13`` is
    represented by every exact/sub-band row whose name is ``q13`` or begins ``q13.``.
    Contig comparison is canonicalized because reference tables may use ``chr5`` while the
    normalized event contract uses ``5``.
    """
    _parse_band_name(designation)
    canonical = canonical_contig(contig)
    return [
        band
        for band in table.bands
        if canonical_contig(band.contig) == canonical
        and (band.name == designation or band.name.startswith(f"{designation}."))
    ]


def _interval_fully_covered(start: int, end: int, intervals: Sequence[tuple[int, int]]) -> bool:
    cursor = start
    for left, right in sorted(intervals):
        if right <= cursor:
            continue
        if left > cursor:
            return False
        cursor = max(cursor, right)
        if cursor >= end:
            return True
    return cursor >= end


def _compact_endpoint(
    *,
    table: CytobandTable,
    contig: str,
    source_band: str,
    affected_intervals: Sequence[tuple[int, int]],
) -> str:
    """Compact to the immediate parent only when all of its reference leaves are covered."""
    parent = _parent_band(source_band)
    if parent is None:
        return source_band
    children = _leaf_bands_under(table, contig, parent)
    if not children:
        return source_band
    if all(
        _interval_fully_covered(child.start, child.end, affected_intervals) for child in children
    ):
        return parent
    return source_band


def _resolve_event(event: CytobandEvent, table: CytobandTable) -> _ResolvedEvent:
    contig = canonical_contig(event.contig)
    if event.state not in _ALLOWED_STATES:
        raise CytobandRegionError(
            f"cytoband region normalization requires an abnormal dosage state, got {event.state}"
        )
    arm, _, _ = _parse_band_name(event.band)
    try:
        start, end = table.band_interval(contig, event.band)
    except (KeyError, ValueError) as error:
        raise CytobandRegionError(
            f"cannot resolve cytoband {contig}{event.band} in resource {table.resource_id}: {error}"
        ) from error
    if event.source_id is not None:
        if not event.source_id.strip():
            raise CytobandRegionError("source_id must be non-empty when supplied")
        if "/" in event.source_id or "\\" in event.source_id:
            raise CytobandRegionError("source_id must not contain a filesystem path")
    return _ResolvedEvent(
        contig=contig,
        arm=arm,
        state=event.state,
        start=start,
        end=end,
        source_band=event.band,
        source_id=event.source_id,
    )


def _contig_sort_key(contig: str) -> tuple[int, str]:
    canonical = canonical_contig(contig)
    if canonical.isdigit():
        return int(canonical), canonical
    return {"X": 23, "Y": 24}.get(canonical, 99), canonical


def _validate_no_overlap(events: Sequence[_ResolvedEvent]) -> None:
    ordered = sorted(
        events,
        key=lambda item: (_contig_sort_key(item.contig), item.start, item.end, item.state.value),
    )
    for previous, current in zip(ordered, ordered[1:], strict=False):
        if previous.contig != current.contig:
            continue
        if current.start < previous.end:
            raise CytobandRegionError(
                "source cytoband events overlap after reference resolution: "
                f"{previous.contig}{previous.source_band} and "
                f"{current.contig}{current.source_band}"
            )


def _are_contiguous(left: _ResolvedEvent, right: _ResolvedEvent) -> bool:
    return (
        left.contig == right.contig
        and left.arm == right.arm
        and left.state == right.state
        and left.end == right.start
    )


def _presentation_endpoints(
    *,
    table: CytobandTable,
    group: Sequence[_ResolvedEvent],
) -> tuple[str, str]:
    intervals = [(item.start, item.end) for item in group]
    genomic_first = _compact_endpoint(
        table=table,
        contig=group[0].contig,
        source_band=group[0].source_band,
        affected_intervals=intervals,
    )
    genomic_last = _compact_endpoint(
        table=table,
        contig=group[-1].contig,
        source_band=group[-1].source_band,
        affected_intervals=intervals,
    )

    # p-arm reference order runs telomere -> centromere.  Cytogenetic interval notation
    # presents the more proximal/lower-numbered boundary first.  Keep genomic coordinates
    # separately so presentation order can never alter interval arithmetic.
    if group[0].arm == "p" and genomic_first != genomic_last:
        return genomic_last, genomic_first
    return genomic_first, genomic_last


def normalize_cytoband_regions(
    events: Iterable[CytobandEvent],
    table: CytobandTable,
) -> list[CytobandRegion]:
    """Normalize cytoband events into deterministic contiguous regions.

    Rules:

    * every source band must resolve against the supplied build-locked resource;
    * overlapping parent/child or duplicate designations fail closed;
    * only directly coordinate-adjacent bands merge;
    * state, chromosome and p/q arm must match, so the centromere is never crossed;
    * one normalization pass may compact an endpoint only to its immediate parent, and only
      when all reference leaf bands underneath that parent are covered by the affected region;
    * source labels and non-identifying source IDs are retained for auditability.
    """
    resolved = [_resolve_event(event, table) for event in events]
    if not resolved:
        return []
    _validate_no_overlap(resolved)
    resolved.sort(key=lambda item: (_contig_sort_key(item.contig), item.start, item.end))

    groups: list[list[_ResolvedEvent]] = []
    current: list[_ResolvedEvent] = []
    for event in resolved:
        if not current or _are_contiguous(current[-1], event):
            current.append(event)
            continue
        groups.append(current)
        current = [event]
    if current:
        groups.append(current)

    regions: list[CytobandRegion] = []
    for group in groups:
        first_label, last_label = _presentation_endpoints(table=table, group=group)
        regions.append(
            CytobandRegion(
                contig=group[0].contig,
                arm=group[0].arm,
                state=group[0].state,
                start=group[0].start,
                end=group[-1].end,
                display_first_band=first_label,
                display_last_band=last_label,
                source_bands=tuple(item.source_band for item in group),
                source_ids=tuple(item.source_id for item in group if item.source_id is not None),
            )
        )
    return regions
