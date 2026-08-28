"""Conservative fusion-evidence projection from independently annotated breakpoints."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Literal

from .breakpoint_annotation import (
    AnnotatedBreakpoint,
    AnnotatedBreakpointPair,
    TranscriptBreakpointHit,
)


class FusionFrameStatus(StrEnum):
    IN_FRAME = "in_frame"
    OUT_OF_FRAME = "out_of_frame"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class FusionPartnerEvidence:
    gene: str | None
    preferred_transcript: str | None
    region: Literal["exon", "intron", "transcript", "unknown"]
    exon_number: int | None
    intron_number: int | None
    strand: str | None


@dataclass(frozen=True)
class FusionEvidence:
    gene_a: FusionPartnerEvidence
    gene_b: FusionPartnerEvidence
    orientation: str | None
    frame_status: FusionFrameStatus = FusionFrameStatus.UNKNOWN


def normalize_orientation(value: str | None) -> str | None:
    """Return only an unambiguous two-breakend strand representation."""
    if value is None:
        return None
    compact = value.strip().replace("/", "").replace(" ", "")
    return compact if compact in {"++", "+-", "-+", "--"} else None


def _preferred_hit(breakpoint: AnnotatedBreakpoint | None) -> TranscriptBreakpointHit | None:
    if breakpoint is None or not breakpoint.transcripts:
        return None
    preferred = [hit for hit in breakpoint.transcripts if hit.preferred]
    candidates = preferred or list(breakpoint.transcripts)
    genes = {hit.gene_name for hit in candidates}
    if len(genes) != 1:
        return None
    return candidates[0]


def _partner(breakpoint: AnnotatedBreakpoint | None) -> FusionPartnerEvidence:
    hit = _preferred_hit(breakpoint)
    if hit is None:
        return FusionPartnerEvidence(None, None, "unknown", None, None, None)
    return FusionPartnerEvidence(
        gene=hit.gene_name,
        preferred_transcript=hit.transcript_id,
        region=hit.region,
        exon_number=hit.exon_number,
        intron_number=hit.intron_number,
        strand=hit.strand,
    )


def fusion_evidence_from_breakpoints(
    pair: AnnotatedBreakpointPair,
    *,
    orientation: str | None = None,
) -> FusionEvidence:
    """Project Gene A/B and transcript context without asserting a fusion frame.

    Frame status stays ``unknown``. GFF CDS phase alone is insufficient to determine which
    coding bases a BND retains; callers may set in/out-of-frame only after an orientation-
    aware transcript-junction implementation supplies that evidence.
    """
    return FusionEvidence(
        gene_a=_partner(pair.primary),
        gene_b=_partner(pair.secondary),
        orientation=normalize_orientation(orientation),
    )
