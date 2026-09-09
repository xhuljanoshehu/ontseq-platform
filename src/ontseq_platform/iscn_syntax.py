"""Pure syntax helpers for the deliberately limited CNV event-fragment subset."""

from __future__ import annotations

import re

_VALID_CYTOBAND = re.compile(r"^[pq](?:ter|[0-9]+(?:\.[0-9]+)?)$")
_SUPPORTED_CNV_EVENT_TYPES = {
    "chromosome_gain",
    "chromosome_loss",
    "deletion",
    "duplication",
}


def render_supported_cnv_fragment(
    *,
    event_type: str,
    chromosome: str,
    copy_number: float | None,
    cytoband_start: str | None,
    cytoband_end: str | None,
    whole_chromosome_span_confirmed: bool = False,
) -> str | None:
    """Return one canonical fragment, or ``None`` outside the supported safe subset."""

    chrom = chromosome.removeprefix("chr")
    if event_type not in _SUPPORTED_CNV_EVENT_TYPES or copy_number is None or chrom in {"X", "Y"}:
        return None
    if event_type == "chromosome_gain" and whole_chromosome_span_confirmed:
        return f"+{chrom}"
    if event_type == "chromosome_loss" and whole_chromosome_span_confirmed:
        return f"-{chrom}"
    if event_type in {"chromosome_gain", "chromosome_loss"}:
        return None
    start = cytoband_start.strip() if cytoband_start is not None else None
    end = cytoband_end.strip() if cytoband_end is not None else start
    if start is None or _VALID_CYTOBAND.fullmatch(start) is None:
        return None
    if end is None or _VALID_CYTOBAND.fullmatch(end) is None:
        return None
    band_range = start if start == end else f"{start}{end}"
    prefix = "del" if event_type == "deletion" else "dup"
    return f"{prefix}({chrom})({band_range})"
