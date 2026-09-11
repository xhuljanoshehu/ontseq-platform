"""Deterministic inline-SVG plots for the portable HTML report.

Presentation layer only. Every value drawn here is copied from an already
normalized report object; nothing is recomputed, rescaled against a policy or
interpreted. The portable report is archived inside the checksummed release
bundle, so each plot is a deterministic pure-standard-library string: no
JavaScript, no external font and no network reference — the HTML stays offline
and self-contained.
"""

from __future__ import annotations

import html
import math
from collections.abc import Sequence
from dataclasses import dataclass

#: Palette mirrors report.py's stylesheet as concrete values, because an SVG
#: must render identically where CSS variables are unavailable (print path,
#: external viewers).
_INK = "#172033"
_MUTED = "#5e687a"
_ACCENT = "#1d5f8a"
_ACCENT_SOFT = "#dbe9f2"
_REFERENCE = "#8a4b08"
_GRID = "#e3e7ee"

_WIDTH = 1000
_HEIGHT = 260
_LEFT = 64
_RIGHT = 16
_TOP = 18
_BOTTOM = 60  # Leaves room for the chromosome group labels below the axis.


@dataclass(frozen=True)
class CoverageBar:
    """One target's already-normalized mean depth, positioned in genome order."""

    region_id: str
    chromosome: str
    start: int
    mean_depth: float


def chromosome_sort_key(name: str) -> tuple[int, str]:
    """Natural chromosome ordering: 1-22, X, Y, M/MT, then anything else by name."""
    core = name.lower().removeprefix("chr")
    if core.isdigit():
        return (int(core), "")
    special = {"x": 23, "y": 24, "m": 25, "mt": 25}
    if core in special:
        return (special[core], "")
    return (100, core)


def _nice_ceiling(value: float) -> float:
    """Smallest 1/2/2.5/5×10^n step ≥ ``value``, so axis maxima stay readable."""
    if not math.isfinite(value) or value <= 0:
        return 1.0
    # 10.0 keeps the exponentiation typed (and valued) as float for negative exponents.
    magnitude = 10.0 ** math.floor(math.log10(value))
    for factor in (1, 2, 2.5, 5, 10):
        if factor * magnitude >= value:
            return factor * magnitude
    return 10.0 * magnitude


def _text(
    x: float,
    y: float,
    content: str,
    *,
    size: int = 10,
    anchor: str = "middle",
    color: str = _MUTED,
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'fill="{color}" text-anchor="{anchor}">{html.escape(content)}</text>'
    )


def coverage_depth_svg(
    bars: Sequence[CoverageBar],
    *,
    title: str,
    reference_depths: Sequence[float] = (),
) -> str:
    """Render per-target mean depth as an inline-SVG bar plot in genome order.

    ``reference_depths`` are drawn as dashed horizontal lines, each labelled
    with its own value; they are the coverage policy's descriptive bins, shown
    so the plot stays consistent with the tables beside it. A measured zero is
    drawn as a visible flat marker, not as a missing bar. An empty input renders
    as an empty string, so the caller omits the figure entirely rather than
    showing an empty frame.
    """
    if not bars:
        return ""
    ordered = sorted(
        bars,
        key=lambda bar: (chromosome_sort_key(bar.chromosome), bar.start, bar.region_id),
    )
    for bar in ordered:
        if not math.isfinite(bar.mean_depth) or bar.mean_depth < 0:
            raise ValueError("A coverage bar's mean depth must be finite and nonnegative")
    maximum = max(bar.mean_depth for bar in ordered)
    y_max = _nice_ceiling(max(maximum, max(reference_depths, default=0.0)))
    plot_right = _WIDTH - _RIGHT
    plot_bottom = _HEIGHT - _BOTTOM
    span = plot_bottom - _TOP
    count = len(ordered)
    slot = (plot_right - _LEFT) / count
    bar_width = max(1.0, min(40.0, slot * 0.8))

    def y_of(value: float) -> float:
        return plot_bottom - span * (value / y_max)

    parts: list[str] = [
        f'<svg viewBox="0 0 {_WIDTH} {_HEIGHT}" role="img" '
        f'aria-label="{html.escape(title)}" xmlns="http://www.w3.org/2000/svg" '
        f'style="width:100%;height:auto;display:block">',
        f"<title>{html.escape(title)}</title>",
    ]
    # Gridlines and y-axis labels at zero, half and full scale.
    for fraction, label_value in ((0.0, "0"), (0.5, f"{y_max / 2:g}"), (1.0, f"{y_max:g}")):
        y = plot_bottom - span * fraction
        parts.append(
            f'<line x1="{_LEFT}" y1="{y:.1f}" x2="{plot_right}" y2="{y:.1f}" '
            f'stroke="{_GRID}" stroke-width="1"/>'
        )
        parts.append(_text(_LEFT - 8, y + 3.5, f"{label_value}×", anchor="end"))
    # Descriptive policy reference lines, each labelled with its own value.
    for depth in reference_depths:
        y = y_of(depth)
        parts.append(
            f'<line x1="{_LEFT}" y1="{y:.1f}" x2="{plot_right}" y2="{y:.1f}" '
            f'stroke="{_REFERENCE}" stroke-width="1" stroke-dasharray="5 4"/>'
        )
        parts.append(
            _text(plot_right - 4, y - 4, f"{depth:g}×", anchor="end", color=_REFERENCE, size=9)
        )
    # Bars in genome order, carrying their exact values as hover/reader text.
    for index, bar in enumerate(ordered):
        center = _LEFT + slot * (index + 0.5)
        top = y_of(bar.mean_depth)
        tooltip = (
            f"{html.escape(bar.region_id)} · {html.escape(bar.chromosome)}:{bar.start} · "
            f"mean {bar.mean_depth:.2f}×"
        )
        if bar.mean_depth == 0:
            # A measured zero stays visible instead of collapsing into nothing.
            parts.append(
                f'<rect x="{center - bar_width / 2:.2f}" y="{plot_bottom - 2:.2f}" '
                f'width="{bar_width:.2f}" height="2" fill="{_ACCENT_SOFT}">'
                f"<title>{tooltip} (measured zero)</title></rect>"
            )
            continue
        parts.append(
            f'<rect x="{center - bar_width / 2:.2f}" y="{top:.2f}" '
            f'width="{bar_width:.2f}" height="{plot_bottom - top:.2f}" fill="{_ACCENT}">'
            f"<title>{tooltip}</title></rect>"
        )
    # Chromosome separators and labels, thinned when a group gets too narrow.
    group_start = 0
    for index in range(1, count + 1):
        boundary = index == count or ordered[index].chromosome != ordered[group_start].chromosome
        if not boundary:
            continue
        left_edge = _LEFT + slot * group_start
        right_edge = _LEFT + slot * index
        if group_start > 0:
            parts.append(
                f'<line x1="{left_edge:.1f}" y1="{_TOP}" x2="{left_edge:.1f}" '
                f'y2="{plot_bottom}" stroke="{_GRID}" stroke-width="1"/>'
            )
        if right_edge - left_edge >= 26:
            parts.append(
                _text(
                    (left_edge + right_edge) / 2,
                    plot_bottom + 16,
                    ordered[group_start].chromosome,
                    size=9,
                )
            )
        group_start = index
    parts.append(_text(_LEFT, _TOP - 6, "mean depth per target (×)", anchor="start", size=9))
    parts.append("</svg>")
    return "".join(parts)
