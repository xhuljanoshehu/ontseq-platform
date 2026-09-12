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


@dataclass(frozen=True)
class ReadLengthBin:
    """One normalized Cramino read-length bin; ``end=None`` marks the open final bin."""

    start: int
    end: int | None
    count: int
    bases: int


def read_length_histogram_svg(
    bins: Sequence[ReadLengthBin],
    *,
    title: str,
    n50_bp: int | None = None,
) -> str:
    """Render the normalized read-length distribution as an inline-SVG histogram.

    Values come from the Cramino histogram sidecar verbatim. The open final bin
    (no end) is drawn in the soft accent colour so it cannot be mistaken for a
    measured closed bin. An empty input renders as an empty string, so the caller
    omits the figure entirely.
    """
    if not bins:
        return ""
    for item in bins:
        if (
            not isinstance(item.start, int)
            or isinstance(item.start, bool)
            or (
                item.end is not None
                and (not isinstance(item.end, int) or isinstance(item.end, bool))
            )
            or not isinstance(item.count, int)
            or isinstance(item.count, bool)
            or not isinstance(item.bases, int)
            or isinstance(item.bases, bool)
            or item.start < 0
            or (item.end is not None and item.end <= item.start)
            or item.count < 0
            or item.bases < 0
        ):
            raise ValueError("A read-length bin is not numeric/valid")
    ordered = sorted(bins, key=lambda item: (item.start, item.end is None, item.end or 0))
    closed_widths = [item.end - item.start for item in ordered if item.end is not None]
    if not closed_widths:
        return ""
    open_width = closed_widths[len(closed_widths) // 2]
    spans = [
        (item.start, item.end if item.end is not None else item.start + open_width, item)
        for item in ordered
    ]
    x_max = max(end for _, end, _ in spans)
    if x_max <= 0:
        return ""
    plot_right = _WIDTH - _RIGHT
    plot_bottom = _HEIGHT - 46
    span_y = plot_bottom - _TOP
    y_max = _nice_ceiling(max(item.count for _, _, item in spans))
    if y_max <= 0:
        return ""

    def x_of(value: float) -> float:
        return _LEFT + (plot_right - _LEFT) * (value / x_max)

    def y_of(value: float) -> float:
        return plot_bottom - span_y * (value / y_max)

    parts: list[str] = [
        f'<svg viewBox="0 0 {_WIDTH} {_HEIGHT}" role="img" '
        f'aria-label="{html.escape(title)}" xmlns="http://www.w3.org/2000/svg" '
        f'style="width:100%;height:auto;display:block">',
        f"<title>{html.escape(title)}</title>",
    ]
    for fraction, label_value in ((0.0, "0"), (0.5, f"{y_max / 2:g}"), (1.0, f"{y_max:g}")):
        y = plot_bottom - span_y * fraction
        parts.append(
            f'<line x1="{_LEFT}" y1="{y:.1f}" x2="{plot_right}" y2="{y:.1f}" '
            f'stroke="{_GRID}" stroke-width="1"/>'
        )
        parts.append(_text(_LEFT - 8, y + 3.5, label_value, anchor="end"))
    if n50_bp is not None and 0 < n50_bp <= x_max:
        x = x_of(n50_bp)
        parts.append(
            f'<line x1="{x:.1f}" y1="{_TOP}" x2="{x:.1f}" y2="{plot_bottom}" '
            f'stroke="{_REFERENCE}" stroke-width="1" stroke-dasharray="5 4"/>'
        )
        parts.append(_text(x + 4, _TOP + 10, "N50", anchor="start", color=_REFERENCE, size=9))
    for start, end, item in spans:
        x0 = x_of(start)
        x1 = x_of(end)
        top = y_of(item.count)
        end_label = "open" if item.end is None else f"{end:,}"
        tooltip = (
            f"{start:,}–{end_label} bp · {item.count:,} reads · {item.bases / 1_000_000:.2f} Mb"
        )
        fill = _ACCENT if item.end is not None else _ACCENT_SOFT
        parts.append(
            f'<rect x="{x0:.2f}" y="{top:.2f}" width="{max(x1 - x0, 0.5):.2f}" '
            f'height="{plot_bottom - top:.2f}" fill="{fill}">'
            f"<title>{tooltip}</title></rect>"
        )
    for value in (0, x_max // 2, x_max):
        parts.append(_text(x_of(value), plot_bottom + 16, f"{value / 1000:g} kb", size=9))
    parts.append(_text(_LEFT, _TOP - 6, "reads per length bin", anchor="start", size=9))
    parts.append("</svg>")
    return "".join(parts)


@dataclass(frozen=True)
class CnvChromosomeBar:
    """One chromosome's normalized multi-bin copy-number consensus."""

    chromosome: str
    median_cn: float
    min_cn: float
    max_cn: float
    rounded_cn: int


def cnv_genome_svg(
    chromosomes: Sequence[CnvChromosomeBar],
    *,
    title: str,
    baseline: float,
) -> str:
    """Median copy number per chromosome with min–max whiskers, as inline SVG.

    Values come verbatim from the normalized QDNAseq/ACE report: the bar is the
    chromosome-level median across bin sizes, the whisker is its min–max range,
    and ``baseline`` is the fitted ACE ploidy drawn as a dashed reference. A
    median of zero is a measured value and stays visible as a flat marker. An
    empty input renders as an empty string, so the caller omits the figure.
    """
    if not chromosomes:
        return ""
    for item in chromosomes:
        values = (item.median_cn, item.min_cn, item.max_cn)
        if (
            any(not isinstance(value, int | float) or isinstance(value, bool) for value in values)
            or any(not math.isfinite(value) for value in values)
            or item.min_cn < 0
            or item.min_cn > item.median_cn
            or item.median_cn > item.max_cn
            or not isinstance(item.rounded_cn, int)
            or isinstance(item.rounded_cn, bool)
            or item.rounded_cn < 0
        ):
            raise ValueError("A CNV chromosome bar is not numeric/valid")
    if not isinstance(baseline, int | float) or isinstance(baseline, bool):
        raise ValueError("The CNV baseline must be numeric")
    if not math.isfinite(baseline) or baseline <= 0:
        raise ValueError("The CNV baseline must be a positive finite value")
    ordered = sorted(chromosomes, key=lambda item: chromosome_sort_key(item.chromosome))
    y_max = _nice_ceiling(max(max(item.max_cn for item in ordered), baseline, 1.0))
    plot_right = _WIDTH - _RIGHT
    plot_bottom = _HEIGHT - 56
    span = plot_bottom - _TOP
    count = len(ordered)
    slot = (plot_right - _LEFT) / count
    bar_width = max(2.0, min(28.0, slot * 0.55))

    def y_of(value: float) -> float:
        return plot_bottom - span * (value / y_max)

    parts: list[str] = [
        f'<svg viewBox="0 0 {_WIDTH} {_HEIGHT}" role="img" '
        f'aria-label="{html.escape(title)}" xmlns="http://www.w3.org/2000/svg" '
        f'style="width:100%;height:auto;display:block">',
        f"<title>{html.escape(title)}</title>",
    ]
    for fraction, label_value in ((0.0, "0"), (0.5, f"{y_max / 2:g}"), (1.0, f"{y_max:g}")):
        y = plot_bottom - span * fraction
        parts.append(
            f'<line x1="{_LEFT}" y1="{y:.1f}" x2="{plot_right}" y2="{y:.1f}" '
            f'stroke="{_GRID}" stroke-width="1"/>'
        )
        parts.append(_text(_LEFT - 8, y + 3.5, label_value, anchor="end"))
    baseline_y = y_of(baseline)
    parts.append(
        f'<line x1="{_LEFT}" y1="{baseline_y:.1f}" x2="{plot_right}" y2="{baseline_y:.1f}" '
        f'stroke="{_REFERENCE}" stroke-width="1" stroke-dasharray="5 4"/>'
    )
    parts.append(
        _text(
            plot_right - 4,
            baseline_y - 4,
            f"fitted ploidy {baseline:g}",
            anchor="end",
            color=_REFERENCE,
            size=9,
        )
    )
    for index, item in enumerate(ordered):
        center = _LEFT + slot * (index + 0.5)
        tooltip = (
            f"{html.escape(item.chromosome)} · median CN {item.median_cn:.3f} · "
            f"rounded {item.rounded_cn} · range {item.min_cn:.3f}–{item.max_cn:.3f}"
        )
        parts.append(
            f'<line x1="{center:.2f}" y1="{y_of(item.max_cn):.2f}" '
            f'x2="{center:.2f}" y2="{y_of(item.min_cn):.2f}" '
            f'stroke="{_INK}" stroke-width="1"/>'
        )
        for cap_value in (item.min_cn, item.max_cn):
            parts.append(
                f'<line x1="{center - 3:.2f}" y1="{y_of(cap_value):.2f}" '
                f'x2="{center + 3:.2f}" y2="{y_of(cap_value):.2f}" '
                f'stroke="{_INK}" stroke-width="1"/>'
            )
        if item.median_cn == 0:
            # A measured zero median stays visible instead of vanishing.
            parts.append(
                f'<rect x="{center - bar_width / 2:.2f}" y="{plot_bottom - 2:.2f}" '
                f'width="{bar_width:.2f}" height="2" fill="{_ACCENT_SOFT}">'
                f"<title>{tooltip} (measured zero)</title></rect>"
            )
        else:
            parts.append(
                f'<rect x="{center - bar_width / 2:.2f}" y="{y_of(item.median_cn):.2f}" '
                f'width="{bar_width:.2f}" '
                f'height="{plot_bottom - y_of(item.median_cn):.2f}" fill="{_ACCENT}">'
                f"<title>{tooltip}</title></rect>"
            )
        parts.append(_text(center, plot_bottom + 14, item.chromosome, size=8))
    parts.append(
        _text(_LEFT, _TOP - 6, "median copy number per chromosome", anchor="start", size=9)
    )
    parts.append("</svg>")
    return "".join(parts)


@dataclass(frozen=True)
class MethylationCell:
    """One already-normalized region × modification measurement."""

    region_id: str
    chromosome: str
    start: int | None
    modification_label: str
    fraction: float | None
    valid_call_count: int
    sites_at_minimum_coverage: int
    sites_total: int


_METH_ZERO = (244, 246, 249)
_METH_ONE = (29, 95, 138)
_CODE_RANK = {"5mC": 0, "5hmC": 1, "6mA": 2}
_NOCALL_FILL = "#d9dde4"
_NOCALL_HATCH = "#9aa3b2"


def _rgb_hex(channels: tuple[int, int, int]) -> str:
    return f"#{channels[0]:02x}{channels[1]:02x}{channels[2]:02x}"


def _methylation_fill(fraction: float) -> str:
    mixed = [
        round(start + (end - start) * fraction)
        for start, end in zip(_METH_ZERO, _METH_ONE, strict=True)
    ]
    return _rgb_hex((mixed[0], mixed[1], mixed[2]))


def methylation_heatmap_svg(cells: Sequence[MethylationCell], *, title: str) -> str:
    """Region × modification heatmap of already-normalized modified fractions.

    ``fraction=None`` is a coverage-floor miss and is drawn hatched, never as a
    measured zero. A measured 0.0 is a pale filled cell. An empty input renders
    as an empty string, so the caller omits the figure entirely.
    """
    if not cells:
        return ""
    seen: set[tuple[str, int | None, str, str]] = set()
    for item in cells:
        counts = (item.valid_call_count, item.sites_at_minimum_coverage, item.sites_total)
        if (
            any(
                not isinstance(value, int) or isinstance(value, bool) or value < 0
                for value in counts
            )
            or item.sites_at_minimum_coverage > item.sites_total
        ):
            raise ValueError("A methylation heatmap cell is not numeric/valid")
        if item.start is not None and (
            not isinstance(item.start, int) or isinstance(item.start, bool) or item.start < 0
        ):
            raise ValueError("A methylation heatmap cell is not numeric/valid")
        if item.fraction is not None and (
            not isinstance(item.fraction, int | float)
            or isinstance(item.fraction, bool)
            or not math.isfinite(item.fraction)
            or item.fraction < 0
            or item.fraction > 1
        ):
            raise ValueError("A methylation heatmap fraction is not in 0-1")
        key = (item.chromosome, item.start, item.region_id, item.modification_label)
        if key in seen:
            raise ValueError("Methylation heatmap cells contain a duplicate region/code pair")
        seen.add(key)

    def column_key(
        column: tuple[str, int | None, str],
    ) -> tuple[tuple[int, str], bool, int, str]:
        chromosome, start, region_id = column
        return (chromosome_sort_key(chromosome), start is None, start or 0, region_id)

    columns = sorted(
        {(item.chromosome, item.start, item.region_id) for item in cells},
        key=column_key,
    )
    rows = sorted(
        {item.modification_label for item in cells},
        key=lambda name: (_CODE_RANK.get(name, 50), name),
    )
    lookup = {
        (item.chromosome, item.start, item.region_id, item.modification_label): item
        for item in cells
    }
    left = 72
    top = 40
    bottom = 52
    cell_height = 28
    height = top + cell_height * len(rows) + bottom
    plot_right = _WIDTH - _RIGHT
    slot = (plot_right - left) / len(columns)
    gap = 1.0
    cell_width = max(2.0, slot - gap)

    parts: list[str] = [
        f'<svg viewBox="0 0 {_WIDTH} {height}" role="img" '
        f'aria-label="{html.escape(title)}" xmlns="http://www.w3.org/2000/svg" '
        f'style="width:100%;height:auto;display:block">',
        f"<title>{html.escape(title)}</title>",
        '<defs><pattern id="meth-nocall" patternUnits="userSpaceOnUse" '
        'width="6" height="6"><rect width="6" height="6" '
        f'fill="{_NOCALL_FILL}"/>'
        f'<path d="M0 6 L6 0" stroke="{_NOCALL_HATCH}" stroke-width="1"/>'
        "</pattern>"
        '<linearGradient id="meth-scale" x1="0" x2="1" y1="0" y2="0">'
        f'<stop offset="0" stop-color="{_rgb_hex(_METH_ZERO)}"/>'
        f'<stop offset="1" stop-color="{_rgb_hex(_METH_ONE)}"/>'
        "</linearGradient></defs>",
    ]
    parts.append(f'<rect x="{left:.1f}" y="10" width="120" height="8" fill="url(#meth-scale)"/>')
    parts.append(_text(left, 8, "0%", anchor="start", size=9))
    parts.append(_text(left + 120, 8, "100%", anchor="end", size=9))
    parts.append(
        f'<rect x="{left + 136:.1f}" y="10" width="12" height="8" fill="url(#meth-nocall)"/>'
    )
    parts.append(_text(left + 152, 17, "not measurable", anchor="start", size=9))
    for row_index, label in enumerate(rows):
        y = top + cell_height * row_index
        parts.append(_text(left - 8, y + cell_height / 2 + 3, label, anchor="end", size=10))
        for column_index, column in enumerate(columns):
            current = lookup.get((*column, label))
            x = left + slot * column_index + gap / 2
            if current is None:
                continue
            region = html.escape(current.region_id)
            code = html.escape(current.modification_label)
            if current.fraction is None:
                tooltip = f"{region} · {code} · not measurable"
                fill = "url(#meth-nocall)"
            else:
                tooltip = (
                    f"{region} · {code} · {current.fraction:.1%} · "
                    f"{current.valid_call_count} valid calls"
                )
                fill = _methylation_fill(float(current.fraction))
            parts.append(
                f'<rect x="{x:.2f}" y="{y + 2:.2f}" width="{cell_width:.2f}" '
                f'height="{cell_height - 4:.2f}" fill="{fill}">'
                f"<title>{tooltip}</title></rect>"
            )
    group_start = 0
    plot_bottom = top + cell_height * len(rows)
    for index in range(1, len(columns) + 1):
        boundary = index == len(columns) or columns[index][0] != columns[group_start][0]
        if not boundary:
            continue
        left_edge = left + slot * group_start
        right_edge = left + slot * index
        if group_start > 0:
            parts.append(
                f'<line x1="{left_edge:.1f}" y1="{top}" x2="{left_edge:.1f}" '
                f'y2="{plot_bottom}" stroke="{_GRID}" stroke-width="1"/>'
            )
        if right_edge - left_edge >= 26:
            parts.append(
                _text(
                    (left_edge + right_edge) / 2,
                    plot_bottom + 14,
                    columns[group_start][0],
                    size=9,
                )
            )
        group_start = index
    parts.append(_text(left, top - 8, "modified fraction by region", anchor="start", size=9))
    parts.append("</svg>")
    return "".join(parts)
