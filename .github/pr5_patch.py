"""One-shot patch for the CNV genome-overview wiring (deleted after success)."""

from pathlib import Path


def patch(path, replacements):
    p = Path(path)
    data = p.read_bytes()
    nl = '\r\n' if b'\r\n' in data else '\n'
    for old_s, new_s in replacements:
        old = old_s.replace('\n', nl).encode('utf-8')
        new = new_s.replace('\n', nl).encode('utf-8')
        if data.count(old) != 1:
            raise SystemExit(f'{path}: expected 1 match, found {data.count(old)}')
        data = data.replace(old, new, 1)
    p.write_bytes(data)
    print(f'patched {path}')


patch('src/ontseq_platform/report_plots.py', [
    (
        '''    parts.append(_text(_LEFT, _TOP - 6, "reads per length bin", anchor="start", size=9))
    parts.append("</svg>")
    return "".join(parts)
''',
        '''    parts.append(_text(_LEFT, _TOP - 6, "reads per length bin", anchor="start", size=9))
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
    """Median copy number per chromosome with min-max whiskers, as inline SVG.

    Values come verbatim from the normalized QDNAseq/ACE report: the bar is the
    chromosome-level median across bin sizes, the whisker is its min-max range,
    and ``baseline`` is the fitted ACE ploidy drawn as a dashed reference. A
    median of zero is a measured value and stays visible as a flat marker. An
    empty input renders as an empty string, so the caller omits the figure.
    """
    if not chromosomes:
        return ""
    for item in chromosomes:
        values = (item.median_cn, item.min_cn, item.max_cn)
        if (
            any(
                not isinstance(value, int | float) or isinstance(value, bool)
                for value in values
            )
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
''',
    ),
])

patch('src/ontseq_platform/cnv/extension.py', [
    (
        '''from ..report import render_html
from ..sidecars import tabular_sidecar
''',
        '''from ..report import render_html
from ..report_plots import CnvChromosomeBar, cnv_genome_svg
from ..sidecars import tabular_sidecar
''',
    ),
    (
        '''    images: list[str] = []
    for label, name in (
        ("ACE purity/ploidy fit landscape", cnv.primary_fit.fit_plot),
        ("Absolute copy-number profile", cnv.primary_fit.copy_number_plot),
    ):
''',
        '''    bars = [
        CnvChromosomeBar(
            chromosome=chromosome.chromosome,
            median_cn=chromosome.median_copy_number,
            min_cn=chromosome.min_copy_number,
            max_cn=chromosome.max_copy_number,
            rounded_cn=chromosome.rounded_copy_number,
        )
        for chromosome in cnv.chromosome_consensus
    ]
    genome_svg = cnv_genome_svg(
        bars,
        title="Genome-wide copy-number overview",
        baseline=cnv.primary_fit.ploidy,
    )
    genome_figure = (
        "<h3>Genome-wide copy-number overview</h3><figure class='plot' "
        "style='margin:14px 0'>"
        + genome_svg
        + "<figcaption style='color:#5e687a;font-size:12px;margin-top:6px'>Median copy "
        "number per chromosome from the multi-bin consensus; whiskers show each "
        "chromosome's min–max range across bin sizes; the dashed line marks the "
        "fitted ACE ploidy. Autosomes only — no X/Y statement is made or implied. "
        "Descriptive technical evidence, not an adequacy or reportability assessment."
        "</figcaption></figure>"
        if genome_svg
        else ""
    )
    images: list[str] = []
    for label, name in (
        ("ACE purity/ploidy fit landscape", cnv.primary_fit.fit_plot),
        ("Absolute copy-number profile", cnv.primary_fit.copy_number_plot),
    ):
''',
    ),
    (
        '''        f"fit error {cnv.primary_fit.fit_error:.6g}.</p>"
        "<h3>Multi-resolution fits</h3><table><thead><tr><th>Bin (kbp)</th>"
''',
        '''        f"fit error {cnv.primary_fit.fit_error:.6g}.</p>"
        + genome_figure
        + "<h3>Multi-resolution fits</h3><table><thead><tr><th>Bin (kbp)</th>"
''',
    ),
])

patch('tests/test_report_plots.py', [
    (
        '''from ontseq_platform.report_plots import (
    CoverageBar,
    ReadLengthBin,
    chromosome_sort_key,
    coverage_depth_svg,
    read_length_histogram_svg,
)
''',
        '''from ontseq_platform.report_plots import (
    CoverageBar,
    CnvChromosomeBar,
    ReadLengthBin,
    chromosome_sort_key,
    cnv_genome_svg,
    coverage_depth_svg,
    read_length_histogram_svg,
)
''',
    ),
    (
        '''if __name__ == "__main__":
    unittest.main()
''',
        '''class CnvGenomeSvgTests(unittest.TestCase):
    def test_empty_renders_nothing(self) -> None:
        self.assertEqual(cnv_genome_svg([], title="empty", baseline=2.0), "")

    def test_output_is_deterministic_and_wellformed(self) -> None:
        bars = [
            CnvChromosomeBar("chr2", 2.01, 1.98, 2.05, 2),
            CnvChromosomeBar("chr1", 1.02, 0.99, 1.04, 1),
            CnvChromosomeBar("chr10", 3.2, 3.0, 3.4, 3),
        ]
        first = cnv_genome_svg(bars, title="t", baseline=2.0)
        self.assertEqual(first, cnv_genome_svg(bars, title="t", baseline=2.0))
        ET.fromstring(first)

    def test_genome_order_and_labels(self) -> None:
        bars = [
            CnvChromosomeBar("chr10", 3.2, 3.0, 3.4, 3),
            CnvChromosomeBar("chr2", 2.01, 1.98, 2.05, 2),
            CnvChromosomeBar("chr1", 1.02, 0.99, 1.04, 1),
        ]
        svg = cnv_genome_svg(bars, title="t", baseline=2.0)
        self.assertLess(svg.index("chr1 · median"), svg.index("chr2 · median"))
        self.assertLess(svg.index("chr2 · median"), svg.index("chr10 · median"))
        self.assertIn("fitted ploidy 2", svg)
        self.assertIn("range 1.990–2.050", svg)

    def test_measured_zero_median_stays_visible(self) -> None:
        svg = cnv_genome_svg(
            [CnvChromosomeBar("chr5", 0.0, 0.0, 0.1, 0)], title="t", baseline=2.0
        )
        self.assertIn("measured zero", svg)

    def test_invalid_values_are_refused(self) -> None:
        with self.assertRaises(ValueError):
            cnv_genome_svg([CnvChromosomeBar("chr1", 2.0, 2.5, 3.0, 2)], title="t", baseline=2.0)
        with self.assertRaises(ValueError):
            cnv_genome_svg(
                [CnvChromosomeBar("chr1", float("nan"), 1.0, 3.0, 2)],
                title="t",
                baseline=2.0,
            )
        with self.assertRaises(ValueError):
            cnv_genome_svg([CnvChromosomeBar("chr1", 2.0, 1.0, 3.0, 2)], title="t", baseline=0.0)


if __name__ == "__main__":
    unittest.main()
''',
    ),
])

patch('CHANGELOG.md', [
    (
        '''  reportability change.

### Fixed
''',
        '''  reportability change.
- Render a deterministic genome-wide copy-number overview in the CNV section: median
  copy number per chromosome from the multi-bin consensus with min–max whiskers and
  the fitted ACE ploidy as the dashed reference, next to the retained ACE PNG panels.
  Pure presentation of the normalized QDNAseq/ACE report; no fit, threshold, contract
  or reportability change.

### Fixed
''',
    ),
])
