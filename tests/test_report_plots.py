"""Regression tests for the portable report's inline-SVG coverage plot."""

from __future__ import annotations

import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import (
    FileFingerprint,
    GenomeBuild,
    ModuleRunStatus,
    ToolRecord,
)
from ontseq_platform.qc import read_length_histogram_from_tsv
from ontseq_platform.report import render_html
from ontseq_platform.report_plots import (
    CoverageBar,
    CnvChromosomeBar,
    ReadLengthBin,
    chromosome_sort_key,
    cnv_genome_svg,
    coverage_depth_svg,
    read_length_histogram_svg,
)
from ontseq_platform.report_sections import coverage_section
from ontseq_platform.target_coverage import (
    TargetCoveragePolicy,
    TargetCoverageRegion,
    TargetCoverageReport,
)


def _region(
    chromosome: str,
    start: int,
    end: int,
    region_id: str,
    mean_depth: float,
    bases: dict[str, int],
) -> TargetCoverageRegion:
    length = end - start
    return TargetCoverageRegion(
        chromosome=chromosome,
        start=start,
        end=end,
        region_id=region_id,
        mean_depth=mean_depth,
        bases_at_threshold=bases,
        fraction_at_threshold={key: value / length for key, value in bases.items()},
    )


def _report() -> TargetCoverageReport:
    policy = TargetCoveragePolicy(
        profile_id="synthetic-coverage-plot",
        status="technical_defaults_only",
        thresholds=[1, 10, 20, 30],
        note="Synthetic plot policy",
    )
    regions = [
        _region(
            "chr1",
            100,
            300,
            "ROI_A",
            25.0,
            {"1x": 200, "10x": 180, "20x": 120, "30x": 40},
        ),
        _region("chr1", 400, 500, "ROI_ZERO", 0.0, {"1x": 0, "10x": 0, "20x": 0, "30x": 0}),
        _region("chr2", 0, 50, "ROI_B", 60.0, {"1x": 50, "10x": 50, "20x": 48, "30x": 30}),
    ]
    interval_bases = sum(region.end - region.start for region in regions)
    summary_metrics: dict[str, float | int] = {
        "region_count": len(regions),
        "interval_bases": interval_bases,
        "interval_weighted_mean_depth": sum(
            region.mean_depth * (region.end - region.start) for region in regions
        )
        / interval_bases,
        "minimum_region_mean_depth": min(region.mean_depth for region in regions),
        "median_region_mean_depth": 25.0,
        "maximum_region_mean_depth": max(region.mean_depth for region in regions),
        "overlapping_interval_count": 0,
    }
    for key in ("1x", "10x", "20x", "30x"):
        summary_metrics[f"interval_bases_at_{key}_fraction"] = (
            sum(region.bases_at_threshold[key] for region in regions) / interval_bases
        )
    return TargetCoverageReport(
        sample_id="SYNTHETIC_001",
        genome_build=GenomeBuild.GRCH38,
        target_bed_version="synthetic-v1",
        status=ModuleRunStatus.COMPLETED,
        policy=policy,
        summary_metrics=summary_metrics,
        regions=regions,
        target_bed_fingerprint=FileFingerprint(size_bytes=42, sha256="c" * 64),
        tool=ToolRecord(name="mosdepth", version="0.3.14", parameters={}),
        warnings=[policy.note],
        limitations=["Technical bins are descriptive only."],
    )


class CoverageDepthSvgTests(unittest.TestCase):
    def test_empty_input_renders_nothing(self) -> None:
        """An absent measurement becomes no figure at all, never an empty frame."""
        self.assertEqual(coverage_depth_svg([], title="empty"), "")

    def test_output_is_deterministic(self) -> None:
        bars = [CoverageBar("A", "chr1", 0, 12.5), CoverageBar("B", "chr2", 10, 30.0)]
        first = coverage_depth_svg(bars, title="t", reference_depths=[20])
        self.assertEqual(first, coverage_depth_svg(bars, title="t", reference_depths=[20]))

    def test_genome_order_not_input_order(self) -> None:
        bars = [
            CoverageBar("TEN", "chr10", 0, 5.0),
            CoverageBar("TWO", "chr2", 0, 5.0),
            CoverageBar("X", "chrX", 0, 5.0),
        ]
        svg = coverage_depth_svg(bars, title="t")
        self.assertLess(svg.index("TWO · chr2:0"), svg.index("TEN · chr10:0"))
        self.assertLess(svg.index("TEN · chr10:0"), svg.index("X · chrX:0"))

    def test_labels_are_escaped(self) -> None:
        svg = coverage_depth_svg(
            [CoverageBar("<script>alert(1)</script>", "chr1", 0, 5.0)], title="t"
        )
        self.assertNotIn("<script>", svg)
        self.assertIn("&lt;script&gt;", svg)

    def test_a_measured_zero_stays_visible(self) -> None:
        """A target measured at zero depth is drawn, not dropped like a missing bar."""
        svg = coverage_depth_svg(
            [CoverageBar("ZERO", "chr1", 0, 0.0), CoverageBar("A", "chr1", 5, 20.0)],
            title="t",
        )
        self.assertIn("measured zero", svg)

    def test_reference_lines_carry_their_values(self) -> None:
        svg = coverage_depth_svg(
            [CoverageBar("A", "chr1", 0, 55.0)], title="t", reference_depths=[20, 30]
        )
        self.assertIn("20×", svg)
        self.assertIn("30×", svg)
        self.assertIn('stroke-dasharray="5 4"', svg)

    def test_nonfinite_or_negative_depth_is_refused(self) -> None:
        for bad in (float("nan"), float("inf"), -1.0):
            with self.assertRaises(ValueError):
                coverage_depth_svg([CoverageBar("A", "chr1", 0, bad)], title="t")

    def test_output_is_wellformed_xml(self) -> None:
        svg = coverage_depth_svg(
            [CoverageBar("A", "chr1", 0, 12.5), CoverageBar("B", "chr2", 10, 30.0)],
            title="wellformed",
            reference_depths=[20],
        )
        ET.fromstring(svg)

    def test_chromosome_sort_key_is_natural(self) -> None:
        self.assertLess(chromosome_sort_key("chr2"), chromosome_sort_key("chr10"))
        self.assertLess(chromosome_sort_key("22"), chromosome_sort_key("chrX"))
        self.assertLess(chromosome_sort_key("chrX"), chromosome_sort_key("chrY"))


class CoverageSectionPlotTests(unittest.TestCase):
    def test_section_embeds_the_plot_when_a_report_is_supplied(self) -> None:
        document = coverage_section(_report(), None)
        self.assertIn("<svg", document)
        self.assertIn("figcaption", document)
        self.assertIn("ROI_A", document)
        self.assertIn("measured zero", document)
        self.assertIn("20×", document)

    def test_section_without_a_report_has_no_plot(self) -> None:
        document = coverage_section(None, None)
        self.assertNotIn("<svg", document)
        self.assertIn("Coverage is not assessed.", document)


class ReadLengthHistogramTsvTests(unittest.TestCase):
    def test_roundtrip_parse(self) -> None:
        text = "start_bp\tend_bp\tread_count\tbase_count\n0\t1000\t50\t40000\n1000\t\t12\t36000\n"
        self.assertEqual(
            read_length_histogram_from_tsv(text),
            [(0, 1000, 50, 40000), (1000, None, 12, 36000)],
        )

    def test_empty_text_parses_to_no_bins(self) -> None:
        self.assertEqual(read_length_histogram_from_tsv(""), [])

    def test_malformed_rows_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            read_length_histogram_from_tsv("wrong\theader\n")
        with self.assertRaises(ValueError):
            read_length_histogram_from_tsv(
                "start_bp\tend_bp\tread_count\tbase_count\n0\t100\t-1\t0\n"
            )
        with self.assertRaises(ValueError):
            read_length_histogram_from_tsv("start_bp\tend_bp\tread_count\tbase_count\n0\tx\t1\t0\n")


class ReadLengthHistogramSvgTests(unittest.TestCase):
    def test_empty_renders_nothing(self) -> None:
        self.assertEqual(read_length_histogram_svg([], title="empty"), "")

    def test_output_is_deterministic_and_wellformed(self) -> None:
        bins = [
            ReadLengthBin(0, 1000, 50, 40000),
            ReadLengthBin(1000, 2000, 120, 180000),
            ReadLengthBin(2000, None, 30, 90000),
        ]
        first = read_length_histogram_svg(bins, title="t", n50_bp=1500)
        self.assertEqual(first, read_length_histogram_svg(bins, title="t", n50_bp=1500))
        ET.fromstring(first)

    def test_n50_marker_and_open_bin_are_drawn(self) -> None:
        bins = [ReadLengthBin(0, 1000, 50, 40000), ReadLengthBin(1000, None, 30, 90000)]
        svg = read_length_histogram_svg(bins, title="t", n50_bp=800)
        self.assertIn("N50", svg)
        self.assertIn("open", svg)

    def test_invalid_bins_are_refused(self) -> None:
        with self.assertRaises(ValueError):
            read_length_histogram_svg([ReadLengthBin(0, 1000, -1, 0)], title="t")
        with self.assertRaises(ValueError):
            read_length_histogram_svg([ReadLengthBin(1000, 500, 1, 10)], title="t")


class ReportHistogramIntegrationTests(unittest.TestCase):
    def test_render_html_includes_the_histogram_when_supplied(self) -> None:
        result = build_demo_result()
        bins = [
            ReadLengthBin(0, 1000, 50, 40000),
            ReadLengthBin(1000, 2000, 120, 180000),
            ReadLengthBin(2000, None, 30, 90000),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = render_html(result, Path(temporary) / "report.html", qc_histogram=bins)
            document = path.read_text(encoding="utf-8")
        self.assertIn("Read length distribution", document)

    def test_render_html_without_histogram_renders_no_figure(self) -> None:
        result = build_demo_result()
        with tempfile.TemporaryDirectory() as temporary:
            path = render_html(result, Path(temporary) / "report.html")
            document = path.read_text(encoding="utf-8")
        self.assertNotIn("Read length distribution", document)


class CnvGenomeSvgTests(unittest.TestCase):
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
        self.assertIn("range 0.990–1.040", svg)

    def test_measured_zero_median_stays_visible(self) -> None:
        svg = cnv_genome_svg(
            [CnvChromosomeBar("chr5", 0.0, 0.0, 0.1, 0)], title="t", baseline=2.0
        )
        self.assertIn("measured zero", svg)

    def test_invalid_values_are_refused(self) -> None:
        with self.assertRaises(ValueError):
            cnv_genome_svg(
                [CnvChromosomeBar("chr1", 2.0, 2.5, 3.0, 2)], title="t", baseline=2.0
            )
        with self.assertRaises(ValueError):
            cnv_genome_svg(
                [CnvChromosomeBar("chr1", float("nan"), 1.0, 3.0, 2)],
                title="t",
                baseline=2.0,
            )
        with self.assertRaises(ValueError):
            cnv_genome_svg(
                [CnvChromosomeBar("chr1", 2.0, 1.0, 3.0, 2)], title="t", baseline=0.0
            )


if __name__ == "__main__":
    unittest.main()
