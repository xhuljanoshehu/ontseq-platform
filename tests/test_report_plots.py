"""Regression tests for the portable report's inline-SVG plots."""

from __future__ import annotations

import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from ontseq_platform.demo import build_demo_result
from ontseq_platform.methylation import (
    MODKIT_MODIFICATION_NAMES,
    MethylationPolicy,
    MethylationRegionSource,
    MethylationRegionSummary,
    MethylationReport,
    ModificationCode,
)
from ontseq_platform.models import (
    FileFingerprint,
    GenomeBuild,
    ModuleRunStatus,
    ToolRecord,
)
from ontseq_platform.qc import read_length_histogram_from_tsv
from ontseq_platform.report import render_html
from ontseq_platform.report_plots import (
    CnvChromosomeBar,
    CoverageBar,
    MethylationCell,
    ReadLengthBin,
    chromosome_sort_key,
    cnv_genome_svg,
    coverage_depth_svg,
    methylation_heatmap_svg,
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


def _methylation_row(
    *,
    region_id: str,
    chromosome: str,
    start: int,
    end: int,
    code: ModificationCode,
    fraction: float | None,
    sites_total: int = 10,
    sites_at_minimum_coverage: int = 8,
    valid_call_count: int = 80,
    modified_call_count: int = 60,
) -> MethylationRegionSummary:
    if fraction is None:
        sites_at_minimum_coverage = 0
        valid_call_count = 0
        modified_call_count = 0
        canonical = 0
        mean_coverage = None
        median = None
    else:
        canonical = valid_call_count - modified_call_count
        mean_coverage = valid_call_count / sites_at_minimum_coverage
        median = fraction
    return MethylationRegionSummary(
        region_id=region_id,
        chromosome=chromosome,
        start=start,
        end=end,
        modification_code=code,
        modification_name=MODKIT_MODIFICATION_NAMES[code],
        sites_total=sites_total,
        sites_at_minimum_coverage=sites_at_minimum_coverage,
        valid_call_count=valid_call_count,
        modified_call_count=modified_call_count,
        canonical_call_count=canonical,
        other_mod_call_count=0,
        mean_modified_fraction=fraction,
        median_site_modified_fraction=median,
        mean_valid_coverage=mean_coverage,
    )


def _methylation_report() -> MethylationReport:
    policy = MethylationPolicy(
        profile_id="synthetic-heatmap",
        status="technical_defaults_only",
        expected_version="0.6.4",
        modification_codes=[ModificationCode.FIVE_MC, ModificationCode.FIVE_HMC],
        note="Synthetic heatmap policy",
        region_source=MethylationRegionSource.TARGET_BED,
    )
    regions = [
        _methylation_row(
            region_id="ROI_A",
            chromosome="chr1",
            start=100,
            end=300,
            code=ModificationCode.FIVE_MC,
            fraction=0.75,
        ),
        _methylation_row(
            region_id="ROI_A",
            chromosome="chr1",
            start=100,
            end=300,
            code=ModificationCode.FIVE_HMC,
            fraction=0.1,
            modified_call_count=8,
        ),
        _methylation_row(
            region_id="ROI_ZERO",
            chromosome="chr2",
            start=0,
            end=50,
            code=ModificationCode.FIVE_MC,
            fraction=0.0,
            modified_call_count=0,
        ),
        _methylation_row(
            region_id="ROI_ZERO",
            chromosome="chr2",
            start=0,
            end=50,
            code=ModificationCode.FIVE_HMC,
            fraction=None,
        ),
        _methylation_row(
            region_id="ROI_C",
            chromosome="chr10",
            start=10,
            end=20,
            code=ModificationCode.FIVE_MC,
            fraction=0.4,
            modified_call_count=32,
        ),
        _methylation_row(
            region_id="ROI_C",
            chromosome="chr10",
            start=10,
            end=20,
            code=ModificationCode.FIVE_HMC,
            fraction=None,
        ),
    ]
    return MethylationReport(
        sample_id="SYNTHETIC_001",
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.COMPLETED,
        policy=policy,
        region_source=MethylationRegionSource.TARGET_BED,
        summary_metrics={
            "region_row_count": len(regions),
            "fail_call_count": 3,
            "nocall_call_count": 7,
        },
        regions=regions,
        bedmethyl_fingerprint=FileFingerprint(size_bytes=42, sha256="c" * 64),
        tool=ToolRecord(name="modkit", version="0.6.4", parameters={}),
        warnings=[policy.note],
        limitations=["Descriptive technical measurements only."],
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
        svg = cnv_genome_svg([CnvChromosomeBar("chr5", 0.0, 0.0, 0.1, 0)], title="t", baseline=2.0)
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


class MethylationHeatmapSvgTests(unittest.TestCase):
    def test_empty_renders_nothing(self) -> None:
        self.assertEqual(methylation_heatmap_svg([], title="empty"), "")

    def test_output_is_deterministic_and_wellformed(self) -> None:
        cells = [
            MethylationCell("A", "chr2", 0, "5mC", 0.75, 80, 8, 10),
            MethylationCell("B", "chr1", 10, "5mC", None, 0, 0, 4),
        ]
        first = methylation_heatmap_svg(cells, title="t")
        self.assertEqual(first, methylation_heatmap_svg(cells, title="t"))
        ET.fromstring(first)

    def test_genome_order_and_unmeasurable_cells(self) -> None:
        cells = [
            MethylationCell("TEN", "chr10", 0, "5mC", 0.4, 40, 4, 4),
            MethylationCell("TWO", "chr2", 0, "5mC", 0.0, 80, 8, 8),
            MethylationCell("ONE", "chr1", 0, "5mC", None, 0, 0, 3),
        ]
        svg = methylation_heatmap_svg(cells, title="t")
        self.assertLess(svg.index("ONE · 5mC"), svg.index("TWO · 5mC"))
        self.assertLess(svg.index("TWO · 5mC"), svg.index("TEN · 5mC"))
        self.assertIn("not measurable", svg)
        self.assertIn("url(#meth-nocall)", svg)
        self.assertIn("0.0%", svg)

    def test_modification_row_order_is_canonical(self) -> None:
        cells = [
            MethylationCell("A", "chr1", 0, "5hmC", 0.1, 10, 2, 2),
            MethylationCell("A", "chr1", 0, "5mC", 0.8, 10, 2, 2),
        ]
        svg = methylation_heatmap_svg(cells, title="t")
        five_mc = svg.index(">5mC<")
        five_hmc = svg.index(">5hmC<")
        self.assertLess(five_mc, five_hmc)

    def test_labels_are_escaped(self) -> None:
        svg = methylation_heatmap_svg(
            [MethylationCell("<script>", "chr1", 0, "5mC", 0.5, 10, 2, 2)],
            title="t",
        )
        self.assertNotIn("<script>", svg)
        self.assertIn("&lt;script&gt;", svg)

    def test_invalid_values_are_refused(self) -> None:
        with self.assertRaises(ValueError):
            methylation_heatmap_svg(
                [MethylationCell("A", "chr1", 0, "5mC", 1.2, 10, 2, 2)], title="t"
            )
        with self.assertRaises(ValueError):
            methylation_heatmap_svg(
                [MethylationCell("A", "chr1", 0, "5mC", float("nan"), 10, 2, 2)],
                title="t",
            )
        with self.assertRaises(ValueError):
            methylation_heatmap_svg(
                [
                    MethylationCell("A", "chr1", 0, "5mC", 0.1, 10, 2, 2),
                    MethylationCell("A", "chr1", 0, "5mC", 0.2, 10, 2, 2),
                ],
                title="t",
            )


class ReportMethylationIntegrationTests(unittest.TestCase):
    def test_nested_targets_with_shared_label_keep_distinct_heatmap_cells(self) -> None:
        report = _methylation_report()
        report.regions = [
            _methylation_row(
                region_id="ROI_A",
                chromosome="chr1",
                start=100,
                end=500,
                code="m",
                fraction=0.25,
                modified_call_count=20,
            ),
            _methylation_row(
                region_id="ROI_A",
                chromosome="chr1",
                start=100,
                end=300,
                code="m",
                fraction=0.75,
                modified_call_count=60,
            ),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = render_html(
                build_demo_result(), Path(temporary) / "report.html", methylation_report=report
            )
            document = path.read_text(encoding="utf-8")
        self.assertIn("ROI_A · 5mC · 75.0%", document)
        self.assertIn("ROI_A · 5mC · 25.0%", document)
        self.assertLess(document.index("75.0%"), document.index("25.0%"))

    def test_render_html_includes_the_heatmap_when_supplied(self) -> None:
        result = build_demo_result()
        with tempfile.TemporaryDirectory() as temporary:
            path = render_html(
                result,
                Path(temporary) / "report.html",
                methylation_report=_methylation_report(),
            )
            document = path.read_text(encoding="utf-8")
        self.assertIn("Modified-base fractions", document)
        self.assertIn("not measurable", document)
        self.assertIn("failed-threshold calls 3", document)
        self.assertIn('href="#methylation"', document)

    def test_render_html_without_methylation_omits_the_section(self) -> None:
        result = build_demo_result()
        with tempfile.TemporaryDirectory() as temporary:
            path = render_html(result, Path(temporary) / "report.html")
            document = path.read_text(encoding="utf-8")
        self.assertNotIn("Modified-base fractions", document)
        self.assertNotIn('href="#methylation"', document)


if __name__ == "__main__":
    unittest.main()
