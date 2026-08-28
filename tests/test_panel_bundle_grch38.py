from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ontseq_platform.panel_bundle import (
    NORMALIZED_INTERVAL_BASES,
    PANEL_BUNDLE_ID,
    SOURCE_BED_SHA256,
    SOURCE_INTERVAL_COUNT,
    SOURCE_REGIONS_SHA256,
    PanelBundleError,
    import_panel_sources,
    read_source_bed,
    read_source_regions,
    sha256_file,
    validate_panel_sources,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
BUNDLE_ROOT = REPO_ROOT / "configs" / "panels" / PANEL_BUNDLE_ID
SOURCE_BED = BUNDLE_ROOT / "source" / "250611_fusion_panel_with_buffer.bed"
SOURCE_REGIONS = BUNDLE_ROOT / "source" / "250611_fusion_panel_with_buffer.interval_list"
NORMALIZED_BED = BUNDLE_ROOT / "derived" / "selection_panel.normalized.bed"


class CommittedPanelSourceTests(unittest.TestCase):
    def test_sources_are_the_locked_original_bytes(self) -> None:
        self.assertEqual(sha256_file(SOURCE_BED), SOURCE_BED_SHA256)
        self.assertEqual(sha256_file(SOURCE_REGIONS), SOURCE_REGIONS_SHA256)

    def test_sources_have_identical_numeric_intervals_in_order(self) -> None:
        records = validate_panel_sources(SOURCE_BED, SOURCE_REGIONS)
        self.assertEqual(len(records), SOURCE_INTERVAL_COUNT)
        self.assertEqual(
            tuple((item.chromosome, item.start, item.end) for item in records),
            read_source_regions(SOURCE_REGIONS),
        )

    def test_every_normalized_start_is_exactly_one_base_lower(self) -> None:
        source = read_source_bed(SOURCE_BED)
        normalized = [
            line.split("\t")
            for line in NORMALIZED_BED.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ]
        self.assertEqual(len(normalized), SOURCE_INTERVAL_COUNT)
        for original, derivative in zip(source, normalized, strict=True):
            self.assertEqual(int(derivative[1]), original.start - 1)
            self.assertEqual(int(derivative[2]), original.end)

    def test_normalized_span_and_unresolved_igh_are_explicit(self) -> None:
        rows = [line.split("\t") for line in NORMALIZED_BED.read_text().splitlines() if line]
        self.assertEqual(sum(int(row[2]) - int(row[1]) for row in rows), NORMALIZED_INTERVAL_BASES)
        labels = {row[3] for row in rows}
        self.assertIn("IGH_REVIEW_REQUIRED", labels)
        self.assertNotIn("IGH", labels)


class PanelImportTests(unittest.TestCase):
    def test_import_preserves_sources_and_separates_derivative(self) -> None:
        with TemporaryDirectory() as raw:
            destination = Path(raw) / PANEL_BUNDLE_ID
            summary = import_panel_sources(SOURCE_BED, SOURCE_REGIONS, destination)
            self.assertEqual(summary.interval_count, 111)
            self.assertEqual(summary.interval_bases, NORMALIZED_INTERVAL_BASES)
            self.assertEqual(summary.source_bed_sha256, SOURCE_BED_SHA256)
            self.assertEqual(summary.source_regions_sha256, SOURCE_REGIONS_SHA256)
            self.assertNotEqual(summary.normalized_bed_sha256, summary.source_bed_sha256)

    def test_import_refuses_to_overwrite_existing_artifacts(self) -> None:
        with TemporaryDirectory() as raw:
            destination = Path(raw) / PANEL_BUNDLE_ID
            import_panel_sources(SOURCE_BED, SOURCE_REGIONS, destination)
            with self.assertRaisesRegex(PanelBundleError, "overwrite"):
                import_panel_sources(SOURCE_BED, SOURCE_REGIONS, destination)

    def test_disagreeing_regions_list_fails(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            bed = root / "source.bed"
            regions = root / "source.interval_list"
            bed.write_text("chr1\t10\t20\tGENE\n", encoding="utf-8")
            regions.write_text("chr1:11-20\n", encoding="utf-8")
            with self.assertRaisesRegex(PanelBundleError, "identical numeric intervals"):
                validate_panel_sources(bed, regions, require_locked_sources=False)


if __name__ == "__main__":
    unittest.main()
