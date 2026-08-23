from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.models import GenomeBuild, ModuleRunStatus, ToolRecord
from ontseq_platform.target_coverage import (
    TargetCoveragePolicy,
    load_target_bed,
    normalize_target_coverage,
)


class TargetCoverageTests(unittest.TestCase):
    def test_normalizes_regions_and_thresholds_without_source_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / "targets.bed"
            regions = root / "sample.regions.bed.gz"
            thresholds = root / "sample.thresholds.bed.gz"
            bed.write_text(
                "chr1\t100\t200\tGENE_A\nchr2\t300\t500\tGENE_B\n",
                encoding="utf-8",
            )
            with gzip.open(regions, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t100\t200\tGENE_A\t25.0\n")
                handle.write("chr2\t300\t500\tGENE_B\t12.0\n")
            with gzip.open(thresholds, "wt", encoding="utf-8") as handle:
                handle.write("#chrom\tstart\tend\tregion\t1X\t10X\t20X\t30X\n")
                handle.write("chr1\t100\t200\tGENE_A\t100\t95\t80\t40\n")
                handle.write("chr2\t300\t500\tGENE_B\t200\t180\t60\t10\n")

            policy = TargetCoveragePolicy(
                profile_id="synthetic",
                status="technical_defaults_only",
                thresholds=[1, 10, 20, 30],
                note="Synthetic technical policy",
            )
            report = normalize_target_coverage(
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                target_bed=bed,
                target_bed_version="synthetic-v1",
                regions_path=regions,
                thresholds_path=thresholds,
                policy=policy,
                tool=ToolRecord(name="mosdepth", version="0.3.14"),
            )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.summary_metrics["region_count"], 2)
        self.assertEqual(report.summary_metrics["interval_bases"], 300)
        self.assertAlmostEqual(report.summary_metrics["interval_weighted_mean_depth"], 16.3333333)
        self.assertAlmostEqual(report.summary_metrics["interval_bases_at_10x_fraction"], 275 / 300)
        self.assertEqual(report.regions[0].region_id, "GENE_A")
        self.assertEqual(report.regions[0].fraction_at_threshold["20x"], 0.8)
        self.assertEqual(report.target_bed_role, "analysis_roi_unbuffered")
        self.assertNotIn(str(bed), report.model_dump_json())

    def test_rejects_duplicate_target_intervals(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / "targets.bed"
            bed.write_text("chr1\t100\t200\tA\nchr1\t100\t200\tB\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                load_target_bed(bed)

    def test_overlapping_intervals_are_explicitly_warned(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / "targets.bed"
            regions = root / "regions.bed.gz"
            thresholds = root / "thresholds.bed.gz"
            bed.write_text("chr1\t100\t200\tA\nchr1\t150\t250\tB\n", encoding="utf-8")
            with gzip.open(regions, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t100\t200\tA\t20\nchr1\t150\t250\tB\t20\n")
            with gzip.open(thresholds, "wt", encoding="utf-8") as handle:
                handle.write("#chrom\tstart\tend\tregion\t1X\n")
                handle.write("chr1\t100\t200\tA\t100\nchr1\t150\t250\tB\t100\n")
            report = normalize_target_coverage(
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                target_bed=bed,
                target_bed_version="v1",
                regions_path=regions,
                thresholds_path=thresholds,
                policy=TargetCoveragePolicy(
                    profile_id="synthetic",
                    status="technical_defaults_only",
                    thresholds=[1],
                    note="Synthetic technical policy",
                ),
                tool=ToolRecord(name="mosdepth", version="0.3.14"),
            )
        self.assertEqual(report.summary_metrics["overlapping_interval_count"], 1)
        self.assertTrue(any("overlap" in warning.lower() for warning in report.warnings))


if __name__ == "__main__":
    unittest.main()
