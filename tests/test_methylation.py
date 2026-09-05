from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.execution import CommandResult
from ontseq_platform.methylation import (
    MethylationPolicy,
    MethylationRegionSource,
    ModificationCode,
    _bed_regions,
    modkit_version,
    normalize_methylation,
    parse_bedmethyl,
    run_methylation,
)
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    ModuleRunStatus,
    SampleManifest,
    ToolRecord,
    Verdict,
)


def _policy(**overrides: object) -> MethylationPolicy:
    defaults: dict[str, object] = {
        "profile_id": "synthetic",
        "status": "technical_defaults_only",
        "expected_version": "0.4.1",
        "minimum_valid_coverage": 5,
        "note": "Synthetic technical policy",
    }
    defaults.update(overrides)
    return MethylationPolicy.model_validate(defaults)


def _row(
    chromosome: str,
    start: int,
    code: str,
    valid: int,
    modified: int,
) -> str:
    """One bedMethyl record. Only columns 1-4, 10 and 12 carry meaning for this adapter."""
    percent = (modified / valid * 100) if valid else 0.0
    canonical = valid - modified
    return "\t".join(
        [
            chromosome,
            str(start),
            str(start + 1),
            code,
            str(min(valid, 1000)),
            "+",
            str(start),
            str(start + 1),
            "255,0,0",
            str(valid),
            f"{percent:.2f}",
            str(modified),
            str(canonical),
            "0",
            "0",
            "0",
            "0",
            "0",
        ]
    )


class ParsingTests(unittest.TestCase):
    def test_reads_tab_separated_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.bedmethyl"
            path.write_text(
                "\n".join([_row("chr1", 100, "m", 20, 15), _row("chr1", 200, "m", 10, 2)]) + "\n",
                encoding="utf-8",
            )
            sites, skipped = parse_bedmethyl(path, allowed_codes=[ModificationCode.FIVE_MC])
        self.assertEqual(len(sites), 2)
        self.assertEqual(skipped, 0)
        self.assertEqual(sites[0].valid_coverage, 20)
        self.assertEqual(sites[0].modified_calls, 15)

    def test_reads_the_space_delimited_layout_older_builds_emit(self) -> None:
        """A file whose trailing count columns are space-separated is still a bedMethyl."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.bedmethyl"
            fields = _row("chr1", 100, "m", 20, 15).split("\t")
            path.write_text(
                "\t".join(fields[:9]) + "\t" + " ".join(fields[9:]) + "\n", encoding="utf-8"
            )
            sites, _ = parse_bedmethyl(path, allowed_codes=[ModificationCode.FIVE_MC])
        self.assertEqual(len(sites), 1)
        self.assertEqual(sites[0].valid_coverage, 20)

    def test_refuses_a_modification_code_the_policy_does_not_list(self) -> None:
        """Dropping an unexpected code would silently change what the fraction means."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.bedmethyl"
            path.write_text(_row("chr1", 100, "h", 20, 15) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as raised:
                parse_bedmethyl(path, allowed_codes=[ModificationCode.FIVE_MC])
        self.assertIn("policy does not list", str(raised.exception))

    def test_non_canonical_contigs_are_counted_not_silently_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.bedmethyl"
            path.write_text(
                "\n".join(
                    [_row("chr1", 100, "m", 20, 15), _row("chrUn_KI270302v1", 5, "m", 20, 15)]
                )
                + "\n",
                encoding="utf-8",
            )
            sites, skipped = parse_bedmethyl(path, allowed_codes=[ModificationCode.FIVE_MC])
        self.assertEqual(len(sites), 1)
        self.assertEqual(skipped, 1)

    def test_rejects_more_modified_calls_than_valid_calls(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.bedmethyl"
            path.write_text(_row("chr1", 100, "m", 5, 9) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                parse_bedmethyl(path, allowed_codes=[ModificationCode.FIVE_MC])

    def test_parses_the_version_banner(self) -> None:
        self.assertEqual(modkit_version("mod_kit 0.4.1"), "0.4.1")


class NormalizationTests(unittest.TestCase):
    def _report(self, rows: list[str], **policy_overrides: object):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.bedmethyl"
            path.write_text("\n".join(rows) + "\n", encoding="utf-8")
            return normalize_methylation(
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                bedmethyl_path=path,
                policy=_policy(**policy_overrides),
                tool=ToolRecord(name="modkit", version="0.4.1"),
            )

    def test_aggregates_per_chromosome_by_call_weight(self) -> None:
        report = self._report(
            [
                _row("chr1", 100, "m", 20, 15),
                _row("chr1", 200, "m", 10, 5),
                _row("chr2", 300, "m", 8, 0),
            ]
        )
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.region_source, MethylationRegionSource.CHROMOSOME)
        self.assertEqual([item.region_id for item in report.regions], ["chr1", "chr2"])
        chr1 = report.regions[0]
        self.assertEqual(chr1.sites_at_minimum_coverage, 2)
        self.assertEqual(chr1.valid_call_count, 30)
        self.assertEqual(chr1.modified_call_count, 20)
        self.assertAlmostEqual(chr1.mean_modified_fraction or 0.0, 20 / 30)
        self.assertAlmostEqual(chr1.median_site_modified_fraction or 0.0, 0.625)
        self.assertEqual(report.regions[1].mean_modified_fraction, 0.0)

    def test_a_region_with_no_qualifying_site_reports_none_not_zero(self) -> None:
        """An unmeasured region is not an unmethylated region."""
        report = self._report([_row("chr1", 100, "m", 20, 15), _row("chr2", 300, "m", 2, 1)])
        by_id = {item.region_id: item for item in report.regions}
        self.assertEqual(by_id["chr2"].sites_total, 1)
        self.assertEqual(by_id["chr2"].sites_at_minimum_coverage, 0)
        self.assertIsNone(by_id["chr2"].mean_modified_fraction)
        self.assertEqual(by_id["chr2"].valid_call_count, 0)

    def test_a_pileup_below_the_coverage_floor_is_no_call(self) -> None:
        report = self._report([_row("chr1", 100, "m", 2, 1)])
        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertTrue(
            any("not say the DNA is" in warning for warning in report.warnings),
            report.warnings,
        )

    def test_a_version_mismatch_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            self._report([_row("chr1", 100, "m", 20, 15)], expected_version="0.5.0")

    def test_the_report_carries_no_filesystem_paths(self) -> None:
        """Reviewer artifacts name checksums, never the operator's directory layout."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.bedmethyl"
            path.write_text(_row("chr1", 100, "m", 20, 15) + "\n", encoding="utf-8")
            report = normalize_methylation(
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                bedmethyl_path=path,
                policy=_policy(),
                tool=ToolRecord(name="modkit", version="0.4.1"),
            )
            self.assertNotIn(directory, report.model_dump_json())
        self.assertIsNotNone(report.bedmethyl_fingerprint.sha256)


class RegionAggregationTests(unittest.TestCase):
    def test_sites_are_counted_in_every_overlapping_target(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / "targets.bed"
            bed.write_text("chr1\t100\t300\tA\nchr1\t200\t400\tB\n", encoding="utf-8")
            path = root / "sample.bedmethyl"
            path.write_text(
                "\n".join(
                    [
                        _row("chr1", 150, "m", 10, 5),
                        _row("chr1", 250, "m", 10, 10),
                        _row("chr1", 350, "m", 10, 0),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            report = normalize_methylation(
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                bedmethyl_path=path,
                policy=_policy(region_source="target_bed"),
                tool=ToolRecord(name="modkit", version="0.4.1"),
                regions=_bed_regions(bed),
                target_bed=bed,
            )
        by_id = {item.region_id: item for item in report.regions}
        self.assertEqual(by_id["A"].sites_at_minimum_coverage, 2)
        self.assertEqual(by_id["B"].sites_at_minimum_coverage, 2)
        self.assertAlmostEqual(by_id["A"].mean_modified_fraction or 0.0, 15 / 20)
        self.assertAlmostEqual(by_id["B"].mean_modified_fraction or 0.0, 10 / 20)
        self.assertIsNotNone(report.target_bed_fingerprint)

    def test_two_intervals_sharing_a_name_stay_separate_rows(self) -> None:
        """A panel names intervals; it does not promise the names are unique.

        Two exons of one gene are routinely both called TP53. Bucketing on the name merged
        them, then emitted one row per interval reading the merged bucket — each row
        counting sites from outside its own range, and the totals double-counted.
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / "targets.bed"
            bed.write_text("chr1\t100\t200\tTP53\nchr1\t300\t400\tTP53\n", encoding="utf-8")
            path = root / "sample.bedmethyl"
            path.write_text(
                "\n".join([_row("chr1", 150, "m", 10, 5), _row("chr1", 350, "m", 10, 9)]) + "\n",
                encoding="utf-8",
            )
            report = normalize_methylation(
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                bedmethyl_path=path,
                policy=_policy(region_source="target_bed"),
                tool=ToolRecord(name="modkit", version="0.4.1"),
                regions=_bed_regions(bed),
                target_bed=bed,
            )
        self.assertEqual(len(report.regions), 2)
        for row in report.regions:
            self.assertEqual(row.region_id, "TP53")
            self.assertEqual(row.sites_total, 1)
        by_start = {row.start: row for row in report.regions}
        self.assertAlmostEqual(by_start[100].mean_modified_fraction or 0.0, 5 / 10)
        self.assertAlmostEqual(by_start[300].mean_modified_fraction or 0.0, 9 / 10)

    def test_a_site_outside_every_target_reaches_no_region(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / "targets.bed"
            bed.write_text("chr1\t100\t200\tA\n", encoding="utf-8")
            path = root / "sample.bedmethyl"
            path.write_text(
                "\n".join([_row("chr1", 150, "m", 10, 5), _row("chr2", 900, "m", 10, 5)]) + "\n",
                encoding="utf-8",
            )
            report = normalize_methylation(
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                bedmethyl_path=path,
                policy=_policy(region_source="target_bed"),
                tool=ToolRecord(name="modkit", version="0.4.1"),
                regions=_bed_regions(bed),
                target_bed=bed,
            )
        self.assertEqual(len(report.regions), 1)
        self.assertEqual(report.regions[0].sites_total, 1)
        self.assertEqual(report.summary_metrics["site_count"], 2)


class _FakeRunner:
    """Stand in for modkit and samtools, writing the pileup a real run would produce."""

    def __init__(self, *, version: str = "0.4.1", tagged_reads: str = "1200") -> None:
        self.version = version
        self.tagged_reads = tagged_reads
        self.commands: list[list[str]] = []
        self.rows: list[str] = [_row("chr1", 100, "m", 20, 15)]

    def run(self, argv, *, timeout_seconds: int = 300) -> CommandResult:
        argv = [str(item) for item in argv]
        self.commands.append(argv)
        if argv[1:2] == ["--version"]:
            return CommandResult(tuple(argv), 0, f"mod_kit {self.version}", "")
        if "view" in argv:
            return CommandResult(tuple(argv), 0, self.tagged_reads, "")
        if "pileup" in argv:
            Path(argv[3]).write_text("\n".join(self.rows) + "\n", encoding="utf-8")
            return CommandResult(tuple(argv), 0, "", "")
        raise AssertionError(f"unexpected command: {argv}")


class AdapterTests(unittest.TestCase):
    """The stage boundary: what run_methylation refuses before modkit ever produces a file."""

    def _fixture(self, root: Path) -> tuple[SampleManifest, AlignedBamIntakeReport]:
        bam = root / "sample.bam"
        bam.write_bytes(b"not a real bam")
        (root / "sample.bam.bai").write_bytes(b"")
        manifest = SampleManifest(
            sample_id="SYNTHETIC_001",
            run_id="RUN_001",
            input=InputSpec(
                kind=InputKind.ALIGNED_BAM, path=str(bam), index_path=str(bam) + ".bai"
            ),
            assay=AssaySpec(
                mode=AssayMode.LOW_COVERAGE_WGS,
                genome_build=GenomeBuild.GRCH38,
                reference_id="synthetic-grch38",
            ),
            analysis=AnalysisSpec(profile="synthetic", modules=[AnalysisModule.METHYLATION]),
        )
        intake = AlignedBamIntakeReport(
            sample_id="SYNTHETIC_001",
            reference_id="synthetic-grch38",
            genome_build=GenomeBuild.GRCH38,
            checks=[],
            verdict=Verdict.PASS,
        )
        return manifest, intake

    def test_a_bam_without_modified_base_tags_is_refused(self) -> None:
        """An empty pileup would read as unmethylated DNA, so the run stops before one."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._fixture(root)
            with self.assertRaises(ValueError) as raised:
                run_methylation(
                    manifest,
                    intake,
                    _policy(cpg_only=False, combine_strands=False),
                    output_dir=root / "out",
                    runner=_FakeRunner(tagged_reads="0"),
                )
        self.assertIn("no MM modified-base tags", str(raised.exception))

    def test_an_enriched_run_aggregating_over_chromosomes_says_the_bed_did_not_apply(
        self,
    ) -> None:
        """The report records the design; it must not let that read as a restriction.

        With region_source=chromosome the BED never reaches modkit, so the pileup covers the
        whole genome — but the report still carries the design's checksum, which reads as
        though the fractions came from the targets.
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._fixture(root)
            bed = root / "targets.bed"
            bed.write_text("chr1\t100\t200\tA\n", encoding="utf-8")
            manifest = manifest.model_copy(
                update={
                    "assay": manifest.assay.model_copy(
                        update={
                            "mode": AssayMode.ADAPTIVE_SAMPLING,
                            "target_bed": str(bed),
                        }
                    )
                }
            )
            report = run_methylation(
                manifest,
                intake,
                _policy(cpg_only=False, combine_strands=False),
                output_dir=root / "out",
                runner=_FakeRunner(),
            )
        self.assertIsNotNone(report.target_bed_fingerprint)
        self.assertTrue(
            any("was not restricted to that design" in item for item in report.warnings),
            report.warnings,
        )

    def test_a_version_outside_the_lock_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._fixture(root)
            with self.assertRaises(ValueError):
                run_methylation(
                    manifest,
                    intake,
                    _policy(cpg_only=False, combine_strands=False),
                    output_dir=root / "out",
                    runner=_FakeRunner(version="0.5.0"),
                )

    def test_a_cpg_restriction_without_a_reference_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._fixture(root)
            with self.assertRaises(ValueError) as raised:
                run_methylation(
                    manifest, intake, _policy(), output_dir=root / "out", runner=_FakeRunner()
                )
        self.assertIn("reference FASTA", str(raised.exception))

    def test_a_completed_pileup_is_normalized_and_records_the_tag_count(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._fixture(root)
            runner = _FakeRunner()
            report = run_methylation(
                manifest,
                intake,
                _policy(cpg_only=False, combine_strands=False),
                output_dir=root / "out",
                runner=runner,
            )
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.reads_with_modified_base_tags, 1200)
        self.assertEqual(report.tool.name, "modkit")
        pileup = next(item for item in runner.commands if "pileup" in item)
        self.assertIn("--only-tabs", pileup)
        self.assertIn("--filter-threshold", pileup)

    def test_an_unanswerable_tag_probe_warns_rather_than_blocking(self) -> None:
        """Not knowing whether tags are present is a distinct answer from knowing there are none."""

        class _NoFilterExpressions(_FakeRunner):
            def run(self, argv, *, timeout_seconds: int = 300) -> CommandResult:
                argv = [str(item) for item in argv]
                if "view" in argv:
                    self.commands.append(argv)
                    return CommandResult(tuple(argv), 1, "", "unrecognised option -e")
                return super().run(argv, timeout_seconds=timeout_seconds)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._fixture(root)
            report = run_methylation(
                manifest,
                intake,
                _policy(cpg_only=False, combine_strands=False),
                output_dir=root / "out",
                runner=_NoFilterExpressions(),
            )
        self.assertIsNone(report.reads_with_modified_base_tags)
        self.assertTrue(
            any("was not verified" in warning for warning in report.warnings), report.warnings
        )


class PolicyTests(unittest.TestCase):
    def test_a_code_cannot_be_both_reported_and_ignored(self) -> None:
        with self.assertRaises(ValueError):
            _policy(modification_codes=["m"], ignored_codes=["m"])

    def test_strand_folding_requires_a_cpg_restriction(self) -> None:
        with self.assertRaises(ValueError):
            _policy(cpg_only=False, combine_strands=True)


if __name__ == "__main__":
    unittest.main()
