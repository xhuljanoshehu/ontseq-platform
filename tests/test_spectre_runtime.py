from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cnv.spectre import (
    SpectrePolicy,
    build_spectre_depth_only_argv,
    normalize_spectre_vcf,
    spectre_version,
)
from ontseq_platform.models import EventType, GenomeBuild, ModuleRunStatus, ToolRecord


def _policy() -> SpectrePolicy:
    return SpectrePolicy(
        profile_id="spectre-depth-only-v1",
        expected_version="0.2.1",
        genome_build=GenomeBuild.GRCH38,
        mosdepth_bin_size_bp=1000,
        minimum_mapping_quality=20,
        minimum_cnv_length_bp=100_000,
        ploidy=2,
        timeout_seconds=3600,
        note="Synthetic Spectre depth-only technical policy.",
    )


def _tool(version: str = "0.2.1") -> ToolRecord:
    return ToolRecord(name="Spectre", version=version, parameters={})


def _write_vcf(path: Path, *records: str) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=Spectre\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n"
        + "".join(f"{record}\n" for record in records),
        encoding="utf-8",
    )


class SpectreDepthOnlyContractTests(unittest.TestCase):
    def test_version_parser_normalizes_current_patch_banner_to_pinned_semver(self) -> None:
        self.assertEqual(
            spectre_version("INFO Spectre version: 0.2.1-patch-july15"),
            "0.2.1",
        )

    def test_depth_only_argv_uses_required_inputs_without_support_callers(self) -> None:
        argv = build_spectre_depth_only_argv(
            spectre="spectre",
            coverage_path=Path("/synthetic/sample.regions.bed.gz"),
            sample_id="SYNTHETIC_001",
            output_dir=Path("/synthetic/out"),
            reference_fasta=Path("/synthetic/reference.fa"),
            policy=_policy(),
            threads=4,
        )

        self.assertEqual(argv[:2], ("spectre", "CNVCaller"))
        self.assertIn("--coverage", argv)
        self.assertIn("/synthetic/sample.regions.bed.gz", argv)
        self.assertIn("--sample-id", argv)
        self.assertIn("SYNTHETIC_001", argv)
        self.assertIn("--output-dir", argv)
        self.assertIn("--reference", argv)
        self.assertIn("--min-cnv-len", argv)
        self.assertIn("100000", argv)
        self.assertIn("--ploidy", argv)
        self.assertIn("2", argv)
        self.assertIn("--threads", argv)
        self.assertIn("4", argv)
        self.assertNotIn("--snfj", argv)
        self.assertNotIn("--snv", argv)
        self.assertNotIn("--population", argv)
        self.assertNotIn("--cancer", argv)

    def test_spectre_native_zero_based_coordinates_are_not_shifted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectre.vcf"
            _write_vcf(
                path,
                (
                    "chr7\t100000\tSpectre.DEL.1\tN\t<DEL>\t35\tPASS\t"
                    "END=300000;SVLEN=200000;SVTYPE=DEL;CN=1;\tGT:HO:GQ:CD\t0/1:0.4:35:10"
                ),
                (
                    "chr8\t500000\tSpectre.DUP.2\tN\t<DUP>\t42\tPASS\t"
                    "END=750000;SVLEN=250000;SVTYPE=DUP;CN=3;\tGT:HO:GQ:CD\t0/1:0.5:42:15"
                ),
            )

            report = normalize_spectre_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=_tool(),
            )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.raw_record_count, 2)
        self.assertEqual(report.accepted_record_count, 2)
        self.assertEqual(report.rejected_record_count, 0)
        self.assertEqual(
            [event.event_type for event in report.events],
            [EventType.DELETION, EventType.DUPLICATION],
        )
        deletion = report.events[0]
        self.assertEqual(deletion.primary.start, 100000)
        self.assertEqual(deletion.primary.end, 300000)
        self.assertEqual(deletion.length_bp, 200000)
        self.assertEqual(deletion.copy_number, 1.0)
        duplication = report.events[1]
        self.assertEqual(duplication.primary.start, 500000)
        self.assertEqual(duplication.primary.end, 750000)
        self.assertEqual(duplication.copy_number, 3.0)

    def test_empty_complete_vcf_is_explicit_technical_no_call(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectre.vcf"
            _write_vcf(path)

            report = normalize_spectre_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=_tool(),
            )

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.raw_record_count, 0)
        self.assertEqual(report.accepted_record_count, 0)
        self.assertTrue(any("not a biological negative" in item.lower() for item in report.warnings))

    def test_nonfinite_copy_number_is_rejected_not_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectre.vcf"
            _write_vcf(
                path,
                (
                    "chr7\t100000\tSpectre.DEL.1\tN\t<DEL>\t35\tPASS\t"
                    "END=300000;SVLEN=200000;SVTYPE=DEL;CN=nan;\tGT:HO:GQ:CD\t0/1:0.4:35:10"
                ),
            )

            report = normalize_spectre_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=_tool(),
            )

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.raw_record_count, 1)
        self.assertEqual(report.accepted_record_count, 0)
        self.assertEqual(report.rejected_record_count, 1)
        self.assertEqual(report.rejection_counts, {"nonfinite_copy_number": 1})

    def test_build_mismatch_fails_closed_before_normalization(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectre.vcf"
            _write_vcf(path)

            with self.assertRaisesRegex(ValueError, "genome build"):
                normalize_spectre_vcf(
                    path,
                    sample_id="SYNTHETIC_001",
                    genome_build=GenomeBuild.GRCH37,
                    policy=_policy(),
                    tool=_tool(),
                )

    def test_runtime_version_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectre.vcf"
            _write_vcf(path)

            with self.assertRaisesRegex(ValueError, "version"):
                normalize_spectre_vcf(
                    path,
                    sample_id="SYNTHETIC_001",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                    tool=_tool("0.2.0"),
                )


if __name__ == "__main__":
    unittest.main()
