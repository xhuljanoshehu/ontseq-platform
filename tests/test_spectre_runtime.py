from __future__ import annotations

import gzip
import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

import yaml

from ontseq_platform.cnv import spectre as spectre_module
from ontseq_platform.cnv.spectre import (
    SpectrePolicy,
    build_spectre_depth_only_argv,
    normalize_spectre_vcf,
    spectre_version,
)
from ontseq_platform.execution import CommandResult
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


def _write_vcf(
    path: Path,
    *records: str,
    sample_id: str = "SYNTHETIC_001",
    source: str = "Spectre",
) -> None:
    payload = (
        "##fileformat=VCFv4.2\n"
        f"##source={source}\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n"
        + "".join(f"{record}\n" for record in records)
    )
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(payload)
    else:
        path.write_text(payload, encoding="utf-8")


def _runtime_inputs(root: Path) -> tuple[Path, Path]:
    coverage = root / "synthetic.regions.bed.gz"
    with gzip.open(coverage, "wt", encoding="utf-8") as handle:
        handle.write("chr7\t0\t1000\t30\nchr7\t1000\t2000\t30\nchr7\t2000\t3000\t15\n")
    Path(f"{coverage}.csi").write_bytes(b"synthetic-csi")
    reference = root / "synthetic.reference.fa"
    reference.write_text(">chr7\n" + "A" * 80 + "\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr7\t80\t6\t80\t81\n", encoding="utf-8")
    return coverage, reference


class _FakeSpectreRunner:
    def __init__(
        self,
        *,
        version: str = "0.2.1",
        records: tuple[str, ...] | None = None,
        output_suffix: str = ".vcf.gz",
        fail_execution: bool = False,
    ) -> None:
        self.version = version
        self.records = (
            records
            if records is not None
            else (
                "chr7\t100000\tSpectre.DEL.1\tN\t<DEL>\t.\t.\t"
                "END=300000;SVLEN=200000;SVTYPE=DEL;CN=1\tGT:HO:GQ:DP\t0/1:0.4:35:10",
            )
        )
        self.output_suffix = output_suffix
        self.fail_execution = fail_execution
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        args = tuple(str(item) for item in argv)
        self.calls.append(args)
        if args == ("spectre", "version"):
            return CommandResult(
                args,
                0,
                f"Spectre input: version\nSpectre version: {self.version}\n",
                "",
            )

        output_dir = Path(args[args.index("--output-dir") + 1])
        sample_id = args[args.index("--sample-id") + 1]
        output_dir.mkdir(parents=True, exist_ok=True)
        if self.fail_execution:
            (output_dir / "partial.tmp").write_text("partial", encoding="utf-8")
            return CommandResult(args, 9, "", "synthetic Spectre failure")

        _write_vcf(
            output_dir / f"{sample_id}{self.output_suffix}",
            *self.records,
            sample_id=sample_id,
        )
        (output_dir / f"{sample_id}_cnv.bed.gz").write_bytes(b"synthetic-bed")
        (output_dir / f"{sample_id}_cnv.bed.gz.tbi").write_bytes(b"synthetic-index")
        (output_dir / f"{sample_id}.spc.gz").write_bytes(b"synthetic-spc")
        (output_dir / "metadata.mdr").write_text("synthetic metadata\n", encoding="utf-8")
        image_dir = output_dir / "img"
        image_dir.mkdir()
        (image_dir / f"{sample_id}_plot_cnv_chr7.png").write_bytes(b"synthetic-png")
        return CommandResult(args, 0, "Spectre finished\n", "")


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

    def test_native_unfiltered_final_call_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectre.vcf"
            _write_vcf(
                path,
                (
                    "chr7\t100000\tSpectre.DEL.1\tN\t<DEL>\t.\t.\t"
                    "END=300000;SVLEN=200000;SVTYPE=DEL;CN=1\tGT:HO:GQ:DP\t0/1:0.4:35:10"
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
        self.assertEqual(report.accepted_record_count, 1)
        self.assertEqual(report.events[0].evidence[0].filters, [])

    def test_malformed_coordinates_are_rejected_without_partial_normalization(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectre.vcf"
            _write_vcf(
                path,
                (
                    "chr7\tnan\tSpectre.DEL.1\tN\t<DEL>\t35\tPASS\t"
                    "END=300000;SVLEN=200000;SVTYPE=DEL;CN=1\tGT:HO:GQ:DP\t0/1:0.4:35:10"
                ),
                (
                    "chr7\t100000\tSpectre.DEL.2\tN\t<DEL>\t35\tPASS\t"
                    "END=100000;SVLEN=0;SVTYPE=DEL;CN=1\tGT:HO:GQ:DP\t0/1:0.4:35:10"
                ),
                (
                    "chr7\t100000\tSpectre.DEL.3\tN\t<DEL>\t35\tPASS\t"
                    "END=300000;SVLEN=199999;SVTYPE=DEL;CN=1\tGT:HO:GQ:DP\t0/1:0.4:35:10"
                ),
                (
                    "chr7\t100000\tSpectre.DEL.4\tN\t<DEL>\tinf\tPASS\t"
                    "END=300000;SVLEN=200000;SVTYPE=DEL;CN=1\tGT:HO:GQ:DP\t0/1:0.4:35:10"
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
        self.assertEqual(report.raw_record_count, 4)
        self.assertEqual(report.accepted_record_count, 0)
        self.assertEqual(
            report.rejection_counts,
            {
                "invalid_native_interval": 1,
                "malformed_quality": 1,
                "malformed_start": 1,
                "native_length_mismatch": 1,
            },
        )

    def test_sample_header_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectre.vcf"
            _write_vcf(path, sample_id="OTHER_SAMPLE")

            with self.assertRaisesRegex(ValueError, "sample"):
                normalize_spectre_vcf(
                    path,
                    sample_id="SYNTHETIC_001",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                    tool=_tool(),
                )

    def test_non_spectre_vcf_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "other.vcf"
            _write_vcf(path, source="OTHER_CALLER")

            with self.assertRaisesRegex(ValueError, "Spectre"):
                normalize_spectre_vcf(
                    path,
                    sample_id="SYNTHETIC_001",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                    tool=_tool(),
                )

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
        self.assertTrue(
            any("not a biological negative" in item.lower() for item in report.warnings)
        )

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

    def test_pinned_technical_policy_is_loadable_and_real_tool_qualified(self) -> None:
        path = Path(__file__).parents[1] / "configs" / "cnv" / "spectre.depth-only.technical.yaml"
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
        policy = SpectrePolicy.model_validate(payload)

        self.assertEqual(policy.expected_version, "0.2.1")
        self.assertEqual(policy.mode, "depth_only")
        self.assertEqual(policy.runtime_package, "spectre-cnv==0.2.1")
        self.assertEqual(
            policy.runtime_wheel_sha256,
            "c50998798a33b22455f3681a69c5c61e82df3fe049224da6d2e54b068bc35510",
        )
        self.assertTrue(policy.real_tool_qualified)
        self.assertEqual(policy.analytical_validation, "not_validated")

    def test_policy_rejects_version_that_disagrees_with_pinned_runtime(self) -> None:
        payload = _policy().model_dump(mode="json")
        payload["expected_version"] = "0.2.0"

        with self.assertRaisesRegex(ValueError, "0.2.1"):
            SpectrePolicy.model_validate(payload)

    def test_runtime_depth_only_preserves_native_artifacts_without_support_callers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage, reference = _runtime_inputs(root)
            output_dir = root / "spectre-output"
            runner = _FakeSpectreRunner()

            report = spectre_module.run_spectre_depth_only(
                coverage_path=coverage,
                sample_id="SYNTHETIC_001",
                output_dir=output_dir,
                reference_fasta=reference,
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                runner=runner,
                threads=4,
            )
            self.assertTrue((output_dir / "SYNTHETIC_001.vcf.gz").is_file())

        self.assertEqual(runner.calls[0], ("spectre", "version"))
        argv = runner.calls[1]
        for forbidden in ("--snfj", "--snv", "--population", "--cancer"):
            self.assertNotIn(forbidden, argv)
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.tool.version, "0.2.1")
        self.assertEqual(report.tool.parameters["mode"], "depth_only")
        self.assertFalse(report.tool.parameters["snfj_input"])
        self.assertFalse(report.tool.parameters["snv_input"])
        self.assertFalse(report.tool.parameters["population_mode"])
        self.assertFalse(report.tool.parameters["cancer_mode"])
        self.assertEqual(
            {artifact.relative_path for artifact in report.native_artifacts},
            {
                "SYNTHETIC_001.spc.gz",
                "SYNTHETIC_001.vcf.gz",
                "SYNTHETIC_001_cnv.bed.gz",
                "SYNTHETIC_001_cnv.bed.gz.tbi",
                "img/SYNTHETIC_001_plot_cnv_chr7.png",
                "metadata.mdr",
            },
        )
        self.assertTrue(all(artifact.fingerprint.sha256 for artifact in report.native_artifacts))

    def test_runtime_preserves_native_empty_vcf_as_technical_no_call(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage, reference = _runtime_inputs(root)
            output_dir = root / "spectre-output"
            runner = _FakeSpectreRunner(records=(), output_suffix=".vcf")

            report = spectre_module.run_spectre_depth_only(
                coverage_path=coverage,
                sample_id="SYNTHETIC_001",
                output_dir=output_dir,
                reference_fasta=reference,
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                runner=runner,
            )
            self.assertTrue((output_dir / "SYNTHETIC_001.vcf").is_file())

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.raw_record_count, 0)
        self.assertTrue(
            any("not a biological negative" in warning.lower() for warning in report.warnings)
        )

    def test_runtime_failure_never_promotes_partial_native_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage, reference = _runtime_inputs(root)
            output_dir = root / "spectre-output"
            runner = _FakeSpectreRunner(fail_execution=True)

            with self.assertRaisesRegex(ValueError, "exit code 9"):
                spectre_module.run_spectre_depth_only(
                    coverage_path=coverage,
                    sample_id="SYNTHETIC_001",
                    output_dir=output_dir,
                    reference_fasta=reference,
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                    runner=runner,
                )

            self.assertFalse(output_dir.exists())

    def test_runtime_version_mismatch_stops_before_caller_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage, reference = _runtime_inputs(root)
            runner = _FakeSpectreRunner(version="0.2.0")

            with self.assertRaisesRegex(ValueError, "version"):
                spectre_module.run_spectre_depth_only(
                    coverage_path=coverage,
                    sample_id="SYNTHETIC_001",
                    output_dir=root / "spectre-output",
                    reference_fasta=reference,
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                    runner=runner,
                )

        self.assertEqual(runner.calls, [("spectre", "version")])


if __name__ == "__main__":
    unittest.main()
