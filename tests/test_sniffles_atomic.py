from __future__ import annotations

import hashlib
import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.execution import CommandResult, ToolExecutionError
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    SampleManifest,
    SnifflesPolicy,
    Verdict,
)
from ontseq_platform.sniffles import run_sniffles

VCF_HEADER = """##fileformat=VCFv4.2
##source=Sniffles2_2.8.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
"""
VALID_VCF = (
    VCF_HEADER
    + "chr1\t1000\t.\tN\t<DEL>\t60\tPASS\t"
    "SVTYPE=DEL;END=1200;SVLEN=-201;SUPPORT=8\tGT:DR:DV\t0/1:8:8\n"
)


def _policy() -> SnifflesPolicy:
    return SnifflesPolicy(
        profile_id="atomic-finalization-test",
        status="technical_defaults_only",
        min_support=5,
        min_sv_length=50,
        mapq=20,
        note="Synthetic technical thresholds only.",
    )


def _manifest(bam: Path, bai: Path) -> SampleManifest:
    return SampleManifest(
        sample_id="SYNTHETIC_001",
        run_id="SYNTHETIC_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path=str(bam),
            index_path=str(bai),
        ),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="SYNTHETIC_REF",
        ),
        analysis=AnalysisSpec(
            profile="synthetic",
            modules=[AnalysisModule.QC, AnalysisModule.SV, AnalysisModule.REPORT],
        ),
    )


def _intake() -> AlignedBamIntakeReport:
    return AlignedBamIntakeReport(
        sample_id="SYNTHETIC_001",
        reference_id="SYNTHETIC_REF",
        genome_build=GenomeBuild.GRCH38,
        checks=[],
        verdict=Verdict.PASS,
    )


class ScenarioRunner:
    def __init__(
        self,
        vcf_text: str,
        *,
        returncode: int = 0,
        stderr: str = "",
        raise_after_write: bool = False,
    ) -> None:
        self.vcf_text = vcf_text
        self.returncode = returncode
        self.stderr = stderr
        self.raise_after_write = raise_after_write
        self.output_paths: list[Path] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(str(item) for item in argv)
        if "--version" in normalized:
            return CommandResult(normalized, 0, "Sniffles2, Version 2.8.0\n", "")
        output = Path(normalized[normalized.index("--vcf") + 1])
        self.output_paths.append(output)
        output.write_text(self.vcf_text, encoding="utf-8")
        if self.raise_after_write:
            raise ToolExecutionError("synthetic timeout after partial VCF write")
        return CommandResult(normalized, self.returncode, "", self.stderr)


class SnifflesAtomicFinalizationTests(unittest.TestCase):
    def _inputs(self, root: Path) -> tuple[SampleManifest, AlignedBamIntakeReport]:
        bam = root / "synthetic.bam"
        bai = root / "synthetic.bam.bai"
        bam.write_bytes(b"synthetic")
        bai.write_bytes(b"synthetic")
        return _manifest(bam, bai), _intake()

    def _assert_no_staged_files(self, root: Path, output: Path) -> None:
        self.assertEqual(list(root.glob(f".{output.name}.*")), [])

    def test_partial_nonzero_exit_leaves_no_final_or_temp_and_retry_succeeds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._inputs(root)
            output = root / "calls.vcf"
            failing = ScenarioRunner(
                VCF_HEADER + "chr1\t1000\tPARTIAL",
                returncode=1,
                stderr="synthetic caller failure",
            )
            with self.assertRaisesRegex(ValueError, "synthetic caller failure"):
                run_sniffles(
                    manifest,
                    intake,
                    _policy(),
                    output_vcf=output,
                    runner=failing,
                )
            self.assertFalse(output.exists())
            self._assert_no_staged_files(root, output)

            succeeding = ScenarioRunner(VALID_VCF)
            report = run_sniffles(
                manifest,
                intake,
                _policy(),
                output_vcf=output,
                runner=succeeding,
            )
            self.assertTrue(output.is_file())
            self.assertEqual(report.accepted_record_count, 1)
            self._assert_no_staged_files(root, output)

    def test_success_exit_with_malformed_vcf_is_never_promoted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._inputs(root)
            output = root / "calls.vcf"
            runner = ScenarioRunner("##fileformat=VCFv4.2\n")
            with self.assertRaisesRegex(ValueError, "complete VCF document"):
                run_sniffles(
                    manifest,
                    intake,
                    _policy(),
                    output_vcf=output,
                    runner=runner,
                )
            self.assertFalse(output.exists())
            self._assert_no_staged_files(root, output)

    def test_exception_after_partial_write_cleans_staging_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._inputs(root)
            output = root / "calls.vcf"
            runner = ScenarioRunner(VCF_HEADER, raise_after_write=True)
            with self.assertRaisesRegex(ToolExecutionError, "synthetic timeout"):
                run_sniffles(
                    manifest,
                    intake,
                    _policy(),
                    output_vcf=output,
                    runner=runner,
                )
            self.assertFalse(output.exists())
            self._assert_no_staged_files(root, output)

    def test_valid_vcf_is_promoted_only_after_normalization_and_fingerprint_matches(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, intake = self._inputs(root)
            output = root / "calls.vcf"
            runner = ScenarioRunner(VALID_VCF)
            report = run_sniffles(
                manifest,
                intake,
                _policy(),
                output_vcf=output,
                runner=runner,
                threads=2,
            )
            self.assertTrue(output.is_file())
            self.assertEqual(output.read_text(encoding="utf-8"), VALID_VCF)
            self.assertEqual(len(runner.output_paths), 1)
            self.assertNotEqual(runner.output_paths[0], output)
            self.assertFalse(runner.output_paths[0].exists())
            self.assertEqual(
                report.vcf_fingerprint.sha256,
                hashlib.sha256(output.read_bytes()).hexdigest(),
            )
            self.assertEqual(report.vcf_fingerprint.size_bytes, output.stat().st_size)
            self._assert_no_staged_files(root, output)


if __name__ == "__main__":
    unittest.main()
