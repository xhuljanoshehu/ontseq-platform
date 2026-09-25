from __future__ import annotations

import tempfile
import unittest
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

from ontseq_platform.cutesv import run_cutesv
from ontseq_platform.cutesv_build import BUILD_ID, FIXED_BODY_SHA256, LAUNCH_CONTRACT
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CuteSvPolicy,
    GenomeBuild,
    InputKind,
    InputSpec,
    SampleManifest,
    Verdict,
)

VALID_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
chr1\t1000\tcute-1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1200;SVLEN=-201;RE=8\tGT:DR:DV\t0/1:8:8
"""

QUALIFIED_CUTESV_IDENTITY = {
    "cutesv_source_sha256": "1" * 64,
    "cutesv_source_body_sha256": FIXED_BODY_SHA256,
    "cutesv_source_shebang_sha256": "2" * 64,
    "cutesv_launch_contract": LAUNCH_CONTRACT,
    "cutesv_launch_interpreter": "python3",
    "cutesv_build_id": BUILD_ID,
    "cutesv_execution_body_sha256": FIXED_BODY_SHA256,
}


@contextmanager
def qualified_cutesv_mock() -> Iterator[None]:
    """Isolate mocked-runner tests from the separately tested executable qualification layer."""
    with (
        patch(
            "ontseq_platform.cutesv.executable_identity",
            return_value=QUALIFIED_CUTESV_IDENTITY.copy(),
        ),
        patch(
            "ontseq_platform.cutesv.prepare_executable",
            side_effect=lambda executable, _directory, _expected: executable,
        ),
        patch(
            "ontseq_platform.pipeline.runner.cutesv_executable_identity",
            return_value=QUALIFIED_CUTESV_IDENTITY.copy(),
        ),
    ):
        yield


class CuteSvRunner:
    def __init__(self, *, returncode: int = 0, vcf: str = VALID_VCF) -> None:
        self.returncode = returncode
        self.vcf = vcf
        self.call_argv: tuple[str, ...] | None = None
        self.staged_vcf: Path | None = None

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(str(item) for item in argv)
        if "--version" in normalized:
            return CommandResult(normalized, 0, "cuteSV 2.1.3\n", "")
        self.call_argv = normalized
        self.staged_vcf = Path(normalized[3])
        self.assert_work_directory(Path(normalized[4]))
        self.staged_vcf.write_text(self.vcf, encoding="utf-8")
        return CommandResult(normalized, self.returncode, "", "synthetic failure")

    @staticmethod
    def assert_work_directory(path: Path) -> None:
        if not path.is_dir():
            raise AssertionError("cuteSV work directory was not created before execution")


def _manifest(root: Path) -> tuple[SampleManifest, AlignedBamIntakeReport, Path]:
    bam = root / "synthetic.bam"
    bai = root / "synthetic.bam.bai"
    reference = root / "synthetic.fasta"
    bam.write_bytes(b"synthetic")
    bai.write_bytes(b"synthetic")
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    manifest = SampleManifest(
        sample_id="SYNTHETIC_CUTESV_001",
        run_id="SYNTHETIC_CUTESV_RUN_001",
        input=InputSpec(kind=InputKind.ALIGNED_BAM, path=str(bam), index_path=str(bai)),
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
    intake = AlignedBamIntakeReport(
        sample_id=manifest.sample_id,
        reference_id=manifest.assay.reference_id,
        genome_build=manifest.assay.genome_build,
        checks=[],
        verdict=Verdict.PASS,
    )
    return manifest, intake, reference


def _policy() -> CuteSvPolicy:
    return CuteSvPolicy(
        profile_id="synthetic-cutesv-atomic",
        status="technical_defaults_only",
        note="Synthetic technical thresholds only.",
    )


class CuteSvAtomicTests(unittest.TestCase):
    def setUp(self) -> None:
        qualification = qualified_cutesv_mock()
        qualification.__enter__()
        self.addCleanup(qualification.__exit__, None, None, None)

    def test_productive_call_uses_locked_parameters_and_promotes_valid_vcf(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, intake, reference = _manifest(root)
            output = root / "calls.vcf"
            runner = CuteSvRunner()
            report = run_cutesv(
                manifest,
                intake,
                _policy(),
                reference_fasta=reference,
                output_vcf=output,
                runner=runner,
                threads=2,
            )
            self.assertTrue(output.is_file())
            self.assertEqual(report.accepted_record_count, 1)
            self.assertIn("--min_support", runner.call_argv or ())
            self.assertIsNotNone(runner.staged_vcf)
            self.assertFalse(runner.staged_vcf.exists() if runner.staged_vcf else True)

    def test_nonzero_exit_never_leaves_a_final_or_staging_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, intake, reference = _manifest(root)
            output = root / "calls.vcf"
            runner = CuteSvRunner(returncode=1)
            with self.assertRaisesRegex(ValueError, "synthetic failure"):
                run_cutesv(
                    manifest,
                    intake,
                    _policy(),
                    reference_fasta=reference,
                    output_vcf=output,
                    runner=runner,
                )
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".cutesv-*")), [])

    def test_large_work_files_use_managed_scratch_and_keep_atomic_output(self) -> None:
        for outcome in ("success", "tool_error", "invalid_vcf", "runner_error"):
            with self.subTest(outcome=outcome), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                output_root = root / "results"
                scratch_root = root / "native-scratch"
                output_root.mkdir()
                scratch_root.mkdir()
                manifest, intake, reference = _manifest(root)
                output = output_root / "calls.vcf"

                class ScratchRunner(CuteSvRunner):
                    @staticmethod
                    def assert_work_directory(
                        path: Path, expected_root: Path = scratch_root, failure: str = outcome
                    ) -> None:
                        assert path.is_dir()
                        assert path.parent == expected_root
                        (path / "synthetic-signatures").write_text("synthetic")
                        if failure == "runner_error":
                            raise OSError("synthetic I/O error")

                runner = ScratchRunner(
                    returncode=1 if outcome == "tool_error" else 0,
                    vcf="invalid" if outcome == "invalid_vcf" else VALID_VCF,
                )
                with patch.dict("os.environ", {"ONTSEQ_CUTESV_SCRATCH_ROOT": str(scratch_root)}):
                    if outcome == "success":
                        report = run_cutesv(
                            manifest,
                            intake,
                            _policy(),
                            reference_fasta=reference,
                            output_vcf=output,
                            runner=runner,
                        )
                        self.assertEqual(report.accepted_record_count, 1)
                        self.assertTrue(output.is_file())
                        self.assertEqual(
                            report.tool.parameters["scratch_storage"], "managed_directory"
                        )
                    else:
                        with self.assertRaises((ValueError, OSError)):
                            run_cutesv(
                                manifest,
                                intake,
                                _policy(),
                                reference_fasta=reference,
                                output_vcf=output,
                                runner=runner,
                            )
                        self.assertFalse(output.exists())
                self.assertEqual(list(scratch_root.iterdir()), [])
                self.assertEqual(list(output_root.glob(".cutesv-*")), [])
                self.assertIsNotNone(runner.staged_vcf)
                if runner.staged_vcf is not None:
                    self.assertEqual(runner.staged_vcf.parent.parent, output_root)

    def test_default_scratch_uses_user_cache_instead_of_ram_temp(self) -> None:
        import os

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, intake, reference = _manifest(root)
            runner = CuteSvRunner()
            with (
                patch.dict(os.environ),
                patch("ontseq_platform.cutesv.Path.home", return_value=root),
            ):
                os.environ.pop("ONTSEQ_CUTESV_SCRATCH_ROOT", None)
                report = run_cutesv(
                    manifest,
                    intake,
                    _policy(),
                    reference_fasta=reference,
                    output_vcf=root / "calls.vcf",
                    runner=runner,
                )
            expected = root / ".cache" / "ontseq" / "cutesv"
            self.assertEqual(report.tool.parameters["scratch_root"], str(expected))
            self.assertEqual(list(expected.iterdir()), [])

    def test_relative_scratch_override_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, intake, reference = _manifest(root)
            with (
                patch.dict("os.environ", {"ONTSEQ_CUTESV_SCRATCH_ROOT": "relative"}),
                self.assertRaisesRegex(ValueError, "absolute"),
            ):
                run_cutesv(
                    manifest,
                    intake,
                    _policy(),
                    reference_fasta=reference,
                    output_vcf=root / "calls.vcf",
                    runner=CuteSvRunner(),
                )
            self.assertFalse((root / "calls.vcf").exists())


if __name__ == "__main__":
    unittest.main()
