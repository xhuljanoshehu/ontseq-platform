from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.cutesv import run_cutesv
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


if __name__ == "__main__":
    unittest.main()
