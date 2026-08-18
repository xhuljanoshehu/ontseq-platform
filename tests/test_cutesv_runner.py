from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.cutesv_runner import CuteSVExecutionPolicy, run_cutesv
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    BamHeaderSummary,
    GenomeBuild,
    InputKind,
    InputSpec,
    SampleManifest,
    Verdict,
)
from ontseq_platform.reference import contig_signature, sha256_file

VCF_TEXT = """##fileformat=VCFv4.2
##source=cuteSV-2.1.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_001
chr1\t1000\tSECRET_ID\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1200;SVLEN=-201;RE=8;AF=0.4\tGT:DR:DV\t0/1:12:8
"""


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
            modules=[AnalysisModule.SV, AnalysisModule.REPORT],
        ),
    )


def _header_signature() -> str:
    return contig_signature([("chr1", 10000), ("chr2", 20000)])


def _intake(*, signature: str | None = None) -> AlignedBamIntakeReport:
    return AlignedBamIntakeReport(
        sample_id="SYNTHETIC_001",
        reference_id="SYNTHETIC_REF",
        genome_build=GenomeBuild.GRCH38,
        header=BamHeaderSummary(
            sort_order="coordinate",
            sequence_count=2,
            total_reference_bases=30000,
            contig_signature_sha256=signature or _header_signature(),
            read_group_count=0,
            sample_tag_count=0,
            program_count=1,
        ),
        checks=[],
        verdict=Verdict.PASS,
    )


def _write_reference(directory: Path) -> tuple[Path, Path]:
    fasta = directory / "synthetic.fa"
    fai = directory / "synthetic.fa.fai"
    fasta.write_text(">chr1\nA\n>chr2\nA\n", encoding="utf-8")
    fai.write_text("chr1\t10000\t0\t80\t81\nchr2\t20000\t0\t80\t81\n", encoding="utf-8")
    return fasta, fai


class RecordingCuteSVRunner:
    def __init__(self, *, version: str = "2.1.4", runtime_returncode: int = 0) -> None:
        self.version = version
        self.runtime_returncode = runtime_returncode
        self.commands: list[tuple[str, ...]] = []
        self.work_dir_existed = False

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(str(item) for item in argv)
        self.commands.append(normalized)
        if "--version" in normalized:
            return CommandResult(normalized, 0, f"cuteSV {self.version}\n", "")

        output = Path(normalized[3])
        work_dir = Path(normalized[4])
        self.work_dir_existed = work_dir.is_dir()
        output.write_text(VCF_TEXT, encoding="utf-8")
        return CommandResult(normalized, self.runtime_returncode, "", "")


class CuteSVRunnerTests(unittest.TestCase):
    def test_shell_free_execution_is_privacy_safe_and_ont_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "synthetic.bam"
            bai = root / "synthetic.bam.bai"
            bam.write_bytes(b"synthetic")
            bai.write_bytes(b"synthetic-index")
            fasta, fai = _write_reference(root)
            output = root / "result.vcf"
            runner = RecordingCuteSVRunner()

            report = run_cutesv(
                _manifest(bam, bai),
                _intake(),
                CuteSVExecutionPolicy(),
                reference_fasta=fasta,
                output_vcf=output,
                runner=runner,
                threads=3,
            )

            self.assertEqual(len(runner.commands), 2)
            self.assertEqual(runner.commands[0], ("cuteSV", "--version"))
            command = runner.commands[1]
            self.assertEqual(command[0], "cuteSV")
            self.assertEqual(command[1], str(bam))
            self.assertEqual(command[2], str(fasta))
            self.assertEqual(command[3], str(output))
            self.assertTrue(runner.work_dir_existed)
            self.assertIn("--ignore_sequence", command)
            self.assertIn("--min_mapq", command)
            self.assertIn("--min_support", command)
            self.assertIn("--min_size", command)
            self.assertIn("--max_size", command)
            self.assertIn("--max_cluster_bias_INS", command)
            self.assertIn("--diff_ratio_merging_INS", command)
            self.assertIn("--max_cluster_bias_DEL", command)
            self.assertIn("--diff_ratio_merging_DEL", command)
            self.assertIn("--max_cluster_bias_TRA", command)
            self.assertIn("--diff_ratio_filtering_TRA", command)
            for forbidden_flag in [
                "--report_readid",
                "--retain_work_dir",
                "--write_old_sigs",
                "--genotype",
            ]:
                self.assertNotIn(forbidden_flag, command)

            self.assertEqual(report.tool.version, "2.1.4")
            self.assertEqual(report.tool.parameters["threads"], 3)
            self.assertEqual(report.tool.parameters["max_size"], -1)
            self.assertEqual(report.tool.parameters["report_read_ids"], False)
            self.assertEqual(report.tool.parameters["ignore_sequence"], True)
            self.assertEqual(report.tool.parameters["reference_fai_sha256"], sha256_file(fai))
            self.assertEqual(
                report.tool.parameters["reference_contig_signature_sha256"],
                _header_signature(),
            )

            serialized = report.model_dump_json()
            for forbidden_path in [str(bam), str(bai), str(fasta), str(fai), str(output), directory]:
                self.assertNotIn(forbidden_path, serialized)
            self.assertNotIn("SECRET_ID", serialized)

    def test_reference_contig_mismatch_fails_before_tool_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "synthetic.bam"
            bai = root / "synthetic.bam.bai"
            bam.write_bytes(b"synthetic")
            bai.write_bytes(b"synthetic-index")
            fasta, _fai = _write_reference(root)
            output = root / "result.vcf"
            runner = RecordingCuteSVRunner()

            with self.assertRaisesRegex(ValueError, "contig signature does not match"):
                run_cutesv(
                    _manifest(bam, bai),
                    _intake(signature=contig_signature([("chr1", 9999)])),
                    CuteSVExecutionPolicy(),
                    reference_fasta=fasta,
                    output_vcf=output,
                    runner=runner,
                )

            self.assertEqual(runner.commands, [])
            self.assertFalse(output.exists())

    def test_version_mismatch_fails_closed_before_runtime_call(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "synthetic.bam"
            bai = root / "synthetic.bam.bai"
            bam.write_bytes(b"synthetic")
            bai.write_bytes(b"synthetic-index")
            fasta, _fai = _write_reference(root)
            output = root / "result.vcf"
            runner = RecordingCuteSVRunner(version="2.1.3")

            with self.assertRaisesRegex(ValueError, "does not match policy lock"):
                run_cutesv(
                    _manifest(bam, bai),
                    _intake(),
                    CuteSVExecutionPolicy(),
                    reference_fasta=fasta,
                    output_vcf=output,
                    runner=runner,
                )

            self.assertEqual(len(runner.commands), 1)
            self.assertFalse(output.exists())

    def test_nonzero_runtime_removes_partial_vcf(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "synthetic.bam"
            bai = root / "synthetic.bam.bai"
            bam.write_bytes(b"synthetic")
            bai.write_bytes(b"synthetic-index")
            fasta, _fai = _write_reference(root)
            output = root / "result.vcf"
            runner = RecordingCuteSVRunner(runtime_returncode=2)

            with self.assertRaisesRegex(ValueError, "non-zero exit code"):
                run_cutesv(
                    _manifest(bam, bai),
                    _intake(),
                    CuteSVExecutionPolicy(),
                    reference_fasta=fasta,
                    output_vcf=output,
                    runner=runner,
                )

            self.assertFalse(output.exists())

    def test_existing_output_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "synthetic.bam"
            bai = root / "synthetic.bam.bai"
            bam.write_bytes(b"synthetic")
            bai.write_bytes(b"synthetic-index")
            fasta, _fai = _write_reference(root)
            output = root / "result.vcf"
            output.write_text("KEEP", encoding="utf-8")
            runner = RecordingCuteSVRunner()

            with self.assertRaisesRegex(ValueError, "Refusing to overwrite"):
                run_cutesv(
                    _manifest(bam, bai),
                    _intake(),
                    CuteSVExecutionPolicy(),
                    reference_fasta=fasta,
                    output_vcf=output,
                    runner=runner,
                )

            self.assertEqual(output.read_text(encoding="utf-8"), "KEEP")
            self.assertEqual(runner.commands, [])


if __name__ == "__main__":
    unittest.main()
