from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from pydantic import ValidationError

from ontseq_platform.execution import CommandResult
from ontseq_platform.models import GenomeBuild, ModuleRunStatus
from ontseq_platform.multicaller_contracts import CallerInputRole
from ontseq_platform.tumor.wakhan import (
    WakhanPhasedCnaPolicy,
    build_wakhan_argv,
    run_wakhan_phased_cna,
)
from ontseq_platform.tumor_inputs import (
    TumorAuxiliaryInputArtifact,
    TumorInputBundle,
)


def _sha_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _write_bam(path: Path, payload: bytes) -> None:
    path.write_bytes(payload)
    Path(f"{path}.bai").write_bytes(b"index:" + payload)


def _write_reference(path: Path) -> str:
    path.write_text(">chr7\n" + "A" * 4000 + "\n", encoding="utf-8")
    Path(f"{path}.fai").write_text("chr7\t4000\t6\t4000\t4001\n", encoding="utf-8")
    return _sha_bytes(path.read_bytes())


def _write_phased_vcf(path: Path) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr7\t500\trs1\tA\tG\t60\tPASS\t.\tGT\t0|1\n",
        encoding="utf-8",
    )


def _write_breakpoints(path: Path) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n##source=Severus\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )


def _artifact(
    *,
    artifact_id: str,
    role: CallerInputRole,
    path: Path,
    analysis_sample_id: str = "TUMOR_001",
    source_sample_id: str = "TUMOR_001",
    reference_sha256: str,
    producer_lane_id: str | None = None,
    producer_lane_sha256: str | None = None,
) -> TumorAuxiliaryInputArtifact:
    return TumorAuxiliaryInputArtifact(
        artifact_id=artifact_id,
        role=role,
        analysis_sample_id=analysis_sample_id,
        source_sample_id=source_sample_id,
        genome_build=GenomeBuild.GRCH38,
        reference_id="GRCh38-test",
        reference_sha256=reference_sha256,
        sha256=_sha_bytes(path.read_bytes()),
        producer_lane_id=producer_lane_id,
        producer_lane_sha256=producer_lane_sha256,
    )


def _bundle(
    *,
    tumor_bam: Path,
    phased_vcf: Path | None,
    reference_sha256: str,
    phased_source_sample_id: str = "NORMAL_001",
    breakpoints_vcf: Path | None = None,
    breakpoint_parent_id: str | None = None,
    breakpoint_parent_sha256: str | None = None,
) -> TumorInputBundle:
    artifacts = [
        _artifact(
            artifact_id="tumor-bam",
            role=CallerInputRole.TUMOR_BAM,
            path=tumor_bam,
            reference_sha256=reference_sha256,
        )
    ]
    if phased_vcf is not None:
        artifacts.append(
            _artifact(
                artifact_id="phased-vcf",
                role=CallerInputRole.PHASED_VARIANTS,
                path=phased_vcf,
                source_sample_id=phased_source_sample_id,
                reference_sha256=reference_sha256,
            )
        )
    if breakpoints_vcf is not None:
        artifacts.append(
            _artifact(
                artifact_id="severus-breakpoints",
                role=CallerInputRole.SEVERUS_BREAKPOINTS,
                path=breakpoints_vcf,
                reference_sha256=reference_sha256,
                producer_lane_id=breakpoint_parent_id,
                producer_lane_sha256=breakpoint_parent_sha256,
            )
        )
    return TumorInputBundle(
        analysis_sample_id="TUMOR_001",
        genome_build=GenomeBuild.GRCH38,
        reference_id="GRCh38-test",
        reference_sha256=reference_sha256,
        artifacts=artifacts,
    )


def _policy(reference_sha256: str, *, mode: str) -> WakhanPhasedCnaPolicy:
    return WakhanPhasedCnaPolicy(
        profile_id=f"wakhan-{mode}-test",
        mode=mode,
        expected_version="0.4.4",
        genome_build=GenomeBuild.GRCH38,
        reference_id="GRCh38-test",
        reference_sha256=reference_sha256,
        timeout_seconds=300,
        note="Synthetic Wakhan phased CNA adapter contract.",
    )


class FakeWakhanRunner:
    def __init__(self) -> None:
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        del timeout_seconds
        args = tuple(str(item) for item in argv)
        self.calls.append(args)
        outdir = Path(args[args.index("--out-dir-plots") + 1])
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "solutions_ranks.tsv").write_text(
            "rank\tploidy\tpurity\tconfidence\n1\t2.0\t0.70\t0.90\n",
            encoding="utf-8",
        )
        solution = outdir / "solution_2.0_0.70_0.90"
        solution.mkdir()
        (solution / "integer_profile.bed").write_text(
            "chrom\tstart\tend\thp1_cn\thp2_cn\nchr7\t0\t1000000\t1\t2\n",
            encoding="utf-8",
        )
        coverage = outdir / "coverage_data"
        coverage.mkdir()
        (coverage / "coverage.csv").write_text(
            "chr,start,end,hp1,hp2\nchr7,0,50000,10,20\n",
            encoding="utf-8",
        )
        phasing = outdir / "phasing_output"
        phasing.mkdir()
        (phasing / "rephased.vcf.gz").write_bytes(b"synthetic-rephased-vcf")
        return CommandResult(args, 0, "", "")


class WakhanRuntimeTests(unittest.TestCase):
    def test_current_published_0_4_4_policy_is_accepted(self) -> None:
        policy = WakhanPhasedCnaPolicy(
            profile_id="wakhan-current-version-contract",
            mode="tumor_only",
            expected_version="0.4.4",
            genome_build=GenomeBuild.GRCH38,
            reference_id="GRCh38-test",
            reference_sha256=_sha_bytes(b"reference"),
            timeout_seconds=300,
            note="Current published Wakhan version contract.",
        )

        self.assertEqual(policy.expected_version, "0.4.4")

    def test_obsolete_0_5_0_policy_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValidationError, "0.4.4"):
            WakhanPhasedCnaPolicy(
                profile_id="wakhan-obsolete-version-contract",
                mode="tumor_only",
                expected_version="0.5.0",
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_sha256=_sha_bytes(b"reference"),
                timeout_seconds=300,
                note="Obsolete Wakhan version contract must fail closed.",
            )

    def test_tumor_normal_command_uses_normal_phasing_and_cpd_without_breakpoints(self) -> None:
        argv = build_wakhan_argv(
            python_executable="python",
            wakhan_script=Path("/opt/wakhan/wakhan.py"),
            tumor_bam=Path("/tmp/tumor.bam"),
            phased_vcf=Path("/tmp/normal.phased.vcf.gz"),
            reference_fasta=Path("/tmp/ref.fa"),
            sample_id="TUMOR_001",
            output_dir=Path("/tmp/wakhan"),
            policy=_policy("a" * 64, mode="tumor_normal"),
            threads=8,
        )

        self.assertEqual(argv[:3], ("python", "/opt/wakhan/wakhan.py", "all"))
        self.assertIn("--normal-phased-vcf", argv)
        self.assertNotIn("--tumor-phased-vcf", argv)
        self.assertIn("--change-point-detection-for-cna", argv)
        self.assertNotIn("--breakpoints", argv)

    def test_tumor_only_command_uses_tumor_phasing_and_explicit_breakpoints(self) -> None:
        argv = build_wakhan_argv(
            python_executable="python",
            wakhan_script=Path("/opt/wakhan/wakhan.py"),
            tumor_bam=Path("/tmp/tumor.bam"),
            phased_vcf=Path("/tmp/tumor.phased.vcf.gz"),
            reference_fasta=Path("/tmp/ref.fa"),
            sample_id="TUMOR_001",
            output_dir=Path("/tmp/wakhan"),
            policy=_policy("a" * 64, mode="tumor_only"),
            threads=8,
            breakpoints_vcf=Path("/tmp/severus.vcf"),
        )

        self.assertIn("--tumor-phased-vcf", argv)
        self.assertNotIn("--normal-phased-vcf", argv)
        self.assertIn("--breakpoints", argv)
        self.assertNotIn("--change-point-detection-for-cna", argv)
        self.assertNotIn("--use-sv-haplotypes", argv)

    def test_missing_registered_phased_input_fails_before_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            phased = root / "normal.phased.vcf"
            reference = root / "ref.fa"
            script = root / "wakhan.py"
            _write_bam(tumor, b"tumor")
            _write_phased_vcf(phased)
            reference_sha256 = _write_reference(reference)
            script.write_text("# synthetic wakhan entry point\n", encoding="utf-8")
            runner = FakeWakhanRunner()

            with self.assertRaisesRegex(ValueError, "phased"):
                run_wakhan_phased_cna(
                    tumor_bam=tumor,
                    phased_vcf=phased,
                    reference_fasta=reference,
                    sample_id="TUMOR_001",
                    output_dir=root / "output",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        phased_vcf=None,
                        reference_sha256=reference_sha256,
                    ),
                    policy=_policy(reference_sha256, mode="tumor_normal"),
                    observed_runtime_version="0.4.4",
                    wakhan_script=script,
                    runner=runner,
                )

        self.assertEqual(runner.calls, [])

    def test_breakpoint_parent_hash_mismatch_fails_before_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            phased = root / "normal.phased.vcf"
            breakpoints = root / "severus.vcf"
            reference = root / "ref.fa"
            script = root / "wakhan.py"
            _write_bam(tumor, b"tumor")
            _write_phased_vcf(phased)
            _write_breakpoints(breakpoints)
            reference_sha256 = _write_reference(reference)
            script.write_text("# synthetic wakhan entry point\n", encoding="utf-8")
            runner = FakeWakhanRunner()

            with self.assertRaisesRegex(ValueError, "parent lane"):
                run_wakhan_phased_cna(
                    tumor_bam=tumor,
                    phased_vcf=phased,
                    breakpoints_vcf=breakpoints,
                    expected_breakpoint_parent_lane_id="lane-severus",
                    expected_breakpoint_parent_lane_sha256="b" * 64,
                    reference_fasta=reference,
                    sample_id="TUMOR_001",
                    output_dir=root / "output",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        phased_vcf=phased,
                        reference_sha256=reference_sha256,
                        breakpoints_vcf=breakpoints,
                        breakpoint_parent_id="lane-severus",
                        breakpoint_parent_sha256="c" * 64,
                    ),
                    policy=_policy(reference_sha256, mode="tumor_normal"),
                    observed_runtime_version="0.4.4",
                    wakhan_script=script,
                    runner=runner,
                )

        self.assertEqual(runner.calls, [])

    def test_successful_breakpoint_run_retains_native_outputs_and_parent_lock(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            phased = root / "normal.phased.vcf"
            breakpoints = root / "severus.vcf"
            reference = root / "ref.fa"
            script = root / "wakhan.py"
            _write_bam(tumor, b"tumor")
            _write_phased_vcf(phased)
            _write_breakpoints(breakpoints)
            reference_sha256 = _write_reference(reference)
            script.write_text("# synthetic wakhan entry point\n", encoding="utf-8")
            parent_sha256 = "d" * 64
            runner = FakeWakhanRunner()

            report = run_wakhan_phased_cna(
                tumor_bam=tumor,
                phased_vcf=phased,
                breakpoints_vcf=breakpoints,
                expected_breakpoint_parent_lane_id="lane-severus",
                expected_breakpoint_parent_lane_sha256=parent_sha256,
                reference_fasta=reference,
                sample_id="TUMOR_001",
                output_dir=root / "output",
                inputs=_bundle(
                    tumor_bam=tumor,
                    phased_vcf=phased,
                    reference_sha256=reference_sha256,
                    breakpoints_vcf=breakpoints,
                    breakpoint_parent_id="lane-severus",
                    breakpoint_parent_sha256=parent_sha256,
                ),
                policy=_policy(reference_sha256, mode="tumor_normal"),
                observed_runtime_version="0.4.4",
                wakhan_script=script,
                runner=runner,
            )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.mode, "tumor_normal")
        self.assertEqual(report.breakpoint_parent_lane_id, "lane-severus")
        self.assertEqual(report.breakpoint_parent_lane_sha256, parent_sha256)
        self.assertEqual(report.tool.version, "0.4.4")
        native_paths = {item.relative_path for item in report.native_artifacts}
        self.assertIn("solutions_ranks.tsv", native_paths)
        self.assertIn(
            "solution_2.0_0.70_0.90/integer_profile.bed",
            native_paths,
        )
        self.assertIn("coverage_data/coverage.csv", native_paths)
        self.assertIn("phasing_output/rephased.vcf.gz", native_paths)
        self.assertFalse(report.real_tool_qualified)
        self.assertTrue(report.research_only)

    def test_runtime_version_drift_fails_before_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            phased = root / "normal.phased.vcf"
            reference = root / "ref.fa"
            script = root / "wakhan.py"
            _write_bam(tumor, b"tumor")
            _write_phased_vcf(phased)
            reference_sha256 = _write_reference(reference)
            script.write_text("# synthetic wakhan entry point\n", encoding="utf-8")
            runner = FakeWakhanRunner()

            with self.assertRaisesRegex(ValueError, "version"):
                run_wakhan_phased_cna(
                    tumor_bam=tumor,
                    phased_vcf=phased,
                    reference_fasta=reference,
                    sample_id="TUMOR_001",
                    output_dir=root / "output",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        phased_vcf=phased,
                        reference_sha256=reference_sha256,
                    ),
                    policy=_policy(reference_sha256, mode="tumor_normal"),
                    observed_runtime_version="0.4.3",
                    wakhan_script=script,
                    runner=runner,
                )

        self.assertEqual(runner.calls, [])


if __name__ == "__main__":
    unittest.main()
