from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.sv.severus import (
    SeverusPairedPolicy,
    build_severus_paired_argv,
    run_severus_paired,
)

from ontseq_platform.execution import CommandResult
from ontseq_platform.models import EventType, GenomeBuild, ModuleRunStatus
from ontseq_platform.multicaller_contracts import CallerInputRole
from ontseq_platform.tumor_inputs import (
    TumorAuxiliaryInputArtifact,
    TumorInputBundle,
)


def _sha_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha_file(path: Path) -> str:
    return _sha_bytes(path.read_bytes())


def _write_bam_fixture(path: Path, payload: bytes) -> None:
    path.write_bytes(payload)
    Path(f"{path}.bai").write_bytes(b"synthetic-bai:" + payload)


def _bundle(
    tumor_bam: Path,
    normal_bam: Path | None,
    *,
    genome_build: GenomeBuild = GenomeBuild.GRCH38,
    reference_id: str = "GRCh38-test",
    reference_sha256: str | None = None,
) -> TumorInputBundle:
    ref_sha = reference_sha256 or _sha_bytes(b"reference")
    artifacts = [
        TumorAuxiliaryInputArtifact(
            artifact_id="tumor-bam",
            role=CallerInputRole.TUMOR_BAM,
            analysis_sample_id="TUMOR_001",
            source_sample_id="TUMOR_001",
            genome_build=genome_build,
            reference_id=reference_id,
            reference_sha256=ref_sha,
            sha256=_sha_file(tumor_bam),
        )
    ]
    if normal_bam is not None:
        artifacts.append(
            TumorAuxiliaryInputArtifact(
                artifact_id="normal-bam",
                role=CallerInputRole.MATCHED_NORMAL_BAM,
                analysis_sample_id="TUMOR_001",
                source_sample_id="NORMAL_001",
                genome_build=genome_build,
                reference_id=reference_id,
                reference_sha256=ref_sha,
                sha256=_sha_file(normal_bam),
            )
        )
    return TumorInputBundle(
        analysis_sample_id="TUMOR_001",
        genome_build=genome_build,
        reference_id=reference_id,
        reference_sha256=ref_sha,
        artifacts=artifacts,
    )


def _policy(
    *,
    genome_build: GenomeBuild = GenomeBuild.GRCH38,
    reference_id: str = "GRCh38-test",
    reference_sha256: str | None = None,
) -> SeverusPairedPolicy:
    return SeverusPairedPolicy(
        profile_id="severus-paired-test",
        expected_version="1.7",
        genome_build=genome_build,
        reference_id=reference_id,
        reference_sha256=reference_sha256 or _sha_bytes(b"reference"),
        min_support=3,
        min_sv_size=50,
        minimum_mapping_quality=10,
        tumor_in_normal_ratio=0.01,
        control_coverage_threshold=3,
        timeout_seconds=300,
        note="Synthetic paired Severus contract for research-only adapter tests.",
    )


class FakeSeverusRunner:
    def __init__(self, *, version: str = "1.7", empty: bool = False) -> None:
        self.version = version
        self.empty = empty
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        del timeout_seconds
        args = tuple(str(item) for item in argv)
        self.calls.append(args)
        if "--version" in args:
            return CommandResult(args, 0, f"{self.version}\n", "")

        out_dir = Path(args[args.index("--out-dir") + 1])
        somatic = out_dir / "somatic_SVs"
        plots = somatic / "plots"
        plots.mkdir(parents=True, exist_ok=True)
        records = "" if self.empty else (
            "chr7\t1001\tseverus_1\tN\t<DEL>\t60\tPASS\t"
            "SVTYPE=DEL;END=2000;DETAILED_TYPE=DEL;CLUSTERID=severus_1;HP=1\n"
        )
        (somatic / "severus_somatic.vcf").write_text(
            "##fileformat=VCFv4.2\n"
            "##source=Severus\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + records,
            encoding="utf-8",
        )
        (out_dir / "breakpoints_double.csv").write_text(
            "sample,breakpoint_1,breakpoint_2\nTUMOR_001,chr7:1000,chr7:2000\n",
            encoding="utf-8",
        )
        (somatic / "breakpoint_clusters_list.tsv").write_text(
            "#cluster_id\ttype\nseverus_1\tsimple\n",
            encoding="utf-8",
        )
        (somatic / "breakpoint_clusters.tsv").write_text(
            "#cluster_id\tadjacency\nseverus_1\tchr7:1000-chr7:2000\n",
            encoding="utf-8",
        )
        (plots / "severus_1.html").write_text("<html>synthetic graph</html>", encoding="utf-8")
        (out_dir / "read_qual.txt").write_text("synthetic quality\n", encoding="utf-8")
        (out_dir / "severus.log").write_text("Starting Severus 1.7\n", encoding="utf-8")
        return CommandResult(args, 0, "", "")


class SeverusPairedRuntimeTests(unittest.TestCase):
    def test_paired_command_requires_control_and_avoids_tumor_only_or_sensitive_flags(self) -> None:
        argv = build_severus_paired_argv(
            severus="severus",
            tumor_bam=Path("/tmp/tumor.bam"),
            normal_bam=Path("/tmp/normal.bam"),
            tumor_sample_id="TUMOR_001",
            normal_sample_id="NORMAL_001",
            output_dir=Path("/tmp/severus"),
            policy=_policy(),
            threads=4,
        )

        self.assertIn("--target-bam", argv)
        self.assertIn("--control-bam", argv)
        self.assertIn("--target-sample", argv)
        self.assertIn("--control-sample", argv)
        self.assertNotIn("--PON", argv)
        self.assertNotIn("--output-read-ids", argv)
        self.assertNotIn("--write-alignments", argv)

    def test_runtime_fails_closed_when_matched_normal_is_not_registered(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            _write_bam_fixture(tumor, b"tumor")
            inputs = _bundle(tumor, None)

            with self.assertRaisesRegex(ValueError, "matched normal"):
                run_severus_paired(
                    tumor_bam=tumor,
                    normal_bam=root / "missing-normal.bam",
                    tumor_sample_id="TUMOR_001",
                    normal_sample_id="NORMAL_001",
                    output_dir=root / "severus",
                    inputs=inputs,
                    policy=_policy(),
                    runner=FakeSeverusRunner(),
                )

    def test_runtime_rejects_bundle_build_or_reference_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            _write_bam_fixture(tumor, b"tumor")
            _write_bam_fixture(normal, b"normal")

            mismatches = (
                (
                    _bundle(tumor, normal, genome_build=GenomeBuild.GRCH37),
                    _policy(),
                    "genome build",
                ),
                (
                    _bundle(tumor, normal, reference_id="other-reference"),
                    _policy(),
                    "reference ID",
                ),
            )
            for inputs, policy, error in mismatches:
                with self.subTest(error=error), self.assertRaisesRegex(ValueError, error):
                    run_severus_paired(
                        tumor_bam=tumor,
                        normal_bam=normal,
                        tumor_sample_id="TUMOR_001",
                        normal_sample_id="NORMAL_001",
                        output_dir=root / f"severus-{error.replace(' ', '-')}",
                        inputs=inputs,
                        policy=policy,
                        runner=FakeSeverusRunner(),
                    )

    def test_runtime_retains_somatic_vcf_breakpoints_clusters_and_graphs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            _write_bam_fixture(tumor, b"tumor")
            _write_bam_fixture(normal, b"normal")
            runner = FakeSeverusRunner()

            report = run_severus_paired(
                tumor_bam=tumor,
                normal_bam=normal,
                tumor_sample_id="TUMOR_001",
                normal_sample_id="NORMAL_001",
                output_dir=root / "severus",
                inputs=_bundle(tumor, normal),
                policy=_policy(),
                runner=runner,
                threads=2,
            )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(len(report.events), 1)
        self.assertEqual(report.events[0].event_type, EventType.DELETION)
        self.assertFalse(report.events[0].reportable)
        paths = {artifact.relative_path for artifact in report.native_artifacts}
        self.assertIn("somatic_SVs/severus_somatic.vcf", paths)
        self.assertIn("breakpoints_double.csv", paths)
        self.assertIn("somatic_SVs/breakpoint_clusters_list.tsv", paths)
        self.assertIn("somatic_SVs/breakpoint_clusters.tsv", paths)
        self.assertIn("somatic_SVs/plots/severus_1.html", paths)
        run_argv = runner.calls[-1]
        self.assertNotIn("--PON", run_argv)
        self.assertNotIn("--output-read-ids", run_argv)

    def test_empty_somatic_vcf_is_no_call_not_biological_negative(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            _write_bam_fixture(tumor, b"tumor")
            _write_bam_fixture(normal, b"normal")

            report = run_severus_paired(
                tumor_bam=tumor,
                normal_bam=normal,
                tumor_sample_id="TUMOR_001",
                normal_sample_id="NORMAL_001",
                output_dir=root / "severus",
                inputs=_bundle(tumor, normal),
                policy=_policy(),
                runner=FakeSeverusRunner(empty=True),
            )

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.events, [])
        self.assertTrue(any("not a biological negative" in item for item in report.warnings))

    def test_runtime_rejects_version_drift_before_caller_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            _write_bam_fixture(tumor, b"tumor")
            _write_bam_fixture(normal, b"normal")
            runner = FakeSeverusRunner(version="1.6")

            with self.assertRaisesRegex(ValueError, "version"):
                run_severus_paired(
                    tumor_bam=tumor,
                    normal_bam=normal,
                    tumor_sample_id="TUMOR_001",
                    normal_sample_id="NORMAL_001",
                    output_dir=root / "severus",
                    inputs=_bundle(tumor, normal),
                    policy=_policy(),
                    runner=runner,
                )

        self.assertEqual(len(runner.calls), 1)
        self.assertIn("--version", runner.calls[0])


if __name__ == "__main__":
    unittest.main()
