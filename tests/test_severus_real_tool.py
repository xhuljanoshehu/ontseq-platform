"""Pinned Severus 1.7 interoperability using generated synthetic tumor/normal BAMs only."""

from __future__ import annotations

import hashlib
import importlib
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.models import GenomeBuild, ModuleRunStatus
from ontseq_platform.multicaller_contracts import CallerInputRole
from ontseq_platform.sv.severus import SeverusPairedPolicy, run_severus_paired
from ontseq_platform.tumor_inputs import (
    TumorAuxiliaryInputArtifact,
    TumorInputBundle,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@unittest.skipUnless(
    os.environ.get("ONTSEQ_SEVERUS_REAL_TOOL") == "1",
    "opt-in pinned Severus 1.7 interoperability lane",
)
class SeverusBinaryTests(unittest.TestCase):
    def setUp(self) -> None:
        severus = shutil.which("severus")
        self.assertIsNotNone(severus)
        assert severus is not None
        self.severus = severus
        self.pysam = importlib.import_module("pysam")

        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

        self.tumor_bam = self.root / "tumor.bam"
        self.normal_bam = self.root / "normal.bam"
        self._write_bam(self.tumor_bam, "SYNTHETIC_TUMOR", tumor=True)
        self._write_bam(self.normal_bam, "SYNTHETIC_NORMAL", tumor=False)

        self.reference_sha256 = hashlib.sha256(b"synthetic-severus-grch38-contract").hexdigest()

    def _write_bam(self, path: Path, sample_id: str, *, tumor: bool) -> None:
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr7", "LN": 250_000}],
            "RG": [{"ID": "rg1", "SM": sample_id}],
        }
        with self.pysam.AlignmentFile(str(path), "wb", header=header) as handle:
            read_number = 0
            for start in range(20_000, 220_000, 2_000):
                copies = 8 if tumor and 90_000 <= start < 120_000 else 6
                for _copy in range(copies):
                    read = self.pysam.AlignedSegment()
                    read.query_name = f"{sample_id}-read-{read_number:06d}"
                    read.query_sequence = "A" * 1_000
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 60
                    read.cigar = ((0, 1_000),)
                    read.query_qualities = self.pysam.qualitystring_to_array("I" * 1_000)
                    read.set_tag("RG", "rg1")
                    read.set_tag("NM", 0)
                    handle.write(read)
                    read_number += 1
        self.pysam.index(str(path))

    def _artifact(
        self,
        role: CallerInputRole,
        path: Path,
        *,
        artifact_id: str,
        source_sample_id: str,
    ) -> TumorAuxiliaryInputArtifact:
        return TumorAuxiliaryInputArtifact(
            artifact_id=artifact_id,
            role=role,
            analysis_sample_id="SYNTHETIC_TUMOR",
            source_sample_id=source_sample_id,
            genome_build=GenomeBuild.GRCH38,
            reference_id="synthetic-grch38",
            reference_sha256=self.reference_sha256,
            sha256=_sha256(path),
        )

    def _bundle(self) -> TumorInputBundle:
        return TumorInputBundle(
            analysis_sample_id="SYNTHETIC_TUMOR",
            genome_build=GenomeBuild.GRCH38,
            reference_id="synthetic-grch38",
            reference_sha256=self.reference_sha256,
            artifacts=[
                self._artifact(
                    CallerInputRole.TUMOR_BAM,
                    self.tumor_bam,
                    artifact_id="tumor-bam",
                    source_sample_id="SYNTHETIC_TUMOR",
                ),
                self._artifact(
                    CallerInputRole.MATCHED_NORMAL_BAM,
                    self.normal_bam,
                    artifact_id="normal-bam",
                    source_sample_id="SYNTHETIC_NORMAL",
                ),
            ],
        )

    def _policy(self) -> SeverusPairedPolicy:
        return SeverusPairedPolicy(
            profile_id="severus-1.7-paired-real-tool-smoke",
            expected_version="1.7",
            real_tool_qualified=True,
            genome_build=GenomeBuild.GRCH38,
            reference_id="synthetic-grch38",
            reference_sha256=self.reference_sha256,
            min_support=3,
            min_sv_size=50,
            minimum_mapping_quality=10,
            tumor_in_normal_ratio=0.01,
            control_coverage_threshold=0,
            timeout_seconds=900,
            note="Synthetic paired Severus 1.7 interoperability smoke.",
        )

    def test_real_severus_1_7_runs_through_paired_adapter(self) -> None:
        output_dir = self.root / "severus-output"
        report = run_severus_paired(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            tumor_sample_id="SYNTHETIC_TUMOR",
            normal_sample_id="SYNTHETIC_NORMAL",
            output_dir=output_dir,
            inputs=self._bundle(),
            policy=self._policy(),
            severus=self.severus,
            threads=2,
        )

        self.assertIn(report.status, {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL})
        self.assertEqual(report.tool.version, "1.7")
        self.assertTrue(report.policy.real_tool_qualified)
        self.assertTrue(report.native_artifacts)
        self.assertTrue(report.somatic_vcf_fingerprint.sha256)
        self.assertTrue(all(event.reportable is False for event in report.events))
        self.assertFalse(report.tool.parameters["output_read_ids"])
        self.assertFalse(report.tool.parameters["write_alignments"])
        native_paths = {item.relative_path for item in report.native_artifacts}
        self.assertTrue(
            {
                "somatic_SVs/severus_somatic.vcf",
                "somatic_SVs/severus_somatic.vcf.gz",
            }.intersection(native_paths)
        )


if __name__ == "__main__":
    unittest.main()
