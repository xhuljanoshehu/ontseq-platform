"""Pinned Severus paired interoperability using generated synthetic long-read BAMs only."""

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


def _sha_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


@unittest.skipUnless(
    os.environ.get("ONTSEQ_SEVERUS_REAL_TOOL") == "1",
    "opt-in pinned Severus 1.7 paired interoperability lane",
)
class SeverusBinaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertIsNotNone(shutil.which("severus"))
        self.pysam = importlib.import_module("pysam")
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

        self.reference_id = "synthetic-grch38"
        self.reference_sha256 = _sha_text(self.reference_id)
        self.tumor_bam = self.root / "tumor.bam"
        self.normal_bam = self.root / "normal.bam"
        self._write_long_read_bam(self.tumor_bam, "SYNTHETIC_TUMOR")
        self._write_long_read_bam(self.normal_bam, "SYNTHETIC_NORMAL")

    def _write_long_read_bam(self, path: Path, sample_id: str) -> None:
        contig_length = 300_000
        read_length = 15_000
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr7", "LN": contig_length}],
            "RG": [{"ID": "rg1", "SM": sample_id}],
        }
        with self.pysam.AlignmentFile(str(path), "wb", header=header) as handle:
            for read_number, start in enumerate(range(0, contig_length - read_length, 1_000)):
                read = self.pysam.AlignedSegment()
                read.query_name = f"{sample_id}-read-{read_number:06d}"
                read.query_sequence = "A" * read_length
                read.flag = 0
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigar = ((0, read_length),)
                read.query_qualities = self.pysam.qualitystring_to_array("I" * read_length)
                read.set_tag("RG", "rg1")
                read.set_tag("NM", 0)
                handle.write(read)
        self.pysam.index(str(path))

    def _artifact(
        self,
        *,
        role: CallerInputRole,
        path: Path,
        artifact_id: str,
        source_sample_id: str,
    ) -> TumorAuxiliaryInputArtifact:
        return TumorAuxiliaryInputArtifact(
            artifact_id=artifact_id,
            role=role,
            analysis_sample_id="SYNTHETIC_TUMOR",
            source_sample_id=source_sample_id,
            genome_build=GenomeBuild.GRCH38,
            reference_id=self.reference_id,
            reference_sha256=self.reference_sha256,
            sha256=_sha256(path),
        )

    def _inputs(self) -> TumorInputBundle:
        return TumorInputBundle(
            analysis_sample_id="SYNTHETIC_TUMOR",
            genome_build=GenomeBuild.GRCH38,
            reference_id=self.reference_id,
            reference_sha256=self.reference_sha256,
            artifacts=[
                self._artifact(
                    role=CallerInputRole.TUMOR_BAM,
                    path=self.tumor_bam,
                    artifact_id="tumor-bam",
                    source_sample_id="SYNTHETIC_TUMOR",
                ),
                self._artifact(
                    role=CallerInputRole.MATCHED_NORMAL_BAM,
                    path=self.normal_bam,
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
            reference_id=self.reference_id,
            reference_sha256=self.reference_sha256,
            min_support=3,
            min_sv_size=50,
            minimum_mapping_quality=10,
            tumor_in_normal_ratio=0.01,
            control_coverage_threshold=1,
            timeout_seconds=900,
            note="Synthetic paired Severus 1.7 interoperability smoke.",
        )

    def test_real_severus_1_7_paired_run_completes_through_adapter(self) -> None:
        report = run_severus_paired(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            tumor_sample_id="SYNTHETIC_TUMOR",
            normal_sample_id="SYNTHETIC_NORMAL",
            output_dir=self.root / "severus-output",
            inputs=self._inputs(),
            policy=self._policy(),
            threads=2,
        )

        self.assertEqual(report.tool.version, "1.7")
        self.assertTrue(report.policy.real_tool_qualified)
        self.assertIn(report.status, {ModuleRunStatus.NO_CALL, ModuleRunStatus.COMPLETED})
        self.assertTrue(report.native_artifacts)
        self.assertTrue(report.somatic_vcf_fingerprint.sha256)
        self.assertTrue(all(event.reportable is False for event in report.events))
        self.assertEqual(report.tumor_bam_fingerprint.sha256, _sha256(self.tumor_bam))
        self.assertEqual(report.normal_bam_fingerprint.sha256, _sha256(self.normal_bam))


if __name__ == "__main__":
    unittest.main()
