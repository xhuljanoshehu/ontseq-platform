"""Pinned Wakhan 0.4.4 interoperability using generated synthetic long-read inputs only."""

from __future__ import annotations

import hashlib
import importlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.models import GenomeBuild, ModuleRunStatus
from ontseq_platform.multicaller_contracts import CallerInputRole
from ontseq_platform.tumor.wakhan import WakhanPhasedCnaPolicy, run_wakhan_phased_cna
from ontseq_platform.tumor_inputs import (
    TumorAuxiliaryInputArtifact,
    TumorInputBundle,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@unittest.skipUnless(
    os.environ.get("ONTSEQ_WAKHAN_REAL_TOOL") == "1",
    "opt-in pinned Wakhan 0.4.4 interoperability lane",
)
class WakhanBinaryTests(unittest.TestCase):
    def setUp(self) -> None:
        wakhan = os.environ.get("ONTSEQ_WAKHAN_EXECUTABLE") or shutil.which("wakhan")
        self.assertIsNotNone(wakhan)
        assert wakhan is not None
        self.wakhan_script = Path(wakhan)
        self.assertTrue(self.wakhan_script.is_file())

        runtime_python = os.environ.get("ONTSEQ_WAKHAN_PYTHON") or sys.executable
        self.runtime_python = str(Path(runtime_python))
        self.assertTrue(Path(self.runtime_python).is_file())
        completed = subprocess.run(
            [
                self.runtime_python,
                "-c",
                "import importlib.metadata; print(importlib.metadata.version('wakhan'))",
            ],
            check=True,
            capture_output=True,
            encoding="utf-8",
        )
        self.runtime_version = completed.stdout.strip()
        self.assertEqual(self.runtime_version, "0.4.4")

        self.pysam = importlib.import_module("pysam")
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

        self.reference = self.root / "synthetic.reference.fa"
        self.reference_length = 500_000
        self.reference_contigs = ("chr7", "chr8")
        self._write_reference()
        self.reference_sha256 = _sha256(self.reference)

        self.variant_positions = tuple(range(20_000, 480_001, 20_000))
        self.phased_vcf = self.root / "synthetic.phased.vcf.gz"
        self._write_phased_vcf()

        self.tumor_bam = self.root / "tumor.bam"
        self._write_tumor_bam()

    def _write_reference(self) -> None:
        line = "A" * 80
        with self.reference.open("w", encoding="ascii") as handle:
            for chromosome in self.reference_contigs:
                handle.write(f">{chromosome}\n")
                full_lines, remainder = divmod(self.reference_length, len(line))
                for _ in range(full_lines):
                    handle.write(f"{line}\n")
                if remainder:
                    handle.write(f"{'A' * remainder}\n")
        self.pysam.faidx(str(self.reference))

    def _write_phased_vcf(self) -> None:
        plain = self.root / "synthetic.phased.vcf"
        rows = [
            f"chr7\t{position + 1}\trs{index}\tA\tG\t60\tPASS\t.\tGT:PS\t0|1:1"
            for index, position in enumerate(self.variant_positions, start=1)
        ]
        plain.write_text(
            "##fileformat=VCFv4.2\n"
            + "".join(
                f"##contig=<ID={chromosome},length={self.reference_length}>\n"
                for chromosome in self.reference_contigs
            )
            + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_TUMOR\n"
            + "\n".join(rows)
            + "\n",
            encoding="utf-8",
        )
        self.pysam.tabix_compress(str(plain), str(self.phased_vcf), force=True)
        self.pysam.tabix_index(str(self.phased_vcf), preset="vcf", force=True)

    def _write_tumor_bam(self) -> None:
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": chromosome, "LN": self.reference_length}
                for chromosome in self.reference_contigs
            ],
            "RG": [{"ID": "rg1", "SM": "SYNTHETIC_TUMOR"}],
        }
        read_length = 20_000
        with self.pysam.AlignmentFile(str(self.tumor_bam), "wb", header=header) as handle:
            read_number = 0
            for start in range(0, self.reference_length - read_length, 2_000):
                if 150_000 <= start < 260_000:
                    copies = 3
                elif 320_000 <= start < 420_000:
                    copies = 9
                else:
                    copies = 6
                overlapping = [
                    position
                    for position in self.variant_positions
                    if start <= position < start + read_length
                ]
                for copy_index in range(copies):
                    sequence = bytearray(b"A" * read_length)
                    if copy_index % 2:
                        for position in overlapping:
                            sequence[position - start] = ord("G")
                    read = self.pysam.AlignedSegment()
                    read.query_name = f"SYNTHETIC_TUMOR-read-{read_number:06d}"
                    read.query_sequence = sequence.decode("ascii")
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 60
                    read.cigar = ((0, read_length),)
                    read.query_qualities = self.pysam.qualitystring_to_array("I" * read_length)
                    read.set_tag("RG", "rg1")
                    read.set_tag("NM", len(overlapping) if copy_index % 2 else 0)
                    handle.write(read)
                    read_number += 1
        self.pysam.index(str(self.tumor_bam))

    def _bundle(self) -> TumorInputBundle:
        return TumorInputBundle(
            analysis_sample_id="SYNTHETIC_TUMOR",
            genome_build=GenomeBuild.GRCH38,
            reference_id="synthetic-grch38",
            reference_sha256=self.reference_sha256,
            artifacts=[
                TumorAuxiliaryInputArtifact(
                    artifact_id="tumor-bam",
                    role=CallerInputRole.TUMOR_BAM,
                    analysis_sample_id="SYNTHETIC_TUMOR",
                    source_sample_id="SYNTHETIC_TUMOR",
                    genome_build=GenomeBuild.GRCH38,
                    reference_id="synthetic-grch38",
                    reference_sha256=self.reference_sha256,
                    sha256=_sha256(self.tumor_bam),
                ),
                TumorAuxiliaryInputArtifact(
                    artifact_id="tumor-phased-vcf",
                    role=CallerInputRole.PHASED_VARIANTS,
                    analysis_sample_id="SYNTHETIC_TUMOR",
                    source_sample_id="SYNTHETIC_TUMOR",
                    genome_build=GenomeBuild.GRCH38,
                    reference_id="synthetic-grch38",
                    reference_sha256=self.reference_sha256,
                    sha256=_sha256(self.phased_vcf),
                ),
            ],
        )

    def _policy(self) -> WakhanPhasedCnaPolicy:
        return WakhanPhasedCnaPolicy(
            profile_id="wakhan-0.4.4-tumor-only-real-tool-smoke",
            mode="tumor_only",
            expected_version="0.4.4",
            real_tool_qualified=True,
            genome_build=GenomeBuild.GRCH38,
            reference_id="synthetic-grch38",
            reference_sha256=self.reference_sha256,
            timeout_seconds=900,
            note="Synthetic Wakhan 0.4.4 tumor-only interoperability smoke.",
        )

    def test_real_wakhan_0_4_4_runs_through_adapter(self) -> None:
        output_dir = self.root / "wakhan-output"
        report = run_wakhan_phased_cna(
            tumor_bam=self.tumor_bam,
            phased_vcf=self.phased_vcf,
            reference_fasta=self.reference,
            sample_id="SYNTHETIC_TUMOR",
            output_dir=output_dir,
            inputs=self._bundle(),
            policy=self._policy(),
            observed_runtime_version=self.runtime_version,
            wakhan_script=self.wakhan_script,
            python_executable=self.runtime_python,
            threads=2,
        )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.tool.version, "0.4.4")
        self.assertTrue(report.real_tool_qualified)
        self.assertTrue(report.policy.real_tool_qualified)
        self.assertTrue(report.native_artifacts)
        roles = {item.role for item in report.native_artifacts}
        self.assertIn("ranked_solutions", roles)
        self.assertIn("integer_copy_number_profile_bed", roles)
        self.assertTrue(all(item.sensitive_output for item in report.native_artifacts))
        self.assertTrue(all(not item.exportable for item in report.native_artifacts))


if __name__ == "__main__":
    unittest.main()
