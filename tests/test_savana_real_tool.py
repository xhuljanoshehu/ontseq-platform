"""Pinned SAVANA interoperability using only generated synthetic inputs."""

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
from ontseq_platform.tumor.savana import (
    SavanaPairedPolicy,
    SavanaTumorOnlyPolicy,
    run_savana_paired,
    run_savana_tumor_only,
)
from ontseq_platform.tumor_inputs import TumorAuxiliaryInputArtifact, TumorInputBundle


_CONTIGS = [*(f"chr{number}" for number in range(1, 23)), "chrX", "chrY"]
_CONTIG_LENGTHS = {chromosome: 10_000 for chromosome in _CONTIGS}
_CONTIG_LENGTHS["chr7"] = 700_000
_SNP_STARTS = tuple(range(50_000, 600_001, 50_000))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@unittest.skipUnless(
    os.environ.get("ONTSEQ_SAVANA_REAL_TOOL") == "1",
    "opt-in pinned SAVANA 1.3.8 interoperability lane",
)
class SavanaBinaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertIsNotNone(shutil.which("savana"))
        self.pysam = importlib.import_module("pysam")
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

        self.reference = self.root / "synthetic.reference.fa"
        with self.reference.open("w", encoding="ascii") as handle:
            sequence_line = "A" * 80
            for chromosome in _CONTIGS:
                handle.write(f">{chromosome}\n")
                length = _CONTIG_LENGTHS[chromosome]
                full_lines, remainder = divmod(length, len(sequence_line))
                for _ in range(full_lines):
                    handle.write(f"{sequence_line}\n")
                if remainder:
                    handle.write(f"{'A' * remainder}\n")
        self.pysam.faidx(str(self.reference))
        self.reference_sha256 = _sha256(self.reference)

        self.tumor_bam = self.root / "tumor.bam"
        self.normal_bam = self.root / "normal.bam"
        self._write_bam(self.tumor_bam, "SYNTHETIC_TUMOR")
        self._write_bam(self.normal_bam, "SYNTHETIC_NORMAL")

        self.snp_vcf = self.root / "synthetic.snps.vcf"
        vcf_rows = [
            f"chr7\t{start + 1}\trs{index}\tA\tG\t60\tPASS\t.\tGT\t0/1"
            for index, start in enumerate(_SNP_STARTS, start=1)
        ]
        self.snp_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            f"##contig=<ID=chr7,length={_CONTIG_LENGTHS['chr7']}>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_NORMAL\n"
            + "\n".join(vcf_rows)
            + "\n",
            encoding="utf-8",
        )

    def _write_bam(self, path: Path, sample_id: str) -> None:
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": chromosome, "LN": _CONTIG_LENGTHS[chromosome]}
                for chromosome in _CONTIGS
            ],
            "RG": [{"ID": "rg1", "SM": sample_id}],
        }
        with self.pysam.AlignmentFile(str(path), "wb", header=header) as handle:
            read_number = 0
            heterozygous_starts = set(_SNP_STARTS)
            is_tumor = sample_id == "SYNTHETIC_TUMOR"
            for start in range(0, _CONTIG_LENGTHS["chr7"] - 1_000, 5_000):
                copies = 24 if start in heterozygous_starts else 4
                alt_copies = 6 if is_tumor else 12
                for copy in range(copies):
                    read = self.pysam.AlignedSegment()
                    read.query_name = f"{sample_id}-read-{read_number:06d}"
                    first_base = (
                        "G"
                        if start in heterozygous_starts and copy < alt_copies
                        else "A"
                    )
                    read.query_sequence = first_base + ("A" * 499)
                    read.flag = 0
                    read.reference_id = _CONTIGS.index("chr7")
                    read.reference_start = start
                    read.mapping_quality = 60
                    read.cigar = ((0, 500),)
                    read.query_qualities = self.pysam.qualitystring_to_array("I" * 500)
                    read.set_tag("RG", "rg1")
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

    def _bundle(self, *, paired: bool) -> TumorInputBundle:
        artifacts = [
            self._artifact(
                CallerInputRole.TUMOR_BAM,
                self.tumor_bam,
                artifact_id="tumor-bam",
                source_sample_id="SYNTHETIC_TUMOR",
            ),
            self._artifact(
                CallerInputRole.SNP_VCF,
                self.snp_vcf,
                artifact_id="snp-vcf",
                source_sample_id="SYNTHETIC_NORMAL",
            ),
        ]
        if paired:
            artifacts.append(
                self._artifact(
                    CallerInputRole.MATCHED_NORMAL_BAM,
                    self.normal_bam,
                    artifact_id="normal-bam",
                    source_sample_id="SYNTHETIC_NORMAL",
                )
            )
        return TumorInputBundle(
            analysis_sample_id="SYNTHETIC_TUMOR",
            genome_build=GenomeBuild.GRCH38,
            reference_id="synthetic-grch38",
            reference_sha256=self.reference_sha256,
            artifacts=artifacts,
        )

    def _paired_policy(self) -> SavanaPairedPolicy:
        return SavanaPairedPolicy(
            profile_id="savana-1.3.8-paired-real-tool-smoke",
            expected_version="1.3.8",
            real_tool_qualified=False,
            genome_build=GenomeBuild.GRCH38,
            reference_id="synthetic-grch38",
            reference_sha256=self.reference_sha256,
            minimum_sv_length=30,
            minimum_mapping_quality=5,
            minimum_support=3,
            minimum_allele_fraction=0.01,
            cn_bin_size_kbp=10,
            timeout_seconds=600,
            note="Synthetic paired SAVANA 1.3.8 interoperability smoke.",
        )

    def _tumor_only_policy(self) -> SavanaTumorOnlyPolicy:
        return SavanaTumorOnlyPolicy(
            profile_id="savana-1.3.8-tumor-only-real-tool-smoke",
            expected_version="1.3.8",
            real_tool_qualified=False,
            genome_build=GenomeBuild.GRCH38,
            reference_id="synthetic-grch38",
            reference_sha256=self.reference_sha256,
            minimum_sv_length=30,
            minimum_mapping_quality=5,
            minimum_support=3,
            minimum_allele_fraction=0.01,
            cn_bin_size_kbp=10,
            timeout_seconds=600,
            note="Synthetic tumor-only SAVANA 1.3.8 interoperability smoke.",
        )

    def _assert_research_report(self, report, *, mode: str) -> None:  # noqa: ANN001
        self.assertEqual(report.mode, mode)
        self.assertEqual(report.tool.version, "1.3.8")
        self.assertIn(report.status, {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL})
        self.assertTrue(report.native_artifacts)
        self.assertTrue(all(event.reportable is False for event in report.events))
        self.assertTrue(
            all(
                not artifact.exportable
                for artifact in report.native_artifacts
                if artifact.sensitive_output
            )
        )

    def test_real_paired_and_tumor_only_modes_complete_through_adapter(self) -> None:
        paired = run_savana_paired(
            tumor_bam=self.tumor_bam,
            normal_bam=self.normal_bam,
            snp_vcf=self.snp_vcf,
            reference_fasta=self.reference,
            sample_id="SYNTHETIC_TUMOR",
            normal_sample_id="SYNTHETIC_NORMAL",
            output_dir=self.root / "paired-output",
            inputs=self._bundle(paired=True),
            policy=self._paired_policy(),
            threads=1,
        )
        self._assert_research_report(paired, mode="paired")
        self.assertIsNotNone(paired.normal_bam_fingerprint)

        tumor_only = run_savana_tumor_only(
            tumor_bam=self.tumor_bam,
            snp_vcf=self.snp_vcf,
            reference_fasta=self.reference,
            sample_id="SYNTHETIC_TUMOR",
            output_dir=self.root / "tumor-only-output",
            inputs=self._bundle(paired=False),
            policy=self._tumor_only_policy(),
            threads=1,
        )
        self._assert_research_report(tumor_only, mode="tumor_only")
        self.assertIsNone(tumor_only.normal_bam_fingerprint)


if __name__ == "__main__":
    unittest.main()
