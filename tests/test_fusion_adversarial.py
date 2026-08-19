from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.fusion import GeneAnnotationIndex
from ontseq_platform.fusion_redundancy import analyze_fusion_redundancy
from ontseq_platform.fusion_workflow import interpret_sniffles_vcf_fusions
from ontseq_platform.models import GenomeBuild, ModuleRunStatus, SnifflesPolicy, ToolRecord
from ontseq_platform.sniffles import normalize_sniffles_vcf


class FusionAdversarialTests(unittest.TestCase):
    def _policy(self) -> SnifflesPolicy:
        return SnifflesPolicy(
            profile_id="synthetic-adversarial",
            status="technical_defaults_only",
            note="Synthetic adversarial software test only.",
        )

    def _tool(self) -> ToolRecord:
        return ToolRecord(name="Sniffles2", version="2.8.0")

    def _annotation(self, root: Path) -> GeneAnnotationIndex:
        path = root / "synthetic.genes.bed"
        path.write_text(
            "chr1\t50\t150\tGENE1\t0\t+\tENST_SYNTH_GENE1\n"
            "chr2\t150\t250\tGENE2\t0\t-\tENST_SYNTH_GENE2\n",
            encoding="utf-8",
        )
        return GeneAnnotationIndex.from_bed(
            path,
            resource_id="synthetic-adversarial-genes",
            resource_version="v1",
            genome_build=GenomeBuild.GRCH38,
        )

    def test_filtered_and_low_support_records_do_not_shift_descriptor_join(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "mixed.synthetic.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t101\tFILTERED\tN\tN]chr2:201]\t60\tLowQual\t"
                "SVTYPE=BND;SUPPORT=12\n"
                "chr1\t101\tACCEPTED\tN\tN]chr2:201]\t60\tPASS\t"
                "SVTYPE=BND;SUPPORT=12\n"
                "chr1\t101\tLOW_SUPPORT\tN\tN]chr2:201]\t60\tPASS\t"
                "SVTYPE=BND;SUPPORT=2\n",
                encoding="utf-8",
            )
            normalized = normalize_sniffles_vcf(
                vcf,
                sample_id="SYNTHETIC_ADV_001",
                genome_build=GenomeBuild.GRCH38,
                policy=self._policy(),
                tool=self._tool(),
            )
            report = interpret_sniffles_vcf_fusions(
                normalized,
                vcf,
                self._annotation(root),
            )

        self.assertEqual(normalized.raw_record_count, 3)
        self.assertEqual(normalized.accepted_record_count, 1)
        self.assertEqual(normalized.rejected_record_count, 2)
        self.assertEqual(normalized.events[0].event_id, "SNIFFLES2-000002")
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(len(report.candidates), 1)
        self.assertEqual(report.candidates[0].source_event_id, "SNIFFLES2-000002")
        self.assertIsNotNone(report.candidates[0].breakend_descriptor)
        self.assertEqual(report.breakend_descriptor_count, 1)

    def test_malformed_alt_can_normalize_via_chr2_end_but_orientation_stays_missing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "fallback.synthetic.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t101\tFALLBACK\tN\t<BND>\t60\tPASS\t"
                "SVTYPE=BND;CHR2=chr2;END=201;SUPPORT=12\n",
                encoding="utf-8",
            )
            normalized = normalize_sniffles_vcf(
                vcf,
                sample_id="SYNTHETIC_ADV_002",
                genome_build=GenomeBuild.GRCH38,
                policy=self._policy(),
                tool=self._tool(),
            )
            report = interpret_sniffles_vcf_fusions(
                normalized,
                vcf,
                self._annotation(root),
            )

        self.assertEqual(normalized.accepted_record_count, 1)
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(len(report.candidates), 1)
        self.assertIsNone(report.candidates[0].breakend_descriptor)
        self.assertEqual(report.breakend_descriptor_count, 0)
        self.assertEqual(
            report.missing_breakend_descriptor_event_ids,
            ["SNIFFLES2-000001"],
        )
        self.assertTrue(
            any("orientation remains unavailable" in warning for warning in report.warnings)
        )

    def test_reciprocal_breakends_are_flagged_without_auto_deduplication(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "reciprocal.synthetic.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t101\tBND_A\tN\tN]chr2:201]\t60\tPASS\t"
                "SVTYPE=BND;SUPPORT=12\n"
                "chr2\t201\tBND_B\tN\tN]chr1:101]\t60\tPASS\t"
                "SVTYPE=BND;SUPPORT=12\n",
                encoding="utf-8",
            )
            normalized = normalize_sniffles_vcf(
                vcf,
                sample_id="SYNTHETIC_ADV_003",
                genome_build=GenomeBuild.GRCH38,
                policy=self._policy(),
                tool=self._tool(),
            )
            report = interpret_sniffles_vcf_fusions(
                normalized,
                vcf,
                self._annotation(root),
            )
            redundancy = analyze_fusion_redundancy(report)

        self.assertEqual(len(report.candidates), 2)
        self.assertEqual(len(redundancy.groups), 1)
        self.assertEqual(
            redundancy.groups[0].source_event_ids,
            ["SNIFFLES2-000001", "SNIFFLES2-000002"],
        )
        self.assertEqual(redundancy.potentially_redundant_candidate_count, 2)
        self.assertFalse(redundancy.auto_deduplicated)
        self.assertFalse(redundancy.groups[0].auto_deduplicated)
        self.assertTrue(redundancy.warnings)

    def test_exact_duplicate_breakends_are_flagged_without_collapsing_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "duplicate.synthetic.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t101\tBND_A\tN\tN]chr2:201]\t60\tPASS\t"
                "SVTYPE=BND;SUPPORT=12\n"
                "chr1\t101\tBND_B\tN\tN]chr2:201]\t60\tPASS\t"
                "SVTYPE=BND;SUPPORT=9\n",
                encoding="utf-8",
            )
            normalized = normalize_sniffles_vcf(
                vcf,
                sample_id="SYNTHETIC_ADV_004",
                genome_build=GenomeBuild.GRCH38,
                policy=self._policy(),
                tool=self._tool(),
            )
            report = interpret_sniffles_vcf_fusions(
                normalized,
                vcf,
                self._annotation(root),
            )
            redundancy = analyze_fusion_redundancy(report)

        self.assertEqual(len(report.candidates), 2)
        self.assertEqual(len(redundancy.groups), 1)
        self.assertEqual(redundancy.candidate_count, 2)
        self.assertEqual(redundancy.potentially_redundant_candidate_count, 2)
        self.assertFalse(redundancy.auto_deduplicated)


if __name__ == "__main__":
    unittest.main()
