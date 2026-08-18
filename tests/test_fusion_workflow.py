from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.breakends import SnifflesJunctionOrientation
from ontseq_platform.fusion import GeneAnnotationIndex
from ontseq_platform.fusion_workflow import interpret_sniffles_vcf_fusions
from ontseq_platform.models import (
    EventType,
    Evidence,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    SnifflesCallReport,
    SnifflesPolicy,
    ToolRecord,
)


class FusionWorkflowTests(unittest.TestCase):
    def test_vcf_descriptor_is_joined_without_raw_vcf_content_leak(self) -> None:
        event = GenomicEvent(
            event_id="SNIFFLES2-000001",
            event_type=EventType.TRANSLOCATION,
            primary=Locus(chromosome="chr9", start=133_729_450, end=133_729_451),
            secondary=Locus(chromosome="chr22", start=23_632_500, end=23_632_501),
            evidence=[
                Evidence(
                    caller="Sniffles2",
                    caller_version="2.8.0",
                    support_reads=12,
                    precise=True,
                )
            ],
            confidence="unclassified",
            reportable=False,
        )
        normalized = SnifflesCallReport(
            sample_id="SYNTHETIC_FUSION_001",
            genome_build=GenomeBuild.GRCH38,
            status=ModuleRunStatus.COMPLETED,
            policy=SnifflesPolicy(
                profile_id="synthetic-fusion",
                status="technical_defaults_only",
                note="Synthetic technical policy only.",
            ),
            events=[event],
            raw_record_count=1,
            accepted_record_count=1,
            rejected_record_count=0,
            rejection_counts={},
            tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            vcf_fingerprint=FileFingerprint(size_bytes=1, sha256="0" * 64),
        )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            annotation_path = root / "genes.bed"
            annotation_path.write_text(
                "chr9\t133729000\t133730000\tABL1\t0\t+\tENST_SYNTH_ABL1\n"
                "chr22\t23632000\t23633000\tBCR\t0\t+\tENST_SYNTH_BCR\n",
                encoding="utf-8",
            )
            annotation = GeneAnnotationIndex.from_bed(
                annotation_path,
                resource_id="synthetic-genes",
                resource_version="v1",
                genome_build=GenomeBuild.GRCH38,
            )
            vcf_path = root / "secret-source-name.vcf"
            raw_alt = "PRIVATESEQ]chr22:23632501]"
            vcf_path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                f"chr9\t133729451\tPRIVATE_VCF_ID\tN\t{raw_alt}\t60\tPASS\t"
                "SVTYPE=BND;SUPPORT=12\n",
                encoding="utf-8",
            )
            fusion_report = interpret_sniffles_vcf_fusions(
                normalized,
                vcf_path,
                annotation,
                known_pairs={("BCR", "ABL1")},
            )

        self.assertEqual(fusion_report.breakend_descriptor_count, 1)
        self.assertEqual(len(fusion_report.candidates), 1)
        candidate = fusion_report.candidates[0]
        self.assertIsNotNone(candidate.breakend_descriptor)
        assert candidate.breakend_descriptor is not None
        self.assertEqual(
            candidate.breakend_descriptor.sniffles_junction_orientation,
            SnifflesJunctionOrientation.PLUS_PLUS,
        )
        self.assertFalse(candidate.gene_pairs[0].orientation_resolved)
        serialized = fusion_report.model_dump_json()
        self.assertNotIn("PRIVATE_VCF_ID", serialized)
        self.assertNotIn(raw_alt, serialized)
        self.assertNotIn("secret-source-name.vcf", serialized)
        self.assertNotIn("PRIVATESEQ", serialized)

    def test_vcf_locus_mismatch_becomes_unresolved_not_silent_join(self) -> None:
        event = GenomicEvent(
            event_id="SNIFFLES2-000001",
            event_type=EventType.TRANSLOCATION,
            primary=Locus(chromosome="chr9", start=100, end=101),
            secondary=Locus(chromosome="chr22", start=200, end=201),
            evidence=[],
            confidence="unclassified",
            reportable=False,
        )
        normalized = SnifflesCallReport(
            sample_id="SYNTHETIC_FUSION_002",
            genome_build=GenomeBuild.GRCH38,
            status=ModuleRunStatus.COMPLETED,
            policy=SnifflesPolicy(
                profile_id="synthetic-fusion",
                status="technical_defaults_only",
                note="Synthetic technical policy only.",
            ),
            events=[event],
            raw_record_count=1,
            accepted_record_count=1,
            rejected_record_count=0,
            rejection_counts={},
            tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            vcf_fingerprint=FileFingerprint(size_bytes=1, sha256="0" * 64),
        )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            annotation_path = root / "genes.bed"
            annotation_path.write_text(
                "chr9\t1\t1000\tGENE9\t0\t+\nchr22\t1\t1000\tGENE22\t0\t+\n",
                encoding="utf-8",
            )
            annotation = GeneAnnotationIndex.from_bed(
                annotation_path,
                resource_id="synthetic-genes",
                resource_version="v1",
                genome_build=GenomeBuild.GRCH38,
            )
            vcf_path = root / "calls.vcf"
            vcf_path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr9\t102\tBND1\tN\tN]chr22:201]\t60\tPASS\tSVTYPE=BND;SUPPORT=8\n",
                encoding="utf-8",
            )
            fusion_report = interpret_sniffles_vcf_fusions(
                normalized,
                vcf_path,
                annotation,
            )

        self.assertEqual(fusion_report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(fusion_report.candidates, [])
        self.assertEqual(
            fusion_report.unresolved_source_event_ids,
            ["SNIFFLES2-000001"],
        )


if __name__ == "__main__":
    unittest.main()
