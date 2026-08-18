from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cutesv import (
    CuteSVPolicy,
    cutesv_observations,
    normalize_cutesv_vcf,
)
from ontseq_platform.models import EventType, GenomeBuild, ModuleRunStatus, ToolRecord

VCF_HEADER = """##fileformat=VCFv4.2
##source=cuteSV_2.1.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
"""


def _record(*fields: str) -> str:
    return "\t".join(fields) + "\n"


def _policy() -> CuteSVPolicy:
    return CuteSVPolicy(min_support=5, min_sv_length=50)


class CuteSVNormalizerTests(unittest.TestCase):
    def test_vcf_normalization_and_privacy(self) -> None:
        vcf = VCF_HEADER
        vcf += _record(
            "chr1",
            "1000",
            "SECRET_ID",
            "N",
            "<DEL>",
            "60",
            "PASS",
            "SVTYPE=DEL;END=1200;SVLEN=-201;RE=8;AF=0.4;RNAMES=SECRET_READ",
            "GT:DR:DV",
            "0/1:12:8",
        )
        vcf += _record(
            "chr1",
            "2000",
            "LOW_SUPPORT",
            "N",
            "<DUP>",
            ".",
            "PASS",
            "SVTYPE=DUP;END=2400;SVLEN=401;RE=2;AF=0.2",
            "GT:DR:DV",
            "0/1:8:2",
        )
        vcf += _record(
            "chr1",
            "3000",
            "BND_SECRET",
            "N",
            "N]chr2:5000]",
            ".",
            "PASS",
            "SVTYPE=BND;RE=7;AF=0.35;RNAMES=PRIVATE_READ",
            "GT:DR:DV",
            "0/1:13:7",
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.cutesv.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_cutesv_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="cuteSV", version="2.1.4"),
            )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.raw_record_count, 3)
        self.assertEqual(report.accepted_record_count, 2)
        self.assertEqual(report.rejected_record_count, 1)
        self.assertEqual(report.rejection_counts, {"support_below_policy": 1})

        deletion, breakend = report.events
        self.assertEqual(deletion.event_id, "CUTESV-000001")
        self.assertEqual(deletion.event_type, EventType.DELETION)
        self.assertEqual((deletion.primary.start, deletion.primary.end), (999, 1200))
        self.assertEqual(deletion.evidence[0].support_reads, 8)
        self.assertEqual(deletion.evidence[0].variant_allele_fraction, 0.4)

        self.assertEqual(breakend.event_id, "CUTESV-000003")
        self.assertEqual(breakend.event_type, EventType.TRANSLOCATION)
        self.assertEqual((breakend.primary.start, breakend.primary.end), (2999, 3000))
        self.assertIsNotNone(breakend.secondary)
        assert breakend.secondary is not None
        self.assertEqual(
            (breakend.secondary.chromosome, breakend.secondary.start, breakend.secondary.end),
            ("chr2", 4999, 5000),
        )

        serialized = report.model_dump_json()
        for forbidden in [
            "SECRET_ID",
            "SECRET_READ",
            "BND_SECRET",
            "PRIVATE_READ",
            "RNAMES",
        ]:
            self.assertNotIn(forbidden, serialized)

    def test_no_call_is_not_negative(self) -> None:
        vcf = VCF_HEADER + _record(
            "chr1",
            "1000",
            "LOW_SUPPORT",
            "N",
            "<DEL>",
            ".",
            "PASS",
            "SVTYPE=DEL;END=1200;SVLEN=-201;RE=1;AF=0.1",
            "GT:DR:DV",
            "0/1:9:1",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.cutesv.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_cutesv_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="cuteSV", version="2.1.4"),
            )

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.events, [])
        self.assertTrue(
            any("not a biological or clinical negative" in item for item in report.warnings)
        )

    def test_malformed_bnd_without_alt_mate_is_rejected(self) -> None:
        vcf = VCF_HEADER + _record(
            "chr1",
            "3000",
            "BAD_BND",
            "N",
            "<BND>",
            ".",
            "PASS",
            "SVTYPE=BND;RE=7;AF=0.35",
            "GT:DR:DV",
            "0/1:13:7",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.cutesv.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_cutesv_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="cuteSV", version="2.1.4"),
            )

        self.assertEqual(report.rejection_counts, {"missing_breakend_mate": 1})
        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)

    def test_version_lock_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.cutesv.vcf"
            path.write_text(VCF_HEADER, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "does not match policy lock"):
                normalize_cutesv_vcf(
                    path,
                    sample_id="SYNTHETIC_001",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                    tool=ToolRecord(name="cuteSV", version="2.1.3"),
                )

    def test_normalized_events_bridge_to_concordance_observations(self) -> None:
        vcf = VCF_HEADER + _record(
            "chr1",
            "3000",
            "BND_SECRET",
            "N",
            "N]chr2:5000]",
            ".",
            "PASS",
            "SVTYPE=BND;RE=7;AF=0.35",
            "GT:DR:DV",
            "0/1:13:7",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.cutesv.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_cutesv_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="cuteSV", version="2.1.4"),
            )

        observations = cutesv_observations(report)
        self.assertEqual(len(observations), 1)
        self.assertEqual(observations[0].caller, "cuteSV")
        self.assertEqual(observations[0].event_type, EventType.TRANSLOCATION)
        self.assertTrue(observations[0].research_only)
        self.assertFalse(observations[0].reportable)


if __name__ == "__main__":
    unittest.main()
