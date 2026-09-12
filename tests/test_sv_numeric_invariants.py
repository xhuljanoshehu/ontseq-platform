from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cutesv import normalize_cutesv_vcf
from ontseq_platform.models import CuteSvPolicy, GenomeBuild, SnifflesPolicy, ToolRecord
from ontseq_platform.sniffles import normalize_sniffles_vcf


def _write(path: Path, body: str) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + body,
        encoding="utf-8",
    )


class SvNumericInvariantTests(unittest.TestCase):
    def test_sniffles_rejects_nonfinite_quality_and_coverage(self) -> None:
        policy = SnifflesPolicy(
            profile_id="test",
            status="technical_defaults_only",
            min_support=5,
            min_sv_length=50,
            mapq=20,
            note="test",
        )
        tool = ToolRecord(name="Sniffles2", version=policy.expected_version)
        rows = {
            "quality": (
                "chr1\t1000\t.\tN\t<DEL>\tinf\tPASS\t"
                "SVTYPE=DEL;END=1200;SVLEN=-201;SUPPORT=8\tGT:DR:DV\t0/1:8:8\n",
                "malformed_quality",
            ),
            "coverage": (
                "chr1\t1000\t.\tN\t<DEL>\t60\tPASS\t"
                "SVTYPE=DEL;END=1200;SVLEN=-201;SUPPORT=8;COVERAGE=10,inf\t"
                "GT:DR:DV\t0/1:8:8\n",
                "malformed_coverage",
            ),
            "vaf": (
                "chr1\t1000\t.\tN\t<DEL>\t60\tPASS\t"
                "SVTYPE=DEL;END=1200;SVLEN=-201;SUPPORT=8;VAF=nan\t"
                "GT:DR:DV\t0/1:8:8\n",
                "malformed_vaf",
            ),
        }
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "sniffles.vcf"
            for name, (row, reason) in rows.items():
                with self.subTest(name=name):
                    _write(path, row)
                    report = normalize_sniffles_vcf(
                        path,
                        sample_id="S",
                        genome_build=GenomeBuild.GRCH38,
                        policy=policy,
                        tool=tool,
                    )
                    self.assertEqual(report.accepted_record_count, 0)
                    self.assertEqual(report.rejection_counts, {reason: 1})

    def test_cutesv_rejects_nonfinite_quality(self) -> None:
        policy = CuteSvPolicy(
            profile_id="test",
            status="technical_defaults_only",
            note="test",
        )
        tool = ToolRecord(name="cuteSV", version=policy.expected_version)
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "cutesv.vcf"
            _write(
                path,
                "chr1\t1000\t.\tN\t<DEL>\tinf\tPASS\t"
                "SVTYPE=DEL;END=1200;SVLEN=-201;RE=8\tGT:DR:DV\t0/1:8:8\n",
            )
            report = normalize_cutesv_vcf(
                path,
                sample_id="S",
                genome_build=GenomeBuild.GRCH38,
                policy=policy,
                tool=tool,
            )
            self.assertEqual(report.accepted_record_count, 0)
            self.assertEqual(report.rejection_counts, {"malformed_quality": 1})


if __name__ == "__main__":
    unittest.main()
