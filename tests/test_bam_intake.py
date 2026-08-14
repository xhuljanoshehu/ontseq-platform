from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.bam_intake import AlignedBamInspector, parse_sam_header
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    ReferenceContig,
    ReferenceLock,
    SampleManifest,
    Verdict,
)

HEADER = "\n".join(
    [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:1000",
        "@SQ\tSN:chr2\tLN:2000",
        "@RG\tID:rg1\tSM:SECRET_INTERNAL_NAME",
        "@PG\tID:SECRET_OPERATOR_VALUE\tPN:minimap2",
    ]
)


class FakeSamtoolsRunner:
    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(argv)
        if normalized[1:] == ("--version",):
            return CommandResult(normalized, 0, "samtools 1.24\n", "")
        if normalized[1:3] == ("quickcheck", "-v"):
            return CommandResult(normalized, 0, "", "")
        if normalized[1:3] == ("view", "-H"):
            return CommandResult(normalized, 0, HEADER, "")
        if normalized[1] == "idxstats":
            return CommandResult(
                normalized,
                0,
                "chr1\t1000\t10\t0\nchr2\t2000\t20\t0\n*\t0\t0\t1\n",
                "",
            )
        raise AssertionError(f"Unexpected command: {normalized}")


def _lock(chr2_length: int = 2000) -> ReferenceLock:
    return ReferenceLock(
        reference_id="SYNTHETIC_REF",
        genome_build=GenomeBuild.GRCH38,
        contigs=[
            ReferenceContig(name="chr1", length=1000),
            ReferenceContig(name="chr2", length=chr2_length),
        ],
        source_fai_sha256="0" * 64,
    )


def _manifest(bam: Path, index: Path) -> SampleManifest:
    return SampleManifest(
        sample_id="SYNTHETIC_001",
        run_id="SYNTHETIC_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path=str(bam),
            index_path=str(index),
        ),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="SYNTHETIC_REF",
        ),
        analysis=AnalysisSpec(profile="lcwgs", modules=[]),
    )


class BamIntakeTests(unittest.TestCase):
    def test_parse_header_extracts_only_non_identifying_summary_fields(self) -> None:
        header = parse_sam_header(HEADER)
        self.assertEqual(header.sort_order, "coordinate")
        self.assertEqual(header.contigs, (("chr1", 1000), ("chr2", 2000)))
        self.assertEqual(header.read_group_count, 1)
        self.assertEqual(header.sample_tag_count, 1)
        self.assertEqual(header.program_count, 1)

    def test_aligned_bam_intake_passes_locked_synthetic_header(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "synthetic.bam"
            index = root / "synthetic.bam.bai"
            bam.write_bytes(b"synthetic-not-a-real-bam")
            index.write_bytes(b"synthetic-not-a-real-index")

            report = AlignedBamInspector(runner=FakeSamtoolsRunner()).inspect(
                _manifest(bam, index), _lock(), include_checksums=True
            )

        self.assertEqual(report.verdict, Verdict.PASS)
        self.assertIsNotNone(report.header)
        self.assertIsNotNone(report.input_fingerprint)
        self.assertIsNotNone(report.input_fingerprint.sha256)
        self.assertNotIn("SECRET_INTERNAL_NAME", report.model_dump_json())
        self.assertNotIn("SECRET_OPERATOR_VALUE", report.model_dump_json())

    def test_reference_length_mismatch_fails_intake(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "synthetic.bam"
            index = root / "synthetic.bam.bai"
            bam.write_bytes(b"synthetic")
            index.write_bytes(b"synthetic")
            report = AlignedBamInspector(runner=FakeSamtoolsRunner()).inspect(
                _manifest(bam, index), _lock(chr2_length=1999)
            )

        self.assertEqual(report.verdict, Verdict.FAIL)
        self.assertTrue(
            any(
                check.name == "sequence_dictionary" and check.status == "FAIL"
                for check in report.checks
            )
        )


if __name__ == "__main__":
    unittest.main()
