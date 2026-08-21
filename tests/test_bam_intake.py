from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.bam_intake import (
    AlignedBamInspector,
    ParsedBamHeader,
    _reference_check,
    parse_sam_header,
)
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CheckStatus,
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
    def __init__(self) -> None:
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(argv)
        self.calls.append(normalized)
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

            runner = FakeSamtoolsRunner()
            report = AlignedBamInspector(runner=runner).inspect(
                _manifest(bam, index), _lock(), include_checksums=True
            )

        self.assertEqual(report.verdict, Verdict.PASS)
        self.assertIsNotNone(report.header)
        self.assertIsNotNone(report.input_fingerprint)
        self.assertIsNotNone(report.input_fingerprint.sha256)
        self.assertNotIn("SECRET_INTERNAL_NAME", report.model_dump_json())
        self.assertNotIn("SECRET_OPERATOR_VALUE", report.model_dump_json())
        idxstats_call = next(call for call in runner.calls if call[1] == "idxstats")
        self.assertEqual(idxstats_call[1:3], ("idxstats", "-X"))
        self.assertEqual(idxstats_call[-2:], (str(bam), str(index)))

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

    def test_exact_sequence_dictionary_passes(self) -> None:
        status, message, details = _reference_check(parse_sam_header(HEADER), _lock())

        self.assertEqual(status, CheckStatus.PASS)
        self.assertIn("matches", message)
        self.assertEqual(details["extra_contigs"], 0)

    def test_missing_locked_contig_fails_with_actionable_counts(self) -> None:
        header = ParsedBamHeader("coordinate", (("chr1", 1000),), 0, 0, 0)

        status, message, details = _reference_check(header, _lock())

        self.assertEqual(status, CheckStatus.FAIL)
        self.assertEqual(details["missing_contigs"], 1)
        self.assertEqual(details["expected_reference_bases"], 3000)
        self.assertEqual(details["observed_reference_bases"], 1000)
        self.assertEqual(
            message,
            "Reference dictionary mismatch: expected 2, observed 1; "
            "1 missing, 0 length mismatches, 0 extra contigs",
        )

    def test_extra_bam_contig_fails_under_strict_policy(self) -> None:
        header = ParsedBamHeader(
            "coordinate", (("chr1", 1000), ("chr2", 2000), ("chr3", 3000)), 0, 0, 0
        )

        status, _message, details = _reference_check(header, _lock())

        self.assertEqual(status, CheckStatus.FAIL)
        self.assertEqual(details["extra_contigs"], 1)
        self.assertEqual(details["observed_reference_bases"], 6000)

    def test_extra_bam_contig_warns_only_under_explicit_policy(self) -> None:
        header = ParsedBamHeader(
            "coordinate", (("chr1", 1000), ("chr2", 2000), ("chr3", 3000)), 0, 0, 0
        )
        permissive_lock = _lock().model_copy(update={"allow_extra_contigs": True})

        status, message, _details = _reference_check(header, permissive_lock)

        self.assertEqual(status, CheckStatus.WARN)
        self.assertIn("explicit extra-contig policy", message)

    def test_reordered_sequence_dictionary_fails(self) -> None:
        header = ParsedBamHeader("coordinate", (("chr2", 2000), ("chr1", 1000)), 0, 0, 0)

        status, message, details = _reference_check(header, _lock())

        self.assertEqual(status, CheckStatus.FAIL)
        self.assertTrue(details["contig_order_mismatch"])
        self.assertIn("contig order differs", message)

    def test_nanorepeat_chr1_fixture_still_fails_against_its_subset_fai(self) -> None:
        """The official fixture kept 86 hs37d5 @SQ lines but bundles only chr1 FASTA."""
        extras = tuple((f"decoy_{number}", 1000 + number) for number in range(85))
        header = ParsedBamHeader("coordinate", (("1", 249_250_621), *extras), 1, 1, 2)
        chr1_only = ReferenceLock(
            reference_id="GRCh37_CHR1_FIXTURE",
            genome_build=GenomeBuild.GRCH37,
            contigs=[ReferenceContig(name="1", length=249_250_621)],
            source_fai_sha256="1" * 64,
        )

        status, message, details = _reference_check(header, chr1_only)

        self.assertEqual(status, CheckStatus.FAIL)
        self.assertEqual(details["expected_contigs"], 1)
        self.assertEqual(details["observed_contigs"], 86)
        self.assertEqual(details["extra_contigs"], 85)
        self.assertIn("85 extra contigs", message)

    def test_region_extraction_passes_when_full_alignment_dictionary_is_locked(self) -> None:
        """Read sparsity does not matter; compatibility is defined by the full BAM header."""
        full_lock = ReferenceLock(
            reference_id="FULL_REFERENCE",
            genome_build=GenomeBuild.GRCH37,
            contigs=[
                ReferenceContig(name="1", length=249_250_621),
                ReferenceContig(name="2", length=243_199_373),
            ],
            source_fai_sha256="2" * 64,
        )
        header = ParsedBamHeader("coordinate", (("1", 249_250_621), ("2", 243_199_373)), 1, 1, 2)

        status, _message, _details = _reference_check(header, full_lock)

        self.assertEqual(status, CheckStatus.PASS)

    def test_mismatch_diagnostics_do_not_copy_custom_sequence_names(self) -> None:
        sensitive_name = "PATIENT_123_PRIVATE"
        header = ParsedBamHeader("coordinate", (("chr1", 1000), (sensitive_name, 20)), 0, 0, 0)

        _status, message, details = _reference_check(header, _lock())

        persisted = f"{message} {details}"
        self.assertNotIn(sensitive_name, persisted)
        self.assertIn("observed_dictionary_sha256", details)


if __name__ == "__main__":
    unittest.main()
