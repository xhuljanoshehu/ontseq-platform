from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    EventType,
    GenomeBuild,
    InputKind,
    InputSpec,
    ModuleRunStatus,
    SampleManifest,
    SnifflesPolicy,
    ToolRecord,
    Verdict,
)
from ontseq_platform.sniffles import normalize_sniffles_vcf, run_sniffles

VCF_HEADER = """##fileformat=VCFv4.2
##source=Sniffles2_2.8.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
"""


def _policy() -> SnifflesPolicy:
    return SnifflesPolicy(
        profile_id="synthetic-conservative",
        status="technical_defaults_only",
        min_support=5,
        min_sv_length=50,
        mapq=20,
        note="Synthetic technical thresholds only.",
    )


def _record(*fields: str) -> str:
    return "\t".join(fields) + "\n"


def _manifest(bam: Path, bai: Path) -> SampleManifest:
    return SampleManifest(
        sample_id="SYNTHETIC_001",
        run_id="SYNTHETIC_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path=str(bam),
            index_path=str(bai),
        ),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="SYNTHETIC_REF",
        ),
        analysis=AnalysisSpec(
            profile="synthetic",
            modules=[AnalysisModule.QC, AnalysisModule.SV, AnalysisModule.REPORT],
        ),
    )


def _intake() -> AlignedBamIntakeReport:
    return AlignedBamIntakeReport(
        sample_id="SYNTHETIC_001",
        reference_id="SYNTHETIC_REF",
        genome_build=GenomeBuild.GRCH38,
        checks=[],
        verdict=Verdict.PASS,
    )


class RecordingRunner:
    def __init__(self, vcf_text: str) -> None:
        self.vcf_text = vcf_text
        self.commands: list[tuple[str, ...]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(str(item) for item in argv)
        self.commands.append(normalized)
        if "--version" in normalized:
            return CommandResult(normalized, 0, "Sniffles2, Version 2.8.0\n", "")
        output = Path(normalized[normalized.index("--vcf") + 1])
        output.write_text(self.vcf_text, encoding="utf-8")
        return CommandResult(normalized, 0, "", "")


class SnifflesAdapterTests(unittest.TestCase):
    def test_vcf_is_normalized_and_rejections_are_explicit(self) -> None:
        vcf = VCF_HEADER
        vcf += _record(
            "chr1",
            "1000",
            "SECRET_RECORD",
            "N",
            "<DEL>",
            "60",
            "PASS",
            "PRECISE;SVTYPE=DEL;END=1200;SVLEN=-201;SUPPORT=8;"
            "COVERAGE=10,12,14;STRAND=+;NM=0.4;STDEV_POS=1.2;STDEV_LEN=0;"
            "RNAMES=SECRET_READ",
            "GT:GQ:DR:DV",
            "0/1:60:8:8",
        )
        vcf += _record(
            "chr1",
            "2000",
            "LOW_SUPPORT",
            "N",
            "<DUP>",
            "60",
            "PASS",
            "SVTYPE=DUP;END=2400;SUPPORT=2",
            "GT:GQ:DR:DV",
            "0/1:30:9:2",
        )
        vcf += _record(
            "chr1",
            "3000",
            "FILTERED",
            "N",
            "<INV>",
            "50",
            "GT",
            "SVTYPE=INV;END=3400;SUPPORT=9",
            "GT:GQ:DR:DV",
            "0/1:30:5:9",
        )
        vcf += _record(
            "chr1",
            "4000",
            "BND_SECRET",
            "N",
            "N]chr2:5000]",
            "45",
            "PASS",
            "PRECISE;SVTYPE=BND;SUPPORT=7;VAF=0.35;STRAND=-",
            "GT:GQ:DR:DV",
            "0/1:40:13:7",
        )
        vcf += _record(
            "chrM",
            "6000",
            "NONCANONICAL",
            "N",
            "<DEL>",
            "40",
            "PASS",
            "SVTYPE=DEL;END=6200;SUPPORT=8",
            "GT:GQ:DR:DV",
            "0/1:40:8:8",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.raw_record_count, 5)
        self.assertEqual(report.accepted_record_count, 2)
        self.assertEqual(report.rejected_record_count, 3)
        self.assertEqual(report.rejection_counts["support_below_policy"], 1)
        self.assertEqual(report.rejection_counts["filter_not_pass:GT"], 1)
        self.assertEqual(report.rejection_counts["invalid_primary_locus"], 1)

        deletion, breakend = report.events
        self.assertEqual(deletion.event_type, EventType.DELETION)
        self.assertEqual((deletion.primary.start, deletion.primary.end), (999, 1200))
        self.assertEqual(deletion.length_bp, 201)
        self.assertEqual(deletion.evidence[0].variant_allele_fraction, 0.5)
        self.assertEqual(deletion.evidence[0].coverage_context, [10.0, 12.0, 14.0])
        self.assertEqual(deletion.evidence[0].supporting_read_strands, "+")
        self.assertEqual(deletion.evidence[0].position_standard_deviation, 1.2)
        self.assertEqual(deletion.evidence[0].length_standard_deviation, 0.0)
        self.assertTrue(deletion.evidence[0].precise)
        self.assertEqual(breakend.event_type, EventType.TRANSLOCATION)
        self.assertIsNotNone(breakend.secondary)
        assert breakend.secondary is not None
        self.assertEqual((breakend.secondary.chromosome, breakend.secondary.start), ("chr2", 4999))

        serialized = report.model_dump_json()
        self.assertNotIn("SECRET_RECORD", serialized)
        self.assertNotIn("SECRET_READ", serialized)
        self.assertNotIn("N]chr2:5000]", serialized)
        self.assertTrue(all(not event.reportable for event in report.events))

    def test_all_four_breakend_forms_resolve_without_persisting_alt_or_form(self) -> None:
        alternates = (
            "AC[chr2:5000[",
            "GT]chr2:5001]",
            "[chr2:5002[CA",
            "]chr2:5003]TG",
        )
        vcf = VCF_HEADER + "".join(
            _record(
                "chr1",
                str(1000 + index),
                f"PRIVATE_{index}",
                "N",
                alternate,
                "45",
                "PASS",
                "PRECISE;SVTYPE=BND;SUPPORT=7",
                "GT:DR:DV",
                "0/1:13:7",
            )
            for index, alternate in enumerate(alternates)
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )

        self.assertEqual(report.accepted_record_count, 4)
        self.assertEqual(
            [event.secondary.start if event.secondary else None for event in report.events],
            [4999, 5000, 5001, 5002],
        )
        serialized = report.model_dump_json()
        for alternate in alternates:
            self.assertNotIn(alternate, serialized)
        self.assertNotIn("PRIVATE_", serialized)
        self.assertNotIn("breakend_alt_form", serialized)
        self.assertNotIn("local_then_", serialized)

    def test_symbolic_breakend_keeps_chr2_end_fallback_without_inventing_alt_form(self) -> None:
        vcf = VCF_HEADER + _record(
            "chr1",
            "1000",
            ".",
            "N",
            "<BND>",
            "45",
            "PASS",
            "PRECISE;SVTYPE=BND;CHR2=chr2;END=5000;SUPPORT=7",
            "GT:DR:DV",
            "0/1:13:7",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )

        event = report.events[0]
        self.assertEqual(event.secondary.chromosome if event.secondary else None, "chr2")
        self.assertEqual(event.secondary.start if event.secondary else None, 4999)
        self.assertNotIn("breakend_alt_form", report.model_dump_json())

    def test_bracket_breakend_info_mate_must_agree_with_alt(self) -> None:
        info_values = (
            "PRECISE;SVTYPE=BND;CHR2=2;END=05000;SUPPORT=7",
            "PRECISE;SVTYPE=BND;CHR2=chr3;END=5000;SUPPORT=7",
            "PRECISE;SVTYPE=BND;CHR2=chr2;END=5001;SUPPORT=7",
            "PRECISE;SVTYPE=BND;CHR2=chr2;SUPPORT=7",
            "PRECISE;SVTYPE=BND;CHR2=2;SUPPORT=7",
            "PRECISE;SVTYPE=BND;CHR2=chr3;SUPPORT=7",
            "PRECISE;SVTYPE=BND;CHR2=chr2,chr3;SUPPORT=7",
        )
        vcf = VCF_HEADER + "".join(
            _record(
                "chr1",
                str(1000 + index),
                ".",
                "N",
                "N[chr2:5000[",
                "45",
                "PASS",
                info,
                "GT:DR:DV",
                "0/1:13:7",
            )
            for index, info in enumerate(info_values)
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )

        self.assertEqual(report.accepted_record_count, 3)
        self.assertEqual(
            report.rejection_counts,
            {"conflicting_breakend_mate": 3, "malformed_breakend_mate": 1},
        )

    def test_chr2_end_fallback_rejects_other_non_bracket_alts(self) -> None:
        vcf = VCF_HEADER + "".join(
            _record(
                "chr1",
                str(1000 + index),
                ".",
                "N",
                alternate,
                "45",
                "PASS",
                "PRECISE;SVTYPE=BND;CHR2=chr2;END=5000;SUPPORT=7",
                "GT:DR:DV",
                "0/1:13:7",
            )
            for index, alternate in enumerate((".", "N", "<DEL>"))
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.rejection_counts, {"unsupported_breakend_alt": 3})

    def test_malformed_bracketed_breakend_is_rejected_instead_of_partially_matched(self) -> None:
        vcf = VCF_HEADER + _record(
            "chr1",
            "1000",
            ".",
            "N",
            "N[chr2:5000[N",
            "45",
            "PASS",
            "PRECISE;SVTYPE=BND;SUPPORT=7",
            "GT:DR:DV",
            "0/1:13:7",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.rejection_counts, {"malformed_breakend_alt": 1})

    def test_extreme_breakend_position_is_counted_not_raised(self) -> None:
        alternate = f"N[chr2:{'9' * 4301}["
        vcf = VCF_HEADER + _record(
            "chr1",
            "1000",
            ".",
            "N",
            alternate,
            "45",
            "PASS",
            "PRECISE;SVTYPE=BND;SUPPORT=7",
            "GT:DR:DV",
            "0/1:13:7",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.rejection_counts, {"malformed_breakend_alt": 1})
        self.assertNotIn(alternate, report.model_dump_json())

    def test_no_accepted_calls_is_no_call_not_negative(self) -> None:
        vcf = VCF_HEADER + _record(
            "chr1",
            "1000",
            ".",
            "N",
            "<DEL>",
            "60",
            "PASS",
            "SVTYPE=DEL;END=1200;SUPPORT=1",
            "GT:DR:DV",
            "0/1:10:1",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )
        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertTrue(any("not a clinical negative" in item for item in report.warnings))

    def test_runner_preserves_filtered_audit_evidence_but_accepts_pass_only(self) -> None:
        vcf = VCF_HEADER + _record(
            "chr1",
            "1000",
            ".",
            "N",
            "<DEL>",
            "60",
            "PASS",
            "SVTYPE=DEL;END=1200;SUPPORT=8",
            "GT:DR:DV",
            "0/1:8:8",
        )
        runner = RecordingRunner(vcf)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "synthetic.bam"
            bai = root / "synthetic.bam.bai"
            bam.write_bytes(b"synthetic")
            bai.write_bytes(b"synthetic")
            output = root / "calls.vcf"
            report = run_sniffles(
                _manifest(bam, bai),
                _intake(),
                _policy(),
                output_vcf=output,
                runner=runner,
                threads=2,
            )

        command = runner.commands[1]
        self.assertNotIn("--pass-only", command)
        self.assertIn("--symbolic", command)
        self.assertIn("--no-progress", command)
        self.assertNotIn("--output-rnames", command)
        self.assertEqual(report.tool.version, "2.8.0")
        self.assertFalse(report.tool.parameters["output_read_names"])
        self.assertFalse(report.tool.parameters["caller_pass_only"])
        self.assertTrue(report.tool.parameters["normalizer_pass_only"])
        self.assertNotIn("--allow-overwrite", command)

    def test_policy_rejects_unlocked_sniffles_version(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(VCF_HEADER, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "does not match policy lock"):
                normalize_sniffles_vcf(
                    path,
                    sample_id="SYNTHETIC_001",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                    tool=ToolRecord(name="Sniffles2", version="2.7.0"),
                )

    def test_normalizer_rechecks_minimum_sv_length(self) -> None:
        vcf = VCF_HEADER + _record(
            "chr1",
            "1000",
            ".",
            "N",
            "<DEL>",
            "60",
            "PASS",
            "SVTYPE=DEL;END=1020;SVLEN=-20;SUPPORT=8",
            "GT:DR:DV",
            "0/1:8:8",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(vcf, encoding="utf-8")
            report = normalize_sniffles_vcf(
                path,
                sample_id="SYNTHETIC_001",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
                tool=ToolRecord(name="Sniffles2", version="2.8.0"),
            )
        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.rejection_counts, {"sv_length_below_policy": 1})


if __name__ == "__main__":
    unittest.main()
