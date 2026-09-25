"""The haplotype lane's contract without real tools: refusals, arithmetic and assessability.

The pinned-binary behaviour these tests assume is exercised separately in
``tests/test_haplotype_methylation_real_tool.py``.
"""

from __future__ import annotations

import tempfile
import unittest
from collections.abc import Mapping, Sequence
from pathlib import Path

from pydantic import ValidationError

from ontseq_platform.execution import CommandResult
from ontseq_platform.io import load_model
from ontseq_platform.methylation import ModificationCode, _Region, _Site
from ontseq_platform.methylation_lanes.haplotype import (
    PHASED_OUTPUT_FILES,
    HaplotypeAssessability,
    HaplotypeMethylationPolicy,
    HaplotypeMethylationReport,
    NotAssessableReason,
    PhasingProvenance,
    RegionPhasing,
    run_haplotype_methylation,
    split_unphased,
)
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    ModuleRunStatus,
    SampleManifest,
    Verdict,
)

ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = [PhasingProvenance(program="whatshap", version="2.3", command_line_sha256="0" * 64)]
REGION_A = _Region(region_id="REGION_A", chromosome="chr1", start=0, end=30)
REGION_B = _Region(region_id="REGION_B", chromosome="chr1", start=60, end=90)


def _policy(**overrides: object) -> HaplotypeMethylationPolicy:
    values: dict[str, object] = {
        "profile_id": "synthetic",
        "status": "technical_defaults_only",
        "note": "Synthetic technical policy",
    }
    values.update(overrides)
    return HaplotypeMethylationPolicy.model_validate(values)


def _row(start: int, valid: int, modified: int, *, code: str = "m", fail: int = 0) -> str:
    percent = modified / valid * 100 if valid else 0.0
    return "\t".join(
        str(item)
        for item in (
            "chr1",
            start,
            start + 1,
            code,
            valid,
            ".",
            start,
            start + 1,
            "255,0,0",
            valid,
            f"{percent:.2f}",
            modified,
            valid - modified,
            0,
            0,
            fail,
            0,
            0,
        )
    )


def _site(start: int, valid: int, modified: int) -> _Site:
    return _Site(
        chromosome="chr1",
        start=start,
        end=start + 1,
        code=ModificationCode.FIVE_MC,
        valid_coverage=valid,
        modified_calls=modified,
        canonical_calls=valid - modified,
        other_mod_calls=0,
        fail_calls=0,
        nocall_calls=0,
        delete_calls=0,
        diff_calls=0,
    )


class _FakeRunner:
    """Stand in for modkit and samtools, writing the layout the pinned binary produces."""

    def __init__(
        self,
        pileups: Mapping[str, Sequence[str]],
        *,
        counts: Mapping[str, str] | None = None,
        version: str = "0.6.4",
        failed_processing: int = 0,
        extra_file: str | None = None,
    ) -> None:
        self.pileups = pileups
        self.counts = {
            "[MM]": "40",
            "exists([HP])": "26",
            "exists([HP]) && [HP] != 1 && [HP] != 2": "0",
            "exists([HP]) && !exists([PS])": "0",
            **(counts or {}),
        }
        self.version = version
        self.failed_processing = failed_processing
        self.extra_file = extra_file
        self.commands: list[list[str]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        items = [str(item) for item in argv]
        self.commands.append(items)
        if items[1:2] == ["--version"]:
            return CommandResult(tuple(items), 0, f"mod_kit {self.version}", "")
        if "view" in items:
            expression = items[items.index("-e") + 1]
            count = "0" if "C[+-]" in expression else self.counts[expression]
            return CommandResult(tuple(items), 0, count, "")
        if "pileup" in items:
            output = Path(items[3])
            output.mkdir()
            for name in PHASED_OUTPUT_FILES:
                rows = self.pileups.get(name.split(".")[0], [])
                (output / name).write_text("".join(f"{row}\n" for row in rows), encoding="utf-8")
            if self.extra_file:
                (output / self.extra_file).write_text("", encoding="utf-8")
            if self.failed_processing:
                log = Path(items[items.index("--log-filepath") + 1])
                log.write_text(f"~{self.failed_processing} failed processing\n", encoding="utf-8")
            return CommandResult(tuple(items), 0, "", "")
        raise AssertionError(f"unexpected command: {items}")


def _scan(phasing: Mapping[str, RegionPhasing]):  # type: ignore[no-untyped-def]
    def scan(bam: Path, regions: Sequence[_Region]) -> dict[_Region, RegionPhasing]:
        return {region: phasing[region.region_id] for region in regions}

    return scan


def _provenance(bam: Path, accepted: Sequence[str]) -> list[PhasingProvenance]:
    return PROVENANCE


#: Region A: HP1 8 reads methylated, HP2 6 reads unmethylated, 6 unphased reads half
#: methylated; three CpGs. Region B: HP1 only at depth, in two phase blocks.
PILEUPS = {
    "combined": [_row(1, 20, 11), _row(4, 20, 11), _row(7, 20, 11), _row(61, 9, 9)],
    "hp1": [_row(1, 8, 8), _row(4, 8, 8), _row(7, 8, 8), _row(61, 6, 6)],
    "hp2": [_row(1, 6, 0), _row(4, 6, 0), _row(7, 6, 0), _row(61, 3, 3)],
}
PHASING = {
    "REGION_A": RegionPhasing(frozenset({100}), hp1_reads=8, hp2_reads=6, unphased_reads=6),
    "REGION_B": RegionPhasing(frozenset({100, 200}), hp1_reads=6, hp2_reads=3),
}


class _Case(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        bam = self.root / "sample.bam"
        bam.write_bytes(b"not a real bam")
        self.reference = self.root / "reference.fa"
        self.reference.write_text(">chr1\nACG\n", encoding="utf-8")
        self.regions = self.root / "regions.bed"
        self.regions.write_text(
            "chr1\t0\t30\tREGION_A\tscore\nchr1\t60\t90\tREGION_B\tscore\n", encoding="utf-8"
        )
        self.manifest = SampleManifest(
            sample_id="SYNTHETIC_HAPLOTYPE",
            run_id="RUN_001",
            input=InputSpec(kind=InputKind.ALIGNED_BAM, path=str(bam), index_path=f"{bam}.bai"),
            assay=AssaySpec(
                mode=AssayMode.LOW_COVERAGE_WGS,
                genome_build=GenomeBuild.GRCH38,
                reference_id="synthetic-grch38",
            ),
            analysis=AnalysisSpec(profile="synthetic", modules=[AnalysisModule.METHYLATION]),
        )
        self.intake = AlignedBamIntakeReport(
            sample_id="SYNTHETIC_HAPLOTYPE",
            reference_id="synthetic-grch38",
            genome_build=GenomeBuild.GRCH38,
            checks=[],
            verdict=Verdict.PASS,
        )

    def _run(
        self,
        runner: _FakeRunner,
        *,
        policy: HaplotypeMethylationPolicy | None = None,
        phasing: Mapping[str, RegionPhasing] = PHASING,
    ) -> HaplotypeMethylationReport:
        return run_haplotype_methylation(
            self.manifest,
            self.intake,
            policy or _policy(),
            region_bed=self.regions,
            output_dir=self.root / "out",
            reference_fasta=self.reference,
            runner=runner,
            modkit="modkit-not-installed-in-unit-tests",
            threads=1,
            scan=_scan(phasing),
            provenance=_provenance,
        )

    @property
    def phased_dir(self) -> Path:
        return self.root / "out/SYNTHETIC_HAPLOTYPE.modkit-phased"


class PolicyTests(unittest.TestCase):
    def test_the_shipped_technical_policy_loads_and_is_not_validated(self) -> None:
        policy = load_model(
            ROOT / "configs/methylation/haplotype.technical.yaml", HaplotypeMethylationPolicy
        )
        self.assertEqual(policy.status, "technical_defaults_only")
        self.assertEqual(policy.expected_version, "0.6.4")
        self.assertEqual(policy.modification_codes, [ModificationCode.FIVE_MC])
        self.assertIn("not validated", policy.note)

    def test_inconsistent_policies_are_refused(self) -> None:
        for overrides in (
            {"modification_codes": ["m", "m"]},
            {"cpg_only": False, "combine_strands": True},
            {"accepted_phasing_programs": ["WhatsHap"]},
            {"accepted_phasing_programs": []},
            {"minimum_informative_reads_per_haplotype": 0},
        ):
            with self.subTest(overrides=overrides), self.assertRaises(ValidationError):
                _policy(**overrides)


class UnphasedArithmeticTests(unittest.TestCase):
    def test_unphased_calls_are_combined_minus_both_haplotypes(self) -> None:
        unphased = split_unphased(
            [_site(1, 20, 11), _site(4, 8, 8)],
            [_site(1, 8, 8), _site(4, 8, 8)],
            [_site(1, 6, 0)],
        )
        # Site 4 is fully explained by HP1 and leaves no unphased row at all.
        self.assertEqual(len(unphased), 1)
        self.assertEqual((unphased[0].valid_coverage, unphased[0].modified_calls), (6, 3))

    def test_a_haplotype_exceeding_the_combined_pileup_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "exceed the combined pileup"):
            split_unphased([_site(1, 5, 5)], [_site(1, 4, 4)], [_site(1, 2, 0)])

    def test_a_haplotype_site_missing_from_the_combined_pileup_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "missing from the combined pileup"):
            split_unphased([_site(1, 5, 5)], [_site(4, 4, 4)], [])


class AdapterTests(_Case):
    def test_a_single_phase_block_region_is_compared_and_a_multi_block_region_is_not(
        self,
    ) -> None:
        runner = _FakeRunner(PILEUPS)
        report = self._run(runner)
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        rows = {row.region_id: row for row in report.rows}

        region_a = rows["REGION_A"]
        self.assertEqual(region_a.hp1.summary.mean_modified_fraction, 1.0)
        self.assertEqual(region_a.hp2.summary.mean_modified_fraction, 0.0)
        self.assertEqual(region_a.haplotype_difference, 1.0)
        self.assertEqual(region_a.unphased.valid_call_count, 18)
        self.assertEqual(region_a.unphased.modified_call_count, 9)

        region_b = rows["REGION_B"]
        self.assertIsNone(region_b.haplotype_difference)
        self.assertEqual(
            region_b.hp1.reasons,
            [NotAssessableReason.MULTIPLE_PHASE_BLOCKS, NotAssessableReason.INSUFFICIENT_SITES],
        )
        self.assertEqual(
            region_b.hp2.reasons,
            [
                NotAssessableReason.MULTIPLE_PHASE_BLOCKS,
                NotAssessableReason.INSUFFICIENT_READS,
                NotAssessableReason.INSUFFICIENT_SITES,
            ],
        )
        self.assertEqual(report.summary_metrics["comparable_row_count"], 1)
        self.assertEqual(report.summary_metrics["haplotype_calls_multiple_phase_blocks"], 2)
        self.assertEqual(set(report.bedmethyl_fingerprints), {"combined", "hp1", "hp2"})
        self.assertTrue(report.research_only)

        pileup = next(command for command in runner.commands if "pileup" in command)
        self.assertIn("--phased", pileup)
        self.assertEqual(pileup[pileup.index("--filter-threshold") + 1], "0.8")
        self.assertEqual(pileup[pileup.index("--modified-bases") + 1], "5mC")
        self.assertEqual(pileup[-2:], ["--cpg", "--combine-strands"])
        # The pinned include-bed parser rejects BED4/5, so the design is projected to BED3.
        include = Path(pileup[pileup.index("--include-bed") + 1])
        self.assertEqual(include.read_text(encoding="utf-8"), "chr1\t0\t30\nchr1\t60\t90\n")

    def test_a_haplotype_below_the_assessability_rule_is_not_assessable(self) -> None:
        report = self._run(
            _FakeRunner(PILEUPS),
            policy=_policy(minimum_informative_reads_per_haplotype=7),
        )
        region_a = next(row for row in report.rows if row.region_id == "REGION_A")
        self.assertIs(region_a.hp1.assessability, HaplotypeAssessability.ASSESSABLE)
        self.assertEqual(region_a.hp2.reasons, [NotAssessableReason.INSUFFICIENT_READS])
        # The measured fraction stays visible; only the comparison is withheld.
        self.assertEqual(region_a.hp2.summary.mean_modified_fraction, 0.0)
        self.assertIsNone(region_a.haplotype_difference)
        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertTrue(any("does not say there is none" in item for item in report.warnings))

    def test_a_region_without_haplotagged_reads_is_not_assessable(self) -> None:
        report = self._run(
            _FakeRunner(PILEUPS),
            phasing={
                "REGION_A": RegionPhasing(unphased_reads=20),
                "REGION_B": PHASING["REGION_B"],
            },
        )
        region_a = next(row for row in report.rows if row.region_id == "REGION_A")
        self.assertIn(NotAssessableReason.NO_PHASED_READS, region_a.hp1.reasons)
        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)

    def test_tag_problems_are_refused_before_the_pileup_runs(self) -> None:
        cases = {
            "no MM modified-base tags": {"[MM]": "0"},
            "no HP haplotype tags": {"exists([HP])": "0"},
            "HP value other than 1 or 2": {"exists([HP]) && [HP] != 1 && [HP] != 2": "2"},
            "carry no PS phase set": {"exists([HP]) && !exists([PS])": "1"},
        }
        for message, counts in cases.items():
            with self.subTest(message=message):
                runner = _FakeRunner(PILEUPS, counts=counts)
                with self.assertRaisesRegex(ValueError, message):
                    self._run(runner)
                self.assertFalse(any("pileup" in command for command in runner.commands))

    def test_an_unevaluable_tag_expression_is_refused_rather_than_assumed_absent(self) -> None:
        runner = _FakeRunner(PILEUPS, counts={"exists([HP]) && !exists([PS])": "error"})
        with self.assertRaisesRegex(ValueError, "could not evaluate"):
            self._run(runner)

    def test_a_version_outside_the_lock_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "does not match policy lock"):
            self._run(_FakeRunner(PILEUPS, version="0.7.0"))

    def test_an_unqualified_output_layout_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "qualified combined/hp1/hp2 layout"):
            self._run(_FakeRunner(PILEUPS, extra_file="hp3.bedmethyl"))

    def test_failed_records_discard_the_partial_output(self) -> None:
        with self.assertRaisesRegex(ValueError, "failed record"):
            self._run(_FakeRunner(PILEUPS, failed_processing=3))
        self.assertFalse(self.phased_dir.exists())

    def test_existing_outputs_are_never_overwritten(self) -> None:
        self.phased_dir.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError, "Refusing to overwrite"):
            self._run(_FakeRunner(PILEUPS))

    def test_the_locked_reference_is_required(self) -> None:
        self.reference.unlink()
        with self.assertRaisesRegex(ValueError, "locked reference FASTA"):
            self._run(_FakeRunner(PILEUPS))


class ReportContractTests(_Case):
    def test_a_difference_across_phase_blocks_cannot_be_constructed(self) -> None:
        payload = self._run(_FakeRunner(PILEUPS)).model_dump(mode="json")
        payload["rows"][0]["phase_sets"] = [100, 200]
        with self.assertRaisesRegex(ValidationError, "single phase block"):
            HaplotypeMethylationReport.model_validate(payload)

    def test_the_status_follows_whether_any_region_was_comparable(self) -> None:
        payload = self._run(_FakeRunner(PILEUPS)).model_dump(mode="json")
        payload["status"] = "NO_CALL"
        with self.assertRaisesRegex(ValidationError, "comparable"):
            HaplotypeMethylationReport.model_validate(payload)

    def test_an_assessable_haplotype_needs_a_measured_fraction(self) -> None:
        payload = self._run(_FakeRunner(PILEUPS)).model_dump(mode="json")
        row = payload["rows"][1]
        row["hp1"]["assessability"] = "ASSESSABLE"
        row["hp1"]["reasons"] = []
        row["hp1"]["summary"]["mean_modified_fraction"] = None
        with self.assertRaises(ValidationError):
            HaplotypeMethylationReport.model_validate(payload)

    def test_no_local_path_reaches_the_report(self) -> None:
        text = self._run(_FakeRunner(PILEUPS)).model_dump_json()
        self.assertNotIn("sample.bam", text)
        self.assertNotIn(str(self.root), text)


if __name__ == "__main__":
    unittest.main()
