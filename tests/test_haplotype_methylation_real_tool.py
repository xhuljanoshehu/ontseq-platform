"""Pinned modkit 0.6.4 ``--phased`` behaviour on synthetic haplotagged MM/ML BAMs.

These tests are the evidence behind the lane's layout and refusal rules. They prove tool
interoperability on generated reads; they say nothing about phasing accuracy or about
allele-specific methylation in biological samples.
"""

from __future__ import annotations

import array
import importlib
import os
import shutil
import tempfile
import unittest
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ontseq_platform.execution import CommandResult, SubprocessRunner
from ontseq_platform.methylation_lanes.haplotype import (
    HaplotypeAssessability,
    HaplotypeMethylationPolicy,
    HaplotypeMethylationReport,
    NotAssessableReason,
    run_haplotype_methylation,
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

HAPLOTAG_PG = {
    "ID": "whatshap",
    "PN": "whatshap",
    "VN": "2.3",
    "CL": "whatshap haplotag --reference synthetic.fa phased.vcf.gz input.bam",
}


@dataclass(frozen=True)
class _Read:
    start: int
    methylated: bool
    hp: int | None = None
    ps: int | None = None
    reverse: bool = False
    #: Soft-clip the first CpG triplet; its MM call must not be placed on the reference.
    clipped: bool = False


class _RecordingRunner(SubprocessRunner):
    def __init__(self) -> None:
        self.commands: list[tuple[str, ...]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        self.commands.append(tuple(argv))
        return super().run(argv, timeout_seconds=timeout_seconds)


@unittest.skipUnless(
    os.environ.get("ONTSEQ_MODKIT_REAL_TOOL") == "1",
    "opt-in pinned modkit 0.6.4 interoperability lane",
)
class PhasedPileupTests(unittest.TestCase):
    def setUp(self) -> None:
        # An explicitly enabled CI lane must fail, not skip, when a tool is missing.
        self.assertIsNotNone(shutil.which("modkit"))
        self.assertIsNotNone(shutil.which("samtools"))
        self.pysam = importlib.import_module("pysam")
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.reference = self.root / "synthetic.fa"
        self.reference.write_text(">chr1\n" + "ACG" * 100 + "\n", encoding="utf-8")
        self.pysam.faidx(str(self.reference))
        self.regions = self.root / "regions.bed"
        self.regions.write_text(
            "chr1\t0\t30\tREGION_A\nchr1\t60\t90\tREGION_B\nchr1\t120\t150\tREGION_C\n",
            encoding="utf-8",
        )
        self.policy = HaplotypeMethylationPolicy(
            profile_id="SYNTHETIC_PHASED",
            status="technical_defaults_only",
            note="Synthetic interoperability only",
        )

    def _bam(
        self, name: str, reads: Sequence[_Read], *, programs: Sequence[dict[str, str]]
    ) -> Path:
        bam = self.root / f"{name}.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 300}],
            "PG": list(programs),
        }
        with self.pysam.AlignmentFile(str(bam), "wb", header=header) as target:
            for number, read in enumerate(
                sorted(reads, key=lambda item: item.start + 3 * item.clipped)
            ):
                record = self.pysam.AlignedSegment()
                record.query_name = f"SYNTHETIC_{number}"
                record.query_sequence = "ACG" * 10
                record.query_qualities = array.array("B", [40] * 30)
                record.flag = 16 if read.reverse else 0
                record.reference_id = 0
                record.reference_start = read.start + (3 if read.clipped else 0)
                record.mapping_quality = 60
                record.cigarstring = "3S27M" if read.clipped else "30M"
                record.set_tag("MN", 30)
                record.set_tag("MM", "C+m?," + ",".join(["0"] * 10) + ";")
                record.set_tag("ML", array.array("B", [255 if read.methylated else 0] * 10))
                if read.hp is not None:
                    record.set_tag("HP", read.hp)
                if read.ps is not None:
                    record.set_tag("PS", read.ps)
                target.write(record)
        self.pysam.index(str(bam))
        return bam

    def _inputs(self, bam: Path) -> tuple[SampleManifest, AlignedBamIntakeReport]:
        manifest = SampleManifest(
            sample_id="SYNTHETIC_PHASED",
            run_id="SYNTHETIC_RUN",
            input=InputSpec(kind=InputKind.ALIGNED_BAM, path=str(bam), index_path=f"{bam}.bai"),
            assay=AssaySpec(
                mode=AssayMode.LOW_COVERAGE_WGS,
                genome_build=GenomeBuild.GRCH38,
                reference_id="synthetic-not-a-human-reference",
            ),
            analysis=AnalysisSpec(profile="synthetic", modules=[AnalysisModule.METHYLATION]),
        )
        intake = AlignedBamIntakeReport(
            sample_id=manifest.sample_id,
            reference_id=manifest.assay.reference_id,
            genome_build=GenomeBuild.GRCH38,
            checks=[],
            verdict=Verdict.PASS,
        )
        return manifest, intake

    def _run(self, bam: Path, runner: _RecordingRunner | None = None) -> HaplotypeMethylationReport:
        manifest, intake = self._inputs(bam)
        return run_haplotype_methylation(
            manifest,
            intake,
            self.policy,
            region_bed=self.regions,
            output_dir=self.root / "out",
            reference_fasta=self.reference,
            runner=runner,
            threads=1,
        )

    def test_phased_pileup_splits_haplotypes_inside_one_phase_block_only(self) -> None:
        reads = [
            # REGION_A: one phase block. HP1 fully methylated on both strands, HP2 not,
            # and six untagged reads, half of them methylated.
            *(_Read(0, True, hp=1, ps=100) for _ in range(6)),
            *(_Read(0, True, hp=1, ps=100, reverse=True) for _ in range(2)),
            *(_Read(0, False, hp=2, ps=100) for _ in range(6)),
            *(_Read(0, False, hp=2, ps=100, clipped=True) for _ in range(2)),
            *(_Read(0, number % 2 == 0) for number in range(6)),
            # REGION_B: HP1 from two phase blocks with opposite states, which the pinned
            # binary pools into one "HP1".
            *(_Read(60, True, hp=1, ps=100) for _ in range(3)),
            *(_Read(60, False, hp=1, ps=200) for _ in range(3)),
            *(_Read(60, True, hp=2, ps=200) for _ in range(6)),
            # REGION_C: untagged reads only.
            *(_Read(120, True) for _ in range(5)),
        ]
        runner = _RecordingRunner()
        report = self._run(self._bam("phased", reads, programs=[HAPLOTAG_PG]), runner)

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.reads_with_modified_base_tags, len(reads))
        self.assertEqual(report.reads_with_haplotype_tags, 28)
        self.assertEqual([item.program for item in report.phasing], ["whatshap"])
        self.assertEqual(report.tool.parameters["phased"], True)
        self.assertEqual(
            sorted(
                path.name for path in (self.root / "out/SYNTHETIC_PHASED.modkit-phased").iterdir()
            ),
            ["combined.bedmethyl", "hp1.bedmethyl", "hp2.bedmethyl"],
        )
        pileups = [command for command in runner.commands if "pileup" in command]
        self.assertEqual(len(pileups), 1)
        self.assertIn("--phased", pileups[0])
        rows = {row.region_id: row for row in report.rows}

        region_a = rows["REGION_A"]
        self.assertEqual(region_a.phase_sets, [100])
        self.assertEqual(region_a.hp1.informative_reads, 8)
        self.assertEqual(region_a.hp1.summary.valid_call_count, 80)
        self.assertEqual(region_a.hp1.summary.mean_modified_fraction, 1.0)
        # Two soft-clipped reads add calls to the nine CpGs they align to, not to the
        # first CpG, whose call lies in the clipped bases.
        self.assertEqual(region_a.hp2.informative_reads, 8)
        self.assertEqual(region_a.hp2.summary.valid_call_count, 78)
        self.assertEqual(region_a.hp2.summary.sites_total, 10)
        self.assertEqual(region_a.hp2.summary.mean_modified_fraction, 0.0)
        self.assertIs(region_a.hp1.assessability, HaplotypeAssessability.ASSESSABLE)
        self.assertIs(region_a.hp2.assessability, HaplotypeAssessability.ASSESSABLE)
        self.assertEqual(region_a.haplotype_difference, 1.0)
        # Unphased calls are combined minus both haplotypes, never dropped.
        self.assertEqual(region_a.unphased_informative_reads, 6)
        self.assertEqual(region_a.unphased.valid_call_count, 60)
        self.assertEqual(region_a.unphased.modified_call_count, 30)

        region_b = rows["REGION_B"]
        self.assertEqual(region_b.phase_sets, [100, 200])
        # The pooled HP1 is still reported, but it is never compared.
        self.assertEqual(region_b.hp1.summary.mean_modified_fraction, 0.5)
        for call in (region_b.hp1, region_b.hp2):
            self.assertIs(call.assessability, HaplotypeAssessability.NOT_ASSESSABLE)
            self.assertIn(NotAssessableReason.MULTIPLE_PHASE_BLOCKS, call.reasons)
        self.assertIsNone(region_b.haplotype_difference)
        self.assertEqual(region_b.unphased.sites_total, 0)

        region_c = rows["REGION_C"]
        self.assertEqual(region_c.phase_sets, [])
        self.assertIn(NotAssessableReason.NO_PHASED_READS, region_c.hp1.reasons)
        self.assertEqual(region_c.hp1.summary.sites_total, 0)
        self.assertIsNone(region_c.hp1.summary.mean_modified_fraction)
        self.assertEqual(region_c.unphased.mean_modified_fraction, 1.0)
        self.assertIsNone(region_c.haplotype_difference)

    def _assert_refused_before_pileup(self, bam: Path, message: str) -> None:
        runner = _RecordingRunner()
        with self.assertRaisesRegex(ValueError, message):
            self._run(bam, runner)
        self.assertFalse(any("pileup" in command for command in runner.commands))
        self.assertFalse((self.root / "out/SYNTHETIC_PHASED.modkit-phased").exists())

    def test_haplotype_labels_other_than_one_and_two_are_refused(self) -> None:
        for value in (3, 0):
            with self.subTest(hp=value):
                bam = self._bam(
                    f"hp{value}",
                    [_Read(0, True, hp=1, ps=100), _Read(0, True, hp=value, ps=100)],
                    programs=[HAPLOTAG_PG],
                )
                self._assert_refused_before_pileup(bam, "HP value other than 1 or 2")

    def test_haplotagged_reads_without_a_phase_set_are_refused(self) -> None:
        bam = self._bam(
            "no-ps",
            [_Read(0, True, hp=1, ps=100), _Read(0, False, hp=2)],
            programs=[HAPLOTAG_PG],
        )
        self._assert_refused_before_pileup(bam, "carry no PS phase set")

    def test_haplotype_tags_of_unknown_provenance_are_refused(self) -> None:
        reads = [_Read(0, True, hp=1, ps=100), _Read(0, False, hp=2, ps=100)]
        not_haplotag = {**HAPLOTAG_PG, "CL": "whatshap phase -o phased.vcf.gz input.vcf.gz"}
        for name, programs in (("no-pg", []), ("phase-only", [not_haplotag])):
            with self.subTest(header=name):
                bam = self._bam(name, reads, programs=programs)
                self._assert_refused_before_pileup(bam, "no haplotagging step")

    def test_a_bam_without_haplotype_tags_is_refused(self) -> None:
        bam = self._bam("untagged", [_Read(0, True)], programs=[HAPLOTAG_PG])
        self._assert_refused_before_pileup(bam, "no HP haplotype tags")


if __name__ == "__main__":
    unittest.main()
