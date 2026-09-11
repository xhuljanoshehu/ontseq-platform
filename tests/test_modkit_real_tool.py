"""Pinned-binary interoperability on generated signals, not patient sensitivity/specificity."""

from __future__ import annotations

import array
import importlib
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.methylation import MethylationPolicy, run_methylation
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


@unittest.skipUnless(
    os.environ.get("ONTSEQ_MODKIT_REAL_TOOL") == "1",
    "opt-in pinned modkit 0.6.4 interoperability lane",
)
class ModkitBinaryTests(unittest.TestCase):
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
        self.bed = self.root / "synthetic.bed"
        self.bed.write_text(
            "chr1\t0\t30\tGENE\nchr1\t60\t90\tGENE\nchr1\t120\t150\tLOW\nchr1\t240\t270\tEMPTY\n",
            encoding="utf-8",
        )

    def _inputs(self, tags: bool = True) -> tuple[SampleManifest, AlignedBamIntakeReport]:
        bam = self.root / ("with-tags.bam" if tags else "without-tags.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 300}]}
        with self.pysam.AlignmentFile(str(bam), "wb", header=header) as target:
            for start, probabilities in [
                (0, [255] * 15 + [0] * 5),
                (60, [0] * 10),
                (120, [255] * 2),
            ]:
                for number, probability in enumerate(probabilities):
                    record = self.pysam.AlignedSegment()
                    record.query_name = f"SYNTHETIC_{start}_{number}"
                    record.query_sequence = "ACG" * 10
                    record.query_qualities = array.array("B", [40] * 30)
                    record.flag = 0
                    record.reference_id = 0
                    record.reference_start = start
                    record.mapping_quality = 60
                    record.cigarstring = "30M"
                    if tags:
                        # MN is the sequence length at tag-writing time; modkit 0.6.x
                        # requires correct MN tags for --combine-strands.
                        record.set_tag("MN", 30)
                        record.set_tag("MM", "C+m?," + ",".join(["0"] * 10) + ";")
                        record.set_tag("ML", array.array("B", [probability] * 10))
                    target.write(record)
        self.pysam.index(str(bam))
        manifest = SampleManifest(
            sample_id="SYNTHETIC_MODKIT",
            run_id="SYNTHETIC_RUN",
            input=InputSpec(
                kind=InputKind.ALIGNED_BAM, path=str(bam), index_path=str(bam) + ".bai"
            ),
            assay=AssaySpec(
                mode=AssayMode.ADAPTIVE_SAMPLING,
                genome_build=GenomeBuild.GRCH38,
                reference_id="synthetic-not-a-human-reference",
                target_bed=str(self.bed),
                target_bed_version="SYNTHETIC_V1",
            ),
            analysis=AnalysisSpec(profile="synthetic", modules=[AnalysisModule.METHYLATION]),
        )
        intake = AlignedBamIntakeReport(
            sample_id=manifest.sample_id,
            reference_id=manifest.assay.reference_id,
            genome_build=GenomeBuild.GRCH38,
            checks=[],
            verdict=Verdict.PASS,
            limitations=["Synthetic adapter input only, not a full-reference intake validation."],
        )
        return manifest, intake

    def _policy(self) -> MethylationPolicy:
        return MethylationPolicy(
            profile_id="SYNTHETIC_MODKIT",
            expected_version="0.6.4",
            status="technical_defaults_only",
            minimum_valid_coverage=5,
            cpg_only=True,
            combine_strands=True,
            filter_threshold=0.8,
            region_source="target_bed",
            note="Synthetic interoperability only",
        )

    def test_real_pileup_preserves_fraction_zero_missing_and_repeated_labels(self) -> None:
        manifest, intake = self._inputs()
        report = run_methylation(
            manifest,
            intake,
            self._policy(),
            output_dir=self.root / "pileup",
            reference_fasta=self.reference,
            threads=1,
        )
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.reads_with_modified_base_tags, 32)
        rows = {row.start: row for row in report.regions}
        self.assertEqual(rows[0].valid_call_count, 200)
        self.assertEqual(rows[0].modified_call_count, 150)
        self.assertEqual(rows[0].mean_modified_fraction, 0.75)
        self.assertEqual(rows[60].valid_call_count, 100)
        self.assertEqual(rows[60].mean_modified_fraction, 0.0)
        self.assertIsNone(rows[120].mean_modified_fraction)
        self.assertIsNone(rows[240].mean_modified_fraction)
        self.assertEqual(report.tool.parameters["ontseq_region_assignment"], "interval-identity-v1")
        self.assertEqual(report.tool.version, "0.6.4")

    def test_real_samtools_detects_independent_same_base_groups_before_modkit(self) -> None:
        reference = self.root / "independent.fa"
        reference.write_text(">chr1\nCC\n", encoding="utf-8")
        self.pysam.faidx(str(reference))
        bam = self.root / "independent.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 2}],
        }
        with self.pysam.AlignmentFile(str(bam), "wb", header=header) as target:
            record = self.pysam.AlignedSegment()
            record.query_name = "INDEPENDENT_GROUPS"
            record.query_sequence = "CC"
            record.query_qualities = array.array("B", [40, 40])
            record.flag = 0
            record.reference_id = 0
            record.reference_start = 0
            record.mapping_quality = 60
            record.cigarstring = "2M"
            record.set_tag("MN", 2)
            record.set_tag("MM", "C+m?,0;C+h?,1;")
            record.set_tag("ML", array.array("B", [200, 250]))
            target.write(record)
        self.pysam.index(str(bam))
        manifest = SampleManifest(
            sample_id="INDEPENDENT_MODKIT",
            run_id="SYNTHETIC_RUN",
            input=InputSpec(
                kind=InputKind.ALIGNED_BAM,
                path=str(bam),
                index_path=str(bam) + ".bai",
            ),
            assay=AssaySpec(
                mode=AssayMode.LOW_COVERAGE_WGS,
                genome_build=GenomeBuild.GRCH38,
                reference_id="synthetic-independent-reference",
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
        policy = MethylationPolicy(
            profile_id="SYNTHETIC_INDEPENDENT_GROUPS",
            expected_version="0.6.4",
            status="technical_defaults_only",
            modification_codes=["m", "h"],
            cpg_only=False,
            combine_strands=False,
            region_source="chromosome",
            note="Regression for the modkit 0.6.4 independent-group safety gate",
        )

        with self.assertRaisesRegex(ValueError, "independent cytosine MM groups"):
            run_methylation(
                manifest,
                intake,
                policy,
                output_dir=self.root / "independent-output",
                reference_fasta=reference,
                threads=1,
            )
        self.assertFalse(
            (self.root / "independent-output" / "INDEPENDENT_MODKIT.modkit.bedmethyl").exists()
        )

    def test_real_bam_without_modification_tags_never_becomes_unmethylated(self) -> None:
        manifest, intake = self._inputs(tags=False)
        with self.assertRaisesRegex(ValueError, "no MM modified-base tags"):
            run_methylation(
                manifest,
                intake,
                self._policy(),
                output_dir=self.root / "missing",
                reference_fasta=self.reference,
                threads=1,
            )
        self.assertFalse((self.root / "missing" / "SYNTHETIC_MODKIT.modkit.bedmethyl").exists())
