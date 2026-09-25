"""Reading HP/PS tags and the haplotagging header with pysam, on generated fixtures only."""

from __future__ import annotations

import array
import tempfile
import unittest
from importlib import import_module
from pathlib import Path
from unittest import SkipTest

from ontseq_platform.methylation import _Region
from ontseq_platform.methylation_lanes.haplotype import (
    RegionPhasing,
    read_phasing_provenance,
    scan_region_phasing,
)

try:
    pysam = import_module("pysam")
except ModuleNotFoundError as error:
    if error.name != "pysam":
        raise
    raise SkipTest("Haplotype tag scanning requires pysam in the analysis runtime") from error

COMMAND_LINE = "whatshap haplotag --reference /private/site/ref.fa calls.vcf.gz in.bam"
REGION = _Region(region_id="REGION", chromosome="chr1", start=0, end=30)


class PhasingScanTests(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)

    def _bam(
        self,
        records: list[tuple[int, dict[str, int]]],
        *,
        programs: list[dict[str, str]] | None = None,
    ) -> Path:
        bam = self.root / "phasing.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 300}],
            "PG": programs
            if programs is not None
            else [{"ID": "whatshap", "PN": "whatshap", "VN": "2.3", "CL": COMMAND_LINE}],
        }
        with pysam.AlignmentFile(str(bam), "wb", header=header) as target:
            for number, (flag, tags) in enumerate(records):
                record = pysam.AlignedSegment()
                record.query_name = f"SYNTHETIC_{number}"
                record.query_sequence = "ACG" * 10
                record.query_qualities = array.array("B", [40] * 30)
                record.flag = flag
                record.reference_id = 0
                record.reference_start = 0
                record.mapping_quality = 60
                record.cigarstring = "30M"
                for name, value in tags.items():
                    record.set_tag(name, value)
                target.write(record)
        pysam.index(str(bam))
        return bam

    def test_only_informative_reads_are_counted_per_haplotype(self) -> None:
        bam = self._bam(
            [
                (0, {"HP": 1, "PS": 100}),
                (16, {"HP": 1, "PS": 100}),
                (0, {"HP": 2, "PS": 100}),
                (0, {}),
                # secondary, supplementary, QC-failed and duplicate records are not
                # informative reads, whatever their tags say
                (0x100, {"HP": 2, "PS": 300}),
                (0x800, {"HP": 2, "PS": 300}),
                (0x200, {"HP": 2, "PS": 300}),
                (0x400, {"HP": 2, "PS": 300}),
            ]
        )
        self.assertEqual(
            scan_region_phasing(bam, [REGION]),
            {REGION: RegionPhasing(frozenset({100}), hp1_reads=2, hp2_reads=1, unphased_reads=1)},
        )

    def test_unsupported_tags_are_refused_while_scanning(self) -> None:
        for tags, message in (
            ({"HP": 3, "PS": 100}, "HP value other than 1 or 2"),
            ({"HP": 1}, "no PS phase set"),
        ):
            with self.subTest(tags=tags):
                bam = self._bam([(0, tags)])
                with self.assertRaisesRegex(ValueError, message):
                    scan_region_phasing(bam, [REGION])

    def test_the_haplotagging_step_is_recorded_without_its_command_line(self) -> None:
        provenance = read_phasing_provenance(self._bam([]), ["whatshap", "longphase"])
        self.assertEqual(len(provenance), 1)
        self.assertEqual((provenance[0].program, provenance[0].version), ("whatshap", "2.3"))
        self.assertNotIn("/private/site", provenance[0].model_dump_json())
        self.assertEqual(len(provenance[0].command_line_sha256), 64)

    def test_haplotype_tags_of_unknown_provenance_are_refused(self) -> None:
        for programs in (
            [],
            [{"ID": "whatshap", "PN": "whatshap", "CL": "whatshap phase -o out.vcf in.vcf"}],
            [{"ID": "custom", "PN": "custom-phaser", "CL": "custom-phaser haplotag in.bam"}],
        ):
            with (
                self.subTest(programs=programs),
                self.assertRaisesRegex(ValueError, "no haplotagging step"),
            ):
                read_phasing_provenance(self._bam([], programs=programs), ["whatshap"])
