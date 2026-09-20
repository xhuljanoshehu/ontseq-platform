"""Pinned Spectre interoperability on generated depth signals, never patient data."""

from __future__ import annotations

import importlib
import os
import shutil
import tempfile
import unittest
from pathlib import Path

import yaml

from ontseq_platform.cnv.spectre import SpectrePolicy, run_spectre_depth_only
from ontseq_platform.models import EventType, GenomeBuild, ModuleRunStatus


@unittest.skipUnless(
    os.environ.get("ONTSEQ_SPECTRE_REAL_TOOL") == "1",
    "opt-in pinned Spectre 0.2.1 depth-only interoperability lane",
)
class SpectreBinaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertIsNotNone(shutil.which("spectre"))
        self.pysam = importlib.import_module("pysam")
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.reference = self.root / "synthetic.reference.fa"
        with self.reference.open("w", encoding="ascii") as handle:
            handle.write(">chr7\n")
            sequence_line = "A" * 80
            for _ in range(2_000_000 // len(sequence_line)):
                handle.write(f"{sequence_line}\n")
        self.pysam.faidx(str(self.reference))

        policy_path = (
            Path(__file__).parents[1] / "configs" / "cnv" / "spectre.depth-only.technical.yaml"
        )
        self.policy = SpectrePolicy.model_validate(
            yaml.safe_load(policy_path.read_text(encoding="utf-8"))
        )

    def _coverage(self, name: str, *, altered: bool) -> Path:
        bed = self.root / f"{name}.regions.bed"
        with bed.open("w", encoding="ascii") as handle:
            for start in range(0, 2_000_000, self.policy.mosdepth_bin_size_bp):
                end = start + self.policy.mosdepth_bin_size_bp
                depth = 30
                if altered and 500_000 <= start < 800_000:
                    depth = 12
                elif altered and 1_100_000 <= start < 1_400_000:
                    depth = 60
                handle.write(f"chr7\t{start}\t{end}\t{depth}\n")
        compressed = bed.with_suffix(".bed.gz")
        self.pysam.tabix_compress(str(bed), str(compressed), force=True)
        self.pysam.tabix_index(str(compressed), preset="bed", force=True, csi=True)
        return compressed

    def test_real_depth_only_calls_and_empty_output_are_preserved(self) -> None:
        altered_output = self.root / "altered-output"
        report = run_spectre_depth_only(
            coverage_path=self._coverage("altered", altered=True),
            sample_id="SYNTHETIC_DEPTH_ONLY",
            output_dir=altered_output,
            reference_fasta=self.reference,
            genome_build=GenomeBuild.GRCH38,
            policy=self.policy,
            threads=1,
        )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.tool.version, "0.2.1")
        self.assertEqual(
            [event.event_type for event in report.events],
            [EventType.DELETION, EventType.DUPLICATION],
        )
        self.assertEqual(
            [(event.primary.start, event.primary.end) for event in report.events],
            [(495_000, 800_000), (1_095_000, 1_400_000)],
        )
        self.assertEqual(
            [event.copy_number for event in report.events],
            [1.0, 4.0],
        )
        for forbidden_parameter in (
            "snfj_input",
            "snv_input",
            "population_mode",
            "cancer_mode",
        ):
            self.assertFalse(report.tool.parameters[forbidden_parameter])
        retained = {artifact.relative_path for artifact in report.native_artifacts}
        self.assertIn("SYNTHETIC_DEPTH_ONLY.vcf.gz", retained)
        self.assertIn("SYNTHETIC_DEPTH_ONLY.vcf.gz.tbi", retained)
        self.assertIn("SYNTHETIC_DEPTH_ONLY.spc.gz", retained)
        self.assertTrue(all(artifact.fingerprint.sha256 for artifact in report.native_artifacts))

        empty_output = self.root / "empty-output"
        empty = run_spectre_depth_only(
            coverage_path=self._coverage("uniform", altered=False),
            sample_id="SYNTHETIC_EMPTY",
            output_dir=empty_output,
            reference_fasta=self.reference,
            genome_build=GenomeBuild.GRCH38,
            policy=self.policy,
            threads=1,
        )
        self.assertEqual(empty.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(empty.events, [])
        self.assertEqual(empty.raw_record_count, 0)
        self.assertTrue((empty_output / "SYNTHETIC_EMPTY.vcf").is_file())
        self.assertTrue(
            any("not a biological negative" in warning.lower() for warning in empty.warnings)
        )


if __name__ == "__main__":
    unittest.main()
