"""Pinned ichorCNA interoperability using generated cfDNA-like depth only."""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cnv.ichorcna import (
    IchorCnaResourceArtifact,
    IchorCnaResourceBundle,
    IchorCnaResourceRole,
    IchorCnaUlpWgsPolicy,
    run_ichorcna_ulp_wgs,
)
from ontseq_platform.models import GenomeBuild, ModuleRunStatus
from ontseq_platform.multicaller_contracts import CallerAssayRegime

_FIXED_STEP_CHROM = re.compile(r"(?:^|\s)chrom=([^\s]+)")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@unittest.skipUnless(
    os.environ.get("ONTSEQ_ICHORCNA_REAL_TOOL") == "1",
    "opt-in pinned ichorCNA 0.5.1 interoperability lane",
)
class IchorCnaBinaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.rscript = "Rscript"

        self.ichor_script = (
            Path(__file__).resolve().parents[1] / "scripts" / "run_ichorcna_0_5_1.R"
        )
        self.assertTrue(self.ichor_script.is_file())

        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

        self.gc_wig = self._system_file("gc_hg38_1000kb.wig")
        self.map_wig = self._system_file("map_hg38_1000kb.wig")
        self.centromere = self._system_file("GRCh38.GCA_000001405.2_centromere_acen.txt")
        self.coverage_wig = self.root / "synthetic.1mb.wig"
        self._write_synthetic_coverage()

    def _system_file(self, name: str) -> Path:
        expression = f'cat(system.file("extdata", "{name}", package="ichorCNA"))'
        completed = subprocess.run(
            [self.rscript, "-e", expression],
            check=True,
            capture_output=True,
            encoding="utf-8",
        )
        path = Path(completed.stdout.strip())
        self.assertTrue(path.is_file(), f"Missing ichorCNA package resource: {name}")
        return path

    def _write_synthetic_coverage(self) -> None:
        current_chromosome = ""
        bin_index = 0
        with (
            self.gc_wig.open("r", encoding="utf-8") as source,
            self.coverage_wig.open("w", encoding="utf-8") as target,
        ):
            for raw_line in source:
                line = raw_line.rstrip("\n")
                if line.startswith("fixedStep"):
                    match = _FIXED_STEP_CHROM.search(line)
                    self.assertIsNotNone(match)
                    assert match is not None
                    current_chromosome = match.group(1)
                    bin_index = 0
                    target.write(f"{line}\n")
                    continue
                if not line or line.startswith(("#", "track", "browser")):
                    target.write(f"{line}\n")
                    continue

                depth = 100 + ((bin_index % 7) - 3) * 2
                if current_chromosome == "chr7" and 40 <= bin_index < 65:
                    depth = 58
                elif current_chromosome == "chr7" and 80 <= bin_index < 105:
                    depth = 152
                target.write(f"{depth}\n")
                bin_index += 1

    def _artifact(
        self,
        *,
        resource_id: str,
        role: IchorCnaResourceRole,
        path: Path,
        bin_size_bp: int | None,
    ) -> IchorCnaResourceArtifact:
        return IchorCnaResourceArtifact(
            resource_id=resource_id,
            role=role,
            genome_build=GenomeBuild.GRCH38,
            bin_size_bp=bin_size_bp,
            size_bytes=path.stat().st_size,
            sha256=_sha256(path),
        )

    def _resources(self) -> IchorCnaResourceBundle:
        return IchorCnaResourceBundle(
            sample_id="SYNTHETIC_CFDNA",
            assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
            genome_build=GenomeBuild.GRCH38,
            bin_size_bp=1_000_000,
            resources=[
                self._artifact(
                    resource_id="coverage",
                    role=IchorCnaResourceRole.COVERAGE_WIG,
                    path=self.coverage_wig,
                    bin_size_bp=1_000_000,
                ),
                self._artifact(
                    resource_id="gc",
                    role=IchorCnaResourceRole.GC_WIG,
                    path=self.gc_wig,
                    bin_size_bp=1_000_000,
                ),
                self._artifact(
                    resource_id="mappability",
                    role=IchorCnaResourceRole.MAPPABILITY_WIG,
                    path=self.map_wig,
                    bin_size_bp=1_000_000,
                ),
                self._artifact(
                    resource_id="centromere",
                    role=IchorCnaResourceRole.CENTROMERE,
                    path=self.centromere,
                    bin_size_bp=None,
                ),
            ],
        )

    def _policy(self) -> IchorCnaUlpWgsPolicy:
        return IchorCnaUlpWgsPolicy(
            profile_id="ichorcna-0.5.1-ulp-wgs-real-tool-smoke",
            expected_version="0.5.1",
            assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
            genome_build=GenomeBuild.GRCH38,
            genome_build_label="hg38",
            genome_style="UCSC",
            bin_size_bp=1_000_000,
            script_sha256=_sha256(self.ichor_script),
            ploidy_candidates=[2, 3],
            normal_fraction_candidates=[0.5, 0.6, 0.7, 0.8, 0.9, 0.95],
            max_copy_number=5,
            minimum_mappability=0.75,
            estimate_normal=True,
            estimate_ploidy=True,
            estimate_subclone_prevalence=True,
            include_homozygous_deletion=False,
            transaction_e=0.9999,
            transaction_strength=10_000,
            timeout_seconds=900,
            real_tool_qualified=True,
            note="Synthetic ichorCNA 0.5.1 ULP-WGS interoperability smoke.",
        )

    def test_real_ichorcna_0_5_1_runs_through_adapter(self) -> None:
        output_dir = self.root / "ichor-output"
        report = run_ichorcna_ulp_wgs(
            coverage_wig=self.coverage_wig,
            gc_wig=self.gc_wig,
            mappability_wig=self.map_wig,
            centromere=self.centromere,
            ichorcna_script=self.ichor_script,
            sample_id="SYNTHETIC_CFDNA",
            output_dir=output_dir,
            resources=self._resources(),
            policy=self._policy(),
            assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
            rscript=self.rscript,
        )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.tool.version, "0.5.1")
        self.assertTrue(report.policy.real_tool_qualified)
        self.assertGreaterEqual(report.tumor_fraction, 0)
        self.assertGreater(report.ploidy, 0)
        self.assertTrue(report.segments)
        self.assertTrue(all(segment.reportable is False for segment in report.segments))
        roles = {artifact.role for artifact in report.native_artifacts}
        self.assertIn("parameters", roles)
        self.assertIn("segments", roles)
        self.assertIn("r_workspace", roles)
        self.assertTrue(all(artifact.fingerprint.sha256 for artifact in report.native_artifacts))


if __name__ == "__main__":
    unittest.main()
