from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path
from unittest.mock import patch

from ontseq_platform.cutesv import run_cutesv
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CuteSvPolicy,
    GenomeBuild,
    InputKind,
    InputSpec,
    SampleManifest,
    Verdict,
)

VALID_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
chr1\t1000\tcute-1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1200;SVLEN=-201;RE=8\tGT:DR:DV\t0/1:8:8
"""


class LocalErrorRunner:
    """Model cuteSV returning zero plus a complete partial VCF after a swallowed worker error."""

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(str(item) for item in argv)
        if "--version" in normalized:
            return CommandResult(normalized, 0, "cuteSV 2.1.3\n", "")
        Path(normalized[3]).write_text(VALID_VCF, encoding="utf-8")
        return CommandResult(normalized, 0, "", "[INFO] LocalError: {e}\n")


def _inputs(root: Path) -> tuple[SampleManifest, AlignedBamIntakeReport, Path]:
    bam = root / "synthetic.bam"
    bai = root / "synthetic.bam.bai"
    reference = root / "synthetic.fasta"
    bam.write_bytes(b"synthetic")
    bai.write_bytes(b"synthetic")
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    manifest = SampleManifest(
        sample_id="SYNTHETIC_CUTESV_LOCALERROR",
        run_id="SYNTHETIC_CUTESV_LOCALERROR_RUN",
        input=InputSpec(kind=InputKind.ALIGNED_BAM, path=str(bam), index_path=str(bai)),
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
    intake = AlignedBamIntakeReport(
        sample_id=manifest.sample_id,
        reference_id=manifest.assay.reference_id,
        genome_build=manifest.assay.genome_build,
        checks=[],
        verdict=Verdict.PASS,
    )
    return manifest, intake, reference


def _policy() -> CuteSvPolicy:
    return CuteSvPolicy(
        profile_id="synthetic-cutesv-localerror",
        status="technical_defaults_only",
        note="Synthetic failure-propagation regression only.",
    )


class CuteSvLocalErrorGuardTests(unittest.TestCase):
    def test_zero_exit_with_upstream_localerror_marker_fails_closed(self) -> None:
        """A swallowed clustering-worker error must not promote a partial VCF."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, intake, reference = _inputs(root)
            output = root / "calls.vcf"
            scratch = root / "scratch"
            scratch.mkdir()

            with (
                patch.dict("os.environ", {"ONTSEQ_CUTESV_SCRATCH_ROOT": str(scratch)}),
                self.assertRaisesRegex(ValueError, "LocalError|worker"),
            ):
                run_cutesv(
                    manifest,
                    intake,
                    _policy(),
                    reference_fasta=reference,
                    output_vcf=output,
                    runner=LocalErrorRunner(),
                )

            self.assertFalse(output.exists())
            self.assertEqual(list(scratch.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
