from __future__ import annotations

import unittest

from pydantic import ValidationError

from ontseq_platform.models import (
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    PrivacySpec,
    SampleManifest,
)


class ManifestValidationTests(unittest.TestCase):
    def test_aligned_bam_requires_index(self) -> None:
        with self.assertRaises(ValidationError):
            InputSpec(kind=InputKind.ALIGNED_BAM, path="sample.bam")

    def test_adaptive_sampling_requires_versioned_bed(self) -> None:
        with self.assertRaises(ValidationError):
            AssaySpec(
                mode=AssayMode.ADAPTIVE_SAMPLING,
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38",
            )

    def test_direct_identifiers_are_blocked(self) -> None:
        with self.assertRaises(ValidationError):
            PrivacySpec(contains_direct_identifiers=True)

    def test_minimal_lcwgs_manifest(self) -> None:
        manifest = SampleManifest(
            sample_id="SYNTHETIC_001",
            run_id="RUN_001",
            input=InputSpec(kind=InputKind.POD5, path="/secure/run"),
            assay=AssaySpec(
                mode=AssayMode.LOW_COVERAGE_WGS,
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-v1",
            ),
            analysis=AnalysisSpec(profile="lcwgs", modules=[]),
        )
        self.assertEqual(manifest.assay.mode, AssayMode.LOW_COVERAGE_WGS)


if __name__ == "__main__":
    unittest.main()
