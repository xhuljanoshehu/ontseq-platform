from __future__ import annotations

import hashlib
import unittest

from ontseq_platform.tumor_inputs import (
    TumorAuxiliaryInputArtifact,
    TumorInputBundle,
)
from pydantic import ValidationError

from ontseq_platform.models import GenomeBuild
from ontseq_platform.multicaller_contracts import CallerInputRole


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _artifact(
    role: CallerInputRole,
    *,
    artifact_id: str | None = None,
    analysis_sample_id: str = "TUMOR_001",
    source_sample_id: str = "TUMOR_001",
    genome_build: GenomeBuild = GenomeBuild.GRCH38,
    reference_id: str = "GRCh38-test",
    reference_sha256: str = _sha("reference"),
    sha256: str | None = None,
    producer_lane_id: str | None = None,
    producer_lane_sha256: str | None = None,
) -> TumorAuxiliaryInputArtifact:
    return TumorAuxiliaryInputArtifact(
        artifact_id=artifact_id or f"artifact-{role.value}",
        role=role,
        analysis_sample_id=analysis_sample_id,
        source_sample_id=source_sample_id,
        genome_build=genome_build,
        reference_id=reference_id,
        reference_sha256=reference_sha256,
        sha256=sha256 or _sha(role.value),
        producer_lane_id=producer_lane_id,
        producer_lane_sha256=producer_lane_sha256,
    )


class TumorAuxiliaryInputContractTests(unittest.TestCase):
    def test_bundle_accepts_explicitly_bound_tumor_normal_and_phasing_inputs(self) -> None:
        bundle = TumorInputBundle(
            analysis_sample_id="TUMOR_001",
            genome_build=GenomeBuild.GRCH38,
            reference_id="GRCh38-test",
            reference_sha256=_sha("reference"),
            artifacts=[
                _artifact(CallerInputRole.TUMOR_BAM),
                _artifact(
                    CallerInputRole.MATCHED_NORMAL_BAM,
                    source_sample_id="NORMAL_001",
                ),
                _artifact(
                    CallerInputRole.PHASED_VARIANTS,
                    producer_lane_id="lane-phasing",
                    producer_lane_sha256=_sha("lane-phasing"),
                ),
            ],
        )

        self.assertEqual(bundle.analysis_sample_id, "TUMOR_001")
        self.assertEqual(
            {item.role for item in bundle.artifacts},
            {
                CallerInputRole.TUMOR_BAM,
                CallerInputRole.MATCHED_NORMAL_BAM,
                CallerInputRole.PHASED_VARIANTS,
            },
        )

    def test_bundle_rejects_analysis_sample_mismatch(self) -> None:
        with self.assertRaisesRegex(ValidationError, "analysis sample"):
            TumorInputBundle(
                analysis_sample_id="TUMOR_001",
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_sha256=_sha("reference"),
                artifacts=[
                    _artifact(
                        CallerInputRole.SNP_VCF,
                        analysis_sample_id="TUMOR_002",
                    )
                ],
            )

    def test_bundle_rejects_genome_build_mismatch(self) -> None:
        with self.assertRaisesRegex(ValidationError, "genome build"):
            TumorInputBundle(
                analysis_sample_id="TUMOR_001",
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_sha256=_sha("reference"),
                artifacts=[
                    _artifact(
                        CallerInputRole.BAF_TABLE,
                        genome_build=GenomeBuild.GRCH37,
                    )
                ],
            )

    def test_bundle_rejects_reference_identity_mismatch(self) -> None:
        for updates in (
            {"reference_id": "different-reference"},
            {"reference_sha256": _sha("different-reference")},
        ):
            with self.subTest(updates=updates), self.assertRaisesRegex(
                ValidationError, "reference"
            ):
                TumorInputBundle(
                    analysis_sample_id="TUMOR_001",
                    genome_build=GenomeBuild.GRCH38,
                    reference_id="GRCh38-test",
                    reference_sha256=_sha("reference"),
                    artifacts=[
                        _artifact(
                            CallerInputRole.PHASED_VARIANTS,
                            **updates,
                        )
                    ],
                )

    def test_artifact_requires_content_fingerprint(self) -> None:
        with self.assertRaises(ValidationError):
            TumorAuxiliaryInputArtifact(
                artifact_id="artifact-snp",
                role=CallerInputRole.SNP_VCF,
                analysis_sample_id="TUMOR_001",
                source_sample_id="TUMOR_001",
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_sha256=_sha("reference"),
                sha256=None,
            )

    def test_parent_lane_binding_requires_id_and_hash_together(self) -> None:
        cases = (
            {"producer_lane_id": "lane-parent"},
            {"producer_lane_sha256": _sha("lane-parent")},
        )
        for updates in cases:
            with self.subTest(updates=updates), self.assertRaises(ValidationError):
                _artifact(CallerInputRole.PHASED_VARIANTS, **updates)

    def test_bundle_rejects_duplicate_roles(self) -> None:
        with self.assertRaisesRegex(ValidationError, "roles"):
            TumorInputBundle(
                analysis_sample_id="TUMOR_001",
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_sha256=_sha("reference"),
                artifacts=[
                    _artifact(
                        CallerInputRole.SNP_VCF,
                        artifact_id="snp-one",
                    ),
                    _artifact(
                        CallerInputRole.SNP_VCF,
                        artifact_id="snp-two",
                    ),
                ],
            )

    def test_bundle_rejects_non_tumor_auxiliary_roles(self) -> None:
        with self.assertRaisesRegex(ValidationError, "not a tumor auxiliary input role"):
            _artifact(CallerInputRole.COVERAGE_PROFILE)


if __name__ == "__main__":
    unittest.main()
