from __future__ import annotations

import unittest

from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import (
    ISCNProposalStatus,
    LegacyResourceContext,
    PipelineResult,
)


class PipelineResultV020Tests(unittest.TestCase):
    def test_pipeline_result_010_and_020_remain_readable_as_legacy_unspecified(self) -> None:
        for schema_version in ("0.1.0", "0.2.0"):
            with self.subTest(schema_version=schema_version):
                payload = build_demo_result().model_dump(mode="json")
                payload["schema_version"] = schema_version
                payload.pop("reference_context", None)
                payload.pop("sidecars", None)
                payload["iscn"] = {
                    "notation": "46,XX,+8",
                    "standard_edition": "ISCN 2024",
                    "conformance_profile": "subset-v0.1-unvalidated",
                    "review_status": "REVIEW_REQUIRED",
                    "source_event_ids": ["CNV-001"],
                    "warnings": ["Legacy unvalidated proposal."],
                }
                restored = PipelineResult.model_validate(payload)
                self.assertIsInstance(restored.reference_context, LegacyResourceContext)
                self.assertEqual(restored.reference_context.status, "legacy_unspecified")
                self.assertEqual(
                    restored.iscn.proposal_status,
                    ISCNProposalStatus.LEGACY_UNSPECIFIED,
                )

    def test_pipeline_result_030_refuses_legacy_iscn_semantics(self) -> None:
        payload = build_demo_result().model_dump(mode="json")
        payload["iscn"] = {
            "notation": "46,XX",
            "conformance_profile": "legacy-unspecified",
        }

        with self.assertRaises(ValueError):
            PipelineResult.model_validate(payload)

    def test_pipeline_result_030_binds_the_proposal_to_the_manifest_build(self) -> None:
        payload = build_demo_result().model_dump(mode="json")
        payload["iscn"]["genome_build"] = "GRCh37"
        payload["iscn"]["resource_provenance"]["genome_build"] = "GRCh37"

        with self.assertRaisesRegex(ValueError, "genome build must match the run manifest"):
            PipelineResult.model_validate(payload)

    def test_pipeline_result_030_rejects_provenance_mismatches(self) -> None:
        mutations = {
            "reference dictionary contract": (
                "reference_dictionary_contract",
                "grch38_canonical_25",
            ),
            "reference bundle ID": ("reference_bundle_id", "SYNTHETIC_OTHER_BUNDLE"),
            "reference bundle version": ("reference_bundle_version", "other-v2"),
            "reference lock SHA256": ("reference_lock_sha256", "3" * 64),
            "annotation cache SHA256": ("annotation_cache_sha256", "4" * 64),
            "cytoband SHA256": ("cytoband_sha256", "5" * 64),
            "cytoband release": ("cytoband_release", "other-release"),
        }
        for expected_label, (field, replacement) in mutations.items():
            with self.subTest(field=field):
                payload = build_demo_result().model_dump(mode="json")
                payload["iscn"]["resource_provenance"][field] = replacement

                with self.assertRaisesRegex(
                    ValueError,
                    rf"provenance does not match.*{expected_label}",
                ):
                    PipelineResult.model_validate(payload)

    def test_pipeline_result_030_binds_provenance_to_reference_context_build(self) -> None:
        payload = build_demo_result().model_dump(mode="json")
        payload["reference_context"]["genome_build"] = "GRCh37"

        with self.assertRaisesRegex(ValueError, "provenance does not match.*genome build"):
            PipelineResult.model_validate(payload)

    def test_partial_iscn_proposal_requires_a_completed_cnv_module(self) -> None:
        for cnv_status in ("NOT_RUN", "FAILED"):
            with self.subTest(cnv_status=cnv_status):
                payload = build_demo_result().model_dump(mode="json")
                cnv_module = next(item for item in payload["modules"] if item["module"] == "cnv")
                cnv_module["status"] = cnv_status

                with self.assertRaises(ValueError):
                    PipelineResult.model_validate(payload)


if __name__ == "__main__":
    unittest.main()
