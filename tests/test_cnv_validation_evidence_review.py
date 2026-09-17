from __future__ import annotations

import hashlib
import unittest

from ontseq_platform.cnv_validation_contracts import CnvDataBasis, CnvRepeatKind
from ontseq_platform.cnv_validation_evidence import (
    CnvContributionStatus,
    CnvEvidenceRecordKind,
    CnvFullEvidenceRecord,
    CnvNativeArtifactReference,
    CnvNormalizedEventRecord,
    CnvRunOutcomeState,
    canonical_evidence_sha256,
    seal_cnv_validation_evidence,
)
from ontseq_platform.models import EventType, GenomeBuild, Locus


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _base_full_payload() -> dict[str, object]:
    parameters: dict[str, object] = {"bin_size_kbp": 500, "policy": "synthetic"}
    return {
        "record_id": "segment-primary",
        "registration_sha256": _sha("registration"),
        "specimen_id": "SYNTHETIC_CNV_REVIEW",
        "biological_specimen_id": "SYNTHETIC_BIO_REVIEW",
        "caller_id": "qdnaseq_ace",
        "caller_version": "synthetic-1",
        "adapter_policy_sha256": _sha("adapter-policy"),
        "execution_identity_sha256": _sha("runtime"),
        "genome_build": GenomeBuild.GRCH38,
        "data_basis": CnvDataBasis.LCWGS_GENOME_WIDE,
        "reference_id": "synthetic-grch38",
        "reference_sha256": _sha("reference"),
        "input_sha256": _sha("input"),
        "coverage_x": 5.0,
        "coverage_definition": "synthetic mean autosomal depth",
        "tumor_fraction": 0.2,
        "tumor_fraction_method": "synthetic orthogonal fraction",
        "tumor_fraction_timepoint": "same synthetic aliquot",
        "bin_size_kbp": 500,
        "repeat_kind": CnvRepeatKind.INDEPENDENT,
        "replicate_id": "replicate-1",
        "record_kind": CnvEvidenceRecordKind.CALLER_SEGMENT,
        "native_record_id": "segment-primary",
        "run_outcome": CnvRunOutcomeState.OBSERVED,
        "contribution_status": CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
        "caller_parameters": parameters,
        "caller_parameters_sha256": canonical_evidence_sha256(parameters),
        "dependency_versions": {"QDNAseq": "synthetic-1", "ACE": "synthetic-1"},
        "native_artifact_ids": ["caller-output"],
        "event_type": EventType.DELETION,
        "primary": Locus(chromosome="7", start=1_000_000, end=6_000_000),
        "copy_number": 1.0,
    }


def _full_record(**updates: object) -> CnvFullEvidenceRecord:
    payload = _base_full_payload()
    payload.update(updates)
    return CnvFullEvidenceRecord.model_validate(payload)


def _normalized(source_ids: list[str], **updates: object) -> CnvNormalizedEventRecord:
    payload: dict[str, object] = {
        "normalized_event_id": "normalized-del-7",
        "registration_sha256": _sha("registration"),
        "specimen_id": "SYNTHETIC_CNV_REVIEW",
        "biological_specimen_id": "SYNTHETIC_BIO_REVIEW",
        "caller_id": "qdnaseq_ace",
        "caller_version": "synthetic-1",
        "adapter_policy_sha256": _sha("adapter-policy"),
        "execution_identity_sha256": _sha("runtime"),
        "normalization_policy_sha256": _sha("normalization-policy"),
        "genome_build": GenomeBuild.GRCH38,
        "data_basis": CnvDataBasis.LCWGS_GENOME_WIDE,
        "reference_id": "synthetic-grch38",
        "reference_sha256": _sha("reference"),
        "input_sha256": _sha("input"),
        "coverage_x": 5.0,
        "coverage_definition": "synthetic mean autosomal depth",
        "tumor_fraction": 0.2,
        "tumor_fraction_method": "synthetic orthogonal fraction",
        "tumor_fraction_timepoint": "same synthetic aliquot",
        "bin_size_kbp": 500,
        "repeat_kind": CnvRepeatKind.INDEPENDENT,
        "replicate_id": "replicate-1",
        "source_full_evidence_ids": source_ids,
        "event_type": EventType.DELETION,
        "primary": Locus(chromosome="7", start=1_000_000, end=6_000_000),
        "normalized_copy_number": 1.0,
        "contribution_status": CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
    }
    payload.update(updates)
    return CnvNormalizedEventRecord.model_validate(payload)


def _artifact() -> CnvNativeArtifactReference:
    return CnvNativeArtifactReference(
        artifact_id="caller-output",
        role="caller_segments",
        relative_path="cnv/synthetic/caller-output.tsv",
        sha256=_sha("caller-output"),
        size_bytes=100,
        media_type="text/tab-separated-values",
    )


def _seal(full_evidence: list[CnvFullEvidenceRecord], normalized: CnvNormalizedEventRecord) -> None:
    seal_cnv_validation_evidence(
        manifest_id="SYNTHETIC_CNV_REVIEW_MANIFEST",
        registration_sha256=_sha("registration"),
        native_artifacts=[_artifact()],
        full_evidence=full_evidence,
        normalized_events=[normalized],
    )


class CnvNormalizedEventReviewTests(unittest.TestCase):
    def test_positive_normalized_event_requires_locus_bearing_caller_source(self) -> None:
        fit_only = _full_record(
            record_id="fit-only",
            native_record_id="fit-only",
            record_kind=CnvEvidenceRecordKind.CALLER_FIT,
            fit_group_id="fit-500",
            selected_fit=True,
            event_type=None,
            primary=None,
            copy_number=None,
        )
        with self.assertRaises(ValueError):
            _seal([fit_only], _normalized(["fit-only"]))

    def test_primary_normalized_event_cannot_be_promoted_from_excluded_only_source(self) -> None:
        excluded_segment = _full_record(
            record_id="segment-excluded",
            native_record_id="segment-excluded",
            contribution_status=CnvContributionStatus.EXCLUDED_FROM_PRIMARY_METRIC,
            contribution_reason="Synthetic preregistered exclusion.",
        )
        with self.assertRaises(ValueError):
            _seal([excluded_segment], _normalized(["segment-excluded"]))

    def test_primary_event_may_keep_secondary_fit_as_additional_source(self) -> None:
        primary_segment = _full_record()
        secondary_fit = _full_record(
            record_id="fit-secondary",
            native_record_id="fit-secondary",
            record_kind=CnvEvidenceRecordKind.CALLER_FIT,
            contribution_status=CnvContributionStatus.SECONDARY,
            fit_group_id="fit-500",
            selected_fit=True,
            event_type=None,
            primary=None,
            copy_number=None,
        )
        _seal(
            [primary_segment, secondary_fit],
            _normalized(["segment-primary", "fit-secondary"]),
        )


if __name__ == "__main__":
    unittest.main()
