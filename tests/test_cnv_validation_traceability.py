from __future__ import annotations

import hashlib
import math
import unittest

from pydantic import ValidationError

from ontseq_platform.cnv_validation_contracts import CnvDataBasis, CnvRepeatKind
from ontseq_platform.cnv_validation_evidence import (
    CnvAggregateEvidenceTrace,
    CnvContributionStatus,
    CnvEvidenceRecordKind,
    CnvFullEvidenceRecord,
    CnvNativeArtifactReference,
    CnvNormalizedEventRecord,
    CnvRunOutcomeState,
    CnvTraceabilityIndex,
    CnvValidationEvidenceManifest,
    canonical_evidence_sha256,
    seal_cnv_traceability_index,
    seal_cnv_validation_evidence,
    verify_cnv_traceability,
)
from ontseq_platform.models import EventType, GenomeBuild, Locus


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _source_record() -> CnvFullEvidenceRecord:
    parameters: dict[str, object] = {"bin_size_kbp": 500, "policy": "synthetic"}
    return CnvFullEvidenceRecord(
        record_id="full-segment-500-1",
        registration_sha256=_sha("registration"),
        specimen_id="SYNTHETIC_CNV_001",
        biological_specimen_id="SYNTHETIC_BIO_001",
        caller_id="qdnaseq_ace",
        caller_version="synthetic-1",
        adapter_policy_sha256=_sha("adapter-policy"),
        execution_identity_sha256=_sha("runtime"),
        genome_build=GenomeBuild.GRCH38,
        data_basis=CnvDataBasis.LCWGS_GENOME_WIDE,
        reference_id="synthetic-grch38",
        reference_sha256=_sha("reference"),
        input_sha256=_sha("input"),
        coverage_x=5.0,
        coverage_definition="synthetic mean autosomal depth",
        tumor_fraction=0.2,
        tumor_fraction_method="synthetic orthogonal fraction",
        tumor_fraction_timepoint="same synthetic aliquot",
        bin_size_kbp=500,
        repeat_kind=CnvRepeatKind.INDEPENDENT,
        replicate_id="replicate-1",
        record_kind=CnvEvidenceRecordKind.CALLER_SEGMENT,
        native_record_id="segment-1",
        run_outcome=CnvRunOutcomeState.OBSERVED,
        contribution_status=CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
        caller_parameters=parameters,
        caller_parameters_sha256=canonical_evidence_sha256(parameters),
        dependency_versions={"QDNAseq": "synthetic-1", "ACE": "synthetic-1"},
        native_artifact_ids=["segments-500"],
        event_type=EventType.DELETION,
        primary=Locus(chromosome="7", start=1_000_000, end=6_000_000),
        copy_number=1.0,
        log2_ratio=-0.5,
    )


def _artifact() -> CnvNativeArtifactReference:
    return CnvNativeArtifactReference(
        artifact_id="segments-500",
        role="caller_segments",
        relative_path="cnv/qdnaseq/500/segments.tsv",
        sha256=_sha("segments-500"),
        size_bytes=123,
        media_type="text/tab-separated-values",
    )


def _normalized_event() -> CnvNormalizedEventRecord:
    return CnvNormalizedEventRecord(
        normalized_event_id="normalized-del-7-1",
        registration_sha256=_sha("registration"),
        specimen_id="SYNTHETIC_CNV_001",
        biological_specimen_id="SYNTHETIC_BIO_001",
        caller_id="qdnaseq_ace",
        caller_version="synthetic-1",
        adapter_policy_sha256=_sha("adapter-policy"),
        execution_identity_sha256=_sha("runtime"),
        normalization_policy_sha256=_sha("normalization-policy"),
        genome_build=GenomeBuild.GRCH38,
        data_basis=CnvDataBasis.LCWGS_GENOME_WIDE,
        reference_id="synthetic-grch38",
        reference_sha256=_sha("reference"),
        input_sha256=_sha("input"),
        coverage_x=5.0,
        coverage_definition="synthetic mean autosomal depth",
        tumor_fraction=0.2,
        tumor_fraction_method="synthetic orthogonal fraction",
        tumor_fraction_timepoint="same synthetic aliquot",
        bin_size_kbp=500,
        repeat_kind=CnvRepeatKind.INDEPENDENT,
        replicate_id="replicate-1",
        source_full_evidence_ids=["full-segment-500-1"],
        event_type=EventType.DELETION,
        primary=Locus(chromosome="7", start=1_000_000, end=6_000_000),
        normalized_copy_number=1.0,
        contribution_status=CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
    )


def _manifest() -> CnvValidationEvidenceManifest:
    return seal_cnv_validation_evidence(
        manifest_id="SYNTHETIC_CNV_EVIDENCE_001",
        registration_sha256=_sha("registration"),
        native_artifacts=[_artifact()],
        full_evidence=[_source_record()],
        normalized_events=[_normalized_event()],
    )


def _trace(**updates: object) -> CnvAggregateEvidenceTrace:
    payload: dict[str, object] = {
        "trace_id": "trace-primary-del-7",
        "aggregate_id": "aggregate-primary-del-7",
        "question_id": "primary_sensitivity",
        "stratum_key": {
            "caller": "qdnaseq_ace",
            "build": "GRCh38",
            "coverage_band": "5x",
            "tumor_fraction_band": "0.2",
            "event_class": "deletion",
            "event_size_band": "5mb",
            "bin_size_kbp": 500,
        },
        "numerator_full_evidence_ids": ["full-segment-500-1"],
        "denominator_full_evidence_ids": ["full-segment-500-1"],
        "excluded_full_evidence_ids": [],
        "normalized_event_ids": ["normalized-del-7-1"],
    }
    payload.update(updates)
    return CnvAggregateEvidenceTrace.model_validate(payload)


def _index(
    evidence: CnvValidationEvidenceManifest,
    *traces: CnvAggregateEvidenceTrace,
) -> CnvTraceabilityIndex:
    selected = list(traces) if traces else [_trace()]
    return seal_cnv_traceability_index(
        index_id="SYNTHETIC_CNV_TRACEABILITY_001",
        evidence_manifest_sha256=evidence.manifest_sha256,
        traces=selected,
    )


class CnvAggregateEvidenceTraceContractTests(unittest.TestCase):
    def test_trace_carries_question_aggregate_and_stratum_identity(self) -> None:
        trace = _trace()
        self.assertEqual(trace.aggregate_id, "aggregate-primary-del-7")
        self.assertEqual(trace.question_id, "primary_sensitivity")
        self.assertEqual(trace.stratum_key["caller"], "qdnaseq_ace")
        self.assertEqual(trace.stratum_key["bin_size_kbp"], 500)

    def test_trace_requires_at_least_one_evidence_reference(self) -> None:
        with self.assertRaises(ValidationError):
            _trace(
                numerator_full_evidence_ids=[],
                denominator_full_evidence_ids=[],
                excluded_full_evidence_ids=[],
                normalized_event_ids=[],
            )

    def test_trace_rejects_duplicate_ids_within_each_membership_list(self) -> None:
        for field in (
            "numerator_full_evidence_ids",
            "denominator_full_evidence_ids",
            "excluded_full_evidence_ids",
            "normalized_event_ids",
        ):
            with self.subTest(field=field), self.assertRaises(ValidationError):
                _trace(**{field: ["duplicate", "duplicate"]})

    def test_trace_keeps_numerator_denominator_and_excluded_memberships_distinct(self) -> None:
        with self.assertRaises(ValidationError):
            _trace(excluded_full_evidence_ids=["full-segment-500-1"])
        with self.assertRaises(ValidationError):
            _trace(
                numerator_full_evidence_ids=[],
                denominator_full_evidence_ids=["full-segment-500-1"],
                excluded_full_evidence_ids=["full-segment-500-1"],
            )

    def test_numerator_must_be_a_subset_of_denominator(self) -> None:
        trace = _trace()
        self.assertEqual(trace.numerator_full_evidence_ids, trace.denominator_full_evidence_ids)
        with self.assertRaises(ValidationError):
            _trace(
                numerator_full_evidence_ids=["full-segment-500-1"],
                denominator_full_evidence_ids=["different-denominator-record"],
            )

    def test_stratum_key_rejects_nonfinite_values(self) -> None:
        with self.assertRaises(ValidationError):
            _trace(stratum_key={"coverage": math.nan})


class CnvTraceabilityIndexContractTests(unittest.TestCase):
    def test_verifier_rejects_wrong_evidence_manifest_sha(self) -> None:
        evidence = _manifest()
        index = seal_cnv_traceability_index(
            index_id="SYNTHETIC_CNV_TRACEABILITY_WRONG_MANIFEST",
            evidence_manifest_sha256=_sha("different-manifest"),
            traces=[_trace()],
        )
        with self.assertRaises(ValueError):
            verify_cnv_traceability(evidence, index)

    def test_verifier_rejects_dangling_full_evidence_id(self) -> None:
        evidence = _manifest()
        index = _index(
            evidence,
            _trace(
                numerator_full_evidence_ids=["missing-full-evidence"],
                denominator_full_evidence_ids=["missing-full-evidence"],
            ),
        )
        with self.assertRaises(ValueError):
            verify_cnv_traceability(evidence, index)

    def test_verifier_rejects_dangling_normalized_event_id(self) -> None:
        evidence = _manifest()
        index = _index(evidence, _trace(normalized_event_ids=["missing-normalized-event"]))
        with self.assertRaises(ValueError):
            verify_cnv_traceability(evidence, index)

    def test_verifier_requires_normalized_sources_in_trace_membership(self) -> None:
        evidence = _manifest()
        index = _index(
            evidence,
            _trace(
                numerator_full_evidence_ids=[],
                denominator_full_evidence_ids=[],
                excluded_full_evidence_ids=[],
            ),
        )
        with self.assertRaises(ValueError):
            verify_cnv_traceability(evidence, index)

    def test_index_rejects_duplicate_trace_ids(self) -> None:
        evidence = _manifest()
        with self.assertRaises(ValidationError):
            _index(evidence, _trace(), _trace(aggregate_id="aggregate-other"))

    def test_index_rejects_duplicate_analytical_trace_addresses(self) -> None:
        evidence = _manifest()
        with self.assertRaises(ValidationError):
            _index(evidence, _trace(), _trace(trace_id="trace-second-copy"))

    def test_index_is_tamper_evident(self) -> None:
        evidence = _manifest()
        index = _index(evidence)
        payload = index.model_dump(mode="json")
        payload["traces"][0]["aggregate_id"] = "tampered-aggregate"
        with self.assertRaises(ValidationError):
            CnvTraceabilityIndex.model_validate(payload)

    def test_trace_walks_aggregate_to_normalized_to_full_to_native_artifact(self) -> None:
        evidence = _manifest()
        index = _index(evidence)
        verify_cnv_traceability(evidence, index)

        trace = index.traces[0]
        normalized_by_id = {item.normalized_event_id: item for item in evidence.normalized_events}
        full_by_id = {item.record_id: item for item in evidence.full_evidence}
        artifact_by_id = {item.artifact_id: item for item in evidence.native_artifacts}

        normalized = normalized_by_id[trace.normalized_event_ids[0]]
        source = full_by_id[normalized.source_full_evidence_ids[0]]
        artifact = artifact_by_id[source.native_artifact_ids[0]]

        self.assertEqual(normalized.normalized_event_id, "normalized-del-7-1")
        self.assertEqual(source.record_id, "full-segment-500-1")
        self.assertEqual(artifact.artifact_id, "segments-500")


if __name__ == "__main__":
    unittest.main()
