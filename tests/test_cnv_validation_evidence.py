from __future__ import annotations

import hashlib
import math
import unittest

from pydantic import ValidationError

from ontseq_platform.cnv_validation_contracts import CnvDataBasis, CnvRepeatKind
from ontseq_platform.cnv_validation_evidence import (
    CnvContributionStatus,
    CnvEvidenceRecordKind,
    CnvFullEvidenceRecord,
    CnvNativeArtifactReference,
    CnvNumericMeasurement,
    CnvRunOutcomeState,
    canonical_evidence_sha256,
)
from ontseq_platform.models import GenomeBuild


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _parameters() -> dict[str, object]:
    return {
        "bin_sizes_kbp": [100, 500, 1000],
        "ace_penalty": 0.6,
        "normalization": {"coordinate_system": "zero_based_half_open"},
    }


def _full_record(**updates: object) -> CnvFullEvidenceRecord:
    parameters = _parameters()
    payload: dict[str, object] = {
        "record_id": "full-fit-500-selected",
        "registration_sha256": _sha("registration"),
        "specimen_id": "SYNTHETIC_CNV_001",
        "biological_specimen_id": "SYNTHETIC_BIO_001",
        "caller_id": "qdnaseq_ace",
        "caller_version": "synthetic-1",
        "adapter_policy_sha256": _sha("qdnaseq-policy"),
        "execution_identity_sha256": _sha("qdnaseq-runtime"),
        "genome_build": GenomeBuild.GRCH38,
        "data_basis": CnvDataBasis.LCWGS_GENOME_WIDE,
        "reference_id": "synthetic-grch38",
        "reference_sha256": _sha("reference"),
        "input_sha256": _sha("input"),
        "coverage_x": 5.25,
        "coverage_definition": "synthetic mean autosomal depth",
        "tumor_fraction": 0.21,
        "tumor_fraction_method": "synthetic orthogonal fraction",
        "tumor_fraction_timepoint": "same synthetic aliquot",
        "bin_size_kbp": 500,
        "repeat_kind": CnvRepeatKind.INDEPENDENT,
        "replicate_id": "replicate-1",
        "record_kind": CnvEvidenceRecordKind.CALLER_FIT,
        "native_record_id": "fit-500-selected",
        "run_outcome": CnvRunOutcomeState.OBSERVED,
        "contribution_status": CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
        "caller_parameters": parameters,
        "caller_parameters_sha256": canonical_evidence_sha256(parameters),
        "dependency_versions": {"QDNAseq": "1.42.0", "ACE": "1.24.0"},
        "native_artifact_ids": ["fit-model-500", "fit-plot-500"],
        "fit_group_id": "fit-500",
        "selected_fit": True,
        "cellularity": 0.23,
        "ploidy": 2.10,
        "fit_error": 0.012,
        "numeric_measurements": [
            CnvNumericMeasurement(
                name="segment_count",
                value=12.0,
                unit="count",
                source_field="segment_count",
            )
        ],
    }
    payload.update(updates)
    return CnvFullEvidenceRecord.model_validate(payload)


class CnvFullEvidenceContractTests(unittest.TestCase):
    def test_continuous_source_values_and_metadata_are_preserved(self) -> None:
        record = _full_record()
        self.assertEqual(record.coverage_x, 5.25)
        self.assertEqual(record.coverage_definition, "synthetic mean autosomal depth")
        self.assertEqual(record.tumor_fraction, 0.21)
        self.assertEqual(record.tumor_fraction_method, "synthetic orthogonal fraction")
        self.assertEqual(record.tumor_fraction_timepoint, "same synthetic aliquot")

    def test_multiresolution_records_remain_independently_addressable(self) -> None:
        records = [
            _full_record(
                record_id=f"full-fit-{size}-selected",
                native_record_id=f"fit-{size}-selected",
                fit_group_id=f"fit-{size}",
                bin_size_kbp=size,
            )
            for size in (100, 500, 1000)
        ]
        self.assertEqual([item.bin_size_kbp for item in records], [100, 500, 1000])
        self.assertEqual(len({item.address_key for item in records}), 3)

    def test_selected_and_alternative_fits_are_separate_retained_records(self) -> None:
        selected = _full_record()
        alternative = _full_record(
            record_id="full-fit-500-alt-1",
            native_record_id="fit-500-alt-1",
            selected_fit=False,
            cellularity=0.31,
            ploidy=2.35,
            fit_error=0.014,
            contribution_status=CnvContributionStatus.SECONDARY,
        )
        self.assertTrue(selected.selected_fit)
        self.assertFalse(alternative.selected_fit)
        self.assertEqual(selected.fit_group_id, alternative.fit_group_id)
        self.assertNotEqual(selected.record_id, alternative.record_id)

    def test_nonfinite_numeric_measurement_is_rejected(self) -> None:
        for value in (math.nan, math.inf, -math.inf):
            with self.subTest(value=value), self.assertRaises(ValidationError):
                CnvNumericMeasurement(
                    name="invalid",
                    value=value,
                    source_field="synthetic",
                )

    def test_nonfinite_caller_parameter_is_rejected_recursively(self) -> None:
        parameters = {"nested": {"bad": math.nan}}
        with self.assertRaises(ValidationError):
            _full_record(
                caller_parameters=parameters,
                caller_parameters_sha256=_sha("not-accepted"),
            )

    def test_caller_parameter_checksum_must_match_snapshot(self) -> None:
        with self.assertRaises(ValidationError):
            _full_record(caller_parameters_sha256=_sha("different-parameters"))

    def test_native_artifact_path_must_be_safe_and_relative(self) -> None:
        valid = CnvNativeArtifactReference(
            artifact_id="bins-500",
            role="caller_bins",
            relative_path="cnv/qdnaseq/500/bins.tsv",
            sha256=_sha("bins"),
            size_bytes=123,
            media_type="text/tab-separated-values",
        )
        self.assertEqual(valid.relative_path, "cnv/qdnaseq/500/bins.tsv")

        for path in ("../secret.tsv", "/tmp/absolute.tsv", "C:\\temp\\bins.tsv", "a\\b.tsv"):
            with self.subTest(path=path), self.assertRaises(ValidationError):
                CnvNativeArtifactReference(
                    artifact_id="unsafe",
                    role="caller_bins",
                    relative_path=path,
                    sha256=_sha("unsafe"),
                    size_bytes=1,
                    media_type="text/tab-separated-values",
                )

    def test_failed_no_call_not_assessable_and_observed_are_distinct(self) -> None:
        observed = _full_record()
        failed = _full_record(
            record_id="full-run-failed",
            native_record_id="run-failed",
            record_kind=CnvEvidenceRecordKind.RUN_SUMMARY,
            run_outcome=CnvRunOutcomeState.FAILED,
            contribution_status=CnvContributionStatus.FAILED,
            outcome_reason="Synthetic execution failure.",
            native_artifact_ids=[],
            fit_group_id=None,
            selected_fit=None,
            cellularity=None,
            ploidy=None,
            fit_error=None,
        )
        no_call = _full_record(
            record_id="full-run-no-call",
            native_record_id="run-no-call",
            record_kind=CnvEvidenceRecordKind.RUN_SUMMARY,
            run_outcome=CnvRunOutcomeState.NO_CALL,
            contribution_status=CnvContributionStatus.NO_CALL,
            outcome_reason="Synthetic completed run with no callable solution.",
            native_artifact_ids=[],
            fit_group_id=None,
            selected_fit=None,
            cellularity=None,
            ploidy=None,
            fit_error=None,
        )
        not_assessable = _full_record(
            record_id="full-run-not-assessable",
            native_record_id="run-not-assessable",
            record_kind=CnvEvidenceRecordKind.RUN_SUMMARY,
            run_outcome=CnvRunOutcomeState.NOT_ASSESSABLE,
            contribution_status=CnvContributionStatus.NOT_ASSESSABLE,
            outcome_reason="Synthetic assessability mask excludes the requested territory.",
            native_artifact_ids=[],
            fit_group_id=None,
            selected_fit=None,
            cellularity=None,
            ploidy=None,
            fit_error=None,
        )
        self.assertEqual(
            {
                observed.run_outcome,
                failed.run_outcome,
                no_call.run_outcome,
                not_assessable.run_outcome,
            },
            {
                CnvRunOutcomeState.OBSERVED,
                CnvRunOutcomeState.FAILED,
                CnvRunOutcomeState.NO_CALL,
                CnvRunOutcomeState.NOT_ASSESSABLE,
            },
        )

    def test_biological_negative_requires_orthogonal_truth_linkage(self) -> None:
        common = {
            "record_id": "full-run-biological-negative",
            "native_record_id": "run-biological-negative",
            "record_kind": CnvEvidenceRecordKind.RUN_SUMMARY,
            "run_outcome": CnvRunOutcomeState.BIOLOGICAL_NEGATIVE_ORTHOGONAL_TRUTH,
            "contribution_status": CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
            "outcome_reason": "Synthetic orthogonal truth establishes an assessable negative.",
            "native_artifact_ids": [],
            "fit_group_id": None,
            "selected_fit": None,
            "cellularity": None,
            "ploidy": None,
            "fit_error": None,
        }
        with self.assertRaises(ValidationError):
            _full_record(**common, truth_support_resource_ids=[])

        biological_negative = _full_record(
            **common,
            truth_support_resource_ids=["synthetic-truth-negative-v1"],
        )
        self.assertEqual(
            biological_negative.run_outcome,
            CnvRunOutcomeState.BIOLOGICAL_NEGATIVE_ORTHOGONAL_TRUTH,
        )

    def test_excluded_record_is_retained_as_full_evidence(self) -> None:
        excluded = _full_record(
            record_id="full-fit-500-excluded",
            native_record_id="fit-500-excluded",
            selected_fit=False,
            contribution_status=CnvContributionStatus.EXCLUDED_FROM_PRIMARY_METRIC,
            contribution_reason="Synthetic preregistered exclusion classification.",
        )
        self.assertEqual(
            excluded.contribution_status,
            CnvContributionStatus.EXCLUDED_FROM_PRIMARY_METRIC,
        )
        self.assertEqual(excluded.run_outcome, CnvRunOutcomeState.OBSERVED)


if __name__ == "__main__":
    unittest.main()
