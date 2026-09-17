from __future__ import annotations

import hashlib
import unittest

from pydantic import ValidationError

from ontseq_platform.cnv_validation_contracts import (
    CnvAcceptanceMetric,
    CnvAcceptanceQuestion,
    CnvCallerLock,
    CnvStratificationDimension,
    CnvStratificationPlan,
    CnvValidationMatrix,
)
from ontseq_platform.models import BenchmarkThresholds, EventType, GenomeBuild


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _matrix() -> CnvValidationMatrix:
    return CnvValidationMatrix(
        matrix_id="synthetic-cnv-validation-v1",
        callers=[
            CnvCallerLock(
                caller_id="qdnaseq_ace",
                caller_version="synthetic-1",
                adapter_policy_sha256=_sha("qdnaseq-policy"),
                execution_identity_sha256=_sha("qdnaseq-runtime"),
            ),
            CnvCallerLock(
                caller_id="spectre",
                caller_version="synthetic-1",
                adapter_policy_sha256=_sha("spectre-policy"),
                execution_identity_sha256=_sha("spectre-runtime"),
            ),
        ],
        genome_builds=[GenomeBuild.GRCH37, GenomeBuild.GRCH38],
        stratification=CnvStratificationPlan(
            primary_dimensions=[
                CnvStratificationDimension.COVERAGE,
                CnvStratificationDimension.TUMOR_FRACTION,
                CnvStratificationDimension.EVENT_CLASS,
                CnvStratificationDimension.EVENT_SIZE,
            ],
            secondary_dimensions=[
                CnvStratificationDimension.CALLER,
                CnvStratificationDimension.GENOME_BUILD,
                CnvStratificationDimension.BIN_SIZE,
            ],
            coverage_cutpoints_x=[2.0, 5.0, 10.0],
            tumor_fraction_cutpoints=[0.1, 0.2, 0.5],
            event_size_cutpoints_bp=[1_000_000, 5_000_000, 10_000_000],
        ),
        event_classes=[
            EventType.CHROMOSOME_GAIN,
            EventType.CHROMOSOME_LOSS,
            EventType.DELETION,
            EventType.DUPLICATION,
        ],
        qdnaseq_bin_sizes_kbp=[100, 500, 1000],
        matching_thresholds=BenchmarkThresholds(minimum_reciprocal_overlap=0.5),
        acceptance=[
            CnvAcceptanceQuestion(
                question_id="primary-sensitivity",
                metric=CnvAcceptanceMetric.SENSITIVITY,
                minimum_evaluable_denominator=5,
                minimum_acceptable=0.9,
            )
        ],
        rationale="Synthetic software contract only; no biological performance claim.",
    )


class CnvValidationContractTests(unittest.TestCase):
    def test_matrix_keeps_dimensions_and_cutpoints_explicit(self) -> None:
        matrix = _matrix()
        self.assertEqual(matrix.stratification.coverage_cutpoints_x, [2.0, 5.0, 10.0])
        self.assertEqual(matrix.stratification.tumor_fraction_cutpoints, [0.1, 0.2, 0.5])
        self.assertEqual(matrix.qdnaseq_bin_sizes_kbp, [100, 500, 1000])
        self.assertTrue(matrix.retain_all_evidence)

    def test_cutpoints_must_be_unique_and_strictly_sorted(self) -> None:
        payload = _matrix().model_dump()
        payload["stratification"]["coverage_cutpoints_x"] = [5.0, 2.0, 5.0]
        with self.assertRaises(ValidationError):
            CnvValidationMatrix.model_validate(payload)

    def test_primary_secondary_exploratory_dimensions_cannot_overlap(self) -> None:
        payload = _matrix().model_dump()
        payload["stratification"]["exploratory_dimensions"] = ["coverage"]
        with self.assertRaises(ValidationError):
            CnvValidationMatrix.model_validate(payload)

    def test_tumor_fraction_cutpoints_stay_inside_probability_domain(self) -> None:
        payload = _matrix().model_dump()
        payload["stratification"]["tumor_fraction_cutpoints"] = [-0.01, 0.5]
        with self.assertRaises(ValidationError):
            CnvValidationMatrix.model_validate(payload)

    def test_event_classes_reject_non_cnv_event_types(self) -> None:
        payload = _matrix().model_dump()
        payload["event_classes"].append(EventType.INVERSION)
        with self.assertRaises(ValidationError):
            CnvValidationMatrix.model_validate(payload)

    def test_acceptance_question_requires_a_bound(self) -> None:
        with self.assertRaises(ValidationError):
            CnvAcceptanceQuestion(
                question_id="invalid-unbounded-question",
                metric=CnvAcceptanceMetric.SENSITIVITY,
                minimum_evaluable_denominator=5,
            )


if __name__ == "__main__":
    unittest.main()
