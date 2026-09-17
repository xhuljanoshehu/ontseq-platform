from __future__ import annotations

import hashlib
import unittest
from datetime import UTC, datetime

from pydantic import ValidationError

from ontseq_platform.cnv_validation_contracts import (
    CnvAcceptanceMetric,
    CnvAcceptanceQuestion,
    CnvAssessabilityMask,
    CnvCallerLock,
    CnvDataBasis,
    CnvNegativeUniverse,
    CnvRepeatKind,
    CnvStratificationDimension,
    CnvStratificationPlan,
    CnvTruthSource,
    CnvValidationCohort,
    CnvValidationMatrix,
    CnvValidationSpecimen,
)
from ontseq_platform.cnv_validation_registration import (
    CnvValidationRegistration,
    canonical_content_sha256,
    cnv_cohort_eligibility_reasons,
    preregister_cnv_validation,
)
from ontseq_platform.models import (
    BenchmarkThresholds,
    EventType,
    GenomeBuild,
    GenomicEvent,
    Locus,
)


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
            )
        ],
        genome_builds=[GenomeBuild.GRCH38],
        data_bases=[CnvDataBasis.LCWGS_GENOME_WIDE],
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
                CnvStratificationDimension.DATA_BASIS,
                CnvStratificationDimension.BIN_SIZE,
            ],
            coverage_cutpoints_x=[2.0, 5.0, 10.0],
            tumor_fraction_cutpoints=[0.1, 0.2, 0.5],
            event_size_cutpoints_bp=[1_000_000, 5_000_000, 10_000_000],
        ),
        event_classes=[EventType.DELETION, EventType.DUPLICATION],
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
        rationale="Synthetic preregistration contract; no biological performance claim.",
    )


def _specimen(*, tumor_fraction: float | None = 0.25) -> CnvValidationSpecimen:
    return CnvValidationSpecimen(
        specimen_id="SYNTHETIC_CNV_001",
        biological_specimen_id="SYNTHETIC_BIO_001",
        repeat_kind=CnvRepeatKind.INDEPENDENT,
        access_basis="public",
        material_type="synthetic-DNA",
        genome_build=GenomeBuild.GRCH38,
        data_basis=CnvDataBasis.LCWGS_GENOME_WIDE,
        reference_id="synthetic-grch38",
        reference_sha256=_sha("reference"),
        input_sha256=_sha("input"),
        coverage_x=5.0,
        coverage_definition="synthetic mean autosomal depth",
        tumor_fraction=tumor_fraction,
        tumor_fraction_method=("synthetic orthogonal fraction" if tumor_fraction is not None else None),
        tumor_fraction_timepoint=("synthetic same aliquot" if tumor_fraction is not None else None),
        truth_sources=[
            CnvTruthSource(
                method_name="synthetic-karyotype",
                method_version="v1",
                resource_id="synthetic-truth-v1",
                resource_sha256=_sha("truth"),
                provenance_reference="synthetic fixture only",
            )
        ],
        assessability_mask=CnvAssessabilityMask(
            resource_id="synthetic-assessability-v1",
            resource_sha256=_sha("assessability"),
            unit="regions",
            definition="Synthetic locked assessable territory for software-contract tests.",
        ),
        truth_events=[
            GenomicEvent(
                event_id="truth-del-1",
                event_type=EventType.DELETION,
                primary=Locus(chromosome="chr7", start=1_000_000, end=6_000_000),
                copy_number=1.0,
            )
        ],
        negative_universe=CnvNegativeUniverse(
            universe_id="synthetic-negative-bins",
            unit="genomic_bins",
            assessable_units=100,
            resource_sha256=_sha("negative-universe"),
            definition="Synthetic non-event bins after the registered assessability mask.",
        ),
    )


def _cohort(*, tumor_fraction: float | None = 0.25) -> CnvValidationCohort:
    return CnvValidationCohort(specimens=[_specimen(tumor_fraction=tumor_fraction)])


def _registration(
    matrix: CnvValidationMatrix | None = None,
    cohort: CnvValidationCohort | None = None,
) -> CnvValidationRegistration:
    return preregister_cnv_validation(
        matrix or _matrix(),
        cohort or _cohort(),
        registration_id="synthetic-cnv-registration",
        registered_at=datetime(2026, 9, 17, tzinfo=UTC),
        code_sha256=_sha("code"),
        software_version="0.8.2",
    )


class CnvValidationRegistrationTests(unittest.TestCase):
    def test_canonical_hash_is_independent_of_mapping_key_order(self) -> None:
        self.assertEqual(
            canonical_content_sha256({"a": 1, "b": 2}),
            canonical_content_sha256({"b": 2, "a": 1}),
        )

    def test_primary_tumor_fraction_dimension_rejects_unknown_fraction(self) -> None:
        cohort = _cohort(tumor_fraction=None)
        reasons = cnv_cohort_eligibility_reasons(_matrix(), cohort)
        self.assertTrue(any("tumour fraction" in reason.lower() for reason in reasons))
        with self.assertRaises(ValueError):
            _registration(cohort=cohort)

    def test_build_outside_registered_matrix_is_ineligible(self) -> None:
        matrix = _matrix().model_copy(update={"genome_builds": [GenomeBuild.GRCH37]})
        reasons = cnv_cohort_eligibility_reasons(matrix, _cohort())
        self.assertTrue(any("genome build" in reason.lower() for reason in reasons))

    def test_data_basis_outside_registered_matrix_is_ineligible(self) -> None:
        matrix = _matrix().model_copy(
            update={"data_bases": [CnvDataBasis.ADAPTIVE_SAMPLING_OFF_TARGET]}
        )
        reasons = cnv_cohort_eligibility_reasons(matrix, _cohort())
        self.assertTrue(any("data basis" in reason.lower() for reason in reasons))

    def test_truth_event_outside_registered_event_classes_is_ineligible(self) -> None:
        matrix = _matrix().model_copy(update={"event_classes": [EventType.DUPLICATION]})
        reasons = cnv_cohort_eligibility_reasons(matrix, _cohort())
        self.assertTrue(any("event class" in reason.lower() for reason in reasons))

    def test_specificity_question_requires_negative_universe(self) -> None:
        specificity = CnvAcceptanceQuestion(
            question_id="specificity",
            metric=CnvAcceptanceMetric.SPECIFICITY,
            minimum_evaluable_denominator=10,
            minimum_acceptable=0.95,
        )
        matrix = _matrix().model_copy(update={"acceptance": [specificity]})
        specimen_payload = _specimen().model_dump()
        specimen_payload["negative_universe"] = None
        cohort = CnvValidationCohort(
            specimens=[CnvValidationSpecimen.model_validate(specimen_payload)]
        )
        reasons = cnv_cohort_eligibility_reasons(matrix, cohort)
        self.assertTrue(any("negative universe" in reason.lower() for reason in reasons))
        with self.assertRaises(ValueError):
            _registration(matrix=matrix, cohort=cohort)

    def test_registration_lock_changes_when_matrix_changes(self) -> None:
        first = _registration()
        changed = _matrix().model_copy(update={"qdnaseq_bin_sizes_kbp": [100, 500]})
        second = _registration(matrix=changed)
        self.assertNotEqual(first.lock_sha256, second.lock_sha256)
        self.assertNotEqual(first.matrix_sha256, second.matrix_sha256)

    def test_registration_rejects_mutated_embedded_content(self) -> None:
        registration = _registration()
        payload = registration.model_dump(mode="json")
        payload["matrix"]["rationale"] = "Mutated after registration and therefore invalid."
        with self.assertRaises(ValidationError):
            CnvValidationRegistration.model_validate(payload)

    def test_registration_timestamp_requires_timezone(self) -> None:
        with self.assertRaises(ValueError):
            preregister_cnv_validation(
                _matrix(),
                _cohort(),
                registration_id="synthetic-cnv-registration",
                registered_at=datetime(2026, 9, 17),
                code_sha256=_sha("code"),
                software_version="0.8.2",
            )

    def test_registration_declares_outcomes_unseen(self) -> None:
        registration = _registration()
        self.assertTrue(registration.outcomes_unseen)
        self.assertTrue(registration.research_only)


if __name__ == "__main__":
    unittest.main()
