from __future__ import annotations

import hashlib
import unittest

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
            ),
            CnvCallerLock(
                caller_id="spectre",
                caller_version="synthetic-1",
                adapter_policy_sha256=_sha("spectre-policy"),
                execution_identity_sha256=_sha("spectre-runtime"),
            ),
        ],
        genome_builds=[GenomeBuild.GRCH37, GenomeBuild.GRCH38],
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


def _truth_event() -> GenomicEvent:
    return GenomicEvent(
        event_id="truth-del-1",
        event_type=EventType.DELETION,
        primary=Locus(chromosome="chr7", start=1_000_000, end=6_000_000),
        copy_number=1.0,
    )


def _specimen(*, tumor_fraction: float | None = 0.25) -> CnvValidationSpecimen:
    assessability_sha256 = _sha("synthetic-assessability")
    return CnvValidationSpecimen(
        specimen_id="SYNTHETIC_CNV_001",
        biological_specimen_id="SYNTHETIC_BIO_001",
        repeat_kind=CnvRepeatKind.INDEPENDENT,
        access_basis="public",
        material_type="synthetic-DNA",
        genome_build=GenomeBuild.GRCH38,
        data_basis=CnvDataBasis.LCWGS_GENOME_WIDE,
        reference_id="synthetic-grch38",
        reference_sha256=_sha("synthetic-reference"),
        input_sha256=_sha("synthetic-input"),
        coverage_x=5.0,
        coverage_definition="synthetic mean autosomal depth",
        tumor_fraction=tumor_fraction,
        tumor_fraction_method=(
            "synthetic orthogonal fraction" if tumor_fraction is not None else None
        ),
        tumor_fraction_timepoint=("synthetic same aliquot" if tumor_fraction is not None else None),
        truth_sources=[
            CnvTruthSource(
                method_name="synthetic-karyotype",
                method_version="v1",
                resource_id="synthetic-truth-v1",
                resource_sha256=_sha("synthetic-truth"),
                provenance_reference="synthetic fixture only",
            )
        ],
        assessability_mask=CnvAssessabilityMask(
            resource_id="synthetic-assessability-v1",
            resource_sha256=assessability_sha256,
            unit="regions",
            definition="Synthetic locked assessable territory for software-contract tests.",
        ),
        truth_events=[_truth_event()],
        negative_universe=CnvNegativeUniverse(
            universe_id="synthetic-negative-bins",
            unit="genomic_bins",
            assessable_units=100,
            resource_sha256=_sha("synthetic-negative-universe"),
            assessability_mask_sha256=assessability_sha256,
            definition="Synthetic non-event bins after a locked mask.",
        ),
    )


class CnvValidationContractTests(unittest.TestCase):
    def test_matrix_keeps_dimensions_and_cutpoints_explicit(self) -> None:
        matrix = _matrix()
        self.assertEqual(matrix.stratification.coverage_cutpoints_x, [2.0, 5.0, 10.0])
        self.assertEqual(matrix.stratification.tumor_fraction_cutpoints, [0.1, 0.2, 0.5])
        self.assertEqual(matrix.qdnaseq_bin_sizes_kbp, [100, 500, 1000])
        self.assertEqual(matrix.data_bases, [CnvDataBasis.LCWGS_GENOME_WIDE])
        self.assertTrue(matrix.retain_all_evidence)

    def test_cutpoints_must_be_unique_and_strictly_sorted(self) -> None:
        payload = _matrix().model_dump()
        payload["stratification"]["coverage_cutpoints_x"] = [5.0, 2.0, 5.0]
        with self.assertRaises(ValidationError):
            CnvValidationMatrix.model_validate(payload)

    def test_registered_numeric_dimension_requires_cutpoints(self) -> None:
        payload = _matrix().model_dump()
        payload["stratification"]["coverage_cutpoints_x"] = []
        with self.assertRaises(ValidationError):
            CnvValidationMatrix.model_validate(payload)

        payload = _matrix().model_dump()
        payload["stratification"]["tumor_fraction_cutpoints"] = []
        with self.assertRaises(ValidationError):
            CnvValidationMatrix.model_validate(payload)

        payload = _matrix().model_dump()
        payload["stratification"]["event_size_cutpoints_bp"] = []
        with self.assertRaises(ValidationError):
            CnvValidationMatrix.model_validate(payload)

    def test_qdnaseq_caller_requires_registered_bin_sizes(self) -> None:
        payload = _matrix().model_dump()
        payload["qdnaseq_bin_sizes_kbp"] = []
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

    def test_unknown_tumor_fraction_is_explicit_and_not_imputed(self) -> None:
        specimen = _specimen(tumor_fraction=None)
        self.assertIsNone(specimen.tumor_fraction)
        self.assertIsNone(specimen.tumor_fraction_method)
        self.assertIsNone(specimen.tumor_fraction_timepoint)

    def test_partial_tumor_fraction_metadata_is_rejected(self) -> None:
        payload = _specimen().model_dump()
        payload["tumor_fraction_method"] = None
        with self.assertRaises(ValidationError):
            CnvValidationSpecimen.model_validate(payload)

    def test_truth_rejects_non_cnv_event_types(self) -> None:
        payload = _specimen().model_dump()
        payload["truth_events"] = [
            GenomicEvent(
                event_id="truth-inv-1",
                event_type=EventType.INVERSION,
                primary=Locus(chromosome="chr7", start=1_000_000, end=6_000_000),
            ).model_dump()
        ]
        with self.assertRaises(ValidationError):
            CnvValidationSpecimen.model_validate(payload)

    def test_negative_universe_is_optional_but_never_invented(self) -> None:
        payload = _specimen().model_dump()
        payload["negative_universe"] = None
        specimen = CnvValidationSpecimen.model_validate(payload)
        self.assertIsNone(specimen.negative_universe)

    def test_negative_universe_is_bound_to_same_assessability_mask(self) -> None:
        payload = _specimen().model_dump()
        payload["negative_universe"]["assessability_mask_sha256"] = _sha("different-mask")
        with self.assertRaises(ValidationError):
            CnvValidationSpecimen.model_validate(payload)

    def test_assessability_mask_is_required(self) -> None:
        payload = _specimen().model_dump()
        payload.pop("assessability_mask")
        with self.assertRaises(ValidationError):
            CnvValidationSpecimen.model_validate(payload)

    def test_data_basis_is_preserved_explicitly(self) -> None:
        specimen = _specimen()
        self.assertEqual(specimen.data_basis, CnvDataBasis.LCWGS_GENOME_WIDE)

    def test_cohort_requires_unique_specimen_ids(self) -> None:
        with self.assertRaises(ValidationError):
            CnvValidationCohort(specimens=[_specimen(), _specimen()])

    def test_repeat_specimen_requires_repeat_group(self) -> None:
        payload = _specimen().model_dump()
        payload["repeat_kind"] = CnvRepeatKind.BETWEEN_RUN
        with self.assertRaises(ValidationError):
            CnvValidationSpecimen.model_validate(payload)


if __name__ == "__main__":
    unittest.main()
