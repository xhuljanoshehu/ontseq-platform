from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from datetime import datetime
from typing import Any, Literal

from pydantic import Field, model_validator

from .cnv_validation_contracts import (
    CnvAcceptanceMetric,
    CnvStratificationDimension,
    CnvValidationCohort,
    CnvValidationMatrix,
)
from .models import StrictModel

SHA256 = r"^[0-9a-f]{64}$"


def canonical_content_sha256(value: StrictModel | Mapping[str, Any]) -> str:
    """Hash canonical JSON content independently of key order and whitespace."""
    payload = value.model_dump(mode="json") if isinstance(value, StrictModel) else dict(value)
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(rendered.encode()).hexdigest()


def cnv_cohort_eligibility_reasons(
    matrix: CnvValidationMatrix,
    cohort: CnvValidationCohort,
) -> list[str]:
    """Return prospective registration blockers without imputing missing study metadata."""
    reasons: list[str] = []
    primary = set(matrix.stratification.primary_dimensions)
    specificity_required = any(
        question.metric == CnvAcceptanceMetric.SPECIFICITY for question in matrix.acceptance
    )

    for specimen in cohort.specimens:
        if specimen.genome_build not in matrix.genome_builds:
            reasons.append(
                f"{specimen.specimen_id}: genome build is outside the registered matrix."
            )
        if specimen.data_basis not in matrix.data_bases:
            reasons.append(f"{specimen.specimen_id}: data basis is outside the registered matrix.")
        if CnvStratificationDimension.COVERAGE in primary and specimen.coverage_x is None:
            reasons.append(
                f"{specimen.specimen_id}: primary coverage stratification requires "
                "measured coverage."
            )
        if CnvStratificationDimension.TUMOR_FRACTION in primary and specimen.tumor_fraction is None:
            reasons.append(
                f"{specimen.specimen_id}: primary tumour fraction stratification requires "
                "an explicit tumour fraction."
            )
        for event in specimen.truth_events:
            if event.event_type not in matrix.event_classes:
                reasons.append(
                    f"{specimen.specimen_id}: truth event {event.event_id} has an event class "
                    "outside the registered matrix."
                )
        if specificity_required and specimen.negative_universe is None:
            reasons.append(
                f"{specimen.specimen_id}: specificity requires a registered negative universe."
            )

    return reasons


class CnvValidationRegistration(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    registration_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    matrix: CnvValidationMatrix
    cohort: CnvValidationCohort
    code_sha256: str = Field(pattern=SHA256)
    software_version: str = Field(min_length=1)
    registered_at: datetime
    outcomes_unseen: Literal[True] = True
    matrix_sha256: str = Field(pattern=SHA256)
    cohort_sha256: str = Field(pattern=SHA256)
    lock_sha256: str = Field(pattern=SHA256)
    timestamp_basis: Literal["operator_declaration_not_trusted_timestamp"] = (
        "operator_declaration_not_trusted_timestamp"
    )
    sensitive_output: Literal[True] = True
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def verify_lock(self) -> CnvValidationRegistration:
        if self.registered_at.utcoffset() is None:
            raise ValueError("Registration timestamp requires a timezone")
        if canonical_content_sha256(self.matrix) != self.matrix_sha256:
            raise ValueError("Matrix content does not match registration checksum")
        if canonical_content_sha256(self.cohort) != self.cohort_sha256:
            raise ValueError("Cohort content does not match registration checksum")
        expected = canonical_content_sha256(self.model_dump(mode="json", exclude={"lock_sha256"}))
        if expected != self.lock_sha256:
            raise ValueError("Registration content does not match its lock")
        return self


def preregister_cnv_validation(
    matrix: CnvValidationMatrix,
    cohort: CnvValidationCohort,
    *,
    registration_id: str,
    registered_at: datetime,
    code_sha256: str,
    software_version: str,
) -> CnvValidationRegistration:
    """Freeze an eligible study declaration before outcome evaluation."""
    reasons = cnv_cohort_eligibility_reasons(matrix, cohort)
    if reasons:
        raise ValueError("CNV validation cohort is ineligible: " + " ".join(reasons))
    if registered_at.utcoffset() is None:
        raise ValueError("Registration timestamp requires a timezone")

    values: dict[str, Any] = {
        "schema_version": "0.1.0",
        "registration_id": registration_id,
        "matrix": matrix.model_dump(mode="json"),
        "cohort": cohort.model_dump(mode="json"),
        "code_sha256": code_sha256,
        "software_version": software_version,
        "registered_at": registered_at.isoformat().replace("+00:00", "Z"),
        "outcomes_unseen": True,
        "matrix_sha256": canonical_content_sha256(matrix),
        "cohort_sha256": canonical_content_sha256(cohort),
        "timestamp_basis": "operator_declaration_not_trusted_timestamp",
        "sensitive_output": True,
        "research_only": True,
    }
    values["lock_sha256"] = canonical_content_sha256(values)
    return CnvValidationRegistration.model_validate(values)
