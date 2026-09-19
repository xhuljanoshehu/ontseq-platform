# CNV Validation Block 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the first CNV-validation block: prospectively defined study/stratification contracts, cohort/truth contracts, fail-closed preregistration locks, schema export, and synthetic tests without executing real biological data or selecting a production caller.

**Architecture:** Add a CNV-specific validation contract layer above the existing `BenchmarkCase`/`benchmark_case()` primitive. Keep study design models separate from preregistration/hash logic. Reuse existing `StrictModel`, `GenomeBuild`, `EventType`, `GenomicEvent`, and `BenchmarkThresholds`; do not alter event matching, QDNAseq/ACE execution, dilution logic, or report generation in this block.

**Tech Stack:** Python 3.11+, Pydantic v2, `unittest`/pytest runner, existing ONTSeq schema exporter, SHA-256 canonical JSON locks.

**Spec:** `docs/superpowers/specs/2026-09-17-cnv-validation-program-design.md`

## Global Constraints

- Research Use Only; no clinical-validity claim.
- No patient genomic payloads, BAM/POD5/FASTQ/VCF content, identifiers, or institutionally restricted truth data in Git.
- Repository tests use synthetic fixtures only.
- Prospective stratification must never delete, overwrite, or hide underlying analytical evidence.
- Continuous measured values remain preserved; categorical strata are derived views only.
- Specificity must never be inferred from event-level TP/FP/FN without a prospectively defined negative denominator.
- No caller becomes a production default in this block.
- Unknown tumour/blast fraction must never be imputed.
- GRCh37 and GRCh38 remain explicit, separate build states.
- Existing benchmark matching and dilution/LoD implementations remain unchanged.
- No merge to `main` and no package-version bump as part of this implementation block.

---

## File structure

- Create `src/ontseq_platform/cnv_validation_contracts.py` — study matrix, strata, caller locks, acceptance-question, specimen/truth/cohort models.
- Create `src/ontseq_platform/cnv_validation_registration.py` — canonical hashing, eligibility checks, registration model and preregistration function.
- Create `tests/test_cnv_validation_contracts.py` — RED/GREEN contract tests.
- Create `tests/test_cnv_validation_registration.py` — RED/GREEN eligibility and lock tests.
- Modify `scripts/export_schemas.py` — export the three new public schemas.
- Create generated `schemas/cnv-validation-matrix.schema.json`.
- Create generated `schemas/cnv-validation-cohort.schema.json`.
- Create generated `schemas/cnv-validation-registration.schema.json`.

No caller adapter, aggregator, CLI, report, or real-data config is created in Block 1.

---

### Task 1: Prospective matrix and stratification contracts

**Files:**
- Create: `tests/test_cnv_validation_contracts.py`
- Create: `src/ontseq_platform/cnv_validation_contracts.py`

**Interfaces:**
- Consumes: `StrictModel`, `GenomeBuild`, `EventType`, `BenchmarkThresholds` from `ontseq_platform.models`.
- Produces:
  - `CnvStratificationRole`
  - `CnvStratificationDimension`
  - `CnvCallerLock`
  - `CnvStratificationPlan`
  - `CnvAcceptanceMetric`
  - `CnvAcceptanceQuestion`
  - `CnvValidationMatrix`

- [ ] **Step 1: Write failing tests for immutable prospective dimensions**

Create `tests/test_cnv_validation_contracts.py` with synthetic-only tests. Start with:

```python
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
        with self.assertRaises(ValidationError):
            _matrix().model_copy(
                update={
                    "stratification": _matrix().stratification.model_copy(
                        update={"coverage_cutpoints_x": [5.0, 2.0, 5.0]}
                    )
                }
            ).model_validate(
                _matrix().model_copy(
                    update={
                        "stratification": _matrix().stratification.model_copy(
                            update={"coverage_cutpoints_x": [5.0, 2.0, 5.0]}
                        )
                    }
                ).model_dump()
            )

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
```

- [ ] **Step 2: Run the contract tests and verify RED**

Run:

```bash
pytest -q tests/test_cnv_validation_contracts.py
```

Expected: collection/import failure because `ontseq_platform.cnv_validation_contracts` does not exist.

- [ ] **Step 3: Implement the minimum contract models**

Create `src/ontseq_platform/cnv_validation_contracts.py`. The public shape must be:

```python
from __future__ import annotations

import math
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from .models import BenchmarkThresholds, EventType, GenomeBuild, StrictModel

SHA256 = r"^[0-9a-f]{64}$"
CNV_EVENT_TYPES = frozenset(
    {
        EventType.CHROMOSOME_GAIN,
        EventType.CHROMOSOME_LOSS,
        EventType.DELETION,
        EventType.DUPLICATION,
    }
)


class CnvStratificationRole(StrEnum):
    PRIMARY = "primary"
    SECONDARY = "secondary"
    EXPLORATORY = "exploratory"


class CnvStratificationDimension(StrEnum):
    CALLER = "caller"
    GENOME_BUILD = "genome_build"
    COVERAGE = "coverage"
    TUMOR_FRACTION = "tumor_fraction"
    EVENT_CLASS = "event_class"
    EVENT_SIZE = "event_size"
    BIN_SIZE = "bin_size"


class CnvCallerLock(StrictModel):
    caller_id: str = Field(pattern=r"^[a-z][a-z0-9_]{1,63}$")
    caller_version: str = Field(min_length=1)
    adapter_policy_sha256: str = Field(pattern=SHA256)
    execution_identity_sha256: str = Field(pattern=SHA256)


class CnvStratificationPlan(StrictModel):
    primary_dimensions: list[CnvStratificationDimension] = Field(min_length=1)
    secondary_dimensions: list[CnvStratificationDimension] = Field(default_factory=list)
    exploratory_dimensions: list[CnvStratificationDimension] = Field(default_factory=list)
    coverage_cutpoints_x: list[float] = Field(default_factory=list)
    tumor_fraction_cutpoints: list[float] = Field(default_factory=list)
    event_size_cutpoints_bp: list[int] = Field(default_factory=list)

    @model_validator(mode="after")
    def coherent_strata(self) -> "CnvStratificationPlan":
        groups = [self.primary_dimensions, self.secondary_dimensions, self.exploratory_dimensions]
        flattened = [item for group in groups for item in group]
        if len(flattened) != len(set(flattened)):
            raise ValueError("A stratification dimension may have exactly one role")
        if self.coverage_cutpoints_x != sorted(set(self.coverage_cutpoints_x)):
            raise ValueError("Coverage cutpoints must be unique and strictly sorted")
        if any(not math.isfinite(x) or x <= 0 for x in self.coverage_cutpoints_x):
            raise ValueError("Coverage cutpoints must be finite and greater than zero")
        if self.tumor_fraction_cutpoints != sorted(set(self.tumor_fraction_cutpoints)):
            raise ValueError("Tumour-fraction cutpoints must be unique and strictly sorted")
        if any(not math.isfinite(x) or x <= 0 or x >= 1 for x in self.tumor_fraction_cutpoints):
            raise ValueError("Tumour-fraction cutpoints must lie strictly between zero and one")
        if self.event_size_cutpoints_bp != sorted(set(self.event_size_cutpoints_bp)):
            raise ValueError("Event-size cutpoints must be unique and strictly sorted")
        if any(x <= 0 for x in self.event_size_cutpoints_bp):
            raise ValueError("Event-size cutpoints must be positive")
        return self


class CnvAcceptanceMetric(StrEnum):
    SENSITIVITY = "sensitivity"
    PRECISION = "precision"
    F1 = "f1"
    FALSE_POSITIVE_BURDEN = "false_positive_burden"
    NO_CALL_RATE = "no_call_rate"
    TECHNICAL_FAILURE_RATE = "technical_failure_rate"
    NOT_ASSESSABLE_RATE = "not_assessable_rate"
    SPECIFICITY = "specificity"
    COPY_NUMBER_ERROR = "copy_number_error"
    CELLULARITY_ERROR = "cellularity_error"
    PLOIDY_ERROR = "ploidy_error"
    REPRODUCIBILITY = "reproducibility"


class CnvAcceptanceQuestion(StrictModel):
    question_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    metric: CnvAcceptanceMetric
    minimum_evaluable_denominator: int = Field(ge=1)
    minimum_acceptable: float | None = None
    maximum_acceptable: float | None = None
    filters: dict[str, str | int | float | bool] = Field(default_factory=dict)

    @model_validator(mode="after")
    def bounded_question(self) -> "CnvAcceptanceQuestion":
        bounds = [self.minimum_acceptable, self.maximum_acceptable]
        if all(value is None for value in bounds):
            raise ValueError("Acceptance question requires a minimum or maximum bound")
        if any(value is not None and not math.isfinite(value) for value in bounds):
            raise ValueError("Acceptance bounds must be finite")
        if (
            self.minimum_acceptable is not None
            and self.maximum_acceptable is not None
            and self.minimum_acceptable > self.maximum_acceptable
        ):
            raise ValueError("Acceptance minimum cannot exceed maximum")
        return self


class CnvValidationMatrix(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    matrix_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    callers: list[CnvCallerLock] = Field(min_length=1)
    genome_builds: list[GenomeBuild] = Field(min_length=1)
    stratification: CnvStratificationPlan
    event_classes: list[EventType] = Field(min_length=1)
    qdnaseq_bin_sizes_kbp: list[int] = Field(default_factory=list)
    matching_thresholds: BenchmarkThresholds = Field(default_factory=BenchmarkThresholds)
    acceptance: list[CnvAcceptanceQuestion] = Field(min_length=1)
    retain_all_evidence: Literal[True] = True
    rationale: str = Field(min_length=20)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_matrix(self) -> "CnvValidationMatrix":
        caller_ids = [item.caller_id for item in self.callers]
        if len(caller_ids) != len(set(caller_ids)):
            raise ValueError("Caller IDs must be unique")
        if self.genome_builds != list(dict.fromkeys(self.genome_builds)):
            raise ValueError("Genome builds must be unique")
        if any(event not in CNV_EVENT_TYPES for event in self.event_classes):
            raise ValueError("CNV validation matrix contains a non-CNV event class")
        if len(self.event_classes) != len(set(self.event_classes)):
            raise ValueError("CNV event classes must be unique")
        if self.qdnaseq_bin_sizes_kbp != sorted(set(self.qdnaseq_bin_sizes_kbp)):
            raise ValueError("QDNAseq bin sizes must be unique and sorted")
        if any(value <= 0 for value in self.qdnaseq_bin_sizes_kbp):
            raise ValueError("QDNAseq bin sizes must be positive")
        question_ids = [item.question_id for item in self.acceptance]
        if len(question_ids) != len(set(question_ids)):
            raise ValueError("Acceptance question IDs must be unique")
        return self
```

Do not add bin-assignment functions yet; continuous-to-categorical derivation belongs to the later evidence/aggregation block.

- [ ] **Step 4: Run the contract tests and verify GREEN**

Run:

```bash
pytest -q tests/test_cnv_validation_contracts.py
```

Expected: PASS.

- [ ] **Step 5: Run targeted existing benchmark/dilution regression tests**

Run:

```bash
pytest -q tests/test_benchmark.py tests/test_dilution.py
```

Expected: PASS with no modifications to benchmark or dilution behavior.

- [ ] **Step 6: Commit Task 1**

```bash
git add src/ontseq_platform/cnv_validation_contracts.py tests/test_cnv_validation_contracts.py
git commit -m "feat(cnv): add validation matrix contracts"
```

---

### Task 2: Cohort, orthogonal truth, and explicit negative-denominator contracts

**Files:**
- Modify: `src/ontseq_platform/cnv_validation_contracts.py`
- Extend: `tests/test_cnv_validation_contracts.py`

**Interfaces:**
- Consumes: Task 1 matrix types plus `GenomicEvent`.
- Produces:
  - `CnvRepeatKind`
  - `CnvTruthSource`
  - `CnvNegativeUniverse`
  - `CnvValidationSpecimen`
  - `CnvValidationCohort`

- [ ] **Step 1: Add failing tests for explicit truth and missing-data semantics**

Append tests that require:

```python
from ontseq_platform.cnv_validation_contracts import (
    CnvNegativeUniverse,
    CnvRepeatKind,
    CnvTruthSource,
    CnvValidationCohort,
    CnvValidationSpecimen,
)
from ontseq_platform.models import GenomicEvent, Locus


def _truth_event() -> GenomicEvent:
    return GenomicEvent(
        event_id="truth-del-1",
        event_type=EventType.DELETION,
        primary=Locus(chromosome="chr7", start=1_000_000, end=6_000_000),
        copy_number=1.0,
    )


def _specimen(*, tumor_fraction: float | None = 0.25) -> CnvValidationSpecimen:
    return CnvValidationSpecimen(
        specimen_id="SYNTHETIC_CNV_001",
        biological_specimen_id="SYNTHETIC_BIO_001",
        repeat_kind=CnvRepeatKind.INDEPENDENT,
        access_basis="public",
        material_type="synthetic-DNA",
        genome_build=GenomeBuild.GRCH38,
        reference_id="synthetic-grch38",
        reference_sha256=_sha("synthetic-reference"),
        input_sha256=_sha("synthetic-input"),
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
                resource_sha256=_sha("synthetic-truth"),
                provenance_reference="synthetic fixture only",
            )
        ],
        truth_events=[_truth_event()],
        negative_universe=CnvNegativeUniverse(
            universe_id="synthetic-negative-bins",
            unit="genomic_bins",
            assessable_units=100,
            resource_sha256=_sha("synthetic-negative-universe"),
            definition="Synthetic non-event bins after a locked mask.",
        ),
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


def test_cohort_requires_unique_specimen_ids(self) -> None:
    with self.assertRaises(ValidationError):
        CnvValidationCohort(specimens=[_specimen(), _specimen()])
```

- [ ] **Step 2: Run the new tests and verify RED**

Run:

```bash
pytest -q tests/test_cnv_validation_contracts.py
```

Expected: FAIL because cohort/truth types are not defined.

- [ ] **Step 3: Implement the minimum cohort/truth models**

Add to `cnv_validation_contracts.py`:

```python
from .models import GenomicEvent


class CnvRepeatKind(StrEnum):
    INDEPENDENT = "independent"
    WITHIN_RUN = "within_run"
    BETWEEN_RUN = "between_run"


class CnvTruthSource(StrictModel):
    method_name: str = Field(min_length=1)
    method_version: str = Field(min_length=1)
    resource_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,127}$")
    resource_sha256: str = Field(pattern=SHA256)
    provenance_reference: str = Field(min_length=3)
    orthogonal_to_evaluated_caller: Literal[True] = True


class CnvNegativeUniverse(StrictModel):
    universe_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,127}$")
    unit: Literal["genomic_bins", "regions"]
    assessable_units: int = Field(ge=1)
    resource_sha256: str = Field(pattern=SHA256)
    definition: str = Field(min_length=20)


class CnvValidationSpecimen(StrictModel):
    specimen_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    biological_specimen_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    repeat_kind: CnvRepeatKind
    repeat_group_id: str | None = Field(default=None, pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    access_basis: Literal["public", "institutionally_authorized"]
    material_type: str = Field(min_length=3)
    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    reference_sha256: str = Field(pattern=SHA256)
    input_sha256: str = Field(pattern=SHA256)
    coverage_x: float | None = Field(default=None, ge=0)
    coverage_definition: str | None = None
    tumor_fraction: float | None = Field(default=None, ge=0, le=1)
    tumor_fraction_method: str | None = None
    tumor_fraction_timepoint: str | None = None
    truth_sources: list[CnvTruthSource] = Field(min_length=1)
    truth_events: list[GenomicEvent] = Field(default_factory=list)
    negative_universe: CnvNegativeUniverse | None = None

    @model_validator(mode="after")
    def coherent_specimen(self) -> "CnvValidationSpecimen":
        if self.coverage_x is None and self.coverage_definition is not None:
            raise ValueError("Coverage definition cannot exist without a measured coverage value")
        if self.coverage_x is not None and not self.coverage_definition:
            raise ValueError("Measured coverage requires its definition")
        fraction_metadata = [self.tumor_fraction_method, self.tumor_fraction_timepoint]
        if self.tumor_fraction is None and any(value is not None for value in fraction_metadata):
            raise ValueError("Tumour-fraction metadata cannot imply a missing fraction")
        if self.tumor_fraction is not None and any(not value for value in fraction_metadata):
            raise ValueError("Measured tumour fraction requires method and timepoint")
        if any(event.event_type not in CNV_EVENT_TYPES for event in self.truth_events):
            raise ValueError("Truth events must be CNV event types")
        event_ids = [event.event_id for event in self.truth_events]
        if len(event_ids) != len(set(event_ids)):
            raise ValueError("Truth event IDs must be unique within a specimen")
        if self.repeat_kind == CnvRepeatKind.INDEPENDENT and self.repeat_group_id is not None:
            raise ValueError("Independent specimens cannot declare a repeat group")
        if self.repeat_kind != CnvRepeatKind.INDEPENDENT and self.repeat_group_id is None:
            raise ValueError("Repeat specimens require repeat_group_id")
        return self


class CnvValidationCohort(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    specimens: list[CnvValidationSpecimen] = Field(min_length=1)
    sensitive_output: Literal[True] = True
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def unique_specimens(self) -> "CnvValidationCohort":
        ids = [item.specimen_id for item in self.specimens]
        if len(ids) != len(set(ids)):
            raise ValueError("Validation specimen IDs must be unique")
        return self
```

- [ ] **Step 4: Run contract tests and verify GREEN**

Run:

```bash
pytest -q tests/test_cnv_validation_contracts.py
```

Expected: PASS.

- [ ] **Step 5: Commit Task 2**

```bash
git add src/ontseq_platform/cnv_validation_contracts.py tests/test_cnv_validation_contracts.py
git commit -m "feat(cnv): add cohort and truth contracts"
```

---

### Task 3: Fail-closed cohort eligibility and preregistration lock

**Files:**
- Create: `tests/test_cnv_validation_registration.py`
- Create: `src/ontseq_platform/cnv_validation_registration.py`

**Interfaces:**
- Consumes: `CnvValidationMatrix`, `CnvValidationCohort`.
- Produces:
  - `canonical_content_sha256(value) -> str`
  - `cnv_cohort_eligibility_reasons(matrix, cohort) -> list[str]`
  - `CnvValidationRegistration`
  - `preregister_cnv_validation(...) -> CnvValidationRegistration`

- [ ] **Step 1: Write failing registration/eligibility tests**

Create `tests/test_cnv_validation_registration.py` and reuse local synthetic helper builders. Tests must assert at least:

```python
def test_primary_tumor_fraction_dimension_rejects_unknown_fraction_for_registration(self) -> None:
    matrix = _matrix()
    cohort = CnvValidationCohort(specimens=[_specimen(tumor_fraction=None)])
    reasons = cnv_cohort_eligibility_reasons(matrix, cohort)
    self.assertTrue(any("tumour fraction" in reason.lower() for reason in reasons))
    with self.assertRaises(ValueError):
        preregister_cnv_validation(
            matrix,
            cohort,
            registration_id="synthetic-cnv-registration",
            registered_at=datetime(2026, 9, 17, tzinfo=UTC),
            code_sha256=_sha("code"),
            software_version="0.8.2",
        )


def test_build_outside_registered_matrix_is_ineligible(self) -> None:
    matrix = _matrix().model_copy(update={"genome_builds": [GenomeBuild.GRCH37]})
    cohort = CnvValidationCohort(specimens=[_specimen()])
    reasons = cnv_cohort_eligibility_reasons(matrix, cohort)
    self.assertTrue(any("genome build" in reason.lower() for reason in reasons))


def test_specificity_question_requires_negative_universe(self) -> None:
    payload = _matrix().model_dump()
    payload["acceptance"] = [
        CnvAcceptanceQuestion(
            question_id="specificity",
            metric=CnvAcceptanceMetric.SPECIFICITY,
            minimum_evaluable_denominator=10,
            minimum_acceptable=0.95,
        ).model_dump()
    ]
    matrix = CnvValidationMatrix.model_validate(payload)
    specimen_payload = _specimen().model_dump()
    specimen_payload["negative_universe"] = None
    cohort = CnvValidationCohort(
        specimens=[CnvValidationSpecimen.model_validate(specimen_payload)]
    )
    reasons = cnv_cohort_eligibility_reasons(matrix, cohort)
    self.assertTrue(any("negative universe" in reason.lower() for reason in reasons))


def test_registration_lock_changes_when_matrix_changes(self) -> None:
    first = preregister_cnv_validation(...)
    changed_matrix = _matrix().model_copy(
        update={"qdnaseq_bin_sizes_kbp": [100, 500]}
    )
    second = preregister_cnv_validation(... changed_matrix ...)
    self.assertNotEqual(first.lock_sha256, second.lock_sha256)


def test_registration_rejects_mutated_embedded_content(self) -> None:
    registration = preregister_cnv_validation(...)
    payload = registration.model_dump(mode="json")
    payload["matrix"]["rationale"] = "Mutated after registration and therefore invalid."
    with self.assertRaises(ValidationError):
        CnvValidationRegistration.model_validate(payload)


def test_registration_timestamp_requires_timezone(self) -> None:
    with self.assertRaises(ValueError):
        preregister_cnv_validation(
            _matrix(),
            CnvValidationCohort(specimens=[_specimen()]),
            registration_id="synthetic-cnv-registration",
            registered_at=datetime(2026, 9, 17),
            code_sha256=_sha("code"),
            software_version="0.8.2",
        )
```

- [ ] **Step 2: Run registration tests and verify RED**

Run:

```bash
pytest -q tests/test_cnv_validation_registration.py
```

Expected: import failure because registration module does not exist.

- [ ] **Step 3: Implement canonical content hashing**

Create `src/ontseq_platform/cnv_validation_registration.py` with canonical JSON hashing independent of whitespace/key order:

```python
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
    payload = value.model_dump(mode="json") if isinstance(value, StrictModel) else dict(value)
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(rendered.encode()).hexdigest()
```

- [ ] **Step 4: Implement eligibility checks without imputation**

Implement:

```python
def cnv_cohort_eligibility_reasons(
    matrix: CnvValidationMatrix, cohort: CnvValidationCohort
) -> list[str]:
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
        if CnvStratificationDimension.COVERAGE in primary and specimen.coverage_x is None:
            reasons.append(
                f"{specimen.specimen_id}: primary coverage stratification requires measured coverage."
            )
        if (
            CnvStratificationDimension.TUMOR_FRACTION in primary
            and specimen.tumor_fraction is None
        ):
            reasons.append(
                f"{specimen.specimen_id}: primary tumour fraction stratification requires an explicit tumour fraction."
            )
        if specificity_required and specimen.negative_universe is None:
            reasons.append(
                f"{specimen.specimen_id}: specificity requires a registered negative universe."
            )
    return reasons
```

Do not invent coverage, tumour fraction, negative regions, ploidy, or truth data.

- [ ] **Step 5: Implement registration model and preregistration function**

Required shape:

```python
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
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def verify_lock(self) -> "CnvValidationRegistration":
        if self.registered_at.utcoffset() is None:
            raise ValueError("Registration timestamp requires a timezone")
        if canonical_content_sha256(self.matrix) != self.matrix_sha256:
            raise ValueError("Matrix content does not match registration checksum")
        if canonical_content_sha256(self.cohort) != self.cohort_sha256:
            raise ValueError("Cohort content does not match registration checksum")
        expected = canonical_content_sha256(
            self.model_dump(mode="json", exclude={"lock_sha256"})
        )
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
    reasons = cnv_cohort_eligibility_reasons(matrix, cohort)
    if reasons:
        raise ValueError("CNV validation cohort is ineligible: " + " ".join(reasons))
    values: dict[str, Any] = {
        "schema_version": "0.1.0",
        "registration_id": registration_id,
        "matrix": matrix.model_dump(mode="json"),
        "cohort": cohort.model_dump(mode="json"),
        "code_sha256": code_sha256,
        "software_version": software_version,
        "registered_at": registered_at.isoformat(),
        "outcomes_unseen": True,
        "matrix_sha256": canonical_content_sha256(matrix),
        "cohort_sha256": canonical_content_sha256(cohort),
        "timestamp_basis": "operator_declaration_not_trusted_timestamp",
        "research_only": True,
    }
    values["lock_sha256"] = canonical_content_sha256(values)
    return CnvValidationRegistration.model_validate(values)
```

- [ ] **Step 6: Run registration tests and verify GREEN**

Run:

```bash
pytest -q tests/test_cnv_validation_registration.py
```

Expected: PASS.

- [ ] **Step 7: Run both new test files together**

Run:

```bash
pytest -q tests/test_cnv_validation_contracts.py tests/test_cnv_validation_registration.py
```

Expected: PASS.

- [ ] **Step 8: Commit Task 3**

```bash
git add src/ontseq_platform/cnv_validation_registration.py tests/test_cnv_validation_registration.py
git commit -m "feat(cnv): add preregistered validation locks"
```

---

### Task 4: Export schemas and verify repository integration

**Files:**
- Modify: `scripts/export_schemas.py`
- Create: `schemas/cnv-validation-matrix.schema.json`
- Create: `schemas/cnv-validation-cohort.schema.json`
- Create: `schemas/cnv-validation-registration.schema.json`
- Extend: `tests/test_cnv_validation_contracts.py`

**Interfaces:**
- Consumes: public contract and registration models from Tasks 1-3.
- Produces: deterministic JSON Schema artifacts checked by the existing `scripts/export_schemas.py --check` gate.

- [ ] **Step 1: Add failing schema-registration test**

Add to `tests/test_cnv_validation_contracts.py`:

```python
def test_public_cnv_validation_models_emit_json_schema(self) -> None:
    for model in (CnvValidationMatrix, CnvValidationCohort):
        schema = model.model_json_schema()
        self.assertEqual(schema["type"], "object")
        self.assertFalse(schema.get("additionalProperties", True))
```

Add an analogous assertion for `CnvValidationRegistration` in `tests/test_cnv_validation_registration.py`.

- [ ] **Step 2: Run tests and verify they are still GREEN before exporter wiring**

Run:

```bash
pytest -q tests/test_cnv_validation_contracts.py tests/test_cnv_validation_registration.py
```

Expected: PASS. This confirms schema generation works before filesystem exporter integration.

- [ ] **Step 3: Wire models into `scripts/export_schemas.py`**

Add imports:

```python
from ontseq_platform.cnv_validation_contracts import CnvValidationCohort, CnvValidationMatrix
from ontseq_platform.cnv_validation_registration import CnvValidationRegistration
```

Add to the existing model dictionary in `_render()`:

```python
"cnv-validation-matrix": CnvValidationMatrix,
"cnv-validation-cohort": CnvValidationCohort,
"cnv-validation-registration": CnvValidationRegistration,
```

- [ ] **Step 4: Generate the schemas**

Run:

```bash
python scripts/export_schemas.py
```

Expected new files:

```text
schemas/cnv-validation-matrix.schema.json
schemas/cnv-validation-cohort.schema.json
schemas/cnv-validation-registration.schema.json
```

- [ ] **Step 5: Verify schema freshness gate**

Run:

```bash
python scripts/export_schemas.py --check
```

Expected: exit code 0 and no stale-schema message.

- [ ] **Step 6: Run targeted quality gates**

Run:

```bash
pytest -q \
  tests/test_cnv_validation_contracts.py \
  tests/test_cnv_validation_registration.py \
  tests/test_benchmark.py \
  tests/test_dilution.py
ruff check src/ontseq_platform/cnv_validation_contracts.py \
  src/ontseq_platform/cnv_validation_registration.py \
  tests/test_cnv_validation_contracts.py \
  tests/test_cnv_validation_registration.py \
  scripts/export_schemas.py
mypy src/ontseq_platform/cnv_validation_contracts.py \
  src/ontseq_platform/cnv_validation_registration.py
```

Expected: all PASS.

- [ ] **Step 7: Run the repository CI-equivalent Python test suite**

Run the same test command used by the repository CI for the Python suite (read `.github/workflows/ci.yml` if it changed). At minimum:

```bash
pytest -q
```

Expected: no regressions outside the new subsystem.

- [ ] **Step 8: Commit Task 4**

```bash
git add \
  scripts/export_schemas.py \
  schemas/cnv-validation-matrix.schema.json \
  schemas/cnv-validation-cohort.schema.json \
  schemas/cnv-validation-registration.schema.json \
  tests/test_cnv_validation_contracts.py \
  tests/test_cnv_validation_registration.py
git commit -m "chore(cnv): export validation contract schemas"
```

---

## Plan self-review

### Spec coverage

Block 1 covers the design spec's requested matrix models, cohort/truth models, prospective strata, acceptance questions, content locking, preregistration, eligibility, schema export, and synthetic unit tests. Full evidence rows, aggregate metrics, QDNAseq/ACE adapter mapping, independent caller execution, and real validation are deliberately deferred to Blocks 2-6.

### Data-retention boundary

This block enforces `retain_all_evidence=True` in the matrix but does not yet create evidence rows. The actual proof that excluded records remain retrievable belongs to Block 2 and must not be falsely claimed after Block 1.

### Scientific boundary

No event-size, coverage, tumour-fraction, sensitivity, specificity, LoD, cellularity, or ploidy value in repository tests is a scientific threshold. Test values are synthetic software fixtures. Real acceptance cut-points require a separately registered study matrix before outcomes are inspected.

### Implementation stop condition

Stop Block 1 if any existing benchmark/dilution behavior has to be changed to make these contracts work. That would indicate a boundary error and requires design review rather than widening this block.
