from __future__ import annotations

import math
from datetime import datetime
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from .models import FileFingerprint, GenomeBuild, ModuleRunStatus, StrictModel

_SHA256 = r"^[0-9a-f]{64}$"
_CANONICAL_CHROMOSOME = r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$"


class MarlinSourceKind(StrEnum):
    PRECOMPUTED_METHYLATION = "PRECOMPUTED_METHYLATION"
    MODKIT_DERIVED = "MODKIT_DERIVED"


class MarlinClassificationDecision(StrEnum):
    HIGH_CONFIDENCE = "HIGH_CONFIDENCE"
    UNKNOWN = "UNKNOWN"


class MarlinDualRuntimeVerdict(StrEnum):
    PASS = "PASS"
    FAIL = "FAIL"


class MarlinProbeObservation(StrictModel):
    chromosome: str = Field(pattern=_CANONICAL_CHROMOSOME)
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    methylation_fraction: float | None = Field(default=None, ge=0, le=1)
    probe_id: str = Field(min_length=1)

    @model_validator(mode="after")
    def validate_interval_and_fraction(self) -> MarlinProbeObservation:
        if self.end <= self.start:
            raise ValueError("MARLIN probe end must be greater than start")
        if self.methylation_fraction is not None and not math.isfinite(self.methylation_fraction):
            raise ValueError("MARLIN methylation fraction must be finite")
        return self


class MarlinFeatureSummary(StrictModel):
    expected_feature_count: Literal[357340] = 357340
    observed_model_feature_count: int = Field(ge=0, le=357340)
    explicit_na_feature_count: int = Field(ge=0, le=357340)
    absent_feature_count: int = Field(ge=0, le=357340)
    non_model_probe_count: int = Field(ge=0)
    observed_fraction: float = Field(ge=0, le=1)
    feature_vector_sha256: str = Field(pattern=_SHA256)
    feature_artifact_sha256: str = Field(pattern=_SHA256)
    preprocessing_contract_version: Literal["marlin-v1-binarize-1"] = "marlin-v1-binarize-1"

    @model_validator(mode="after")
    def validate_partition_and_fraction(self) -> MarlinFeatureSummary:
        partition = (
            self.observed_model_feature_count
            + self.explicit_na_feature_count
            + self.absent_feature_count
        )
        if partition != self.expected_feature_count:
            raise ValueError("MARLIN feature counts must partition the 357340 model features")
        expected_fraction = self.observed_model_feature_count / self.expected_feature_count
        if not math.isclose(self.observed_fraction, expected_fraction, rel_tol=0, abs_tol=1e-12):
            raise ValueError("observed_fraction does not match observed model-feature count")
        return self


class MarlinModelUnitScore(StrictModel):
    model_id: int = Field(ge=1, le=42)
    label: str = Field(min_length=1)
    score: float = Field(ge=0, le=1)

    @model_validator(mode="after")
    def score_is_finite(self) -> MarlinModelUnitScore:
        if not math.isfinite(self.score):
            raise ValueError("MARLIN model-unit score must be finite")
        return self


class MarlinGroupedScore(StrictModel):
    label: str = Field(min_length=1)
    score: float = Field(ge=0, le=1.00001)

    @model_validator(mode="after")
    def score_is_finite(self) -> MarlinGroupedScore:
        if not math.isfinite(self.score):
            raise ValueError("MARLIN grouped score must be finite")
        return self


class MarlinArtifactLock(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    lock_id: str = Field(min_length=3)
    model_version: str = Field(min_length=1)
    model_source_uri: str = Field(min_length=1)
    model_sha256: str = Field(pattern=_SHA256)
    code_source_uri: str = Field(min_length=1)
    code_version_or_commit: str = Field(min_length=1)
    code_manifest_sha256: str = Field(pattern=_SHA256)
    feature_source_uri: str = Field(min_length=1)
    feature_sha256: str = Field(pattern=_SHA256)
    canonical_feature_list_sha256: str = Field(pattern=_SHA256)
    class_annotation_source_uri: str = Field(min_length=1)
    class_annotation_sha256: str = Field(pattern=_SHA256)
    probe_resource_source_uri: str = Field(min_length=1)
    probe_resource_sha256: str = Field(pattern=_SHA256)
    genome_build: GenomeBuild
    expected_feature_count: Literal[357340] = 357340
    expected_model_unit_count: Literal[42] = 42
    preprocessing_contract_version: Literal["marlin-v1-binarize-1"] = "marlin-v1-binarize-1"
    runtime_contract_version: Literal["marlin-runtime-v1"] = "marlin-runtime-v1"
    R_version: str = Field(min_length=1)
    keras_version: str = Field(min_length=1)
    tensorflow_version: str = Field(min_length=1)
    python_version_if_used: str | None = Field(default=None, min_length=1)
    execution_backend: str = Field(min_length=1)
    created_at: datetime
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def lock_is_v1_grch37(self) -> MarlinArtifactLock:
        if self.genome_build is not GenomeBuild.GRCH37:
            raise ValueError("MARLIN v1 ONTSeq artifact lock currently requires GRCh37/hg19")
        if self.created_at.utcoffset() is None:
            raise ValueError("MARLIN artifact-lock timestamp requires a timezone")
        return self


class MarlinArtifactSetIdentity(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    artifact_set_sha256: str = Field(pattern=_SHA256)
    model_version: str = Field(min_length=1)
    model_sha256: str = Field(pattern=_SHA256)
    code_version_or_commit: str = Field(min_length=1)
    code_manifest_sha256: str = Field(pattern=_SHA256)
    feature_sha256: str = Field(pattern=_SHA256)
    canonical_feature_list_sha256: str = Field(pattern=_SHA256)
    class_annotation_sha256: str = Field(pattern=_SHA256)
    probe_resource_sha256: str = Field(pattern=_SHA256)
    genome_build: GenomeBuild
    expected_feature_count: Literal[357340] = 357340
    expected_model_unit_count: Literal[42] = 42
    preprocessing_contract_version: Literal["marlin-v1-binarize-1"] = "marlin-v1-binarize-1"
    runtime_contract_version: Literal["marlin-runtime-v1"] = "marlin-runtime-v1"
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def identity_is_marlin_v1_grch37(self) -> MarlinArtifactSetIdentity:
        if self.genome_build is not GenomeBuild.GRCH37:
            raise ValueError("MARLIN v1 artifact-set identity currently requires GRCh37/hg19")
        return self


class MarlinRuntimeIdentity(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    runtime_id: str = Field(min_length=3)
    R_version: str = Field(min_length=1)
    keras_version: str = Field(min_length=1)
    tensorflow_version: str = Field(min_length=1)
    python_version: str = Field(pattern=r"^\d+\.\d+\.\d+$")
    execution_backend: Literal["cpu", "gpu"]
    runtime_contract_version: Literal["marlin-runtime-v1"] = "marlin-runtime-v1"
    created_at: datetime
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def identity_is_normalized(self) -> MarlinRuntimeIdentity:
        if self.created_at.utcoffset() is None:
            raise ValueError("MARLIN runtime-identity timestamp requires a timezone")
        text_fields = (
            self.runtime_id,
            self.R_version,
            self.keras_version,
            self.tensorflow_version,
            self.python_version,
        )
        if any(value != value.strip() for value in text_fields):
            raise ValueError("MARLIN runtime-identity text fields must be whitespace-normalized")
        return self


class MarlinBridgeLock(StrictModel):
    """Evidence lock required before the native modkit-to-MARLIN bridge can be used.

    The hashes under ``validation_*`` identify the same-specimen comparison that justified
    enabling this adapter. They are evidence provenance, not identities that every future
    specimen must share. Runtime inputs are instead constrained by build, adapter, modkit,
    pileup semantics, probe-resource and feature-resource identities.
    """

    schema_version: Literal["0.1.0"] = "0.1.0"
    bridge_id: str = Field(min_length=3)
    status: Literal["validated_same_specimen_bridge"] = "validated_same_specimen_bridge"
    genome_build: GenomeBuild
    adapter_version: Literal["marlin-modkit-bridge-v1"] = "marlin-modkit-bridge-v1"
    modkit_version: Literal["0.6.4"] = "0.6.4"
    pileup_semantics: Literal["5mC+5hmC-combine-mods-v1"] = "5mC+5hmC-combine-mods-v1"
    probe_resource_sha256: str = Field(pattern=_SHA256)
    feature_artifact_sha256: str = Field(pattern=_SHA256)
    validation_precomputed_input_sha256: str = Field(pattern=_SHA256)
    validation_modkit_bedmethyl_sha256: str = Field(pattern=_SHA256)
    validation_precomputed_feature_vector_sha256: str = Field(pattern=_SHA256)
    validation_modkit_feature_vector_sha256: str = Field(pattern=_SHA256)
    compared_feature_count: int = Field(ge=1, le=357340)
    concordant_feature_count: int = Field(ge=0, le=357340)
    feature_agreement_fraction: float = Field(ge=0, le=1)
    evidence_reference: str = Field(min_length=3)
    validated_at: datetime
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def bridge_evidence_is_consistent(self) -> MarlinBridgeLock:
        if self.genome_build is not GenomeBuild.GRCH37:
            raise ValueError("MARLIN native bridge currently requires GRCh37/hg19")
        if self.validated_at.utcoffset() is None:
            raise ValueError("MARLIN bridge-lock timestamp requires a timezone")
        if self.concordant_feature_count > self.compared_feature_count:
            raise ValueError("MARLIN bridge concordant count cannot exceed compared count")
        expected = self.concordant_feature_count / self.compared_feature_count
        if not math.isclose(self.feature_agreement_fraction, expected, rel_tol=0, abs_tol=1e-12):
            raise ValueError("MARLIN bridge feature agreement fraction is inconsistent with counts")
        return self


class MarlinRuntimeCompatibilityProfile(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    reference_runtime_lock_id: str = Field(min_length=1)
    feature_vector_sha256: str = Field(pattern=_SHA256)
    reference_scores: list[float] = Field(min_length=42, max_length=42)
    absolute_score_tolerance: float = Field(gt=0, le=1)
    score_sum_tolerance: float = Field(gt=0, le=1)
    top_model_unit_index: int = Field(ge=0, le=41)
    execution_backend: str = Field(min_length=1)
    created_at: datetime
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def profile_is_valid_softmax_reference(self) -> MarlinRuntimeCompatibilityProfile:
        if self.created_at.utcoffset() is None:
            raise ValueError("MARLIN runtime-profile timestamp requires a timezone")
        if any(
            not math.isfinite(score) or score < 0 or score > 1 for score in self.reference_scores
        ):
            raise ValueError("MARLIN reference scores must be finite probabilities")
        if abs(sum(self.reference_scores) - 1.0) > self.score_sum_tolerance:
            raise ValueError("MARLIN reference scores violate the declared softmax-sum tolerance")
        observed_top = max(range(42), key=self.reference_scores.__getitem__)
        if observed_top != self.top_model_unit_index:
            raise ValueError("top_model_unit_index does not match the reference score vector")
        return self


class MarlinDualRuntimeCompatibilityReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    comparison_id: str = Field(min_length=3)
    artifact_set_identity: MarlinArtifactSetIdentity
    reference_artifact_lock_id: str = Field(min_length=3)
    candidate_artifact_lock_id: str = Field(min_length=3)
    feature_vector_sha256: str = Field(pattern=_SHA256)
    reference_profile_id: str = Field(min_length=1)
    reference_runtime_identity: MarlinRuntimeIdentity
    candidate_runtime_identity: MarlinRuntimeIdentity
    absolute_score_tolerance: float = Field(gt=0, le=1)
    score_sum_tolerance: float = Field(gt=0, le=1)
    reference_scores: list[float] = Field(min_length=42, max_length=42)
    candidate_scores: list[float] = Field(min_length=42, max_length=42)
    absolute_differences: list[float] = Field(min_length=42, max_length=42)
    max_absolute_difference: float = Field(ge=0)
    reference_top_model_unit_index: int = Field(ge=0, le=41)
    candidate_top_model_unit_index: int = Field(ge=0, le=41)
    all_scores_within_tolerance: bool
    top_model_unit_matches: bool
    softmax_invariants_pass: bool
    verdict: MarlinDualRuntimeVerdict
    created_at: datetime
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def report_is_internally_consistent(self) -> MarlinDualRuntimeCompatibilityReport:
        if self.created_at.utcoffset() is None:
            raise ValueError("MARLIN dual-runtime report timestamp requires a timezone")
        if any(
            not math.isfinite(score) or score < 0 or score > 1
            for score in (*self.reference_scores, *self.candidate_scores)
        ):
            raise ValueError("MARLIN dual-runtime scores must be finite probabilities")
        if any(not math.isfinite(delta) or delta < 0 for delta in self.absolute_differences):
            raise ValueError(
                "MARLIN dual-runtime score differences must be finite and non-negative"
            )
        if not math.isfinite(self.max_absolute_difference):
            raise ValueError("MARLIN dual-runtime max difference must be finite")
        expected_differences = [
            abs(reference - candidate)
            for reference, candidate in zip(
                self.reference_scores, self.candidate_scores, strict=True
            )
        ]
        if any(
            not math.isclose(observed, expected, rel_tol=0, abs_tol=1e-15)
            for observed, expected in zip(
                self.absolute_differences, expected_differences, strict=True
            )
        ):
            raise ValueError("MARLIN dual-runtime score differences are inconsistent")
        expected_max = max(self.absolute_differences)
        if not math.isclose(
            self.max_absolute_difference, expected_max, rel_tol=0, abs_tol=1e-15
        ):
            raise ValueError("MARLIN dual-runtime max difference is inconsistent")
        reference_top = max(range(42), key=self.reference_scores.__getitem__)
        candidate_top = max(range(42), key=self.candidate_scores.__getitem__)
        if reference_top != self.reference_top_model_unit_index:
            raise ValueError("MARLIN dual-runtime reference top model unit is inconsistent")
        if candidate_top != self.candidate_top_model_unit_index:
            raise ValueError("MARLIN dual-runtime candidate top model unit is inconsistent")
        expected_top_match = reference_top == candidate_top
        if self.top_model_unit_matches != expected_top_match:
            raise ValueError("MARLIN dual-runtime top-model-unit match flag is inconsistent")
        expected_within = all(
            delta <= self.absolute_score_tolerance for delta in self.absolute_differences
        )
        if self.all_scores_within_tolerance != expected_within:
            raise ValueError("MARLIN dual-runtime score-tolerance flag is inconsistent")
        expected_softmax = (
            abs(sum(self.reference_scores) - 1.0) <= self.score_sum_tolerance
            and abs(sum(self.candidate_scores) - 1.0) <= self.score_sum_tolerance
        )
        if self.softmax_invariants_pass != expected_softmax:
            raise ValueError("MARLIN dual-runtime softmax flag is inconsistent")
        expected_pass = expected_within and expected_top_match and expected_softmax
        if (self.verdict is MarlinDualRuntimeVerdict.PASS) != expected_pass:
            raise ValueError("MARLIN dual-runtime PASS verdict is inconsistent with gate results")
        return self


class MarlinPredictionReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str = Field(min_length=1)
    status: ModuleRunStatus
    decision: MarlinClassificationDecision
    source_kind: MarlinSourceKind
    genome_build: GenomeBuild
    input_fingerprint: FileFingerprint
    feature_summary: MarlinFeatureSummary
    raw_model_scores: list[MarlinModelUnitScore]
    class_scores: list[MarlinGroupedScore]
    family_scores: list[MarlinGroupedScore]
    lineage_scores: list[MarlinGroupedScore]
    top_class: str | None = Field(default=None, min_length=1)
    top_class_score: float | None = Field(default=None, ge=0, le=1.00001)
    confidence_threshold: float = Field(default=0.8, ge=0.8, le=0.8)
    artifact_lock_id: str = Field(min_length=1)
    runtime_profile_id: str = Field(min_length=1)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def report_is_internally_consistent(self) -> MarlinPredictionReport:
        if self.genome_build is not GenomeBuild.GRCH37:
            raise ValueError("MARLIN v1 prediction reports require GRCh37/hg19")
        if self.status not in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL}:
            raise ValueError("MARLIN normalized report status is COMPLETED or NO_CALL")
        if self.top_class_score is not None and not math.isfinite(self.top_class_score):
            raise ValueError("MARLIN top-class score must be finite")

        if self.status is ModuleRunStatus.NO_CALL:
            if self.decision is not MarlinClassificationDecision.UNKNOWN:
                raise ValueError("MARLIN NO_CALL must carry UNKNOWN decision")
            if self.top_class is not None or self.top_class_score is not None:
                raise ValueError("MARLIN NO_CALL cannot carry a top classification")
            if (
                self.raw_model_scores
                or self.class_scores
                or self.family_scores
                or self.lineage_scores
            ):
                raise ValueError("MARLIN NO_CALL cannot carry model or grouped score payloads")
            return self

        if self.feature_summary.observed_model_feature_count == 0:
            raise ValueError("MARLIN COMPLETED report requires at least one observed model feature")
        if len(self.raw_model_scores) != 42:
            raise ValueError("MARLIN COMPLETED report requires exactly 42 raw model scores")
        ids = [item.model_id for item in self.raw_model_scores]
        labels = [item.label for item in self.raw_model_scores]
        if ids != list(range(1, 43)) or len(set(labels)) != 42:
            raise ValueError(
                "MARLIN raw model scores require model IDs in locked order 1..42 and unique labels"
            )
        if abs(sum(item.score for item in self.raw_model_scores) - 1.0) > 1e-5:
            raise ValueError("MARLIN raw model scores violate softmax-sum tolerance")
        if not self.class_scores:
            raise ValueError("MARLIN COMPLETED report requires grouped current-class scores")
        if self.top_class is None or self.top_class_score is None:
            raise ValueError("MARLIN COMPLETED report requires top current class and score")
        top = min(self.class_scores, key=lambda item: (-item.score, item.label))
        if self.top_class != top.label or not math.isclose(
            self.top_class_score, top.score, rel_tol=0, abs_tol=1e-12
        ):
            raise ValueError("MARLIN top classification must match grouped current-class scores")
        expected_decision = (
            MarlinClassificationDecision.HIGH_CONFIDENCE
            if self.top_class_score >= self.confidence_threshold
            else MarlinClassificationDecision.UNKNOWN
        )
        if self.decision is not expected_decision:
            raise ValueError("MARLIN decision is inconsistent with the fixed confidence threshold")
        return self
