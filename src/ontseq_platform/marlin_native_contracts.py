"""Independent native research contracts; these confer no legacy bridge validation."""

from __future__ import annotations

import math
from typing import Literal

from pydantic import Field, model_validator

from .marlin_contracts import (
    MarlinClassificationDecision,
    MarlinFeatureSummary,
    MarlinGroupedScore,
    MarlinModelUnitScore,
)
from .models import FileFingerprint, GenomeBuild, ModuleRunStatus, StrictModel, ToolRecord

ADAPTER_VERSION = "marlin-native-research-v1"
MODEL_SHA256 = "6210a674a0e7690b7b6c03184f5732a65037cf00ff4383e1892f26e90fdcb217"
SHA256_PATTERN = r"^[0-9a-f]{64}$"


class NativeMarlinAsset(StrictModel):
    path: str = Field(min_length=1)
    sha256: str = Field(pattern=SHA256_PATTERN)
    source_uri: str = Field(min_length=1)


class NativeMarlinProbeMap(NativeMarlinAsset):
    genome_build: GenomeBuild


class NativeMarlinInstallation(StrictModel):
    schema_version: Literal["1.0.0"] = "1.0.0"
    adapter_version: Literal["marlin-native-research-v1"] = "marlin-native-research-v1"
    model: NativeMarlinAsset
    features: NativeMarlinAsset
    classes: NativeMarlinAsset
    probe_maps: dict[GenomeBuild, NativeMarlinProbeMap]
    python_executable: str = Field(min_length=1)
    python_version: str = Field(pattern=r"^3\.10\.\d+$")
    tensorflow_version: Literal["2.13.1"] = "2.13.1"
    keras_version: Literal["2.13.1"] = "2.13.1"
    runtime_manifest: NativeMarlinAsset
    runtime_archive_sha256: str = Field(pattern=SHA256_PATTERN)
    confinement: Literal["linux-seccomp-network-deny-v1"] = "linux-seccomp-network-deny-v1"

    @model_validator(mode="after")
    def build_keys_match(self) -> NativeMarlinInstallation:
        if not self.probe_maps or any(k != v.genome_build for k, v in self.probe_maps.items()):
            raise ValueError("probe map build differs from installation key")
        return self


class NativeMarlinReadiness(StrictModel):
    ready: bool
    reason: str = Field(min_length=1)


class NativeMarlinReport(StrictModel):
    schema_version: Literal["1.0.0"] = "1.0.0"
    adapter_version: Literal["marlin-native-research-v1"] = "marlin-native-research-v1"
    run_id: str = Field(min_length=1)
    sample_id: str = Field(min_length=1)
    genome_build: GenomeBuild
    status: ModuleRunStatus
    reason: str = Field(min_length=1)
    decision: MarlinClassificationDecision | None = None
    feature_summary: MarlinFeatureSummary | None = None
    raw_model_scores: list[MarlinModelUnitScore] = Field(default_factory=list)
    class_scores: list[MarlinGroupedScore] = Field(default_factory=list)
    family_scores: list[MarlinGroupedScore] = Field(default_factory=list)
    lineage_scores: list[MarlinGroupedScore] = Field(default_factory=list)
    top_class: str | None = None
    top_class_score: float | None = Field(default=None, ge=0, le=1.00001)
    confidence_threshold: float = Field(default=0.8, ge=0.8, le=0.8)
    validation_status: Literal["UNVALIDATED_RESEARCH"] = "UNVALIDATED_RESEARCH"
    research_only: Literal[True] = True
    tools: list[ToolRecord] = Field(default_factory=list)
    parameters: dict[str, object] = Field(default_factory=dict)
    input_fingerprints: dict[str, FileFingerprint] = Field(default_factory=dict)
    installation_signature: dict[str, str] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(
        default_factory=lambda: [
            "UNVALIDATED_RESEARCH: synthetic technical acceptance is not analytical validation.",
            "MARLIN scores are model ranks, not clinical diagnosis probabilities.",
        ]
    )

    @model_validator(mode="after")
    def evidence_matches_status(self) -> NativeMarlinReport:
        scores = (self.raw_model_scores, self.class_scores, self.family_scores, self.lineage_scores)
        if self.status in {ModuleRunStatus.NOT_RUN, ModuleRunStatus.FAILED}:
            if any(scores) or any(
                v is not None
                for v in (self.decision, self.feature_summary, self.top_class, self.top_class_score)
            ):
                raise ValueError("blocked/failed MARLIN reports cannot carry prediction evidence")
            return self
        if self.feature_summary is None or not self.tools or not self.installation_signature:
            raise ValueError("MARLIN outcome requires feature, tool and installation evidence")
        if not {"bam", "reference", "bedmethyl", "tensor"} <= self.input_fingerprints.keys():
            raise ValueError("MARLIN outcome requires exact input fingerprints")
        if any(v.sha256 is None for v in self.input_fingerprints.values()):
            raise ValueError("MARLIN input fingerprint requires SHA256")
        if self.status == ModuleRunStatus.NO_CALL:
            if self.feature_summary.observed_model_feature_count or any(scores):
                raise ValueError("NO_CALL requires zero observed features and no inference scores")
            if self.decision != MarlinClassificationDecision.UNKNOWN:
                raise ValueError("NO_CALL decision must be UNKNOWN")
            if self.top_class is not None or self.top_class_score is not None:
                raise ValueError("NO_CALL cannot carry a top class")
            return self
        if self.feature_summary.observed_model_feature_count == 0:
            raise ValueError("completed MARLIN requires observed features")
        if [s.model_id for s in self.raw_model_scores] != list(range(1, 43)):
            raise ValueError("MARLIN requires 42 scores in model ID order")
        for group in scores:
            if not group or not math.isclose(math.fsum(s.score for s in group), 1, abs_tol=1e-5):
                raise ValueError("MARLIN score groups must sum to one")
        for group in scores[1:]:
            if len({s.label for s in group}) != len(group):
                raise ValueError("MARLIN grouped labels must be unique")
            if group != sorted(group, key=lambda s: (-s.score, s.label)):
                raise ValueError("MARLIN grouped scores must be sorted")
        top = self.class_scores[0]
        if self.top_class != top.label or self.top_class_score != top.score:
            raise ValueError("MARLIN top class must equal the leading grouped class")
        expected = (
            MarlinClassificationDecision.HIGH_CONFIDENCE
            if top.score >= 0.8
            else MarlinClassificationDecision.UNKNOWN
        )
        if self.decision != expected:
            raise ValueError("MARLIN decision contradicts fixed confidence threshold")
        if "worker_output" not in self.input_fingerprints:
            raise ValueError("completed MARLIN requires worker output fingerprint")
        return self
