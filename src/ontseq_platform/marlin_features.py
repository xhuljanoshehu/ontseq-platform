from __future__ import annotations

import hashlib
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from .marlin_contracts import MarlinFeatureSummary, MarlinProbeObservation
from .marlin_input import MarlinPrecomputedInput
from .models import GenomeBuild

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_FEATURE_TO_BYTE = {-1.0: b"\xff", 0.0: b"\x00", 1.0: b"\x01"}
_EXPECTED_FEATURE_COUNT = 357340


@dataclass(frozen=True)
class _FeatureSummaryDraft:
    expected_feature_count: int
    observed_model_feature_count: int
    explicit_na_feature_count: int
    absent_feature_count: int
    non_model_probe_count: int
    observed_fraction: float
    feature_vector_sha256: str
    feature_artifact_sha256: str


@dataclass(frozen=True)
class _FeatureVectorDraft:
    values: tuple[float, ...]
    summary: _FeatureSummaryDraft


@dataclass(frozen=True)
class MarlinFeatureVector:
    values: tuple[float, ...]
    summary: MarlinFeatureSummary


def _validate_feature_ids(feature_ids: Sequence[str]) -> tuple[str, ...]:
    normalized = tuple(feature_ids)
    if not normalized:
        raise ValueError("MARLIN feature list cannot be empty")
    if any(not feature_id or feature_id != feature_id.strip() for feature_id in normalized):
        raise ValueError("MARLIN feature identifiers must be non-empty and whitespace-free")
    if len(set(normalized)) != len(normalized):
        raise ValueError("MARLIN feature list contains duplicate identifiers")
    return normalized


def _build_feature_vector_for_ids(
    feature_ids: Sequence[str],
    observations: Mapping[str, MarlinProbeObservation],
    *,
    feature_artifact_sha256: str,
) -> _FeatureVectorDraft:
    """Build the published -1/0/+1 representation for an explicit ordered feature list."""
    ordered = _validate_feature_ids(feature_ids)
    if _SHA256.fullmatch(feature_artifact_sha256) is None:
        raise ValueError("feature_artifact_sha256 must be lowercase SHA-256")
    for key, observation in observations.items():
        if key != observation.probe_id:
            raise ValueError("MARLIN observation mapping key must equal observation.probe_id")

    feature_set = set(ordered)
    values: list[float] = []
    observed_count = 0
    explicit_na_count = 0
    absent_count = 0

    for feature_id in ordered:
        observation = observations.get(feature_id)
        if observation is None:
            values.append(0.0)
            absent_count += 1
            continue
        if observation.methylation_fraction is None:
            values.append(0.0)
            explicit_na_count += 1
            continue
        values.append(1.0 if observation.methylation_fraction >= 0.5 else -1.0)
        observed_count += 1

    non_model_count = sum(1 for probe_id in observations if probe_id not in feature_set)
    payload = b"".join(_FEATURE_TO_BYTE[value] for value in values)
    digest = hashlib.sha256(payload).hexdigest()
    expected_count = len(ordered)
    summary = _FeatureSummaryDraft(
        expected_feature_count=expected_count,
        observed_model_feature_count=observed_count,
        explicit_na_feature_count=explicit_na_count,
        absent_feature_count=absent_count,
        non_model_probe_count=non_model_count,
        observed_fraction=observed_count / expected_count,
        feature_vector_sha256=digest,
        feature_artifact_sha256=feature_artifact_sha256,
    )
    return _FeatureVectorDraft(values=tuple(values), summary=summary)


def build_marlin_feature_vector(
    source: MarlinPrecomputedInput,
    feature_ids: Sequence[str],
    *,
    feature_artifact_sha256: str,
) -> MarlinFeatureVector:
    """Build the production MARLIN v1 vector under the fixed 357,340-feature contract."""
    if source.genome_build is not GenomeBuild.GRCH37:
        raise ValueError("MARLIN v1 feature construction currently requires GRCh37/hg19")
    ordered = _validate_feature_ids(feature_ids)
    if len(ordered) != _EXPECTED_FEATURE_COUNT:
        raise ValueError(
            f"MARLIN v1 feature construction requires exactly {_EXPECTED_FEATURE_COUNT} features; "
            f"found {len(ordered)}"
        )
    observations = {item.probe_id: item for item in source.observations}
    if len(observations) != len(source.observations):
        raise ValueError("MARLIN input contains duplicate probe identifiers")
    draft = _build_feature_vector_for_ids(
        ordered,
        observations,
        feature_artifact_sha256=feature_artifact_sha256,
    )
    summary = MarlinFeatureSummary(
        observed_model_feature_count=draft.summary.observed_model_feature_count,
        explicit_na_feature_count=draft.summary.explicit_na_feature_count,
        absent_feature_count=draft.summary.absent_feature_count,
        non_model_probe_count=draft.summary.non_model_probe_count,
        observed_fraction=draft.summary.observed_fraction,
        feature_vector_sha256=draft.summary.feature_vector_sha256,
        feature_artifact_sha256=draft.summary.feature_artifact_sha256,
    )
    return MarlinFeatureVector(values=draft.values, summary=summary)
