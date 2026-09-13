from __future__ import annotations

from pathlib import Path

import pytest
from openpyxl import Workbook

from ontseq_platform.marlin_classification import (
    MarlinClassAnnotation,
    classify_marlin_scores,
    load_marlin_class_annotations,
)
from ontseq_platform.marlin_contracts import (
    MarlinClassificationDecision,
    MarlinFeatureSummary,
    MarlinModelUnitScore,
    MarlinSourceKind,
)
from ontseq_platform.marlin_runtime import MarlinRuntimeResult
from ontseq_platform.models import FileFingerprint, GenomeBuild, ModuleRunStatus


def _annotations() -> tuple[MarlinClassAnnotation, ...]:
    rows: list[MarlinClassAnnotation] = []
    for model_id in range(1, 43):
        if model_id in {1, 2}:
            cls, family, lineage = "Class A", "Family A", "Myeloid"
        elif model_id == 3:
            cls, family, lineage = "Class B", "Family B", "Lymphoid"
        else:
            cls, family, lineage = f"Class {model_id}", f"Family {model_id}", "Other"
        rows.append(
            MarlinClassAnnotation(
                model_id=model_id,
                class_name_current=cls,
                methylation_class_family=family,
                lineage=lineage,
            )
        )
    return tuple(rows)


def _feature_summary(observed: int = 1000) -> MarlinFeatureSummary:
    return MarlinFeatureSummary(
        observed_model_feature_count=observed,
        explicit_na_feature_count=0,
        absent_feature_count=357340 - observed,
        non_model_probe_count=0,
        observed_fraction=observed / 357340,
        feature_vector_sha256="a" * 64,
        feature_artifact_sha256="b" * 64,
    )


def _runtime(scores: tuple[float, ...]) -> MarlinRuntimeResult:
    return MarlinRuntimeResult(
        artifact_lock_id="LOCK",
        feature_vector_sha256="a" * 64,
        model_scores=[
            MarlinModelUnitScore(
                model_id=index + 1,
                label=f"model-unit-{index + 1}",
                score=score,
            )
            for index, score in enumerate(scores)
        ],
        raw_score_sha256="c" * 64,
        execution_backend="cpu",
    )


def _classify(runtime: MarlinRuntimeResult | None, observed: int = 1000):
    return classify_marlin_scores(
        sample_id="SAMPLE_001",
        runtime_result=runtime,
        feature_summary=_feature_summary(observed),
        input_fingerprint=FileFingerprint(size_bytes=10, sha256="d" * 64),
        source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
        genome_build=GenomeBuild.GRCH37,
        artifact_lock_id="LOCK",
        runtime_profile_id="PROFILE",
        annotations=_annotations(),
    )


def test_valid_low_grouped_score_is_completed_unknown() -> None:
    scores = tuple([0.40, 0.29, 0.31] + [0.0] * 39)
    report = _classify(_runtime(scores))
    assert report.status is ModuleRunStatus.COMPLETED
    assert report.decision is MarlinClassificationDecision.UNKNOWN
    assert report.top_class == "Class A"
    assert report.top_class_score == pytest.approx(0.69)
    assert report.runtime_profile_id == "PROFILE"


def test_grouped_score_at_exact_threshold_is_high_confidence() -> None:
    scores = tuple([0.50, 0.30, 0.20] + [0.0] * 39)
    report = _classify(_runtime(scores))
    assert report.top_class_score == pytest.approx(0.8)
    assert report.decision is MarlinClassificationDecision.HIGH_CONFIDENCE


def test_zero_observed_features_is_no_call_without_runtime() -> None:
    report = _classify(None, observed=0)
    assert report.status is ModuleRunStatus.NO_CALL
    assert report.decision is MarlinClassificationDecision.UNKNOWN
    assert report.raw_model_scores == []
    assert report.top_class is None
    assert report.runtime_profile_id == "PROFILE"


def test_nonzero_observed_features_require_runtime_result() -> None:
    with pytest.raises(ValueError, match="runtime result"):
        _classify(None, observed=1)


def test_runtime_feature_digest_and_lock_must_match_classification_inputs() -> None:
    runtime = _runtime(tuple([1.0] + [0.0] * 41))
    with pytest.raises(ValueError, match="feature vector"):
        classify_marlin_scores(
            sample_id="SAMPLE_001",
            runtime_result=runtime.model_copy(update={"feature_vector_sha256": "e" * 64}),
            feature_summary=_feature_summary(),
            input_fingerprint=FileFingerprint(size_bytes=10, sha256="d" * 64),
            source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
            genome_build=GenomeBuild.GRCH37,
            artifact_lock_id="LOCK",
            runtime_profile_id="PROFILE",
            annotations=_annotations(),
        )
    with pytest.raises(ValueError, match="artifact lock"):
        classify_marlin_scores(
            sample_id="SAMPLE_001",
            runtime_result=runtime.model_copy(update={"artifact_lock_id": "OTHER"}),
            feature_summary=_feature_summary(),
            input_fingerprint=FileFingerprint(size_bytes=10, sha256="d" * 64),
            source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
            genome_build=GenomeBuild.GRCH37,
            artifact_lock_id="LOCK",
            runtime_profile_id="PROFILE",
            annotations=_annotations(),
        )


def test_raw_scores_are_preserved_and_grouping_is_deterministic() -> None:
    runtime = _runtime(tuple([0.4, 0.3, 0.2, 0.1] + [0.0] * 38))
    report = _classify(runtime)
    assert report.raw_model_scores == runtime.model_scores
    class_scores = {item.label: item.score for item in report.class_scores}
    assert class_scores["Class A"] == pytest.approx(0.7)
    family_scores = {item.label: item.score for item in report.family_scores}
    assert family_scores["Family A"] == pytest.approx(0.7)
    lineage_scores = {item.label: item.score for item in report.lineage_scores}
    assert lineage_scores["Myeloid"] == pytest.approx(0.7)
    assert lineage_scores["Lymphoid"] == pytest.approx(0.2)
    assert lineage_scores["Other"] == pytest.approx(0.1)


def _write_annotation_xlsx(path: Path, rows: list[tuple[object, ...]]) -> None:
    workbook = Workbook()
    sheet = workbook.active
    sheet.append(["model_id", "class_name_current", "mcf", "lineage"])
    for row in rows:
        sheet.append(row)
    workbook.save(path)


def test_annotation_loader_reads_exact_42_model_ids(tmp_path: Path) -> None:
    path = tmp_path / "annotations.xlsx"
    rows = [(i, f"Class {i}", f"Family {i}", "Lineage") for i in range(1, 43)]
    _write_annotation_xlsx(path, rows)
    annotations = load_marlin_class_annotations(path)
    assert [item.model_id for item in annotations] == list(range(1, 43))
    assert annotations[0].class_name_current == "Class 1"


def test_annotation_loader_rejects_missing_columns_or_duplicate_ids(tmp_path: Path) -> None:
    missing = tmp_path / "missing.xlsx"
    workbook = Workbook()
    sheet = workbook.active
    sheet.append(["model_id", "class_name_current", "lineage"])
    workbook.save(missing)
    with pytest.raises(ValueError, match="mcf"):
        load_marlin_class_annotations(missing)

    duplicate = tmp_path / "duplicate.xlsx"
    rows = [(i, f"Class {i}", f"Family {i}", "Lineage") for i in range(1, 42)]
    rows.append((41, "Duplicate", "Family", "Lineage"))
    _write_annotation_xlsx(duplicate, rows)
    with pytest.raises(ValueError, match="model IDs"):
        load_marlin_class_annotations(duplicate)


def test_annotation_loader_sorts_rows_by_model_id(tmp_path: Path) -> None:
    path = tmp_path / "unsorted.xlsx"
    rows = [(i, f"Class {i}", f"Family {i}", "Lineage") for i in range(42, 0, -1)]
    _write_annotation_xlsx(path, rows)
    annotations = load_marlin_class_annotations(path)
    assert [item.model_id for item in annotations] == list(range(1, 43))


def test_grouped_scores_preserve_allowed_softmax_roundoff() -> None:
    annotations = tuple(
        MarlinClassAnnotation(
            model_id=i,
            class_name_current="Class A",
            methylation_class_family="Family A",
            lineage="Lineage A",
        )
        for i in range(1, 43)
    )
    scores = tuple([1.000005 / 42] * 42)
    runtime = _runtime(scores)
    report = classify_marlin_scores(
        sample_id="SAMPLE_001",
        runtime_result=runtime,
        feature_summary=_feature_summary(),
        input_fingerprint=FileFingerprint(size_bytes=10, sha256="d" * 64),
        source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
        genome_build=GenomeBuild.GRCH37,
        artifact_lock_id="LOCK",
        runtime_profile_id="PROFILE",
        annotations=annotations,
    )
    assert report.top_class_score == pytest.approx(sum(scores))
    assert report.top_class_score > 1.0
