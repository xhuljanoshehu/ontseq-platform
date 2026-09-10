"""Adversarial regressions: descriptive coverage must not invent observability."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
from pydantic import ValidationError

from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import FileFingerprint, GenomeBuild, ModuleRunStatus, ToolRecord
from ontseq_platform.report import render_html
from ontseq_platform.target_coverage import (
    TargetCoveragePolicy,
    TargetCoverageRegion,
    TargetCoverageReport,
    _parse_float,
)
from ontseq_platform.workbook import render_workbook


def region_payload() -> dict:
    return dict(
        chromosome="chr1",
        start=0,
        end=100,
        region_id="SYNTHETIC_ROI",
        mean_depth=12.5,
        bases_at_threshold={"1x": 100, "10x": 50, "20x": 0, "30x": 0},
        fraction_at_threshold={"1x": 1.0, "10x": 0.5, "20x": 0.0, "30x": 0.0},
    )


def coverage() -> TargetCoverageReport:
    return TargetCoverageReport(
        sample_id=build_demo_result().manifest.sample_id,
        genome_build=GenomeBuild.GRCH38,
        target_bed_version="synthetic-v1",
        status=ModuleRunStatus.COMPLETED,
        policy=TargetCoveragePolicy(
            profile_id="synthetic", status="technical_defaults_only", note="test"
        ),
        summary_metrics={"region_count": 1, "interval_bases": 100},
        regions=[TargetCoverageRegion(**region_payload())],
        target_bed_fingerprint=FileFingerprint(size_bytes=1, sha256="0" * 64),
        tool=ToolRecord(name="mosdepth", version="0.3.14"),
    )


@pytest.mark.parametrize("depth", [math.inf, -math.inf, math.nan])
def test_nonfinite_depth_is_not_valid_coverage(depth: float) -> None:
    payload = region_payload()
    payload["mean_depth"] = depth
    with pytest.raises(ValidationError):
        TargetCoverageRegion(**payload)


@pytest.mark.parametrize("text", ["nan", "inf", "-inf", "1e999"])
def test_parser_refuses_nonfinite_depth(text: str) -> None:
    with pytest.raises(ValueError):
        _parse_float(text, field="mean depth")


def test_fraction_must_represent_the_recorded_base_count() -> None:
    payload = region_payload()
    payload["fraction_at_threshold"]["10x"] = 0.9
    with pytest.raises(ValidationError, match="fraction"):
        TargetCoverageRegion(**payload)


def test_higher_depth_threshold_cannot_have_more_bases() -> None:
    payload = region_payload()
    payload["bases_at_threshold"]["20x"] = 75
    payload["fraction_at_threshold"]["20x"] = 0.75
    with pytest.raises(ValidationError, match="non-increasing|monotonic"):
        TargetCoverageRegion(**payload)


@pytest.mark.parametrize(
    "key,value",
    [
        ("interval_weighted_mean_depth", 123.0),
        ("minimum_region_mean_depth", 123.0),
        ("median_region_mean_depth", 123.0),
        ("maximum_region_mean_depth", 123.0),
        ("interval_bases_at_10x_fraction", 0.9),
        ("interval_weighted_mean_depth", math.inf),
    ],
)
def test_supplied_summary_must_agree_with_regions(key: str, value: float) -> None:
    payload = coverage().model_dump()
    payload["summary_metrics"][key] = value
    with pytest.raises(ValidationError):
        TargetCoverageReport.model_validate(payload)


def test_missing_optional_summaries_stay_missing() -> None:
    assert "interval_weighted_mean_depth" not in coverage().summary_metrics


@pytest.mark.parametrize("renderer,suffix", [(render_html, ".html"), (render_workbook, ".xlsx")])
@pytest.mark.parametrize(
    "field,value", [("sample_id", "ANOTHER_SAMPLE"), ("genome_build", GenomeBuild.GRCH37)]
)
def test_report_rejects_foreign_coverage_before_writing(
    tmp_path: Path, renderer, suffix, field, value
) -> None:
    sidecar = coverage().model_copy(update={field: value})
    output = tmp_path / ("report" + suffix)
    with pytest.raises(ValueError, match="sample|build"):
        renderer(build_demo_result(), output, target_coverage=sidecar)
    assert not output.exists()
