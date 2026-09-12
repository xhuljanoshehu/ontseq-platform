"""Explicit profile integration keeps the selected assay, resources and tool pin intact."""

from __future__ import annotations

from pathlib import Path

import pytest
from test_profile_analysis import (
    CONFIGS,
    _bam,
    _full_fixture_contigs,
    _header_from_contigs,
    _resource_root,
)

from ontseq_platform.models import AnalysisModule
from ontseq_platform.profile_analysis import AnalyzeSettings, build_profile_run_configuration


@pytest.mark.parametrize("profile", ["AML_LCWGS_GRCh38", "AML_AS_111_GRCh38"])
@pytest.mark.parametrize("include", [False, True])
def test_profile_methylation_is_opt_in_with_assay_specific_policy(
    tmp_path: Path, profile: str, include: bool
) -> None:
    root = _resource_root(tmp_path / "resources")
    bam = _bam(tmp_path)
    config = build_profile_run_configuration(
        AnalyzeSettings(
            bam=bam,
            profile_id=profile,
            resource_root=root,
            configuration_root=CONFIGS,
            output_dir=tmp_path / "results",
            include_methylation=include,
            executables={"modkit": "trusted-modkit"},
        ),
        header_text=_header_from_contigs(_full_fixture_contigs()),
    )
    assert (AnalysisModule.METHYLATION in config.manifest.analysis.modules) is include
    assert config.manifest.analysis.profile == profile
    assert config.executable("modkit") == "trusted-modkit"
    assert config.reference_fasta is not None
    if include:
        assert config.methylation_policy is not None
        assert config.methylation_policy.expected_version == "0.6.4"
        expected = "target_bed" if "AS_111" in profile else "chromosome"
        assert config.methylation_policy.region_source.value == expected
        assert config.methylation_policy.status == "technical_defaults_only"
    else:
        assert config.methylation_policy is None
