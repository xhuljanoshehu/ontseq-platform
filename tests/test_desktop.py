from pathlib import Path

import pytest
from pydantic import ValidationError

from ontseq_platform.desktop import (
    DesktopAnalysisRequest,
    DesktopBackend,
    DesktopConfig,
    DesktopReferenceProfile,
    load_desktop_config,
    locate_bam_index,
    sanitize_sample_id,
    save_desktop_config,
)
from ontseq_platform.models import AssayMode, GenomeBuild


def test_sanitize_sample_id_removes_unsafe_characters() -> None:
    assert sanitize_sample_id(" AML Patient / 01 ") == "AML_Patient_01"
    assert sanitize_sample_id("x") == "S_x"
    assert len(sanitize_sample_id("A" * 100)) == 64


def test_locate_bam_index_supports_bam_bai_and_bai(tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM")
    bam_bai = tmp_path / "sample.bam.bai"
    bam_bai.write_bytes(b"BAI")
    assert locate_bam_index(bam) == bam_bai

    bam_bai.unlink()
    short_bai = tmp_path / "sample.bai"
    short_bai.write_bytes(b"BAI")
    assert locate_bam_index(bam) == short_bai


def test_desktop_config_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "desktop.json"
    config = DesktopConfig(
        backend_mode="wsl",
        wsl_distribution="Ubuntu",
        output_root=str(tmp_path / "results"),
        reference_profiles=[
            DesktopReferenceProfile(
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_lock_path=str(tmp_path / "grch38.lock.json"),
            )
        ],
    )
    assert save_desktop_config(config, path) == path
    assert load_desktop_config(path) == config


def test_duplicate_reference_builds_are_rejected(tmp_path: Path) -> None:
    profile = DesktopReferenceProfile(
        genome_build=GenomeBuild.GRCH38,
        reference_id="GRCh38-test",
        reference_lock_path=str(tmp_path / "lock.json"),
    )
    with pytest.raises(ValidationError):
        DesktopConfig(reference_profiles=[profile, profile])


def test_adaptive_sampling_profile_requires_bed_path_and_version_together(tmp_path: Path) -> None:
    with pytest.raises(ValidationError):
        DesktopReferenceProfile(
            genome_build=GenomeBuild.GRCH38,
            reference_id="GRCh38-test",
            reference_lock_path=str(tmp_path / "lock.json"),
            adaptive_sampling_target_bed_path=str(tmp_path / "targets.bed"),
        )


def test_manifest_is_pseudonymized_and_requests_research_modules(tmp_path: Path) -> None:
    bam = tmp_path / "SAMPLE_001.bam"
    bai = tmp_path / "SAMPLE_001.bam.bai"
    bed = tmp_path / "targets.bed"
    bam.write_bytes(b"BAM")
    bai.write_bytes(b"BAI")
    bed.write_text("chr1\t1\t100\n", encoding="utf-8")
    profile = DesktopReferenceProfile(
        genome_build=GenomeBuild.GRCH38,
        reference_id="GRCh38-test",
        reference_lock_path=str(tmp_path / "lock.json"),
        adaptive_sampling_target_bed_path=str(bed),
        adaptive_sampling_target_bed_version="test-v1",
    )
    backend = DesktopBackend(
        DesktopConfig(
            backend_mode="local",
            wsl_project_root=str(tmp_path),
            output_root=str(tmp_path / "results"),
            reference_profiles=[profile],
        )
    )
    request = DesktopAnalysisRequest(
        bam_path=bam,
        sample_id="SAMPLE_001",
        genome_build=GenomeBuild.GRCH38,
        assay_mode=AssayMode.ADAPTIVE_SAMPLING,
        output_dir=tmp_path / "results" / "SAMPLE_001",
        run_id="RUN_TEST_001",
    )

    manifest = backend._build_manifest(  # noqa: SLF001 - explicit unit test of deterministic builder
        request,
        sample_id="SAMPLE_001",
        bam=bam,
        bai=bai,
        profile=profile,
    )

    assert manifest["sample_id"] == "SAMPLE_001"
    assert manifest["privacy"] == {
        "pseudonymized": True,
        "contains_direct_identifiers": False,
        "cloud_upload_approved": False,
    }
    analysis = manifest["analysis"]
    assert isinstance(analysis, dict)
    assert analysis["modules"] == ["qc", "cnv", "sv", "fusion", "iscn", "report"]
    assay = manifest["assay"]
    assert isinstance(assay, dict)
    assert assay["target_bed_version"] == "test-v1"


def test_project_root_expression_supports_home_shortcut() -> None:
    assert DesktopBackend._project_root_expression("~") == '"$HOME"'  # noqa: SLF001
    assert DesktopBackend._project_root_expression("~/ontseq-platform").startswith('"$HOME"/')  # noqa: SLF001
