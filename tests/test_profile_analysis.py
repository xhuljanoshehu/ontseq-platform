from __future__ import annotations

import argparse
import io
import json
import shutil
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest
import yaml

from ontseq_platform.io import load_model
from ontseq_platform.models import (
    GenomeBuild,
    QCPolicy,
    ReferenceContig,
    ReferenceDictionaryContract,
    ReferenceLock,
    SnifflesPolicy,
    TargetBedRole,
)
from ontseq_platform.pipeline.components import RunComponents
from ontseq_platform.profile_analysis import (
    AnalyzeSettings,
    ProfileRuntimeSettings,
    build_profile_run_configuration,
)
from ontseq_platform.reference import canonical_contigs, grch38_canonical_25_contigs, sha256_file
from ontseq_platform.reference_catalog import ReferenceBundleInstaller
from ontseq_platform.resource_bootstrap import (
    REFERENCE_BUNDLE_ID,
    GRCh38ResourceBootstrapper,
)
from ontseq_platform.resource_commands import handle_references_command
from ontseq_platform.resource_registry import ResourceRegistry
from ontseq_platform.service.app import ServiceConfig, _resolvable_profile_ids
from ontseq_platform.target_coverage import TargetCoveragePolicy

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs"
REFERENCE_FIXTURE = ROOT / "tests" / "fixtures" / "reference_catalog" / "GRCh38_FIXTURE_v1"


def _header(build: GenomeBuild) -> str:
    return _header_from_contigs(canonical_contigs(build))


def _header_from_contigs(contigs: tuple[tuple[str, int], ...]) -> str:
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines.extend(f"@SQ\tSN:{name}\tLN:{length}" for name, length in contigs)
    return "\n".join(lines) + "\n"


def _full_fixture_contigs() -> tuple[tuple[str, int], ...]:
    return (*grch38_canonical_25_contigs(), ("GL000008.2", 209709))


def _resource_root(root: Path) -> Path:
    source = root / "fixture-source"
    shutil.copytree(REFERENCE_FIXTURE, source)
    recipe_path = source / "bundle.recipe.yaml"
    recipe = yaml.safe_load(recipe_path.read_text(encoding="utf-8"))
    recipe["bundle_id"] = REFERENCE_BUNDLE_ID
    recipe_path.write_text(yaml.safe_dump(recipe, sort_keys=False), encoding="utf-8")

    # Transaction/profile tests stay byte-small. Canonical dictionary enforcement itself is
    # covered in test_bam_resolution.py and test_reference_canonical.py; there is no runtime flag
    # that can bypass it in the product.
    def fixture_reference_lock(fai_path: Path, **kwargs: object) -> ReferenceLock:
        contigs = _full_fixture_contigs()
        return ReferenceLock(
            reference_id=str(kwargs["reference_id"]),
            genome_build=GenomeBuild.GRCH38,
            contigs=[ReferenceContig(name=name, length=length) for name, length in contigs],
            allow_extra_contigs=False,
            source_fai_sha256=sha256_file(fai_path),
        )

    with (
        patch("ontseq_platform.reference.validate_canonical_reference"),
        patch("ontseq_platform.reference_catalog.validate_canonical_reference"),
        patch(
            "ontseq_platform.reference_catalog.reference_lock_from_fai",
            side_effect=fixture_reference_lock,
        ),
    ):
        ReferenceBundleInstaller(root).import_bundle(source)
        shutil.rmtree(source)
        GRCh38ResourceBootstrapper(root, packaged_config_root=CONFIGS).activate()
    return root


def _bam(root: Path) -> Path:
    bam = root / "AML sample.bam"
    bam.write_bytes(b"synthetic BAM placeholder")
    Path(f"{bam}.bai").write_bytes(b"synthetic BAI placeholder")
    return bam


@pytest.mark.parametrize(
    ("profile_id", "adaptive", "dictionary_contract"),
    (
        ("AML_LCWGS_GRCh38", False, ReferenceDictionaryContract.EXACT_FULL),
        ("AML_AS_111_GRCh38", True, ReferenceDictionaryContract.EXACT_FULL),
        (
            "AML_LCWGS_GRCh38_CANONICAL25",
            False,
            ReferenceDictionaryContract.GRCH38_CANONICAL_25,
        ),
        (
            "AML_AS_111_GRCh38_CANONICAL25",
            True,
            ReferenceDictionaryContract.GRCH38_CANONICAL_25,
        ),
    ),
)
def test_profile_configuration_resolves_every_resource_without_manual_paths(
    profile_id: str,
    adaptive: bool,
    dictionary_contract: ReferenceDictionaryContract,
) -> None:
    with TemporaryDirectory() as raw:
        root = Path(raw).resolve()
        resources = _resource_root(root / "resources")
        bam = _bam(root)

        header_contigs = (
            grch38_canonical_25_contigs()
            if dictionary_contract == ReferenceDictionaryContract.GRCH38_CANONICAL_25
            else _full_fixture_contigs()
        )
        config = build_profile_run_configuration(
            AnalyzeSettings(
                bam=bam,
                profile_id=profile_id,
                resource_root=resources,
                output_dir=root / "results",
                configuration_root=CONFIGS,
                pipeline_version="test",
            ),
            header_text=_header_from_contigs(header_contigs),
        )

        assert config.manifest.assay.genome_build == GenomeBuild.GRCH38
        assert config.manifest.analysis.profile == profile_id
        assert config.manifest.run_id.startswith("AML_sample-")
        assert config.manifest.run_id.endswith("Z")
        assert config.reference_fasta is not None
        assert config.annotation_cache is not None
        assert config.resource_context is not None
        assert config.resource_context.reference_bundle_id == REFERENCE_BUNDLE_ID
        assert config.resource_context.reference_dictionary_contract == dictionary_contract
        assert tuple((item.name, item.length) for item in config.reference_lock.contigs) == (
            header_contigs
        )
        assert config.reference_lock.allow_extra_contigs is False
        if adaptive:
            assert config.selection_target_bed is not None
            assert config.manifest.assay.target_bed is not None
            assert Path(config.manifest.assay.target_bed) != config.selection_target_bed
            assert config.manifest.assay.target_bed_role == TargetBedRole.ANALYSIS_ROI_UNBUFFERED
        else:
            assert config.selection_target_bed is None
            assert config.manifest.assay.target_bed is None
            assert config.target_coverage_policy is None


def test_full_profile_does_not_fall_back_to_canonical25_header() -> None:
    with TemporaryDirectory() as raw:
        root = Path(raw).resolve()
        resources = _resource_root(root / "resources")
        bam = _bam(root)

        with pytest.raises(ValueError, match="1 missing"):
            build_profile_run_configuration(
                AnalyzeSettings(
                    bam=bam,
                    profile_id="AML_LCWGS_GRCh38",
                    resource_root=resources,
                    output_dir=root / "results",
                    configuration_root=CONFIGS,
                ),
                header_text=_header_from_contigs(grch38_canonical_25_contigs()),
            )


def test_profile_configuration_rejects_grch37_header_before_pipeline() -> None:
    with TemporaryDirectory() as raw:
        root = Path(raw).resolve()
        resources = _resource_root(root / "resources")
        bam = _bam(root)

        with pytest.raises(ValueError, match="profile requires GRCh38"):
            build_profile_run_configuration(
                AnalyzeSettings(
                    bam=bam,
                    profile_id="AML_LCWGS_GRCh38",
                    resource_root=resources,
                    output_dir=root / "results",
                    configuration_root=CONFIGS,
                ),
                header_text=_header(GenomeBuild.GRCH37),
            )


def test_profile_configuration_can_use_fast_pinned_resource_preflight() -> None:
    with TemporaryDirectory() as raw:
        root = Path(raw).resolve()
        resources = _resource_root(root / "resources")
        bam = _bam(root)

        with (
            patch(
                "ontseq_platform.resource_registry.sha256_file",
                side_effect=AssertionError("interactive profile preflight must not compute SHA256"),
            ),
            patch("ontseq_platform.bam_resolution.validate_full_dictionary"),
        ):
            config = build_profile_run_configuration(
                AnalyzeSettings(
                    bam=bam,
                    profile_id="AML_AS_111_GRCh38",
                    resource_root=resources,
                    output_dir=root / "results",
                    configuration_root=CONFIGS,
                    verify_resource_checksums=False,
                ),
                header_text=_header(GenomeBuild.GRCH38),
            )

        assert config.resource_context is not None
        assert config.resource_context.panel_bundle_id == "AML_AS_111_GRCh38_v1"


def test_profile_configuration_preserves_trusted_runtime_overrides() -> None:
    with TemporaryDirectory() as raw:
        root = Path(raw).resolve()
        resources = _resource_root(root / "resources")
        bam = _bam(root)
        qc_policy = load_model(CONFIGS / "qc" / "defaults.yaml", QCPolicy).model_copy(
            update={"note": "service-selected QC policy"}
        )
        sniffles_policy = load_model(
            CONFIGS / "sv" / "sniffles2.conservative.technical.yaml", SnifflesPolicy
        ).model_copy(update={"profile_id": "service-selected-sniffles"})
        target_policy = load_model(
            CONFIGS / "qc" / "adaptive_target_coverage.technical.yaml",
            TargetCoveragePolicy,
        ).model_copy(update={"profile_id": "service-selected-coverage"})
        components = RunComponents.model_validate(
            {
                "selection_id": "service-components",
                "status": "technical_defaults_only",
                "components": {"cnv": {"provider": "qdnaseq_ace", "enabled": False}},
            }
        )
        runtime = ProfileRuntimeSettings(
            qc_policy=qc_policy,
            sniffles_policy=sniffles_policy,
            cutesv_policy=None,
            sv_consensus_policy=None,
            sv_evidence_policy=None,
            target_coverage_policy=target_policy,
            sv_minimum_mean_depth=27.5,
            components=components,
        )

        with patch("ontseq_platform.bam_resolution.validate_full_dictionary"):
            as_config = build_profile_run_configuration(
                AnalyzeSettings(
                    bam=bam,
                    profile_id="AML_AS_111_GRCh38",
                    resource_root=resources,
                    output_dir=root / "results",
                    configuration_root=CONFIGS,
                    runtime_settings=runtime,
                ),
                header_text=_header(GenomeBuild.GRCH38),
            )
            lcwgs_config = build_profile_run_configuration(
                AnalyzeSettings(
                    bam=bam,
                    profile_id="AML_LCWGS_GRCh38",
                    resource_root=resources,
                    output_dir=root / "results",
                    configuration_root=CONFIGS,
                    runtime_settings=runtime,
                ),
                header_text=_header(GenomeBuild.GRCH38),
            )

        assert as_config.qc_policy.note == "service-selected QC policy"
        assert as_config.sniffles_policy is not None
        assert as_config.sniffles_policy.profile_id == "service-selected-sniffles"
        assert as_config.target_coverage_policy is not None
        assert as_config.target_coverage_policy.profile_id == "service-selected-coverage"
        assert as_config.cutesv_policy is None
        assert as_config.sv_consensus_policy is None
        assert as_config.sv_evidence_policy is None
        assert as_config.sv_minimum_mean_depth == 27.5
        assert as_config.components is components
        assert lcwgs_config.target_coverage_policy is None


def test_service_advertises_only_profiles_with_a_resolvable_pinned_context() -> None:
    with TemporaryDirectory() as raw:
        root = Path(raw).resolve()
        resources = _resource_root(root / "resources")
        config = ServiceConfig(
            reference_lock=None,
            output_dir=root / "results",
            allowed_roots=[root],
            qc_policy=root / "qc.yaml",
            sniffles_policy=root / "sniffles.yaml",
            target_coverage_policy=root / "coverage.yaml",
            resource_root=resources,
        )

        assert _resolvable_profile_ids(config) == [
            "AML_LCWGS_GRCh38",
            "AML_AS_111_GRCh38",
            "AML_LCWGS_GRCh38_CANONICAL25",
            "AML_AS_111_GRCh38_CANONICAL25",
        ]

        registry = ResourceRegistry(resources)
        context = registry.resolve_profile("AML_AS_111_GRCh38", verify_files=False)
        selection = Path(context.resource_paths["panel.selection_panel_buffered"])
        selection.write_bytes(selection.read_bytes() + b"corrupt-size")

        assert _resolvable_profile_ids(config) == [
            "AML_LCWGS_GRCh38",
            "AML_LCWGS_GRCh38_CANONICAL25",
        ]


def test_reference_status_reports_only_profiles_with_ready_panel_and_knowledge() -> None:
    with TemporaryDirectory() as raw:
        root = Path(raw).resolve()
        resources = _resource_root(root / "resources")

        def status_payload() -> dict[str, object]:
            output = io.StringIO()
            args = argparse.Namespace(
                command="references",
                references_command="status",
                resource_root=resources,
                config_root=CONFIGS,
                as_json=True,
            )
            with redirect_stdout(output):
                assert handle_references_command(args)
            payload = json.loads(output.getvalue())
            assert isinstance(payload, dict)
            return payload

        ready = status_payload()
        assert ready["profiles"] == [
            "AML_AS_111_GRCh38",
            "AML_AS_111_GRCh38_CANONICAL25",
            "AML_LCWGS_GRCh38",
            "AML_LCWGS_GRCh38_CANONICAL25",
        ]
        statuses = {item["profile_id"]: item for item in ready["profile_status"]}
        assert statuses["AML_AS_111_GRCh38"]["panel_bundle"] == "AML_AS_111_GRCh38_v1"
        assert statuses["AML_AS_111_GRCh38"]["knowledge_bundle"] == "HEMATOLOGY_v3"

        registry = ResourceRegistry(resources)
        context = registry.resolve_profile("AML_LCWGS_GRCh38", verify_files=False)
        knowledge = Path(context.resource_paths["knowledge.rearrangement_knowledge"])
        knowledge.write_bytes(knowledge.read_bytes() + b"corrupt-size")

        broken = status_payload()
        assert broken["profiles"] == []
        assert all(item["valid"] is False for item in broken["profile_status"])
