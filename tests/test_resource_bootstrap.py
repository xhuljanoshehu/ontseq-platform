from __future__ import annotations

import argparse
import io
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import MagicMock, patch

import yaml  # type: ignore[import-untyped]

from ontseq_platform.models import PanelBundle
from ontseq_platform.reference_catalog import ReferenceBundleInstaller
from ontseq_platform.resource_bootstrap import (
    KNOWLEDGE_BUNDLE_ID,
    PANEL_BUNDLE_ID,
    PROFILE_IDS,
    REFERENCE_BUNDLE_ID,
    GRCh38ResourceBootstrapper,
    ResourceBootstrapError,
)
from ontseq_platform.resource_registry import ResourceRegistry

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs"
REFERENCE_FIXTURE = ROOT / "tests" / "fixtures" / "reference_catalog" / "GRCh38_FIXTURE_v1"


def _install_official_id_reference(resource_root: Path) -> Path:
    """Install the tiny GRCh38 cache fixture under the ID pinned by the curated profiles."""

    source = resource_root / "fixture-source"
    shutil.copytree(REFERENCE_FIXTURE, source)
    recipe_path = source / "bundle.recipe.yaml"
    recipe = yaml.safe_load(recipe_path.read_text(encoding="utf-8"))
    recipe["bundle_id"] = REFERENCE_BUNDLE_ID
    recipe_path.write_text(
        yaml.safe_dump(recipe, sort_keys=False),
        encoding="utf-8",
        newline="\n",
    )
    installed = ReferenceBundleInstaller(resource_root).import_bundle(source)
    shutil.rmtree(source)
    return installed.path


class Grch38ResourceBootstrapTests(unittest.TestCase):
    def setUp(self) -> None:
        """The bootstrap contract uses a real tiny FASTA, not a multi-gigabyte CI genome."""

        self._canonical_patches = (
            patch("ontseq_platform.reference.validate_canonical_reference"),
            patch("ontseq_platform.reference_catalog.validate_canonical_reference"),
        )
        for canonical_patch in self._canonical_patches:
            canonical_patch.start()
            self.addCleanup(canonical_patch.stop)

    def test_curated_resources_compile_activate_and_resolve_both_profiles(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)

            result = GRCh38ResourceBootstrapper(
                resource_root, packaged_config_root=CONFIGS
            ).activate()

            self.assertEqual(result.reference_bundle_id, REFERENCE_BUNDLE_ID)
            self.assertEqual(result.panel_bundle_id, PANEL_BUNDLE_ID)
            self.assertEqual(result.knowledge_bundle_id, KNOWLEDGE_BUNDLE_ID)
            self.assertEqual(result.profile_ids, PROFILE_IDS)
            self.assertGreaterEqual(result.panel_summary.compilation.resolved_target_count, 2)
            self.assertIn(
                "IGH_REVIEW_REQUIRED",
                result.panel_summary.compilation.unresolved_targets,
            )
            self.assertEqual(len(result.activated_paths), 4)
            self.assertEqual(result.repaired_paths, ())
            self.assertEqual(result.already_active_paths, ())

            installed_source = (
                resource_root
                / "panels"
                / PANEL_BUNDLE_ID
                / "source"
                / "250611_fusion_panel_with_buffer.bed"
            )
            curated_source = (
                CONFIGS
                / "panels"
                / PANEL_BUNDLE_ID
                / "source"
                / "250611_fusion_panel_with_buffer.bed"
            )
            self.assertEqual(installed_source.read_bytes(), curated_source.read_bytes())

            panel_directory = resource_root / "panels" / PANEL_BUNDLE_ID
            panel = PanelBundle.model_validate(
                yaml.safe_load((panel_directory / "bundle.yaml").read_text(encoding="utf-8"))
            )
            for resource_id in (
                panel.analysis_roi_resource_id,
                panel.transcript_cache_resource_id,
            ):
                resource = panel.resource(resource_id)
                self.assertIsNotNone(resource.sha256)
                self.assertIsNotNone(resource.size_bytes)
                self.assertTrue((panel_directory / resource.path).is_file())

            registry = ResourceRegistry(resource_root)
            for profile_id in PROFILE_IDS:
                with self.subTest(profile=profile_id):
                    context = registry.resolve_profile(profile_id)
                    self.assertEqual(context.reference_bundle_id, REFERENCE_BUNDLE_ID)
                    self.assertEqual(context.knowledge_bundle_id, KNOWLEDGE_BUNDLE_ID)
            self.assertIsNone(registry.resolve_profile("AML_LCWGS_GRCh38").panel_bundle_id)
            self.assertEqual(
                registry.resolve_profile("AML_AS_111_GRCh38").panel_bundle_id,
                PANEL_BUNDLE_ID,
            )

    def test_repeated_activation_is_idempotent_and_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)
            bootstrapper = GRCh38ResourceBootstrapper(resource_root, packaged_config_root=CONFIGS)
            first = bootstrapper.activate()
            manifest = resource_root / "panels" / PANEL_BUNDLE_ID / "bundle.yaml"
            first_manifest = manifest.read_bytes()

            second = bootstrapper.activate()

            self.assertEqual(manifest.read_bytes(), first_manifest)
            self.assertEqual(second.activated_paths, ())
            self.assertEqual(second.repaired_paths, ())
            self.assertEqual(len(second.already_active_paths), 4)
            self.assertEqual(
                second.panel_summary.manifest_sha256,
                first.panel_summary.manifest_sha256,
            )

    def test_missing_reference_never_activates_dependent_resources(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            resource_root = Path(raw).resolve()
            bootstrapper = GRCh38ResourceBootstrapper(resource_root, packaged_config_root=CONFIGS)

            with self.assertRaisesRegex(ResourceBootstrapError, "must be installed and valid"):
                bootstrapper.activate()

            self.assertFalse((resource_root / "knowledge").exists())
            self.assertFalse((resource_root / "panels").exists())
            self.assertFalse((resource_root / "profiles").exists())

    def test_checksum_failure_in_staging_exposes_nothing(self) -> None:
        with (
            tempfile.TemporaryDirectory() as raw,
            tempfile.TemporaryDirectory() as config_raw,
        ):
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)
            configs = Path(config_raw) / "configs"
            shutil.copytree(CONFIGS, configs)
            knowledge = (
                configs / "knowledge_bundles" / KNOWLEDGE_BUNDLE_ID / "aml_rearrangements.v0.1.json"
            )
            knowledge.write_bytes(knowledge.read_bytes() + b"\n")

            with self.assertRaisesRegex(ResourceBootstrapError, "checksum mismatch"):
                GRCh38ResourceBootstrapper(resource_root, packaged_config_root=configs).activate()

            self.assertFalse((resource_root / "knowledge" / KNOWLEDGE_BUNDLE_ID).exists())
            self.assertFalse((resource_root / "panels" / PANEL_BUNDLE_ID).exists())
            self.assertFalse((resource_root / "profiles").exists())
            staging = resource_root / ".bootstrap-staging"
            self.assertFalse(any(staging.iterdir()))

    def test_grch37_profile_is_rejected_without_fallback_or_activation(self) -> None:
        with (
            tempfile.TemporaryDirectory() as raw,
            tempfile.TemporaryDirectory() as config_raw,
        ):
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)
            configs = Path(config_raw) / "configs"
            shutil.copytree(CONFIGS, configs)
            profile_path = configs / "profiles" / "AML_LCWGS_GRCh38.yaml"
            profile = yaml.safe_load(profile_path.read_text(encoding="utf-8"))
            profile["genome_build"] = "GRCh37"
            profile_path.write_text(
                yaml.safe_dump(profile, sort_keys=False),
                encoding="utf-8",
                newline="\n",
            )

            with self.assertRaisesRegex(ResourceBootstrapError, "is not GRCh38"):
                GRCh38ResourceBootstrapper(resource_root, packaged_config_root=configs).activate()

            self.assertFalse((resource_root / "knowledge" / KNOWLEDGE_BUNDLE_ID).exists())
            self.assertFalse((resource_root / "panels" / PANEL_BUNDLE_ID).exists())
            self.assertFalse((resource_root / "profiles").exists())

    def test_different_existing_bundle_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)
            bootstrapper = GRCh38ResourceBootstrapper(resource_root, packaged_config_root=CONFIGS)
            bootstrapper.activate()
            knowledge = (
                resource_root / "knowledge" / KNOWLEDGE_BUNDLE_ID / "aml_rearrangements.v0.1.json"
            )
            knowledge.write_bytes(b"locally-divergent\n")

            with self.assertRaisesRegex(FileExistsError, "different bundle"):
                bootstrapper.activate()

            self.assertEqual(knowledge.read_bytes(), b"locally-divergent\n")

    def test_repair_restores_knowledge_panel_and_profile_without_manual_deletion(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)
            bootstrapper = GRCh38ResourceBootstrapper(resource_root, packaged_config_root=CONFIGS)
            bootstrapper.activate()
            knowledge = (
                resource_root / "knowledge" / KNOWLEDGE_BUNDLE_ID / "aml_rearrangements.v0.1.json"
            )
            panel = (
                resource_root
                / "panels"
                / PANEL_BUNDLE_ID
                / "source"
                / "250611_fusion_panel_with_buffer.bed"
            )
            profile = resource_root / "profiles" / "AML_AS_111_GRCh38.yaml"
            expected_knowledge = knowledge.read_bytes()
            expected_panel = panel.read_bytes()
            expected_profile = profile.read_bytes()
            knowledge.write_bytes(b"corrupt-knowledge\n")
            panel.write_bytes(b"corrupt-panel\n")
            profile.write_bytes(b"corrupt-profile\n")

            result = bootstrapper.repair()

            self.assertEqual(result.activated_paths, ())
            self.assertEqual(
                set(result.repaired_paths),
                {
                    resource_root / "knowledge" / KNOWLEDGE_BUNDLE_ID,
                    resource_root / "panels" / PANEL_BUNDLE_ID,
                    profile,
                },
            )
            self.assertEqual(len(result.already_active_paths), 1)
            self.assertEqual(knowledge.read_bytes(), expected_knowledge)
            self.assertEqual(panel.read_bytes(), expected_panel)
            self.assertEqual(profile.read_bytes(), expected_profile)
            registry = ResourceRegistry(resource_root)
            for profile_id in PROFILE_IDS:
                with self.subTest(profile=profile_id):
                    registry.resolve_profile(profile_id, verify_files=True)

    def test_failed_family_repair_rolls_back_an_earlier_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)
            bootstrapper = GRCh38ResourceBootstrapper(resource_root, packaged_config_root=CONFIGS)
            bootstrapper.activate()
            knowledge = (
                resource_root / "knowledge" / KNOWLEDGE_BUNDLE_ID / "aml_rearrangements.v0.1.json"
            )
            panel_directory = resource_root / "panels" / PANEL_BUNDLE_ID
            panel = panel_directory / "source" / "250611_fusion_panel_with_buffer.bed"
            corrupt_knowledge = b"retain-this-corrupt-knowledge\n"
            corrupt_panel = b"retain-this-corrupt-panel\n"
            knowledge.write_bytes(corrupt_knowledge)
            panel.write_bytes(corrupt_panel)
            real_replace = os.replace

            def fail_panel_publish(source: str | Path, destination: str | Path) -> None:
                if Path(destination) == panel_directory and Path(source).parent.name == "panels":
                    raise OSError("synthetic panel publish failure")
                real_replace(source, destination)

            with (
                patch("ontseq_platform.resource_bootstrap.os.replace", fail_panel_publish),
                self.assertRaisesRegex(OSError, "synthetic panel publish failure"),
            ):
                bootstrapper.repair()

            self.assertEqual(knowledge.read_bytes(), corrupt_knowledge)
            self.assertEqual(panel.read_bytes(), corrupt_panel)

    def test_cli_reference_repair_routes_to_full_profile_family_repair(self) -> None:
        try:
            from ontseq_platform.resource_commands import _bootstrap_if_official
        except ModuleNotFoundError as error:
            if error.name != "ontseq_platform.resource_commands":
                raise
            self.skipTest("the references CLI is introduced by the profile/pipeline stack stage")
        bootstrap_result = MagicMock()
        bootstrap_result.profile_ids = PROFILE_IDS
        bootstrap_result.panel_summary.compilation.roi_interval_count = 111
        bootstrap_result.repaired_paths = (Path("knowledge") / KNOWLEDGE_BUNDLE_ID,)
        bootstrapper = MagicMock()
        bootstrapper.repair.return_value = bootstrap_result
        args = argparse.Namespace(resource_root=Path("resources"), config_root=CONFIGS)

        with (
            patch(
                "ontseq_platform.resource_commands.GRCh38ResourceBootstrapper",
                return_value=bootstrapper,
            ),
            redirect_stdout(io.StringIO()),
        ):
            _bootstrap_if_official(args, REFERENCE_BUNDLE_ID, repair_existing=True)

        bootstrapper.repair.assert_called_once_with()
        bootstrapper.activate.assert_not_called()

    def test_namespace_symlink_is_rejected_without_touching_external_tree(self) -> None:
        with (
            tempfile.TemporaryDirectory() as raw,
            tempfile.TemporaryDirectory() as external_raw,
        ):
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)
            external = Path(external_raw).resolve()
            sentinel = external / "do-not-touch.txt"
            sentinel.write_text("external data\n", encoding="utf-8")
            try:
                (resource_root / "knowledge").symlink_to(
                    external,
                    target_is_directory=True,
                )
            except OSError as exc:
                self.skipTest(f"directory symlinks are unavailable: {exc}")

            bootstrapper = GRCh38ResourceBootstrapper(
                resource_root,
                packaged_config_root=CONFIGS,
            )
            with self.assertRaisesRegex(
                ResourceBootstrapError,
                "symbolic link or junction",
            ):
                bootstrapper.repair()

            self.assertEqual(sentinel.read_text(encoding="utf-8"), "external data\n")
            self.assertFalse((external / KNOWLEDGE_BUNDLE_ID).exists())

    def test_incomplete_rollback_retains_transaction_backups_without_a_marker(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            resource_root = Path(raw).resolve()
            _install_official_id_reference(resource_root)
            bootstrapper = GRCh38ResourceBootstrapper(
                resource_root,
                packaged_config_root=CONFIGS,
            )
            bootstrapper.activate()
            knowledge_directory = resource_root / "knowledge" / KNOWLEDGE_BUNDLE_ID
            panel_directory = resource_root / "panels" / PANEL_BUNDLE_ID
            knowledge = knowledge_directory / "aml_rearrangements.v0.1.json"
            panel = panel_directory / "source" / "250611_fusion_panel_with_buffer.bed"
            knowledge.write_bytes(b"corrupt-knowledge\n")
            panel.write_bytes(b"corrupt-panel\n")
            real_replace = os.replace

            def fail_publish_and_one_rollback(
                source: str | Path,
                destination: str | Path,
            ) -> None:
                source_path = Path(source)
                destination_path = Path(destination)
                if destination_path == panel_directory and source_path.parent.name == "panels":
                    raise OSError("synthetic panel publish failure")
                if destination_path == knowledge_directory and source_path.parent.name == "backups":
                    raise OSError("synthetic rollback restore failure")
                real_replace(source, destination)

            with (
                patch(
                    "ontseq_platform.resource_bootstrap.os.replace",
                    fail_publish_and_one_rollback,
                ),
                self.assertRaisesRegex(
                    ResourceBootstrapError,
                    "rollback was incomplete; recovery data was retained",
                ),
            ):
                bootstrapper.repair()

            transactions = list((resource_root / ".bootstrap-staging").glob("grch38.*"))
            self.assertEqual(len(transactions), 1)
            self.assertTrue((transactions[0] / "backups").is_dir())
            self.assertTrue(any((transactions[0] / "backups").iterdir()))


if __name__ == "__main__":
    unittest.main()
