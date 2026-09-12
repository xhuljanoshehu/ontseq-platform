from __future__ import annotations

import argparse
import io
import json
import shutil
import sqlite3
import tempfile
import unittest
from contextlib import closing, redirect_stdout
from pathlib import Path
from unittest.mock import MagicMock, patch

import yaml  # type: ignore[import-untyped]

from ontseq_platform.aml_rearrangements import load_aml_knowledge
from ontseq_platform.models import (
    AmlKnowledgeLock,
    AssayMode,
    GenomeBuild,
    KnowledgeBundle,
    ReferenceDictionaryContract,
)
from ontseq_platform.reference import sha256_file
from ontseq_platform.reference_catalog import ReferenceCatalog
from ontseq_platform.resource_bootstrap import (
    GRCH37_KNOWLEDGE_BUNDLE_ID,
    GRCH37_PANEL_BUNDLE_ID,
    GRCH37_PROFILE_IDS,
    GRCH37_REFERENCE_BUNDLE_ID,
    KNOWLEDGE_BUNDLE_ID,
    GRCh37ResourceBootstrapper,
    ResourceBootstrapError,
)
from ontseq_platform.resource_commands import (
    _authority_catalog_for_import,
    _bootstrap_if_official,
    handle_references_command,
)

CONFIGS = Path(__file__).resolve().parents[1] / "configs"
KNOWLEDGE = CONFIGS / "knowledge_bundles" / GRCH37_KNOWLEDGE_BUNDLE_ID
PANEL = CONFIGS / "panels" / GRCH37_PANEL_BUNDLE_ID


def _synthetic_grch37_annotation_cache(path: Path) -> None:
    mappings = []
    for line in (
        (PANEL / "source" / "target_gene_map.grch37.tsv").read_text(encoding="utf-8").splitlines()
    ):
        if not line or line.startswith("#") or line.startswith("source_target_label\t"):
            continue
        fields = line.split("\t")
        if fields[3] != ".":
            mappings.append(fields)
    selection_coordinates = {
        fields[3]: (int(fields[1]), int(fields[2]))
        for line in (PANEL / "derived" / "selection_panel.normalized.bed")
        .read_text(encoding="utf-8")
        .splitlines()
        if line and not line.startswith("#")
        for fields in (line.split("\t"),)
    }
    with closing(sqlite3.connect(path)) as connection:
        connection.executescript(
            """
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE genes (
              gene_id TEXT PRIMARY KEY, gene_name TEXT NOT NULL, chrom TEXT NOT NULL,
              start INTEGER NOT NULL, end INTEGER NOT NULL, strand TEXT NOT NULL,
              gene_type TEXT
            );
            CREATE TABLE transcripts (
              transcript_id TEXT PRIMARY KEY, gene_id TEXT NOT NULL, transcript_name TEXT,
              chrom TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL,
              strand TEXT NOT NULL, transcript_type TEXT, tags_json TEXT, mane_status TEXT,
              mane_refseq_id TEXT, appris TEXT, is_canonical INTEGER NOT NULL,
              is_basic INTEGER NOT NULL, cds_length INTEGER NOT NULL,
              transcript_length INTEGER NOT NULL
            );
            """
        )
        connection.execute("INSERT INTO metadata VALUES ('genome_build', 'GRCh37')")
        for index, fields in enumerate(mappings, start=1):
            _source_label, _target_label, chromosome, gene_id, native_name, _resolution = fields
            selection_start, selection_end = selection_coordinates[_target_label]
            start = selection_start + 1
            end = selection_end - 1
            connection.execute(
                "INSERT INTO genes VALUES (?, ?, ?, ?, ?, '+', 'protein_coding')",
                (gene_id, native_name, chromosome, start, end),
            )
            connection.execute(
                """
                INSERT INTO transcripts VALUES
                (?, ?, ?, ?, ?, ?, '+', 'protein_coding', '[]', NULL, NULL, NULL,
                 1, 1, 300, 500)
                """,
                (
                    f"ENST_SYNTHETIC_{index:03d}",
                    gene_id,
                    f"{native_name}-synthetic",
                    chromosome,
                    start,
                    end,
                ),
            )
        connection.commit()


class Grch37BootstrapTests(unittest.TestCase):
    def test_cli_status_and_full_validation_include_both_builds(self) -> None:
        for command in ("status", "validate"):
            with self.subTest(command=command):
                registries = {}
                for build in GenomeBuild:
                    registry = MagicMock()
                    registry.profiles = {f"AML_LCWGS_{build.value}": object()}
                    registry.diagnostics = ()
                    context = registry.resolve_profile.return_value
                    context.reference_bundle_id = f"{build.value}_TEST"
                    context.knowledge_bundle_id = f"KNOWLEDGE_{build.value}"
                    context.panel_bundle_id = None
                    registries[build] = registry
                installer = MagicMock()
                report = MagicMock()
                report.bundle_id = "SYNTHETIC"
                report.valid = True
                installer.status.return_value = (report,)
                installer.validate.return_value = report
                args = argparse.Namespace(
                    command="references",
                    references_command=command,
                    resource_root=Path("synthetic-resources"),
                    config_root=CONFIGS,
                    as_json=False,
                    bundle_id=None,
                )
                stream = io.StringIO()
                with (
                    patch(
                        "ontseq_platform.resource_commands.ReferenceBundleInstaller",
                        return_value=installer,
                    ),
                    patch(
                        "ontseq_platform.resource_commands.ResourceRegistry",
                        side_effect=lambda _root, *, active_build, selected=registries: selected[
                            active_build
                        ],
                    ),
                    patch("ontseq_platform.resource_commands._print_reports"),
                    redirect_stdout(stream),
                ):
                    self.assertTrue(handle_references_command(args))
                for build, registry in registries.items():
                    registry.resolve_profile.assert_called_once_with(
                        f"AML_LCWGS_{build.value}", verify_files=command == "validate"
                    )
                    self.assertIn(f"AML_LCWGS_{build.value}", stream.getvalue())

    def test_official_native_recipe_is_complete_and_has_no_mane_or_grch38_sources(self) -> None:
        catalog = ReferenceCatalog.discover(CONFIGS / "reference_bundles")
        recipe = catalog.get(GRCH37_REFERENCE_BUNDLE_ID)
        self.assertEqual(recipe.genome_build, GenomeBuild.GRCH37)
        sources = [resource for resource in recipe.resources if not resource.generated]
        self.assertEqual(len(sources), 9)
        self.assertLess(sum(resource.size_bytes or 0 for resource in sources), 2_000_000_000)
        self.assertEqual(
            recipe.resource("grch37_full_assembly_fasta_gz").sha256,
            "a37049f50cd181ffd057e8edde5f8d5a39afd3fed89ef17eb423323507ae2dd6",
        )
        for resource in sources:
            self.assertRegex(resource.sha256 or "", "^[0-9a-f]{64}$")
            self.assertGreater(resource.size_bytes or 0, 0)
            self.assertNotIn("hg38", resource.source_url or "")
            self.assertNotIn("GRCh38", resource.source_url or "")
            self.assertNotEqual(resource.role, "mane_gff3")
        self.assertEqual(
            recipe.resource("annotation_cache").derived_from,
            ["gencode_v19_all_gtf", "ucsc_hg19_cytobands"],
        )
        hg19_chrm = recipe.resource("ucsc_hg19_chrm_fasta_gz")
        self.assertEqual(hg19_chrm.role, "ucsc_hg19_chrm_fasta_archive")
        self.assertEqual(hg19_chrm.size_bytes, 5_537)
        self.assertEqual(
            hg19_chrm.sha256,
            "d2368da512728003230b13a394b22b3197ffe368c6219cf62e5b6dfc4a5d4891",
        )
        self.assertEqual(
            recipe.resource("ucsc_hg19_canonical25_fasta").derived_from,
            ["genome_fasta", "fasta_index", "ucsc_hg19_chrm_fasta_gz"],
        )
        self.assertEqual(
            recipe.resource("ucsc_hg19_canonical25_fasta").sha256,
            "1a814c7523f47e88d4fc1707a3bf9f18d51f2856a60477f56db13a833db3ca99",
        )
        self.assertEqual(
            recipe.resource("ucsc_hg19_canonical25_fasta_index").derived_from,
            ["ucsc_hg19_canonical25_fasta"],
        )
        self.assertEqual(
            recipe.resource("ucsc_hg19_canonical25_fasta_index").sha256,
            "1fb875bbf510ef0c7bd5b9b92bd430afadccdd0b2bc2106d6ffce1d9f0f00339",
        )
        self.assertTrue(recipe.resource("blacklist_bed").generated)

    def test_explicit_grch37_dictionaries_and_assays_are_published(self) -> None:
        bootstrap = GRCh37ResourceBootstrapper(packaged_config_root=CONFIGS)
        profiles = bootstrap._load_curated_profiles()
        self.assertEqual(
            tuple(
                (
                    profile.profile_id,
                    profile.assay_mode,
                    profile.reference_dictionary_contract,
                    profile.panel_bundle,
                    profile.adaptive_sampling,
                )
                for profile in profiles
            ),
            (
                (
                    "AML_LCWGS_GRCh37",
                    AssayMode.LOW_COVERAGE_WGS,
                    ReferenceDictionaryContract.EXACT_FULL,
                    None,
                    "disabled",
                ),
                (
                    "AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25",
                    AssayMode.LOW_COVERAGE_WGS,
                    ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25,
                    None,
                    "disabled",
                ),
                (
                    "AML_AS_111_GRCh37",
                    AssayMode.ADAPTIVE_SAMPLING,
                    ReferenceDictionaryContract.EXACT_FULL,
                    GRCH37_PANEL_BUNDLE_ID,
                    "enabled",
                ),
                (
                    "AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25",
                    AssayMode.ADAPTIVE_SAMPLING,
                    ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25,
                    GRCH37_PANEL_BUNDLE_ID,
                    "enabled",
                ),
            ),
        )
        for profile in profiles:
            self.assertEqual(profile.reference_bundle, GRCH37_REFERENCE_BUNDLE_ID)
            self.assertEqual(profile.genome_build, GenomeBuild.GRCH37)
        self.assertEqual(
            GRCH37_PROFILE_IDS,
            (
                "AML_LCWGS_GRCh37",
                "AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25",
                "AML_AS_111_GRCh37",
                "AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25",
            ),
        )

    def test_knowledge_is_checksum_locked_coordinate_free_and_source_attributed(self) -> None:
        bundle = KnowledgeBundle.model_validate(
            yaml.safe_load((KNOWLEDGE / "bundle.yaml").read_text(encoding="utf-8"))
        )
        self.assertEqual(bundle.genome_build, GenomeBuild.GRCH37)
        self.assertFalse(bundle.coordinate_bearing)
        for resource in bundle.resources:
            path = KNOWLEDGE / resource.path
            self.assertEqual(resource.sha256, sha256_file(path))
            self.assertEqual(resource.size_bytes, path.stat().st_size)
            self.assertIsNone(resource.coordinate_system)
        resource = KNOWLEDGE / "hematology_rearrangements.v0.3.grch37.json"
        lock = AmlKnowledgeLock.model_validate(
            json.loads(
                (KNOWLEDGE / "hematology_rearrangements.v0.3.grch37.lock.json").read_text(
                    encoding="utf-8"
                )
            )
        )
        records = load_aml_knowledge(resource, lock)
        original_path = (
            CONFIGS
            / "knowledge_bundles"
            / KNOWLEDGE_BUNDLE_ID
            / "hematology_rearrangements.v0.3.json"
        )
        original = json.loads(original_path.read_text(encoding="utf-8"))
        derived = json.loads(resource.read_text(encoding="utf-8"))
        self.assertEqual(derived["records"], original["records"])
        self.assertEqual(len(records), 38)
        self.assertEqual(sum(len(record.pathologies) for record in records), 74)
        self.assertNotIn("panel", derived["scope"])
        self.assertEqual(
            derived["scope"]["provenance"]["source_sha256"], sha256_file(original_path)
        )

    def test_activation_requires_a_valid_native_reference(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            bootstrap = GRCh37ResourceBootstrapper(raw, packaged_config_root=CONFIGS)
            with self.assertRaisesRegex(ResourceBootstrapError, "must be installed and valid"):
                bootstrap.activate()
            self.assertFalse((Path(raw) / "profiles").exists())

    def test_grch37_activation_is_idempotent_and_never_touches_grch38(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            annotation_cache = root / "cache"
            _synthetic_grch37_annotation_cache(annotation_cache)
            sentinel = root / "knowledge" / KNOWLEDGE_BUNDLE_ID / "existing.txt"
            sentinel.parent.mkdir(parents=True)
            sentinel.write_bytes(b"existing GRCh38 resources remain unchanged\n")
            bootstrap = GRCh37ResourceBootstrapper(root, packaged_config_root=CONFIGS)
            # This transaction test isolates already-validated cache activation; native cache
            # compilation/build rejection is exercised by reference-catalog integration tests.
            with patch.object(
                bootstrap, "_reference_annotation_cache", return_value=annotation_cache
            ):
                first = bootstrap.activate()
                second = bootstrap.activate()
            self.assertIsNotNone(first.panel_summary)
            assert first.panel_summary is not None
            self.assertEqual(first.panel_bundle_id, GRCH37_PANEL_BUNDLE_ID)
            self.assertEqual(first.panel_summary.compilation.target_count, 111)
            self.assertEqual(first.panel_summary.compilation.resolved_target_count, 108)
            self.assertEqual(first.panel_summary.compilation.roi_interval_count, 108)
            self.assertEqual(
                first.panel_summary.compilation.unresolved_targets,
                (
                    "CT45A2_REVIEW_REQUIRED",
                    "IGH_REVIEW_REQUIRED",
                    "GPR128_REVIEW_REQUIRED",
                ),
            )
            self.assertEqual(len(first.activated_paths), 6)
            self.assertEqual(second.activated_paths, ())
            self.assertEqual(len(second.already_active_paths), 6)
            self.assertTrue((root / "panels" / GRCH37_PANEL_BUNDLE_ID).is_dir())
            self.assertEqual(sentinel.read_bytes(), b"existing GRCh38 resources remain unchanged\n")

    def test_grch37_repair_restores_pinned_knowledge(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            annotation_cache = root / "cache"
            _synthetic_grch37_annotation_cache(annotation_cache)
            bootstrap = GRCh37ResourceBootstrapper(root, packaged_config_root=CONFIGS)
            with patch.object(
                bootstrap, "_reference_annotation_cache", return_value=annotation_cache
            ):
                bootstrap.activate()
                target = (
                    root
                    / "knowledge"
                    / GRCH37_KNOWLEDGE_BUNDLE_ID
                    / "hematology_rearrangements.v0.3.grch37.json"
                )
                original = target.read_bytes()
                target.write_bytes(b"corrupt synthetic fixture\n")
                with self.assertRaisesRegex(FileExistsError, "overwrite a different"):
                    bootstrap.activate()
                repaired = bootstrap.repair()
            self.assertEqual(target.read_bytes(), original)
            self.assertEqual(len(repaired.repaired_paths), 1)

    def test_cross_build_curated_profile_is_rejected_before_activation(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            configs = root / "configs"
            (configs / "profiles").mkdir(parents=True)
            profile = yaml.safe_load(
                (CONFIGS / "profiles" / "AML_LCWGS_GRCh37.yaml").read_text(encoding="utf-8")
            )
            profile["genome_build"] = "GRCh38"
            (configs / "profiles" / "AML_LCWGS_GRCh37.yaml").write_text(
                yaml.safe_dump(profile), encoding="utf-8"
            )
            bootstrap = GRCh37ResourceBootstrapper(root / "resources", packaged_config_root=configs)
            with self.assertRaisesRegex(ResourceBootstrapError, "is not GRCh37"):
                bootstrap._load_curated_profiles()

    def test_cross_build_curated_knowledge_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            source = root / "knowledge"
            shutil.copytree(KNOWLEDGE, source)
            manifest = source / "bundle.yaml"
            payload = yaml.safe_load(manifest.read_text(encoding="utf-8"))
            payload["genome_build"] = "GRCh38"
            manifest.write_text(yaml.safe_dump(payload), encoding="utf-8")
            bootstrap = GRCh37ResourceBootstrapper(root, packaged_config_root=CONFIGS)
            with self.assertRaisesRegex(ResourceBootstrapError, "is not GRCh37"):
                bootstrap._require_curated_bundle(
                    source, KnowledgeBundle, GRCH37_KNOWLEDGE_BUNDLE_ID
                )

    def test_cli_activates_grch37_family_and_reports_native_panel(self) -> None:
        result = MagicMock()
        result.profile_ids = GRCH37_PROFILE_IDS
        result.panel_summary.compilation.roi_interval_count = 108
        result.repaired_paths = ()
        bootstrap = MagicMock()
        bootstrap.activate.return_value = result
        stream = io.StringIO()
        args = argparse.Namespace(resource_root=Path("resources"), config_root=CONFIGS)
        with (
            patch(
                "ontseq_platform.resource_commands.GRCh37ResourceBootstrapper",
                return_value=bootstrap,
            ),
            patch("ontseq_platform.resource_commands.GRCh38ResourceBootstrapper") as grch38,
            redirect_stdout(stream),
        ):
            _bootstrap_if_official(args, GRCH37_REFERENCE_BUNDLE_ID)
        bootstrap.activate.assert_called_once_with()
        grch38.assert_not_called()
        self.assertIn("AML_LCWGS_GRCh37", stream.getvalue())
        self.assertIn("AML_AS_111_GRCh37", stream.getvalue())
        self.assertIn("panel intervals: 108 ROI", stream.getvalue())

    def test_official_grch37_import_requires_packaged_authority(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            source = Path(raw)
            (source / "bundle.recipe.yaml").write_text(
                f"bundle_id: {GRCH37_REFERENCE_BUNDLE_ID}\n", encoding="utf-8"
            )
            catalog = MagicMock()
            with patch(
                "ontseq_platform.resource_commands._packaged_authority_catalog",
                return_value=catalog,
            ):
                self.assertIs(_authority_catalog_for_import(source), catalog)
            catalog.get.assert_called_once_with(GRCH37_REFERENCE_BUNDLE_ID)


if __name__ == "__main__":
    unittest.main()
