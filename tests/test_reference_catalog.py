from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import sqlite3
import tempfile
import unittest
from contextlib import closing, contextmanager, redirect_stdout
from pathlib import Path
from unittest.mock import patch

import yaml

from ontseq_platform.models import (
    CoordinateSystem,
    GenomeBuild,
    ReferenceBundle,
    ResourceFile,
)
from ontseq_platform.reference_catalog import (
    ReferenceBundleInstaller,
    ReferenceCatalog,
    ResourceValidationState,
    _fasta_index,
    _gunzip,
    _normalize_ucsc_context_table,
    validate_fasta_fai_consistency,
    validate_reference_bundle_directory,
)

FIXTURES = Path(__file__).parent / "fixtures" / "reference_cache"
FIXTURE_CATALOG = Path(__file__).parent / "fixtures" / "reference_catalog"


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


TINY_FASTA = b">chr1\nACGT\n"
TINY_FAI = b"chr1\t4\t6\t4\t5\n"


@contextmanager
def _allow_tiny_transaction_fixture():
    """Keep transaction tests byte-small without adding a production validation escape hatch."""

    with (
        patch("ontseq_platform.reference.validate_canonical_reference"),
        patch("ontseq_platform.reference_catalog.validate_canonical_reference"),
    ):
        yield


def _source_resource(
    resource_id: str,
    role: str,
    path: str,
    content: bytes,
    *,
    coordinate_system: CoordinateSystem | None = None,
) -> ResourceFile:
    return ResourceFile(
        resource_id=resource_id,
        role=role,
        path=path,
        sha256=_sha256(content),
        size_bytes=len(content),
        source_url=f"https://fixtures.invalid/{resource_id}",
        release="fixture-v1",
        coordinate_system=coordinate_system,
    )


def _recipe() -> tuple[ReferenceBundle, dict[str, bytes]]:
    contents = {
        "genome": TINY_FASTA,
        "fai": TINY_FAI,
        "gencode": (FIXTURES / "gencode.v50.fragment.gtf").read_bytes(),
        "mane": (FIXTURES / "MANE.GRCh38.v1.5.fragment.gff3").read_bytes(),
        "cytobands": (FIXTURES / "cytoBand.hg38.fragment.tsv").read_bytes(),
    }
    resources = [
        _source_resource("genome", "genome_fasta", "genome.fa", contents["genome"]),
        _source_resource("fai", "fasta_index", "genome.fa.fai", contents["fai"]),
        _source_resource(
            "gencode",
            "gencode_gtf",
            "sources/gencode.gtf",
            contents["gencode"],
            coordinate_system=CoordinateSystem.ONE_BASED_INCLUSIVE,
        ),
        _source_resource(
            "mane",
            "mane_gff3",
            "sources/mane.gff3",
            contents["mane"],
            coordinate_system=CoordinateSystem.ONE_BASED_INCLUSIVE,
        ),
        _source_resource(
            "cytobands",
            "cytobands",
            "sources/cytobands.tsv",
            contents["cytobands"],
            coordinate_system=CoordinateSystem.ZERO_BASED_HALF_OPEN,
        ),
        ResourceFile(
            resource_id="reference_lock",
            role="reference_lock",
            path="reference.lock.json",
            generated=True,
            derived_from=["fai"],
        ),
        ResourceFile(
            resource_id="chromosome_sizes",
            role="chromosome_sizes",
            path="chromosome-sizes.tsv",
            generated=True,
            derived_from=["fai"],
            coordinate_system=CoordinateSystem.ZERO_BASED_HALF_OPEN,
        ),
        ResourceFile(
            resource_id="annotation_cache",
            role="annotation_cache",
            path="annotation.sqlite",
            generated=True,
            derived_from=["gencode", "mane", "cytobands"],
            coordinate_system=CoordinateSystem.ZERO_BASED_HALF_OPEN,
        ),
    ]
    return (
        ReferenceBundle(
            bundle_id="GRCh38_FIXTURE_v1",
            version="fixture-v1",
            genome_build=GenomeBuild.GRCH38,
            resources=resources,
            reference_lock_resource_id="reference_lock",
            fasta_resource_id="genome",
            fai_resource_id="fai",
            annotation_cache_resource_id="annotation_cache",
        ),
        contents,
    )


def _change_resource(
    bundle: ReferenceBundle, resource_id: str, **updates: object
) -> ReferenceBundle:
    return bundle.model_copy(
        update={
            "resources": [
                resource.model_copy(update=updates)
                if resource.resource_id == resource_id
                else resource
                for resource in bundle.resources
            ]
        }
    )


class _MemoryOpener:
    def __init__(self, contents: dict[str, bytes]) -> None:
        self.contents = contents
        self.requests: list[str] = []

    def __call__(self, url: str) -> io.BytesIO:
        self.requests.append(url)
        resource_id = url.rsplit("/", 1)[-1]
        return io.BytesIO(self.contents[resource_id])


class DeterministicReferenceDerivationTests(unittest.TestCase):
    def test_gzip_fasta_is_decompressed_and_indexed_without_external_tools(self) -> None:
        fasta = b">chr1 fixture\nACGT\nAC\n>chrM\nAAA\n"
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            archive = root / "genome.fa.gz"
            with (
                archive.open("wb") as raw_handle,
                gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as compressed,
            ):
                compressed.write(fasta)
            uncompressed = root / "genome.fa"
            index = root / "genome.fa.fai"

            _gunzip(archive, uncompressed)
            _fasta_index(uncompressed, index)

            self.assertEqual(uncompressed.read_bytes(), fasta)
            self.assertEqual(
                index.read_text(encoding="utf-8").splitlines(),
                ["chr1\t6\t14\t4\t5", "chrM\t3\t28\t3\t4"],
            )

    def test_fasta_fai_consistency_checks_lengths_offsets_and_wrapping(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            fasta = root / "genome.fa"
            fai = root / "genome.fa.fai"
            fasta.write_bytes(b">chr1 fixture\nACGT\nAC\n>chrM\nAAA\n")
            _fasta_index(fasta, fai)

            validate_fasta_fai_consistency(fasta, fai)

            fai.write_text("chr1\t6\t0\t4\t5\nchrM\t3\t28\t3\t4\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not derived from the installed FASTA"):
                validate_fasta_fai_consistency(fasta, fai)


class ReferenceBundleInstallerTests(unittest.TestCase):
    def setUp(self) -> None:
        self._fixture_gate = _allow_tiny_transaction_fixture()
        self._fixture_gate.__enter__()
        self.addCleanup(self._fixture_gate.__exit__, None, None, None)

    def test_install_stages_validates_compiles_and_activates_bundle(self) -> None:
        recipe, contents = _recipe()
        opener = _MemoryOpener(contents)
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            installed = ReferenceBundleInstaller(root, opener=opener).install(recipe)

            self.assertTrue(installed.validation.valid)
            self.assertEqual(installed.path, root / "references" / recipe.bundle_id)
            self.assertTrue((installed.path / "bundle.yaml").is_file())
            self.assertTrue((installed.path / "annotation.sqlite").is_file())
            self.assertTrue((installed.path / "reference.lock.json").is_file())
            active = ReferenceBundle.model_validate(
                yaml.safe_load((installed.path / "bundle.yaml").read_text(encoding="utf-8"))
            )
            self.assertTrue(all(resource.sha256 for resource in active.resources))
            self.assertEqual(len(opener.requests), 5)
            self.assertFalse(any((root / "references" / ".staging").iterdir()))

    def test_failed_checksum_never_activates_partial_bundle(self) -> None:
        recipe, contents = _recipe()
        contents["mane"] = b"not the pinned MANE source"
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            installer = ReferenceBundleInstaller(root, opener=_MemoryOpener(contents))
            with self.assertRaisesRegex(ValueError, "expected (?:SHA256|378 bytes)"):
                installer.install(recipe)

            self.assertFalse((root / "references" / recipe.bundle_id).exists())
            self.assertFalse(any((root / "references" / ".staging").iterdir()))

    def test_install_rejects_a_release_pinned_fai_for_different_fasta_bytes(self) -> None:
        recipe, contents = _recipe()
        mismatched_fai = b"chr1\t5\t6\t5\t6\n"
        contents["fai"] = mismatched_fai
        recipe = recipe.model_copy(
            update={
                "resources": [
                    item.model_copy(
                        update={
                            "sha256": _sha256(mismatched_fai),
                            "size_bytes": len(mismatched_fai),
                        }
                    )
                    if item.resource_id == "fai"
                    else item
                    for item in recipe.resources
                ]
            }
        )
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            installer = ReferenceBundleInstaller(root, opener=_MemoryOpener(contents))

            with self.assertRaisesRegex(ValueError, "not derived from the installed FASTA"):
                installer.install(recipe)

            self.assertFalse((root / "references" / recipe.bundle_id).exists())

    def test_repair_replaces_only_bad_source_and_rebuilds_its_derivative(self) -> None:
        recipe, contents = _recipe()
        opener = _MemoryOpener(contents)
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            installer = ReferenceBundleInstaller(root, opener=opener)
            installed = installer.install(recipe)
            lock_before = (installed.path / "reference.lock.json").read_bytes()
            cache_before = (installed.path / "annotation.sqlite").read_bytes()
            (installed.path / "sources" / "gencode.gtf").write_text("corrupted", encoding="utf-8")
            opener.requests.clear()

            repaired = installer.repair(recipe)

            self.assertTrue(repaired.validation.valid)
            self.assertEqual(opener.requests, ["https://fixtures.invalid/gencode"])
            self.assertEqual((installed.path / "reference.lock.json").read_bytes(), lock_before)
            self.assertEqual((installed.path / "annotation.sqlite").read_bytes(), cache_before)

    def test_offline_repair_never_attempts_a_network_fallback(self) -> None:
        recipe, contents = _recipe()
        opener = _MemoryOpener(contents)
        with tempfile.TemporaryDirectory() as raw:
            installer = ReferenceBundleInstaller(Path(raw), opener=opener)
            installed = installer.install(recipe)
            opener.requests.clear()
            (installed.path / "sources" / "mane.gff3").unlink()

            with self.assertRaisesRegex(RuntimeError, "offline mode"):
                installer.repair(recipe, offline=True)
            self.assertEqual(opener.requests, [])

    def test_repair_cannot_silently_upgrade_an_installed_reference_release(self) -> None:
        recipe, contents = _recipe()
        with tempfile.TemporaryDirectory() as raw:
            installer = ReferenceBundleInstaller(Path(raw), opener=_MemoryOpener(contents))
            installer.install(recipe)
            changed_release = recipe.model_copy(update={"version": "fixture-v2"})

            with self.assertRaisesRegex(ValueError, "exact installed bundle"):
                installer.repair(changed_release)

    def test_repair_rejects_changed_source_and_generator_contracts_at_same_version(self) -> None:
        recipe, contents = _recipe()
        changed_contracts = (
            _change_resource(
                recipe,
                "gencode",
                source_url="https://untrusted-mirror.invalid/gencode",
            ),
            _change_resource(
                recipe,
                "annotation_cache",
                derived_from=["gencode", "mane"],
            ),
        )
        for changed_recipe in changed_contracts:
            with (
                self.subTest(changed=changed_recipe.resources),
                tempfile.TemporaryDirectory() as raw,
            ):
                installer = ReferenceBundleInstaller(Path(raw), opener=_MemoryOpener(contents))
                installer.install(recipe)

                with self.assertRaisesRegex(
                    ValueError,
                    "immutable catalog source/generator contract.*new bundle ID/version",
                ):
                    installer.repair(changed_recipe)

    def test_import_is_local_checksum_validated_and_atomic(self) -> None:
        recipe, contents = _recipe()
        with (
            tempfile.TemporaryDirectory() as source_raw,
            tempfile.TemporaryDirectory() as target_raw,
        ):
            source = ReferenceBundleInstaller(
                Path(source_raw), opener=_MemoryOpener(contents)
            ).install(recipe)
            target_installer = ReferenceBundleInstaller(Path(target_raw))

            imported = target_installer.import_bundle(source.path)

            self.assertTrue(imported.validation.valid)
            self.assertEqual(imported.bundle.bundle_id, recipe.bundle_id)
            self.assertTrue((imported.path / "annotation.sqlite").is_file())

    def test_validate_and_import_reject_a_re_pinned_cross_file_mismatch(self) -> None:
        recipe, contents = _recipe()
        with (
            tempfile.TemporaryDirectory() as source_raw,
            tempfile.TemporaryDirectory() as target_raw,
        ):
            source = ReferenceBundleInstaller(
                Path(source_raw), opener=_MemoryOpener(contents)
            ).install(recipe)
            mismatched_fai = b"chr1\t5\t6\t5\t6\n"
            fai_path = source.path / "genome.fa.fai"
            fai_path.write_bytes(mismatched_fai)
            manifest_path = source.path / "bundle.yaml"
            manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
            for resource in manifest["resources"]:
                if resource["resource_id"] == "fai":
                    resource["sha256"] = _sha256(mismatched_fai)
                    resource["size_bytes"] = len(mismatched_fai)
            manifest_path.write_text(
                yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8", newline="\n"
            )

            report = validate_reference_bundle_directory(source.path)
            self.assertFalse(report.valid)
            self.assertIn("not derived from the installed FASTA", " ".join(report.errors))
            with self.assertRaisesRegex(ValueError, "not a valid activated reference bundle"):
                ReferenceBundleInstaller(Path(target_raw)).import_bundle(source.path)

    def test_pinned_recipe_source_tree_import_compiles_without_network(self) -> None:
        source = FIXTURE_CATALOG / "GRCh38_FIXTURE_v1"
        with tempfile.TemporaryDirectory() as raw:
            opener = _MemoryOpener({})
            installer = ReferenceBundleInstaller(Path(raw), opener=opener)

            imported = installer.import_bundle(source)

            self.assertTrue(imported.validation.valid)
            self.assertEqual(imported.bundle.bundle_id, "GRCh38_FIXTURE_v1")
            self.assertTrue((imported.path / "annotation.sqlite").is_file())
            self.assertEqual(opener.requests, [])

    def test_catalog_authority_rejects_self_repinned_import_for_known_bundle_id(self) -> None:
        authority, contents = _recipe()
        self_repinned = _change_resource(
            authority,
            "gencode",
            source_url="https://untrusted-mirror.invalid/gencode",
        )
        catalog = ReferenceCatalog({authority.bundle_id: authority})
        with (
            tempfile.TemporaryDirectory() as source_raw,
            tempfile.TemporaryDirectory() as target_raw,
        ):
            source = ReferenceBundleInstaller(
                Path(source_raw), opener=_MemoryOpener(contents)
            ).install(self_repinned)

            with self.assertRaisesRegex(ValueError, "immutable catalog source/generator contract"):
                ReferenceBundleInstaller(Path(target_raw)).import_bundle(
                    source.path,
                    authority_catalog=catalog,
                )
            self.assertFalse((Path(target_raw) / "references" / authority.bundle_id).exists())

    def test_catalog_authority_rejects_repinned_generated_cache_content(self) -> None:
        authority, contents = _recipe()
        catalog = ReferenceCatalog({authority.bundle_id: authority})
        with (
            tempfile.TemporaryDirectory() as source_raw,
            tempfile.TemporaryDirectory() as target_raw,
        ):
            source = ReferenceBundleInstaller(
                Path(source_raw), opener=_MemoryOpener(contents)
            ).install(authority)
            cache_path = source.path / "annotation.sqlite"
            with closing(sqlite3.connect(cache_path)) as connection:
                connection.execute("UPDATE genes SET gene_name = 'TAMPERED'")
                connection.commit()
            active = ReferenceBundle.model_validate(
                yaml.safe_load((source.path / "bundle.yaml").read_text(encoding="utf-8"))
            )
            repinned = _change_resource(
                active,
                "annotation_cache",
                sha256=_sha256(cache_path.read_bytes()),
                size_bytes=cache_path.stat().st_size,
            )
            (source.path / "bundle.yaml").write_text(
                yaml.safe_dump(
                    repinned.model_dump(mode="json", exclude_none=True),
                    sort_keys=False,
                ),
                encoding="utf-8",
                newline="\n",
            )
            self.assertTrue(validate_reference_bundle_directory(source.path).valid)

            with self.assertRaisesRegex(ValueError, "deterministic catalog derivation"):
                ReferenceBundleInstaller(Path(target_raw)).import_bundle(
                    source.path,
                    authority_catalog=catalog,
                )
            self.assertFalse((Path(target_raw) / "references" / authority.bundle_id).exists())

    def test_catalog_authority_does_not_claim_custom_bundle_ids(self) -> None:
        authority, contents = _recipe()
        custom = authority.model_copy(update={"bundle_id": "CUSTOM_GRCh38_FIXTURE_v1"})
        catalog = ReferenceCatalog({authority.bundle_id: authority})
        with (
            tempfile.TemporaryDirectory() as source_raw,
            tempfile.TemporaryDirectory() as target_raw,
        ):
            source = ReferenceBundleInstaller(
                Path(source_raw), opener=_MemoryOpener(contents)
            ).install(custom)

            imported = ReferenceBundleInstaller(Path(target_raw)).import_bundle(
                source.path,
                authority_catalog=catalog,
            )

            self.assertEqual(imported.bundle.bundle_id, custom.bundle_id)
            self.assertTrue(imported.validation.valid)

    def test_cli_custom_import_does_not_require_a_reference_catalog(self) -> None:
        try:
            from ontseq_platform.resource_commands import handle_references_command
        except ModuleNotFoundError as error:
            if error.name != "ontseq_platform.resource_commands":
                raise
            self.skipTest("the references CLI is introduced by the profile/pipeline stack stage")
        recipe, contents = _recipe()
        custom = recipe.model_copy(update={"bundle_id": "CUSTOM_NO_CATALOG_v1"})
        with (
            tempfile.TemporaryDirectory() as source_raw,
            tempfile.TemporaryDirectory() as target_raw,
            tempfile.TemporaryDirectory() as config_raw,
            redirect_stdout(io.StringIO()),
        ):
            source = ReferenceBundleInstaller(
                Path(source_raw), opener=_MemoryOpener(contents)
            ).install(custom)
            args = argparse.Namespace(
                command="references",
                references_command="import",
                resource_root=Path(target_raw),
                config_root=Path(config_raw),
                path=source.path,
            )

            self.assertTrue(handle_references_command(args))
            imported = Path(target_raw) / "references" / custom.bundle_id
            self.assertTrue(validate_reference_bundle_directory(imported).valid)

    def test_reference_namespace_symlinks_are_rejected_before_repair(self) -> None:
        recipe, contents = _recipe()
        with (
            tempfile.TemporaryDirectory() as external_raw,
            tempfile.TemporaryDirectory() as target_raw,
        ):
            external_root = Path(external_raw)
            external = ReferenceBundleInstaller(
                external_root, opener=_MemoryOpener(contents)
            ).install(recipe)
            target_root = Path(target_raw)
            try:
                (target_root / "references").symlink_to(
                    external_root / "references", target_is_directory=True
                )
            except OSError as exc:
                self.skipTest(f"directory symlinks are unavailable: {exc}")
            before = (external.path / "bundle.yaml").read_bytes()

            with self.assertRaisesRegex(ValueError, "symbolic link or junction"):
                ReferenceBundleInstaller(target_root).repair(recipe)

            self.assertEqual((external.path / "bundle.yaml").read_bytes(), before)


class ReferenceCatalogTests(unittest.TestCase):
    def test_ucsc_context_tables_compile_to_standard_half_open_bed(self) -> None:
        fixtures = {
            "repeatmasker": (
                "0\t100\t10\t0\t0\tchr1\t10\t20\t0\t+\tAluY\tSINE\tAlu\n",
                "chr1\t10\t20\tAluY|SINE|Alu\n",
            ),
            "simple_repeats": (
                "0\tchr2\t30\t40\ttrf\t2\t5.0\n",
                "chr2\t30\t40\ttrf\n",
            ),
            "segmental_duplication": (
                "0\tchr3\t50\t75\tdupA\textra\n",
                "chr3\t50\t75\tdupA\n",
            ),
        }
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            for role, (source_line, expected) in fixtures.items():
                with self.subTest(role=role):
                    source = root / f"{role}.txt.gz"
                    output = root / f"{role}.bed"
                    with gzip.open(source, "wt", encoding="utf-8") as handle:
                        handle.write(source_line)
                    _normalize_ucsc_context_table(source, output, role=role)
                    self.assertEqual(output.read_text(encoding="utf-8"), expected)

    def test_committed_fixture_catalog_pins_every_source_byte(self) -> None:
        catalog = ReferenceCatalog.discover(FIXTURE_CATALOG)
        bundle = catalog.get("GRCh38_FIXTURE_v1")
        bundle_root = FIXTURE_CATALOG / bundle.bundle_id

        for resource in bundle.resources:
            if resource.generated:
                self.assertIsNone(resource.sha256)
                continue
            source = bundle_root.joinpath(*resource.path.split("/"))
            self.assertEqual(resource.size_bytes, source.stat().st_size)
            self.assertEqual(resource.sha256, _sha256(source.read_bytes()))

    def test_tiny_ci_fixture_has_no_productive_grch38_activation_path(self) -> None:
        source = FIXTURE_CATALOG / "GRCh38_FIXTURE_v1"
        with (
            tempfile.TemporaryDirectory() as raw,
            self.assertRaisesRegex(ValueError, "canonical assembly validation failed"),
        ):
            ReferenceBundleInstaller(Path(raw)).import_bundle(source)

    def test_official_release_lock_is_discoverable_and_every_source_is_pinned(self) -> None:
        inventory_root = (
            Path(__file__).parents[1]
            / "configs"
            / "reference_bundles"
            / "GRCh38_GENCODE50_MANE1.5_v1"
        )
        catalog = ReferenceCatalog.discover(inventory_root)

        bundle = catalog.get("GRCh38_GENCODE50_MANE1.5_v1")
        self.assertEqual(bundle.genome_build.value, "GRCh38")
        sources = [resource for resource in bundle.resources if not resource.generated]
        self.assertEqual(len(sources), 9)
        self.assertTrue(all(resource.sha256 is not None for resource in sources))
        self.assertTrue(all(resource.size_bytes is not None for resource in sources))
        self.assertTrue(all(resource.source_url for resource in sources))
        self.assertTrue(all(resource.source_date for resource in sources))
        self.assertEqual(
            bundle.resource("hoffman_umap_hg38_k100").role,
            "mappability",
        )
        mappability = bundle.resource("hoffman_umap_hg38_k100")
        self.assertEqual(mappability.resource_id, "hoffman_umap_hg38_k100")
        self.assertEqual(
            mappability.source_url,
            "https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz",
        )
        self.assertEqual(mappability.coordinate_system.value, "zero_based_half_open")

    def test_catalog_discovers_only_named_recipe_or_bundle_manifests(self) -> None:
        recipe, _ = _recipe()
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            bundle_dir = root / recipe.bundle_id
            bundle_dir.mkdir()
            (bundle_dir / "bundle.recipe.yaml").write_text(
                yaml.safe_dump(recipe.model_dump(mode="json"), sort_keys=False),
                encoding="utf-8",
            )
            (root / "loose.yaml").write_text("not: a bundle\n", encoding="utf-8")

            catalog = ReferenceCatalog.discover(root)

        self.assertEqual([bundle.bundle_id for bundle in catalog.list()], [recipe.bundle_id])

    def test_active_grch38_installer_rejects_grch37_catalog_recipes(self) -> None:
        recipe, _ = _recipe()
        grch37 = recipe.model_copy(update={"genome_build": GenomeBuild.GRCH37})
        with tempfile.TemporaryDirectory() as raw:
            manifest = Path(raw) / "bundle.yaml"
            manifest.write_text(
                yaml.safe_dump(grch37.model_dump(mode="json"), sort_keys=False),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "only GRCh38"):
                ReferenceCatalog.from_manifests([manifest])

    def test_status_reports_checksum_mismatch_without_changing_files(self) -> None:
        recipe, contents = _recipe()
        with tempfile.TemporaryDirectory() as raw, _allow_tiny_transaction_fixture():
            installer = ReferenceBundleInstaller(Path(raw), opener=_MemoryOpener(contents))
            installed = installer.install(recipe)
            gencode = installed.path / "sources" / "gencode.gtf"
            gencode.write_text("corrupt", encoding="utf-8")

            status = installer.status()

        self.assertEqual(len(status), 1)
        states = {item.resource_id: item.state for item in status[0].resources}
        self.assertEqual(states["gencode"], ResourceValidationState.SIZE_MISMATCH)


if __name__ == "__main__":
    unittest.main()
