from __future__ import annotations

import os
import unittest
from pathlib import Path

from ontseq_platform.reference_catalog import ReferenceBundleInstaller, ReferenceCatalog


@unittest.skipUnless(
    os.environ.get("ONTSEQ_FULL_GRCH38_BUNDLE_SMOKE") == "1",
    "set ONTSEQ_FULL_GRCH38_BUNDLE_SMOKE=1 for the opt-in real-resource smoke",
)
class FullGrch38ReferenceBundleSmokeTests(unittest.TestCase):
    """Explicitly opt-in because this test downloads and installs multi-gigabyte resources."""

    def test_install_or_validate_release_locked_bundle(self) -> None:
        root_value = os.environ.get("ONTSEQ_FULL_GRCH38_RESOURCE_ROOT")
        if not root_value:
            self.fail("the full smoke requires ONTSEQ_FULL_GRCH38_RESOURCE_ROOT")
        manifest = Path(
            os.environ.get("ONTSEQ_FULL_GRCH38_BUNDLE_MANIFEST")
            or Path(__file__).parents[1]
            / "configs"
            / "reference_bundles"
            / "GRCh38_GENCODE50_MANE1.5_v1"
            / "bundle.recipe.yaml"
        )
        bundle = ReferenceCatalog.from_manifests([manifest]).get("GRCh38_GENCODE50_MANE1.5_v1")
        installer = ReferenceBundleInstaller(Path(root_value))
        destination = installer.references_root / bundle.bundle_id
        installed = installer.repair(bundle) if destination.is_dir() else installer.install(bundle)
        self.assertTrue(installed.validation.valid)


if __name__ == "__main__":
    unittest.main()
