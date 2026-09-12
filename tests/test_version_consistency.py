"""Synthetic release metadata rejects stale Desktop resource defaults."""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from scripts import check_version_consistency as guard


class VersionConsistencyTests(unittest.TestCase):
    def setUp(self) -> None:
        root = Path(self.enterContext(TemporaryDirectory()))
        paths = {
            name: root / value.relative_to(guard.ROOT)
            for name, value in vars(guard).items()
            if name.isupper() and isinstance(value, Path)
        }
        self.enterContext(patch.multiple(guard, **paths))
        self.files = {
            "PYPROJECT": '[project]\nversion = "9.8.7"\n',
            "UV_LOCK": (
                '[[package]]\nname = "ontseq-platform"\nversion = "9.8.7"\n'
                'source = { editable = "." }\n'
            ),
            "PACKAGE_INIT": '__version__ = "9.8.7"\n',
            "CITATION": "version: 9.8.7\n",
            "DESKTOP_PROJECT": "<Version>9.8.7</Version>",
            "DESKTOP_VERSION": (
                'public const string CoreValue = "9.8.7";\n'
                'public const string Value = CoreValue + "-engineering";\n'
            ),
            "DESKTOP_LAUNCHER": "private const string ReleaseVersion = DesktopVersion.CoreValue;",
            "DESKTOP_MODELS": (
                'public const string DefaultResourceRootWsl = "~/.local/share/ontseq/resources-v"'
                " + DesktopVersion.CoreValue;"
            ),
            "DESKTOP_SETTINGS_EXAMPLE": (
                '{"resourceRootWsl": "~/.local/share/ontseq/resources-v9.8.7"}'
            ),
            "DESKTOP_CHANGELOG": "## 9.8.7-engineering\n## 0.6.2-engineering\n",
            "DESKTOP_README": (
                "Current engineering build: Desktop/Core v9.8.7.\n"
                "## v9.8.7 user path\nruntime-v9.8.7\n"
                "~/.local/share/ontseq/resources-v9.8.7\n"
                "Historical default: ~/.local/share/ontseq/resources\n"
            ),
            "DESKTOP_FIRST_RUN": (
                "# ONTSeq Desktop v9.8.7\n"
                "ontseq-desktop-v9.8.7-win-x64-setup-engineering\n"
                "~/.local/share/ontseq/resources-v9.8.7\n"
            ),
            "DESKTOP_ISOLATED_TESTING": "~/.local/share/ontseq/resources-v9.8.7\n",
            "DESKTOP_PACKAGED_README": "~/.local/share/ontseq/resources-v9.8.7\n",
            "OPERATOR_FIRST_RUN": "~/.local/share/ontseq/resources-v9.8.7\n",
            "DESKTOP_API_CONTRACT": "~/.local/share/ontseq/resources-v9.8.7\n",
            "REFERENCE_SYSTEM": (
                "~/.local/share/ontseq/resources-v9.8.7\n"
                "GRCh37_GENCODE19_HG19_v2\nAML_AS_111_GRCh37_v1\nPipelineResult 0.3.0\n"
            ),
            "DESKTOP_WORKFLOW": (
                'ONTSEQ_VERSION: "9.8.7"\n'
                "name: ontseq-linux-runtime-v${{ env.ONTSEQ_VERSION }}\n"
                "name: ontseq-desktop-v${{ env.ONTSEQ_VERSION }}-win-x64-setup-engineering\n"
            ),
            "CHANGELOG": "## Unreleased\n## 9.8.7\n## 0.6.2\n",
            "README": guard.STATUS_HEADING + "\nDesktop/Core: 9.8.7\n",
        }
        for name, content in self.files.items():
            path = getattr(guard, name)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")

    def test_next_release_derives_resource_root_and_preserves_independent_versions(self) -> None:
        self.assertEqual(guard._desktop_resource_root(), "~/.local/share/ontseq/resources-v9.8.7")
        self.assertEqual(guard._problems(), [])

    def test_literal_duplicates_cannot_drift_from_shared_core_version(self) -> None:
        cases = (
            ("DESKTOP_VERSION", 'public const string Value = "9.8.7-engineering";', "Value"),
            (
                "DESKTOP_LAUNCHER",
                'private const string ReleaseVersion = "9.8.7";',
                "ReleaseVersion",
            ),
            (
                "DESKTOP_MODELS",
                "public const string DefaultResourceRootWsl = "
                '"~/.local/share/ontseq/resources-v9.8.7";',
                "DefaultResourceRootWsl",
            ),
        )
        for name, content, label in cases:
            with self.subTest(name=name):
                path = getattr(guard, name)
                if name == "DESKTOP_VERSION":
                    content += '\npublic const string CoreValue = "9.8.7";'
                path.write_text(content, encoding="utf-8")
                try:
                    with self.assertRaisesRegex(guard.Mismatch, label + " must derive"):
                        guard._problems()
                finally:
                    path.write_text(self.files[name], encoding="utf-8")

    def test_stale_resource_example_is_rejected(self) -> None:
        guard.DESKTOP_SETTINGS_EXAMPLE.write_text(
            '{"resourceRootWsl": "~/.local/share/ontseq/resources-v0.6.2"}', encoding="utf-8"
        )
        problems = guard._problems()
        self.assertEqual(len(problems), 1)
        self.assertIn("desktop.settings.example.json resourceRootWsl", problems[0])

    def test_stale_operator_resource_root_is_rejected(self) -> None:
        guard.OPERATOR_FIRST_RUN.write_text(
            "~/.local/share/ontseq/resources-v0.6.2", encoding="utf-8"
        )
        problems = guard._problems()
        self.assertEqual(len(problems), 1)
        self.assertIn("DESKTOP_FIRST_RUN.md", problems[0])


if __name__ == "__main__":
    unittest.main()
