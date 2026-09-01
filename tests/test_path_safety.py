from __future__ import annotations

import stat
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from ontseq_platform.path_safety import is_link_like


class PathSafetyTests(unittest.TestCase):
    def test_python311_windows_reparse_point_is_link_like(self) -> None:
        attributes = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
        fake_stat = SimpleNamespace(st_file_attributes=attributes)
        candidate = Path("junction")

        with (
            patch("ontseq_platform.path_safety.path_lexists", return_value=True),
            patch("ontseq_platform.path_safety.os.name", "nt"),
            patch.object(Path, "is_symlink", return_value=False),
            patch.object(Path, "lstat", return_value=fake_stat),
        ):
            self.assertTrue(is_link_like(candidate))

    def test_python311_windows_plain_directory_is_not_link_like(self) -> None:
        fake_stat = SimpleNamespace(st_file_attributes=0)
        candidate = Path("directory")

        with (
            patch("ontseq_platform.path_safety.path_lexists", return_value=True),
            patch("ontseq_platform.path_safety.os.name", "nt"),
            patch.object(Path, "is_symlink", return_value=False),
            patch.object(Path, "lstat", return_value=fake_stat),
        ):
            self.assertFalse(is_link_like(candidate))


if __name__ == "__main__":
    unittest.main()
