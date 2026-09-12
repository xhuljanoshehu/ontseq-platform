"""Real filesystem capability gates for integration tests on restricted hosts."""

from __future__ import annotations

import errno
import unittest
from pathlib import Path


def symlink_or_skip(case: unittest.TestCase, link: Path, target: Path) -> None:
    """Keep real symlink tests enabled wherever the host grants this capability.

    Windows requires a privilege or Developer Mode even for temporary synthetic files.
    Other filesystem errors must still fail the test instead of being hidden.
    """
    try:
        link.symlink_to(target)
    except OSError as error:
        if getattr(error, "winerror", None) == 1314 or error.errno in {
            errno.ENOSYS,
            errno.ENOTSUP,
        }:
            case.skipTest(f"real symlinks unavailable on this host: {error}")
        raise
