"""Identity of explicitly qualified modkit builds; a version banner is not a patch proof."""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass
from pathlib import Path

PR709_BUILD_ID = "ontseq-modkit-pr709-v1"
PR709_SOURCE_COMMIT = "9fb9aea763aa1b78ac1172fa6d6734e7955ccc70"
# Executables added only after passing the recorded real-tool qualification matrix.
# No environment variable, policy file or adjacent receipt can grant this exception.
PR709_BINARY_SHA256: frozenset[str] = frozenset(
    {
        "f71542ee41ec964202a751919c3f47e0d54398ba0e5c5c1ccb0990c8a65769aa",
    }
)


@dataclass(frozen=True)
class ModkitBinaryIdentity:
    executable: str
    sha256: str | None
    build_id: str | None = None
    source_commit: str | None = None

    def parameters(self) -> dict[str, str | None]:
        return {
            "modkit_binary_sha256": self.sha256,
            "modkit_build_id": self.build_id,
            "modkit_source_commit": self.source_commit,
        }


def identify_modkit_binary(executable: str) -> ModkitBinaryIdentity:
    """Resolve exactly as a local subprocess does, then identify the executable's bytes."""
    located = shutil.which(executable)
    if located is None:
        # Test/remote CommandRunner implementations may not expose a local executable.
        # An unresolved tool never qualifies for the scientific-correctness exception.
        return ModkitBinaryIdentity(executable, None)
    path = Path(located).resolve(strict=True)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    sha256 = digest.hexdigest()
    reviewed = sha256 in PR709_BINARY_SHA256
    return ModkitBinaryIdentity(
        str(path),
        sha256,
        PR709_BUILD_ID if reviewed else None,
        PR709_SOURCE_COMMIT if reviewed else None,
    )
