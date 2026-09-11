"""Run-scoped content digests for large immutable analysis inputs.

Stage plans are rebuilt during every invocation so content-addressed resume can decide
whether a stage is reusable. Several plans fingerprint the same multi-gigabyte BAM. This
cache removes duplicate reads *within one RunContext* while preserving the deliberate
re-verification reads used by intake and release verification.
"""

from __future__ import annotations

from pathlib import Path

from .envelope import sha256_file

FileIdentity = tuple[str, int, int, int, int, int]


class RunInputDigestCache:
    """Memoize stable input digests for one pipeline run invocation only."""

    def __init__(self) -> None:
        self._digests: dict[FileIdentity, str] = {}

    @staticmethod
    def _identity(path: Path) -> FileIdentity:
        resolved = path.resolve(strict=True)
        stat = resolved.stat()
        return (
            str(resolved),
            int(stat.st_dev),
            int(stat.st_ino),
            int(stat.st_size),
            int(stat.st_mtime_ns),
            int(stat.st_ctime_ns),
        )

    def digest(self, path: Path) -> tuple[str, bool]:
        """Return SHA-256 plus whether the file stayed unchanged while it was read.

        The cache key includes the resolved path and filesystem identity/state. A changed
        file therefore misses the cache. A digest is inserted only after a second stat
        proves that the same target, size, mtime and ctime survived the full read.
        """

        before = self._identity(path)
        cached = self._digests.get(before)
        if cached is not None:
            return cached, True

        digest = sha256_file(Path(before[0]))
        try:
            after = self._identity(path)
        except OSError:
            return digest, False
        stable = before == after
        if stable:
            self._digests[before] = digest
        return digest, stable
