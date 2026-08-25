from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.pipeline import envelope as envelope_module
from ontseq_platform.pipeline.envelope import (
    RunEnvelope,
    forget_input_digests,
    sha256_file,
    sha256_input,
)


class _DigestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.base = Path(self._tmp.name)
        forget_input_digests()

    def tearDown(self) -> None:
        forget_input_digests()
        self._tmp.cleanup()

    def _write(self, name: str, content: bytes) -> Path:
        path = self.base / name
        path.write_bytes(content)
        return path


class InputDigestTests(_DigestCase):
    def test_the_digest_is_the_one_sha256_file_would_have_recorded(self) -> None:
        """Provenance must not change. This only stops the number being recomputed."""
        path = self._write("sample.bam", b"aligned reads")
        digest, stable = sha256_input(path)
        self.assertEqual(digest, sha256_file(path))
        self.assertTrue(stable)

    def test_a_file_is_read_once_however_many_stages_ask_for_it(self) -> None:
        path = self._write("sample.bam", b"aligned reads" * 1000)
        reads = 0
        original = envelope_module.sha256_file

        def counting(target: Path, **kwargs: object) -> str:
            nonlocal reads
            reads += 1
            return original(target, **kwargs)  # type: ignore[arg-type]

        envelope_module.sha256_file = counting  # type: ignore[assignment]
        try:
            digests = {sha256_input(path)[0] for _ in range(6)}
        finally:
            envelope_module.sha256_file = original  # type: ignore[assignment]
        self.assertEqual(reads, 1, "six stage plans should cost one read, not six")
        self.assertEqual(len(digests), 1)

    def test_two_names_for_one_file_share_the_read(self) -> None:
        path = self._write("sample.bam", b"aligned reads")
        link = self.base / "hardlink.bam"
        os.link(path, link)
        self.assertEqual(sha256_input(path)[0], sha256_input(link)[0])

    def test_a_replaced_file_is_read_again(self) -> None:
        path = self._write("sample.bam", b"first")
        first = sha256_input(path)[0]
        path.write_bytes(b"second and longer")
        self.assertNotEqual(sha256_input(path)[0], first)
        self.assertEqual(sha256_input(path)[0], sha256_file(path))

    def test_different_files_do_not_share_a_digest(self) -> None:
        one = self._write("a.bam", b"one")
        two = self._write("b.bam", b"two")
        self.assertNotEqual(sha256_input(one)[0], sha256_input(two)[0])

    def test_the_map_stays_bounded_for_a_long_lived_daemon(self) -> None:
        limit = envelope_module._INPUT_DIGEST_LIMIT
        for index in range(limit + 25):
            sha256_input(self._write(f"s{index}.bam", f"sample {index}".encode()))
        self.assertLessEqual(len(envelope_module._INPUT_DIGESTS), limit)


class UnstableReadTests(_DigestCase):
    """A file that moves while it is being read is reported, and is not remembered.

    The stability check came from the intake hardening on `main`; memoising must not
    weaken it. An unstable read produces a digest of something that no longer exists, so
    caching it would let a later stage fingerprint a state the file was never in.
    """

    def _read_and_disturb(self, path: Path, content: bytes) -> None:
        original = envelope_module.sha256_file

        def disturbing(target: Path, **kwargs: object) -> str:
            digest = original(target, **kwargs)  # type: ignore[arg-type]
            target.write_bytes(content)  # the file moves under the reader
            return digest

        envelope_module.sha256_file = disturbing  # type: ignore[assignment]
        try:
            self.last = sha256_input(path)
        finally:
            envelope_module.sha256_file = original  # type: ignore[assignment]

    def test_an_unstable_read_is_reported_as_unstable(self) -> None:
        path = self._write("sample.bam", b"first")
        self._read_and_disturb(path, b"changed while being read")
        self.assertFalse(self.last[1])

    def test_an_unstable_read_is_not_remembered(self) -> None:
        path = self._write("sample.bam", b"first")
        self._read_and_disturb(path, b"changed while being read")
        self.assertEqual(envelope_module._INPUT_DIGESTS, {})
        # The next caller sees the file as it now is, not the state that was mid-read.
        self.assertEqual(sha256_input(path)[0], sha256_file(path))


class VerificationIsNotMemoisedTests(_DigestCase):
    """The reason `sha256_input` is a separate function rather than a faster `sha256_file`.

    Two callers re-read a file precisely to notice that it changed: `verify()`, against an
    artifact's recorded checksum, and the intake stage, against the digest its own plan
    recorded. Answering either from a map populated when the file was first read would
    reduce the check to comparing a value with itself. Both call `sha256_file`.
    """

    def test_verify_still_catches_a_change_the_cache_cannot_see(self) -> None:
        envelope = RunEnvelope.create(self.base, run_id="RUN_001", sample_id="SAMPLE_001")
        artifact = envelope.atomic_write_text("qc/report.json", '{"verdict": "PASS"}')
        target = envelope.path(artifact.relative_path)

        # Populate the map for this exact file, then tamper in the way a naive cache cannot
        # see: identical byte count, and the modification time put back afterwards.
        self.assertEqual(sha256_input(target)[0], artifact.sha256)
        stat = target.stat()
        target.write_text('{"verdict": "FAIL"}', encoding="utf-8")
        os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        self.assertEqual(target.stat().st_size, artifact.size_bytes, "test needs equal sizes")

        # The memoised digest is stale by construction — that is what makes the point.
        self.assertEqual(sha256_input(target)[0], artifact.sha256)
        # Reading the file is not fooled, and verification reads the file.
        self.assertNotEqual(sha256_file(target), artifact.sha256)
        self.assertEqual(
            envelope.verify([artifact]),
            [f"{artifact.relative_path}: checksum changed"],
        )


if __name__ == "__main__":
    unittest.main()
