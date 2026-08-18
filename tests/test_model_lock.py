from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.model_lock import (
    ModelLockError,
    exit_code,
    fingerprint,
    human_size,
    render,
)


def _historic_signature(path: Path) -> str:
    """The digest ``basecall.model_signature`` computed before it was extracted.

    Reimplemented rather than imported, deliberately. Importing the current implementation
    would only prove it equals itself; this pins the value against the algorithm that
    already produced every ``model_sha256`` a site may have recorded. If the two ever
    disagree, previously locked models start failing preflight for no stated reason.
    """
    digest = hashlib.sha256()
    for item in sorted(path.rglob("*")):
        if item.is_file():
            digest.update(item.relative_to(path).as_posix().encode("utf-8"))
            digest.update(hashlib.sha256(item.read_bytes()).hexdigest().encode("ascii"))
    return digest.hexdigest()


def _model(root: Path, files: dict[str, bytes] | None = None) -> Path:
    directory = root / "model"
    contents = {"config.toml": b"[model]\n", "weights/0.tensor": b"\x01\x02\x03"}
    for name, payload in (files if files is not None else contents).items():
        target = directory / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


class DigestTests(unittest.TestCase):
    def test_the_digest_matches_the_one_already_recorded_in_policies(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary))
            self.assertEqual(fingerprint(model).signature, _historic_signature(model))

    def test_the_digest_is_stable_across_calls(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary))
            self.assertEqual(fingerprint(model).signature, fingerprint(model).signature)

    def test_changing_one_byte_changes_the_digest(self) -> None:
        """Otherwise the lock detects nothing it exists to detect."""
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary))
            before = fingerprint(model).signature
            (model / "config.toml").write_bytes(b"[model]\n\n")
            self.assertNotEqual(fingerprint(model).signature, before)

    def test_renaming_a_file_changes_the_digest(self) -> None:
        """Names enter the digest, so a rearranged model is a different model."""
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary))
            before = fingerprint(model).signature
            (model / "config.toml").rename(model / "config2.toml")
            self.assertNotEqual(fingerprint(model).signature, before)

    def test_the_containing_directory_name_does_not_enter_the_digest(self) -> None:
        """A site that moves its models must not have to re-lock every policy."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = _model(root / "a")
            second = _model(root / "b")
            self.assertEqual(fingerprint(first).signature, fingerprint(second).signature)


class RefusalTests(unittest.TestCase):
    def test_a_missing_path_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, self.assertRaises(ModelLockError):
            fingerprint(Path(temporary) / "absent")

    def test_a_bare_model_name_is_refused_rather_than_resolved(self) -> None:
        """Guessing at Dorado's cache would produce a lock on a directory nobody chose."""
        with self.assertRaises(ModelLockError):
            fingerprint(Path("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"))

    def test_a_file_is_refused_with_a_message_naming_the_reason(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "model.tar"
            path.write_bytes(b"x")
            with self.assertRaises(ModelLockError) as caught:
                fingerprint(path)
            self.assertIn("not a directory", str(caught.exception))


class ConcernTests(unittest.TestCase):
    """Every directory yields a valid-looking digest, including a broken one."""

    def test_a_healthy_model_has_no_concerns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = fingerprint(_model(Path(temporary)))
            self.assertEqual(result.concerns, ())
            self.assertEqual(exit_code(result), 0)

    def test_an_empty_directory_is_refused_as_a_lock_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            empty = Path(temporary) / "model"
            empty.mkdir()
            result = fingerprint(empty)
            self.assertEqual(result.file_count, 0)
            self.assertTrue(result.concerns)
            self.assertEqual(exit_code(result), 2)

    def test_an_empty_file_is_named(self) -> None:
        """The usual shape of an interrupted download, and it hashes perfectly happily."""
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary), {"config.toml": b"x", "weights/0.tensor": b""})
            result = fingerprint(model)
            self.assertEqual(result.empty_files, ("weights/0.tensor",))
            self.assertIn("weights/0.tensor", " ".join(result.concerns))
            self.assertEqual(exit_code(result), 2)

    def test_a_broken_symlink_is_reported_rather_than_skipped(self) -> None:
        """It contributes nothing to the digest, so the digest cannot reveal it."""
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary))
            (model / "dangling.tensor").symlink_to(model / "does-not-exist")
            result = fingerprint(model)
            self.assertEqual(result.broken_links, ("dangling.tensor",))
            self.assertEqual(exit_code(result), 2)

    def test_a_broken_symlink_does_not_change_the_digest(self) -> None:
        """It is reported, not silently folded in: the bytes on disk are unchanged."""
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary))
            before = fingerprint(model).signature
            (model / "dangling.tensor").symlink_to(model / "does-not-exist")
            self.assertEqual(fingerprint(model).signature, before)

    def test_a_symlink_to_a_real_file_is_hashed_like_any_other(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary))
            before = fingerprint(model).signature
            (model / "alias.toml").symlink_to(model / "config.toml")
            self.assertNotEqual(fingerprint(model).signature, before)


class AccountingTests(unittest.TestCase):
    def test_the_file_count_and_size_describe_what_was_hashed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model = _model(Path(temporary), {"a": b"xx", "b/c": b"yyy"})
            result = fingerprint(model)
            self.assertEqual(result.file_count, 2)
            self.assertEqual(result.total_bytes, 5)
            self.assertEqual([item.relative_path for item in result.files], ["a", "b/c"])

    def test_human_size_keeps_the_exact_figure_beside_the_rounded_one(self) -> None:
        self.assertEqual(human_size(0), "0 B (0 bytes)")
        self.assertIn("(1536 bytes)", human_size(1536))
        self.assertIn("1.5 KiB", human_size(1536))
        self.assertIn("GiB", human_size(2 * 1024**3))


class RenderTests(unittest.TestCase):
    def test_a_healthy_model_renders_both_lines_to_record(self) -> None:
        """The path and the checksum together: one without the other locks nothing."""
        with tempfile.TemporaryDirectory() as temporary:
            result = fingerprint(_model(Path(temporary)))
            text = render(result)
            self.assertIn(f"model_sha256: {result.signature}", text)
            self.assertIn(f"model: {result.path}", text)

    def test_a_broken_model_is_not_offered_as_something_to_paste(self) -> None:
        """Printing the line beside the warning is how the warning gets ignored."""
        with tempfile.TemporaryDirectory() as temporary:
            empty = Path(temporary) / "model"
            empty.mkdir()
            text = render(fingerprint(empty))
            self.assertIn("DO NOT LOCK", text)
            self.assertNotIn("model_sha256:", text)

    def test_the_file_listing_appears_only_when_asked(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = fingerprint(_model(Path(temporary)))
            self.assertNotIn("config.toml", render(result))
            self.assertIn("config.toml", render(result, list_files=True))


if __name__ == "__main__":
    unittest.main()
