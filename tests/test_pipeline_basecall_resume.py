from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ontseq_platform.pipeline.envelope import stage_signature
from ontseq_platform.pipeline.input_digest import RunInputDigestCache
from ontseq_platform.pipeline.runner import StageFailure, _pod5_external_inputs


class _DigestContext:
    def __init__(self) -> None:
        self.input_digests = RunInputDigestCache()

    def fingerprint_external_input(
        self, path: Path, *, label: str | None = None
    ) -> tuple[str, str]:
        digest, stable = self.input_digests.digest(path)
        if not stable:
            raise StageFailure("required external input changed while it was being fingerprinted")
        return (label or path.name, digest)


class _MutatingContext(_DigestContext):
    def __init__(self, directory: Path) -> None:
        super().__init__()
        self.directory = directory
        self.mutated = False

    def fingerprint_external_input(
        self, path: Path, *, label: str | None = None
    ) -> tuple[str, str]:
        result = super().fingerprint_external_input(path, label=label)
        if not self.mutated:
            (self.directory / "late.pod5").write_bytes(b"late synthetic signal")
            self.mutated = True
        return result


def _signature(external_inputs: tuple[tuple[str, str], ...]) -> str:
    return stage_signature(
        stage="basecall",
        upstream=(),
        parameters={"model": "synthetic-model"},
        tool_versions={"dorado": "0.9.0"},
        external_inputs=external_inputs,
    )


class Pod5ResumeInputTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name) / "pod5"
        (self.directory / "nested").mkdir(parents=True)
        self.first = self.directory / "a.pod5"
        self.second = self.directory / "nested" / "b.pod5"
        self.first.write_bytes(b"AAAAAA")
        self.second.write_bytes(b"BBBBBB")

    def test_unchanged_inventory_has_stable_signature_and_relative_labels(self) -> None:
        context = _DigestContext()
        before = _pod5_external_inputs(context, self.directory)
        after = _pod5_external_inputs(context, self.directory)

        self.assertEqual(before, after)
        self.assertEqual(
            [label for label, _digest in before],
            ["pod5:a.pod5", "pod5:nested/b.pod5"],
        )
        self.assertEqual(_signature(before), _signature(after))

    def test_content_change_invalidates_basecall_resume_signature(self) -> None:
        context = _DigestContext()
        before = _pod5_external_inputs(context, self.directory)
        self.first.write_bytes(b"CCCCCC")
        after = _pod5_external_inputs(context, self.directory)

        self.assertNotEqual(dict(before)["pod5:a.pod5"], dict(after)["pod5:a.pod5"])
        self.assertNotEqual(_signature(before), _signature(after))

    def test_membership_change_invalidates_basecall_resume_signature(self) -> None:
        context = _DigestContext()
        before = _pod5_external_inputs(context, self.directory)
        (self.directory / "new.pod5").write_bytes(b"new synthetic signal")
        after = _pod5_external_inputs(context, self.directory)

        self.assertNotEqual(_signature(before), _signature(after))
        self.assertIn("pod5:new.pod5", dict(after))

    def test_inventory_drift_during_fingerprinting_fails_closed(self) -> None:
        with self.assertRaisesRegex(StageFailure, "input set changed"):
            _pod5_external_inputs(_MutatingContext(self.directory), self.directory)

    def test_empty_directory_is_rejected(self) -> None:
        empty = Path(self.temp.name) / "empty"
        empty.mkdir()
        with self.assertRaisesRegex(StageFailure, "contains no .pod5"):
            _pod5_external_inputs(_DigestContext(), empty)


if __name__ == "__main__":
    unittest.main()
