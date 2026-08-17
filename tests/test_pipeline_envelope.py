from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.pipeline.envelope import (
    SUBDIRECTORIES,
    Artifact,
    EnvelopeError,
    RunEnvelope,
    canonical_signature,
    sha256_file,
    stage_signature,
)


class _EnvelopeCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.base = Path(self._tmp.name)
        self.envelope = RunEnvelope.create(self.base, run_id="RUN_001", sample_id="SAMPLE_001")

    def tearDown(self) -> None:
        self._tmp.cleanup()


class LayoutTests(_EnvelopeCase):
    def test_every_declared_subdirectory_exists(self) -> None:
        for name in SUBDIRECTORIES:
            self.assertTrue((self.envelope.root / name).is_dir(), name)

    def test_root_is_run_then_sample(self) -> None:
        self.assertEqual(self.envelope.root, self.base / "RUN_001" / "SAMPLE_001")

    def test_creating_twice_is_idempotent(self) -> None:
        again = RunEnvelope.create(self.base, run_id="RUN_001", sample_id="SAMPLE_001")
        self.assertEqual(again.root, self.envelope.root)

    def test_unsafe_identifiers_are_rejected(self) -> None:
        for run_id in ("", "..", "a/b"):
            with self.assertRaises(EnvelopeError):
                RunEnvelope.create(self.base, run_id=run_id, sample_id="SAMPLE_001")


class PathSafetyTests(_EnvelopeCase):
    """Reviewer artifacts must never carry a path from outside the envelope."""

    def test_traversal_is_refused(self) -> None:
        with self.assertRaises(EnvelopeError):
            self.envelope.path("../../etc/passwd")

    def test_absolute_paths_are_refused(self) -> None:
        with self.assertRaises(EnvelopeError):
            self.envelope.path("/etc/passwd")

    def test_relative_refuses_a_path_outside_the_envelope(self) -> None:
        with self.assertRaises(EnvelopeError):
            self.envelope.relative(Path("/etc/passwd"))

    def test_relative_returns_posix_form_inside(self) -> None:
        target = self.envelope.path("qc/report.json")
        target.write_text("{}", encoding="utf-8")
        self.assertEqual(self.envelope.relative(target), "qc/report.json")


class AtomicWriteTests(_EnvelopeCase):
    def test_write_produces_a_fingerprinted_artifact(self) -> None:
        artifact = self.envelope.atomic_write_text("qc/report.json", '{"a": 1}')
        self.assertEqual(artifact.relative_path, "qc/report.json")
        self.assertEqual(artifact.size_bytes, 8)
        self.assertEqual(artifact.sha256, sha256_file(self.envelope.path("qc/report.json")))

    def test_no_temporary_files_are_left_behind(self) -> None:
        self.envelope.atomic_write_text("qc/report.json", "{}")
        leftovers = [p.name for p in self.envelope.path("qc").iterdir() if p.name.startswith(".")]
        self.assertEqual(leftovers, [])

    def test_rewriting_replaces_atomically(self) -> None:
        self.envelope.atomic_write_text("qc/report.json", "old")
        artifact = self.envelope.atomic_write_text("qc/report.json", "new")
        self.assertEqual(self.envelope.path("qc/report.json").read_text(), "new")
        self.assertEqual(artifact.size_bytes, 3)

    def test_adopt_moves_an_external_file_into_the_envelope(self) -> None:
        source = self.base / "scratch.txt"
        source.write_text("payload", encoding="utf-8")
        artifact = self.envelope.adopt(source, "evidence/sv/candidates.vcf")
        self.assertFalse(source.exists())
        self.assertEqual(artifact.relative_path, "evidence/sv/candidates.vcf")
        self.assertEqual(
            self.envelope.path("evidence/sv/candidates.vcf").read_text(), "payload"
        )

    def test_intermediate_directories_are_not_exportable(self) -> None:
        exported = self.envelope.atomic_write_text("reports/report.html", "<html></html>")
        internal = self.envelope.atomic_write_text("work/tmp.bam", "x")
        alignment = self.envelope.atomic_write_text("alignment/sample.bam", "x")
        self.assertTrue(exported.exportable)
        self.assertFalse(internal.exportable)
        self.assertFalse(alignment.exportable)

    def test_raw_genomic_formats_are_never_exportable(self) -> None:
        """The rule that keeps a VCF out of Git must keep it out of a bundle too."""
        for relative_path in (
            "evidence/sv/candidates.vcf",
            "evidence/sv/candidates.vcf.gz",
            "qc/reads.fastq",
            "normalized/depth.bigwig",
        ):
            artifact = self.envelope.atomic_write_text(relative_path, "x")
            self.assertFalse(artifact.exportable, relative_path)

    def test_normalized_json_next_to_a_vcf_stays_exportable(self) -> None:
        artifact = self.envelope.atomic_write_text("evidence/sv/sniffles.json", "{}")
        self.assertTrue(artifact.exportable)

    def test_fingerprinting_a_missing_artifact_raises(self) -> None:
        with self.assertRaises(EnvelopeError):
            self.envelope.fingerprint("qc/absent.json")


class VerificationTests(_EnvelopeCase):
    def test_untouched_artifacts_verify(self) -> None:
        artifact = self.envelope.atomic_write_text("qc/report.json", "{}")
        self.assertEqual(self.envelope.verify([artifact]), [])

    def test_missing_artifact_is_reported(self) -> None:
        artifact = self.envelope.atomic_write_text("qc/report.json", "{}")
        self.envelope.path("qc/report.json").unlink()
        self.assertEqual(self.envelope.verify([artifact]), ["qc/report.json: missing"])

    def test_modified_artifact_is_reported(self) -> None:
        artifact = self.envelope.atomic_write_text("qc/report.json", "{}")
        self.envelope.path("qc/report.json").write_text("{ }", encoding="utf-8")
        self.assertEqual(self.envelope.verify([artifact]), ["qc/report.json: size changed"])

    def test_same_size_different_content_is_still_caught(self) -> None:
        artifact = self.envelope.atomic_write_text("qc/report.json", "ab")
        self.envelope.path("qc/report.json").write_text("cd", encoding="utf-8")
        self.assertEqual(self.envelope.verify([artifact]), ["qc/report.json: checksum changed"])


class SignatureTests(unittest.TestCase):
    """Resume must be content-addressed, so every input has to move the signature."""

    def _artifact(self, sha: str) -> Artifact:
        return Artifact(relative_path="qc/report.json", size_bytes=2, sha256=sha, exportable=True)

    def _signature(self, **overrides: object) -> str:
        kwargs: dict[str, object] = {
            "stage": "qc",
            "upstream": [self._artifact("a" * 64)],
            "parameters": {"threads": 4},
            "tool_versions": {"cramino": "1.3.0"},
            "external_inputs": [("input.bam", "b" * 64)],
        }
        kwargs.update(overrides)
        return stage_signature(**kwargs)  # type: ignore[arg-type]

    def test_identical_inputs_give_an_identical_signature(self) -> None:
        self.assertEqual(self._signature(), self._signature())

    def test_parameter_change_moves_the_signature(self) -> None:
        self.assertNotEqual(self._signature(), self._signature(parameters={"threads": 8}))

    def test_tool_version_change_moves_the_signature(self) -> None:
        self.assertNotEqual(
            self._signature(), self._signature(tool_versions={"cramino": "1.4.0"})
        )

    def test_upstream_content_change_moves_the_signature(self) -> None:
        self.assertNotEqual(self._signature(), self._signature(upstream=[self._artifact("c" * 64)]))

    def test_external_input_change_moves_the_signature(self) -> None:
        self.assertNotEqual(
            self._signature(), self._signature(external_inputs=[("input.bam", "d" * 64)])
        )

    def test_stage_identity_moves_the_signature(self) -> None:
        self.assertNotEqual(self._signature(), self._signature(stage="sv"))

    def test_key_order_does_not_matter(self) -> None:
        first = canonical_signature({"a": 1, "b": 2})
        second = canonical_signature({"b": 2, "a": 1})
        self.assertEqual(first, second)

    def test_signature_is_total_for_unusual_value_types(self) -> None:
        # A mid-run crash while hashing parameters would be the worst possible failure.
        self.assertTrue(canonical_signature({"path": Path("/tmp/x"), "when": object()}))


if __name__ == "__main__":
    unittest.main()
