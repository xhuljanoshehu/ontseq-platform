from __future__ import annotations

import gzip
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import check_repository_safety as safety


class RepositorySafetyTests(unittest.TestCase):
    def test_renamed_nanopolish_inputs_and_reports_are_detected(self) -> None:
        # All payloads are synthetic; they live only in a temporary directory.
        artifacts = {
            "calls.tsv": (
                "chromosome\tstart\tend\tread_name\tlog_lik_ratio\t"
                "log_lik_methylated\tlog_lik_unmethylated\n"
            ),
            "renamed.json": json.dumps({"sensitive_output": True}),
            "renamed.csv": '"source_a_id","source_b_id","estimated_source_a_fraction"\n',
            "renamed.html": "<!doctype html><title>Methylation source-fraction experiment</title>",
            "per_read.tsv": "read_id\tref_position\tmod_qual\tmod_code\n",
            "validation.csv": ("read_group_budget,replicate_seed,estimated_source_a_fraction\n"),
            "validation.html": '<meta name="ontseq-sensitive-output" content="true">',
        }
        with tempfile.TemporaryDirectory() as directory:
            for name, content in artifacts.items():
                for compressed in (False, True):
                    with self.subTest(name=name, compressed=compressed):
                        path = Path(directory) / (name + (".gz" if compressed else ""))
                        raw = ("\ufeff" + content).encode("utf-8")
                        path.write_bytes(gzip.compress(raw) if compressed else raw)
                        self.assertIsNotNone(safety.sensitive_content_reason(path))

    def test_legitimate_configuration_and_schemas_remain_allowed(self) -> None:
        paths = [
            Path("configs/qc/target_coverage_expectations.grch38.tsv"),
            Path("schemas/methylation-mixture-report.schema.json"),
            Path("schemas/methylation-mixture-policy.schema.json"),
        ]
        for path in paths:
            with self.subTest(path=path):
                self.assertTrue(path.is_file())
                self.assertIsNone(safety.sensitive_content_reason(path))
                ignored = subprocess.run(
                    ["git", "check-ignore", "--no-index", str(path)],
                    capture_output=True,
                    check=False,
                )
                self.assertEqual(ignored.returncode, 1, ignored.stdout)

    def test_default_artifact_names_are_ignored(self) -> None:
        for name in (
            "sample.nanopolish.tsv",
            "sample.nanopolish.tsv.gz",
            "sample.methylation_calls.tsv",
            "sample.methylation-mixture.json",
            "sample.methylation-mixture.csv",
            "sample.methylation-mixture.html",
            "sample.modkit.tsv",
            "sample.modkit.tsv.gz",
            "sample.methylation-validation.json",
            "sample.methylation-validation.csv",
            "sample.methylation-validation.html",
            "sample.methylation-holdout.json",
            "sample.methylation-holdout.csv",
            "sample.methylation-holdout.html",
        ):
            with self.subTest(name=name):
                result = subprocess.run(
                    ["git", "check-ignore", "--no-index", name],
                    capture_output=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_forced_tracked_files_are_still_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.methylation-mixture.json"
            path.write_text("{}", encoding="utf-8")
            with (
                patch.object(safety, "_candidate_files", return_value=[path]),
                self.assertRaisesRegex(SystemExit, "prohibited methylation artifact filename"),
            ):
                safety.main()

    def test_unreadable_gzip_and_oversized_expansion_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.tsv.gz"
            path.write_bytes(b"not a gzip file")
            self.assertEqual(
                safety.sensitive_content_reason(path), "cannot inspect structured data file"
            )
            path.write_bytes(gzip.compress(b"x" * 101))
            with patch.object(safety, "MAX_BYTES", 100):
                self.assertEqual(
                    safety.sensitive_content_reason(path),
                    "structured data exceeds 5 MiB inspection limit",
                )

    def test_git_paths_with_newlines_are_not_split(self) -> None:
        completed = subprocess.CompletedProcess([], 0, "a\nfile.tsv\0other.json\0", "")
        with patch.object(safety.subprocess, "run", return_value=completed):
            self.assertEqual(safety._candidate_files(), [Path("a\nfile.tsv"), Path("other.json")])


if __name__ == "__main__":
    unittest.main()
