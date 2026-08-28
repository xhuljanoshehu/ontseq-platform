from __future__ import annotations

import gzip
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ontseq_platform.sidecars import tabular_sidecar


class SidecarTests(unittest.TestCase):
    def test_tsv_metadata_is_streamed_and_relative(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            path = root / "sidecars" / "cnv_segments.tsv"
            path.parent.mkdir()
            path.write_text("chrom\tstart\tend\nchr1\t0\t10\nchr2\t0\t20\n", encoding="utf-8")
            artifact = tabular_sidecar(
                artifact_id="cnv_segments",
                envelope_root=root,
                relative_path="sidecars/cnv_segments.tsv",
            )
            self.assertEqual(artifact.row_count, 2)
            self.assertEqual(artifact.columns, ["chrom", "start", "end"])
            self.assertEqual(artifact.size_bytes, path.stat().st_size)

    def test_gzip_table_is_supported(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            path = root / "sidecars.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write("value\n1\n")
            artifact = tabular_sidecar(
                artifact_id="read_length_histogram",
                envelope_root=root,
                relative_path="sidecars.tsv.gz",
            )
            self.assertEqual(artifact.row_count, 1)
            self.assertEqual(artifact.media_type, "application/gzip")


if __name__ == "__main__":
    unittest.main()
