from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.service.app import _build_manifest
from ontseq_platform.service.guard import GuardError


def _payload(bam: Path) -> dict[str, str]:
    return {
        "bam": str(bam),
        "sample_id": "SAMPLE_001",
        "run_id": "RUN_001",
        "assay": "lcwgs",
        "genome_build": "GRCh38",
    }


class ServiceManifestIndexTests(unittest.TestCase):
    def test_manifest_uses_bam_dot_bai_when_both_exist(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            preferred = Path(f"{bam}.bai")
            alternative = bam.with_suffix(".bai")
            for path in (bam, preferred, alternative):
                path.write_bytes(b"x")

            manifest = _build_manifest(
                _payload(bam),
                reference_id="GRCh38_TEST",
                allowed_roots=[Path(temporary)],
            )

            self.assertEqual(manifest.input.index_path, str(preferred))

    def test_manifest_accepts_short_bai(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            alternative = bam.with_suffix(".bai")
            bam.write_bytes(b"BAM")
            alternative.write_bytes(b"BAI")

            manifest = _build_manifest(
                _payload(bam),
                reference_id="GRCh38_TEST",
                allowed_roots=[Path(temporary)],
            )

            self.assertEqual(manifest.input.index_path, str(alternative))

    def test_manifest_refuses_a_missing_index(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "sample.bam"
            bam.write_bytes(b"BAM")

            with self.assertRaises(GuardError):
                _build_manifest(
                    _payload(bam),
                    reference_id="GRCh38_TEST",
                    allowed_roots=[Path(temporary)],
                )

    def test_outside_path_is_refused_before_index_discovery(self) -> None:
        with (
            tempfile.TemporaryDirectory() as allowed,
            tempfile.TemporaryDirectory() as outside,
        ):
            bam = Path(outside) / "sample.bam"
            bam.write_bytes(b"BAM")

            with self.assertRaises(GuardError) as caught:
                _build_manifest(
                    _payload(bam),
                    reference_id="GRCh38_TEST",
                    allowed_roots=[Path(allowed)],
                )

            self.assertIn("outside the allowed directories", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
