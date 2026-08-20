from __future__ import annotations

import io
import json
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.models import AnalysisModule
from ontseq_platform.service.app import _build_manifest, _read_chunked_body
from ontseq_platform.service.guard import GuardError


def _payload(bam: Path, *, genome_build: str = "GRCh37") -> dict[str, str]:
    return {
        "bam": str(bam),
        "sample_id": "SYNTHETIC_001",
        "run_id": "DESKTOP_SYNTHETIC_001",
        "assay": "lcwgs",
        "genome_build": genome_build,
        "target_bed": "",
        "target_bed_version": "",
    }


class DesktopManifestTests(unittest.TestCase):
    def test_short_bai_is_recorded_and_grch37_requests_cnv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "sample.bam"
            bai = root / "sample.bai"
            bam.write_bytes(b"BAM")
            bai.write_bytes(b"BAI")
            manifest = _build_manifest(
                _payload(bam),
                reference_id="SYNTHETIC_REF",
                allowed_roots=[root],
            )
            self.assertEqual(manifest.input.index_path, str(bai))
            self.assertIn(AnalysisModule.CNV, manifest.analysis.modules)

    def test_bam_dot_bai_has_precedence_in_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "sample.bam"
            preferred = Path(f"{bam}.bai")
            alternative = bam.with_suffix(".bai")
            for path in (bam, preferred, alternative):
                path.write_bytes(b"x")
            manifest = _build_manifest(
                _payload(bam),
                reference_id="SYNTHETIC_REF",
                allowed_roots=[root],
            )
            self.assertEqual(manifest.input.index_path, str(preferred))

    def test_grch38_does_not_request_unbundled_qdnaseq_annotations(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "sample.bam"
            Path(f"{bam}.bai").write_bytes(b"BAI")
            bam.write_bytes(b"BAM")
            manifest = _build_manifest(
                _payload(bam, genome_build="GRCh38"),
                reference_id="SYNTHETIC_REF",
                allowed_roots=[root],
            )
            self.assertNotIn(AnalysisModule.CNV, manifest.analysis.modules)

    def test_bam_outside_allowed_root_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as outside:
            root = Path(temporary)
            bam = Path(outside) / "sample.bam"
            bam.write_bytes(b"BAM")
            Path(f"{bam}.bai").write_bytes(b"BAI")
            with self.assertRaises(GuardError):
                _build_manifest(
                    _payload(bam),
                    reference_id="SYNTHETIC_REF",
                    allowed_roots=[root],
                )


class ChunkedRequestTests(unittest.TestCase):
    def test_valid_chunked_json_is_reassembled(self) -> None:
        payload = json.dumps({"sample_id": "SYNTHETIC_001"}).encode()
        first = payload[:10]
        second = payload[10:]
        framed = (
            f"{len(first):X}\r\n".encode()
            + first
            + b"\r\n"
            + f"{len(second):X}\r\n".encode()
            + second
            + b"\r\n0\r\n\r\n"
        )
        self.assertEqual(_read_chunked_body(io.BytesIO(framed)), payload)

    def test_truncated_chunk_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "truncated"):
            _read_chunked_body(io.BytesIO(b"5\r\nabc\r\n0\r\n\r\n"))

    def test_invalid_chunk_size_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "chunk size"):
            _read_chunked_body(io.BytesIO(b"XYZ\r\nabc\r\n0\r\n\r\n"))


if __name__ == "__main__":
    unittest.main()
