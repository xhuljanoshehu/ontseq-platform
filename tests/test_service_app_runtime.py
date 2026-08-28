from __future__ import annotations

import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ontseq_platform.models import AnalysisModule
from ontseq_platform.pipeline.components import RunComponents
from ontseq_platform.service.app import (
    ServiceConfig,
    _build_manifest,
    _build_profile_configuration,
    _read_chunked_body,
)
from ontseq_platform.service.guard import GuardError

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs"


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

    def test_grch38_requests_real_tool_tested_qdnaseq_cnv(self) -> None:
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
            self.assertIn(AnalysisModule.CNV, manifest.analysis.modules)

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

    def test_profile_post_selects_fast_pinned_resource_preflight_and_core_run_id(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "sample.bam"
            bam.write_bytes(b"BAM")
            Path(f"{bam}.bai").write_bytes(b"BAI")
            resources = root / "resources"
            resources.mkdir()
            config = ServiceConfig(
                reference_lock=None,
                output_dir=root / "results",
                allowed_roots=[root],
                qc_policy=CONFIGS / "qc" / "defaults.yaml",
                sniffles_policy=CONFIGS / "sv" / "sniffles2.conservative.technical.yaml",
                target_coverage_policy=(CONFIGS / "qc" / "adaptive_target_coverage.technical.yaml"),
                cutesv_policy=CONFIGS / "sv" / "cutesv.conservative.technical.yaml",
                sv_consensus_policy=(CONFIGS / "sv" / "sniffles2_cutesv.consensus.technical.yaml"),
                sv_evidence_policy=CONFIGS / "sv" / "evidence-priority.technical.yaml",
                sv_minimum_mean_depth=23.5,
                components=RunComponents.model_validate(
                    {
                        "selection_id": "desktop-service-selection",
                        "status": "technical_defaults_only",
                        "components": {"cnv": {"provider": "qdnaseq_ace", "enabled": False}},
                    }
                ),
                resource_root=resources,
            )
            sentinel = object()

            with (
                patch("ontseq_platform.service.app.ResourceRegistry") as registry_type,
                patch(
                    "ontseq_platform.service.app.build_profile_run_configuration",
                    return_value=sentinel,
                ) as build,
                patch(
                    "ontseq_platform.service.app.windows_to_wsl",
                    side_effect=lambda value: value,
                ),
            ):
                registry_type.return_value.profiles = {"AML_LCWGS_GRCh38": object()}
                result = _build_profile_configuration(
                    {
                        "bam": str(bam),
                        "sample_id": "SAMPLE_001",
                        "profile": "AML_LCWGS_GRCh38",
                        "genome_build": "GRCh38",
                        "assay": "lcwgs",
                    },
                    config,
                )

            self.assertIs(result, sentinel)
            settings = build.call_args.args[0]
            self.assertFalse(settings.verify_resource_checksums)
            self.assertIsNone(settings.run_id)
            runtime = settings.runtime_settings
            self.assertIsNotNone(runtime)
            assert runtime is not None
            self.assertIn("Numeric clinical gates remain null", runtime.qc_policy.note)
            self.assertEqual(
                runtime.sniffles_policy.profile_id,
                "sniffles2-conservative-technical-v1",
            )
            self.assertEqual(
                runtime.target_coverage_policy.profile_id,
                "adaptive_target_coverage_technical_v1",
            )
            self.assertEqual(runtime.cutesv_policy.profile_id, "cutesv-conservative-technical-v1")
            self.assertEqual(
                runtime.sv_consensus_policy.profile_id,
                "sniffles2-cutesv-consensus-technical-v1",
            )
            self.assertEqual(
                runtime.sv_evidence_policy.profile_id,
                "sv-evidence-priority-technical-v1",
            )
            self.assertEqual(runtime.sv_minimum_mean_depth, 23.5)
            self.assertIs(runtime.components, config.components)


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
