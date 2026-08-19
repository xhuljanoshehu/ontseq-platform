import tempfile
import unittest
from pathlib import Path

from pydantic import ValidationError

from ontseq_platform.desktop import (
    DesktopAnalysisRequest,
    DesktopBackend,
    DesktopConfig,
    DesktopReferenceProfile,
    load_desktop_config,
    locate_bam_index,
    sanitize_sample_id,
    save_desktop_config,
)
from ontseq_platform.models import AssayMode, GenomeBuild


class DesktopTests(unittest.TestCase):
    def test_sanitize_sample_id_removes_unsafe_characters(self) -> None:
        self.assertEqual(sanitize_sample_id(" AML Patient / 01 "), "AML_Patient_01")
        self.assertEqual(sanitize_sample_id("x"), "S_x")
        self.assertEqual(len(sanitize_sample_id("A" * 100)), 64)

    def test_locate_bam_index_supports_bam_bai_and_bai(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bam = root / "sample.bam"
            bam.write_bytes(b"BAM")
            bam_bai = root / "sample.bam.bai"
            bam_bai.write_bytes(b"BAI")
            self.assertEqual(locate_bam_index(bam), bam_bai)

            bam_bai.unlink()
            short_bai = root / "sample.bai"
            short_bai.write_bytes(b"BAI")
            self.assertEqual(locate_bam_index(bam), short_bai)

    def test_desktop_config_roundtrip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "desktop.json"
            config = DesktopConfig(
                backend_mode="wsl",
                wsl_distribution="Ubuntu",
                output_root=str(root / "results"),
                reference_profiles=[
                    DesktopReferenceProfile(
                        genome_build=GenomeBuild.GRCH38,
                        reference_id="GRCh38-test",
                        reference_lock_path=str(root / "grch38.lock.json"),
                    )
                ],
            )
            self.assertEqual(save_desktop_config(config, path), path)
            self.assertEqual(load_desktop_config(path), config)

    def test_duplicate_reference_builds_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = DesktopReferenceProfile(
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_lock_path=str(root / "lock.json"),
            )
            with self.assertRaises(ValidationError):
                DesktopConfig(reference_profiles=[profile, profile])

    def test_adaptive_sampling_profile_requires_bed_path_and_version_together(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValidationError):
                DesktopReferenceProfile(
                    genome_build=GenomeBuild.GRCH38,
                    reference_id="GRCh38-test",
                    reference_lock_path=str(root / "lock.json"),
                    adaptive_sampling_target_bed_path=str(root / "targets.bed"),
                )

    def test_manifest_is_pseudonymized_and_requests_research_modules(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bam = root / "SAMPLE_001.bam"
            bai = root / "SAMPLE_001.bam.bai"
            bed = root / "targets.bed"
            bam.write_bytes(b"BAM")
            bai.write_bytes(b"BAI")
            bed.write_text("chr1\t1\t100\n", encoding="utf-8")
            profile = DesktopReferenceProfile(
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_lock_path=str(root / "lock.json"),
                adaptive_sampling_target_bed_path=str(bed),
                adaptive_sampling_target_bed_version="test-v1",
            )
            backend = DesktopBackend(
                DesktopConfig(
                    backend_mode="local",
                    wsl_project_root=str(root),
                    output_root=str(root / "results"),
                    reference_profiles=[profile],
                )
            )
            request = DesktopAnalysisRequest(
                bam_path=bam,
                sample_id="SAMPLE_001",
                genome_build=GenomeBuild.GRCH38,
                assay_mode=AssayMode.ADAPTIVE_SAMPLING,
                output_dir=root / "results" / "SAMPLE_001",
                run_id="RUN_TEST_001",
            )

            manifest = backend._build_manifest(
                request,
                sample_id="SAMPLE_001",
                bam=bam,
                bai=bai,
                profile=profile,
            )

            self.assertEqual(manifest["sample_id"], "SAMPLE_001")
            self.assertEqual(
                manifest["privacy"],
                {
                    "pseudonymized": True,
                    "contains_direct_identifiers": False,
                    "cloud_upload_approved": False,
                },
            )
            analysis = manifest["analysis"]
            self.assertIsInstance(analysis, dict)
            assert isinstance(analysis, dict)
            self.assertEqual(
                analysis["modules"], ["qc", "cnv", "sv", "fusion", "iscn", "report"]
            )
            assay = manifest["assay"]
            self.assertIsInstance(assay, dict)
            assert isinstance(assay, dict)
            self.assertEqual(assay["target_bed_version"], "test-v1")

    def test_project_root_expression_supports_home_shortcut(self) -> None:
        self.assertEqual(DesktopBackend._project_root_expression("~"), '"$HOME"')
        self.assertTrue(
            DesktopBackend._project_root_expression("~/ontseq-platform").startswith('"$HOME"/')
        )


if __name__ == "__main__":
    unittest.main()
