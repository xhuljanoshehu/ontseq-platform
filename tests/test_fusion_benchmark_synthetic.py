from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.fusion_benchmark import load_synthetic_fusion_benchmark
from ontseq_platform.fusion_benchmark_synthetic import execute_synthetic_fusion_fixture


class SyntheticFusionBenchmarkExecutionTests(unittest.TestCase):
    def test_all_repository_fixtures_execute_through_real_software_path(self) -> None:
        suite = load_synthetic_fusion_benchmark(
            Path("configs/sv/fusion.synthetic_benchmark.yaml")
        )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            results = []
            for fixture in suite.fixtures:
                with self.subTest(fixture_id=fixture.fixture_id):
                    result = execute_synthetic_fusion_fixture(
                        fixture,
                        root / fixture.fixture_id,
                    )
                    results.append(result)
                    self.assertTrue(result.benchmark.passed, result.benchmark.failures)
                    self.assertEqual(
                        result.fusion_candidate_count,
                        fixture.expected_candidate_count,
                    )
                    self.assertEqual(
                        result.normalized_event_count,
                        fixture.expected_candidate_count,
                    )
                    self.assertEqual(
                        result.privacy_profile,
                        "synthetic_nonbiological_local_files_only",
                    )

        self.assertEqual(len(results), len(suite.fixtures))

    def test_generated_files_are_explicitly_synthetic_and_local(self) -> None:
        suite = load_synthetic_fusion_benchmark(
            Path("configs/sv/fusion.synthetic_benchmark.yaml")
        )
        fixture = suite.fixtures[0]

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = execute_synthetic_fusion_fixture(fixture, root)
            generated_names = sorted(path.name for path in root.iterdir())
            vcf_text = (root / f"{fixture.fixture_id}.synthetic.vcf").read_text(
                encoding="utf-8"
            )
            bed_text = (root / f"{fixture.fixture_id}.synthetic.genes.bed").read_text(
                encoding="utf-8"
            )

        self.assertTrue(result.benchmark.passed)
        self.assertEqual(
            generated_names,
            [
                f"{fixture.fixture_id}.synthetic.genes.bed",
                f"{fixture.fixture_id}.synthetic.vcf",
            ],
        )
        self.assertIn("SYNTHETIC_BND", vcf_text)
        self.assertIn("ENST_SYNTH_", bed_text)
        self.assertNotIn("PRIVATE_VCF_ID", vcf_text)


if __name__ == "__main__":
    unittest.main()
