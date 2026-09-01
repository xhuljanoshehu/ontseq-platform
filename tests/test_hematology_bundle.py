from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
BUNDLE_ROOT = REPO_ROOT / "configs" / "knowledge_bundles" / "HEMATOLOGY_v3"
MANIFEST = BUNDLE_ROOT / "bundle.yaml"


class HematologyBundleTests(unittest.TestCase):
    def test_existing_aml_knowledge_is_integrated_with_the_grch38_profile(self) -> None:
        manifest = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(manifest["bundle_id"], "HEMATOLOGY_v3")
        self.assertTrue(manifest["coordinate_bearing"])
        self.assertEqual(manifest["genome_build"], "GRCh38")
        for resource in manifest["resources"]:
            path = BUNDLE_ROOT / resource["path"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), resource["sha256"])
            self.assertEqual(path.stat().st_size, resource["size_bytes"])

    def test_recurrent_aml_patterns_remain_candidates_with_source_ids(self) -> None:
        payload = json.loads(
            (BUNDLE_ROOT / "hematology_rearrangements.v0.3.json").read_text(encoding="utf-8")
        )
        by_name = {record["display_name"]: record for record in payload["records"]}
        for expected in ("RUNX1::RUNX1T1", "CBFB::MYH11", "PML::RARA", "BCR::ABL1"):
            self.assertIn(expected, by_name)
            self.assertTrue(by_name[expected]["source_ids"])
            self.assertIn("caveat", by_name[expected])

        self.assertIn("PICALM::MLLT10", by_name)
        self.assertEqual(by_name["PICALM::MLLT10"]["relevance"], "aml_relevant_pattern")

    def test_public_records_are_hematology_filtered_and_pathology_attributed(self) -> None:
        payload = json.loads(
            (BUNDLE_ROOT / "hematology_rearrangements.v0.3.json").read_text(encoding="utf-8")
        )
        self.assertEqual(payload["scope"]["record_count"], 38)
        self.assertEqual(payload["scope"]["pathology_association_count"], 74)
        by_name = {record["display_name"]: record for record in payload["records"]}
        self.assertNotIn("EML4::ALK", by_name)
        self.assertEqual(
            {item["disease_id"] for item in by_name["PML::RARA"]["pathologies"]},
            {"DOID:0060318", "DOID:0081081"},
        )
        self.assertIn(
            "DOID:9119",
            {item["disease_id"] for item in by_name["PICALM::MLLT10"]["pathologies"]},
        )


if __name__ == "__main__":
    unittest.main()
