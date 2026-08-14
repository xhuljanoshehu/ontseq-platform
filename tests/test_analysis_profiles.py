import unittest
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
PROFILE_DIR = ROOT / "configs" / "analysis_profiles"


def _profile(name: str) -> dict[str, Any]:
    with (PROFILE_DIR / name).open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    assert isinstance(loaded, dict)
    return loaded


class AnalysisProfileTests(unittest.TestCase):
    def test_all_profiles_remain_explicitly_unvalidated_and_proposal_only(self) -> None:
        for path in PROFILE_DIR.glob("*.yaml"):
            with self.subTest(profile=path.name):
                profile = _profile(path.name)
                self.assertEqual(profile["status"], "proposed_unvalidated")
                self.assertEqual(profile["modules"]["iscn"], {"mode": "proposal_only"})

    def test_lcwgs_cnv_has_candidates_but_no_unvalidated_primary(self) -> None:
        profile = _profile("lcwgs.yaml")
        cnv = profile["modules"]["cnv"]

        self.assertEqual(cnv["selection_status"], "benchmark_required")
        self.assertNotIn("primary", cnv)
        self.assertEqual(
            {candidate["provider"] for candidate in cnv["candidates"]},
            {"ichorcna", "qdnaseq_ace", "spectre"},
        )

    def test_adaptive_sv_has_candidates_but_no_unvalidated_primary(self) -> None:
        profile = _profile("adaptive_sampling.yaml")
        sv = profile["modules"]["sv"]

        self.assertEqual(sv["selection_status"], "benchmark_required")
        self.assertNotIn("primary", sv)
        self.assertEqual(set(sv["candidates"]), {"sniffles2", "cutesv", "savana", "severus"})
