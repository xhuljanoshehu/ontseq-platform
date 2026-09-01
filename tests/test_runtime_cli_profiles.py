from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from ontseq_platform.models import CuteSvPolicy
from ontseq_platform.runtime_cli import _cutesv_policy_for_run


class OptionalCuteSvPolicyTests(unittest.TestCase):
    def test_legacy_run_without_reference_does_not_implicitly_enable_cutesv(self) -> None:
        with patch("ontseq_platform.runtime_cli._cutesv_policy") as load_policy:
            policy = _cutesv_policy_for_run(Path("default-cutesv.yaml"), None)

        self.assertIsNone(policy)
        load_policy.assert_not_called()

    def test_run_with_reference_keeps_the_dual_caller_policy(self) -> None:
        expected = CuteSvPolicy(
            profile_id="synthetic-cutesv",
            status="technical_defaults_only",
            note="test",
        )
        with patch(
            "ontseq_platform.runtime_cli._cutesv_policy",
            return_value=expected,
        ) as load_policy:
            policy = _cutesv_policy_for_run(
                Path("default-cutesv.yaml"),
                Path("reference.fa"),
            )

        self.assertIs(policy, expected)
        load_policy.assert_called_once_with(Path("default-cutesv.yaml"))


if __name__ == "__main__":
    unittest.main()
