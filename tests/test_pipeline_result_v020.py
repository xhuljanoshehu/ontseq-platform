from __future__ import annotations

import unittest

from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import (
    LegacyResourceContext,
    PipelineResult,
)


class PipelineResultV020Tests(unittest.TestCase):
    def test_pipeline_result_010_remains_readable_as_legacy_unspecified(self) -> None:
        payload = build_demo_result().model_dump(mode="json")
        payload["schema_version"] = "0.1.0"
        payload.pop("reference_context", None)
        payload.pop("sidecars", None)
        restored = PipelineResult.model_validate(payload)
        self.assertIsInstance(restored.reference_context, LegacyResourceContext)
        self.assertEqual(restored.reference_context.status, "legacy_unspecified")


if __name__ == "__main__":
    unittest.main()
