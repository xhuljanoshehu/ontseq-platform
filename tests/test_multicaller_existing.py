from __future__ import annotations

import importlib
import unittest

from ontseq_platform.multicaller_existing import (
    EXISTING_CALLER_BRIDGES,
    get_existing_caller_bridge,
)

from ontseq_platform.multicaller_contracts import CallerMode


class ExistingMultiCallerBridgeTests(unittest.TestCase):
    def test_existing_runtime_bridges_resolve_the_real_adapter_symbols(self) -> None:
        expected = {
            CallerMode.QDNASEQ_ACE_MULTIBIN: (
                "ontseq_platform.cnv.qdnaseq",
                "run_qdnaseq_ace",
            ),
            CallerMode.SNIFFLES2_STANDARD: (
                "ontseq_platform.sniffles",
                "run_sniffles",
            ),
            CallerMode.CUTESV_STANDARD: (
                "ontseq_platform.cutesv",
                "run_cutesv",
            ),
        }

        for mode, (module_name, symbol_name) in expected.items():
            with self.subTest(mode=mode.value):
                bridge = get_existing_caller_bridge(mode)
                self.assertTrue(bridge.runtime_qualified)
                self.assertEqual(bridge.execution_module, module_name)
                self.assertEqual(bridge.execution_symbol, symbol_name)
                module = importlib.import_module(bridge.execution_module)
                self.assertTrue(callable(getattr(module, bridge.execution_symbol)))

    def test_qdnaseq_bridge_reuses_block4_evidence_adapter(self) -> None:
        bridge = get_existing_caller_bridge(CallerMode.QDNASEQ_ACE_MULTIBIN)

        self.assertEqual(
            bridge.evidence_module,
            "ontseq_platform.cnv_validation_qdnaseq",
        )
        self.assertEqual(
            bridge.evidence_symbol,
            "build_qdnaseq_ace_validation_manifest",
        )
        assert bridge.evidence_module is not None
        assert bridge.evidence_symbol is not None
        module = importlib.import_module(bridge.evidence_module)
        self.assertTrue(callable(getattr(module, bridge.evidence_symbol)))

    def test_bridge_contract_never_treats_caller_agreement_as_truth(self) -> None:
        for mode, bridge in EXISTING_CALLER_BRIDGES.items():
            with self.subTest(mode=mode.value):
                self.assertFalse(bridge.caller_agreement_is_truth)
                self.assertTrue(bridge.retain_native_output)

    def test_sniffles_mosaic_is_catalogued_but_not_runtime_qualified(self) -> None:
        bridge = get_existing_caller_bridge(CallerMode.SNIFFLES2_MOSAIC)

        self.assertFalse(bridge.runtime_qualified)
        self.assertIsNone(bridge.execution_symbol)
        self.assertIn("separate", bridge.qualification_reason.lower())

    def test_new_unimplemented_callers_have_no_existing_bridge(self) -> None:
        for mode in (
            CallerMode.SPECTRE_DEPTH_ONLY,
            CallerMode.SEVERUS_PAIRED,
            CallerMode.SAVANA_PAIRED,
            CallerMode.WAKHAN_PHASED_CNA,
            CallerMode.ICHORCNA_ULP_WGS,
        ):
            with self.subTest(mode=mode.value), self.assertRaises(ValueError):
                get_existing_caller_bridge(mode)


if __name__ == "__main__":
    unittest.main()
