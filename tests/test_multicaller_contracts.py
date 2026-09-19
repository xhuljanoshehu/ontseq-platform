from __future__ import annotations

import unittest

from ontseq_platform.multicaller_contracts import (
    CALLER_CATALOG,
    CallerAnalyticalDomain,
    CallerAssayRegime,
    CallerInputRole,
    CallerMode,
    CallerProvider,
    get_caller_catalog_entry,
)


class MultiCallerCatalogTests(unittest.TestCase):
    def test_catalog_contains_exact_selected_provider_modes(self) -> None:
        self.assertEqual(
            set(CALLER_CATALOG),
            {
                CallerMode.QDNASEQ_ACE_MULTIBIN,
                CallerMode.SPECTRE_DEPTH_ONLY,
                CallerMode.SPECTRE_SNIFFLES_SUPPORTED,
                CallerMode.SNIFFLES2_STANDARD,
                CallerMode.SNIFFLES2_MOSAIC,
                CallerMode.CUTESV_STANDARD,
                CallerMode.SEVERUS_PAIRED,
                CallerMode.SEVERUS_SINGLE_SAMPLE,
                CallerMode.SAVANA_PAIRED,
                CallerMode.SAVANA_TUMOR_ONLY,
                CallerMode.WAKHAN_PHASED_CNA,
                CallerMode.ICHORCNA_ULP_WGS,
            },
        )

    def test_savana_paired_and_tumor_only_are_distinct_contracts(self) -> None:
        paired = get_caller_catalog_entry(CallerMode.SAVANA_PAIRED)
        tumor_only = get_caller_catalog_entry(CallerMode.SAVANA_TUMOR_ONLY)

        self.assertEqual(paired.provider, CallerProvider.SAVANA)
        self.assertEqual(tumor_only.provider, CallerProvider.SAVANA)
        self.assertNotEqual(paired.mode_id, tumor_only.mode_id)
        self.assertIn(CallerInputRole.MATCHED_NORMAL_BAM, paired.required_input_roles)
        self.assertNotIn(CallerInputRole.MATCHED_NORMAL_BAM, tumor_only.required_input_roles)
        self.assertIn(CallerInputRole.TUMOR_BAM, paired.required_input_roles)
        self.assertIn(CallerInputRole.TUMOR_BAM, tumor_only.required_input_roles)

    def test_spectre_depth_only_has_no_caller_parent_dependency(self) -> None:
        depth_only = get_caller_catalog_entry(CallerMode.SPECTRE_DEPTH_ONLY)
        supported = get_caller_catalog_entry(CallerMode.SPECTRE_SNIFFLES_SUPPORTED)

        self.assertEqual(depth_only.provider, CallerProvider.SPECTRE)
        self.assertEqual(depth_only.required_parent_modes, [])
        self.assertEqual(
            supported.required_parent_modes,
            [CallerMode.SNIFFLES2_STANDARD],
        )

    def test_every_catalog_entry_declares_domain_regime_and_inputs(self) -> None:
        for mode, entry in CALLER_CATALOG.items():
            with self.subTest(mode=mode.value):
                self.assertEqual(entry.mode_id, mode)
                self.assertGreaterEqual(len(entry.analytical_domains), 1)
                self.assertGreaterEqual(len(entry.compatible_regimes), 1)
                self.assertGreaterEqual(len(entry.required_input_roles), 1)
                self.assertTrue(
                    all(
                        isinstance(item, CallerAnalyticalDomain)
                        for item in entry.analytical_domains
                    )
                )
                self.assertTrue(
                    all(isinstance(item, CallerAssayRegime) for item in entry.compatible_regimes)
                )
                self.assertTrue(
                    all(isinstance(item, CallerInputRole) for item in entry.required_input_roles)
                )

    def test_catalog_entries_have_unique_mode_ids_and_provider_prefixes(self) -> None:
        mode_ids = [entry.mode_id for entry in CALLER_CATALOG.values()]
        self.assertEqual(len(mode_ids), len(set(mode_ids)))
        for entry in CALLER_CATALOG.values():
            self.assertTrue(entry.mode_id.value.startswith(f"{entry.provider.value}:"))

    def test_qdnaseq_and_ichorcna_live_in_different_assay_regimes(self) -> None:
        qdnaseq = get_caller_catalog_entry(CallerMode.QDNASEQ_ACE_MULTIBIN)
        ichorcna = get_caller_catalog_entry(CallerMode.ICHORCNA_ULP_WGS)

        self.assertIn(CallerAssayRegime.LCWGS, qdnaseq.compatible_regimes)
        self.assertNotIn(CallerAssayRegime.CFDNA_ULP_WGS, qdnaseq.compatible_regimes)
        self.assertEqual(ichorcna.compatible_regimes, [CallerAssayRegime.CFDNA_ULP_WGS])


if __name__ == "__main__":
    unittest.main()
