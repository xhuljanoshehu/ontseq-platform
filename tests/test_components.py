"""Component selection: what a run may ask for, and what it may not get away with."""

from __future__ import annotations

import unittest

from pydantic import ValidationError

from ontseq_platform.pipeline.components import (
    SUPPORTED_PROVIDERS,
    ComponentChoice,
    ComponentVersionMismatch,
    RunComponents,
)
from ontseq_platform.pipeline.stages import StageId


def _selection(**components: object) -> RunComponents:
    return RunComponents.model_validate(
        {
            "selection_id": "test-selection",
            "status": "technical_defaults_only",
            "components": components,
        }
    )


class SelectionValidationTests(unittest.TestCase):
    def test_a_provider_without_an_adapter_is_refused_when_the_file_is_read(self) -> None:
        """Refusing late would cost a run that had already spent hours on alignment."""
        with self.assertRaises(ValidationError) as caught:
            _selection(sv={"provider": "cutesv", "version": "2.1.0"})
        self.assertIn("no adapter in this repository", str(caught.exception))
        self.assertIn("sniffles2", str(caught.exception))

    def test_a_stage_that_runs_no_tool_cannot_be_selected(self) -> None:
        with self.assertRaises(ValidationError) as caught:
            _selection(report={"provider": "anything"})
        self.assertIn("runs no external tool", str(caught.exception))

    def test_release_and_assemble_are_not_selectable(self) -> None:
        """A run must not be able to switch off its own report or checksum bundle."""
        for stage in (StageId.ASSEMBLE, StageId.REPORT, StageId.RELEASE):
            self.assertNotIn(stage, SUPPORTED_PROVIDERS)

    def test_a_version_must_look_like_a_version(self) -> None:
        with self.assertRaises(ValidationError):
            _selection(sv={"provider": "sniffles2", "version": "latest"})


class VersionEnforcementTests(unittest.TestCase):
    def test_the_selected_version_must_be_the_installed_one(self) -> None:
        choice = ComponentChoice(provider="sniffles2", version="2.4.0")
        with self.assertRaises(ComponentVersionMismatch) as caught:
            choice.verify({"sniffles": "2.8.0"}, stage=StageId.SV)
        message = str(caught.exception)
        self.assertIn("2.4.0", message)
        self.assertIn("2.8.0", message)

    def test_a_matching_version_passes(self) -> None:
        ComponentChoice(provider="sniffles2", version="2.8.0").verify(
            {"sniffles": "2.8.0"}, stage=StageId.SV
        )

    def test_a_missing_probe_is_a_mismatch_not_a_pass(self) -> None:
        """Silently accepting an absent version would defeat the point of pinning."""
        choice = ComponentChoice(provider="mosdepth", version="0.3.14")
        with self.assertRaises(ComponentVersionMismatch):
            choice.verify({}, stage=StageId.TARGET_COVERAGE)

    def test_an_unpinned_component_accepts_whatever_is_installed(self) -> None:
        ComponentChoice(provider="cramino").verify({"cramino": "1.3.0"}, stage=StageId.QC)

    def test_the_provider_name_and_the_probe_key_may_differ(self) -> None:
        ComponentChoice(provider="sniffles2", version="2.8.0").verify(
            {"sniffles": "2.8.0"}, stage=StageId.SV
        )


class DeselectionTests(unittest.TestCase):
    def test_without_disables_an_existing_entry_and_keeps_its_provider(self) -> None:
        selection = _selection(cnv={"provider": "qdnaseq_ace"}).without(StageId.CNV)
        choice = selection.choice_for(StageId.CNV)
        assert choice is not None
        self.assertFalse(choice.enabled)
        self.assertEqual(choice.provider, "qdnaseq_ace")

    def test_without_can_disable_a_stage_the_file_never_mentioned(self) -> None:
        selection = _selection().without(StageId.SV)
        self.assertEqual(selection.disabled_stages(), (StageId.SV,))

    def test_a_non_selectable_stage_cannot_be_disabled(self) -> None:
        with self.assertRaises(ValueError):
            _selection().without(StageId.RELEASE)

    def test_unpinned_stages_are_reported_so_the_run_can_warn(self) -> None:
        selection = _selection(
            qc={"provider": "cramino"},
            sv={"provider": "sniffles2", "version": "2.8.0"},
        )
        self.assertEqual(selection.unpinned_stages(), (StageId.QC,))

    def test_the_summary_names_state_provider_and_version(self) -> None:
        lines = _selection(sv={"provider": "sniffles2", "version": "2.4.0"}).summary()
        self.assertEqual(len(lines), 1)
        self.assertIn("sniffles2", lines[0])
        self.assertIn("2.4.0", lines[0])
        self.assertTrue(lines[0].startswith("on "))


if __name__ == "__main__":
    unittest.main()
