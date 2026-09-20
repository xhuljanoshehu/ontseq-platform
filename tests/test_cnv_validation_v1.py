from __future__ import annotations

import hashlib
import unittest
from datetime import UTC, datetime

from ontseq_platform.cnv_validation_v1 import (
    CnvValidationV1Design,
    build_cnv_validation_v1_matrix,
    build_cnv_validation_v1_registration,
    default_cnv_validation_v1_design,
    dilution_policy_from_validation_v1,
    expand_cnv_validation_v1_cells,
    materialize_synthetic_validation_v1_cohort,
    summarize_cnv_validation_v1,
)
from pydantic import ValidationError

from ontseq_platform.cnv_validation_contracts import CnvCallerLock
from ontseq_platform.models import EventType


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _callers() -> list[CnvCallerLock]:
    return [
        CnvCallerLock(
            caller_id="qdnaseq_ace",
            caller_version="synthetic-qdnaseq-v1",
            adapter_policy_sha256=_sha("qdnaseq-policy"),
            execution_identity_sha256=_sha("qdnaseq-runtime"),
        ),
        CnvCallerLock(
            caller_id="spectre",
            caller_version="synthetic-spectre-v1",
            adapter_policy_sha256=_sha("spectre-policy"),
            execution_identity_sha256=_sha("spectre-runtime"),
        ),
    ]


class CnvValidationV1Tests(unittest.TestCase):
    def test_default_design_has_preregistered_primary_factor_levels(self) -> None:
        design = default_cnv_validation_v1_design()

        self.assertEqual(design.coverage_levels_x, [0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
        self.assertEqual(design.tumor_fractions, [1.0, 0.5, 0.3, 0.2, 0.1, 0.05])
        self.assertEqual(
            design.subchromosomal_sizes_bp,
            [100_000, 250_000, 500_000, 1_000_000, 5_000_000, 10_000_000],
        )
        self.assertEqual(design.replicates, 3)
        self.assertTrue(design.research_only)

    def test_invalid_factor_levels_fail_closed(self) -> None:
        payload = default_cnv_validation_v1_design().model_dump()
        payload["coverage_levels_x"] = [1.0, 0.5, 1.0]
        with self.assertRaises(ValidationError):
            CnvValidationV1Design.model_validate(payload)

        payload = default_cnv_validation_v1_design().model_dump()
        payload["tumor_fractions"] = [0.5, 1.0]
        with self.assertRaises(ValidationError):
            CnvValidationV1Design.model_validate(payload)

        payload = default_cnv_validation_v1_design().model_dump()
        payload["subchromosomal_sizes_bp"] = [1_000_000, 500_000]
        with self.assertRaises(ValidationError):
            CnvValidationV1Design.model_validate(payload)

    def test_cell_expansion_is_deterministic_and_whole_chromosome_has_no_numeric_size(self) -> None:
        design = default_cnv_validation_v1_design()
        cells = expand_cnv_validation_v1_cells(design)

        expected_subchromosomal = 6 * 6 * 2 * 6 * 3
        expected_whole_chromosome = 6 * 6 * 2 * 3
        self.assertEqual(len(cells), expected_subchromosomal + expected_whole_chromosome)
        self.assertEqual(len({cell.cell_id for cell in cells}), len(cells))

        whole = [
            cell
            for cell in cells
            if cell.event_type in {EventType.CHROMOSOME_GAIN, EventType.CHROMOSOME_LOSS}
        ]
        self.assertTrue(whole)
        self.assertTrue(all(cell.event_size_bp is None for cell in whole))

        first = expand_cnv_validation_v1_cells(design)
        self.assertEqual(first, cells)

    def test_subchromosomal_truth_event_exactly_matches_requested_size(self) -> None:
        design = default_cnv_validation_v1_design()
        cohort = materialize_synthetic_validation_v1_cohort(design)
        event_by_specimen = {
            specimen.specimen_id: specimen.truth_events[0] for specimen in cohort.specimens
        }

        for cell in expand_cnv_validation_v1_cells(design):
            if cell.event_size_bp is None:
                continue
            event = event_by_specimen[cell.specimen_id]
            self.assertEqual(event.primary.end - event.primary.start, cell.event_size_bp)

    def test_synthetic_cohort_preserves_continuous_factors_and_repeat_groups(self) -> None:
        design = default_cnv_validation_v1_design()
        cells = expand_cnv_validation_v1_cells(design)
        cohort = materialize_synthetic_validation_v1_cohort(design)
        specimen_by_id = {item.specimen_id: item for item in cohort.specimens}

        self.assertEqual(len(specimen_by_id), len(cells))
        self.assertEqual(len({item.input_sha256 for item in cohort.specimens}), len(cells))

        for cell in cells[:50]:
            specimen = specimen_by_id[cell.specimen_id]
            self.assertEqual(specimen.coverage_x, cell.coverage_x)
            self.assertEqual(specimen.tumor_fraction, cell.tumor_fraction)
            self.assertEqual(specimen.repeat_group_id, cell.repeat_group_id)
            self.assertEqual(len(specimen.truth_events), 1)

    def test_matrix_derives_unique_strata_cutpoints_from_planned_levels(self) -> None:
        design = default_cnv_validation_v1_design()
        matrix = build_cnv_validation_v1_matrix(design, _callers())

        self.assertEqual(matrix.stratification.coverage_cutpoints_x, [0.75, 1.5, 3.5, 7.5, 15.0])
        self.assertEqual(
            matrix.stratification.tumor_fraction_cutpoints,
            [0.075, 0.15, 0.25, 0.4, 0.75],
        )
        self.assertEqual(
            matrix.stratification.event_size_cutpoints_bp,
            [175_000, 375_000, 750_000, 3_000_000, 7_500_000],
        )
        self.assertEqual(matrix.qdnaseq_bin_sizes_kbp, [100, 500, 1000])

    def test_lcwgs_v1_rejects_incompatible_ichorcna_caller_lock(self) -> None:
        design = default_cnv_validation_v1_design()
        callers = _callers() + [
            CnvCallerLock(
                caller_id="ichorcna",
                caller_version="0.5.1",
                adapter_policy_sha256=_sha("ichor-policy"),
                execution_identity_sha256=_sha("ichor-runtime"),
            )
        ]

        with self.assertRaisesRegex(ValueError, "lcWGS|compatible"):
            build_cnv_validation_v1_matrix(design, callers)

    def test_registration_is_locked_and_changes_when_factor_or_caller_changes(self) -> None:
        design = default_cnv_validation_v1_design()
        first = build_cnv_validation_v1_registration(
            design,
            _callers(),
            registration_id="analytical-validation-v1-registration",
            registered_at=datetime(2026, 9, 20, tzinfo=UTC),
            code_sha256=_sha("code-v1"),
            software_version="0.8.2",
        )
        again = build_cnv_validation_v1_registration(
            design,
            _callers(),
            registration_id="analytical-validation-v1-registration",
            registered_at=datetime(2026, 9, 20, tzinfo=UTC),
            code_sha256=_sha("code-v1"),
            software_version="0.8.2",
        )
        self.assertEqual(first.lock_sha256, again.lock_sha256)
        self.assertTrue(first.outcomes_unseen)

        changed_payload = design.model_dump()
        changed_payload["coverage_levels_x"] = [0.5, 1.0, 2.0, 5.0, 10.0]
        changed = CnvValidationV1Design.model_validate(changed_payload)
        changed_registration = build_cnv_validation_v1_registration(
            changed,
            _callers(),
            registration_id="analytical-validation-v1-registration",
            registered_at=datetime(2026, 9, 20, tzinfo=UTC),
            code_sha256=_sha("code-v1"),
            software_version="0.8.2",
        )
        self.assertNotEqual(first.lock_sha256, changed_registration.lock_sha256)

        caller_changed = _callers()
        caller_changed[1] = caller_changed[1].model_copy(
            update={"execution_identity_sha256": _sha("spectre-runtime-v2")}
        )
        caller_registration = build_cnv_validation_v1_registration(
            design,
            caller_changed,
            registration_id="analytical-validation-v1-registration",
            registered_at=datetime(2026, 9, 20, tzinfo=UTC),
            code_sha256=_sha("code-v1"),
            software_version="0.8.2",
        )
        self.assertNotEqual(first.lock_sha256, caller_registration.lock_sha256)

    def test_dilution_projection_reuses_existing_fraction_order_and_replicates(self) -> None:
        design = default_cnv_validation_v1_design()
        policy = dilution_policy_from_validation_v1(design, expected_samtools_version="1.24")

        self.assertEqual(policy.tumor_fractions, design.tumor_fractions)
        self.assertEqual(policy.replicates, design.replicates)
        self.assertTrue(policy.include_normal_only_control)
        self.assertEqual(policy.status, "technical_defaults_only")

    def test_summary_reports_truth_cells_not_false_caller_votes(self) -> None:
        summary = summarize_cnv_validation_v1(default_cnv_validation_v1_design(), _callers())

        self.assertEqual(summary["positive_truth_cells"], 864)
        self.assertEqual(summary["caller_count"], 2)
        self.assertEqual(summary["qdnaseq_bin_sizes_kbp"], [100, 500, 1000])
        self.assertNotIn("winner", summary)
        self.assertNotIn("majority", summary)


if __name__ == "__main__":
    unittest.main()
