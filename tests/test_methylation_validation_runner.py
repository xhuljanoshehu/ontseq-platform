from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from datetime import UTC, datetime, timedelta
from pathlib import Path
from unittest.mock import patch

from ontseq_platform import __version__, entrypoint, methylation_validation_runner
from ontseq_platform.methylation_mixture import MethylationMixturePolicy, NanopolishSourceMetadata
from ontseq_platform.methylation_validation import (
    SAMPLE_ROLES,
    MethylationValidationEvidence,
    MethylationValidationMatrix,
    MethylationValidationRegistration,
    MethylationValidationReport,
    RecoveryAcceptance,
    ValidationCohort,
    ValidationDecision,
    ValidationSample,
    ValidationTechnicalMetadata,
    content_sha256,
    evaluate_methylation_validation,
    preregister_validation,
)
from ontseq_platform.methylation_validation_cli import render_validation_report
from ontseq_platform.methylation_validation_runner import (
    ValidationLocalInput,
    ValidationLocalInputs,
    execute_registered_validation,
    validation_software_sha256,
)
from ontseq_platform.modbam import (
    MODBAM_ADAPTER_VERSION,
    ModbamAdapterPolicy,
    ModbamSourceMetadata,
)
from ontseq_platform.models import GenomeBuild
from ontseq_platform.reference import sha256_file


def _write_synthetic_table(path: Path, role: str, *, reads: int = 30, ratio: float) -> None:
    lines = [
        "chromosome\tstrand\tstart\tend\tread_name\tlog_lik_ratio\t"
        "log_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_motifs\tsequence"
    ]
    for read_index in range(reads):
        for marker in range(6):
            position = 1000 + 10 * marker
            lines.append(
                f"chr1\t+\t{position}\t{position}\tSYN_{role}_{read_index:03d}\t"
                f"{ratio}\t{ratio}\t0\t1\t1\tACGTA"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class MethylationValidationRunnerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.policy = MethylationMixturePolicy(
            profile_id="synthetic-validation-runner-v1",
            source_a_fractions=[0.0, 0.5, 1.0],
            replicates=10,
            minimum_markers=6,
            maximum_markers=6,
            minimum_reference_valid_calls=5,
            minimum_mixture_valid_calls=3,
            minimum_absolute_delta_beta=0.5,
            uncertainty_draws=100,
            maximum_confidence_interval_width=1.0,
            note="Synthetic execution regression fixture; no biological validation claim.",
        )
        self.matrix = MethylationValidationMatrix(
            matrix_id="synthetic-runner-matrix-v1",
            source_a_fractions=[0.0, 0.5, 1.0],
            read_group_budgets=[20],
            replicate_seeds=list(range(11, 21)),
            acceptance=[
                RecoveryAcceptance(
                    read_group_budget=20,
                    maximum_absolute_bias=0.01,
                    maximum_mae=0.01,
                    maximum_rmse=0.01,
                    minimum_conditional_interval_coverage=0.9,
                    maximum_no_call_rate=0.0,
                    maximum_per_fraction_absolute_bias=0.01,
                    minimum_per_fraction_conditional_interval_coverage=0.9,
                    maximum_per_fraction_no_call_rate=0.0,
                )
            ],
            estimator_policy=self.policy,
            rationale="Synthetic regression controls for an idealised two-source experiment.",
        )
        self.metadata = NanopolishSourceMetadata(
            caller_name="nanopolish",
            caller_version="0.14.0",
            reference_id="SYNTHETIC_REFERENCE",
            reference_genome_build=GenomeBuild.GRCH38,
            reference_sha256="e" * 64,
            sequencing_platform="SYNTHETIC_ONT",
            flow_cell_product_code="SYNTHETIC_FLOWCELL",
            library_kit="SYNTHETIC_KIT",
            basecaller_name="SYNTHETIC_BASECALLER",
            basecaller_version="1.0",
            basecaller_model="SYNTHETIC_MODEL",
            modification_context="CpG_5mC_vs_unmethylated",
        )
        technical = ValidationTechnicalMetadata(
            **self.metadata.model_dump(exclude={"schema_version"}),
            adapter_name="nanopolish_call_methylation_tsv",
            adapter_version="nanopolish-call-table-v2",
            adapter_policy_sha256=content_sha256(self.policy),
        )
        bindings = {}
        samples = {}
        for role in SAMPLE_ROLES:
            path = self.root / f"{role}.tsv"
            _write_synthetic_table(path, role, ratio=5.0 if role.endswith("_a") else -5.0)
            metadata_path = self.root / f"{role}.metadata.json"
            metadata_path.write_text(self.metadata.model_dump_json(), encoding="utf-8")
            bindings[role] = ValidationLocalInput(
                calls_path=path,
                input_format="nanopolish",
                metadata_path=metadata_path,
            )
            source = "SYNTHETIC_SOURCE_A" if role.endswith("_a") else "SYNTHETIC_SOURCE_B"
            samples[role] = ValidationSample(
                accession=f"SYNTHETIC_ACCESSION_{role}",
                dataset_version="synthetic-v1",
                biological_sample_id=f"SYNTHETIC_SAMPLE_{role}",
                donor_id=f"SYNTHETIC_DONOR_{role}",
                biological_source=source,
                cell_type=source,
                material_type="SYNTHETIC_MATERIAL",
                input_sha256=sha256_file(path),
                provenance_reference="Synthetic regression fixture; generated within this test.",
                technical=technical,
            )
        self.inputs = ValidationLocalInputs(samples=bindings)
        self.cohort = ValidationCohort(
            samples=samples,
            independence_evidence="Synthetic distinct labels exercise declarations only.",
            access_basis="public",
        )
        self.registration = self._register()

    def _register(self) -> MethylationValidationRegistration:
        return preregister_validation(
            self.matrix,
            self.cohort,
            registration_id="SYNTHETIC_RUNNER_VALIDATION",
            registered_at=datetime.now(UTC) - timedelta(minutes=1),
            test_outcomes_unseen=True,
            code_sha256=validation_software_sha256(),
            software_version=__version__,
        )

    def _refresh_input_lock(self, role: str) -> None:
        self.cohort.samples[role].input_sha256 = sha256_file(self.inputs.samples[role].calls_path)
        self.registration = self._register()

    def _cli(self, *arguments: str) -> str:
        output = io.StringIO()
        with patch("sys.argv", ["ontseq", *arguments]), contextlib.redirect_stdout(output):
            entrypoint.main()
        return output.getvalue()

    def test_ideal_synthetic_execution_is_deterministic_and_bound(self) -> None:
        first = execute_registered_validation(self.registration, self.inputs)
        second = execute_registered_validation(self.registration, self.inputs)
        self.assertEqual(first.decision, ValidationDecision.PASS)
        self.assertEqual(first.overall_metrics.completed_levels, 30)
        self.assertAlmostEqual(first.overall_metrics.mean_absolute_error, 0.0)
        self.assertEqual(first.overall_metrics.no_call_rate, 0.0)
        self.assertIsNotNone(first.evidence)
        self.assertIsNotNone(second.evidence)
        self.assertEqual(first.evidence.levels, second.evidence.levels)
        self.assertEqual(first.evidence.reference_markers, second.evidence.reference_markers)
        self.assertEqual(first.evidence.registration_sha256, self.registration.lock_sha256)
        self.assertEqual(
            first.evidence.observed_input_sha256,
            {role: sample.input_sha256 for role, sample in self.cohort.samples.items()},
        )
        self.assertTrue(first.research_only)
        self.assertIn("retained whole read groups", " ".join(first.limitations))
        self.assertEqual(
            MethylationValidationReport.model_validate_json(first.model_dump_json()), first
        )

    def test_control_suffix_cannot_disguise_reused_donor_in_registered_runner(self) -> None:
        # Previously this full fixture could register and execute with PASS because
        # the NUL suffix distinguished two visually identical donor declarations.
        for role, sample in self.cohort.samples.items():
            sample.donor_id = "SYNTHETIC_SHARED_DONOR" + ("\x00" if role.startswith("test") else "")
        with self.assertRaisesRegex(ValueError, "control characters"):
            registration = self._register()
            execute_registered_validation(registration, self.inputs)

    def test_dorado_sam_adapter_executes_the_four_sample_matrix(self) -> None:
        # Six synthetic CpGs on each whole read; no external biological input.
        reference = self.root / "synthetic.fa"
        reference.write_bytes(b">chr1\nCGCGCGCGCGCG\n")
        adapter_policy = ModbamAdapterPolicy()
        policy_path = self.root / "adapter-policy.json"
        policy_path.write_text(adapter_policy.model_dump_json(), encoding="utf-8")
        metadata = ModbamSourceMetadata(
            caller_version="synthetic-1",
            reference_id="synthetic-reference-v1",
            reference_genome_build=GenomeBuild.GRCH38,
            reference_sha256=sha256_file(reference),
            sequencing_platform="ONT",
            flow_cell_product_code="synthetic",
            library_kit="synthetic",
            basecaller_version="synthetic-1",
            basecaller_model="synthetic-cpg-v1",
            modification_model="synthetic-5mc-v1",
            cytosine_modifications="5mC",
        )
        technical = ValidationTechnicalMetadata(
            **metadata.model_dump(exclude={"schema_version", "source_modbam_sha256"}),
            adapter_name="sam_mm_ml",
            adapter_version=MODBAM_ADAPTER_VERSION,
            adapter_policy_sha256=content_sha256(adapter_policy),
        )
        for role in SAMPLE_ROLES:
            path = self.root / f"{role}.sam"
            probability = "255" if role.endswith("_a") else "0"
            ml = ",".join([probability] * 6)
            lines = ["@SQ\tSN:chr1\tLN:12", "@PG\tID:dorado\tPN:dorado\tVN:synthetic-1"]
            lines.extend(
                f"SYN_{role}_{i}\t0\tchr1\t1\t60\t12M\t*\t0\t0\tCGCGCGCGCGCG\t*"
                f"\tMM:Z:C+m?,0,0,0,0,0,0;\tML:B:C,{ml}"
                for i in range(30)
            )
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            metadata_path = self.root / f"{role}.dorado.json"
            metadata_path.write_text(metadata.model_dump_json(), encoding="utf-8")
            self.inputs.samples[role] = ValidationLocalInput(
                calls_path=path,
                input_format="modbam",
                metadata_path=metadata_path,
                reference_path=reference,
                adapter_policy_path=policy_path,
            )
            self.cohort.samples[role].input_sha256 = sha256_file(path)
            self.cohort.samples[role].technical = technical
        self.registration = self._register()
        report = execute_registered_validation(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.PASS)
        self.assertEqual(report.overall_metrics.completed_levels, 30)
        self.assertAlmostEqual(report.overall_metrics.mean_absolute_error, 0.0)
        MethylationValidationReport.model_validate_json(report.model_dump_json())

    def test_input_modified_after_registration_is_rejected(self) -> None:
        with self.inputs.samples["test_a"].calls_path.open("a", encoding="utf-8") as handle:
            handle.write("\n")
        with self.assertRaisesRegex(ValueError, "registered SHA-256"):
            execute_registered_validation(self.registration, self.inputs)

    def test_local_metadata_mismatch_is_rejected(self) -> None:
        metadata = self.metadata.model_dump(mode="json")
        metadata["basecaller_version"] = "2.0"
        self.inputs.samples["test_a"].metadata_path.write_text(
            json.dumps(metadata), encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "metadata differs at basecaller_version"):
            execute_registered_validation(self.registration, self.inputs)

    def test_metadata_cannot_silently_drop_reference_provenance(self) -> None:
        metadata = self.metadata.model_dump(mode="json")
        del metadata["reference_sha256"]
        self.inputs.samples["test_b"].metadata_path.write_text(
            json.dumps(metadata), encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "complete Nanopolish metadata"):
            execute_registered_validation(self.registration, self.inputs)

    def test_same_donor_calibration_and_test_is_no_call_before_input_access(self) -> None:
        self.cohort.samples["test_a"].donor_id = self.cohort.samples["calibration_a"].donor_id
        self.registration = self._register()
        self.inputs.samples["test_a"].calls_path = self.root / "does-not-exist.tsv"
        report = execute_registered_validation(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertIsNone(report.evidence)
        self.assertEqual(report.overall_metrics.no_call_rate, 1.0)
        self.assertTrue(any("donors overlap" in reason for reason in report.reasons))

    def test_reused_read_ids_are_rejected_even_with_distinct_file_checksums(self) -> None:
        path = self.inputs.samples["test_a"].calls_path
        _write_synthetic_table(path, "calibration_a", ratio=6.0)
        self._refresh_input_lock("test_a")
        with self.assertRaisesRegex(ValueError, "Read identifiers overlap"):
            execute_registered_validation(self.registration, self.inputs)

    def test_short_test_pool_cannot_reduce_registered_budget(self) -> None:
        _write_synthetic_table(
            self.inputs.samples["test_b"].calls_path, "test_b", reads=12, ratio=-5
        )
        self._refresh_input_lock("test_b")
        report = execute_registered_validation(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.no_call_rate, 1.0)
        self.assertIn("largest registered constant", report.execution_no_call_reason or "")
        MethylationValidationReport.model_validate_json(report.model_dump_json())

    def test_software_change_after_registration_is_rejected(self) -> None:
        with (
            patch.object(
                methylation_validation_runner, "validation_software_sha256", return_value="f" * 64
            ),
            self.assertRaisesRegex(ValueError, "Software/runtime changed"),
        ):
            execute_registered_validation(self.registration, self.inputs)

    def test_source_change_during_execution_is_rejected(self) -> None:
        original = methylation_validation_runner._run_independent_sources

        def mutate_after_execution(*args: object, **kwargs: object) -> MethylationValidationReport:
            report = original(*args, **kwargs)
            with self.inputs.samples["test_b"].calls_path.open("a", encoding="utf-8") as handle:
                handle.write("\n")
            return report

        with (
            patch.object(
                methylation_validation_runner,
                "_run_independent_sources",
                side_effect=mutate_after_execution,
            ),
            self.assertRaisesRegex(ValueError, "input changed during"),
        ):
            execute_registered_validation(self.registration, self.inputs)

    def test_uninformative_calibration_produces_explicit_complete_grid_no_call(self) -> None:
        _write_synthetic_table(
            self.inputs.samples["calibration_b"].calls_path, "calibration_b", ratio=5
        )
        self._refresh_input_lock("calibration_b")
        report = execute_registered_validation(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.observed_levels, 30)
        self.assertEqual(report.overall_metrics.explicit_no_call_levels, 30)
        self.assertEqual(report.overall_metrics.missing_levels, 0)
        self.assertIsNone(report.overall_metrics.mean_absolute_error)

    def test_cli_plan_is_no_call_and_cannot_overwrite_existing_lock(self) -> None:
        matrix_path = self.root / "matrix.json"
        matrix_path.write_text(self.matrix.model_dump_json(), encoding="utf-8")
        output = self.root / "readiness.json"
        args = (
            "methylation-validation-plan",
            "--matrix",
            str(matrix_path),
            "--output",
            str(output),
        )
        self.assertIn("NO_CALL", self._cli(*args))
        before = output.read_bytes()
        readiness = json.loads(before)
        self.assertEqual(readiness["decision"], "NO_CALL")
        self.assertEqual(readiness["expected_levels"], 30)
        self.assertEqual(readiness["matrix_sha256"], content_sha256(self.matrix))
        with self.assertRaisesRegex(SystemExit, "ERROR"):
            self._cli(*args)
        self.assertEqual(output.read_bytes(), before)

    def test_cli_evaluate_missing_evidence_preserves_unavailable_metrics(self) -> None:
        registration_path = self.root / "registration.json"
        registration_path.write_text(self.registration.model_dump_json(), encoding="utf-8")
        output_dir = self.root / "empty-evaluation"
        self.assertIn(
            "NO_CALL",
            self._cli(
                "methylation-validation-evaluate",
                "--registration",
                str(registration_path),
                "--output-dir",
                str(output_dir),
            ),
        )
        json_path = next(output_dir.glob("*.json"))
        report = MethylationValidationReport.model_validate_json(json_path.read_text())
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.missing_levels, 30)
        self.assertIsNone(report.overall_metrics.mean_absolute_error)
        self.assertIn("Research use only", next(output_dir.glob("*.html")).read_text())
        self.assertEqual(len(next(output_dir.glob("*.csv")).read_text().splitlines()), 1)

    def test_cli_register_run_evaluate_and_no_overwrite(self) -> None:
        matrix_path, cohort_path, inputs_path = [
            self.root / name for name in ("matrix.json", "cohort.json", "inputs.json")
        ]
        for path, value in (
            (matrix_path, self.matrix),
            (cohort_path, self.cohort),
            (inputs_path, self.inputs),
        ):
            path.write_text(value.model_dump_json(), encoding="utf-8")
        registration_path = self.root / "cli-registration.json"
        self._cli(
            "methylation-validation-register",
            "--matrix",
            str(matrix_path),
            "--cohort",
            str(cohort_path),
            "--registration-id",
            "SYNTHETIC_CLI",
            "--confirm-test-outcomes-unseen",
            "--output",
            str(registration_path),
        )
        run_output = self.root / "run"
        args = (
            "methylation-validation-run",
            "--registration",
            str(registration_path),
            "--inputs",
            str(inputs_path),
            "--output-dir",
            str(run_output),
        )
        self.assertIn("PASS", self._cli(*args))
        files_before = {path.name: path.read_bytes() for path in run_output.iterdir()}
        self.assertEqual(len(files_before), 3)
        with self.assertRaisesRegex(SystemExit, "already exists"):
            self._cli(*args)
        self.assertEqual(
            files_before, {path.name: path.read_bytes() for path in run_output.iterdir()}
        )
        report = MethylationValidationReport.model_validate_json(
            next(run_output.glob("*.json")).read_text()
        )
        evidence_path = self.root / "evidence.json"
        evidence_path.write_text(report.evidence.model_dump_json(), encoding="utf-8")
        evaluation_output = self.root / "evaluation"
        self.assertIn(
            "PASS",
            self._cli(
                "methylation-validation-evaluate",
                "--registration",
                str(registration_path),
                "--evidence",
                str(evidence_path),
                "--output-dir",
                str(evaluation_output),
            ),
        )
        evaluated = MethylationValidationReport.model_validate_json(
            next(evaluation_output.glob("*.json")).read_text()
        )
        self.assertEqual(evaluated, report)

    def test_cli_tampered_registration_checksum_cannot_run(self) -> None:
        payload = self.registration.model_dump(mode="json")
        payload["matrix"]["replicate_seeds"][0] = 109
        path = self.root / "tampered-registration.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "Matrix content does not match"):
            self._cli(
                "methylation-validation-evaluate",
                "--registration",
                str(path),
                "--output-dir",
                str(self.root / "tampered"),
            )
        self.assertFalse((self.root / "tampered").exists())

    def test_existing_partial_output_is_not_overwritten(self) -> None:
        report = evaluate_methylation_validation(self.registration)
        output = self.root / "partial"
        output.mkdir()
        sentinel = output / f"{self.registration.registration_id}.methylation-validation.html"
        sentinel.write_text("preserve existing output", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "already exists"):
            render_validation_report(report, output)
        self.assertEqual(sentinel.read_text(), "preserve existing output")
        self.assertEqual(list(output.iterdir()), [sentinel])

    def test_imported_evidence_cannot_change_point_estimate_without_replay(self) -> None:
        report = execute_registered_validation(self.registration, self.inputs)
        evidence_data = report.evidence.model_dump(mode="json")
        evidence_data["levels"][10]["result"]["estimated_source_a_fraction"] = 0.55
        with self.assertRaises(ValueError):
            evidence = MethylationValidationEvidence.model_validate(evidence_data)
            evaluate_methylation_validation(self.registration, evidence)


if __name__ == "__main__":
    unittest.main()
