from __future__ import annotations

import csv
import hashlib
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from datetime import UTC, datetime, timedelta
from pathlib import Path
from unittest.mock import patch

from ontseq_platform import __version__, entrypoint, methylation_holdout_runner
from ontseq_platform.methylation_holdout import (
    HOLDOUT_SOURCE_ROLES,
    PairedHoldoutCohort,
    PairedHoldoutRecoveryReport,
    PairedHoldoutRegistration,
    PairedHoldoutSplitPolicy,
    holdout_split_manifest_sha256,
    preregister_paired_holdout,
)
from ontseq_platform.methylation_holdout_runner import (
    PairedHoldoutLocalInputs,
    execute_registered_paired_holdout,
)
from ontseq_platform.methylation_mixture import (
    MethylationMixturePolicy,
    NanopolishSourceMetadata,
    _selection_digest,
    _split_reads,
)
from ontseq_platform.methylation_validation import (
    MethylationValidationMatrix,
    RecoveryAcceptance,
    ValidationDecision,
    ValidationSample,
    ValidationTechnicalMetadata,
    content_sha256,
)
from ontseq_platform.methylation_validation_runner import (
    ValidationLocalInput,
    validation_software_sha256,
)
from ontseq_platform.modbam import (
    MODBAM_ADAPTER_VERSION,
    ModbamAdapterPolicy,
    ModbamSourceMetadata,
)
from ontseq_platform.models import GenomeBuild
from ontseq_platform.reference import sha256_file


def _names(role: str, reads: int = 40) -> tuple[str, ...]:
    return tuple(f"SYN_READ_{role}_{index:03d}" for index in range(reads))


class PairedHoldoutRunnerTests(unittest.TestCase):
    """Exercise actual files, input adapters and recovery; all inputs are synthetic."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.policy = MethylationMixturePolicy(
            profile_id="synthetic-paired-holdout-v1",
            source_a_fractions=[0.0, 0.5, 1.0],
            replicates=10,
            calibration_read_group_fraction=0.25,
            minimum_markers=6,
            maximum_markers=6,
            minimum_reference_valid_calls=5,
            minimum_mixture_valid_calls=3,
            minimum_absolute_delta_beta=0.5,
            uncertainty_draws=100,
            maximum_confidence_interval_width=1.0,
            note="Synthetic software regression only; no clinical or biological validation.",
        )
        self.matrix = MethylationValidationMatrix(
            matrix_id="synthetic-paired-holdout-matrix-v1",
            source_a_fractions=[0.0, 0.5, 1.0],
            read_group_budgets=[20],
            replicate_seeds=list(range(11, 21)),
            estimator_policy=self.policy,
            acceptance=[self._acceptance(20)],
            rationale="Idealized synthetic whole-read recovery with prospective thresholds.",
        )
        self.split_policy = PairedHoldoutSplitPolicy(
            calibration_read_group_fraction=0.25,
            source_a_split_seed=101,
            source_b_split_seed=202,
        )
        self.metadata = NanopolishSourceMetadata(
            caller_name="nanopolish",
            caller_version="0.14.0",
            reference_id="SYNTHETIC_REFERENCE",
            reference_genome_build=GenomeBuild.GRCH38,
            reference_sha256=hashlib.sha256(b"SYNTHETIC_NANOPOLISH_REFERENCE").hexdigest(),
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
        for role in HOLDOUT_SOURCE_ROLES:
            path = self.directory / f"{role}.tsv"
            self._write_nanopolish(path, role)
            metadata_path = self.directory / f"{role}.metadata.json"
            metadata_path.write_text(self.metadata.model_dump_json(), encoding="utf-8")
            bindings[role] = ValidationLocalInput(
                calls_path=path,
                input_format="nanopolish",
                metadata_path=metadata_path,
            )
            samples[role] = ValidationSample(
                accession=f"SYNTHETIC_ACCESSION_{role}",
                dataset_version="synthetic-v1",
                biological_sample_id=f"SYNTHETIC_SAMPLE_{role}",
                # Same donor is allowed; these are two material sources, not four donors.
                donor_id="SYNTHETIC_SHARED_DONOR",
                biological_source=f"SYNTHETIC_BIOLOGICAL_{role}",
                cell_type=f"SYNTHETIC_CELL_{role}",
                material_type="SYNTHETIC_DNA",
                input_sha256=sha256_file(path),
                provenance_reference="Generated entirely within this synthetic software test.",
                technical=technical,
            )
        self.inputs = PairedHoldoutLocalInputs(samples=bindings)
        self.cohort = PairedHoldoutCohort(
            samples=samples,
            source_distinction_evidence="Two synthetic biological source labels.",
            access_basis="public",
        )
        self.registration = self._register()

    @staticmethod
    def _acceptance(budget: int) -> RecoveryAcceptance:
        return RecoveryAcceptance(
            read_group_budget=budget,
            maximum_absolute_bias=0.01,
            maximum_mae=0.01,
            maximum_rmse=0.01,
            minimum_conditional_interval_coverage=0.9,
            maximum_no_call_rate=0.0,
            maximum_per_fraction_absolute_bias=0.01,
            minimum_per_fraction_conditional_interval_coverage=0.9,
            maximum_per_fraction_no_call_rate=0.0,
        )

    def _register(self, *, unseen: bool = True) -> PairedHoldoutRegistration:
        return preregister_paired_holdout(
            self.matrix,
            self.cohort,
            self.split_policy,
            registration_id="SYNTHETIC_HOLDOUT_REGISTRATION",
            registered_at=datetime.now(UTC) - timedelta(minutes=1),
            test_outcomes_unseen=unseen,
            code_sha256=validation_software_sha256(),
            software_version=__version__,
        )

    def _cli_files(self) -> dict[str, Path]:
        models = {
            "matrix": self.matrix,
            "split-policy": self.split_policy,
            "cohort": self.cohort,
            "inputs": self.inputs,
            "registration": self.registration,
        }
        paths = {}
        for name, model in models.items():
            path = self.directory / f"cli-{name}.json"
            path.write_text(model.model_dump_json(indent=2), encoding="utf-8")
            paths[name] = path
        return paths

    @staticmethod
    def _cli(*arguments: str | Path) -> str:
        output = io.StringIO()
        with (
            patch("sys.argv", ["ontseq", *(str(argument) for argument in arguments)]),
            redirect_stdout(output),
        ):
            entrypoint.main()
        return output.getvalue()

    def _write_nanopolish(
        self,
        path: Path,
        role: str,
        *,
        reverse_test: bool = False,
        reads: int = 40,
    ) -> None:
        names = _names(role, reads)
        seed = 101 if role == "source_a" else 202
        test_names = set(_split_reads(names, fraction=0.25, seed=seed)[1]) if reads > 1 else set()
        lines = [
            "chromosome\tstrand\tstart\tend\tread_name\tlog_lik_ratio\t"
            "log_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_motifs\tsequence"
        ]
        for name in names:
            ratio = 5.0 if role == "source_a" else -5.0
            if reverse_test and name in test_names:
                ratio = -ratio
            for marker in range(6):
                position = 100 + marker * 10
                lines.append(
                    f"chr1\t+\t{position}\t{position}\t{name}\t{ratio}\t{ratio}\t0\t1\t1\tACGTA"
                )
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def _refresh_input_locks(self) -> None:
        for role in HOLDOUT_SOURCE_ROLES:
            self.cohort.samples[role].input_sha256 = sha256_file(
                self.inputs.samples[role].calls_path
            )
        self.registration = self._register()

    def _set_budget(self, budget: int) -> None:
        self.matrix.read_group_budgets = [budget]
        self.matrix.acceptance = [self._acceptance(budget)]
        self.registration = self._register()

    def _use_dorado(self, *, modkit: bool = False) -> None:
        reference = self.directory / "synthetic.fa"
        reference.write_bytes(b">chr1\nCGCGCGCGCGCG\n")
        adapter_policy = ModbamAdapterPolicy()
        policy_path = self.directory / "synthetic.adapter.json"
        policy_path.write_text(adapter_policy.model_dump_json(), encoding="utf-8")
        for role in HOLDOUT_SOURCE_ROLES:
            # Synthetic lineage fingerprints identify fixture parents, not external genomes.
            parent = (
                hashlib.sha256(f"SYNTHETIC_PARENT_{role}".encode()).hexdigest() if modkit else None
            )
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
                modkit_version="0.6.1" if modkit else None,
                source_modbam_sha256=parent,
            )
            technical = ValidationTechnicalMetadata(
                **metadata.model_dump(exclude={"schema_version", "source_modbam_sha256"}),
                adapter_name="modkit_extract_full_tsv" if modkit else "sam_mm_ml",
                adapter_version=MODBAM_ADAPTER_VERSION,
                adapter_policy_sha256=content_sha256(adapter_policy),
            )
            path = self.directory / f"{role}.{'modkit.tsv' if modkit else 'sam'}"
            if modkit:
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(
                        [
                            "read_id",
                            "forward_read_position",
                            "ref_position",
                            "chrom",
                            "mod_strand",
                            "ref_strand",
                            "ref_mod_strand",
                            "read_length",
                            "mod_qual",
                            "mod_code",
                            "canonical_base",
                            "modified_primary_base",
                            "inferred",
                            "flag",
                            "alignment_start",
                            "alignment_end",
                        ]
                    )
                    for name in _names(role):
                        for position in range(0, 12, 2):
                            writer.writerow(
                                [
                                    name,
                                    position,
                                    position,
                                    "chr1",
                                    "+",
                                    "+",
                                    "+",
                                    12,
                                    "0.9980469" if role == "source_a" else "0.001953125",
                                    "m",
                                    "C",
                                    "C",
                                    "false",
                                    0,
                                    0,
                                    12,
                                ]
                            )
            else:
                ml = ",".join(["255" if role == "source_a" else "0"] * 6)
                lines = ["@SQ\tSN:chr1\tLN:12", "@PG\tID:basecaller\tPN:dorado\tVN:synthetic-1"]
                lines.extend(
                    f"{name}\t0\tchr1\t1\t60\t12M\t*\t0\t0\tCGCGCGCGCGCG\t*"
                    f"\tMM:Z:C+m?,0,0,0,0,0,0;\tML:B:C,{ml}"
                    for name in _names(role)
                )
                path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            metadata_path = self.directory / f"{role}.dorado.json"
            metadata_path.write_text(metadata.model_dump_json(), encoding="utf-8")
            self.inputs.samples[role] = ValidationLocalInput(
                calls_path=path,
                input_format="modkit" if modkit else "modbam",
                metadata_path=metadata_path,
                reference_path=reference,
                adapter_policy_path=policy_path,
            )
            self.cohort.samples[role].technical = technical
            self.cohort.samples[role].source_modbam_sha256 = parent
        self._refresh_input_locks()

    def test_real_nanopolish_files_recover_30_levels_with_bounded_scope(self) -> None:
        first = execute_registered_paired_holdout(self.registration, self.inputs)
        second = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(first.decision, ValidationDecision.PASS)
        self.assertEqual(first.overall_metrics.completed_levels, 30)
        self.assertEqual(first.overall_metrics.no_call_rate, 0)
        self.assertAlmostEqual(first.overall_metrics.mean_absolute_error, 0)
        self.assertAlmostEqual(first.overall_metrics.root_mean_squared_error, 0)
        self.assertEqual(first.validation_scope, "held_out_read_group_recovery_within_two_sources")
        self.assertEqual(first.donor_independent_validation, "NOT_EVALUATED")
        self.assertEqual(first.estimand, "fraction_of_retained_whole_read_groups_from_source_a")
        self.assertTrue(first.research_only)
        self.assertIsNotNone(first.evidence)
        self.assertIsNotNone(second.evidence)
        self.assertEqual(first.evidence.levels, second.evidence.levels)
        self.assertEqual(
            first.evidence.split_manifest_sha256, second.evidence.split_manifest_sha256
        )
        self.assertEqual(first.evidence.source_read_group_counts, {"source_a": 40, "source_b": 40})
        self.assertEqual(
            first.evidence.calibration_read_group_counts, {"source_a": 10, "source_b": 10}
        )
        self.assertEqual(first.evidence.test_read_group_counts, {"source_a": 30, "source_b": 30})
        self.assertEqual(len(first.evidence.reference_markers), 6)
        self.assertEqual(
            PairedHoldoutRecoveryReport.model_validate_json(first.model_dump_json()), first
        )
        self.assertNotIn("SYN_READ_source_a_000", first.model_dump_json())

    def test_actual_estimator_pools_are_disjoint_whole_reads(self) -> None:
        calibration = []
        tests = []
        actual_aggregate = methylation_holdout_runner._aggregate
        actual_generate = methylation_holdout_runner.generate_validation_levels

        def observe_aggregate(source, read_names):
            calibration.append(dict(source.calls_by_read))
            self.assertEqual(set(read_names), set(source.calls_by_read))
            return actual_aggregate(source, read_names)

        def observe_generate(matrix, source_a, source_b, markers):
            tests.extend([dict(source_a.calls_by_read), dict(source_b.calls_by_read)])
            return actual_generate(matrix, source_a, source_b, markers)

        with (
            patch.object(methylation_holdout_runner, "_aggregate", side_effect=observe_aggregate),
            patch.object(
                methylation_holdout_runner,
                "generate_validation_levels",
                side_effect=observe_generate,
            ),
        ):
            report = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.PASS)
        self.assertEqual(len(calibration), 2)
        self.assertEqual(len(tests), 2)
        for index, role in enumerate(HOLDOUT_SOURCE_ROLES):
            seed = 101 if role == "source_a" else 202
            expected_calibration, expected_test = _split_reads(
                _names(role), fraction=0.25, seed=seed
            )
            self.assertEqual(set(calibration[index]), set(expected_calibration))
            self.assertEqual(set(tests[index]), set(expected_test))
            self.assertFalse(set(calibration[index]) & set(tests[index]))
            self.assertEqual(set(calibration[index]) | set(tests[index]), set(_names(role)))
            self.assertTrue(
                all(
                    len(calls) == 6
                    for calls in (*calibration[index].values(), *tests[index].values())
                )
            )
            suffix = role[-1]
            for pool_role, names in (
                (f"calibration_{suffix}", expected_calibration),
                (f"test_{suffix}", expected_test),
            ):
                self.assertEqual(
                    report.evidence.pool_selection_sha256[pool_role],
                    _selection_digest(names, (), context=f"paired-holdout:{pool_role}"),
                )

    def test_test_only_biology_does_not_change_calibration_markers(self) -> None:
        baseline = execute_registered_paired_holdout(self.registration, self.inputs)
        for role in HOLDOUT_SOURCE_ROLES:
            self._write_nanopolish(self.inputs.samples[role].calls_path, role, reverse_test=True)
        self._refresh_input_locks()
        reversed_report = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(
            reversed_report.evidence.reference_markers, baseline.evidence.reference_markers
        )
        self.assertEqual(
            reversed_report.evidence.reference_marker_sha256,
            baseline.evidence.reference_marker_sha256,
        )
        self.assertEqual(
            reversed_report.evidence.pool_selection_sha256, baseline.evidence.pool_selection_sha256
        )
        self.assertEqual(reversed_report.decision, ValidationDecision.FAIL)
        self.assertGreater(reversed_report.overall_metrics.mean_absolute_error, 0.6)

    def test_sam_mm_ml_executes_the_two_source_holdout_matrix(self) -> None:
        self._use_dorado()
        report = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.PASS)
        self.assertEqual(report.overall_metrics.completed_levels, 30)
        self.assertEqual(report.evidence.calibration_read_group_counts["source_a"], 10)
        self.assertEqual(report.evidence.test_read_group_counts["source_a"], 30)

    def test_modkit_061_rows_execute_the_two_source_holdout_matrix(self) -> None:
        self._use_dorado(modkit=True)
        report = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.PASS)
        self.assertEqual(report.overall_metrics.completed_levels, 30)
        self.assertEqual(report.evidence.source_read_group_counts, {"source_a": 40, "source_b": 40})

    def test_low_read_budget_has_explicit_no_call_at_every_matrix_cell(self) -> None:
        self._set_budget(2)
        report = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.explicit_no_call_levels, 30)
        self.assertEqual(report.overall_metrics.no_call_rate, 1)

    def test_budget_larger_than_fixed_test_pool_cannot_silently_resample(self) -> None:
        self._set_budget(32)
        report = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertIn("holdout pools", report.execution_no_call_reason)
        self.assertIsNone(report.evidence)

    def test_input_and_metadata_changes_after_registration_are_rejected(self) -> None:
        path = self.inputs.samples["source_a"].calls_path
        original = path.read_bytes()
        path.write_bytes(original + b"\n")
        with self.assertRaisesRegex(ValueError, "registered SHA-256"):
            execute_registered_paired_holdout(self.registration, self.inputs)
        path.write_bytes(original)
        altered = self.metadata.model_dump()
        altered["basecaller_version"] = "2.0"
        self.inputs.samples["source_b"].metadata_path.write_text(
            json.dumps(altered), encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "metadata differs at basecaller_version"):
            execute_registered_paired_holdout(self.registration, self.inputs)

    def test_runtime_change_before_or_during_execution_is_rejected(self) -> None:
        with (
            patch.object(
                methylation_holdout_runner, "validation_software_sha256", return_value="f" * 64
            ),
            self.assertRaisesRegex(ValueError, "Software/runtime changed after"),
        ):
            execute_registered_paired_holdout(self.registration, self.inputs)
        with (
            patch.object(
                methylation_holdout_runner,
                "validation_software_sha256",
                side_effect=[
                    self.registration.code_sha256,
                    self.registration.code_sha256,
                    "f" * 64,
                ],
            ),
            self.assertRaisesRegex(ValueError, "Software/runtime changed during"),
        ):
            execute_registered_paired_holdout(self.registration, self.inputs)

    def test_source_change_during_estimation_is_detected_by_final_hash(self) -> None:
        actual_generate = methylation_holdout_runner.generate_validation_levels

        def modify_source_after_estimation(*args, **kwargs):
            result = actual_generate(*args, **kwargs)
            path = self.inputs.samples["source_a"].calls_path
            path.write_bytes(path.read_bytes() + b"\n")
            return result

        with (
            patch.object(
                methylation_holdout_runner,
                "generate_validation_levels",
                side_effect=modify_source_after_estimation,
            ),
            self.assertRaisesRegex(ValueError, "input changed during"),
        ):
            execute_registered_paired_holdout(self.registration, self.inputs)

    def test_reused_read_identifiers_across_distinct_files_are_rejected(self) -> None:
        path = self.inputs.samples["source_b"].calls_path
        path.write_text(
            path.read_text(encoding="utf-8").replace("SYN_READ_source_b_", "SYN_READ_source_a_"),
            encoding="utf-8",
        )
        self._refresh_input_locks()
        with self.assertRaisesRegex(ValueError, "Read identifiers overlap"):
            execute_registered_paired_holdout(self.registration, self.inputs)

    def test_split_manifest_is_bound_to_actual_source_and_pool_counts(self) -> None:
        report = execute_registered_paired_holdout(self.registration, self.inputs)
        evidence = report.evidence
        self.assertEqual(
            evidence.split_manifest_sha256,
            holdout_split_manifest_sha256(
                registration_sha256=self.registration.lock_sha256,
                observed_input_sha256=evidence.observed_input_sha256,
                source_read_group_counts=evidence.source_read_group_counts,
                calibration_read_group_counts=evidence.calibration_read_group_counts,
                test_read_group_counts=evidence.test_read_group_counts,
                pool_selection_sha256=evidence.pool_selection_sha256,
            ),
        )
        payload = report.model_dump(mode="json")
        payload["evidence"]["pool_selection_sha256"]["test_a"] = "f" * 64
        with self.assertRaisesRegex(ValueError, "split manifest"):
            PairedHoldoutRecoveryReport.model_validate(payload)

    def test_registration_split_mutation_is_rejected(self) -> None:
        self.registration.split_policy.source_a_split_seed += 1
        with self.assertRaises(ValueError):
            execute_registered_paired_holdout(self.registration, self.inputs)

    def test_registration_of_previous_unbounded_nanopolish_adapter_is_rejected(self) -> None:
        self.cohort.samples["source_a"].technical.adapter_version = "nanopolish-call-table-v1"
        self.cohort.samples["source_b"].technical.adapter_version = "nanopolish-call-table-v1"
        self.registration = self._register()
        with self.assertRaisesRegex(ValueError, "registered adapter version"):
            execute_registered_paired_holdout(self.registration, self.inputs)

    def test_seen_test_outcomes_cannot_execute_or_pass(self) -> None:
        self.registration = self._register(unseen=False)
        with patch.object(methylation_holdout_runner, "parse_validation_source") as parser:
            report = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        parser.assert_not_called()

    def test_identical_biological_source_declarations_cannot_pass(self) -> None:
        self.cohort.samples["source_b"].biological_source = self.cohort.samples[
            "source_a"
        ].biological_source
        self.registration = self._register()
        report = execute_registered_paired_holdout(self.registration, self.inputs)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertIn("biologically distinct", " ".join(report.reasons))

    def test_cli_plan_is_no_call_and_refuses_to_overwrite_existing_output(self) -> None:
        paths = self._cli_files()
        output_path = self.directory / "readiness.json"
        arguments = (
            "methylation-holdout-plan",
            "--matrix",
            paths["matrix"],
            "--split-policy",
            paths["split-policy"],
            "--output",
            output_path,
        )
        console = self._cli(*arguments)
        original = output_path.read_bytes()
        readiness = json.loads(original)
        self.assertIn("NO_CALL", console)
        self.assertEqual(readiness["decision"], "NO_CALL")
        self.assertEqual(readiness["status"], "NOT_REGISTERED")
        self.assertEqual(readiness["expected_levels"], 30)
        self.assertEqual(readiness["matrix_sha256"], content_sha256(self.matrix))
        self.assertEqual(readiness["split_policy_sha256"], content_sha256(self.split_policy))
        self.assertEqual(readiness["donor_independent_validation"], "NOT_EVALUATED")
        with self.assertRaisesRegex(SystemExit, "ERROR:"):
            self._cli(*arguments)
        self.assertEqual(output_path.read_bytes(), original)

    def test_cli_fingerprint_catalogues_two_roles_without_local_paths(self) -> None:
        paths = self._cli_files()
        output_path = self.directory / "fingerprints.json"
        self._cli(
            "methylation-holdout-fingerprint",
            "--inputs",
            paths["inputs"],
            "--output",
            output_path,
        )
        encoded = output_path.read_text(encoding="utf-8")
        catalogue = json.loads(encoded)
        self.assertEqual(set(catalogue["samples"]), set(HOLDOUT_SOURCE_ROLES))
        self.assertEqual(catalogue["code_sha256"], validation_software_sha256())
        self.assertTrue(catalogue["sensitive_output"])
        self.assertNotIn(self.directory.name, encoded)
        self.assertNotIn("calls_path", encoded)
        self.assertNotIn("metadata_path", encoded)
        for role, binding in self.inputs.samples.items():
            sample = catalogue["samples"][role]
            self.assertEqual(sample["input_format"], "nanopolish")
            self.assertEqual(set(sample["files"]), {"calls", "metadata"})
            for label, path in (("calls", binding.calls_path), ("metadata", binding.metadata_path)):
                self.assertEqual(
                    sample["files"][label],
                    {"sha256": sha256_file(path), "size_bytes": path.stat().st_size},
                )

    def test_cli_registration_and_execution_export_all_levels_with_bounded_scope(self) -> None:
        paths = self._cli_files()
        registration_path = self.directory / "cli-prospective-registration.json"
        registration_id = "SYNTHETIC_CLI_HOLDOUT"
        self._cli(
            "methylation-holdout-register",
            "--matrix",
            paths["matrix"],
            "--split-policy",
            paths["split-policy"],
            "--cohort",
            paths["cohort"],
            "--registration-id",
            registration_id,
            "--confirm-test-outcomes-unseen",
            "--output",
            registration_path,
        )
        registered = PairedHoldoutRegistration.model_validate_json(
            registration_path.read_text(encoding="utf-8")
        )
        self.assertEqual(registered.code_sha256, validation_software_sha256())
        self.assertEqual(registered.split_policy, self.split_policy)
        self.assertTrue(registered.test_outcomes_unseen)
        output_directory = self.directory / "cli-results"
        console = self._cli(
            "methylation-holdout-run",
            "--registration",
            registration_path,
            "--inputs",
            paths["inputs"],
            "--output-dir",
            output_directory,
        )
        stem = registration_id + ".methylation-holdout"
        outputs = {
            extension: output_directory / f"{stem}.{extension}"
            for extension in ("json", "csv", "html")
        }
        self.assertEqual(set(output_directory.iterdir()), set(outputs.values()))
        report = PairedHoldoutRecoveryReport.model_validate_json(
            outputs["json"].read_text(encoding="utf-8")
        )
        self.assertEqual(report.decision, ValidationDecision.PASS)
        self.assertEqual(len(report.evidence.levels), 30)
        with outputs["csv"].open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 30)
        self.assertEqual({int(row["replicate_seed"]) for row in rows}, set(range(11, 21)))
        self.assertEqual({float(row["target_source_a_fraction"]) for row in rows}, {0, 0.5, 1})
        for row in rows:
            self.assertEqual(row["validation_decision"], "PASS")
            self.assertEqual(row["validation_scope"], report.validation_scope)
            self.assertEqual(row["donor_independent_validation"], "NOT_EVALUATED")
        page = outputs["html"].read_text(encoding="utf-8")
        self.assertIn("Research use only.", page)
        self.assertIn("Fraction of retained whole read groups from source A.", page)
        self.assertIn("Donor-independent validation: NOT_EVALUATED", page)
        self.assertIn("Only retained reads from these two sources.", page)
        self.assertIn("Donor-independent validation: NOT_EVALUATED", console)

    def test_cli_evaluation_without_evidence_reports_all_missing_as_no_call(self) -> None:
        paths = self._cli_files()
        output_directory = self.directory / "cli-missing-evidence"
        console = self._cli(
            "methylation-holdout-evaluate",
            "--registration",
            paths["registration"],
            "--output-dir",
            output_directory,
        )
        stem = self.registration.registration_id + ".methylation-holdout"
        report = PairedHoldoutRecoveryReport.model_validate_json(
            (output_directory / f"{stem}.json").read_text(encoding="utf-8")
        )
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertIsNone(report.evidence)
        self.assertEqual(report.overall_metrics.expected_levels, 30)
        self.assertEqual(report.overall_metrics.completed_levels, 0)
        self.assertEqual(report.overall_metrics.missing_levels, 30)
        self.assertEqual(report.overall_metrics.no_call_rate, 1)
        self.assertEqual(report.donor_independent_validation, "NOT_EVALUATED")
        with (output_directory / f"{stem}.csv").open(encoding="utf-8", newline="") as handle:
            self.assertEqual(list(csv.DictReader(handle)), [])
        self.assertIn("NO_CALL", console)
        self.assertIn("Donor-independent validation: NOT_EVALUATED", console)


if __name__ == "__main__":
    unittest.main()
