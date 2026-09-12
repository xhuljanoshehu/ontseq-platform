"""CLI for a prospective two-source holdout study, distinct from donor validation."""

from __future__ import annotations

import argparse
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import Field

from . import __version__
from .methylation_holdout import (
    PairedHoldoutCohort,
    PairedHoldoutEvidence,
    PairedHoldoutRegistration,
    PairedHoldoutSplitPolicy,
    evaluate_paired_holdout,
    preregister_paired_holdout,
)
from .methylation_holdout_runner import (
    PairedHoldoutLocalInputs,
    execute_registered_paired_holdout,
)
from .methylation_validation import MethylationValidationMatrix, content_sha256
from .methylation_validation_cli import (
    _load,
    _write_new,
    fingerprint_local_inputs,
    render_validation_report,
)
from .methylation_validation_runner import validation_software_sha256
from .models import StrictModel

COMMANDS = (
    "methylation-holdout-plan",
    "methylation-holdout-fingerprint",
    "methylation-holdout-register",
    "methylation-holdout-run",
    "methylation-holdout-evaluate",
)


class PairedHoldoutReadiness(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    study_design: Literal["paired_source_read_holdout"] = "paired_source_read_holdout"
    decision: Literal["NO_CALL"] = "NO_CALL"
    status: Literal["NOT_REGISTERED"] = "NOT_REGISTERED"
    matrix_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    split_policy_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    code_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    expected_levels: int = Field(ge=1)
    software_version: str = __version__
    donor_independent_validation: Literal["NOT_EVALUATED"] = "NOT_EVALUATED"
    reasons: list[str] = Field(
        default_factory=lambda: [
            "No checksummed two-source cohort and unseen test experiment have been registered.",
            "No real-data recovery results are available; this is a prospective study template.",
        ]
    )
    research_only: Literal[True] = True
    sensitive_output: Literal[True] = True


def add_subparsers(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    plan = subparsers.add_parser(COMMANDS[0])
    plan.add_argument("--matrix", type=Path, required=True)
    plan.add_argument("--split-policy", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    fingerprint = subparsers.add_parser(COMMANDS[1])
    fingerprint.add_argument("--inputs", type=Path, required=True)
    fingerprint.add_argument("--output", type=Path, required=True)
    register = subparsers.add_parser(COMMANDS[2])
    register.add_argument("--matrix", type=Path, required=True)
    register.add_argument("--split-policy", type=Path, required=True)
    register.add_argument("--cohort", type=Path, required=True)
    register.add_argument("--registration-id", required=True)
    register.add_argument("--confirm-test-outcomes-unseen", action="store_true", required=True)
    register.add_argument("--output", type=Path, required=True)
    run = subparsers.add_parser(COMMANDS[3])
    run.add_argument("--registration", type=Path, required=True)
    run.add_argument("--inputs", type=Path, required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    evaluate = subparsers.add_parser(COMMANDS[4])
    evaluate.add_argument("--registration", type=Path, required=True)
    evaluate.add_argument("--evidence", type=Path)
    evaluate.add_argument("--output-dir", type=Path, required=True)


def run_command(args: argparse.Namespace) -> None:
    try:
        if args.command in COMMANDS[:3]:
            if args.command == COMMANDS[1]:
                inputs = PairedHoldoutLocalInputs.model_validate(_load(args.inputs))
                print(_write_new(fingerprint_local_inputs(inputs.samples.items()), args.output))
                return
            matrix = MethylationValidationMatrix.model_validate(_load(args.matrix))
            split = PairedHoldoutSplitPolicy.model_validate(_load(args.split_policy))
            if split.calibration_read_group_fraction != (
                matrix.estimator_policy.calibration_read_group_fraction
            ):
                raise ValueError("Split and matrix calibration fractions disagree")
            if args.command == COMMANDS[0]:
                readiness = PairedHoldoutReadiness(
                    matrix_sha256=content_sha256(matrix),
                    split_policy_sha256=content_sha256(split),
                    code_sha256=validation_software_sha256(),
                    expected_levels=matrix.expected_levels,
                )
                print(_write_new(readiness, args.output))
                print("NO_CALL: prospective two-source design; no experimental evidence")
                return
            cohort = PairedHoldoutCohort.model_validate(_load(args.cohort))
            registration = preregister_paired_holdout(
                matrix,
                cohort,
                split,
                registration_id=args.registration_id,
                registered_at=datetime.now(UTC),
                test_outcomes_unseen=True,
                code_sha256=validation_software_sha256(),
                software_version=__version__,
            )
            print(_write_new(registration, args.output))
            return
        registration = PairedHoldoutRegistration.model_validate(_load(args.registration))
        if args.command == COMMANDS[3]:
            inputs = PairedHoldoutLocalInputs.model_validate(_load(args.inputs))
            report = execute_registered_paired_holdout(registration, inputs)
        else:
            evidence = (
                PairedHoldoutEvidence.model_validate(_load(args.evidence))
                if (args.evidence is not None)
                else None
            )
            report = evaluate_paired_holdout(registration, evidence)
        for output in render_validation_report(report, args.output_dir):
            print(output)
        print(f"{report.decision.value}: held-out read-group recovery in these two sources only")
        print("Donor-independent validation: NOT_EVALUATED")
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc


def main() -> None:
    parser = argparse.ArgumentParser(description="Prospective two-source held-out read recovery")
    add_subparsers(parser.add_subparsers(dest="command", required=True))
    run_command(parser.parse_args())


if __name__ == "__main__":
    main()
