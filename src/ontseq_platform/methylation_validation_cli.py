"""Explicit local commands for prospective methylation validation."""

from __future__ import annotations

import argparse
import csv
import html
import json
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from . import __version__
from .methylation_holdout import PairedHoldoutRecoveryReport
from .methylation_validation import (
    MethylationValidationEvidence,
    MethylationValidationMatrix,
    MethylationValidationRegistration,
    MethylationValidationReport,
    ValidationCohort,
    content_sha256,
    evaluate_methylation_validation,
    preregister_validation,
)
from .methylation_validation_runner import (
    ValidationLocalInput,
    ValidationLocalInputs,
    ValidationReadiness,
    execute_registered_validation,
    validation_software_sha256,
)
from .models import StrictModel
from .reference import sha256_file

COMMANDS = (
    "methylation-validation-plan",
    "methylation-validation-fingerprint",
    "methylation-validation-register",
    "methylation-validation-run",
    "methylation-validation-evaluate",
)


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_new(payload: StrictModel | dict[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = (
        payload.model_dump_json(indent=2)
        if isinstance(payload, StrictModel)
        else json.dumps(payload, indent=2, allow_nan=False)
    )
    with path.open("x", encoding="utf-8") as handle:
        handle.write(text + "\n")
    return path


def render_validation_report(
    report: MethylationValidationReport | PairedHoldoutRecoveryReport,
    output_dir: Path,
) -> list[Path]:
    report = type(report).model_validate_json(report.model_dump_json())
    paired = isinstance(report, PairedHoldoutRecoveryReport)
    suffix = ".methylation-holdout" if paired else ".methylation-validation"
    stem = report.registration.registration_id + suffix
    paths = [output_dir / f"{stem}.{extension}" for extension in ("json", "csv", "html")]
    if any(path.exists() for path in paths):
        raise ValueError("Validation output already exists; choose a new output directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_new(report, paths[0])
    rows: list[dict[str, Any]] = []
    if report.evidence is not None:
        for item in report.evidence.levels:
            level = item.result
            rows.append(
                {
                    "read_group_budget": item.read_group_budget,
                    "replicate_seed": item.replicate_seed,
                    "target_source_a_fraction": level.target_source_a_fraction,
                    "realized_source_a_read_group_fraction": (
                        level.realized_source_a_read_group_fraction
                    ),
                    "status": level.status.value,
                    "estimated_source_a_fraction": level.estimated_source_a_fraction,
                    "confidence_interval_lower": level.confidence_interval_lower,
                    "confidence_interval_upper": level.confidence_interval_upper,
                    "no_call_reason": level.no_call_reason,
                    "validation_decision": report.decision.value,
                    "validation_scope": report.validation_scope,
                    "study_matrix_id": report.registration.matrix.matrix_id,
                    "donor_independent_validation": (
                        "NOT_EVALUATED" if paired else "source_design_gate_applied"
                    ),
                }
            )
    fields = [
        "read_group_budget",
        "replicate_seed",
        "target_source_a_fraction",
        "realized_source_a_read_group_fraction",
        "status",
        "estimated_source_a_fraction",
        "confidence_interval_lower",
        "confidence_interval_upper",
        "no_call_reason",
        "validation_decision",
        "validation_scope",
        "study_matrix_id",
        "donor_independent_validation",
    ]
    with paths[1].open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    metrics = report.overall_metrics.model_dump()
    metric_rows = "".join(
        f"<tr><th>{html.escape(key)}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in metrics.items()
    )
    stratum_rows = "".join(
        f"<tr><td>{stratum.read_group_budget}</td>"
        f"<td>{stratum.source_a_fraction if stratum.source_a_fraction is not None else 'all'}</td>"
        f"<td>{stratum.decision.value}</td>"
        f"<td>{stratum.metrics.completed_levels}/{stratum.metrics.expected_levels}</td>"
        f"<td>{stratum.metrics.no_call_rate:.3f}</td>"
        + "".join(
            f"<td>{value:.4f}</td>" if value is not None else "<td>unavailable</td>"
            for value in (
                stratum.metrics.mean_bias,
                stratum.metrics.mean_absolute_error,
                stratum.metrics.root_mean_squared_error,
                stratum.metrics.conditional_interval_coverage,
            )
        )
        + f"<td>{html.escape('; '.join(stratum.reasons))}</td></tr>"
        for stratum in report.strata
    )
    reasons = "".join(f"<li>{html.escape(value)}</li>" for value in report.reasons)
    limitations = "".join(f"<li>{html.escape(value)}</li>" for value in report.limitations)
    title = (
        "Paired-source held-out read recovery" if paired else "Independent methylation validation"
    )
    scope_note = (
        "Donor-independent validation: NOT_EVALUATED. Only retained reads from these two sources."
        if paired
        else "Donor-independent study; in-silico retained-read recovery only."
    )
    page = (
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="ontseq-sensitive-output" content="true">'
        f"<title>{title}</title>"
        "<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:16px}"
        "td,th{padding:8px;text-align:left;border-bottom:1px solid #ccc}"
        "table{border-collapse:collapse;width:100%}</style>"
        f"<h1>{title}</h1><p><strong>Research use only.</strong> "
        f"Fraction of retained whole read groups from source A.</p><p>{scope_note}</p>"
        f"<p>Study: {html.escape(report.registration.matrix.matrix_id)}</p>"
        f"<p>{html.escape(report.registration.matrix.rationale)}</p>"
        f"<h2>{report.decision.value}</h2><ul>{reasons}</ul><table>{metric_rows}</table>"
        "<h2>Registered strata</h2><table><tr><th>Read groups</th><th>Source A</th>"
        "<th>Decision</th><th>Completed/planned</th><th>NO_CALL rate</th><th>Bias</th>"
        "<th>MAE</th><th>RMSE</th><th>Coverage</th>"
        f"<th>Reasons</th></tr>{stratum_rows}</table>"
        f"<h2>Interpretation limits</h2><ul>{limitations}</ul></html>"
    )
    with paths[2].open("x", encoding="utf-8") as handle:
        handle.write(page)
    return paths


def fingerprint_local_inputs(
    samples: Iterable[tuple[str, ValidationLocalInput]],
) -> dict[str, Any]:
    """Catalogue actual local bytes without exposing paths or parsing biological outcomes."""
    catalogue: dict[str, Any] = {
        "schema_version": "0.1.0",
        "sensitive_output": True,
        "created_at": datetime.now(UTC).isoformat(),
        "samples": {},
        "code_sha256": validation_software_sha256(),
        "software_version": __version__,
    }
    for role, item in samples:
        files: dict[str, Any] = {}
        for label in ("calls_path", "metadata_path", "reference_path", "adapter_policy_path"):
            path = getattr(item, label)
            if path is not None:
                files[label.removesuffix("_path")] = {
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
        catalogue["samples"][role] = {"input_format": item.input_format, "files": files}
    return catalogue


def add_subparsers(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    plan = subparsers.add_parser(COMMANDS[0])
    plan.add_argument("--matrix", required=True, type=Path)
    plan.add_argument("--output", required=True, type=Path)
    fingerprint = subparsers.add_parser(COMMANDS[1])
    fingerprint.add_argument("--inputs", required=True, type=Path)
    fingerprint.add_argument("--output", required=True, type=Path)
    register = subparsers.add_parser(COMMANDS[2])
    register.add_argument("--matrix", required=True, type=Path)
    register.add_argument("--cohort", required=True, type=Path)
    register.add_argument("--registration-id", required=True)
    register.add_argument("--confirm-test-outcomes-unseen", action="store_true", required=True)
    register.add_argument("--output", required=True, type=Path)
    run = subparsers.add_parser(COMMANDS[3])
    run.add_argument("--registration", required=True, type=Path)
    run.add_argument("--inputs", required=True, type=Path)
    run.add_argument("--output-dir", required=True, type=Path)
    evaluate = subparsers.add_parser(COMMANDS[4])
    evaluate.add_argument("--registration", required=True, type=Path)
    evaluate.add_argument("--evidence", type=Path)
    evaluate.add_argument("--output-dir", required=True, type=Path)


def run_command(args: argparse.Namespace) -> None:
    try:
        if args.command == COMMANDS[0]:
            matrix = MethylationValidationMatrix.model_validate(_load(args.matrix))
            readiness = ValidationReadiness(
                matrix_sha256=content_sha256(matrix),
                code_sha256=validation_software_sha256(),
                expected_levels=matrix.expected_levels,
            )
            print(_write_new(readiness, args.output))
            print("NO_CALL: study design only; no independent experimental evidence")
        elif args.command == COMMANDS[1]:
            inputs = ValidationLocalInputs.model_validate(_load(args.inputs))
            catalogue = fingerprint_local_inputs(inputs.samples.items())
            print(_write_new(catalogue, args.output))
        elif args.command == COMMANDS[2]:
            matrix = MethylationValidationMatrix.model_validate(_load(args.matrix))
            cohort = ValidationCohort.model_validate(_load(args.cohort))
            registration = preregister_validation(
                matrix,
                cohort,
                registration_id=args.registration_id,
                registered_at=datetime.now(UTC),
                test_outcomes_unseen=True,
                code_sha256=validation_software_sha256(),
                software_version=__version__,
            )
            print(_write_new(registration, args.output))
        else:
            registration = MethylationValidationRegistration.model_validate(
                _load(args.registration)
            )
            if args.command == COMMANDS[3]:
                inputs = ValidationLocalInputs.model_validate(_load(args.inputs))
                report = execute_registered_validation(registration, inputs)
            else:
                evidence = (
                    MethylationValidationEvidence.model_validate(_load(args.evidence))
                    if (args.evidence is not None)
                    else None
                )
                report = evaluate_methylation_validation(registration, evidence)
            for output in render_validation_report(report, args.output_dir):
                print(output)
            print(f"{report.decision.value}: retained-read-group recovery only")
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc


def main() -> None:
    parser = argparse.ArgumentParser(description="Research-only independent methylation recovery")
    add_subparsers(parser.add_subparsers(dest="command", required=True))
    run_command(parser.parse_args())


if __name__ == "__main__":
    main()
