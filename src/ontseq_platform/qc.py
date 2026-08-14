from __future__ import annotations

import json
from typing import Any

from .execution import CommandRunner, SubprocessRunner
from .models import (
    CraminoQCReport,
    InputKind,
    QCMetrics,
    QCPolicy,
    SampleManifest,
    ToolRecord,
    Verdict,
)


def _mapping(value: object, field: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"Cramino JSON field {field!r} must be an object")
    return value


def _number(mapping: dict[str, Any], field: str) -> int | float:
    value = mapping.get(field)
    if not isinstance(value, int | float) or isinstance(value, bool):
        raise ValueError(f"Cramino JSON field {field!r} must be numeric")
    return value


def parse_cramino_json(text: str) -> dict[str, float | int | str | None]:
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError("Cramino did not return valid JSON") from exc
    root = _mapping(payload, "root")
    alignment = _mapping(root.get("alignment_stats"), "alignment_stats")
    read_stats = _mapping(root.get("read_stats"), "read_stats")
    metrics: dict[str, float | int | str | None] = {
        "number_of_reads": int(_number(alignment, "num_reads")),
        "number_of_alignments": int(_number(alignment, "num_alignments")),
        "aligned_percent": float(_number(alignment, "percent_from_total")),
        "total_yield_gb": float(_number(read_stats, "yield_gb")),
        "mean_coverage_x": float(_number(read_stats, "mean_coverage")),
        "n50_bp": int(_number(read_stats, "n50")),
        "n75_bp": int(_number(read_stats, "n75")),
        "median_length_bp": float(_number(read_stats, "median_length")),
        "mean_length_bp": float(_number(read_stats, "mean_length")),
    }
    identity_value = root.get("identity_stats")
    if identity_value is not None:
        identity = _mapping(identity_value, "identity_stats")
        metrics.update(
            {
                "median_identity_percent": float(_number(identity, "median_identity")),
                "mean_identity_percent": float(_number(identity, "mean_identity")),
                "modal_identity_percent": float(_number(identity, "modal_identity")),
                "identity_estimated": str(bool(identity.get("is_estimated", False))).lower(),
            }
        )
    return metrics


def evaluate_qc_metrics(
    metrics: dict[str, float | int | str | None], policy: QCPolicy
) -> QCMetrics:
    failed_gates: list[str] = []
    warnings: list[str] = []
    configured = {name: value for name, value in policy.numeric_gates.items() if value is not None}
    for name, minimum in configured.items():
        value = metrics.get(name)
        if not isinstance(value, int | float) or isinstance(value, bool):
            failed_gates.append(name)
            warnings.append(f"Configured QC gate {name!r} could not be evaluated")
        elif value < minimum:
            failed_gates.append(name)
            warnings.append(f"QC metric {name!r} is below its configured minimum")
    if failed_gates:
        verdict = Verdict.FAIL
    elif not configured:
        verdict = Verdict.WARN
        warnings.append(
            "No validated numeric QC gates are configured; metrics are descriptive only."
        )
    elif policy.status != "validated":
        verdict = Verdict.WARN
    else:
        verdict = Verdict.PASS
    if policy.status != "validated":
        warnings.append(policy.note)
    return QCMetrics(
        verdict=verdict,
        metrics=metrics,
        warnings=warnings,
        failed_gates=failed_gates,
    )


def run_cramino_qc(
    manifest: SampleManifest,
    policy: QCPolicy,
    *,
    runner: CommandRunner | None = None,
    cramino: str = "cramino",
    threads: int = 4,
) -> CraminoQCReport:
    if manifest.input.kind != InputKind.ALIGNED_BAM:
        raise ValueError("Cramino aligned-BAM QC requires input.kind=aligned_bam")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    command_runner = runner or SubprocessRunner()
    version_result = command_runner.run([cramino, "--version"], timeout_seconds=30)
    if version_result.returncode != 0:
        raise ValueError("Cramino version probe returned a non-zero exit code")
    version_line = (
        version_result.stdout.splitlines()[0] if version_result.stdout else "cramino unknown"
    )
    version = version_line.removeprefix("cramino ").strip()
    result = command_runner.run(
        [cramino, "--threads", str(threads), "--format", "json", manifest.input.path],
        timeout_seconds=3600,
    )
    if result.returncode != 0:
        raise ValueError(f"Cramino failed with exit code {result.returncode}")
    metrics = parse_cramino_json(result.stdout)
    qc = evaluate_qc_metrics(metrics, policy)
    return CraminoQCReport(
        sample_id=manifest.sample_id,
        qc=qc,
        tool=ToolRecord(
            name="cramino",
            version=version,
            parameters={"threads": threads, "format": "json"},
        ),
        limitations=[
            "Cramino metrics are descriptive until assay-specific QC gates are validated.",
            "aligned_percent follows Cramino percent_from_total semantics for default "
            "aligned mode.",
            "No read-level data or source path is copied into this normalized QC report.",
        ],
    )
