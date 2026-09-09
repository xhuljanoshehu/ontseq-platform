from __future__ import annotations

import csv
import json
import os
import tempfile
from pathlib import Path
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


def read_length_histogram(text: str) -> list[tuple[int, int | None, int, int]]:
    """Extract numeric Cramino read-length bins when histogram output is present."""

    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError("Cramino did not return valid JSON") from exc
    root = _mapping(payload, "root")
    raw_histograms = root.get("histograms")
    if raw_histograms is None:
        return []
    histograms = _mapping(raw_histograms, "histograms")
    raw_read_length = histograms.get("read_length")
    if raw_read_length is None:
        return []
    # Cramino 1.3.0 serializes a Histogram object with step/max_value/bins. Retain the
    # earlier direct-array shape as an import-compatible form for already captured fixtures.
    raw_bins = raw_read_length.get("bins") if isinstance(raw_read_length, dict) else raw_read_length
    if not isinstance(raw_bins, list):
        raise ValueError("Cramino JSON field 'histograms.read_length.bins' must be an array")
    bins: list[tuple[int, int | None, int, int]] = []
    for index, raw_bin in enumerate(raw_bins):
        if not isinstance(raw_bin, dict):
            raise ValueError(f"Cramino read-length histogram bin {index} must be an object")
        start = raw_bin.get("start")
        end = raw_bin.get("end")
        count = raw_bin.get("count")
        bases = raw_bin.get("bases")
        if (
            not isinstance(start, int)
            or isinstance(start, bool)
            or (end is not None and (not isinstance(end, int) or isinstance(end, bool)))
            or not isinstance(count, int)
            or isinstance(count, bool)
            or not isinstance(bases, int)
            or isinstance(bases, bool)
            or start < 0
            or (end is not None and end <= start)
            or count < 0
            or bases < 0
        ):
            raise ValueError(f"Cramino read-length histogram bin {index} is not numeric/valid")
        bins.append((start, end, count, bases))
    return bins


def write_read_length_histogram(
    bins: list[tuple[int, int | None, int, int]], output_path: Path
) -> Path | None:
    """Write a compact numeric TSV sidecar; absence remains explicit as no artifact."""

    if not bins:
        # This helper may be called while re-running QC in an existing envelope.  Empty
        # current evidence must remove an earlier histogram instead of preserving stale
        # numeric data that the new QC result no longer supports.
        output_path.unlink(missing_ok=True)
        return None
    output_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["start_bp", "end_bp", "read_count", "base_count"])
            writer.writerows(
                (start, "" if end is None else end, count, bases)
                for start, end, count, bases in bins
            )
        os.replace(temporary, output_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return output_path


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
    else:
        verdict = Verdict.PASS
    return QCMetrics(
        verdict=verdict,
        metrics=metrics,
        warnings=warnings,
        failed_gates=failed_gates,
    )


def cramino_version(stdout: str) -> str:
    """Parse a Cramino version from its probe output.

    Public because preflight has to reach the same answer this module will. A preflight
    that parsed versions differently from the run it precedes could report a tool as
    identified that the run then cannot identify.
    """
    version_line = stdout.splitlines()[0] if stdout else "cramino unknown"
    return version_line.removeprefix("cramino ").strip()


def _stderr_tail(stderr: str, *, limit: int = 4000) -> str:
    """Return a bounded diagnostic tail without flooding logs on tool failure."""
    return stderr.strip()[-limit:]


def run_cramino_qc(
    manifest: SampleManifest,
    policy: QCPolicy,
    *,
    runner: CommandRunner | None = None,
    cramino: str = "cramino",
    threads: int = 4,
    histogram_output: Path | None = None,
) -> CraminoQCReport:
    if manifest.input.kind != InputKind.ALIGNED_BAM:
        raise ValueError("Cramino aligned-BAM QC requires input.kind=aligned_bam")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    command_runner = runner or SubprocessRunner()
    version_result = command_runner.run([cramino, "--version"], timeout_seconds=30)
    if version_result.returncode != 0:
        diagnostic = _stderr_tail(version_result.stderr)
        message = "Cramino version probe returned a non-zero exit code"
        if diagnostic:
            message += f": {diagnostic}"
        raise ValueError(message)
    version = cramino_version(version_result.stdout)
    argv = [cramino, "--threads", str(threads)]
    raw_histogram: Path | None = None
    if histogram_output is not None:
        histogram_output.parent.mkdir(parents=True, exist_ok=True)
        raw_histogram = histogram_output.with_name(f".{histogram_output.name}.cramino.tsv")
        raw_histogram.unlink(missing_ok=True)
        # Without an explicit FILE, cramino 1.3.0 writes the --hist-count TSV to stdout and
        # interleaves it with --format json. Give the optional argument its own file so stdout
        # remains one valid JSON document; the normalized sidecar below still comes from the
        # structured JSON histogram bins.
        argv.extend(["--hist-count", str(raw_histogram)])
    argv.extend(["--format", "json", manifest.input.path])
    try:
        result = command_runner.run(argv, timeout_seconds=3600)
        if result.returncode != 0:
            diagnostic = _stderr_tail(result.stderr)
            message = f"Cramino failed with exit code {result.returncode}"
            if diagnostic:
                message += f": {diagnostic}"
            raise ValueError(message)
        metrics = parse_cramino_json(result.stdout)
        if histogram_output is not None:
            write_read_length_histogram(read_length_histogram(result.stdout), histogram_output)
    finally:
        if raw_histogram is not None:
            raw_histogram.unlink(missing_ok=True)
    qc = evaluate_qc_metrics(metrics, policy)
    return CraminoQCReport(
        sample_id=manifest.sample_id,
        qc=qc,
        tool=ToolRecord(
            name="cramino",
            version=version,
            parameters={
                "threads": threads,
                "format": "json",
                "read_length_histogram_requested": histogram_output is not None,
            },
        ),
        limitations=[
            "Cramino metrics are descriptive until assay-specific QC gates are validated.",
            "aligned_percent follows Cramino percent_from_total semantics for default "
            "aligned mode.",
            "No read-level data or source path is copied into this normalized QC report.",
        ],
    )
