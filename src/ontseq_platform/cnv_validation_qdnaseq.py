from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path, PurePosixPath
from typing import Any

from .cnv.qdnaseq import CnvFit, QDNAseqCallReport, QDNAseqPolicy
from .cnv_validation_contracts import CnvCallerLock, CnvValidationSpecimen
from .cnv_validation_evidence import (
    CnvContributionStatus,
    CnvEvidenceRecordKind,
    CnvFullEvidenceRecord,
    CnvNativeArtifactReference,
    CnvNormalizedEventRecord,
    CnvNumericMeasurement,
    CnvRunOutcomeState,
    CnvValidationEvidenceManifest,
    canonical_evidence_sha256,
    seal_cnv_validation_evidence,
)
from .models import Locus, ModuleRunStatus
from .pipeline.envelope import sha256_file


def _stable_id(prefix: str, *parts: object) -> str:
    rendered = "\x1f".join(str(part) for part in parts)
    digest = hashlib.sha256(rendered.encode("utf-8")).hexdigest()[:24]
    return f"{prefix}-{digest}"


def _artifact_role_and_media_type(relative_path: str) -> tuple[str, str]:
    name = PurePosixPath(relative_path).name
    if name.endswith(".qdnaseq-ace.summary.json"):
        return "run_summary", "application/json"
    if name.endswith(".consensus.chromosomes.tsv"):
        return "chromosome_consensus", "text/tab-separated-values"
    if name.endswith(".segments.tsv"):
        return "caller_segments", "text/tab-separated-values"
    if name.endswith(".chromosomes.tsv"):
        return "chromosome_summary", "text/tab-separated-values"
    if name.endswith(".bins.tsv"):
        return "caller_bins", "text/tab-separated-values"
    if name.endswith(".ace-models.tsv"):
        return "ace_models", "text/tab-separated-values"
    if name.endswith(".ace-fit.png"):
        return "ace_fit_plot", "image/png"
    if name.endswith(".copy-number.png"):
        return "copy_number_plot", "image/png"
    if name.endswith(".segmented.rds"):
        return "caller_rds", "application/octet-stream"
    return "caller_artifact", "application/octet-stream"


def _artifact_path(output_dir: Path, relative_path: str) -> Path:
    if "\\" in relative_path:
        raise ValueError("QDNAseq evidence artifact paths must use POSIX separators")
    posix = PurePosixPath(relative_path)
    if not relative_path or posix.is_absolute() or ".." in posix.parts:
        raise ValueError(f"Unsafe QDNAseq evidence artifact path: {relative_path!r}")

    root = output_dir.resolve()
    path = (output_dir / Path(*posix.parts)).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"QDNAseq evidence artifact escapes output directory: {relative_path!r}")
    if not path.is_file():
        raise ValueError(f"QDNAseq evidence artifact is missing: {relative_path}")
    return path


def _native_artifacts(
    report: QDNAseqCallReport,
    output_dir: Path,
) -> list[CnvNativeArtifactReference]:
    if len(report.output_files) != len(set(report.output_files)):
        raise ValueError("QDNAseq report output_files must be unique")

    artifacts: list[CnvNativeArtifactReference] = []
    for relative_path in report.output_files:
        path = _artifact_path(output_dir, relative_path)
        digest = sha256_file(path)
        role, media_type = _artifact_role_and_media_type(relative_path)
        artifacts.append(
            CnvNativeArtifactReference(
                artifact_id=_stable_id("qdnaseq-artifact", relative_path, digest),
                role=role,
                relative_path=relative_path,
                sha256=digest,
                size_bytes=path.stat().st_size,
                media_type=media_type,
                description="Retained caller-native QDNAseq+ACE validation artifact.",
            )
        )
    return artifacts


def _dependency_versions(report: QDNAseqCallReport) -> dict[str, str]:
    names = [tool.name for tool in report.tools]
    if len(names) != len(set(names)):
        raise ValueError("QDNAseq report contains duplicate tool-version records")
    versions = {tool.name: tool.version for tool in report.tools}
    for required in ("QDNAseq", "ACE"):
        if not versions.get(required):
            raise ValueError(f"QDNAseq validation evidence requires {required} version provenance")
    return versions


def _caller_version(versions: dict[str, str]) -> str:
    return f"QDNAseq={versions['QDNAseq']};ACE={versions['ACE']}"


def _fit_artifact_paths(fit: CnvFit) -> list[str]:
    result = [
        fit.segment_file,
        fit.chromosome_file,
        fit.fit_plot,
        fit.copy_number_plot,
        fit.rds_file,
    ]
    if fit.bins_file is not None:
        result.append(fit.bins_file)
    if fit.model_file is not None:
        result.append(fit.model_file)
    return result


def _fit_artifact_ids(
    fit: CnvFit,
    artifact_id_by_path: dict[str, str],
) -> list[str]:
    missing = [path for path in _fit_artifact_paths(fit) if path not in artifact_id_by_path]
    if missing:
        raise ValueError(
            "QDNAseq fit references native artifacts absent from output_files: "
            + ", ".join(sorted(missing))
        )
    return [artifact_id_by_path[path] for path in _fit_artifact_paths(fit)]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"QDNAseq validation TSV is missing or empty: {path.name}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"QDNAseq validation TSV has no header: {path.name}")
        return [dict(row) for row in reader]


def _finite_float(row: dict[str, str], key: str) -> float:
    raw = row.get(key, "")
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid QDNAseq numeric value for {key}: {raw!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite QDNAseq numeric value for {key}: {raw!r}")
    return value


def _finite_int(row: dict[str, str], key: str) -> int:
    value = _finite_float(row, key)
    if not value.is_integer():
        raise ValueError(f"Invalid QDNAseq integer value for {key}: {row.get(key)!r}")
    return int(value)


def _optional_measurement(
    row: dict[str, str],
    key: str,
    *,
    unit: str | None = None,
) -> CnvNumericMeasurement | None:
    raw = row.get(key)
    if raw is None or raw == "":
        return None
    return CnvNumericMeasurement(
        name=key,
        value=_finite_float(row, key),
        unit=unit,
        source_field=key,
    )


def _chromosome(row: dict[str, str]) -> str:
    value = row.get("chromosome") or row.get("chr")
    if not value:
        raise ValueError("QDNAseq validation row is missing chromosome")
    return value if value.startswith("chr") else f"chr{value}"


def _locus(row: dict[str, str]) -> Locus:
    if row.get("coordinate_system") != "zero_based_half_open":
        raise ValueError("QDNAseq validation row must declare zero_based_half_open coordinates")
    start = _finite_int(row, "start")
    end = _finite_int(row, "end")
    if start < 0 or end <= start:
        raise ValueError("QDNAseq validation row contains an invalid half-open interval")
    return Locus(chromosome=_chromosome(row), start=start, end=end)


def _resolution_contribution(
    fit: CnvFit,
    policy: QDNAseqPolicy,
) -> CnvContributionStatus:
    if fit.bin_size_kbp == policy.primary_bin_size_kbp:
        return CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS
    return CnvContributionStatus.SECONDARY


def _bin_records(
    *,
    fit: CnvFit,
    policy: QDNAseqPolicy,
    common: dict[str, Any],
    output_dir: Path,
    artifact_id_by_path: dict[str, str],
) -> list[CnvFullEvidenceRecord]:
    if fit.bins_file is None:
        return []
    artifact_id = artifact_id_by_path.get(fit.bins_file)
    if artifact_id is None:
        raise ValueError("QDNAseq bin table is absent from retained native artifacts")

    rows = _read_tsv(_artifact_path(output_dir, fit.bins_file))
    result: list[CnvFullEvidenceRecord] = []
    measurement_fields = (
        ("reads", "count"),
        ("bases", "count"),
        ("gc", None),
        ("mappability", None),
        ("blacklist", None),
        ("residual", None),
        ("loess", None),
        ("corrected", None),
        ("copynumber", None),
        ("segmented", None),
    )
    for index, row in enumerate(rows, start=1):
        measurements = [
            measurement
            for key, unit in measurement_fields
            if (measurement := _optional_measurement(row, key, unit=unit)) is not None
        ]
        result.append(
            CnvFullEvidenceRecord.model_validate(
                {
                    **common,
                    "record_id": f"qdnaseq-bin-{fit.bin_size_kbp}-{index}",
                    "bin_size_kbp": fit.bin_size_kbp,
                    "record_kind": CnvEvidenceRecordKind.CALLER_BIN,
                    "native_record_id": f"bin-{fit.bin_size_kbp}-{index}",
                    "run_outcome": CnvRunOutcomeState.OBSERVED,
                    "contribution_status": _resolution_contribution(fit, policy),
                    "native_artifact_ids": [artifact_id],
                    "fit_group_id": f"ace-fit-{fit.bin_size_kbp}",
                    "primary": _locus(row),
                    "numeric_measurements": measurements,
                }
            )
        )
    return result


def _segment_records(
    *,
    fit: CnvFit,
    policy: QDNAseqPolicy,
    common: dict[str, Any],
    output_dir: Path,
    artifact_id_by_path: dict[str, str],
) -> list[CnvFullEvidenceRecord]:
    artifact_id = artifact_id_by_path.get(fit.segment_file)
    if artifact_id is None:
        raise ValueError("QDNAseq segment table is absent from retained native artifacts")

    rows = _read_tsv(_artifact_path(output_dir, fit.segment_file))
    result: list[CnvFullEvidenceRecord] = []
    for index, row in enumerate(rows, start=1):
        measurements = [
            CnvNumericMeasurement(
                name="bin_count",
                value=float(_finite_int(row, "bin_count")),
                unit="count",
                source_field="bin_count",
            ),
            CnvNumericMeasurement(
                name="call",
                value=_finite_float(row, "call"),
                source_field="call",
            ),
        ]
        raw_qnorm = row.get("qnorm_log10")
        if raw_qnorm not in {None, "", "-Inf"}:
            measurements.append(
                CnvNumericMeasurement(
                    name="qnorm_log10",
                    value=_finite_float(row, "qnorm_log10"),
                    source_field="qnorm_log10",
                )
            )

        result.append(
            CnvFullEvidenceRecord.model_validate(
                {
                    **common,
                    "record_id": f"qdnaseq-segment-{fit.bin_size_kbp}-{index}",
                    "bin_size_kbp": fit.bin_size_kbp,
                    "record_kind": CnvEvidenceRecordKind.CALLER_SEGMENT,
                    "native_record_id": f"segment-{fit.bin_size_kbp}-{index}",
                    "run_outcome": CnvRunOutcomeState.OBSERVED,
                    "contribution_status": _resolution_contribution(fit, policy),
                    "native_artifact_ids": [artifact_id],
                    "fit_group_id": f"ace-fit-{fit.bin_size_kbp}",
                    "primary": _locus(row),
                    "copy_number": _finite_float(row, "absolute_copy_number"),
                    "numeric_measurements": measurements,
                }
            )
        )
    return result


def _chromosome_records(
    *,
    fit: CnvFit,
    policy: QDNAseqPolicy,
    common: dict[str, Any],
    output_dir: Path,
    artifact_id_by_path: dict[str, str],
) -> list[CnvFullEvidenceRecord]:
    artifact_id = artifact_id_by_path.get(fit.chromosome_file)
    if artifact_id is None:
        raise ValueError("QDNAseq chromosome table is absent from retained native artifacts")

    rows = _read_tsv(_artifact_path(output_dir, fit.chromosome_file))
    seen: set[str] = set()
    result: list[CnvFullEvidenceRecord] = []
    for row in rows:
        chromosome = _chromosome(row)
        if chromosome in seen:
            raise ValueError(f"QDNAseq chromosome table contains duplicate row for {chromosome}")
        seen.add(chromosome)
        result.append(
            CnvFullEvidenceRecord.model_validate(
                {
                    **common,
                    "record_id": f"qdnaseq-chromosome-{fit.bin_size_kbp}-{chromosome}",
                    "bin_size_kbp": fit.bin_size_kbp,
                    "record_kind": CnvEvidenceRecordKind.CHROMOSOME_SUMMARY,
                    "native_record_id": f"chromosome-{fit.bin_size_kbp}-{chromosome}",
                    "run_outcome": CnvRunOutcomeState.OBSERVED,
                    "contribution_status": _resolution_contribution(fit, policy),
                    "native_artifact_ids": [artifact_id],
                    "fit_group_id": f"ace-fit-{fit.bin_size_kbp}",
                    "copy_number": _finite_float(row, "copy_number"),
                }
            )
        )
    return result


def _consensus_records(
    *,
    report: QDNAseqCallReport,
    common: dict[str, Any],
    artifact_id_by_path: dict[str, str],
) -> list[CnvFullEvidenceRecord]:
    if not report.chromosome_consensus:
        return []

    paths = [path for path in artifact_id_by_path if path.endswith(".consensus.chromosomes.tsv")]
    if len(paths) != 1:
        raise ValueError("QDNAseq chromosome consensus requires exactly one retained consensus TSV")
    artifact_id = artifact_id_by_path[paths[0]]

    result: list[CnvFullEvidenceRecord] = []
    for item in report.chromosome_consensus:
        result.append(
            CnvFullEvidenceRecord.model_validate(
                {
                    **common,
                    "record_id": f"qdnaseq-consensus-{item.chromosome}",
                    "bin_size_kbp": None,
                    "record_kind": CnvEvidenceRecordKind.CHROMOSOME_SUMMARY,
                    "native_record_id": f"consensus-{item.chromosome}",
                    "run_outcome": CnvRunOutcomeState.OBSERVED,
                    "contribution_status": CnvContributionStatus.SECONDARY,
                    "native_artifact_ids": [artifact_id],
                    "copy_number": item.median_copy_number,
                    "numeric_measurements": [
                        CnvNumericMeasurement(
                            name="rounded_copy_number",
                            value=float(item.rounded_copy_number),
                            source_field="rounded_copy_number",
                        ),
                        CnvNumericMeasurement(
                            name="agreeing_bins",
                            value=float(item.agreeing_bins),
                            unit="count",
                            source_field="agreeing_bins",
                        ),
                        CnvNumericMeasurement(
                            name="contributing_bins",
                            value=float(item.contributing_bins),
                            unit="count",
                            source_field="contributing_bins",
                        ),
                        CnvNumericMeasurement(
                            name="min_copy_number",
                            value=item.min_copy_number,
                            source_field="min_copy_number",
                        ),
                        CnvNumericMeasurement(
                            name="max_copy_number",
                            value=item.max_copy_number,
                            source_field="max_copy_number",
                        ),
                    ],
                }
            )
        )
    return result


def _common_record_payload(
    *,
    registration_sha256: str,
    specimen: CnvValidationSpecimen,
    caller_lock: CnvCallerLock,
    replicate_id: str,
    caller_parameters: dict[str, Any],
    caller_parameters_sha256: str,
    dependency_versions: dict[str, str],
) -> dict[str, Any]:
    return {
        "registration_sha256": registration_sha256,
        "specimen_id": specimen.specimen_id,
        "biological_specimen_id": specimen.biological_specimen_id,
        "caller_id": caller_lock.caller_id,
        "caller_version": caller_lock.caller_version,
        "adapter_policy_sha256": caller_lock.adapter_policy_sha256,
        "execution_identity_sha256": caller_lock.execution_identity_sha256,
        "genome_build": specimen.genome_build,
        "data_basis": specimen.data_basis,
        "reference_id": specimen.reference_id,
        "reference_sha256": specimen.reference_sha256,
        "input_sha256": specimen.input_sha256,
        "coverage_x": specimen.coverage_x,
        "coverage_definition": specimen.coverage_definition,
        "tumor_fraction": specimen.tumor_fraction,
        "tumor_fraction_method": specimen.tumor_fraction_method,
        "tumor_fraction_timepoint": specimen.tumor_fraction_timepoint,
        "repeat_kind": specimen.repeat_kind,
        "replicate_id": replicate_id,
        "repeat_group_id": specimen.repeat_group_id,
        "caller_parameters": caller_parameters,
        "caller_parameters_sha256": caller_parameters_sha256,
        "dependency_versions": dependency_versions,
    }


def _validate_adapter_identity(
    *,
    report: QDNAseqCallReport,
    policy: QDNAseqPolicy,
    specimen: CnvValidationSpecimen,
    caller_lock: CnvCallerLock,
    dependency_versions: dict[str, str],
) -> None:
    if caller_lock.caller_id != "qdnaseq_ace":
        raise ValueError("QDNAseq validation adapter requires caller_id='qdnaseq_ace'")
    if report.sample_id != specimen.specimen_id:
        raise ValueError("QDNAseq report sample does not match validation specimen")
    if report.genome_build != specimen.genome_build:
        raise ValueError("QDNAseq report build does not match validation specimen")

    observed_bins = [fit.bin_size_kbp for fit in report.fits]
    if len(observed_bins) != len(set(observed_bins)):
        raise ValueError("QDNAseq report contains duplicate bin-size fits")
    if set(observed_bins) != set(policy.bin_sizes_kbp):
        raise ValueError("QDNAseq report bin sizes do not match the registered adapter policy")
    if report.primary_fit.bin_size_kbp != policy.primary_bin_size_kbp:
        raise ValueError("QDNAseq report primary fit does not match the adapter policy")

    observed_caller_version = _caller_version(dependency_versions)
    if caller_lock.caller_version != observed_caller_version:
        raise ValueError(
            "QDNAseq caller lock version does not match report provenance: "
            f"{observed_caller_version}"
        )


def _selected_fit_record(
    *,
    fit: CnvFit,
    policy: QDNAseqPolicy,
    common: dict[str, Any],
    native_artifact_ids: list[str],
) -> CnvFullEvidenceRecord:
    primary = fit.bin_size_kbp == policy.primary_bin_size_kbp
    return CnvFullEvidenceRecord.model_validate(
        {
            **common,
            "record_id": f"qdnaseq-fit-{fit.bin_size_kbp}-selected",
            "bin_size_kbp": fit.bin_size_kbp,
            "record_kind": CnvEvidenceRecordKind.CALLER_FIT,
            "native_record_id": f"ace-fit-{fit.bin_size_kbp}-selected",
            "run_outcome": CnvRunOutcomeState.OBSERVED,
            "contribution_status": (
                CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS
                if primary
                else CnvContributionStatus.SECONDARY
            ),
            "native_artifact_ids": native_artifact_ids,
            "fit_group_id": f"ace-fit-{fit.bin_size_kbp}",
            "selected_fit": True,
            "cellularity": fit.cellularity,
            "ploidy": fit.ploidy,
            "fit_error": fit.fit_error,
            "numeric_measurements": [
                CnvNumericMeasurement(
                    name="candidate_count",
                    value=float(fit.candidate_count),
                    unit="count",
                    source_field="candidate_count",
                ),
                CnvNumericMeasurement(
                    name="segment_count",
                    value=float(fit.segment_count),
                    unit="count",
                    source_field="segment_count",
                ),
            ],
        }
    )


def _alternative_fit_records(
    *,
    fit: CnvFit,
    common: dict[str, Any],
    native_artifact_ids: list[str],
) -> list[CnvFullEvidenceRecord]:
    required = {"cellularity", "ploidy", "fit_error"}
    result: list[CnvFullEvidenceRecord] = []
    for index, alternative in enumerate(fit.alternatives, start=1):
        if set(alternative) != required:
            raise ValueError(
                "ACE alternative fit must contain exactly cellularity, ploidy and fit_error"
            )
        result.append(
            CnvFullEvidenceRecord.model_validate(
                {
                    **common,
                    "record_id": f"qdnaseq-fit-{fit.bin_size_kbp}-alternative-{index}",
                    "bin_size_kbp": fit.bin_size_kbp,
                    "record_kind": CnvEvidenceRecordKind.CALLER_FIT,
                    "native_record_id": f"ace-fit-{fit.bin_size_kbp}-alternative-{index}",
                    "run_outcome": CnvRunOutcomeState.OBSERVED,
                    "contribution_status": CnvContributionStatus.SECONDARY,
                    "native_artifact_ids": native_artifact_ids,
                    "fit_group_id": f"ace-fit-{fit.bin_size_kbp}",
                    "selected_fit": False,
                    "cellularity": alternative["cellularity"],
                    "ploidy": alternative["ploidy"],
                    "fit_error": alternative["fit_error"],
                }
            )
        )
    return result


def map_qdnaseq_ace_full_evidence(
    *,
    report: QDNAseqCallReport,
    policy: QDNAseqPolicy,
    output_dir: Path,
    registration_sha256: str,
    specimen: CnvValidationSpecimen,
    caller_lock: CnvCallerLock,
    replicate_id: str,
) -> tuple[list[CnvNativeArtifactReference], list[CnvFullEvidenceRecord]]:
    """Map retained QDNAseq+ACE multi-resolution results into Block-2 full evidence.

    This adapter is additive: it fingerprints caller-native files and emits immutable
    evidence rows. It does not alter QDNAseq/ACE execution, choose a new primary
    resolution, normalize events, or discard non-primary fits.
    """
    artifacts = _native_artifacts(report, output_dir)
    artifact_id_by_path = {item.relative_path: item.artifact_id for item in artifacts}
    dependency_versions = _dependency_versions(report)
    _validate_adapter_identity(
        report=report,
        policy=policy,
        specimen=specimen,
        caller_lock=caller_lock,
        dependency_versions=dependency_versions,
    )

    caller_parameters = policy.model_dump(mode="json")
    caller_parameters_sha256 = canonical_evidence_sha256(caller_parameters)
    common = _common_record_payload(
        registration_sha256=registration_sha256,
        specimen=specimen,
        caller_lock=caller_lock,
        replicate_id=replicate_id,
        caller_parameters=caller_parameters,
        caller_parameters_sha256=caller_parameters_sha256,
        dependency_versions=dependency_versions,
    )

    if report.status == ModuleRunStatus.COMPLETED:
        run_outcome = CnvRunOutcomeState.OBSERVED
        contribution_status = CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS
        outcome_reason = None
    elif report.status == ModuleRunStatus.NO_CALL:
        run_outcome = CnvRunOutcomeState.NO_CALL
        contribution_status = CnvContributionStatus.NO_CALL
        outcome_reason = "QDNAseq+ACE completed without normalized CNV events."
    else:
        raise ValueError(
            f"Unsupported QDNAseq report status for validation evidence: {report.status}"
        )

    full_evidence: list[CnvFullEvidenceRecord] = [
        CnvFullEvidenceRecord.model_validate(
            {
                **common,
                "record_id": "qdnaseq-run-summary",
                "bin_size_kbp": None,
                "record_kind": CnvEvidenceRecordKind.RUN_SUMMARY,
                "native_record_id": "qdnaseq-ace-run-summary",
                "run_outcome": run_outcome,
                "contribution_status": contribution_status,
                "outcome_reason": outcome_reason,
                "native_artifact_ids": [item.artifact_id for item in artifacts],
                "numeric_measurements": [
                    CnvNumericMeasurement(
                        name="configured_bin_count",
                        value=float(len(policy.bin_sizes_kbp)),
                        unit="count",
                        source_field="policy.bin_sizes_kbp",
                    ),
                    CnvNumericMeasurement(
                        name="normalized_event_count",
                        value=float(len(report.events)),
                        unit="count",
                        source_field="report.events",
                    ),
                ],
            }
        )
    ]

    for fit in sorted(report.fits, key=lambda item: item.bin_size_kbp):
        fit_artifact_ids = _fit_artifact_ids(fit, artifact_id_by_path)
        full_evidence.append(
            _selected_fit_record(
                fit=fit,
                policy=policy,
                common=common,
                native_artifact_ids=fit_artifact_ids,
            )
        )
        full_evidence.extend(
            _alternative_fit_records(
                fit=fit,
                common=common,
                native_artifact_ids=fit_artifact_ids,
            )
        )
        full_evidence.extend(
            _bin_records(
                fit=fit,
                policy=policy,
                common=common,
                output_dir=output_dir,
                artifact_id_by_path=artifact_id_by_path,
            )
        )
        full_evidence.extend(
            _segment_records(
                fit=fit,
                policy=policy,
                common=common,
                output_dir=output_dir,
                artifact_id_by_path=artifact_id_by_path,
            )
        )
        full_evidence.extend(
            _chromosome_records(
                fit=fit,
                policy=policy,
                common=common,
                output_dir=output_dir,
                artifact_id_by_path=artifact_id_by_path,
            )
        )

    full_evidence.extend(
        _consensus_records(
            report=report,
            common=common,
            artifact_id_by_path=artifact_id_by_path,
        )
    )
    return artifacts, full_evidence


def _normalized_common_payload(
    *,
    registration_sha256: str,
    specimen: CnvValidationSpecimen,
    caller_lock: CnvCallerLock,
    replicate_id: str,
) -> dict[str, Any]:
    return {
        "registration_sha256": registration_sha256,
        "specimen_id": specimen.specimen_id,
        "biological_specimen_id": specimen.biological_specimen_id,
        "caller_id": caller_lock.caller_id,
        "caller_version": caller_lock.caller_version,
        "adapter_policy_sha256": caller_lock.adapter_policy_sha256,
        "execution_identity_sha256": caller_lock.execution_identity_sha256,
        "normalization_policy_sha256": caller_lock.adapter_policy_sha256,
        "genome_build": specimen.genome_build,
        "data_basis": specimen.data_basis,
        "reference_id": specimen.reference_id,
        "reference_sha256": specimen.reference_sha256,
        "input_sha256": specimen.input_sha256,
        "coverage_x": specimen.coverage_x,
        "coverage_definition": specimen.coverage_definition,
        "tumor_fraction": specimen.tumor_fraction,
        "tumor_fraction_method": specimen.tumor_fraction_method,
        "tumor_fraction_timepoint": specimen.tumor_fraction_timepoint,
        "repeat_kind": specimen.repeat_kind,
        "replicate_id": replicate_id,
        "repeat_group_id": specimen.repeat_group_id,
    }


def _normalized_event_records(
    *,
    report: QDNAseqCallReport,
    policy: QDNAseqPolicy,
    registration_sha256: str,
    specimen: CnvValidationSpecimen,
    caller_lock: CnvCallerLock,
    replicate_id: str,
    full_evidence: list[CnvFullEvidenceRecord],
) -> list[CnvNormalizedEventRecord]:
    primary_segments = [
        record
        for record in full_evidence
        if record.record_kind == CnvEvidenceRecordKind.CALLER_SEGMENT
        and record.bin_size_kbp == policy.primary_bin_size_kbp
        and record.run_outcome == CnvRunOutcomeState.OBSERVED
    ]
    common = _normalized_common_payload(
        registration_sha256=registration_sha256,
        specimen=specimen,
        caller_lock=caller_lock,
        replicate_id=replicate_id,
    )

    result: list[CnvNormalizedEventRecord] = []
    for event in report.events:
        sources = [
            record
            for record in primary_segments
            if record.primary == event.primary
            and (
                event.copy_number is None
                or (
                    record.copy_number is not None
                    and math.isclose(
                        record.copy_number,
                        event.copy_number,
                        rel_tol=0.0,
                        abs_tol=1e-9,
                    )
                )
            )
        ]
        if len(sources) != 1:
            raise ValueError(
                "Each normalized QDNAseq event must map to exactly one primary "
                f"caller segment; {event.event_id} mapped to {len(sources)}"
            )
        source = sources[0]
        result.append(
            CnvNormalizedEventRecord.model_validate(
                {
                    **common,
                    "normalized_event_id": _stable_id(
                        "qdnaseq-normalized",
                        event.event_id,
                        event.primary.chromosome,
                        event.primary.start,
                        event.primary.end,
                    ),
                    "bin_size_kbp": policy.primary_bin_size_kbp,
                    "source_full_evidence_ids": [source.record_id],
                    "event_type": event.event_type,
                    "primary": event.primary,
                    "normalized_copy_number": event.copy_number,
                    "contribution_status": CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
                }
            )
        )
    return result


def build_qdnaseq_ace_validation_manifest(
    *,
    manifest_id: str,
    report: QDNAseqCallReport,
    policy: QDNAseqPolicy,
    output_dir: Path,
    registration_sha256: str,
    specimen: CnvValidationSpecimen,
    caller_lock: CnvCallerLock,
    replicate_id: str,
) -> CnvValidationEvidenceManifest:
    """Build and seal complete retained QDNAseq+ACE validation evidence.

    Native artifacts remain immutable sources. Normalized events are additive primary-view
    records linked back to exactly one locus-bearing primary-resolution caller segment.
    """
    artifacts, full_evidence = map_qdnaseq_ace_full_evidence(
        report=report,
        policy=policy,
        output_dir=output_dir,
        registration_sha256=registration_sha256,
        specimen=specimen,
        caller_lock=caller_lock,
        replicate_id=replicate_id,
    )
    normalized_events = _normalized_event_records(
        report=report,
        policy=policy,
        registration_sha256=registration_sha256,
        specimen=specimen,
        caller_lock=caller_lock,
        replicate_id=replicate_id,
        full_evidence=full_evidence,
    )
    return seal_cnv_validation_evidence(
        manifest_id=manifest_id,
        registration_sha256=registration_sha256,
        native_artifacts=artifacts,
        full_evidence=full_evidence,
        normalized_events=normalized_events,
    )
