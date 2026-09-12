"""Paired-source methylation mixtures and source-fraction quantification.

This module is deliberately separate from :mod:`ontseq_platform.dilution`.  That module
characterises CNV/SV detection after mixing BAMs labelled tumour and normal.  Here the
labels are neutral ``source_a`` and ``source_b`` because methylation alone does not turn a
read fraction into tumour purity, blast fraction or cell fraction.

The first input adapter reads Nanopolish ``call-methylation`` TSV files.  Calls are kept
together by read name, split into calibration and held-out mixture pools, and sampled as
whole read groups.  Read names never enter the report; SHA-256 selection digests prove
that a seeded draw is reproducible without exporting identifiers.

For every exact marker shared by the two calibration pools, methylated and
unmethylated call rates per read group follow

    rate_mix = f * rate_a + (1 - f) * rate_b

and ``f`` is estimated by constrained weighted least squares. A conditional Dirichlet Monte
Carlo interval propagates mixture-pool sampling uncertainty around the fitted four-state model
while the calibration marker lock remains fixed. It is not a clinical confidence interval and
does not account for calibration-source or biological population variation.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import html
import math
import random
import re
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import BinaryIO, Literal, Protocol, cast

from pydantic import Field, model_validator

from . import __version__
from .models import FileFingerprint, GenomeBuild, ModuleRunStatus, StrictModel
from .reference import sha256_file

_IDENTIFIER = r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$"
_ANALYSIS_IDENTIFIER = r"^[A-Za-z0-9][A-Za-z0-9._-]{2,69}$"
_PROVENANCE_VALUE = r"^[A-Za-z0-9][A-Za-z0-9 ._:+()\-]{0,199}$"
_CANONICAL_CHROMOSOME = r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$"
NANOPOLISH_ADAPTER_VERSION = "nanopolish-call-table-v2"
# Versioned intake limits; independent of assay/estimator thresholds.
_NANOPOLISH_MAXIMUM_PHYSICAL_LINE_BYTES = 65_536
_NANOPOLISH_MAXIMUM_HEADER_BYTES = 16_384
_NANOPOLISH_MAXIMUM_COLUMNS = 64
_READ_NAME_CONTROL = re.compile(r"[\x00-\x1f\x7f-\x9f]")
_REQUIRED_NANOPOLISH_COLUMNS = frozenset(
    {
        "chromosome",
        "strand",
        "start",
        "end",
        "read_name",
        "log_lik_ratio",
        "log_lik_methylated",
        "log_lik_unmethylated",
        "num_calling_strands",
        "num_motifs",
        "sequence",
    }
)


class NanopolishCallState(StrEnum):
    METHYLATED = "methylated"
    UNMETHYLATED = "unmethylated"
    AMBIGUOUS = "ambiguous"


class MethylationMixtureReportStatus(StrEnum):
    """Aggregate state of a complete requested dilution grid."""

    COMPLETED = "COMPLETED"
    PARTIAL = "PARTIAL"
    NO_CALL = "NO_CALL"


class TechnicalRecoveryAssessment(StrEnum):
    """Technical agreement with known in-silico fractions, never a clinical verdict."""

    PASS = "PASS"
    FAIL = "FAIL"
    NOT_EVALUABLE = "NOT_EVALUABLE"


class NanopolishSourceMetadata(StrictModel):
    """Declared upstream protocol metadata for one call table.

    ``unknown`` is explicit rather than inferred from a filename.  Critical fields that are
    known for both sources are compatibility gates in the paired-source experiment.
    """

    schema_version: Literal["0.1.0"] = "0.1.0"
    caller_name: Literal["nanopolish", "unknown"] = "unknown"
    caller_version: str = Field(default="unknown", pattern=_PROVENANCE_VALUE)
    reference_id: str = Field(default="unknown", pattern=_PROVENANCE_VALUE)
    reference_genome_build: GenomeBuild | Literal["unknown"] = "unknown"
    reference_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    sequencing_platform: str = Field(default="unknown", pattern=_PROVENANCE_VALUE)
    flow_cell_product_code: str = Field(default="unknown", pattern=_PROVENANCE_VALUE)
    library_kit: str = Field(default="unknown", pattern=_PROVENANCE_VALUE)
    basecaller_name: str = Field(default="unknown", pattern=_PROVENANCE_VALUE)
    basecaller_version: str = Field(default="unknown", pattern=_PROVENANCE_VALUE)
    basecaller_model: str = Field(default="unknown", pattern=_PROVENANCE_VALUE)
    modification_context: Literal["CpG_5mC_vs_unmethylated", "unknown"] = "unknown"

    def unknown_fields(self) -> tuple[str, ...]:
        values = self.model_dump()
        return tuple(
            name
            for name, value in values.items()
            if name != "schema_version" and (value == "unknown" or value is None)
        )


_COMPATIBILITY_METADATA_FIELDS: tuple[str, ...] = (
    "caller_name",
    "caller_version",
    "reference_id",
    "reference_genome_build",
    "reference_sha256",
    "sequencing_platform",
    "flow_cell_product_code",
    "library_kit",
    "basecaller_name",
    "basecaller_version",
    "basecaller_model",
    "modification_context",
)


def _known_metadata_mismatches(
    source_a: NanopolishSourceMetadata, source_b: NanopolishSourceMetadata
) -> tuple[str, ...]:
    mismatches: list[str] = []
    for field_name in _COMPATIBILITY_METADATA_FIELDS:
        value_a = getattr(source_a, field_name)
        value_b = getattr(source_b, field_name)
        if value_a in {None, "unknown"} or value_b in {None, "unknown"}:
            continue
        if value_a != value_b:
            mismatches.append(field_name)
    return tuple(mismatches)


def _meets_minimum_absolute_delta_beta(value: float, minimum: float) -> bool:
    magnitude = abs(value)
    return magnitude > minimum or math.isclose(magnitude, minimum, rel_tol=0.0, abs_tol=1e-12)


def _calibration_read_group_count(read_group_count: int, fraction: float) -> int:
    reference_count = round(read_group_count * fraction)
    return min(max(reference_count, 1), read_group_count - 1)


def _mixture_level_seed(policy_seed: int, fraction_index: int, replicate: int) -> int:
    return policy_seed + 10_000_019 + fraction_index * 10_007 + replicate * 101


class MethylationMixturePolicy(StrictModel):
    """Versioned technical choices for one paired-source experiment."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only"] = "technical_defaults_only"
    log_likelihood_ratio_threshold: float = Field(default=2.5, gt=0)
    source_a_fractions: list[float] = Field(
        default_factory=lambda: [0.0, 0.05, 0.1, 0.2, 0.5, 1.0], min_length=2
    )
    replicates: int = Field(default=3, ge=1)
    seed: int = Field(default=20260904, ge=0, le=2_147_483_647)
    calibration_read_group_fraction: float = Field(default=0.25, gt=0, lt=1)
    total_mixture_read_groups: int | None = Field(default=None, ge=2)
    minimum_reference_valid_calls: int = Field(default=5, ge=1)
    minimum_mixture_valid_calls: int = Field(default=3, ge=1)
    minimum_absolute_delta_beta: float = Field(default=0.2, gt=0, le=1)
    minimum_markers: int = Field(default=20, ge=1)
    maximum_markers: int | None = Field(default=2000, ge=1)
    uncertainty_draws: int = Field(default=500, ge=100)
    confidence_level: float = Field(default=0.95, gt=0.5, lt=1)
    minimum_successful_uncertainty_fraction: float = Field(default=0.8, gt=0, le=1)
    maximum_confidence_interval_width: float = Field(default=0.5, gt=0, le=1)
    maximum_standardized_model_fit_rmse: float = Field(default=3.0, gt=0)
    maximum_extrapolation: float = Field(default=0.1, ge=0, le=1)
    maximum_input_bytes_per_source: int = Field(
        default=100_000_000,
        ge=1,
        description=(
            "Nanopolish adapter v2 applies this limit independently to stored bytes and "
            "decoded bytes, including headers and blank lines."
        ),
    )
    maximum_rows_per_source: int = Field(default=2_000_000, ge=1)
    minimum_quantifiable_fraction_for_recovery: float = Field(default=0.8, gt=0, le=1)
    maximum_absolute_mean_bias: float = Field(default=0.1, ge=0, le=1)
    maximum_mean_absolute_error: float = Field(default=0.15, ge=0, le=1)
    maximum_root_mean_squared_error: float = Field(default=0.2, ge=0, le=1)
    minimum_confidence_interval_coverage_fraction: float = Field(default=0.8, ge=0, le=1)
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def policy_is_coherent(self) -> MethylationMixturePolicy:
        fractions = self.source_a_fractions
        if fractions != sorted(set(fractions)):
            raise ValueError("source_a_fractions must be unique and ordered from low to high")
        if any(value < 0 or value > 1 for value in fractions):
            raise ValueError("source_a_fractions must lie between zero and one")
        if self.maximum_markers is not None and self.maximum_markers < self.minimum_markers:
            raise ValueError("maximum_markers cannot be smaller than minimum_markers")
        return self


class MethylationSourceSummary(StrictModel):
    source_id: str = Field(pattern=_IDENTIFIER)
    genome_build: GenomeBuild
    adapter_version: Literal["nanopolish-call-table-v2"] = "nanopolish-call-table-v2"
    input_security_profile: Literal["bounded-tsv-gzip-v1"] = "bounded-tsv-gzip-v1"
    decoded_size_bytes: int = Field(ge=1)
    input_format: Literal["nanopolish_call_methylation_tsv"] = "nanopolish_call_methylation_tsv"
    coordinate_system: Literal["nanopolish_start_end_as_reported"] = (
        "nanopolish_start_end_as_reported"
    )
    chromosome_normalization: Literal["canonical_chr_prefix"] = "canonical_chr_prefix"
    upstream: NanopolishSourceMetadata = Field(default_factory=NanopolishSourceMetadata)
    fingerprint: FileFingerprint
    rows_total: int = Field(ge=0)
    canonical_call_rows: int = Field(ge=0)
    informative_call_rows: int = Field(ge=0)
    ambiguous_call_rows: int = Field(ge=0)
    skipped_noncanonical_rows: int = Field(ge=0)
    read_group_count: int = Field(ge=1)

    @model_validator(mode="after")
    def row_counts_are_consistent(self) -> MethylationSourceSummary:
        if self.rows_total != self.canonical_call_rows + self.skipped_noncanonical_rows:
            raise ValueError("Nanopolish source row counts are inconsistent")
        if self.canonical_call_rows != self.informative_call_rows + self.ambiguous_call_rows:
            raise ValueError("Nanopolish informative and ambiguous counts are inconsistent")
        if self.fingerprint.sha256 is None:
            raise ValueError("A methylation source must carry a SHA-256 fingerprint")
        if self.upstream.reference_genome_build not in {"unknown", self.genome_build}:
            raise ValueError("Declared reference genome build contradicts the source genome build")
        return self


class MethylationReferenceMarker(StrictModel):
    marker_id: str = Field(min_length=1)
    chromosome: str = Field(pattern=r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    strand: Literal["+", "-"]
    start: int = Field(ge=0)
    end: int = Field(ge=0)
    motif_count: int = Field(ge=1)
    sequence_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    source_a_modified_calls: int = Field(ge=0)
    source_a_unmethylated_calls: int = Field(ge=0)
    source_a_ambiguous_calls: int = Field(ge=0)
    source_a_valid_calls: int = Field(ge=1)
    source_a_total_read_groups: int = Field(ge=1)
    source_b_modified_calls: int = Field(ge=0)
    source_b_unmethylated_calls: int = Field(ge=0)
    source_b_ambiguous_calls: int = Field(ge=0)
    source_b_valid_calls: int = Field(ge=1)
    source_b_total_read_groups: int = Field(ge=1)
    source_a_beta: float = Field(ge=0, le=1)
    source_b_beta: float = Field(ge=0, le=1)
    delta_beta: float = Field(ge=-1, le=1)

    @model_validator(mode="after")
    def marker_is_consistent(self) -> MethylationReferenceMarker:
        expected_id = (
            f"{self.chromosome}:{self.start}-{self.end}:{self.strand}:motifs={self.motif_count}:"
            f"sequence_sha256={self.sequence_sha256}"
        )
        if self.marker_id != expected_id:
            raise ValueError("Methylation marker ID is inconsistent with its exact coordinates")
        if self.end < self.start:
            raise ValueError("Methylation marker end cannot be smaller than start")
        if self.source_a_modified_calls > self.source_a_valid_calls:
            raise ValueError("Source A has more modified than valid calls")
        if self.source_b_modified_calls > self.source_b_valid_calls:
            raise ValueError("Source B has more modified than valid calls")
        if (
            self.source_a_modified_calls + self.source_a_unmethylated_calls
            != self.source_a_valid_calls
        ):
            raise ValueError("Source A valid calls do not equal its classified calls")
        if (
            self.source_b_modified_calls + self.source_b_unmethylated_calls
            != self.source_b_valid_calls
        ):
            raise ValueError("Source B valid calls do not equal its classified calls")
        if (
            self.source_a_valid_calls + self.source_a_ambiguous_calls
            > self.source_a_total_read_groups
        ):
            raise ValueError("Source A marker calls exceed its calibration read groups")
        if (
            self.source_b_valid_calls + self.source_b_ambiguous_calls
            > self.source_b_total_read_groups
        ):
            raise ValueError("Source B marker calls exceed its calibration read groups")
        if not math.isclose(
            self.source_a_beta,
            self.source_a_modified_calls / self.source_a_valid_calls,
            abs_tol=1e-12,
        ):
            raise ValueError("Source A beta is inconsistent with its counts")
        if not math.isclose(
            self.source_b_beta,
            self.source_b_modified_calls / self.source_b_valid_calls,
            abs_tol=1e-12,
        ):
            raise ValueError("Source B beta is inconsistent with its counts")
        if not math.isclose(
            self.delta_beta, self.source_a_beta - self.source_b_beta, abs_tol=1e-12
        ):
            raise ValueError("Marker delta_beta is inconsistent")
        return self


class MethylationMixtureMarkerObservation(StrictModel):
    marker_id: str = Field(min_length=1)
    modified_calls: int = Field(ge=0)
    unmethylated_calls: int = Field(ge=0)
    valid_calls: int = Field(ge=0)
    ambiguous_calls: int = Field(ge=0)
    beta: float | None = Field(default=None, ge=0, le=1)

    @model_validator(mode="after")
    def observation_is_consistent(self) -> MethylationMixtureMarkerObservation:
        if self.modified_calls > self.valid_calls:
            raise ValueError("Mixture has more modified than valid calls")
        if self.modified_calls + self.unmethylated_calls != self.valid_calls:
            raise ValueError("Mixture valid calls do not equal its classified calls")
        if (self.valid_calls > 0) != (self.beta is not None):
            raise ValueError("Mixture beta must exist exactly when valid calls exist")
        if self.beta is not None and not math.isclose(
            self.beta, self.modified_calls / self.valid_calls, abs_tol=1e-12
        ):
            raise ValueError("Mixture beta is inconsistent with its counts")
        return self


class MethylationMixtureLevelResult(StrictModel):
    level_id: str = Field(pattern=_IDENTIFIER)
    replicate: int = Field(ge=1)
    seed: int = Field(ge=0)
    target_source_a_fraction: float = Field(ge=0, le=1)
    source_a_read_groups: int = Field(ge=0)
    source_b_read_groups: int = Field(ge=0)
    total_read_groups: int = Field(ge=1)
    realized_source_a_read_group_fraction: float = Field(ge=0, le=1)
    selection_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    status: ModuleRunStatus
    estimated_source_a_fraction: float | None = Field(default=None, ge=0, le=1)
    unconstrained_source_a_fraction: float | None = None
    confidence_interval_lower: float | None = Field(default=None, ge=0, le=1)
    confidence_interval_upper: float | None = Field(default=None, ge=0, le=1)
    confidence_level: float = Field(gt=0.5, lt=1)
    markers_available: int = Field(ge=0)
    markers_used: int = Field(ge=0)
    valid_mixture_calls: int = Field(ge=0)
    model_fit_standardized_rmse: float | None = Field(default=None, ge=0)
    successful_uncertainty_draws: int = Field(ge=0)
    no_call_reason: str | None = None
    warnings: list[str] = Field(default_factory=list)
    marker_observations: list[MethylationMixtureMarkerObservation] = Field(default_factory=list)

    @model_validator(mode="after")
    def level_is_consistent(self) -> MethylationMixtureLevelResult:
        if self.source_a_read_groups + self.source_b_read_groups != self.total_read_groups:
            raise ValueError("Mixture read-group counts do not add up")
        expected = self.source_a_read_groups / self.total_read_groups
        if not math.isclose(self.realized_source_a_read_group_fraction, expected, abs_tol=1e-12):
            raise ValueError("Realized source fraction is inconsistent with read-group counts")
        if self.markers_used > self.markers_available:
            raise ValueError("A mixture cannot use more markers than are available")
        if any(
            item.valid_calls + item.ambiguous_calls > self.total_read_groups
            for item in self.marker_observations
        ):
            raise ValueError("A marker has more calls than mixture read groups")
        if self.status not in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL}:
            raise ValueError("A mixture estimate must be COMPLETED or NO_CALL")
        estimated = self.estimated_source_a_fraction is not None
        bounded = (
            self.confidence_interval_lower is not None
            and self.confidence_interval_upper is not None
        )
        if (self.confidence_interval_lower is None) != (self.confidence_interval_upper is None):
            raise ValueError("Confidence interval bounds must be both present or both absent")
        if self.status == ModuleRunStatus.COMPLETED:
            if (
                not estimated
                or self.unconstrained_source_a_fraction is None
                or not bounded
                or self.model_fit_standardized_rmse is None
            ):
                raise ValueError(
                    "A completed estimate needs constrained and raw points, interval and fit"
                )
            if self.no_call_reason is not None:
                raise ValueError("A completed estimate cannot carry a no-call reason")
            assert self.confidence_interval_lower is not None
            assert self.confidence_interval_upper is not None
            if self.confidence_interval_upper < self.confidence_interval_lower:
                raise ValueError("Confidence interval bounds are reversed")
        else:
            if (
                estimated
                or self.unconstrained_source_a_fraction is not None
                or bounded
                or self.model_fit_standardized_rmse is not None
            ):
                raise ValueError("A no-call must not expose a quantitative estimate")
            if not self.no_call_reason:
                raise ValueError("A no-call must state its reason")
        return self


def _technical_recovery_assessment(
    *,
    policy: MethylationMixturePolicy,
    quantifiable_fraction: float,
    mean_bias: float | None,
    mean_absolute_error: float | None,
    root_mean_squared_error: float | None,
    confidence_interval_coverage_fraction: float | None,
) -> tuple[TechnicalRecoveryAssessment, list[str]]:
    if (
        mean_bias is None
        or mean_absolute_error is None
        or root_mean_squared_error is None
        or confidence_interval_coverage_fraction is None
        or quantifiable_fraction < policy.minimum_quantifiable_fraction_for_recovery
    ):
        return (
            TechnicalRecoveryAssessment.NOT_EVALUABLE,
            [
                "The quantifiable share of the requested grid is below the locked minimum "
                "for technical recovery assessment."
            ],
        )

    reasons: list[str] = []
    if abs(mean_bias) > policy.maximum_absolute_mean_bias:
        reasons.append(
            f"Absolute mean bias {abs(mean_bias):.4f} exceeds "
            f"{policy.maximum_absolute_mean_bias:.4f}."
        )
    if mean_absolute_error > policy.maximum_mean_absolute_error:
        reasons.append(
            f"Mean absolute error {mean_absolute_error:.4f} exceeds "
            f"{policy.maximum_mean_absolute_error:.4f}."
        )
    if root_mean_squared_error > policy.maximum_root_mean_squared_error:
        reasons.append(
            f"Root mean squared error {root_mean_squared_error:.4f} exceeds "
            f"{policy.maximum_root_mean_squared_error:.4f}."
        )
    if confidence_interval_coverage_fraction < policy.minimum_confidence_interval_coverage_fraction:
        reasons.append(
            "Empirical confidence-interval coverage "
            f"{confidence_interval_coverage_fraction:.4f} is below "
            f"{policy.minimum_confidence_interval_coverage_fraction:.4f}."
        )
    return (
        TechnicalRecoveryAssessment.FAIL if reasons else TechnicalRecoveryAssessment.PASS,
        reasons,
    )


class MethylationMixtureReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    analysis_id: str = Field(pattern=_ANALYSIS_IDENTIFIER)
    genome_build: GenomeBuild
    source_relationship: Literal["declared_biologically_distinct"] = (
        "declared_biologically_distinct"
    )
    source_fraction_definition: Literal["fraction_of_held_out_read_groups_from_source_a"] = (
        "fraction_of_held_out_read_groups_from_source_a"
    )
    estimand: Literal["methylation_signal_mixture_coefficient_for_source_a"] = (
        "methylation_signal_mixture_coefficient_for_source_a"
    )
    known_fraction_basis: Literal["nanopolish_read_group_count"] = "nanopolish_read_group_count"
    algorithm_version: Literal["paired-source-call-rate-wls-conditional-dirichlet-v2"] = (
        "paired-source-call-rate-wls-conditional-dirichlet-v2"
    )
    status: MethylationMixtureReportStatus
    software_version: str = Field(pattern=_PROVENANCE_VALUE)
    git_commit: str = Field(pattern=r"^(?:UNKNOWN|[0-9a-f]{40})$")
    policy: MethylationMixturePolicy
    source_a: MethylationSourceSummary
    source_b: MethylationSourceSummary
    source_a_calibration_read_groups: int = Field(ge=1)
    source_b_calibration_read_groups: int = Field(ge=1)
    source_a_mixture_pool_read_groups: int = Field(ge=1)
    source_b_mixture_pool_read_groups: int = Field(ge=1)
    calibration_selection_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    reference_marker_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    reference_markers: list[MethylationReferenceMarker] = Field(default_factory=list)
    levels: list[MethylationMixtureLevelResult] = Field(min_length=1)
    completed_levels: int = Field(ge=0)
    no_call_levels: int = Field(ge=0)
    quantifiable_fraction: float = Field(ge=0, le=1)
    grid_fully_quantified: bool
    accuracy_metrics_scope: Literal["completed_levels_only"] = "completed_levels_only"
    mean_bias: float | None = None
    mean_absolute_error: float | None = Field(default=None, ge=0)
    root_mean_squared_error: float | None = Field(default=None, ge=0)
    confidence_interval_coverage_fraction: float | None = Field(default=None, ge=0, le=1)
    technical_recovery_assessment: TechnicalRecoveryAssessment
    technical_recovery_reasons: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    sensitive_output: Literal[True] = True
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def report_is_consistent(self) -> MethylationMixtureReport:
        if self.source_a.source_id == self.source_b.source_id:
            raise ValueError("Paired sources must have different source IDs")
        if self.source_a.fingerprint.sha256 == self.source_b.fingerprint.sha256:
            raise ValueError("Paired sources must not have byte-identical fingerprints")
        if self.source_a.genome_build != self.source_b.genome_build:
            raise ValueError("Paired sources must use the same genome build")
        if self.genome_build != self.source_a.genome_build:
            raise ValueError("Report genome build does not match its sources")
        metadata_mismatches = _known_metadata_mismatches(
            self.source_a.upstream, self.source_b.upstream
        )
        if metadata_mismatches:
            raise ValueError(
                "Paired sources have incompatible known upstream metadata: "
                + ", ".join(metadata_mismatches)
            )
        for source in (self.source_a, self.source_b):
            if source.fingerprint.size_bytes > self.policy.maximum_input_bytes_per_source:
                raise ValueError("Source exceeds the locked input-byte safety limit")
            if source.decoded_size_bytes > self.policy.maximum_input_bytes_per_source:
                raise ValueError("Source exceeds the locked decoded-byte safety limit")
            if source.rows_total > self.policy.maximum_rows_per_source:
                raise ValueError("Source exceeds the locked input-row safety limit")
        if (
            self.source_a_calibration_read_groups + self.source_a_mixture_pool_read_groups
            != self.source_a.read_group_count
            or self.source_b_calibration_read_groups + self.source_b_mixture_pool_read_groups
            != self.source_b.read_group_count
        ):
            raise ValueError("Calibration and mixture pools do not partition both sources")
        expected_calibration_counts = (
            _calibration_read_group_count(
                self.source_a.read_group_count,
                self.policy.calibration_read_group_fraction,
            ),
            _calibration_read_group_count(
                self.source_b.read_group_count,
                self.policy.calibration_read_group_fraction,
            ),
        )
        if (
            self.source_a_calibration_read_groups,
            self.source_b_calibration_read_groups,
        ) != expected_calibration_counts:
            raise ValueError("Calibration pool sizes do not match the locked policy")
        marker_ids = [marker.marker_id for marker in self.reference_markers]
        if marker_ids != sorted(set(marker_ids)):
            raise ValueError("Reference marker IDs must be unique and canonically ordered")
        if any(
            marker.source_a_total_read_groups != self.source_a_calibration_read_groups
            or marker.source_b_total_read_groups != self.source_b_calibration_read_groups
            for marker in self.reference_markers
        ):
            raise ValueError("Reference marker totals do not match the calibration pools")
        if any(
            marker.source_a_valid_calls < self.policy.minimum_reference_valid_calls
            or marker.source_b_valid_calls < self.policy.minimum_reference_valid_calls
            or not _meets_minimum_absolute_delta_beta(
                marker.delta_beta, self.policy.minimum_absolute_delta_beta
            )
            for marker in self.reference_markers
        ):
            raise ValueError("Reference marker lock violates the marker-selection policy")
        if (
            self.policy.maximum_markers is not None
            and len(self.reference_markers) > self.policy.maximum_markers
        ):
            raise ValueError("Reference marker lock exceeds the policy maximum")
        if self.reference_marker_sha256 != _reference_marker_digest(self.reference_markers):
            raise ValueError("Reference marker lock checksum is inconsistent")

        expected_grid = [
            (fraction, replicate, f"{self.analysis_id}.L{index + 1:03d}.R{replicate:03d}")
            for index, fraction in enumerate(self.policy.source_a_fractions)
            for replicate in range(1, self.policy.replicates + 1)
        ]
        actual_grid = [
            (level.target_source_a_fraction, level.replicate, level.level_id)
            for level in self.levels
        ]
        if actual_grid != expected_grid:
            raise ValueError("Mixture levels do not match the complete policy fraction grid")
        budgets = {level.total_read_groups for level in self.levels}
        if len(budgets) != 1:
            raise ValueError("Every mixture level must use the same read-group budget")
        if self.policy.total_mixture_read_groups is not None and budgets != {
            self.policy.total_mixture_read_groups
        }:
            raise ValueError("Mixture budget does not match the pinned policy value")
        for level_index, level in enumerate(self.levels):
            fraction_index = level_index // self.policy.replicates
            expected_seed = _mixture_level_seed(self.policy.seed, fraction_index, level.replicate)
            if level.seed != expected_seed:
                raise ValueError("A mixture level seed does not match the locked policy")
            if level.source_a_read_groups > self.source_a_mixture_pool_read_groups or (
                level.source_b_read_groups > self.source_b_mixture_pool_read_groups
            ):
                raise ValueError("A mixture level exceeds an available held-out source pool")
            if level.source_a_read_groups != round(
                level.total_read_groups * level.target_source_a_fraction
            ):
                raise ValueError("A mixture level does not implement its target fraction")
            if not math.isclose(
                level.confidence_level, self.policy.confidence_level, abs_tol=1e-12
            ):
                raise ValueError("A level confidence value differs from the locked policy")
            if level.markers_available != len(self.reference_markers):
                raise ValueError("A level marker count differs from the reference lock")
            observations = level.marker_observations
            observation_ids = [item.marker_id for item in observations]
            if observation_ids != marker_ids:
                raise ValueError("Level observations must exactly match the ordered marker lock")
            expected_used = sum(
                item.valid_calls >= self.policy.minimum_mixture_valid_calls for item in observations
            )
            if level.markers_used != expected_used:
                raise ValueError("Level markers_used is inconsistent with its observations")
            if level.valid_mixture_calls != sum(item.valid_calls for item in observations):
                raise ValueError("Level valid-call total is inconsistent with its observations")
            if level.successful_uncertainty_draws > self.policy.uncertainty_draws:
                raise ValueError("Level uncertainty draws exceed the locked policy")
            observation_counts = {
                item.marker_id: _Counts(
                    modified=item.modified_calls,
                    unmethylated=item.unmethylated_calls,
                    ambiguous=item.ambiguous_calls,
                )
                for item in observations
            }
            expected_raw, expected_point, expected_fit, _expected_observations = _point_estimate(
                self.reference_markers,
                observation_counts,
                minimum_mixture_valid_calls=self.policy.minimum_mixture_valid_calls,
                mixture_total_read_groups=level.total_read_groups,
            )
            if level.status == ModuleRunStatus.COMPLETED:
                assert level.estimated_source_a_fraction is not None
                assert level.unconstrained_source_a_fraction is not None
                assert level.confidence_interval_lower is not None
                assert level.confidence_interval_upper is not None
                assert level.model_fit_standardized_rmse is not None
                if expected_raw is None or expected_point is None or expected_fit is None:
                    raise ValueError("Completed level cannot be recomputed from its observations")
                if not math.isclose(
                    level.unconstrained_source_a_fraction,
                    expected_raw,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ):
                    raise ValueError("Completed raw estimate is inconsistent with marker counts")
                if not math.isclose(
                    level.estimated_source_a_fraction,
                    expected_point,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ):
                    raise ValueError("Completed estimate is inconsistent with marker counts")
                if not math.isclose(
                    level.model_fit_standardized_rmse,
                    expected_fit,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ):
                    raise ValueError("Completed fit is inconsistent with marker counts")
                if level.markers_used < self.policy.minimum_markers:
                    raise ValueError("Completed level violates the minimum-marker policy")
                if (
                    level.unconstrained_source_a_fraction < -self.policy.maximum_extrapolation
                    or level.unconstrained_source_a_fraction > 1 + self.policy.maximum_extrapolation
                ):
                    raise ValueError("Completed level violates the extrapolation policy")
                if (
                    level.model_fit_standardized_rmse
                    > self.policy.maximum_standardized_model_fit_rmse
                ):
                    raise ValueError("Completed level violates the model-fit policy")
                if (
                    level.confidence_interval_upper - level.confidence_interval_lower
                    > self.policy.maximum_confidence_interval_width
                ):
                    raise ValueError("Completed level violates the interval-width policy")
                required_draws = math.ceil(
                    self.policy.uncertainty_draws
                    * self.policy.minimum_successful_uncertainty_fraction
                )
                if level.successful_uncertainty_draws < required_draws:
                    raise ValueError("Completed level violates the uncertainty-draw policy")
                rng = random.Random(level.seed + 9_000_001)
                recomputed_draws = [
                    _draw_estimate(
                        self.reference_markers,
                        observation_counts,
                        minimum_mixture_valid_calls=self.policy.minimum_mixture_valid_calls,
                        mixture_total_read_groups=level.total_read_groups,
                        fitted_source_a_fraction=expected_point,
                        rng=rng,
                    )
                    for _ in range(self.policy.uncertainty_draws)
                ]
                successful = sorted(value for value in recomputed_draws if value is not None)
                if level.successful_uncertainty_draws != len(successful):
                    raise ValueError("Completed uncertainty-draw count is inconsistent")
                tail = (1 - self.policy.confidence_level) / 2
                expected_lower = _quantile(successful, tail)
                expected_upper = _quantile(successful, 1 - tail)
                if not math.isclose(
                    level.confidence_interval_lower,
                    expected_lower,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ) or not math.isclose(
                    level.confidence_interval_upper,
                    expected_upper,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ):
                    raise ValueError(
                        "Completed interval is inconsistent with marker counts and seed"
                    )
        completed = sum(level.status == ModuleRunStatus.COMPLETED for level in self.levels)
        no_call = sum(level.status == ModuleRunStatus.NO_CALL for level in self.levels)
        if completed != self.completed_levels or no_call != self.no_call_levels:
            raise ValueError("Mixture report level counts are inconsistent")
        if completed == len(self.levels):
            expected_status = MethylationMixtureReportStatus.COMPLETED
        elif completed:
            expected_status = MethylationMixtureReportStatus.PARTIAL
        else:
            expected_status = MethylationMixtureReportStatus.NO_CALL
        if self.status != expected_status:
            raise ValueError("Mixture report status is inconsistent with its levels")
        expected_quantifiable = completed / len(self.levels)
        if not math.isclose(self.quantifiable_fraction, expected_quantifiable, abs_tol=1e-12):
            raise ValueError("Report quantifiable_fraction is inconsistent with its levels")
        if self.grid_fully_quantified != (completed == len(self.levels)):
            raise ValueError("Report complete-grid indicator is inconsistent with its levels")
        metrics = (self.mean_bias, self.mean_absolute_error, self.root_mean_squared_error)
        metrics_present = all(value is not None for value in metrics)
        if any(value is not None for value in metrics) and not metrics_present:
            raise ValueError("Accuracy metrics must be all present or all absent")
        if metrics_present != (completed > 0):
            raise ValueError("Accuracy metrics must exist exactly when an estimate completed")
        if (self.confidence_interval_coverage_fraction is not None) != (completed > 0):
            raise ValueError("Interval coverage must exist exactly when an estimate completed")
        if completed:
            completed_levels = [
                level for level in self.levels if level.status == ModuleRunStatus.COMPLETED
            ]
            errors = [
                cast(float, level.estimated_source_a_fraction)
                - level.realized_source_a_read_group_fraction
                for level in completed_levels
            ]
            expected_metrics = (
                sum(errors) / len(errors),
                sum(abs(value) for value in errors) / len(errors),
                math.sqrt(sum(value**2 for value in errors) / len(errors)),
            )
            if not all(
                math.isclose(cast(float, actual), expected, abs_tol=1e-12)
                for actual, expected in zip(metrics, expected_metrics, strict=True)
            ):
                raise ValueError("Accuracy metrics are inconsistent with completed levels")
            interval_hits = sum(
                cast(float, level.confidence_interval_lower)
                <= level.realized_source_a_read_group_fraction
                <= cast(float, level.confidence_interval_upper)
                for level in completed_levels
            )
            expected_coverage = interval_hits / len(completed_levels)
            if not math.isclose(
                cast(float, self.confidence_interval_coverage_fraction),
                expected_coverage,
                abs_tol=1e-12,
            ):
                raise ValueError("Confidence-interval coverage is inconsistent with levels")
        expected_assessment, expected_reasons = _technical_recovery_assessment(
            policy=self.policy,
            quantifiable_fraction=self.quantifiable_fraction,
            mean_bias=self.mean_bias,
            mean_absolute_error=self.mean_absolute_error,
            root_mean_squared_error=self.root_mean_squared_error,
            confidence_interval_coverage_fraction=self.confidence_interval_coverage_fraction,
        )
        if (
            self.technical_recovery_assessment != expected_assessment
            or self.technical_recovery_reasons != expected_reasons
        ):
            raise ValueError("Technical recovery assessment is inconsistent with report metrics")
        if self.limitations != list(_LIMITATIONS):
            raise ValueError("Report limitations differ from the locked software limitations")
        expected_warnings = _report_warnings(
            policy=self.policy,
            source_a=self.source_a,
            source_b=self.source_b,
            total_mixture_read_groups=next(iter(budgets)),
            levels=self.levels,
        )
        if self.warnings != expected_warnings:
            raise ValueError("Report warnings are inconsistent with its inputs and level results")
        return self


@dataclass(frozen=True, order=True)
class _MarkerKey:
    chromosome: str
    strand: Literal["+", "-"]
    start: int
    end: int
    motif_count: int
    sequence_sha256: str

    @property
    def marker_id(self) -> str:
        return (
            f"{self.chromosome}:{self.start}-{self.end}:{self.strand}:motifs={self.motif_count}:"
            f"sequence_sha256={self.sequence_sha256}"
        )


@dataclass(frozen=True)
class _NanopolishCall:
    marker: _MarkerKey
    state: NanopolishCallState


@dataclass
class _Counts:
    modified: int = 0
    unmethylated: int = 0
    ambiguous: int = 0

    @property
    def valid(self) -> int:
        return self.modified + self.unmethylated

    @property
    def beta(self) -> float | None:
        return self.modified / self.valid if self.valid else None

    def add(self, state: NanopolishCallState) -> None:
        if state == NanopolishCallState.METHYLATED:
            self.modified += 1
        elif state == NanopolishCallState.UNMETHYLATED:
            self.unmethylated += 1
        else:
            self.ambiguous += 1


@dataclass(frozen=True)
class _ParsedSource:
    summary: MethylationSourceSummary
    calls_by_read: Mapping[str, tuple[_NanopolishCall, ...]]


class ReadCallSource(Protocol):
    """Adapter-neutral intact read groups accepted by the call-rate estimator."""

    @property
    def calls_by_read(self) -> Mapping[str, tuple[_NanopolishCall, ...]]: ...


@dataclass(frozen=True)
class _LevelContext:
    level_id: str
    replicate: int
    seed: int
    target_fraction: float
    source_a_reads: Sequence[str]
    source_b_reads: Sequence[str]
    selection_sha256: str
    confidence_level: float
    markers_available: int


@contextmanager
def _open_nanopolish(path: Path) -> Iterator[BinaryIO | gzip.GzipFile]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rb") as handle:
            yield handle
    else:
        with path.open("rb") as handle:
            yield handle


@dataclass
class _BoundedNanopolishLines(Iterator[str]):
    """Bound decoded bytes before UTF-8/CSV allocation, including skipped blank lines."""

    handle: BinaryIO | gzip.GzipFile
    maximum_bytes: int
    decoded_size_bytes: int = 0
    physical_lines: int = 0

    def __next__(self) -> str:
        header = self.physical_lines == 0
        line_limit = (
            _NANOPOLISH_MAXIMUM_HEADER_BYTES if header else _NANOPOLISH_MAXIMUM_PHYSICAL_LINE_BYTES
        )
        remaining = self.maximum_bytes - self.decoded_size_bytes
        raw = self.handle.readline(min(line_limit, remaining) + 1)
        if not raw:
            raise StopIteration
        self.decoded_size_bytes += len(raw)
        if self.decoded_size_bytes > self.maximum_bytes:
            raise ValueError("Nanopolish source exceeds the locked decoded-byte safety limit")
        if len(raw) > line_limit:
            kind = "header" if header else "physical line"
            raise ValueError(f"Nanopolish {kind} exceeds the {line_limit}-byte safety limit")
        if raw.count(b"\t") >= _NANOPOLISH_MAXIMUM_COLUMNS:
            raise ValueError("Nanopolish physical line exceeds the column safety limit")
        self.physical_lines += 1
        return raw.decode("utf-8-sig" if header else "utf-8")


def _reject_read_name_controls(read_name: str) -> None:
    if _READ_NAME_CONTROL.search(read_name):
        raise ValueError("Nanopolish read_name contains a forbidden control character")


def _canonical_chromosome(raw: str) -> str | None:
    value = raw.strip()
    if re.fullmatch(_CANONICAL_CHROMOSOME, value) is None:
        return None
    return value if value.startswith("chr") else f"chr{value}"


def _integer(raw: str | None, *, field: str, line_number: int, minimum: int = 0) -> int:
    try:
        value = int(raw or "")
    except ValueError as exc:
        raise ValueError(f"Nanopolish line {line_number} has an invalid {field}") from exc
    if value < minimum:
        raise ValueError(f"Nanopolish line {line_number} has {field} below {minimum}")
    return value


def _finite_float(raw: str | None, *, field: str, line_number: int) -> float:
    try:
        value = float(raw or "")
    except ValueError as exc:
        raise ValueError(f"Nanopolish line {line_number} has an invalid {field}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Nanopolish line {line_number} has a non-finite {field}")
    return value


def parse_nanopolish_source(
    path: Path,
    *,
    source_id: str,
    genome_build: GenomeBuild,
    policy: MethylationMixturePolicy,
    upstream: NanopolishSourceMetadata | None = None,
) -> _ParsedSource:
    """Parse a Nanopolish call table without exporting read identifiers."""

    if re.fullmatch(_IDENTIFIER, source_id) is None:
        raise ValueError("source_id must be a safe pseudonymous identifier")
    if not path.is_file():
        raise ValueError("Nanopolish methylation source is missing or unreadable")
    initial_size = path.stat().st_size
    if initial_size > policy.maximum_input_bytes_per_source:
        raise ValueError(
            f"Nanopolish source has {initial_size} bytes; the locked in-memory safety limit is "
            f"{policy.maximum_input_bytes_per_source} bytes"
        )
    initial_sha256 = sha256_file(path)
    calls: dict[str, list[_NanopolishCall]] = {}
    marker_cache: dict[tuple[str, Literal["+", "-"], int, int, int, str], _MarkerKey] = {}
    seen: set[tuple[str, _MarkerKey]] = set()
    rows_total = 0
    canonical_rows = 0
    informative_rows = 0
    ambiguous_rows = 0
    skipped_noncanonical = 0

    with _open_nanopolish(path) as handle:
        lines = _BoundedNanopolishLines(handle, policy.maximum_input_bytes_per_source)
        header_line = next(lines, None)
        if header_line is None:
            raise ValueError("Nanopolish input has no header")
        try:
            # A header must be one bounded physical line, never a multiline CSV record.
            fieldnames = next(csv.reader([header_line], delimiter="\t", strict=True))
        except (csv.Error, StopIteration) as exc:
            raise ValueError("Nanopolish input has an invalid single-line header") from exc
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError("Nanopolish input contains duplicate column names")
        missing = sorted(_REQUIRED_NANOPOLISH_COLUMNS - set(fieldnames))
        if missing:
            raise ValueError("Nanopolish input is missing columns: " + ", ".join(missing))
        reader = csv.DictReader(lines, fieldnames=fieldnames, delimiter="\t")

        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Nanopolish line {line_number} has unexpected extra columns")
            if not any((value or "").strip() for value in row.values()):
                continue
            rows_total += 1
            if rows_total > policy.maximum_rows_per_source:
                raise ValueError(
                    "Nanopolish source exceeds the locked in-memory row safety limit of "
                    f"{policy.maximum_rows_per_source}"
                )
            chromosome = _canonical_chromosome(row.get("chromosome") or "")
            if chromosome is None:
                skipped_noncanonical += 1
                continue
            start = _integer(row.get("start"), field="start", line_number=line_number)
            end = _integer(row.get("end"), field="end", line_number=line_number)
            if end < start:
                raise ValueError(f"Nanopolish line {line_number} has invalid coordinates")
            motif_count = _integer(
                row.get("num_motifs"),
                field="num_motifs",
                line_number=line_number,
                minimum=1,
            )
            _integer(
                row.get("num_calling_strands"),
                field="num_calling_strands",
                line_number=line_number,
                minimum=1,
            )
            _finite_float(
                row.get("log_lik_methylated"),
                field="log_lik_methylated",
                line_number=line_number,
            )
            _finite_float(
                row.get("log_lik_unmethylated"),
                field="log_lik_unmethylated",
                line_number=line_number,
            )
            llr = _finite_float(
                row.get("log_lik_ratio"), field="log_lik_ratio", line_number=line_number
            )
            strand = (row.get("strand") or "").strip()
            if strand not in {"+", "-"}:
                raise ValueError(f"Nanopolish line {line_number} has invalid strand")
            typed_strand = cast(Literal["+", "-"], strand)
            raw_read_name = row.get("read_name") or ""
            _reject_read_name_controls(raw_read_name)
            read_name = raw_read_name.strip()
            if not read_name:
                raise ValueError(f"Nanopolish line {line_number} has an empty read_name")
            sequence = (row.get("sequence") or "").strip().upper()
            if not sequence:
                raise ValueError(f"Nanopolish line {line_number} has an empty sequence")
            sequence_sha256 = hashlib.sha256(sequence.encode("utf-8")).hexdigest()

            event_threshold = policy.log_likelihood_ratio_threshold * motif_count
            if llr >= event_threshold:
                state = NanopolishCallState.METHYLATED
                informative_rows += 1
            elif llr <= -event_threshold:
                state = NanopolishCallState.UNMETHYLATED
                informative_rows += 1
            else:
                state = NanopolishCallState.AMBIGUOUS
                ambiguous_rows += 1
            canonical_rows += 1
            marker_tuple = (
                chromosome,
                typed_strand,
                start,
                end,
                motif_count,
                sequence_sha256,
            )
            marker = marker_cache.setdefault(marker_tuple, _MarkerKey(*marker_tuple))
            duplicate_key = (read_name, marker)
            if duplicate_key in seen:
                raise ValueError(
                    f"Nanopolish line {line_number} duplicates a marker within one read group"
                )
            seen.add(duplicate_key)
            calls.setdefault(read_name, []).append(
                _NanopolishCall(
                    marker=marker,
                    state=state,
                )
            )

    if not calls:
        raise ValueError("Nanopolish source contains no canonical methylation calls")
    final_size = path.stat().st_size
    final_sha256 = sha256_file(path)
    if final_size != initial_size or final_sha256 != initial_sha256:
        raise ValueError("Nanopolish source changed while it was being parsed")
    fingerprint = FileFingerprint(size_bytes=initial_size, sha256=initial_sha256)
    summary = MethylationSourceSummary(
        source_id=source_id,
        genome_build=genome_build,
        upstream=upstream or NanopolishSourceMetadata(),
        fingerprint=fingerprint,
        decoded_size_bytes=lines.decoded_size_bytes,
        rows_total=rows_total,
        canonical_call_rows=canonical_rows,
        informative_call_rows=informative_rows,
        ambiguous_call_rows=ambiguous_rows,
        skipped_noncanonical_rows=skipped_noncanonical,
        read_group_count=len(calls),
    )
    return _ParsedSource(
        summary=summary,
        calls_by_read={name: tuple(items) for name, items in calls.items()},
    )


def _split_reads(
    read_names: Sequence[str], *, fraction: float, seed: int
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    if len(read_names) < 2:
        raise ValueError("Each source needs at least two read groups for calibration and mixture")
    ordered = sorted(read_names)
    random.Random(seed).shuffle(ordered)
    reference_count = _calibration_read_group_count(len(ordered), fraction)
    return tuple(ordered[:reference_count]), tuple(ordered[reference_count:])


def _selection_digest(
    source_a_reads: Sequence[str], source_b_reads: Sequence[str], *, context: str
) -> str:
    digest = hashlib.sha256()
    digest.update(f"context={context}\n".encode())
    for label, names in (("A", source_a_reads), ("B", source_b_reads)):
        for name in sorted(names):
            _reject_read_name_controls(name)
            digest.update(label.encode())
            digest.update(b"\0")
            digest.update(name.encode("utf-8"))
            digest.update(b"\n")
    return digest.hexdigest()


def _aggregate(source: ReadCallSource, read_names: Sequence[str]) -> dict[_MarkerKey, _Counts]:
    counts: dict[_MarkerKey, _Counts] = {}
    for read_name in read_names:
        for call in source.calls_by_read[read_name]:
            counts.setdefault(call.marker, _Counts()).add(call.state)
    return counts


def _reference_markers(
    source_a_counts: Mapping[_MarkerKey, _Counts],
    source_b_counts: Mapping[_MarkerKey, _Counts],
    *,
    source_a_total_read_groups: int,
    source_b_total_read_groups: int,
    policy: MethylationMixturePolicy,
) -> list[MethylationReferenceMarker]:
    selected: list[MethylationReferenceMarker] = []
    for key in set(source_a_counts) & set(source_b_counts):
        a = source_a_counts[key]
        b = source_b_counts[key]
        if (
            a.valid < policy.minimum_reference_valid_calls
            or b.valid < policy.minimum_reference_valid_calls
        ):
            continue
        assert a.beta is not None
        assert b.beta is not None
        delta = a.beta - b.beta
        if not _meets_minimum_absolute_delta_beta(delta, policy.minimum_absolute_delta_beta):
            continue
        selected.append(
            MethylationReferenceMarker(
                marker_id=key.marker_id,
                chromosome=key.chromosome,
                strand=key.strand,
                start=key.start,
                end=key.end,
                motif_count=key.motif_count,
                sequence_sha256=key.sequence_sha256,
                source_a_modified_calls=a.modified,
                source_a_unmethylated_calls=a.unmethylated,
                source_a_ambiguous_calls=a.ambiguous,
                source_a_valid_calls=a.valid,
                source_a_total_read_groups=source_a_total_read_groups,
                source_b_modified_calls=b.modified,
                source_b_unmethylated_calls=b.unmethylated,
                source_b_ambiguous_calls=b.ambiguous,
                source_b_valid_calls=b.valid,
                source_b_total_read_groups=source_b_total_read_groups,
                source_a_beta=a.beta,
                source_b_beta=b.beta,
                delta_beta=delta,
            )
        )
    selected.sort(key=lambda marker: (-abs(marker.delta_beta), marker.marker_id))
    if policy.maximum_markers is not None:
        selected = selected[: policy.maximum_markers]
    return sorted(selected, key=lambda marker: marker.marker_id)


def _reference_marker_digest(markers: Sequence[MethylationReferenceMarker]) -> str:
    digest = hashlib.sha256()
    for marker in markers:
        values = (
            marker.marker_id,
            str(marker.source_a_modified_calls),
            str(marker.source_a_unmethylated_calls),
            str(marker.source_a_ambiguous_calls),
            str(marker.source_a_valid_calls),
            str(marker.source_a_total_read_groups),
            str(marker.source_b_modified_calls),
            str(marker.source_b_unmethylated_calls),
            str(marker.source_b_ambiguous_calls),
            str(marker.source_b_valid_calls),
            str(marker.source_b_total_read_groups),
        )
        digest.update("\t".join(values).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _binomial_posterior_variance(successes: int, total: int) -> float:
    alpha = successes + 0.5
    beta = total - successes + 0.5
    return alpha * beta / ((alpha + beta) ** 2 * (alpha + beta + 1))


def _channel_weight(
    source_a_count: int,
    source_a_total: int,
    source_b_count: int,
    source_b_total: int,
    mixture_count: int,
    mixture_total: int,
) -> float:
    variance = (
        _binomial_posterior_variance(source_a_count, source_a_total)
        + _binomial_posterior_variance(source_b_count, source_b_total)
        + _binomial_posterior_variance(mixture_count, mixture_total)
    )
    return 1.0 / max(variance, 1e-12)


def _rate_components(
    marker: MethylationReferenceMarker,
    counts: _Counts,
    *,
    mixture_total_read_groups: int,
) -> tuple[tuple[float, float, float, float], tuple[float, float, float, float]]:
    raw_counts = (
        (
            marker.source_a_modified_calls,
            marker.source_b_modified_calls,
            counts.modified,
        ),
        (
            marker.source_a_unmethylated_calls,
            marker.source_b_unmethylated_calls,
            counts.unmethylated,
        ),
    )
    components: list[tuple[float, float, float, float]] = []
    for source_a_count, source_b_count, mixture_count in raw_counts:
        components.append(
            (
                source_a_count / marker.source_a_total_read_groups,
                source_b_count / marker.source_b_total_read_groups,
                mixture_count / mixture_total_read_groups,
                _channel_weight(
                    source_a_count,
                    marker.source_a_total_read_groups,
                    source_b_count,
                    marker.source_b_total_read_groups,
                    mixture_count,
                    mixture_total_read_groups,
                ),
            )
        )
    return components[0], components[1]


def _point_estimate(
    markers: Sequence[MethylationReferenceMarker],
    observations: Mapping[str, _Counts],
    *,
    minimum_mixture_valid_calls: int,
    mixture_total_read_groups: int,
) -> tuple[
    float | None,
    float | None,
    float | None,
    list[MethylationMixtureMarkerObservation],
]:
    numerator = 0.0
    denominator = 0.0
    residual_terms: list[tuple[float, float]] = []
    serialized: list[MethylationMixtureMarkerObservation] = []
    usable: list[
        tuple[
            MethylationReferenceMarker,
            tuple[tuple[float, float, float, float], tuple[float, float, float, float]],
        ]
    ] = []
    for marker in markers:
        counts = observations.get(marker.marker_id, _Counts())
        serialized.append(
            MethylationMixtureMarkerObservation(
                marker_id=marker.marker_id,
                modified_calls=counts.modified,
                unmethylated_calls=counts.unmethylated,
                valid_calls=counts.valid,
                ambiguous_calls=counts.ambiguous,
                beta=counts.beta,
            )
        )
        if counts.valid < minimum_mixture_valid_calls or counts.beta is None:
            continue
        components = _rate_components(
            marker,
            counts,
            mixture_total_read_groups=mixture_total_read_groups,
        )
        for source_a_rate, source_b_rate, mixture_rate, weight in components:
            delta = source_a_rate - source_b_rate
            numerator += weight * delta * (mixture_rate - source_b_rate)
            denominator += weight * delta**2
        usable.append((marker, components))
    if denominator <= 0:
        return None, None, None, serialized
    unconstrained = numerator / denominator
    estimate = min(max(unconstrained, 0.0), 1.0)
    for _marker, components in usable:
        for source_a_rate, source_b_rate, mixture_rate, weight in components:
            expected = source_b_rate + estimate * (source_a_rate - source_b_rate)
            residual_terms.append((weight, (mixture_rate - expected) ** 2))
    degrees_of_freedom = max(len(residual_terms) - 1, 1)
    standardized_rmse = math.sqrt(
        sum(weight * error for weight, error in residual_terms) / degrees_of_freedom
    )
    return unconstrained, estimate, standardized_rmse, serialized


def _state_rates(
    *, modified: int, unmethylated: int, ambiguous: int, total: int
) -> tuple[float, float, float, float]:
    missing = total - modified - unmethylated - ambiguous
    if missing < 0:
        raise ValueError("Methylation call-state counts exceed their read-group total")
    return tuple(count / total for count in (modified, unmethylated, ambiguous, missing))  # type: ignore[return-value]


def _draw_model_state_rates(
    probabilities: Sequence[float], *, concentration: int, rng: random.Random
) -> tuple[float, float]:
    """Draw M/U rates around one fitted four-state mixture distribution.

    Calibration rates remain fixed. This avoids treating independently perturbed noisy
    calibration endpoints as if they were exact regressors, which systematically pulls the
    coefficient toward the centre. Zero-probability states remain zero in this conditional
    model; positive states use a Dirichlet concentration equal to the mixture read-group count.
    """

    if concentration < 1 or len(probabilities) != 4:
        raise ValueError("Conditional Dirichlet draws need four states and positive concentration")
    if any(value < 0 or value > 1 for value in probabilities) or not math.isclose(
        sum(probabilities), 1.0, abs_tol=1e-12
    ):
        raise ValueError("Conditional Dirichlet probabilities must form a distribution")
    draws = [
        rng.gammavariate(probability * concentration, 1.0) if probability > 0 else 0.0
        for probability in probabilities
    ]
    normalizer = sum(draws)
    if normalizer <= 0:
        raise ValueError("Conditional Dirichlet draw has no positive state")
    return draws[0] / normalizer, draws[1] / normalizer


def _draw_estimate(
    markers: Sequence[MethylationReferenceMarker],
    observations: Mapping[str, _Counts],
    *,
    minimum_mixture_valid_calls: int,
    mixture_total_read_groups: int,
    fitted_source_a_fraction: float,
    rng: random.Random,
) -> float | None:
    numerator = 0.0
    denominator = 0.0
    used = 0
    for marker in markers:
        counts = observations.get(marker.marker_id)
        if counts is None or counts.valid < minimum_mixture_valid_calls:
            continue
        source_a_state_rates = _state_rates(
            modified=marker.source_a_modified_calls,
            unmethylated=marker.source_a_unmethylated_calls,
            ambiguous=marker.source_a_ambiguous_calls,
            total=marker.source_a_total_read_groups,
        )
        source_b_state_rates = _state_rates(
            modified=marker.source_b_modified_calls,
            unmethylated=marker.source_b_unmethylated_calls,
            ambiguous=marker.source_b_ambiguous_calls,
            total=marker.source_b_total_read_groups,
        )
        fitted_state_rates = tuple(
            source_b_rate + fitted_source_a_fraction * (source_a_rate - source_b_rate)
            for source_a_rate, source_b_rate in zip(
                source_a_state_rates, source_b_state_rates, strict=True
            )
        )
        mixture_rates = _draw_model_state_rates(
            fitted_state_rates,
            concentration=mixture_total_read_groups,
            rng=rng,
        )
        observed_components = _rate_components(
            marker,
            counts,
            mixture_total_read_groups=mixture_total_read_groups,
        )
        for mixture_rate, observed in zip(mixture_rates, observed_components, strict=True):
            source_a_rate, source_b_rate, _observed_mixture_rate, weight = observed
            delta = source_a_rate - source_b_rate
            numerator += weight * delta * (mixture_rate - source_b_rate)
            denominator += weight * delta**2
        used += 1
    if used == 0 or denominator <= 0:
        return None
    return min(max(numerator / denominator, 0.0), 1.0)


def _quantile(sorted_values: Sequence[float], probability: float) -> float:
    position = (len(sorted_values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    fraction = position - lower
    return sorted_values[lower] * (1 - fraction) + sorted_values[upper] * fraction


def _no_call_level(
    context: _LevelContext,
    *,
    marker_observations: list[MethylationMixtureMarkerObservation],
    markers_used: int,
    successful_draws: int,
    reason: str,
) -> MethylationMixtureLevelResult:
    total = len(context.source_a_reads) + len(context.source_b_reads)
    return MethylationMixtureLevelResult(
        level_id=context.level_id,
        replicate=context.replicate,
        seed=context.seed,
        target_source_a_fraction=context.target_fraction,
        source_a_read_groups=len(context.source_a_reads),
        source_b_read_groups=len(context.source_b_reads),
        total_read_groups=total,
        realized_source_a_read_group_fraction=len(context.source_a_reads) / total,
        selection_sha256=context.selection_sha256,
        status=ModuleRunStatus.NO_CALL,
        confidence_level=context.confidence_level,
        markers_available=context.markers_available,
        markers_used=markers_used,
        valid_mixture_calls=sum(item.valid_calls for item in marker_observations),
        successful_uncertainty_draws=successful_draws,
        no_call_reason=reason,
        marker_observations=marker_observations,
    )


def _estimate_level(
    *,
    level_id: str,
    replicate: int,
    seed: int,
    target_fraction: float,
    source_a_reads: Sequence[str],
    source_b_reads: Sequence[str],
    source_a: ReadCallSource,
    source_b: ReadCallSource,
    reference_markers: Sequence[MethylationReferenceMarker],
    policy: MethylationMixturePolicy,
) -> MethylationMixtureLevelResult:
    selected = _selection_digest(source_a_reads, source_b_reads, context=level_id)
    a_counts = _aggregate(source_a, source_a_reads)
    b_counts = _aggregate(source_b, source_b_reads)
    by_id = {marker.marker_id: marker for marker in reference_markers}
    observations: dict[str, _Counts] = {}
    for key in set(a_counts) | set(b_counts):
        marker_id = key.marker_id
        if marker_id not in by_id:
            continue
        combined = _Counts()
        for source_counts in (a_counts.get(key), b_counts.get(key)):
            if source_counts is not None:
                combined.modified += source_counts.modified
                combined.unmethylated += source_counts.unmethylated
                combined.ambiguous += source_counts.ambiguous
        observations[marker_id] = combined

    unconstrained, point, fit_standardized_rmse, serialized = _point_estimate(
        reference_markers,
        observations,
        minimum_mixture_valid_calls=policy.minimum_mixture_valid_calls,
        mixture_total_read_groups=len(source_a_reads) + len(source_b_reads),
    )
    used = sum(item.valid_calls >= policy.minimum_mixture_valid_calls for item in serialized)
    context = _LevelContext(
        level_id=level_id,
        replicate=replicate,
        seed=seed,
        target_fraction=target_fraction,
        source_a_reads=source_a_reads,
        source_b_reads=source_b_reads,
        selection_sha256=selected,
        confidence_level=policy.confidence_level,
        markers_available=len(reference_markers),
    )
    if len(reference_markers) < policy.minimum_markers:
        return _no_call_level(
            context,
            marker_observations=serialized,
            markers_used=used,
            successful_draws=0,
            reason=(
                f"Only {len(reference_markers)} reference marker(s) passed the locked filters; "
                f"at least {policy.minimum_markers} are required."
            ),
        )
    if (
        used < policy.minimum_markers
        or unconstrained is None
        or point is None
        or fit_standardized_rmse is None
    ):
        return _no_call_level(
            context,
            marker_observations=serialized,
            markers_used=used,
            successful_draws=0,
            reason=(
                f"Only {used} mixture marker(s) reached the valid-call floor; at least "
                f"{policy.minimum_markers} are required."
            ),
        )
    if (
        unconstrained < -policy.maximum_extrapolation
        or unconstrained > 1 + policy.maximum_extrapolation
    ):
        return _no_call_level(
            context,
            marker_observations=serialized,
            markers_used=used,
            successful_draws=0,
            reason=(
                f"The unconstrained source-A coefficient {unconstrained:.4f} lies beyond "
                f"the allowed extrapolation of {policy.maximum_extrapolation:.4f}."
            ),
        )
    if fit_standardized_rmse > policy.maximum_standardized_model_fit_rmse:
        return _no_call_level(
            context,
            marker_observations=serialized,
            markers_used=used,
            successful_draws=0,
            reason=(
                "Marker signals disagree with the paired-source linear model: standardized "
                f"fit RMSE {fit_standardized_rmse:.4f} exceeds "
                f"{policy.maximum_standardized_model_fit_rmse:.4f}."
            ),
        )

    level_warnings: list[str] = []
    if unconstrained < 0 or unconstrained > 1:
        level_warnings.append(
            f"The unconstrained coefficient {unconstrained:.4f} was bounded to {point:.4f}."
        )

    rng = random.Random(seed + 9_000_001)
    draws = [
        _draw_estimate(
            reference_markers,
            observations,
            minimum_mixture_valid_calls=policy.minimum_mixture_valid_calls,
            mixture_total_read_groups=len(source_a_reads) + len(source_b_reads),
            fitted_source_a_fraction=point,
            rng=rng,
        )
        for _ in range(policy.uncertainty_draws)
    ]
    successful = sorted(value for value in draws if value is not None)
    required = math.ceil(policy.uncertainty_draws * policy.minimum_successful_uncertainty_fraction)
    if len(successful) < required:
        return _no_call_level(
            context,
            marker_observations=serialized,
            markers_used=used,
            successful_draws=len(successful),
            reason=(
                f"Only {len(successful)} of {policy.uncertainty_draws} uncertainty draws "
                f"were estimable; at least {required} are required."
            ),
        )
    tail = (1 - policy.confidence_level) / 2
    lower = _quantile(successful, tail)
    upper = _quantile(successful, 1 - tail)
    if upper - lower > policy.maximum_confidence_interval_width:
        return _no_call_level(
            context,
            marker_observations=serialized,
            markers_used=used,
            successful_draws=len(successful),
            reason=(
                f"The {policy.confidence_level:.1%} interval width {upper - lower:.4f} "
                f"exceeds {policy.maximum_confidence_interval_width:.4f}."
            ),
        )

    total = len(source_a_reads) + len(source_b_reads)
    return MethylationMixtureLevelResult(
        level_id=level_id,
        replicate=replicate,
        seed=seed,
        target_source_a_fraction=target_fraction,
        source_a_read_groups=len(source_a_reads),
        source_b_read_groups=len(source_b_reads),
        total_read_groups=total,
        realized_source_a_read_group_fraction=len(source_a_reads) / total,
        selection_sha256=selected,
        status=ModuleRunStatus.COMPLETED,
        estimated_source_a_fraction=point,
        unconstrained_source_a_fraction=unconstrained,
        confidence_interval_lower=lower,
        confidence_interval_upper=upper,
        confidence_level=policy.confidence_level,
        markers_available=len(reference_markers),
        markers_used=used,
        valid_mixture_calls=sum(item.valid_calls for item in serialized),
        model_fit_standardized_rmse=fit_standardized_rmse,
        successful_uncertainty_draws=len(successful),
        warnings=level_warnings,
        marker_observations=serialized,
    )


def _feasible_total(
    fractions: Sequence[float], source_a_available: int, source_b_available: int
) -> int:
    bounds: list[float] = []
    for fraction in fractions:
        if fraction > 0:
            bounds.append(source_a_available / fraction)
        if fraction < 1:
            bounds.append(source_b_available / (1 - fraction))
    total = math.floor(min(bounds))
    while total > 1:
        if all(
            round(total * fraction) <= source_a_available
            and total - round(total * fraction) <= source_b_available
            for fraction in fractions
        ):
            return total
        total -= 1
    raise ValueError("Held-out source pools cannot fund every requested mixture level")


_LIMITATIONS: tuple[str, ...] = (
    "The estimated value is the fraction of held-out read groups drawn from source A. It is "
    "not automatically a DNA-mass fraction, cell fraction, tumour purity or blast fraction.",
    "Calibration and mixture pools are disjoint read groups from the same two biological "
    "sources. Independent donors and orthogonal truth are still required for generalisation.",
    "The conditional Dirichlet Monte Carlo interval holds the calibration marker rates fixed "
    "and propagates only mixture-pool sampling under the fitted four-state linear model. It "
    "does not include calibration-source uncertainty, marker correlation or model "
    "misspecification.",
    "The standardized residual RMSE is a technical disagreement heuristic, not a calibrated "
    "chi-square statistic; methylated and unmethylated channels are correlated.",
    "Nanopolish coordinates and the declared genome build are accepted under the input "
    "contract; unknown upstream metadata remains an explicit provenance limitation.",
    "Copy number, ploidy, allelic imbalance, cell composition and platform transfer are not "
    "corrected by this paired-source technical model.",
    "The current parser is deliberately bounded by the policy byte and row limits and keeps "
    "accepted calls in memory. Multi-gigabyte inputs require a validated streaming or on-disk "
    "backend.",
    "Technical recovery thresholds and the resulting PASS/FAIL/NOT_EVALUABLE assessment are "
    "engineering controls only; they are not validated analytical acceptance criteria.",
)


def _report_warnings(
    *,
    policy: MethylationMixturePolicy,
    source_a: MethylationSourceSummary,
    source_b: MethylationSourceSummary,
    total_mixture_read_groups: int,
    levels: Sequence[MethylationMixtureLevelResult],
) -> list[str]:
    warnings = [policy.note]
    for label, metadata in (("source A", source_a.upstream), ("source B", source_b.upstream)):
        unknown = metadata.unknown_fields()
        if unknown:
            warnings.append(
                f"Upstream provenance is incomplete for {label}: " + ", ".join(unknown) + "."
            )
    if policy.total_mixture_read_groups is None:
        warnings.append(
            "The constant mixture budget was derived from the held-out pools: "
            f"{total_mixture_read_groups} read groups per level. Pin it before comparing "
            "separate experiments."
        )
    if source_a.skipped_noncanonical_rows or source_b.skipped_noncanonical_rows:
        warnings.append(
            "Non-canonical-contig rows were excluded and counted in each source summary."
        )
    no_call_count = sum(level.status == ModuleRunStatus.NO_CALL for level in levels)
    if no_call_count:
        warnings.append(
            f"{no_call_count} of {len(levels)} mixture level(s) returned NO_CALL; inspect the "
            "level-specific reasons rather than treating them as zero source contribution."
        )
    if no_call_count < len(levels):
        warnings.append(
            "Bias, MAE and RMSE summarize completed levels only; use quantifiable_fraction "
            "to assess coverage of the requested grid."
        )
    return warnings


def run_nanopolish_mixture_analysis(
    source_a_path: Path,
    source_b_path: Path,
    *,
    source_a_id: str,
    source_b_id: str,
    source_a_genome_build: GenomeBuild,
    source_b_genome_build: GenomeBuild,
    analysis_id: str,
    policy: MethylationMixturePolicy,
    sources_declared_biologically_distinct: bool,
    source_a_metadata: NanopolishSourceMetadata | None = None,
    source_b_metadata: NanopolishSourceMetadata | None = None,
    software_version: str = __version__,
    git_commit: str = "UNKNOWN",
) -> MethylationMixtureReport:
    """Run the complete paired-source call-table experiment."""

    if re.fullmatch(_ANALYSIS_IDENTIFIER, analysis_id) is None:
        raise ValueError("analysis_id must be a safe pseudonymous identifier")
    if not sources_declared_biologically_distinct:
        raise ValueError(
            "The operator must explicitly declare that source A and source B are biologically "
            "distinct; file names and hashes cannot establish this relationship"
        )
    if (
        re.fullmatch(_IDENTIFIER, source_a_id) is None
        or re.fullmatch(_IDENTIFIER, source_b_id) is None
    ):
        raise ValueError("Both source IDs must be safe pseudonymous identifiers")
    if source_a_id == source_b_id:
        raise ValueError("The two biological sources must have different source IDs")
    if source_a_genome_build != source_b_genome_build:
        raise ValueError("The two methylation sources use different genome builds")
    if re.fullmatch(_PROVENANCE_VALUE, software_version) is None:
        raise ValueError("software_version must be a single safe provenance value")
    if re.fullmatch(r"(?:UNKNOWN|[0-9a-f]{40})", git_commit) is None:
        raise ValueError("git_commit must be UNKNOWN or a full lowercase Git commit")
    metadata_a = source_a_metadata or NanopolishSourceMetadata()
    metadata_b = source_b_metadata or NanopolishSourceMetadata()
    metadata_mismatches = _known_metadata_mismatches(metadata_a, metadata_b)
    if metadata_mismatches:
        raise ValueError(
            "The two sources have incompatible known upstream metadata: "
            + ", ".join(metadata_mismatches)
        )
    source_a = parse_nanopolish_source(
        source_a_path,
        source_id=source_a_id,
        genome_build=source_a_genome_build,
        policy=policy,
        upstream=metadata_a,
    )
    source_b = parse_nanopolish_source(
        source_b_path,
        source_id=source_b_id,
        genome_build=source_b_genome_build,
        policy=policy,
        upstream=metadata_b,
    )
    if source_a.summary.fingerprint.sha256 == source_b.summary.fingerprint.sha256:
        raise ValueError("The two source files are byte-identical, not distinct inputs")

    reference_a, pool_a = _split_reads(
        tuple(source_a.calls_by_read),
        fraction=policy.calibration_read_group_fraction,
        seed=policy.seed + 1_000_003,
    )
    reference_b, pool_b = _split_reads(
        tuple(source_b.calls_by_read),
        fraction=policy.calibration_read_group_fraction,
        seed=policy.seed + 2_000_003,
    )
    calibration_digest = _selection_digest(
        reference_a, reference_b, context=f"{analysis_id}:calibration"
    )
    markers = _reference_markers(
        _aggregate(source_a, reference_a),
        _aggregate(source_b, reference_b),
        source_a_total_read_groups=len(reference_a),
        source_b_total_read_groups=len(reference_b),
        policy=policy,
    )

    feasible = _feasible_total(policy.source_a_fractions, len(pool_a), len(pool_b))
    total = policy.total_mixture_read_groups or feasible
    if total > feasible:
        raise ValueError(
            f"total_mixture_read_groups {total} exceeds the feasible constant budget of {feasible}"
        )
    for fraction in policy.source_a_fractions:
        a_count = round(total * fraction)
        b_count = total - a_count
        if 0 < fraction < 1 and (a_count == 0 or b_count == 0):
            raise ValueError(
                f"A budget of {total} cannot represent interior source fraction {fraction}; "
                "increase total_mixture_read_groups"
            )

    levels: list[MethylationMixtureLevelResult] = []
    for fraction_index, fraction in enumerate(policy.source_a_fractions):
        for replicate in range(1, policy.replicates + 1):
            level_seed = _mixture_level_seed(policy.seed, fraction_index, replicate)
            a_count = round(total * fraction)
            b_count = total - a_count
            rng_a = random.Random(level_seed)
            rng_b = random.Random(level_seed + 1)
            selected_a = tuple(rng_a.sample(pool_a, a_count)) if a_count else ()
            selected_b = tuple(rng_b.sample(pool_b, b_count)) if b_count else ()
            level_id = f"{analysis_id}.L{fraction_index + 1:03d}.R{replicate:03d}"
            levels.append(
                _estimate_level(
                    level_id=level_id,
                    replicate=replicate,
                    seed=level_seed,
                    target_fraction=fraction,
                    source_a_reads=selected_a,
                    source_b_reads=selected_b,
                    source_a=source_a,
                    source_b=source_b,
                    reference_markers=markers,
                    policy=policy,
                )
            )

    completed = [level for level in levels if level.status == ModuleRunStatus.COMPLETED]
    errors = [
        (level.estimated_source_a_fraction or 0.0) - level.realized_source_a_read_group_fraction
        for level in completed
    ]
    mean_bias = sum(errors) / len(errors) if errors else None
    mean_absolute_error = sum(abs(value) for value in errors) / len(errors) if errors else None
    root_mean_squared_error = (
        math.sqrt(sum(value**2 for value in errors) / len(errors)) if errors else None
    )
    confidence_interval_coverage_fraction = (
        sum(
            cast(float, level.confidence_interval_lower)
            <= level.realized_source_a_read_group_fraction
            <= cast(float, level.confidence_interval_upper)
            for level in completed
        )
        / len(completed)
        if completed
        else None
    )
    no_call_count = len(levels) - len(completed)
    quantifiable_fraction = len(completed) / len(levels)
    recovery_assessment, recovery_reasons = _technical_recovery_assessment(
        policy=policy,
        quantifiable_fraction=quantifiable_fraction,
        mean_bias=mean_bias,
        mean_absolute_error=mean_absolute_error,
        root_mean_squared_error=root_mean_squared_error,
        confidence_interval_coverage_fraction=confidence_interval_coverage_fraction,
    )
    warnings = _report_warnings(
        policy=policy,
        source_a=source_a.summary,
        source_b=source_b.summary,
        total_mixture_read_groups=total,
        levels=levels,
    )

    return MethylationMixtureReport(
        analysis_id=analysis_id,
        genome_build=source_a_genome_build,
        status=(
            MethylationMixtureReportStatus.COMPLETED
            if len(completed) == len(levels)
            else (
                MethylationMixtureReportStatus.PARTIAL
                if completed
                else MethylationMixtureReportStatus.NO_CALL
            )
        ),
        software_version=software_version,
        git_commit=git_commit,
        policy=policy,
        source_a=source_a.summary,
        source_b=source_b.summary,
        source_a_calibration_read_groups=len(reference_a),
        source_b_calibration_read_groups=len(reference_b),
        source_a_mixture_pool_read_groups=len(pool_a),
        source_b_mixture_pool_read_groups=len(pool_b),
        calibration_selection_sha256=calibration_digest,
        reference_marker_sha256=_reference_marker_digest(markers),
        reference_markers=markers,
        levels=levels,
        completed_levels=len(completed),
        no_call_levels=no_call_count,
        quantifiable_fraction=quantifiable_fraction,
        grid_fully_quantified=len(completed) == len(levels),
        mean_bias=mean_bias,
        mean_absolute_error=mean_absolute_error,
        root_mean_squared_error=root_mean_squared_error,
        confidence_interval_coverage_fraction=confidence_interval_coverage_fraction,
        technical_recovery_assessment=recovery_assessment,
        technical_recovery_reasons=recovery_reasons,
        warnings=warnings,
        limitations=list(_LIMITATIONS),
    )


def render_methylation_mixture_csv(report: MethylationMixtureReport, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "report_status",
                "technical_recovery_assessment",
                "quantifiable_fraction",
                "accuracy_metrics_scope",
                "software_version",
                "git_commit",
                "sensitive_output",
                "source_a_id",
                "source_a_sha256",
                "source_b_id",
                "source_b_sha256",
                "level_id",
                "replicate",
                "target_source_a_fraction",
                "realized_source_a_read_group_fraction",
                "estimated_source_a_fraction",
                "unconstrained_source_a_fraction",
                "confidence_interval_lower",
                "confidence_interval_upper",
                "status",
                "markers_used",
                "valid_mixture_calls",
                "model_fit_standardized_rmse",
                "no_call_reason",
            ]
        )
        for level in report.levels:
            writer.writerow(
                [
                    report.status.value,
                    report.technical_recovery_assessment.value,
                    report.quantifiable_fraction,
                    report.accuracy_metrics_scope,
                    report.software_version,
                    report.git_commit,
                    report.sensitive_output,
                    report.source_a.source_id,
                    report.source_a.fingerprint.sha256,
                    report.source_b.source_id,
                    report.source_b.fingerprint.sha256,
                    level.level_id,
                    level.replicate,
                    level.target_source_a_fraction,
                    level.realized_source_a_read_group_fraction,
                    level.estimated_source_a_fraction,
                    level.unconstrained_source_a_fraction,
                    level.confidence_interval_lower,
                    level.confidence_interval_upper,
                    level.status.value,
                    level.markers_used,
                    level.valid_mixture_calls,
                    level.model_fit_standardized_rmse,
                    level.no_call_reason,
                ]
            )
    return output_path


def render_methylation_mixture_html(report: MethylationMixtureReport, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def cell(value: object) -> str:
        return html.escape("" if value is None else str(value))

    def rounded(value: float | None) -> float | None:
        return None if value is None else round(value, 4)

    rows: list[str] = []
    points: list[str] = []
    for level in report.levels:
        note = level.no_call_reason or "; ".join(level.warnings)
        rows.append(
            "".join(
                [
                    "<tr>",
                    f"<td>{cell(level.level_id)}</td>",
                    f"<td>{level.target_source_a_fraction:.3f}</td>",
                    f"<td>{level.realized_source_a_read_group_fraction:.3f}</td>",
                    f"<td>{cell(rounded(level.estimated_source_a_fraction))}</td>",
                    f"<td>{cell(rounded(level.confidence_interval_lower))}</td>",
                    f"<td>{cell(rounded(level.confidence_interval_upper))}</td>",
                    f"<td>{cell(level.status.value)}</td>",
                    f"<td>{level.markers_used}</td>",
                    f"<td>{cell(note)}</td>",
                    "</tr>",
                ]
            )
        )
        if level.estimated_source_a_fraction is not None:
            x = level.realized_source_a_read_group_fraction * 100
            y = (1 - level.estimated_source_a_fraction) * 100
            points.append(
                f"<circle cx='{x:.3f}' cy='{y:.3f}' r='1.8'><title>"
                f"target={level.target_source_a_fraction:.4f}; "
                f"realized={level.realized_source_a_read_group_fraction:.4f}; "
                f"estimated={level.estimated_source_a_fraction:.4f}</title></circle>"
            )
    warnings = "".join(f"<li>{cell(item)}</li>" for item in report.warnings)
    limitations = "".join(f"<li>{cell(item)}</li>" for item in report.limitations)
    recovery_reasons = "".join(
        f"<li>{cell(item)}</li>" for item in report.technical_recovery_reasons
    )
    source_a_upstream = cell(
        "; ".join(f"{key}={value}" for key, value in report.source_a.upstream.model_dump().items())
    )
    source_b_upstream = cell(
        "; ".join(f"{key}={value}" for key, value in report.source_b.upstream.model_dump().items())
    )
    table_headers = "".join(
        f"<th>{cell(label)}</th>"
        for label in (
            "Level",
            "Target",
            "Realized",
            "Estimated",
            "CI low",
            "CI high",
            "Status",
            "Markers",
            "Note",
        )
    )
    rmse = rounded(report.root_mean_squared_error)
    recovery = cell(report.technical_recovery_assessment.value)
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>Methylation source-fraction experiment</title>
<style>
body{{font:15px system-ui,sans-serif;max-width:1180px;margin:2rem auto;
padding:0 1rem;color:#172033}}
.notice{{background:#fff3cd;border:1px solid #e8c65a;padding:.8rem;font-weight:700}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));
gap:.7rem;margin:1rem 0}}
.card{{border:1px solid #ccd4df;border-radius:8px;padding:.8rem}}
.card strong{{display:block;font-size:1.5rem}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border:1px solid #d7dde6;padding:.45rem;text-align:left}}
th{{background:#075985;color:white}}
svg{{width:min(100%,560px);border:1px solid #ccd4df;background:#fbfdff}}
circle{{fill:#0f766e;opacity:.8}}
.diag{{stroke:#64748b;stroke-dasharray:3 2}} .axis{{stroke:#172033}}
.muted{{color:#596579}} code{{background:#eef2f6;padding:.1rem .25rem}}
</style></head><body>
<h1>Paired-source methylation mixture</h1>
<p class="notice">RESEARCH USE ONLY — sensitive derived genomic output — controlled storage
required — not analytically or clinically validated</p>
<p>The reported value is a methylation-signal mixture coefficient for source A.
The known comparison value is the share of held-out Nanopolish read groups drawn from source A.
Neither value automatically represents tumour purity, cell fraction or DNA-mass fraction.</p>
<div class="cards">
<div class="card">Status<strong>{cell(report.status.value)}</strong></div>
<div class="card">Technical recovery<strong>{recovery}</strong></div>
<div class="card">Reference markers<strong>{len(report.reference_markers)}</strong></div>
<div class="card">Completed levels<strong>{report.completed_levels}</strong></div>
<div class="card">NO_CALL levels<strong>{report.no_call_levels}</strong></div>
<div class="card">Quantifiable<strong>{report.quantifiable_fraction:.1%}</strong></div>
<div class="card">RMSE<strong>{cell(rmse)}</strong></div>
</div>
<p class="muted">Accuracy metrics scope: {cell(report.accuracy_metrics_scope)}.</p>
<p class="muted">Empirical interval coverage over completed levels:
{cell(rounded(report.confidence_interval_coverage_fraction))}.</p>
<h2>Technical recovery assessment</h2>
<p>This checks the estimates against the known in-silico read-group fractions under locked,
unvalidated engineering thresholds. It is not an analytical or clinical validation verdict.</p>
<ul>{recovery_reasons or "<li>No locked recovery threshold was exceeded.</li>"}</ul>
<h2>Calibration</h2>
<svg viewBox="-10 -8 120 120" role="img"
aria-label="Known versus estimated source A fractions">
<line class="axis" x1="0" y1="100" x2="100" y2="100"/>
<line class="axis" x1="0" y1="0" x2="0" y2="100"/>
<line class="diag" x1="0" y1="100" x2="100" y2="0"/>{"".join(points)}
<text x="42" y="110" font-size="4">known fraction A</text>
<text x="-8" y="-2" font-size="4">estimated A</text>
</svg>
<h2>Mixture levels</h2>
<table><thead><tr>{table_headers}</tr></thead>
<tbody>{"".join(rows)}</tbody></table>
<h2>Provenance</h2>
<p>Analysis <code>{cell(report.analysis_id)}</code>;
build <code>{cell(report.genome_build.value)}</code>;
algorithm <code>{cell(report.algorithm_version)}</code>.</p>
<p>ONTSeq <code>{cell(report.software_version)}</code>;
Git commit <code>{cell(report.git_commit)}</code>.</p>
<p>Source A: <code>{cell(report.source_a.source_id)}</code>,
SHA-256 <code>{cell(report.source_a.fingerprint.sha256)}</code><br>
Source B: <code>{cell(report.source_b.source_id)}</code>,
SHA-256 <code>{cell(report.source_b.fingerprint.sha256)}</code></p>
<p>Source A upstream: <code>{source_a_upstream}</code><br>
Source B upstream: <code>{source_b_upstream}</code></p>
<p>Reference marker lock SHA-256:
<code>{cell(report.reference_marker_sha256)}</code></p>
<h2>Warnings</h2><ul>{warnings}</ul><h2>Limitations</h2><ul>{limitations}</ul>
</body></html>"""
    output_path.write_text(document, encoding="utf-8", newline="\n")
    return output_path
