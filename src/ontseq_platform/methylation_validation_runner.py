"""Local execution of a prospectively locked four-source read-group experiment."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import random
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import Field

from . import __version__
from .methylation_mixture import (
    NANOPOLISH_ADAPTER_VERSION,
    MethylationMixturePolicy,
    MethylationReferenceMarker,
    NanopolishSourceMetadata,
    ReadCallSource,
    _aggregate,
    _estimate_level,
    _reference_marker_digest,
    _reference_markers,
    _selection_digest,
    parse_nanopolish_source,
)
from .methylation_validation import (
    SAMPLE_ROLES,
    MethylationValidationEvidence,
    MethylationValidationMatrix,
    MethylationValidationRegistration,
    MethylationValidationReport,
    SampleRole,
    ValidationLevelEvidence,
    ValidationSample,
    cohort_eligibility_reasons,
    content_sha256,
    evaluate_methylation_validation,
    validation_level_seed,
)
from .models import StrictModel
from .reference import sha256_file


class ValidationLocalInput(StrictModel):
    """Sensitive local bindings; paths never enter a public catalogue."""

    calls_path: Path
    input_format: Literal["nanopolish", "modkit", "modbam"]
    metadata_path: Path
    reference_path: Path | None = None
    adapter_policy_path: Path | None = None


class ValidationLocalInputs(StrictModel):
    samples: dict[SampleRole, ValidationLocalInput]
    sensitive_output: Literal[True] = True


def validation_software_sha256() -> str:
    """Fingerprint actual Python source and runtime, including an uncommitted checkout."""
    package = Path(__file__).parent
    try:
        pysam_version: str | None = importlib.metadata.version("pysam")
    except importlib.metadata.PackageNotFoundError:
        pysam_version = None
    manifest = {
        "source_files": {
            path.relative_to(package).as_posix(): sha256_file(path)
            for path in sorted(package.rglob("*.py"))
        },
        "python": platform.python_version(),
        "pysam": pysam_version,
        "dependencies": {name: importlib.metadata.version(name) for name in ("pydantic", "PyYAML")},
    }
    return hashlib.sha256(
        json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


class ValidationReadiness(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    decision: Literal["NO_CALL"] = "NO_CALL"
    status: Literal["NOT_REGISTERED"] = "NOT_REGISTERED"
    matrix_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    code_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    software_version: str = __version__
    expected_levels: int = Field(ge=1)
    reasons: list[str] = Field(
        default_factory=lambda: [
            "No eligible checksummed four-sample cohort has been registered.",
            "No independent experiment has been executed; accuracy and coverage are unavailable.",
        ]
    )
    research_only: Literal[True] = True
    sensitive_output: Literal[True] = True


def _parse_inputs(
    registration: MethylationValidationRegistration,
    inputs: ValidationLocalInputs,
) -> dict[SampleRole, ReadCallSource]:
    if set(inputs.samples) != set(SAMPLE_ROLES):
        raise ValueError("Local inputs must bind exactly the four registered sample roles")
    if len({item.input_format for item in inputs.samples.values()}) != 1:
        raise ValueError("An independent experiment cannot mix adapter formats")
    sources: dict[SampleRole, ReadCallSource] = {}
    for role in SAMPLE_ROLES:
        sources[role] = parse_validation_source(
            role,
            registration.cohort.samples[role],
            inputs.samples[role],
            registration.matrix.estimator_policy,
        )
    # Different files or specimen declarations cannot override an observed reused read ID.
    seen: set[str] = set()
    for role in SAMPLE_ROLES:
        names = set(sources[role].calls_by_read)
        if seen & names:
            raise ValueError("Read identifiers overlap between registered biological inputs")
        seen.update(names)
    return sources


def parse_validation_source(
    role: str,
    sample: ValidationSample,
    local: ValidationLocalInput,
    policy: MethylationMixturePolicy,
) -> ReadCallSource:
    """Parse one checksummed source identically for both prospective study designs."""
    if sha256_file(local.calls_path) != sample.input_sha256:
        raise ValueError(f"{role}: local input differs from the registered SHA-256")
    metadata_data = json.loads(local.metadata_path.read_text(encoding="utf-8-sig"))
    if local.input_format == "nanopolish":
        if local.adapter_policy_path is not None:
            raise ValueError("Nanopolish uses the registered estimator policy as adapter policy")
        upstream = NanopolishSourceMetadata.model_validate(metadata_data)
        if upstream.unknown_fields():
            raise ValueError(f"{role}: complete Nanopolish metadata is required")
        if sample.technical.adapter_name != "nanopolish_call_methylation_tsv":
            raise ValueError(f"{role}: registered adapter name does not match Nanopolish")
        if sample.technical.adapter_version != NANOPOLISH_ADAPTER_VERSION:
            raise ValueError(f"{role}: registered adapter version does not match implementation")
        if content_sha256(policy) != sample.technical.adapter_policy_sha256:
            raise ValueError(f"{role}: Nanopolish adapter policy checksum differs")
        for field in NanopolishSourceMetadata.model_fields:
            if field != "schema_version" and getattr(upstream, field) != getattr(
                sample.technical, field
            ):
                raise ValueError(f"{role}: metadata differs at {field}")
        return parse_nanopolish_source(
            local.calls_path,
            source_id=role,
            genome_build=sample.technical.reference_genome_build,
            policy=policy,
            upstream=upstream,
        )
    return _parse_modification_input(sample, role, local, policy, metadata_data)


def _parse_modification_input(
    sample: ValidationSample,
    role: str,
    local: ValidationLocalInput,
    policy: MethylationMixturePolicy,
    metadata_data: object,
) -> ReadCallSource:
    from .modbam import (
        MODBAM_COMPATIBILITY_FIELDS,
        ModbamAdapterPolicy,
        ModbamSourceMetadata,
        parse_modbam_source,
        parse_modkit_source,
    )

    if local.reference_path is None or local.adapter_policy_path is None:
        raise ValueError("modBAM/modkit input requires reference_path and adapter_policy_path")
    upstream = ModbamSourceMetadata.model_validate(metadata_data)
    adapter_policy = ModbamAdapterPolicy.model_validate_json(
        local.adapter_policy_path.read_text(encoding="utf-8-sig")
    )
    if content_sha256(adapter_policy) != sample.technical.adapter_policy_sha256:
        raise ValueError(f"{role}: modification adapter policy checksum differs")
    for field in MODBAM_COMPATIBILITY_FIELDS:
        if getattr(upstream, field) != getattr(sample.technical, field):
            raise ValueError(f"{role}: metadata differs at {field}")
    if upstream.source_modbam_sha256 != sample.source_modbam_sha256:
        raise ValueError(f"{role}: parent modBAM fingerprint differs from registration")
    parser = parse_modkit_source if local.input_format == "modkit" else parse_modbam_source
    source = parser(
        local.calls_path,
        source_id=role,
        genome_build=sample.technical.reference_genome_build,
        policy=policy,
        upstream=upstream,
        adapter_policy=adapter_policy,
        reference_path=local.reference_path,
    )
    if source.summary.input_format != sample.technical.adapter_name:
        raise ValueError(f"{role}: registered adapter name does not match parsed format")
    if source.summary.parser_version != sample.technical.adapter_version:
        raise ValueError(f"{role}: registered parser version does not match local runtime")
    return source


def execute_registered_validation(
    registration: MethylationValidationRegistration,
    inputs: ValidationLocalInputs,
) -> MethylationValidationReport:
    """Execute only eligible inputs; missing biological eligibility remains an explicit NO_CALL."""
    registration = MethylationValidationRegistration.model_validate_json(
        registration.model_dump_json()
    )
    inputs = ValidationLocalInputs.model_validate_json(inputs.model_dump_json())
    if cohort_eligibility_reasons(registration.cohort) or not registration.test_outcomes_unseen:
        return evaluate_methylation_validation(registration, None)
    if validation_software_sha256() != registration.code_sha256:
        raise ValueError(
            "Software/runtime changed after registration; a new study lock is required"
        )
    started = datetime.now(UTC)
    if started <= registration.registered_at:
        raise ValueError("Experiment cannot begin before its declared registration")
    sources = _parse_inputs(registration, inputs)
    report = _run_independent_sources(registration, sources, experiment_started_at=started)
    for role in SAMPLE_ROLES:
        if (
            sha256_file(inputs.samples[role].calls_path)
            != registration.cohort.samples[role].input_sha256
        ):
            raise ValueError(f"{role}: input changed during the registered experiment")
    if validation_software_sha256() != registration.code_sha256:
        raise ValueError("Software/runtime changed during the registered experiment")
    return report


def _run_independent_sources(
    registration: MethylationValidationRegistration,
    sources: Mapping[SampleRole, ReadCallSource],
    *,
    experiment_started_at: datetime,
) -> MethylationValidationReport:
    """Internal execution after local input identity and adapter checks."""
    matrix = registration.matrix
    policy = matrix.estimator_policy
    names = {role: tuple(sorted(sources[role].calls_by_read)) for role in SAMPLE_ROLES}
    if any(len(names[role]) < max(matrix.read_group_budgets) for role in ("test_a", "test_b")):
        return evaluate_methylation_validation(
            registration,
            execution_no_call_reason=(
                "Held-out sources cannot supply the largest registered constant read-group budget."
            ),
        )
    if not names["calibration_a"] or not names["calibration_b"]:
        return evaluate_methylation_validation(
            registration,
            execution_no_call_reason="Both calibration samples require retained read groups.",
        )
    markers = _reference_markers(
        _aggregate(sources["calibration_a"], names["calibration_a"]),
        _aggregate(sources["calibration_b"], names["calibration_b"]),
        source_a_total_read_groups=len(names["calibration_a"]),
        source_b_total_read_groups=len(names["calibration_b"]),
        policy=policy,
    )
    levels = generate_validation_levels(matrix, sources["test_a"], sources["test_b"], markers)
    # Rehash after parsing/execution in the public runner before serialising input attestations.
    evidence = MethylationValidationEvidence(
        registration_sha256=registration.lock_sha256,
        code_sha256=validation_software_sha256(),
        software_version=__version__,
        experiment_started_at=experiment_started_at,
        experiment_completed_at=datetime.now(UTC),
        observed_input_sha256={
            role: registration.cohort.samples[role].input_sha256 for role in SAMPLE_ROLES
        },
        read_group_counts={role: len(names[role]) for role in SAMPLE_ROLES},
        reference_markers=markers,
        reference_marker_sha256=_reference_marker_digest(markers),
        calibration_selection_sha256=_selection_digest(
            names["calibration_a"], names["calibration_b"], context="independent-calibration"
        ),
        levels=levels,
    )
    return evaluate_methylation_validation(registration, evidence)


def generate_validation_levels(
    matrix: MethylationValidationMatrix,
    test_a: ReadCallSource,
    test_b: ReadCallSource,
    markers: Sequence[MethylationReferenceMarker],
) -> list[ValidationLevelEvidence]:
    """Shared deterministic grid execution from already separated whole-read test pools."""
    names_a = tuple(sorted(test_a.calls_by_read))
    names_b = tuple(sorted(test_b.calls_by_read))
    levels: list[ValidationLevelEvidence] = []
    for depth_index, budget in enumerate(matrix.read_group_budgets):
        for fraction_index, fraction in enumerate(matrix.source_a_fractions):
            for replicate, replicate_seed in enumerate(matrix.replicate_seeds, 1):
                seed = validation_level_seed(matrix, budget, fraction, replicate)
                count_a = round(budget * fraction)
                selected_a = random.Random(seed).sample(names_a, count_a)
                selected_b = random.Random(seed + 1).sample(names_b, budget - count_a)
                # Fixed bounded ID independent of free-text accession/sample metadata.
                level_id = f"VALIDATION.D{depth_index + 1}.F{fraction_index + 1}.R{replicate}"
                result = _estimate_level(
                    level_id=level_id,
                    replicate=replicate,
                    seed=seed,
                    target_fraction=fraction,
                    source_a_reads=selected_a,
                    source_b_reads=selected_b,
                    source_a=test_a,
                    source_b=test_b,
                    reference_markers=markers,
                    policy=matrix.estimator_policy,
                )
                levels.append(
                    ValidationLevelEvidence(
                        read_group_budget=budget, replicate_seed=replicate_seed, result=result
                    )
                )
    return levels
