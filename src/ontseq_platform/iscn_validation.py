from __future__ import annotations

import re
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from enum import StrEnum
from importlib import import_module
from typing import Protocol, cast


class ISCNValidationStatus(StrEnum):
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"
    NOT_RUN = "NOT_RUN"


@dataclass(frozen=True, slots=True)
class ISCNValidationResult:
    status: ISCNValidationStatus
    engine: str
    messages: tuple[str, ...] = ()


class _ExternalValidationResult(Protocol):
    valid: bool
    errors: Iterable[object]


_CHROMOSOME = r"(?:[1-9]|1[0-9]|2[0-2]|X|Y)"
_BAND = r"[pq][1-9][0-9]*(?:\.[0-9]+)?"
_BASE_RE = re.compile(r"^(?P<count>[0-9]{2,3}),(?P<sex>XX|XY|X|XXY|XYY)(?:,(?P<body>.+))?$")
_NUMERIC_RE = re.compile(rf"^[+-]{_CHROMOSOME}$")
_SINGLE_CHR_STRUCTURAL_RE = re.compile(
    rf"^(?P<event>del|dup|inv)\((?P<chrom>{_CHROMOSOME})\)"
    rf"\((?P<start>{_BAND})(?P<end>{_BAND})?\)$"
)
_TRANSLOCATION_RE = re.compile(
    rf"^t\((?P<chrom_a>{_CHROMOSOME});(?P<chrom_b>{_CHROMOSOME})\)"
    rf"\((?P<band_a>{_BAND});(?P<band_b>{_BAND})\)$"
)


def validate_subset(notation: str) -> ISCNValidationResult:
    """Validate the deliberately small ISCN subset emitted by ONTSeq.

    This validator is intentionally conservative and does not claim to implement the full
    ISCN 2024 grammar. It checks only forms the ONTSeq renderer itself can emit. A full
    external validator is preferred when the optional dependency is installed.
    """

    messages: list[str] = []
    if any(character.isspace() for character in notation):
        return ISCNValidationResult(
            status=ISCNValidationStatus.FAIL,
            engine="ontseq-subset-v0.2",
            messages=("Whitespace is not allowed in the emitted ONTSeq ISCN subset.",),
        )

    match = _BASE_RE.fullmatch(notation)
    if match is None:
        return ISCNValidationResult(
            status=ISCNValidationStatus.FAIL,
            engine="ontseq-subset-v0.2",
            messages=("Karyotype header is outside the implemented ONTSeq subset.",),
        )

    body = match.group("body")
    if not body:
        return ISCNValidationResult(
            status=ISCNValidationStatus.PASS,
            engine="ontseq-subset-v0.2",
        )

    fragments = body.split(",")
    for fragment in fragments:
        if _NUMERIC_RE.fullmatch(fragment):
            continue
        structural = _SINGLE_CHR_STRUCTURAL_RE.fullmatch(fragment)
        if structural:
            start = structural.group("start")
            end = structural.group("end")
            if end and start[0] != end[0] and structural.group("event") in {"del", "dup"}:
                messages.append(
                    f"{fragment}: simple deletion/duplication spans the centromere; "
                    "the ONTSeq subset does not render this automatically."
                )
            continue
        if _TRANSLOCATION_RE.fullmatch(fragment):
            continue
        messages.append(f"Unsupported or malformed fragment: {fragment}")

    if messages:
        return ISCNValidationResult(
            status=ISCNValidationStatus.FAIL,
            engine="ontseq-subset-v0.2",
            messages=tuple(messages),
        )
    return ISCNValidationResult(
        status=ISCNValidationStatus.PASS,
        engine="ontseq-subset-v0.2",
    )


def validate_with_iscn_authenticator(notation: str) -> ISCNValidationResult:
    """Validate using the optional MIT-licensed ``iscn-authenticator`` package.

    The dependency is kept optional so that a network/package-manager issue can never make
    the core analysis unavailable. When absent, callers can fall back to ``validate_subset``.
    """

    try:
        module = import_module("iscn_authenticator")
    except ModuleNotFoundError:
        return ISCNValidationResult(
            status=ISCNValidationStatus.NOT_RUN,
            engine="iscn-authenticator",
            messages=("Optional package iscn-authenticator is not installed.",),
        )

    raw_validator = module.__dict__.get("validate_karyotype")
    if not callable(raw_validator):
        return ISCNValidationResult(
            status=ISCNValidationStatus.WARN,
            engine="iscn-authenticator",
            messages=("Installed iscn-authenticator exposes no validate_karyotype callable.",),
        )
    validator = cast(Callable[[str], _ExternalValidationResult], raw_validator)

    try:
        result = validator(notation)
    except Exception as exc:  # pragma: no cover - defensive boundary around optional package
        return ISCNValidationResult(
            status=ISCNValidationStatus.WARN,
            engine="iscn-authenticator",
            messages=(f"External validator failed unexpectedly: {exc}",),
        )

    errors = tuple(str(error) for error in result.errors or ())
    if result.valid:
        return ISCNValidationResult(
            status=ISCNValidationStatus.PASS,
            engine="iscn-authenticator",
            messages=errors,
        )
    return ISCNValidationResult(
        status=ISCNValidationStatus.FAIL,
        engine="iscn-authenticator",
        messages=errors or ("External validator rejected the karyotype.",),
    )


def validate_iscn(notation: str, *, prefer_external: bool = True) -> ISCNValidationResult:
    if prefer_external:
        external = validate_with_iscn_authenticator(notation)
        if external.status != ISCNValidationStatus.NOT_RUN:
            return external
    return validate_subset(notation)
