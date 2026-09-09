"""Privacy-safe parsing of the four VCF breakend replacement forms."""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum

_LOCAL_SEQUENCE = r"[ACGTNacgtn]+"
_CONTIG = (
    r"[0-9A-Za-z!#$%&+./:;?@^_|~-]"
    r"[0-9A-Za-z!#$%&*+./:;=?@^_|~-]*"
)
_LOCAL_THEN_MATE = re.compile(
    rf"^(?P<local>{_LOCAL_SEQUENCE})(?P<bracket>[\[\]])"
    rf"(?P<chromosome>{_CONTIG}):(?P<position>[0-9]+)"
    rf"(?P=bracket)$",
    flags=re.ASCII,
)
_MATE_THEN_LOCAL = re.compile(
    rf"^(?P<bracket>[\[\]])(?P<chromosome>{_CONTIG}):"
    rf"(?P<position>[0-9]+)(?P=bracket)(?P<local>{_LOCAL_SEQUENCE})$",
    flags=re.ASCII,
)
_CONTIG_ONLY = re.compile(rf"^{_CONTIG}$", flags=re.ASCII)
_POSITION_ONLY = re.compile(r"^[0-9]+$", flags=re.ASCII)
_SYMBOLIC_BREAKEND_ALTS = frozenset({"<BND>", "<TRA>"})
_MAX_VCF_INTEGER = 2_147_483_647


class _BreakendAltForm(StrEnum):
    """Internal bracket placement used only to prove all four VCF forms were parsed."""

    LOCAL_THEN_OPEN = "local_then_open"
    LOCAL_THEN_CLOSE = "local_then_close"
    OPEN_THEN_LOCAL = "open_then_local"
    CLOSE_THEN_LOCAL = "close_then_local"


class BreakendParseError(ValueError):
    """A stable caller-normalization rejection reason for invalid BND representation."""

    def __init__(self, reason: str, message: str) -> None:
        super().__init__(message)
        self.reason = reason


@dataclass(frozen=True)
class _ParsedBreakend:
    """The mate locus and VCF syntax, with local/inserted sequence discarded."""

    mate_chromosome: str
    mate_position_0based: int
    alt_form: _BreakendAltForm


@dataclass(frozen=True)
class ResolvedBreakend:
    """A BND mate resolved from strict bracket syntax or a supported symbolic ALT."""

    mate_chromosome: str
    mate_position_0based: int


def _parse_breakend_alt(alternate: str) -> _ParsedBreakend:
    """Parse one VCF BND ALT without retaining sequence from either side of the mate.

    VCF combines whether the mate appears before or after the local replacement string
    with either ``[`` or ``]`` brackets. Those four forms are preserved as syntax only;
    no transcript or gene orientation is inferred.
    """

    match = _LOCAL_THEN_MATE.fullmatch(alternate)
    local_first = match is not None
    if match is None:
        match = _MATE_THEN_LOCAL.fullmatch(alternate)
    if match is None:
        raise BreakendParseError(
            "malformed_breakend_alt",
            "VCF breakend ALT must be one ASCII allele in one of the four bracket forms",
        )
    position = _positive_vcf_position(
        match.group("position"),
        reason="malformed_breakend_alt",
        field="breakend ALT mate position",
    )
    bracket = match.group("bracket")
    if local_first:
        alt_form = (
            _BreakendAltForm.LOCAL_THEN_OPEN
            if bracket == "["
            else _BreakendAltForm.LOCAL_THEN_CLOSE
        )
    else:
        alt_form = (
            _BreakendAltForm.OPEN_THEN_LOCAL
            if bracket == "["
            else _BreakendAltForm.CLOSE_THEN_LOCAL
        )
    return _ParsedBreakend(
        mate_chromosome=match.group("chromosome"),
        mate_position_0based=position - 1,
        alt_form=alt_form,
    )


def _info_text(value: str | bool | None, *, field: str) -> str:
    if value is None:
        raise BreakendParseError("missing_breakend_mate", f"symbolic breakend requires {field}")
    if not isinstance(value, str) or not value:
        raise BreakendParseError(
            "malformed_breakend_mate", f"breakend {field} must be one scalar value"
        )
    return value


def _positive_vcf_position(raw: str, *, reason: str, field: str) -> int:
    """Parse a positive VCF Integer without exposing Python's digit-limit exceptions."""

    if _POSITION_ONLY.fullmatch(raw) is None or len(raw) > 10:
        raise BreakendParseError(reason, f"{field} is not one positive 32-bit ASCII integer")
    try:
        position = int(raw)
    except (ValueError, OverflowError) as exc:
        raise BreakendParseError(
            reason, f"{field} is not one positive 32-bit ASCII integer"
        ) from exc
    if not 1 <= position <= _MAX_VCF_INTEGER:
        raise BreakendParseError(reason, f"{field} is not one positive 32-bit ASCII integer")
    return position


def _declared_chromosome(value: str | bool | None, *, required: bool) -> str | None:
    if value is None and not required:
        return None
    raw = _info_text(value, field="CHR2")
    if _CONTIG_ONLY.fullmatch(raw) is None:
        raise BreakendParseError(
            "malformed_breakend_mate", "breakend CHR2 is not one ASCII VCF contig name"
        )
    return raw


def _declared_position(value: str | bool | None, *, required: bool) -> int | None:
    if value is None and not required:
        return None
    raw = _info_text(value, field="END")
    return _positive_vcf_position(
        raw,
        reason="malformed_breakend_mate",
        field="breakend END",
    )


def _normalized_chromosome(chromosome: str) -> str:
    """Compare the two canonical naming styles without rewriting recorded coordinates."""

    return chromosome[3:] if chromosome.startswith("chr") else chromosome


def resolve_breakend(
    alternate: str,
    *,
    declared_chromosome: str | bool | None,
    declared_position: str | bool | None,
) -> ResolvedBreakend:
    """Resolve one BND ALT and fail closed on contradictory INFO mate coordinates."""

    if "[" in alternate or "]" in alternate:
        parsed = _parse_breakend_alt(alternate)
        # CHR2 explicitly names the mate chromosome and must agree whenever present.
        # END alone can have caller-specific local-event semantics, so compare END only
        # as part of a CHR2/END mate declaration.
        if declared_chromosome is not None:
            info_chromosome = _declared_chromosome(declared_chromosome, required=True)
            assert info_chromosome is not None
            if _normalized_chromosome(info_chromosome) != _normalized_chromosome(
                parsed.mate_chromosome
            ):
                raise BreakendParseError(
                    "conflicting_breakend_mate", "breakend ALT and CHR2 name different mates"
                )
            if declared_position is not None:
                info_position = _declared_position(declared_position, required=True)
                assert info_position is not None
                if info_position - 1 != parsed.mate_position_0based:
                    raise BreakendParseError(
                        "conflicting_breakend_mate", "breakend ALT and END name different mates"
                    )
        return ResolvedBreakend(
            mate_chromosome=parsed.mate_chromosome,
            mate_position_0based=parsed.mate_position_0based,
        )

    if alternate not in _SYMBOLIC_BREAKEND_ALTS:
        raise BreakendParseError(
            "unsupported_breakend_alt",
            "BND/TRA requires bracket syntax or the symbolic ALT <BND>/<TRA>",
        )
    chromosome = _declared_chromosome(declared_chromosome, required=True)
    position = _declared_position(declared_position, required=True)
    assert chromosome is not None and position is not None
    return ResolvedBreakend(
        mate_chromosome=chromosome,
        mate_position_0based=position - 1,
    )
