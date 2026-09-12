"""Small-variant policy for a germline caller run on a tumour sample.

Clair3 is a germline small-variant caller. Running it on a leukaemia sample produces
variants; it does not produce *somatic* variants, and the difference is not a labelling
detail. A germline heterozygous SNP and a clonal somatic mutation in a pure tumour both sit
near allele fraction 0.5 and are indistinguishable from read counts alone. Deciding between
them needs a matched normal, a population-frequency filter, or a caller designed for the
tumour-only question -- none of which this assay currently has.

So the somatic label is not withheld by convention here, it is absent from the type system.
:class:`SomaticStatus` has exactly one member, and a test pins that it has exactly one. There
is no code path that produces a somatic call, because a boundary that depends on every future
caller remembering to respect it is not a boundary. This is the same trap ADR-022 records for
ClinVar's germline vocabulary, arriving from the other direction: there, germline
classifications were attached to somatic questions; here, a germline caller would answer one.

Two further constraints are carried in the data rather than in documentation.

*Depth.* An allele fraction is only meaningful against the depth that measured it. Every
accepted call is checked against the detection floor computed by
:mod:`ontseq_platform.quantitation` at its own observed depth, so a call sitting below what
its coverage can resolve is rejected with that reason rather than passed on with a quiet
caveat. This is per variant, not per run: an adaptive-sampling BAM carries ~80x inside the
panel and ~9x outside it, and one global threshold would be wrong in both places.

*Indels.* Long reads measurably struggle with them: 66% recall and 42% precision at 52x
(Abel et al., J Mol Diagn 2025), and 17.6% of small indels recovered at 21x (Kato et al.,
ASCO 2024). NPM1's canonical alteration is a 4 bp insertion and therefore falls in exactly
that class. Indels accordingly carry ``requires_orthogonal_confirmation``, which travels with
the call rather than living in a document nobody reads at the bench.

Nothing here executes Clair3. This is the policy and the boundary; the subprocess adapter and
its version probe are separate work, and no verification status is claimed for either until
CI has run the real binary (ADR-015).

Research use only. No threshold here is clinically validated, and no call produced under this
policy is reportable.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum

from ontseq_platform.quantitation import minimum_detectable_vaf

SNV = "snv"
INSERTION = "insertion"
DELETION = "deletion"
MNV = "mnv"

COMPLETED = "COMPLETED"
NO_CALL = "NO_CALL"

REJECT_NOT_PASS = "filter_is_not_pass"
REJECT_LOW_QUALITY = "quality_below_policy_minimum"
REJECT_LOW_DEPTH = "depth_below_policy_minimum"
REJECT_FEW_VARIANT_READS = "variant_reads_below_policy_minimum"
REJECT_BELOW_DETECTION_FLOOR = "allele_fraction_below_the_detection_floor_at_this_depth"
REJECT_UNRESOLVABLE_DEPTH = "depth_supports_no_allele_fraction_at_all"


class SmallVariantError(ValueError):
    """Raised when a record or policy cannot describe a variant call at all."""


class SomaticStatus(StrEnum):
    """Deliberately a single member.

    Adding a second one is a decision about the assay, not a refactor: it requires a matched
    normal, a population-frequency filter, or a tumour-only somatic caller, and it requires
    the validation that makes any of those trustworthy.
    """

    NOT_DETERMINED = "somatic_status_not_determined"


SOMATIC_REASON = (
    "This assay runs a germline caller without a matched normal, so a variant's somatic "
    "status is undetermined. A germline heterozygous SNP and a clonal somatic mutation are "
    "not distinguishable from read counts alone. Resolving this needs a matched normal, a "
    "population-frequency filter, or a tumour-only somatic caller."
)


@dataclass(frozen=True)
class SmallVariant:
    """One record as read from a caller's VCF, before any policy is applied."""

    chrom: str
    position: int
    reference: str
    alternate: str
    depth: int
    variant_reads: int
    quality: float
    filter_status: str

    def __post_init__(self) -> None:
        if self.position < 1:
            raise SmallVariantError(f"position must be one-based and positive: {self.position}")
        if not self.reference or not self.alternate:
            raise SmallVariantError("reference and alternate alleles must both be present")
        if self.depth <= 0:
            raise SmallVariantError(f"depth must be positive: {self.depth}")
        if not 0 <= self.variant_reads <= self.depth:
            raise SmallVariantError(
                f"variant_reads ({self.variant_reads}) must lie within depth ({self.depth})"
            )

    @property
    def vaf(self) -> float:
        return self.variant_reads / self.depth

    @property
    def variant_class(self) -> str:
        if len(self.reference) == len(self.alternate):
            return SNV if len(self.reference) == 1 else MNV
        return INSERTION if len(self.alternate) > len(self.reference) else DELETION

    @property
    def length_change(self) -> int:
        return len(self.alternate) - len(self.reference)

    @property
    def is_indel(self) -> bool:
        return self.variant_class in (INSERTION, DELETION)


def _as_number(raw: Mapping[str, object], key: str, default: float) -> float:
    """Read a numeric config value, or say which key is wrong.

    A config carrying ``min_depth: "thirty"`` should fail by name rather than as a TypeError
    three frames away, and a bool is not a number here even though Python says it is.
    """
    value = raw.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise SmallVariantError(f"clair3 profile field {key!r} must be a number, got {value!r}")
    return float(value)


def _as_flag(raw: Mapping[str, object], key: str, default: bool) -> bool:
    value = raw.get(key, default)
    if not isinstance(value, bool):
        raise SmallVariantError(
            f"clair3 profile field {key!r} must be true or false, got {value!r}"
        )
    return value


@dataclass(frozen=True)
class Clair3Policy:
    """Acceptance thresholds. Engineering defaults, never clinical ones."""

    profile_id: str
    expected_version: str
    pass_only: bool = True
    min_quality: float = 10.0
    min_depth: int = 30
    min_variant_reads: int = 4
    indels_require_orthogonal_confirmation: bool = True
    error_rate: float = 0.01
    alpha: float = 0.01
    power: float = 0.95
    required_model_id: str | None = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, object]) -> Clair3Policy:
        """Build from a parsed config.

        Takes a mapping rather than a path so this module stays dependency-free and can be
        tested wherever the criteria table can. A caller that has yaml reads the file.
        """
        try:
            profile_id = str(raw["profile_id"])
            expected_version = str(raw["expected_version"])
        except KeyError as error:
            raise SmallVariantError(f"clair3 profile is missing {error}") from error
        if str(raw.get("caller_vocabulary", "germline")) != "germline":
            raise SmallVariantError(
                "clair3 is a germline caller; a profile claiming another vocabulary for it "
                "is wrong about the tool, not about the assay"
            )
        return cls(
            profile_id=profile_id,
            expected_version=expected_version,
            pass_only=_as_flag(raw, "pass_only", True),
            min_quality=_as_number(raw, "min_quality", 10.0),
            min_depth=int(_as_number(raw, "min_depth", 30)),
            min_variant_reads=int(_as_number(raw, "min_variant_reads", 4)),
            indels_require_orthogonal_confirmation=_as_flag(
                raw, "indels_require_orthogonal_confirmation", True
            ),
            error_rate=_as_number(raw, "error_rate", 0.01),
            alpha=_as_number(raw, "alpha", 0.01),
            power=_as_number(raw, "power", 0.95),
            required_model_id=(
                str(raw["required_model_id"]) if raw.get("required_model_id") else None
            ),
        )


@dataclass(frozen=True)
class AcceptedCall:
    """A variant that passed policy. Never somatic, never reportable, always caveated."""

    variant: SmallVariant
    somatic_status: SomaticStatus
    requires_orthogonal_confirmation: bool
    detection_floor: float

    #: Fixed under the technical policy, exactly as ADR-007 fixes it for Sniffles2.
    reportable: bool = False

    def caveat(self) -> str:
        parts = [SOMATIC_REASON]
        if self.requires_orthogonal_confirmation:
            parts.append(
                f"This is a {self.variant.variant_class} of {abs(self.variant.length_change)} bp. "
                "Long-read indel precision is measurably poor, so this call is a hypothesis "
                "requiring orthogonal confirmation, not a result."
            )
        parts.append(
            f"Observed allele fraction {self.variant.vaf:.3f} against a detection floor of "
            f"{self.detection_floor:.3f} at {self.variant.depth}x."
        )
        return " ".join(parts)


@dataclass(frozen=True)
class RejectedCall:
    variant: SmallVariant
    reason: str


def evaluate(variant: SmallVariant, policy: Clair3Policy) -> AcceptedCall | RejectedCall:
    """Apply the policy to one variant, and check it against its own depth.

    The detection floor is computed per variant rather than per run. On an adaptive-sampling
    BAM a single global threshold is wrong twice over: too lax outside the panel and too
    strict inside it.
    """
    if policy.pass_only and variant.filter_status.upper() != "PASS":
        return RejectedCall(variant=variant, reason=REJECT_NOT_PASS)
    if variant.quality < policy.min_quality:
        return RejectedCall(variant=variant, reason=REJECT_LOW_QUALITY)
    if variant.depth < policy.min_depth:
        return RejectedCall(variant=variant, reason=REJECT_LOW_DEPTH)
    if variant.variant_reads < policy.min_variant_reads:
        return RejectedCall(variant=variant, reason=REJECT_FEW_VARIANT_READS)

    floor = minimum_detectable_vaf(
        variant.depth,
        error_rate=policy.error_rate,
        alpha=policy.alpha,
        power=policy.power,
    )
    if floor is None:
        return RejectedCall(variant=variant, reason=REJECT_UNRESOLVABLE_DEPTH)
    if variant.vaf < floor:
        return RejectedCall(variant=variant, reason=REJECT_BELOW_DETECTION_FLOOR)

    return AcceptedCall(
        variant=variant,
        somatic_status=SomaticStatus.NOT_DETERMINED,
        requires_orthogonal_confirmation=(
            variant.is_indel and policy.indels_require_orthogonal_confirmation
        ),
        detection_floor=floor,
    )


@dataclass(frozen=True)
class SmallVariantCallSet:
    """The outcome of applying the policy to every record the caller emitted."""

    profile_id: str
    accepted: tuple[AcceptedCall, ...]
    rejected: tuple[RejectedCall, ...]

    @property
    def status(self) -> str:
        """``NO_CALL`` when nothing survived policy, mirroring ADR-007 for Sniffles2.

        An empty accepted set means the stage ran and found nothing it could stand behind. It
        does not mean the sample carries no small variants, and it is not ``COMPLETED`` with
        an empty result, which a reader would take as a negative finding.
        """
        return COMPLETED if self.accepted else NO_CALL

    def rejection_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for call in self.rejected:
            counts[call.reason] = counts.get(call.reason, 0) + 1
        return counts

    @property
    def indels_needing_confirmation(self) -> tuple[AcceptedCall, ...]:
        return tuple(call for call in self.accepted if call.requires_orthogonal_confirmation)


def apply_policy(variants: Sequence[SmallVariant], policy: Clair3Policy) -> SmallVariantCallSet:
    """Split a caller's records into what policy accepts and what it refuses, counting both.

    Rejected records are kept rather than dropped. A call set that silently discarded most of
    its input looks identical to one that had little input, and the two mean opposite things.
    """
    accepted: list[AcceptedCall] = []
    rejected: list[RejectedCall] = []
    for variant in variants:
        outcome = evaluate(variant, policy)
        if isinstance(outcome, AcceptedCall):
            accepted.append(outcome)
        else:
            rejected.append(outcome)
    return SmallVariantCallSet(
        profile_id=policy.profile_id,
        accepted=tuple(accepted),
        rejected=tuple(rejected),
    )
