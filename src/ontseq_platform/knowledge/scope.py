"""What a knowledge base can and cannot say about a finding in *this* sample.

A knowledge base answers a question about a *record*. Reading its answer as a verdict about
the sample in front of you is the single most consequential mistake this layer can make, and
it is easy to make silently, so the vocabulary here exists to keep the two apart.

Three separations carry the design.

**An assertion belongs to its source's vocabulary, and stays there.** ClinVar classifies
under ACMG germline rules: *Pathogenic*, *Likely pathogenic*, *Uncertain significance*. Those
words mean "this variant causes this inherited condition". A somatic driver in AML is not
"pathogenic" in that sense, and a somatic finding annotated *Pathogenic* would read as though
it were. So an assertion is carried verbatim together with the vocabulary it was written in,
and nothing in this package translates one vocabulary into another.

**Origin is checked, not assumed.** ClinVar records their own origin — germline, somatic, or
not provided — and an assay is looking for one or the other. Pairing a germline assertion
with a somatic question is not wrong in itself; it is *informative in a different way*, and
the annotation says which. Where either side is silent the alignment is ``UNKNOWN``, never
quietly treated as agreement.

**A match is a measurement, not a fact.** Region-based records match copy-number findings by
overlap, and an overlap of 0.51 and one of 0.99 are different evidence. The match type and
the reciprocal overlap are both recorded, so a reader can see whether the annotation
describes the same event or merely the same neighbourhood.

Nothing here decides reportability. That decision needs somatic criteria — ELN, ICC, a local
gene list — that this repository does not have and must not invent.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class Origin(StrEnum):
    """Whether an assertion, or a question, concerns inherited or acquired variation."""

    GERMLINE = "germline"
    SOMATIC = "somatic"
    #: The source or the assay did not say. Distinct from either answer.
    UNKNOWN = "unknown"


class ScopeAlignment(StrEnum):
    """How a source assertion's origin relates to what the assay was asking."""

    #: The assertion is about the kind of variation this assay is looking for.
    ALIGNED = "aligned"
    #: A real annotation about a different kind of variation. Still worth showing — a
    #: germline pathogenic CNV under a somatic finding is a secondary finding, not noise —
    #: but never as though it answered the question that was asked.
    MISMATCHED = "mismatched"
    #: One side did not declare. Not evidence of agreement.
    UNKNOWN = "unknown"


class MatchType(StrEnum):
    """How closely a source record corresponds to the finding it was matched to.

    Ordered from strongest to weakest. The distinction matters because a knowledge base
    entry for a 3 Mb recurrent deletion says something quite specific, and the same entry
    matched to a 90 Mb chromosome-arm loss that happens to contain it says much less.
    """

    #: Same coordinates within tolerance.
    EXACT = "exact"
    #: The source record lies wholly inside the finding.
    RECORD_WITHIN_FINDING = "record_within_finding"
    #: The finding lies wholly inside the source record.
    FINDING_WITHIN_RECORD = "finding_within_record"
    #: Partial overlap meeting the configured reciprocal threshold.
    OVERLAP = "overlap"


#: ClinVar's own review-status strings, mapped to the star rating NCBI publishes for them.
#: Kept as data rather than parsed heuristically: the strings are a closed set NCBI defines,
#: and an unrecognised one must surface as unknown rather than be scored by guesswork.
CLINVAR_REVIEW_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no assertion provided": 0,
    "no classifications from unflagged records": 0,
    "no classification for the single variant": 0,
}


def review_stars(review_status: str) -> int | None:
    """Star rating for a ClinVar review status, or ``None`` when the string is unknown.

    ``None`` rather than ``0``: an unrecognised status means this code has not kept up with
    NCBI's vocabulary, which is a different situation from a record NCBI itself rates as
    having no assertion criteria. Scoring the first as the second would hide a stale mapping
    behind a plausible number.
    """
    return CLINVAR_REVIEW_STARS.get(review_status.strip().lower())


def align(source_origin: Origin, assay_intent: Origin) -> tuple[ScopeAlignment, str]:
    """Compare an assertion's origin with what the assay was looking for."""
    if source_origin is Origin.UNKNOWN or assay_intent is Origin.UNKNOWN:
        missing = "the source record" if source_origin is Origin.UNKNOWN else "the assay"
        return ScopeAlignment.UNKNOWN, (
            f"{missing} did not declare whether this concerns germline or somatic "
            "variation, so the classification's scope cannot be checked against the "
            "question that was asked"
        )
    if source_origin is assay_intent:
        return ScopeAlignment.ALIGNED, f"both concern {source_origin.value} variation"
    return ScopeAlignment.MISMATCHED, (
        f"the record classifies {source_origin.value} variation while the assay is "
        f"investigating {assay_intent.value} findings; this may still matter as a "
        "secondary finding, but it does not answer the question that was asked"
    )


@dataclass(frozen=True)
class Interval:
    """A half-open genomic interval on one contig."""

    contig: str
    start: int
    end: int

    def __post_init__(self) -> None:
        # ClinVar writes -1 for records whose placement is unknown, and -1 > -2, so an
        # end-after-start check alone admits a record at chr7:-2--1. It would match nothing,
        # forever, while being counted as usable — a silent hole in the annotation coverage
        # rather than a loud failure.
        if self.start < 0:
            raise ValueError(
                f"genomic coordinates are not negative: {self.contig}:{self.start}-{self.end}"
            )
        if self.end <= self.start:
            raise ValueError(
                f"interval end must exceed start: {self.contig}:{self.start}-{self.end}"
            )

    @property
    def length(self) -> int:
        return self.end - self.start


def canonical_contig(name: str) -> str:
    """Strip a ``chr`` prefix so ``chr7`` and ``7`` compare equal.

    ClinVar publishes bare names and BAM headers usually carry the prefix. A mismatch here
    produces zero annotations everywhere with no error, which is the worst possible failure:
    it looks exactly like a sample with nothing known about it.
    """
    text = name.strip()
    lowered = text.lower()
    if lowered.startswith("chr"):
        text = text[3:]
    return text.upper() if text.lower() in {"x", "y", "mt", "m"} else text


def reciprocal_overlap(left: Interval, right: Interval) -> float:
    """Smaller of the two overlap fractions, which is the conservative reading.

    Using the larger fraction would let a 90 Mb finding "match" a 3 Mb record at 1.0 simply
    by containing it.
    """
    if canonical_contig(left.contig) != canonical_contig(right.contig):
        return 0.0
    shared = min(left.end, right.end) - max(left.start, right.start)
    if shared <= 0:
        return 0.0
    return min(shared / left.length, shared / right.length)


def classify_match(
    finding: Interval,
    record: Interval,
    *,
    minimum_reciprocal_overlap: float,
    exact_tolerance_bp: int,
) -> tuple[MatchType, float] | None:
    """Decide whether a record matches a finding, and how closely. ``None`` when it does not.

    ``exact_tolerance_bp`` exists because breakpoints from a read-depth caller and from a
    clinical array are never bit-identical, and demanding that they be would leave every real
    match classified as a mere overlap.
    """
    if not 0.0 < minimum_reciprocal_overlap <= 1.0:
        raise ValueError("minimum_reciprocal_overlap must lie in (0, 1]")
    if exact_tolerance_bp < 0:
        raise ValueError("exact_tolerance_bp must not be negative")

    overlap = reciprocal_overlap(finding, record)
    if overlap <= 0.0:
        return None

    if (
        abs(finding.start - record.start) <= exact_tolerance_bp
        and abs(finding.end - record.end) <= exact_tolerance_bp
    ):
        return MatchType.EXACT, overlap
    if record.start >= finding.start and record.end <= finding.end:
        return MatchType.RECORD_WITHIN_FINDING, overlap
    if finding.start >= record.start and finding.end <= record.end:
        return MatchType.FINDING_WITHIN_RECORD, overlap
    if overlap >= minimum_reciprocal_overlap:
        return MatchType.OVERLAP, overlap
    return None
