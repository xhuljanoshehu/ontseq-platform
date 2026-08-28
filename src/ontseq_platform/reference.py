from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from .models import (
    GenomeBuild,
    ReferenceContig,
    ReferenceDictionaryContract,
    ReferenceLock,
)

# Nuclear chromosome lengths are stable assembly identifiers. Requiring all 24 nuclear
# chromosomes lets the Desktop reject a tiny test or region-only dictionary without making
# assumptions about optional mitochondrial, decoy, ALT or random contigs.
_CANONICAL_NUCLEAR_LENGTHS: dict[GenomeBuild, tuple[int, ...]] = {
    GenomeBuild.GRCH37: (
        249250621,
        243199373,
        198022430,
        191154276,
        180915260,
        171115067,
        159138663,
        146364022,
        141213431,
        135534747,
        135006516,
        133851895,
        115169878,
        107349540,
        102531392,
        90354753,
        81195210,
        78077248,
        59128983,
        63025520,
        48129895,
        51304566,
        155270560,
        59373566,
    ),
    GenomeBuild.GRCH38: (
        248956422,
        242193529,
        198295559,
        190214555,
        181538259,
        170805979,
        159345973,
        145138636,
        138394717,
        133797422,
        135086622,
        133275309,
        114364328,
        107043718,
        101991189,
        90338345,
        83257441,
        80373285,
        58617616,
        64444167,
        46709983,
        50818468,
        156040895,
        57227415,
    ),
}


@dataclass(frozen=True)
class CanonicalReferenceSummary:
    genome_build: GenomeBuild
    naming_style: str
    contig_count: int
    total_reference_bases: int


def _canonical_profile(genome_build: GenomeBuild, prefix: str) -> dict[str, int]:
    labels = [*(str(number) for number in range(1, 23)), "X", "Y"]
    return {
        f"{prefix}{label}": length
        for label, length in zip(labels, _CANONICAL_NUCLEAR_LENGTHS[genome_build], strict=True)
    }


def canonical_contigs(
    genome_build: GenomeBuild, *, chr_prefix: bool = True
) -> tuple[tuple[str, int], ...]:
    """Return the ordered canonical nuclear dictionary used for build detection."""

    prefix = "chr" if chr_prefix else ""
    return tuple(_canonical_profile(genome_build, prefix).items())


def grch38_canonical_25_contigs() -> tuple[tuple[str, int], ...]:
    """Return chr1-22, chrX, chrY and chrM in the only supported 25-contig order."""

    return (*canonical_contigs(GenomeBuild.GRCH38), ("chrM", 16569))


def validate_grch38_canonical_25(
    contigs: Iterable[tuple[str, int]],
) -> CanonicalReferenceSummary:
    """Require exactly the chr-prefixed GRCh38 Canonical-25 dictionary.

    Build detection deliberately tolerates additional contigs.  A profile that names this
    contract does not: decoys, ALT loci, unplaced scaffolds, missing chrM and reordered
    dictionaries are all distinct alignment references and therefore fail closed.
    """

    observed = tuple(contigs)
    expected = grch38_canonical_25_contigs()
    if observed != expected:
        observed_map = dict(observed)
        expected_map = dict(expected)
        missing = [name for name in expected_map if name not in observed_map]
        extras = [name for name in observed_map if name not in expected_map]
        mismatched = [
            name
            for name, length in expected_map.items()
            if name in observed_map and observed_map[name] != length
        ]
        order_mismatch = not missing and not extras and not mismatched
        raise ValueError(
            "GRCh38 Canonical-25 dictionary must be exactly chr1-22, chrX, chrY, chrM "
            "with standard lengths and order: "
            f"{len(missing)} missing, {len(extras)} extra, "
            f"{len(mismatched)} length mismatches, order_mismatch={order_mismatch}"
        )
    return CanonicalReferenceSummary(
        genome_build=GenomeBuild.GRCH38,
        naming_style="chr-prefixed",
        contig_count=len(expected),
        total_reference_bases=sum(length for _, length in expected),
    )


def reference_lock_for_dictionary_contract(
    reference_lock: ReferenceLock,
    contract: ReferenceDictionaryContract,
) -> ReferenceLock:
    """Derive the exact run lock selected by a profile's explicit dictionary contract."""

    if contract == ReferenceDictionaryContract.EXACT_FULL:
        return reference_lock
    if contract != ReferenceDictionaryContract.GRCH38_CANONICAL_25:
        raise ValueError(f"unsupported reference dictionary contract: {contract.value}")
    if reference_lock.genome_build != GenomeBuild.GRCH38:
        raise ValueError("grch38_canonical_25 cannot be derived from a non-GRCh38 ReferenceLock")

    expected = grch38_canonical_25_contigs()
    source = tuple((item.name, item.length) for item in reference_lock.contigs)
    source_positions = {record: index for index, record in enumerate(source)}
    try:
        positions = [source_positions[record] for record in expected]
    except KeyError as exc:
        raise ValueError(
            "pinned GRCh38 ReferenceLock does not contain the complete Canonical-25 dictionary"
        ) from exc
    if positions != sorted(positions):
        raise ValueError("pinned GRCh38 ReferenceLock orders Canonical-25 contigs incompatibly")

    validate_grch38_canonical_25(expected)
    return ReferenceLock(
        reference_id=reference_lock.reference_id,
        genome_build=reference_lock.genome_build,
        contigs=[ReferenceContig(name=name, length=length) for name, length in expected],
        allow_extra_contigs=False,
        source_fai_sha256=reference_lock.source_fai_sha256,
    )


def validate_canonical_reference(
    contigs: Iterable[tuple[str, int]], genome_build: GenomeBuild
) -> CanonicalReferenceSummary:
    """Require a complete GRCh37/GRCh38 nuclear dictionary for a named Desktop build.

    Both common sequence-name styles are accepted (``chr1`` and ``1``). Additional contigs
    are deliberately ignored here: their exact compatibility remains the BAM intake gate's
    job. A partial, wrong-build or mixed-style dictionary cannot masquerade as a complete
    named build. This checks the FAI dictionary, not the FASTA base content.
    """

    materialized = tuple(contigs)
    observed = dict(materialized)
    if len(observed) != len(materialized):
        raise ValueError("reference contains duplicate contig names")

    candidates: list[tuple[int, int, str, list[str], list[str]]] = []
    matches: list[tuple[str, str]] = []
    for prefix, style in (("chr", "chr-prefixed"), ("", "unprefixed")):
        expected = _canonical_profile(genome_build, prefix)
        missing = [name for name in expected if name not in observed]
        mismatched = [
            name
            for name, length in expected.items()
            if name in observed and observed[name] != length
        ]
        candidates.append(
            (len(missing) + len(mismatched), len(missing), style, missing, mismatched)
        )
        if not missing and not mismatched:
            matches.append((prefix, style))

    if matches:
        prefix, style = matches[0]
        opposite_prefix = "" if prefix else "chr"
        opposite_alias_count = sum(
            name in observed for name in _canonical_profile(genome_build, opposite_prefix)
        )
        if opposite_alias_count:
            raise ValueError(
                f"{genome_build.value} canonical assembly validation failed: the dictionary "
                f"mixes canonical naming styles ({opposite_alias_count} opposite-style "
                "chromosome aliases found). Use one consistent 1-22, X and Y naming style."
            )
        return CanonicalReferenceSummary(
            genome_build=genome_build,
            naming_style=style,
            contig_count=len(materialized),
            total_reference_bases=sum(length for _, length in materialized),
        )

    _, _, style, missing, mismatched = min(candidates)
    examples = [*missing[:3], *mismatched[:3]]
    example_text = f" Examples: {', '.join(examples)}." if examples else ""
    raise ValueError(
        f"{genome_build.value} canonical assembly validation failed: expected complete "
        "chromosomes 1-22, X and Y with standard lengths in one naming style; "
        f"found {len(materialized)} contigs, with {len(missing)} missing and "
        f"{len(mismatched)} length mismatches in the closest ({style}) profile."
        f"{example_text} Partial or wrong-build dictionaries cannot be configured as a "
        f"full {genome_build.value} Desktop reference."
    )


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def contig_signature(contigs: Iterable[tuple[str, int]]) -> str:
    digest = hashlib.sha256()
    for name, length in contigs:
        digest.update(f"{name}\t{length}\n".encode())
    return digest.hexdigest()


def reference_lock_signature(reference_lock: ReferenceLock) -> str:
    """Fingerprint every lock field using deterministic JSON for resume decisions."""

    canonical = json.dumps(
        reference_lock.model_dump(mode="json"),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return hashlib.sha256(canonical).hexdigest()


def reference_lock_from_fai(
    fai_path: Path,
    *,
    reference_id: str,
    genome_build: GenomeBuild,
    allow_extra_contigs: bool = False,
    require_canonical_assembly: bool = False,
) -> ReferenceLock:
    contigs: list[ReferenceContig] = []
    with fai_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Invalid FASTA index line {line_number}: expected at least 2 fields"
                )
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid FASTA index length on line {line_number}: {fields[1]!r}"
                ) from exc
            contigs.append(ReferenceContig(name=fields[0], length=length))
    if not contigs:
        raise ValueError(f"FASTA index contains no contigs: {fai_path}")
    lock = ReferenceLock(
        reference_id=reference_id,
        genome_build=genome_build,
        contigs=contigs,
        allow_extra_contigs=allow_extra_contigs,
        source_fai_sha256=sha256_file(fai_path),
    )
    if require_canonical_assembly:
        validate_canonical_reference(
            ((item.name, item.length) for item in lock.contigs), lock.genome_build
        )
    return lock
