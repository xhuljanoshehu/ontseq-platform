"""Local-file boundary for the Lea historical CNV comparator.

This module is deliberately separate from :mod:`lea_compat`: parsing and scientific
normalization can be unit-tested without filesystem access, while this boundary adds the
checksums and locked reference resources needed for a reproducible local import.

It is runnable directly with ``python -m ontseq_platform.cnv.lea_io``.  The command writes
only a canonical research-only ``CnvCallSet`` JSON.  Source paths are never copied into the
reviewer-facing artifact.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from ..models import ReferenceLock
from .cytobands import load_cytoband_file
from .lea_compat import LEA_ACE_2026_HG19, lea_ace_call_set_from_outputs
from .models import CnvCallSet, CnvDataBasis


def sha256_file(path: Path) -> str:
    """Return the SHA-256 of a local file without exposing its path in result contracts."""
    if not path.is_file():
        raise ValueError(f"required historical artifact is missing or unreadable: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_lea_ace_call_set(
    *,
    cn_csv: Path,
    dels_dups_csv: Path,
    cytoband_file: Path,
    expected_cytoband_sha256: str,
    cytoband_resource_id: str,
    reference_lock: ReferenceLock,
    call_set_id: str,
    sample_id: str,
    data_basis: CnvDataBasis,
) -> CnvCallSet:
    """Load, fingerprint and normalize one historical QDNAseq+ACE result pair.

    The historical profile is frozen to GRCh37.  A reference lock for another build is
    refused; the compatibility lane never performs lift-over.  The cytoband file must match
    an expected SHA-256 supplied by deployment configuration, so a same-named but different
    resource cannot silently change band coordinates.
    """
    profile = LEA_ACE_2026_HG19
    profile.validate()
    if reference_lock.genome_build != profile.genome_build:
        raise ValueError(
            f"reference lock build {reference_lock.genome_build.value} does not match "
            f"historical comparator build {profile.genome_build.value}"
        )
    observed_cytoband_sha256 = sha256_file(cytoband_file)
    if observed_cytoband_sha256 != expected_cytoband_sha256.lower():
        raise ValueError(
            "cytoband resource SHA-256 does not match the deployment lock; historical "
            "compatibility import refused"
        )

    cn_sha256 = sha256_file(cn_csv)
    dels_dups_sha256 = sha256_file(dels_dups_csv)
    cytobands = load_cytoband_file(
        cytoband_file,
        genome_build=profile.genome_build,
        resource_id=cytoband_resource_id,
        source_sha256=observed_cytoband_sha256,
    )
    call_set = lea_ace_call_set_from_outputs(
        cn_lines=cn_csv.read_text(encoding="utf-8-sig").splitlines(),
        dels_dups_lines=dels_dups_csv.read_text(encoding="utf-8-sig").splitlines(),
        cytobands=cytobands,
        contig_lengths={item.name: item.length for item in reference_lock.contigs},
        call_set_id=call_set_id,
        sample_id=sample_id,
        genome_build=profile.genome_build,
        data_basis=data_basis,
    )

    if call_set.tool is None:  # pragma: no cover - enforced by lea_compat
        raise RuntimeError("historical comparator returned no tool provenance")
    tool = call_set.tool.model_copy(
        update={
            "parameters": {
                **call_set.tool.parameters,
                "cn_csv_sha256": cn_sha256,
                "dels_dups_csv_sha256": dels_dups_sha256,
                "cytoband_sha256": observed_cytoband_sha256,
                "reference_fai_sha256": reference_lock.source_fai_sha256,
                "reference_id": reference_lock.reference_id,
            }
        }
    )
    return call_set.model_copy(update={"tool": tool})


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m ontseq_platform.cnv.lea_io",
        description=(
            "Import one historical Lea QDNAseq+ACE result pair into the canonical "
            "research-only CNV call-set contract. No lift-over and no clinical release."
        ),
    )
    parser.add_argument("--cn-csv", type=Path, required=True)
    parser.add_argument("--dels-dups-csv", type=Path, required=True)
    parser.add_argument("--cytobands", type=Path, required=True)
    parser.add_argument("--cytoband-sha256", required=True)
    parser.add_argument("--cytoband-resource-id", required=True)
    parser.add_argument("--reference-lock", type=Path, required=True)
    parser.add_argument("--call-set-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument(
        "--data-basis",
        choices=[item.value for item in CnvDataBasis if item != CnvDataBasis.SIMULATED],
        required=True,
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    args = _parser().parse_args()
    lock = ReferenceLock.model_validate_json(args.reference_lock.read_text(encoding="utf-8"))
    call_set = load_lea_ace_call_set(
        cn_csv=args.cn_csv,
        dels_dups_csv=args.dels_dups_csv,
        cytoband_file=args.cytobands,
        expected_cytoband_sha256=args.cytoband_sha256,
        cytoband_resource_id=args.cytoband_resource_id,
        reference_lock=lock,
        call_set_id=args.call_set_id,
        sample_id=args.sample_id,
        data_basis=CnvDataBasis(args.data_basis),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(call_set.model_dump_json(indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
