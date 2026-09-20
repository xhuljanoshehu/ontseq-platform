"""Narrow Mutect2 VCF reader. Native alleles/filters are evidence, not clinical calls."""

from __future__ import annotations

import gzip
import math
import re
from pathlib import Path
from typing import Any, TextIO


def _number(value: str, *, integer: bool = False) -> int | float | None:
    if value == ".":
        return None
    parsed = int(value) if integer else float(value)
    if not math.isfinite(parsed) or parsed < 0:
        raise ValueError("VCF contains a negative or non-finite number")
    return parsed


def read_candidates(
    path: Path, tumor_name: str, expected_samples: tuple[str, ...]
) -> list[dict[str, Any]]:
    """Preserve every literal SNV/indel/MNV allele, including failed FILTER records.

    Positions are native 1-based VCF positions. Alleles are NOT left-normalized, so
    this output must not be used as an inter-caller join key without normalization.
    AF is the caller estimate; AD fraction uses sum(AD), NOT FORMAT/DP. Missing
    values remain null. An absent/malformed sample or record fails the whole import.
    """
    stream: TextIO = gzip.open(path, "rt") if path.suffix == ".gz" else path.open()
    rows: list[dict[str, Any]] = []
    samples: list[str] | None = None
    saw_fileformat = False
    with stream:
        for line in stream:
            line = line.rstrip("\r\n")
            if line.startswith("##fileformat=VCFv4."):
                saw_fileformat = True
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM\t"):
                if samples is not None:
                    raise ValueError("Duplicate VCF column header")
                columns = line.split("\t")
                if columns[:9] != [
                    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
                ]:
                    raise ValueError("Invalid VCF column header")
                samples = columns[9:]
                if (
                    len(samples) != len(set(samples))
                    or set(samples) != set(expected_samples)
                    or tumor_name not in samples
                ):
                    raise ValueError("VCF sample identity does not match the locked inputs")
                continue
            if not line or line.startswith("#"):
                raise ValueError("Unexpected blank/comment VCF line")
            if not saw_fileformat or samples is None:
                raise ValueError("VCF data before required headers")
            fields = line.split("\t")
            if len(fields) != 9 + len(samples):
                raise ValueError("Wrong number of VCF columns")
            chrom, pos, _, ref, alternate, qual, filters, info, fmt = fields[:9]
            alts = alternate.split(",")
            position = int(pos)
            if (
                not chrom or position < 1 or not filters
                or not re.fullmatch(r"[ACGTNacgtn]+", ref)
                or len(alts) != len(set(alts))
                or any(not re.fullmatch(r"[ACGTNacgtn]+", alt) or alt == ref for alt in alts)
            ):
                raise ValueError("Invalid or unsupported literal small variant")
            info_values = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
            allele_filters = (
                info_values["AS_FilterStatus"].split("|") if "AS_FilterStatus" in info_values else None
            )
            if allele_filters is not None and len(allele_filters) != len(alts):
                raise ValueError("AS_FilterStatus must describe every ALT")
            keys = fmt.split(":")
            values = fields[9 + samples.index(tumor_name)].split(":")
            if len(keys) != len(set(keys)) or len(keys) != len(values):
                raise ValueError("Invalid VCF FORMAT values")
            sample = dict(zip(keys, values, strict=True))
            ad_text = sample.get("AD", ".")
            ad = None if ad_text == "." else [_number(x, integer=True) for x in ad_text.split(",")]
            af_text = sample.get("AF", ".")
            af = None if af_text == "." else [_number(x) for x in af_text.split(",")]
            if ad is not None and (len(ad) != len(alts) + 1 or None in ad):
                raise ValueError("AD must have one count for REF and every ALT")
            if af is not None and (len(af) != len(alts) or any(x is not None and x > 1 for x in af)):
                raise ValueError("AF must have one valid fraction per ALT")
            total = sum(x for x in ad if x is not None) if ad is not None else None
            depth = _number(sample.get("DP", "."), integer=True)
            quality = _number(qual)
            for index, alt in enumerate(alts, start=1):
                alt_count = ad[index] if ad is not None else None
                allele_filter = allele_filters[index - 1] if allele_filters is not None else None
                technical_pass = filters == "PASS" and (
                    allele_filter in ("SITE", "PASS")
                    if allele_filters is not None else len(alts) == 1
                )
                rows.append({
                    "chromosome": chrom, "position": position, "reference": ref,
                    "alternate": alt, "native_alt_index": index, "quality": quality,
                    "depth": depth, "reference_reads": ad[0] if ad is not None else None,
                    "alternate_reads": alt_count,
                    "ad_fraction": alt_count / total if alt_count is not None and total else None,
                    "caller_af": af[index - 1] if af is not None else None,
                    "native_filter": filters, "native_info": info,
                    "technical_pass": technical_pass, "allele_filter": allele_filter,
                    "somatic_status": "NOT_VALIDATED", "reportable": False,
                    "requires_expert_review": True,
                    "requires_orthogonal_confirmation": len(ref) != len(alt),
                    "representation": "native_vcf_not_left_normalized",
                })
    if not saw_fileformat or samples is None:
        raise ValueError("Missing required VCF headers")
    return rows
