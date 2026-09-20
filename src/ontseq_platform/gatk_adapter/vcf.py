"""Loss-conscious extraction of literal alleles; NOT left-normalization or adjudication."""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from pathlib import Path
from typing import TextIO

from .contracts import SmallVariantCandidate


class VcfContractError(ValueError):
    """Malformed or unsupported native evidence; never silently drop a record."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _numbers(value: str | None, length: int, *, integer: bool = False) -> tuple:
    if value is None or value == '.':
        return (None,) * length
    items = value.split(',')
    if len(items) != length:
        raise VcfContractError('FORMAT allele cardinality does not match the ALT field')
    parsed = []
    for item in items:
        if item == '.':
            parsed.append(None)
            continue
        if integer and not re.fullmatch(r'\d+', item):
            raise VcfContractError('Read depth must be a non-negative integer')
        try:
            number = int(item) if integer else float(item)
        except ValueError as exc:
            raise VcfContractError('Invalid numeric FORMAT value') from exc
        if not math.isfinite(number) or number < 0 or (not integer and number > 1):
            raise VcfContractError('AF must be finite in [0,1]; read counts must be non-negative')
        parsed.append(number)
    return tuple(parsed)


def _sample_fields(keys: list[str], value: str, n_alt: int) -> dict:
    if value == '.':
        values = ['.'] * len(keys)
    else:
        values = value.split(':')
        # VCF permits omitted trailing FORMAT values; do not invent zeroes.
        if len(values) > len(keys):
            raise VcfContractError('Too many sample FORMAT fields')
        values += ['.'] * (len(keys) - len(values))
    fields = dict(zip(keys, values, strict=True))
    return {
        'af': _numbers(fields.get('AF'), n_alt),
        'ad': _numbers(fields.get('AD'), n_alt + 1, integer=True),
        'dp': _numbers(fields.get('DP'), 1, integer=True)[0],
    }


def _variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt):
        return 'SNV' if len(ref) == 1 else 'MNV'
    if len(ref) == 1 and alt.startswith(ref):
        return 'INS'
    if len(alt) == 1 and ref.startswith(alt):
        return 'DEL'
    return 'COMPLEX'


def extract_candidates(
    path: Path, *, run_id: str, sample_id: str, normal_sample_id: str | None = None,
    max_records: int = 2_000_000,
) -> tuple[SmallVariantCandidate, ...]:
    """Retain every literal ALT and every FILTER state in a separately versioned schema.

    Coordinates describe the original VCF reference span (including the anchor).
    AF is the caller's FORMAT/AF, not AD/DP, purity, or a validated somatic probability.
    An empty tuple from a valid header is not a biological negative.
    """
    source_digest = sha256_file(path)
    opener = gzip.open if path.suffix == '.gz' else open
    expected = {sample_id} | ({normal_sample_id} if normal_sample_id else set())
    records: list[SmallVariantCandidate] = []
    header: list[str] | None = None
    version_seen = False
    handle: TextIO
    with opener(path, 'rt', encoding='utf-8', errors='strict') as handle:
        for line_number, raw in enumerate(handle, 1):
            if len(raw) > 4 * 1024 * 1024:
                raise VcfContractError('Native VCF line exceeds the bounded parser limit')
            line = raw.rstrip('\r\n')
            if line.startswith('##'):
                if header is not None:
                    raise VcfContractError('Metadata found after VCF header')
                if line.startswith('##fileformat=VCFv4.'):
                    version_seen = True
                continue
            if line.startswith('#CHROM\t'):
                if header is not None:
                    raise VcfContractError('Duplicate VCF header')
                header = line.split('\t')
                if header[:9] != [
                    '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
                ] or not version_seen:
                    raise VcfContractError('Unsupported or missing VCF header')
                if set(header[9:]) != expected or len(header[9:]) != len(expected):
                    raise VcfContractError('VCF sample names do not match the requested sample set')
                continue
            if not line:
                continue
            if header is None:
                raise VcfContractError('VCF data found without a valid header')
            columns = line.split('\t')
            if len(columns) != len(header):
                raise VcfContractError('VCF record column count does not match the header')
            chrom, pos_text, _, ref, alt_text, _, filters, info_text, format_text = columns[:9]
            if not re.fullmatch(r'[1-9]\d*', pos_text):
                raise VcfContractError('VCF POS must be a positive integer')
            pos = int(pos_text)
            alts = alt_text.split(',')
            if not chrom or len(set(alts)) != len(alts):
                raise VcfContractError('Invalid chromosome or duplicate ALT allele')
            if any(not re.fullmatch(r'[ACGTN]+', allele) for allele in [ref, *alts]):
                raise VcfContractError('Unsupported or symbolic allele; native VCF is retained')
            if ref in alts:
                raise VcfContractError('ALT allele cannot equal REF')
            if not filters:
                raise VcfContractError('Missing native FILTER value')
            info = {}
            if info_text != '.':
                for entry in info_text.split(';'):
                    key, separator, value = entry.partition('=')
                    if not key or key in info:
                        raise VcfContractError('Invalid or duplicated INFO field')
                    info[key] = value if separator else True
            keys = format_text.split(':')
            if len(keys) != len(set(keys)) or any(not key for key in keys):
                raise VcfContractError('Invalid FORMAT keys')
            tumor = _sample_fields(keys, columns[header.index(sample_id, 9)], len(alts))
            normal = (
                _sample_fields(keys, columns[header.index(normal_sample_id, 9)], len(alts))
                if normal_sample_id else None
            )
            record_digest = hashlib.sha256(line.encode()).hexdigest()
            for allele_index, alt in enumerate(alts, 1):
                if len(records) >= max_records:
                    raise VcfContractError('Candidate limit exceeded; no truncation is permitted')
                identity = [run_id, sample_id, source_digest, line_number, allele_index]
                allele_filters = info.get('AS_FilterStatus')
                if allele_filters is True:
                    raise VcfContractError('AS_FilterStatus must be a value, not a flag')
                records.append(SmallVariantCandidate(
                    evidence_id=hashlib.sha256(json.dumps(identity).encode()).hexdigest(),
                    run_id=run_id, sample_id=sample_id, chromosome=chrom,
                    start0=pos - 1, end0=pos - 1 + len(ref), native_pos1=pos,
                    ref=ref, alt=alt, allele_index=allele_index,
                    variant_type=_variant_type(ref, alt),
                    tumor_af=tumor['af'][allele_index - 1],
                    normal_af=normal['af'][allele_index - 1] if normal else None,
                    ref_depth=tumor['ad'][0], alt_depth=tumor['ad'][allele_index],
                    depth=tumor['dp'], native_site_filter=filters,
                    native_allele_filter_status=allele_filters,
                    native_record_sha256=record_digest, source_vcf_sha256=source_digest,
                    source_line_number=line_number,
                ))
    if header is None:
        raise VcfContractError('Missing VCF header')
    return tuple(records)
