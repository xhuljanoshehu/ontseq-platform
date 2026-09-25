"""Synthetic regressions for GATK 4.6.2.0 AS_FilterStatus allele mapping."""

from __future__ import annotations

from pathlib import Path

import pytest

from ontseq_platform.gatk_adapter.vcf import VcfContractError, extract_candidates


def _vcf(tmp_path: Path, info: str, *, alts: str = "A,T") -> Path:
    path = tmp_path / "synthetic.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=Mutect2\n"
        "##contig=<ID=chrSynthetic,length=8>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_T\n"
        f"chrSynthetic\t3\t.\tG\t{alts}\t.\tPASS\t{info}\tGT:AD:AF:DP\t1/2:70,20,10:0.2,0.1:100\n"
    )
    return path


def _parse(path: Path):
    return extract_candidates(path, run_id="synthetic-run-001", sample_id="SYNTHETIC_T")


def test_as_filter_status_is_decoded_per_alt(tmp_path):
    records = _parse(_vcf(tmp_path, "AS_FilterStatus=SITE|weak_evidence"))
    assert [record.alt for record in records] == ["A", "T"]
    assert [record.native_allele_filter_status for record in records] == [
        "SITE",
        "weak_evidence",
    ]


def test_as_filter_status_preserves_multiple_filters_within_one_allele(tmp_path):
    records = _parse(_vcf(tmp_path, "AS_FilterStatus=clustered_events,strand_bias|SITE"))
    assert [record.native_allele_filter_status for record in records] == [
        "clustered_events,strand_bias",
        "SITE",
    ]


def test_as_filter_status_cardinality_mismatch_fails_closed(tmp_path):
    with pytest.raises(VcfContractError, match="AS_FilterStatus allele cardinality"):
        _parse(_vcf(tmp_path, "AS_FilterStatus=SITE"))


def test_as_filter_status_flag_fails_closed(tmp_path):
    with pytest.raises(VcfContractError, match="must be a value"):
        _parse(_vcf(tmp_path, "AS_FilterStatus"))


def test_missing_as_filter_status_remains_unknown(tmp_path):
    records = _parse(_vcf(tmp_path, "."))
    assert [record.native_allele_filter_status for record in records] == [None, None]
