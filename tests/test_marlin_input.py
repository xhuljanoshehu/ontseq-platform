from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from ontseq_platform.marlin_contracts import MarlinSourceKind
from ontseq_platform.marlin_input import parse_marlin_probe_bed
from ontseq_platform.models import GenomeBuild


def test_parse_probe_bed_gzip_preserves_na_and_fingerprint(tmp_path: Path) -> None:
    path = tmp_path / "sample.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("chr1\t10\t11\t0.75\tcg0001\n")
        handle.write("chr1\t20\t21\tNA\tcg0002\n")
    parsed = parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)
    assert parsed.source_kind is MarlinSourceKind.PRECOMPUTED_METHYLATION
    assert parsed.genome_build is GenomeBuild.GRCH37
    assert parsed.adapter_id == "marlin_probe_bed_v1"
    assert parsed.input_fingerprint.size_bytes == path.stat().st_size
    assert parsed.input_fingerprint.sha256 is not None
    assert [item.probe_id for item in parsed.observations] == ["cg0001", "cg0002"]
    assert parsed.observations[0].methylation_fraction == 0.75
    assert parsed.observations[1].methylation_fraction is None


def test_parse_plain_text_uses_exact_tab_schema(tmp_path: Path) -> None:
    path = tmp_path / "sample.bed"
    path.write_text("chr2\t0\t1\t0.5\tcgA\n", encoding="utf-8")
    parsed = parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)
    assert parsed.observations[0].start == 0
    assert parsed.observations[0].methylation_fraction == 0.5


@pytest.mark.parametrize(
    "payload, message",
    [
        ("chr1\t1\t2\t0.5\tcg1\textra\n", "exactly 5"),
        ("chr1 1 2 0.5 cg1\n", "exactly 5"),
        ("chr1\t-1\t2\t0.5\tcg1\n", "start"),
        ("chr1\t2\t2\t0.5\tcg1\n", "greater than start"),
        ("chr1\t3\t2\t0.5\tcg1\n", "greater than start"),
        ("chr1\t1\t2\tnan\tcg1\n", "finite"),
        ("chr1\t1\t2\tinf\tcg1\n", "finite"),
        ("chr1\t1\t2\t-0.01\tcg1\n", "between 0 and 1"),
        ("chr1\t1\t2\t1.01\tcg1\n", "between 0 and 1"),
        ("chr1\t1\t2\t0.5\t\n", "probe"),
        ("chrM\t1\t2\t0.5\tcg1\n", "chromosome"),
        ("chr1\t1\t2\t0.5\tcg1\n\n", "blank"),
        ("chr1\t 1\t2\t0.5\tcg1\n", "whitespace"),
    ],
)
def test_malformed_probe_rows_fail_closed(tmp_path: Path, payload: str, message: str) -> None:
    path = tmp_path / "bad.bed"
    path.write_text(payload, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)


def test_duplicate_probe_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "dup.bed"
    path.write_text(
        "chr1\t10\t11\t0.1\tcg0001\nchr1\t12\t13\t0.9\tcg0001\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Duplicate MARLIN probe"):
        parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)


def test_empty_file_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "empty.bed"
    path.write_bytes(b"")
    with pytest.raises(ValueError, match="no MARLIN probe rows"):
        parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)


def test_grch38_is_rejected_for_v1(tmp_path: Path) -> None:
    path = tmp_path / "sample.bed"
    path.write_text("chr1\t10\t11\t0.5\tcg1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="GRCh37"):
        parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH38)


def test_input_is_fingerprinted_before_reading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "sample.bed"
    path.write_text("chr1\t10\t11\t0.5\tcg1\n", encoding="utf-8")
    calls: list[str] = []

    def fake_sha256(_path: Path) -> str:
        calls.append("fingerprint")
        return "a" * 64

    def exploding_open(*args: object, **kwargs: object) -> object:
        assert calls == ["fingerprint"]
        raise RuntimeError("stop-after-order-check")

    monkeypatch.setattr("ontseq_platform.marlin_input.sha256_file", fake_sha256)
    monkeypatch.setattr(Path, "open", exploding_open)
    with pytest.raises(RuntimeError, match="stop-after-order-check"):
        parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)
