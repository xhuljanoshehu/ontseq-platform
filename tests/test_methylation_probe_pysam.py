"""Real direct-BAM interoperability, using only freshly generated synthetic fixtures."""

from __future__ import annotations

import array
import os
import threading
from importlib import import_module
from pathlib import Path
from typing import Any
from unittest import SkipTest

import pytest

from ontseq_platform._methylation_worker import _indexed_offsets
from ontseq_platform.methylation_probe import probe_bam_methylation

try:
    pysam = import_module("pysam")
except ModuleNotFoundError as error:
    if error.name != "pysam":
        raise
    raise SkipTest("Direct BAM interoperability requires pysam in the analysis runtime") from error


def _record(
    position: int,
    *,
    mm: str | None = None,
    ml: list[int] | None = None,
    length: int = 20,
    reference_id: int = 0,
) -> Any:
    record = pysam.AlignedSegment()
    record.query_name = "SYNTHETIC_PRIVATE_READ"
    record.query_sequence = "CG" * (length // 2)
    record.flag = 0
    record.reference_id = reference_id
    record.reference_start = position
    record.mapping_quality = 60
    record.cigarstring = f"{length}M"
    record.query_qualities = array.array("B", [30] * length)
    if mm is not None:
        record.set_tag("MM", mm)
    if ml is not None:
        record.set_tag("ML", array.array("B", ml))
    return record


def _bam(tmp_path: Path, records: list[Any], *, indexed: str | None = None) -> Path:
    path = tmp_path / "SYNTHETIC.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "synthetic1", "LN": 20_000_000}, {"SN": "synthetic2", "LN": 20_000_000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as target:
        for record in records:
            target.write(record)
    if indexed == "bai":
        pysam.index(str(path))
    elif indexed == "csi":
        pysam.index("-c", str(path))
    return path


@pytest.mark.parametrize(
    ("mm", "ml", "reason"),
    [
        (None, None, "complete_no_tags"),
        ("C+m?,0;", [240], "detected_5mc"),
        ("C+m?,0;", None, "invalid_tags"),
        (None, [240], "invalid_tags"),
        ("C+m?,99;", [240], "invalid_tags"),
        ("C+m?,0;", [240, 240], "invalid_tags"),
        ("C+h?,0;", [240], "unsupported_modification"),
    ],
)
def test_supported_pairs_and_invalid_or_other_modifications(
    tmp_path: Path,
    mm: str | None,
    ml: list[int] | None,
    reason: str,
) -> None:
    path = _bam(tmp_path, [_record(0, mm=mm, ml=ml)])
    result = probe_bam_methylation(path, mode="thorough", timeout_seconds=15)
    assert result.reason_code == reason
    assert "SYNTHETIC_PRIVATE_READ" not in repr(result)
    if reason == "detected_5mc":
        assert result.confirmation_offset is not None
        assert result.detected_modifications == ("C+m:5mC",)
        with pysam.AlignmentFile(str(path), "rb") as reader:
            reader.seek(result.confirmation_offset)
            assert next(reader).has_tag("MM")


@pytest.mark.parametrize("indexed", ["bai", "csi"])
def test_indexed_sampling_finds_distant_tags_and_yields_exact_witness(
    tmp_path: Path,
    indexed: str,
) -> None:
    records = [_record(index * 30) for index in range(1200)]
    records.append(_record(10_000_000, mm="C+m?,0;", ml=[240], reference_id=1))
    path = _bam(tmp_path, records, indexed=indexed)
    assert len(_indexed_offsets(path)) > 1
    result = probe_bam_methylation(path, maximum_reads=128, timeout_seconds=15)
    assert result.status == "detected" and result.checked_reads <= 128
    assert result.confirmation_offset is not None
    confirmed = probe_bam_methylation(
        path,
        maximum_reads=1,
        confirmation_offset=result.confirmation_offset,
        timeout_seconds=15,
    )
    assert confirmed.status == "detected" and confirmed.checked_reads == 1


def test_partial_sample_never_asserts_absence_but_thorough_ignores_read_limit(
    tmp_path: Path,
) -> None:
    path = _bam(tmp_path, [_record(index * 30) for index in range(20)])
    quick = probe_bam_methylation(path, maximum_reads=2, timeout_seconds=15)
    assert quick.reason_code == "sample_incomplete" and not quick.complete
    thorough = probe_bam_methylation(path, maximum_reads=2, mode="thorough", timeout_seconds=15)
    assert thorough.status == "not_detected" and thorough.complete and thorough.checked_reads == 20


def test_indexed_negative_stays_unknown(tmp_path: Path) -> None:
    path = _bam(tmp_path, [_record(index * 30) for index in range(20)], indexed="bai")
    result = probe_bam_methylation(path, timeout_seconds=15)
    assert result.reason_code == "sample_incomplete" and not result.complete


def test_ultralong_record_avoids_former_sam_line_limit(tmp_path: Path) -> None:
    length = 200_000
    path = _bam(tmp_path, [_record(0, mm="C+m?,0;", ml=[240], length=length)])
    result = probe_bam_methylation(path, timeout_seconds=15)
    assert result.status == "detected"


def test_truncated_bam_never_returns_negative(tmp_path: Path) -> None:
    path = _bam(tmp_path, [_record(0)])
    path.write_bytes(path.read_bytes()[:-28])
    result = probe_bam_methylation(path, mode="thorough", timeout_seconds=15)
    assert result.status == "unknown" and not result.complete


def test_cancel_and_file_change_are_observed_during_real_stream(tmp_path: Path) -> None:
    path = _bam(tmp_path, [_record(index * 30) for index in range(500)])
    cancelled = threading.Event()
    result = probe_bam_methylation(
        path,
        mode="thorough",
        timeout_seconds=15,
        cancel_event=cancelled,
        progress_callback=lambda _progress: cancelled.set(),
    )
    assert result.reason_code == "cancelled"
    before = path.stat()
    changed = probe_bam_methylation(
        path,
        mode="thorough",
        timeout_seconds=15,
        progress_callback=lambda _progress: os.utime(
            path,
            ns=(before.st_atime_ns, before.st_mtime_ns + 1000000),
        ),
    )
    assert changed.reason_code == "file_changed" and changed.bam_identity is None


def test_stale_witness_does_not_preserve_positive(tmp_path: Path) -> None:
    path = _bam(tmp_path, [_record(0, mm="C+m?,0;", ml=[240])])
    first = probe_bam_methylation(path, timeout_seconds=15)
    path = _bam(tmp_path, [_record(0)])
    second = probe_bam_methylation(
        path,
        confirmation_offset=first.confirmation_offset,
        timeout_seconds=15,
    )
    assert second.status == "not_detected" and second.complete


def test_legacy_aliases_are_supported(tmp_path: Path) -> None:
    record = _record(0)
    record.set_tag("Mm", "C+m?,0;")
    record.set_tag("Ml", array.array("B", [240]))
    path = _bam(tmp_path, [record])
    assert probe_bam_methylation(path, timeout_seconds=15).status == "detected"


@pytest.mark.parametrize("case", ["duplicate_mm", "wrong_ml_type", "wrong_mn", "numeric_5mc"])
def test_modification_contract_rejects_ambiguous_or_unsupported_tags(
    tmp_path: Path,
    case: str,
) -> None:
    record = _record(0, mm="C+m?,0;", ml=[240])
    expected = "invalid_tags"
    if case == "duplicate_mm":
        record.set_tag("MM", "C+m?,0;", replace=False)
    elif case == "wrong_ml_type":
        record.set_tag("ML", array.array("H", [240]))
    elif case == "wrong_mn":
        record.set_tag("MN", 200)
    else:
        record.set_tag("MM", "C+27551?,0;")
        expected = "unsupported_modification"
    path = _bam(tmp_path, [record])
    result = probe_bam_methylation(path, mode="thorough", timeout_seconds=15)
    assert result.status == "unknown" and result.reason_code == expected


def test_record_size_guard_is_independent_of_sam_line_length(tmp_path: Path) -> None:
    path = _bam(tmp_path, [_record(0, length=4_000_002)])
    result = probe_bam_methylation(path, timeout_seconds=15)
    assert result.status == "unknown" and result.reason_code == "record_limit"


def test_malformed_index_falls_back_to_exhaustive_bam_stream(tmp_path: Path) -> None:
    path = _bam(tmp_path, [_record(0)])
    Path(str(path) + ".bai").write_bytes(b"BAI\x01\xff\xff\xff\x7f")
    result = probe_bam_methylation(path, timeout_seconds=15)
    assert result.status == "not_detected" and result.complete


def test_orphaned_mn_tag_is_incomplete_modified_base_evidence(tmp_path: Path) -> None:
    record = _record(0)
    record.set_tag("MN", 20)
    path = _bam(tmp_path, [record])
    result = probe_bam_methylation(path, mode="thorough", timeout_seconds=15)
    assert result.reason_code == "invalid_tags" and result.status == "unknown"


def test_corrupt_compressed_index_uses_safe_sequential_fallback(tmp_path: Path) -> None:
    path = _bam(tmp_path, [_record(0)])
    # Valid gzip header followed by an invalid DEFLATE block; no external fixture.
    Path(str(path) + ".csi").write_bytes(b"\x1f\x8b\x08\x00" + b"\x00" * 6 + b"\x07" * 20)
    result = probe_bam_methylation(path, timeout_seconds=15)
    assert result.status == "not_detected" and result.complete
