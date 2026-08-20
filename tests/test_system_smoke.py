from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cnv.qdnaseq import (
    CnvChromosomeConsensus,
    CnvFit,
    QDNAseqCallReport,
    QDNAseqPolicy,
)
from ontseq_platform.models import (
    CheckStatus,
    EventType,
    Evidence,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    ToolRecord,
)
from ontseq_platform.system_smoke import (
    _verify_release_checksums,
    cnv_truth_checks,
    synthetic_cnv_copy_count,
)


def _fit(bin_size_kbp: int) -> CnvFit:
    return CnvFit(
        bin_size_kbp=bin_size_kbp,
        cellularity=1.0,
        ploidy=2.0,
        fit_error=0.01,
        candidate_count=2,
        segment_count=2,
        alternatives=[],
        segment_file=f"SYNTH.{bin_size_kbp}kbp.segments.tsv",
        chromosome_file=f"SYNTH.{bin_size_kbp}kbp.chromosomes.tsv",
        fit_plot=f"SYNTH.{bin_size_kbp}kbp.fit.png",
        copy_number_plot=f"SYNTH.{bin_size_kbp}kbp.cn.png",
        rds_file=f"SYNTH.{bin_size_kbp}kbp.rds",
    )


def _event(
    event_id: str,
    chromosome: str,
    event_type: EventType,
    copy_number: float,
) -> GenomicEvent:
    return GenomicEvent(
        event_id=event_id,
        event_type=event_type,
        primary=Locus(chromosome=chromosome, start=1, end=100_000_000),
        length_bp=99_999_999,
        copy_number=copy_number,
        evidence=[Evidence(caller="QDNAseq+ACE", caller_version="1.42.0/1.24.0")],
    )


def _report(
    *,
    include_gain8: bool = True,
    chr8_agreeing_bins: int = 3,
    primary_bin: int = 500,
) -> QDNAseqCallReport:
    fits = [_fit(100), _fit(500), _fit(1000)]
    events = [
        _event("CNV_SYNTH_0001", "chr7", EventType.CHROMOSOME_LOSS, 1.0),
    ]
    if include_gain8:
        events.append(
            _event("CNV_SYNTH_0002", "chr8", EventType.CHROMOSOME_GAIN, 3.0)
        )
    return QDNAseqCallReport(
        sample_id="SYNTH",
        genome_build=GenomeBuild.GRCH37,
        status=ModuleRunStatus.COMPLETED,
        primary_fit=next(item for item in fits if item.bin_size_kbp == primary_bin),
        fits=fits,
        chromosome_consensus=[
            CnvChromosomeConsensus(
                chromosome="chr7",
                median_copy_number=1.0,
                rounded_copy_number=1,
                agreeing_bins=3,
                contributing_bins=3,
                min_copy_number=1.0,
                max_copy_number=1.0,
            ),
            CnvChromosomeConsensus(
                chromosome="chr8",
                median_copy_number=3.0,
                rounded_copy_number=3,
                agreeing_bins=chr8_agreeing_bins,
                contributing_bins=3,
                min_copy_number=3.0,
                max_copy_number=3.0,
            ),
        ],
        events=events,
        tools=[
            ToolRecord(name="QDNAseq", version="1.42.0"),
            ToolRecord(name="ACE", version="1.24.0"),
        ],
        output_files=["SYNTH.qdnaseq-ace.summary.json"],
    )


class SystemSmokeTruthTests(unittest.TestCase):
    def setUp(self) -> None:
        self.policy = QDNAseqPolicy(
            profile_id="system-smoke-test",
            bin_sizes_kbp=[100, 500, 1000],
            primary_bin_size_kbp=500,
            note="test policy",
        )

    def test_depth_fixture_encodes_expected_chr7_and_chr8_truth(self) -> None:
        self.assertEqual(synthetic_cnv_copy_count(1), 2)
        self.assertEqual(synthetic_cnv_copy_count(7), 1)
        self.assertEqual(synthetic_cnv_copy_count(8), 3)
        self.assertEqual(synthetic_cnv_copy_count(22), 2)
        with self.assertRaises(ValueError):
            synthetic_cnv_copy_count(23)

    def test_expected_qdnaseq_truth_passes_every_check(self) -> None:
        checks = cnv_truth_checks(_report(), self.policy)
        self.assertGreaterEqual(len(checks), 6)
        self.assertTrue(all(item.status == CheckStatus.PASS for item in checks), checks)

    def test_missing_chr8_gain_is_detected(self) -> None:
        checks = {item.name: item for item in cnv_truth_checks(
            _report(include_gain8=False), self.policy
        )}
        self.assertEqual(checks["expected_chr8_gain"].status, CheckStatus.FAIL)
        self.assertEqual(checks["expected_chr7_loss"].status, CheckStatus.PASS)

    def test_multibin_disagreement_is_detected(self) -> None:
        checks = {item.name: item for item in cnv_truth_checks(
            _report(chr8_agreeing_bins=2), self.policy
        )}
        self.assertEqual(checks["cnv_multibin_consensus"].status, CheckStatus.FAIL)

    def test_wrong_primary_bin_is_detected(self) -> None:
        checks = {item.name: item for item in cnv_truth_checks(
            _report(primary_bin=100), self.policy
        )}
        self.assertEqual(checks["qdnaseq_primary_bin"].status, CheckStatus.FAIL)


class ReleaseChecksumVerificationTests(unittest.TestCase):
    def test_valid_manifest_passes_and_tampering_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            artifact = root / "reports" / "result.txt"
            artifact.parent.mkdir(parents=True)
            artifact.write_text("stable\n", encoding="utf-8")
            release = root / "release"
            release.mkdir()
            digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
            (release / "checksums.sha256").write_text(
                f"{digest}  reports/result.txt\n", encoding="utf-8"
            )

            self.assertEqual(_verify_release_checksums(root), (True, 1))
            artifact.write_text("changed\n", encoding="utf-8")
            self.assertEqual(_verify_release_checksums(root), (False, 0))

    def test_manifest_cannot_escape_envelope_root(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            parent = Path(temp)
            root = parent / "envelope"
            release = root / "release"
            release.mkdir(parents=True)
            outside = parent / "outside.txt"
            outside.write_text("outside\n", encoding="utf-8")
            digest = hashlib.sha256(outside.read_bytes()).hexdigest()
            (release / "checksums.sha256").write_text(
                f"{digest}  ../outside.txt\n", encoding="utf-8"
            )

            self.assertEqual(_verify_release_checksums(root), (False, 0))


if __name__ == "__main__":
    unittest.main()
