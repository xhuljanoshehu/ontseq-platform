from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.clair3_vcf import (
    MALFORMED_INCONSISTENT_COUNTS,
    MALFORMED_MISSING_ALLELE_DEPTH,
    MALFORMED_MISSING_DEPTH,
    MALFORMED_MULTI_ALLELIC,
    MALFORMED_NO_FORMAT,
    MALFORMED_UNPARSABLE_COUNTS,
    MALFORMED_ZERO_DEPTH,
    Clair3ReadError,
    PreconditionFailure,
    check_preconditions,
    clair3_version,
    read_clair3_vcf,
)
from ontseq_platform.execution import CommandResult
from ontseq_platform.small_variants import INSERTION, Clair3Policy

HEADER = "\n".join(
    (
        "##fileformat=VCFv4.2",
        "##source=Clair3",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    )
)


def _record(
    *,
    chrom: str = "chr5",
    position: str = "171390000",
    reference: str = "A",
    alternate: str = "G",
    quality: str = "30.0",
    filter_status: str = "PASS",
    fmt: str = "GT:GQ:DP:AD:AF",
    sample: str = "0/1:30:80:50,30:0.375",
) -> str:
    return "\t".join(
        (chrom, position, ".", reference, alternate, quality, filter_status, ".", fmt, sample)
    )


class _StubRunner:
    def __init__(self, *, returncode: int = 0, stdout: str = "", stderr: str = "") -> None:
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv, *, timeout_seconds: int = 300) -> CommandResult:  # type: ignore[no-untyped-def]
        self.calls.append(tuple(argv))
        return CommandResult(
            argv=tuple(argv),
            returncode=self.returncode,
            stdout=self.stdout,
            stderr=self.stderr,
        )


class ReadVcfTests(unittest.TestCase):
    def _write(self, body: str, *, gzipped: bool = False) -> Path:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        suffix = ".vcf.gz" if gzipped else ".vcf"
        path = Path(directory.name) / f"calls{suffix}"
        text = f"{HEADER}\n{body}\n" if body else f"{HEADER}\n"
        if gzipped:
            path.write_bytes(gzip.compress(text.encode("utf-8")))
        else:
            path.write_text(text, encoding="utf-8")
        return path

    def test_reads_depth_and_alternate_support(self) -> None:
        contents = read_clair3_vcf(self._write(_record()))
        self.assertEqual(len(contents.variants), 1)
        variant = contents.variants[0]
        self.assertEqual(variant.depth, 80)
        self.assertEqual(variant.variant_reads, 30)
        self.assertAlmostEqual(variant.vaf, 0.375)
        self.assertEqual(variant.filter_status, "PASS")

    def test_reads_a_gzipped_vcf(self) -> None:
        contents = read_clair3_vcf(self._write(_record(), gzipped=True))
        self.assertEqual(len(contents.variants), 1)

    def test_alternate_support_is_taken_not_derived(self) -> None:
        """At a site with a third allele, depth minus reference overstates support.

        AD here is ref=50, alt=20, and 10 reads carry something else. Deriving support as
        80 - 50 would report 30.
        """
        contents = read_clair3_vcf(
            self._write(_record(sample="0/1:30:80:50,20,10:0.25", fmt="GT:GQ:DP:AD:AF"))
        )
        self.assertEqual(contents.variants[0].variant_reads, 20)

    def test_an_npm1_style_insertion_survives_the_round_trip(self) -> None:
        contents = read_clair3_vcf(self._write(_record(reference="C", alternate="CTCTG")))
        variant = contents.variants[0]
        self.assertEqual(variant.variant_class, INSERTION)
        self.assertEqual(variant.length_change, 4)

    def test_a_missing_quality_is_read_as_zero_not_dropped(self) -> None:
        contents = read_clair3_vcf(self._write(_record(quality=".")))
        self.assertEqual(contents.variants[0].quality, 0.0)

    def test_non_pass_records_are_kept_for_policy_to_judge(self) -> None:
        """Reading is not filtering. The policy decides what PASS means, not the reader."""
        contents = read_clair3_vcf(self._write(_record(filter_status="LowQual")))
        self.assertEqual(len(contents.variants), 1)
        self.assertEqual(contents.variants[0].filter_status, "LowQual")

    def test_an_empty_vcf_yields_nothing_without_error(self) -> None:
        contents = read_clair3_vcf(self._write(""))
        self.assertEqual(contents.variants, ())
        self.assertEqual(contents.malformed, ())
        self.assertEqual(contents.total_records, 0)

    def test_an_unreadable_path_is_an_error(self) -> None:
        with self.assertRaises(Clair3ReadError):
            read_clair3_vcf(Path("/nonexistent/calls.vcf"))


class MalformedRecordTests(unittest.TestCase):
    def _write(self, body: str) -> Path:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "calls.vcf"
        path.write_text(f"{HEADER}\n{body}\n", encoding="utf-8")
        return path

    def _reason(self, body: str) -> str:
        contents = read_clair3_vcf(self._write(body))
        self.assertEqual(contents.variants, ())
        self.assertEqual(len(contents.malformed), 1)
        return contents.malformed[0][1]

    def test_a_record_without_a_sample_column(self) -> None:
        self.assertEqual(
            self._reason("chr5\t171390000\t.\tA\tG\t30.0\tPASS\t."), MALFORMED_NO_FORMAT
        )

    def test_a_record_without_dp(self) -> None:
        self.assertEqual(
            self._reason(_record(fmt="GT:GQ:AD", sample="0/1:30:50,30")),
            MALFORMED_MISSING_DEPTH,
        )

    def test_a_record_without_ad(self) -> None:
        self.assertEqual(
            self._reason(_record(fmt="GT:GQ:DP:AF", sample="0/1:30:80:0.375")),
            MALFORMED_MISSING_ALLELE_DEPTH,
        )

    def test_a_multi_allelic_record_is_refused_rather_than_split(self) -> None:
        """One AD pair cannot be attributed to two alternates without guessing."""
        self.assertEqual(self._reason(_record(alternate="G,T")), MALFORMED_MULTI_ALLELIC)

    def test_non_numeric_counts(self) -> None:
        self.assertEqual(
            self._reason(_record(sample="0/1:30:deep:50,30:0.375")), MALFORMED_UNPARSABLE_COUNTS
        )

    def test_an_ad_without_an_alternate_count(self) -> None:
        self.assertEqual(
            self._reason(_record(sample="0/1:30:80:50:0.375")), MALFORMED_UNPARSABLE_COUNTS
        )

    def test_zero_depth(self) -> None:
        self.assertEqual(self._reason(_record(sample="0/1:30:0:0,0:0.0")), MALFORMED_ZERO_DEPTH)

    def test_support_exceeding_depth(self) -> None:
        self.assertEqual(
            self._reason(_record(sample="0/1:30:10:5,40:0.9")), MALFORMED_INCONSISTENT_COUNTS
        )

    def test_malformed_records_are_counted_by_reason_not_discarded(self) -> None:
        body = "\n".join(
            (
                _record(),
                _record(alternate="G,T"),
                _record(alternate="G,C"),
                _record(fmt="GT:GQ:AD", sample="0/1:30:50,30"),
            )
        )
        contents = read_clair3_vcf(self._write(body))
        self.assertEqual(len(contents.variants), 1)
        self.assertEqual(contents.total_records, 4)
        self.assertEqual(
            contents.malformed_counts,
            {MALFORMED_MULTI_ALLELIC: 2, MALFORMED_MISSING_DEPTH: 1},
        )


class VersionParsingTests(unittest.TestCase):
    def test_parses_a_bare_version(self) -> None:
        self.assertEqual(clair3_version("2.0.2"), "2.0.2")

    def test_parses_a_version_from_a_banner(self) -> None:
        self.assertEqual(clair3_version("Clair3 v2.0.2\nsomething else"), "2.0.2")

    def test_unparsable_output_returns_the_first_line_rather_than_a_guess(self) -> None:
        self.assertEqual(clair3_version("command not found"), "command not found")

    def test_empty_output_is_unknown(self) -> None:
        self.assertEqual(clair3_version("   "), "unknown")


class PreconditionTests(unittest.TestCase):
    PINNED = Clair3Policy(
        profile_id="t", expected_version="2.0.2", required_model_id="r1041_e82_400bps_sup_v500"
    )

    def test_an_unpinned_model_stops_the_run_before_the_binary_is_probed(self) -> None:
        """The sharper of the two gates: a wrong model does not error, it degrades quietly."""
        runner = _StubRunner(stdout="2.0.2")
        policy = Clair3Policy(profile_id="t", expected_version="2.0.2")
        failure = check_preconditions(runner, binary="clair3", policy=policy)
        assert isinstance(failure, PreconditionFailure)
        self.assertEqual(failure.reason, "model_not_pinned")
        self.assertEqual(runner.calls, [], "the binary must not be probed once already blocked")

    def test_a_matching_version_and_pinned_model_clears_the_run(self) -> None:
        runner = _StubRunner(stdout="Clair3 v2.0.2")
        self.assertIsNone(check_preconditions(runner, binary="clair3", policy=self.PINNED))
        self.assertEqual(runner.calls, [("clair3", "--version")])

    def test_a_version_mismatch_blocks_the_run(self) -> None:
        runner = _StubRunner(stdout="Clair3 v1.1.2")
        failure = check_preconditions(runner, binary="clair3", policy=self.PINNED)
        assert isinstance(failure, PreconditionFailure)
        self.assertEqual(failure.reason, "version_mismatch")
        self.assertIn("1.1.2", failure.detail)
        self.assertIn("2.0.2", failure.detail)

    def test_a_failed_probe_blocks_the_run(self) -> None:
        runner = _StubRunner(returncode=127, stderr="clair3: not found")
        failure = check_preconditions(runner, binary="clair3", policy=self.PINNED)
        assert isinstance(failure, PreconditionFailure)
        self.assertEqual(failure.reason, "version_probe_failed")

    def test_the_version_is_read_from_stderr_too(self) -> None:
        runner = _StubRunner(stdout="", stderr="Clair3 v2.0.2")
        self.assertIsNone(check_preconditions(runner, binary="clair3", policy=self.PINNED))

    def test_a_failure_is_returned_not_raised(self) -> None:
        """Preflight wants every blocked stage in one pass, not a stop at the first."""
        runner = _StubRunner(stdout="Clair3 v1.1.2")
        failure = check_preconditions(runner, binary="clair3", policy=self.PINNED)
        self.assertIsInstance(failure, PreconditionFailure)


if __name__ == "__main__":
    unittest.main()
