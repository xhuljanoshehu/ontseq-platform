from __future__ import annotations

import contextlib
import csv
import gzip
import hashlib
import io
import random
import tempfile
import unittest
import unittest.mock
from pathlib import Path

import yaml

from ontseq_platform import cli
from ontseq_platform.methylation_mixture import (
    NANOPOLISH_ADAPTER_VERSION,
    MethylationMixturePolicy,
    MethylationMixtureReport,
    MethylationMixtureReportStatus,
    MethylationReferenceMarker,
    NanopolishCallState,
    NanopolishSourceMetadata,
    TechnicalRecoveryAssessment,
    _Counts,
    _draw_estimate,
    _MarkerKey,
    _point_estimate,
    _quantile,
    _reference_markers,
    _selection_digest,
    parse_nanopolish_source,
    render_methylation_mixture_csv,
    render_methylation_mixture_html,
    run_nanopolish_mixture_analysis,
)
from ontseq_platform.models import GenomeBuild, ModuleRunStatus

_HEADER = (
    "chromosome",
    "strand",
    "start",
    "end",
    "read_name",
    "log_lik_ratio",
    "log_lik_methylated",
    "log_lik_unmethylated",
    "num_calling_strands",
    "num_motifs",
    "sequence",
)


def _row(
    *,
    chromosome: str = "chr1",
    start: int | str = 100,
    end: int | str | None = None,
    read_name: str = "SYN_READ_001",
    log_likelihood_ratio: float | str = 5.0,
    strand: str = "+",
    motif_count: int | str = 1,
    sequence: str = "ACGTA",
) -> dict[str, str]:
    start_text = str(start)
    if end is None:
        # Nanopolish call-methylation reports the inclusive motif-group endpoint. A
        # single-CpG group therefore commonly has identical start and end values.
        end = int(start_text)
    ratio_text = str(log_likelihood_ratio)
    return {
        "chromosome": chromosome,
        "strand": strand,
        "start": start_text,
        "end": str(end),
        "read_name": read_name,
        "log_lik_ratio": ratio_text,
        "log_lik_methylated": ratio_text,
        "log_lik_unmethylated": "0.0",
        "num_calling_strands": "1",
        "num_motifs": str(motif_count),
        "sequence": sequence,
    }


def _write_table(
    path: Path,
    rows: list[dict[str, str]],
    *,
    header: tuple[str, ...] = _HEADER,
) -> Path:
    lines = ["\t".join(header)]
    lines.extend("\t".join(row.get(column, "") for column in header) for row in rows)
    contents = "\n".join(lines) + "\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            handle.write(contents)
    else:
        path.write_text(contents, encoding="utf-8", newline="\n")
    return path


def _policy(**overrides: object) -> MethylationMixturePolicy:
    values: dict[str, object] = {
        "profile_id": "synthetic-paired-source",
        "status": "technical_defaults_only",
        "source_a_fractions": [0.0, 0.25, 0.5, 0.75, 1.0],
        "replicates": 2,
        "seed": 41,
        "calibration_read_group_fraction": 0.25,
        "total_mixture_read_groups": 20,
        "minimum_reference_valid_calls": 5,
        "minimum_mixture_valid_calls": 3,
        "minimum_absolute_delta_beta": 0.5,
        "minimum_markers": 6,
        "maximum_markers": 8,
        "uncertainty_draws": 100,
        "confidence_level": 0.95,
        "minimum_successful_uncertainty_fraction": 0.8,
        "maximum_confidence_interval_width": 1.0,
        "maximum_standardized_model_fit_rmse": 3.0,
        "minimum_quantifiable_fraction_for_recovery": 0.8,
        "maximum_absolute_mean_bias": 0.1,
        "maximum_mean_absolute_error": 0.15,
        "maximum_root_mean_squared_error": 0.2,
        "minimum_confidence_interval_coverage_fraction": 0.8,
        "note": "Synthetic technical test policy; not a validated threshold.",
    }
    values.update(overrides)
    return MethylationMixturePolicy.model_validate(values)


def _write_distinct_sources(
    root: Path,
    *,
    markers: int = 8,
    reads: int = 40,
    source_b_ratio: float = -5.0,
    source_b_sequence: str = "ACGTA",
    read_name_token: str = "SYN_READ",
) -> tuple[Path, Path]:
    rows_a: list[dict[str, str]] = []
    rows_b: list[dict[str, str]] = []
    for read_index in range(reads):
        for marker_index in range(markers):
            start = 1_000 + marker_index * 10
            rows_a.append(
                _row(
                    start=start,
                    read_name=f"{read_name_token}_A_{read_index:03d}",
                    log_likelihood_ratio=5.0,
                    sequence="ACGTA",
                )
            )
            rows_b.append(
                _row(
                    start=start,
                    read_name=f"{read_name_token}_B_{read_index:03d}",
                    log_likelihood_ratio=source_b_ratio,
                    sequence=source_b_sequence,
                )
            )
    return (
        _write_table(root / "source-a.tsv", rows_a),
        _write_table(root / "source-b.tsv", rows_b),
    )


def _write_extrapolating_sources(root: Path, policy: MethylationMixturePolicy) -> tuple[Path, Path]:
    """Create held-out extremes outside deliberately narrower calibration profiles."""
    names_a = [f"EXTRAP_A_{index:02d}" for index in range(8)]
    names_b = [f"EXTRAP_B_{index:02d}" for index in range(8)]
    shuffled_a = sorted(names_a)
    shuffled_b = sorted(names_b)
    random.Random(policy.seed + 1_000_003).shuffle(shuffled_a)
    random.Random(policy.seed + 2_000_003).shuffle(shuffled_b)
    calibration_a = set(shuffled_a[:4])
    calibration_b = set(shuffled_b[:4])
    methylated_calibration_a = set(sorted(calibration_a)[:3])
    methylated_calibration_b = {sorted(calibration_b)[0]}

    rows_a: list[dict[str, str]] = []
    rows_b: list[dict[str, str]] = []
    for marker_index in range(6):
        start = 2_000 + marker_index * 10
        for read_name in names_a:
            methylated = read_name not in calibration_a or read_name in methylated_calibration_a
            rows_a.append(
                _row(
                    start=start,
                    read_name=read_name,
                    log_likelihood_ratio=5.0 if methylated else -5.0,
                )
            )
        for read_name in names_b:
            methylated = read_name in methylated_calibration_b
            rows_b.append(
                _row(
                    start=start,
                    read_name=read_name,
                    log_likelihood_ratio=5.0 if methylated else -5.0,
                    sequence="ACGTA",
                )
            )
    return (
        _write_table(root / "extrapolation-a.tsv", rows_a),
        _write_table(root / "extrapolation-b.tsv", rows_b),
    )


def _write_differential_callability_sources(
    root: Path, policy: MethylationMixturePolicy
) -> tuple[Path, Path]:
    """Source A calls every marker; source B calls the same markers in 20% of reads."""
    names_a = [f"CALLABILITY_A_{index:02d}" for index in range(40)]
    names_b = [f"CALLABILITY_B_{index:02d}" for index in range(40)]
    shuffled_b = sorted(names_b)
    random.Random(policy.seed + 2_000_003).shuffle(shuffled_b)
    calibration_count = round(len(names_b) * policy.calibration_read_group_fraction)
    calibration_b = shuffled_b[:calibration_count]
    pool_b = shuffled_b[calibration_count:]
    callable_b = set(calibration_b[:2] + pool_b[:6])

    rows_a: list[dict[str, str]] = []
    rows_b: list[dict[str, str]] = []
    for read_name in names_a:
        # The ambiguous row registers all source read groups without becoming a marker.
        rows_a.append(_row(start=9_000, read_name=read_name, log_likelihood_ratio=0.0))
        for marker_index in range(6):
            rows_a.append(
                _row(
                    start=3_000 + marker_index * 10,
                    read_name=read_name,
                    log_likelihood_ratio=5.0,
                )
            )
    for read_name in names_b:
        rows_b.append(
            _row(
                start=9_000,
                read_name=read_name,
                log_likelihood_ratio=0.0,
                sequence="ACGTA",
            )
        )
        if read_name in callable_b:
            for marker_index in range(6):
                rows_b.append(
                    _row(
                        start=3_000 + marker_index * 10,
                        read_name=read_name,
                        log_likelihood_ratio=-5.0,
                        sequence="ACGTA",
                    )
                )
    return (
        _write_table(root / "callability-a.tsv", rows_a),
        _write_table(root / "callability-b.tsv", rows_b),
    )


def _write_sparse_marker_disagreement_sources(
    root: Path, policy: MethylationMixturePolicy
) -> tuple[Path, Path]:
    """Make low-rate held-out marker channels alternate between f=1 and f=0."""
    names_a = [f"DISAGREE_A_{index:03d}" for index in range(400)]
    names_b = [f"DISAGREE_B_{index:03d}" for index in range(400)]
    shuffled_a = sorted(names_a)
    shuffled_b = sorted(names_b)
    random.Random(policy.seed + 1_000_003).shuffle(shuffled_a)
    random.Random(policy.seed + 2_000_003).shuffle(shuffled_b)
    calibration_count = round(len(names_a) * policy.calibration_read_group_fraction)
    calibration_a = shuffled_a[:calibration_count]
    calibration_b = shuffled_b[:calibration_count]
    pool_a = shuffled_a[calibration_count:]
    pool_b = shuffled_b[calibration_count:]

    half_fraction_index = policy.source_a_fractions.index(0.5)
    half_seed = policy.seed + 10_000_019 + half_fraction_index * 10_007 + 101
    selected_a = random.Random(half_seed).sample(pool_a, 100)
    selected_b = random.Random(half_seed + 1).sample(pool_b, 100)

    calibration_signal_a = set(calibration_a[:5])
    calibration_signal_b = set(calibration_b[:5])
    held_out_signal_a = set(selected_a[:10])
    held_out_signal_b = set(selected_b[:10])
    rows_a: list[dict[str, str]] = []
    rows_b: list[dict[str, str]] = []
    for read_name in names_a:
        rows_a.append(_row(start=9_500, read_name=read_name, log_likelihood_ratio=0.0))
        for marker_index in range(6):
            if read_name in calibration_signal_a or (
                marker_index < 3 and read_name in held_out_signal_a
            ):
                rows_a.append(
                    _row(
                        start=4_000 + marker_index * 10,
                        read_name=read_name,
                        log_likelihood_ratio=5.0,
                    )
                )
    for read_name in names_b:
        rows_b.append(
            _row(
                start=9_500,
                read_name=read_name,
                log_likelihood_ratio=0.0,
                sequence="ACGTA",
            )
        )
        for marker_index in range(6):
            if read_name in calibration_signal_b or (
                marker_index >= 3 and read_name in held_out_signal_b
            ):
                rows_b.append(
                    _row(
                        start=4_000 + marker_index * 10,
                        read_name=read_name,
                        log_likelihood_ratio=-5.0,
                        sequence="ACGTA",
                    )
                )
    return (
        _write_table(root / "disagreement-a.tsv", rows_a),
        _write_table(root / "disagreement-b.tsv", rows_b),
    )


def _write_biased_recovery_sources(
    root: Path, policy: MethylationMixturePolicy
) -> tuple[Path, Path]:
    """Keep B low in calibration but make both held-out source pools A-like."""

    names_a = [f"BIASED_A_{index:03d}" for index in range(40)]
    names_b = [f"BIASED_B_{index:03d}" for index in range(40)]
    shuffled_b = sorted(names_b)
    random.Random(policy.seed + 2_000_003).shuffle(shuffled_b)
    calibration_count = round(len(names_b) * policy.calibration_read_group_fraction)
    calibration_b = set(shuffled_b[:calibration_count])
    rows_a: list[dict[str, str]] = []
    rows_b: list[dict[str, str]] = []
    for marker_index in range(8):
        start = 5_000 + marker_index * 10
        for read_name in names_a:
            rows_a.append(_row(start=start, read_name=read_name, log_likelihood_ratio=5.0))
        for read_name in names_b:
            rows_b.append(
                _row(
                    start=start,
                    read_name=read_name,
                    log_likelihood_ratio=-5.0 if read_name in calibration_b else 5.0,
                )
            )
    return (
        _write_table(root / "biased-recovery-a.tsv", rows_a),
        _write_table(root / "biased-recovery-b.tsv", rows_b),
    )


def _run(
    source_a: Path,
    source_b: Path,
    *,
    policy: MethylationMixturePolicy | None = None,
    source_a_id: str = "SOURCE_A",
    source_b_id: str = "SOURCE_B",
    source_a_build: GenomeBuild = GenomeBuild.GRCH38,
    source_b_build: GenomeBuild = GenomeBuild.GRCH38,
    analysis_id: str = "SYNTHETIC_MIXTURE",
    source_a_metadata: NanopolishSourceMetadata | None = None,
    source_b_metadata: NanopolishSourceMetadata | None = None,
    sources_declared_biologically_distinct: bool = True,
    software_version: str = "0.4.1-test",
    git_commit: str = "UNKNOWN",
) -> MethylationMixtureReport:
    return run_nanopolish_mixture_analysis(
        source_a,
        source_b,
        source_a_id=source_a_id,
        source_b_id=source_b_id,
        source_a_genome_build=source_a_build,
        source_b_genome_build=source_b_build,
        analysis_id=analysis_id,
        policy=policy or _policy(),
        sources_declared_biologically_distinct=sources_declared_biologically_distinct,
        source_a_metadata=source_a_metadata,
        source_b_metadata=source_b_metadata,
        software_version=software_version,
        git_commit=git_commit,
    )


def _known_metadata(**overrides: object) -> NanopolishSourceMetadata:
    values: dict[str, object] = {
        "caller_name": "nanopolish",
        "caller_version": "0.13.3",
        "reference_id": "GRCh38-primary",
        "reference_genome_build": "GRCh38",
        "reference_sha256": "b" * 64,
        "sequencing_platform": "Oxford Nanopore",
        "flow_cell_product_code": "FLO-MIN106",
        "library_kit": "SQK-LSK109",
        "basecaller_name": "Albacore",
        "basecaller_version": "2.3.4",
        "basecaller_model": "r9.4.1_450bps",
        "modification_context": "CpG_5mC_vs_unmethylated",
    }
    values.update(overrides)
    return NanopolishSourceMetadata.model_validate(values)


class NanopolishSourceParsingTests(unittest.TestCase):
    def test_plain_and_gzip_inputs_classify_calls_and_count_exclusions(self) -> None:
        rows = [
            _row(read_name="METHYLATED", log_likelihood_ratio=3.0),
            _row(
                chromosome="2",
                start=200,
                read_name="UNMETHYLATED",
                log_likelihood_ratio=-3.0,
            ),
            _row(
                chromosome="chrX",
                start=300,
                read_name="AMBIGUOUS",
                log_likelihood_ratio=0.0,
            ),
            _row(
                chromosome="chrUn_KI270302v1",
                start=400,
                read_name="EXCLUDED",
                log_likelihood_ratio=4.0,
            ),
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for filename in ("calls.tsv", "calls.tsv.gz"):
                with self.subTest(filename=filename):
                    parsed = parse_nanopolish_source(
                        _write_table(root / filename, rows),
                        source_id="SYN_SOURCE",
                        genome_build=GenomeBuild.GRCH38,
                        policy=_policy(),
                    )
                    self.assertEqual(parsed.summary.rows_total, 4)
                    self.assertEqual(parsed.summary.canonical_call_rows, 3)
                    self.assertEqual(parsed.summary.informative_call_rows, 2)
                    self.assertEqual(parsed.summary.ambiguous_call_rows, 1)
                    self.assertEqual(parsed.summary.skipped_noncanonical_rows, 1)
                    self.assertEqual(parsed.summary.read_group_count, 3)
                    self.assertEqual(
                        parsed.calls_by_read["METHYLATED"][0].state,
                        NanopolishCallState.METHYLATED,
                    )
                    self.assertEqual(
                        parsed.calls_by_read["UNMETHYLATED"][0].state,
                        NanopolishCallState.UNMETHYLATED,
                    )
                    self.assertEqual(
                        parsed.calls_by_read["AMBIGUOUS"][0].state,
                        NanopolishCallState.AMBIGUOUS,
                    )
                    self.assertEqual(
                        parsed.calls_by_read["UNMETHYLATED"][0].marker.chromosome,
                        "chr2",
                    )
                    self.assertIsNotNone(parsed.summary.fingerprint.sha256)

    def test_missing_required_column_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_table(
                Path(directory) / "missing.tsv",
                [_row()],
                header=tuple(column for column in _HEADER if column != "sequence"),
            )
            with self.assertRaisesRegex(ValueError, "missing columns: sequence"):
                parse_nanopolish_source(
                    path,
                    source_id="SYN_SOURCE",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                )

    def test_likelihood_threshold_is_multiplied_by_number_of_motifs(self) -> None:
        rows = [
            _row(
                start=100,
                read_name="BELOW_POSITIVE",
                log_likelihood_ratio=4.99,
                motif_count=2,
            ),
            _row(
                start=200,
                read_name="AT_POSITIVE",
                log_likelihood_ratio=5.0,
                motif_count=2,
            ),
            _row(
                start=300,
                read_name="BELOW_NEGATIVE",
                log_likelihood_ratio=-4.99,
                motif_count=2,
            ),
            _row(
                start=400,
                read_name="AT_NEGATIVE",
                log_likelihood_ratio=-5.0,
                motif_count=2,
            ),
        ]
        with tempfile.TemporaryDirectory() as directory:
            parsed = parse_nanopolish_source(
                _write_table(Path(directory) / "motif-threshold.tsv", rows),
                source_id="SYN_SOURCE",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(log_likelihood_ratio_threshold=2.5),
            )

        state_by_read = {
            read_name: calls[0].state for read_name, calls in parsed.calls_by_read.items()
        }
        self.assertEqual(state_by_read["BELOW_POSITIVE"], NanopolishCallState.AMBIGUOUS)
        self.assertEqual(state_by_read["AT_POSITIVE"], NanopolishCallState.METHYLATED)
        self.assertEqual(state_by_read["BELOW_NEGATIVE"], NanopolishCallState.AMBIGUOUS)
        self.assertEqual(state_by_read["AT_NEGATIVE"], NanopolishCallState.UNMETHYLATED)
        self.assertEqual(parsed.summary.informative_call_rows, 2)
        self.assertEqual(parsed.summary.ambiguous_call_rows, 2)

    def test_malformed_rows_are_refused_instead_of_silently_coerced(self) -> None:
        cases = {
            "invalid coordinate": (_row(start="not-an-integer", end=200), "invalid start"),
            "reversed interval": (_row(start=200, end=100), "invalid coordinates"),
            "invalid strand": (_row(strand="*"), "invalid strand"),
            "missing read": (_row(read_name=""), "empty read_name"),
            "missing sequence": (_row(sequence=""), "empty sequence"),
            "non-finite likelihood": (
                _row(log_likelihood_ratio="nan"),
                "non-finite log_lik_",
            ),
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index, (name, (row, message)) in enumerate(cases.items()):
                with self.subTest(case=name):
                    path = _write_table(root / f"malformed-{index}.tsv", [row])
                    with self.assertRaisesRegex(ValueError, message):
                        parse_nanopolish_source(
                            path,
                            source_id="SYN_SOURCE",
                            genome_build=GenomeBuild.GRCH38,
                            policy=_policy(),
                        )

    def test_duplicate_marker_within_one_read_group_is_refused(self) -> None:
        duplicate = _row(read_name="SYN_DUPLICATE", start=100)
        with tempfile.TemporaryDirectory() as directory:
            path = _write_table(Path(directory) / "duplicate.tsv", [duplicate, duplicate])
            with self.assertRaisesRegex(ValueError, "duplicates a marker"):
                parse_nanopolish_source(
                    path,
                    source_id="SYN_SOURCE",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                )

    def test_same_coordinates_on_opposite_strands_are_distinct_exact_markers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_table(
                Path(directory) / "opposite-strands.tsv",
                [
                    _row(read_name="SYN_BISTRAND", start=100, strand="+"),
                    _row(read_name="SYN_BISTRAND", start=100, strand="-"),
                ],
            )
            parsed = parse_nanopolish_source(
                path,
                source_id="SYN_SOURCE",
                genome_build=GenomeBuild.GRCH38,
                policy=_policy(),
            )

        calls = parsed.calls_by_read["SYN_BISTRAND"]
        self.assertEqual(len(calls), 2)
        self.assertEqual({call.marker.strand for call in calls}, {"+", "-"})
        self.assertEqual(len({call.marker.marker_id for call in calls}), 2)

    def test_file_without_canonical_calls_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_table(
                Path(directory) / "noncanonical.tsv",
                [_row(chromosome="chrUn_KI270302v1")],
            )
            with self.assertRaisesRegex(ValueError, "no canonical methylation calls"):
                parse_nanopolish_source(
                    path,
                    source_id="SYN_SOURCE",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(),
                )

    def test_input_byte_limit_fails_before_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_table(Path(directory) / "too-large.tsv", [_row()])
            with self.assertRaisesRegex(ValueError, "in-memory safety limit"):
                parse_nanopolish_source(
                    path,
                    source_id="SYN_SOURCE",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(maximum_input_bytes_per_source=1),
                )

    def test_input_row_limit_fails_at_the_first_excess_row(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_table(
                Path(directory) / "too-many-rows.tsv",
                [
                    _row(start=100, read_name="SYN_ROW_1"),
                    _row(start=200, read_name="SYN_ROW_2"),
                ],
            )
            with self.assertRaisesRegex(ValueError, "row safety limit of 1"):
                parse_nanopolish_source(
                    path,
                    source_id="SYN_SOURCE",
                    genome_build=GenomeBuild.GRCH38,
                    policy=_policy(maximum_rows_per_source=1),
                )


class NanopolishBoundedIntakeTests(unittest.TestCase):
    @staticmethod
    def _parse(path: Path, *, maximum_bytes: int = 100_000_000):
        return parse_nanopolish_source(
            path,
            source_id="SYN_SECURITY_CHECK",
            genome_build=GenomeBuild.GRCH38,
            policy=_policy(maximum_input_bytes_per_source=maximum_bytes),
        )

    def test_compressed_under_limit_but_decoded_over_limit_is_rejected(self) -> None:
        rows = [_row(read_name=f"SYN_READ_{index:03d}") for index in range(200)]
        decoded = (
            "\t".join(_HEADER)
            + "\n"
            + "".join("\t".join(row[column] for column in _HEADER) + "\n" for row in rows)
        ).encode()
        compressed = gzip.compress(decoded, mtime=0)
        self.assertLess(len(compressed), 1024)
        self.assertGreater(len(decoded), 1024)
        self.assertLess(len(decoded), 65_536)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bounded.tsv.gz"
            path.write_bytes(compressed)
            with self.assertRaisesRegex(ValueError, "decoded-byte safety limit"):
                self._parse(path, maximum_bytes=1024)

    def test_decoded_budget_counts_blank_lines_and_concatenated_gzip_members(self) -> None:
        first = ("\t".join(_HEADER) + "\n").encode()
        compressed = gzip.compress(first, mtime=0) + gzip.compress(b"\n" * 2048, mtime=0)
        self.assertLess(len(compressed), 1024)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "members.tsv.gz"
            path.write_bytes(compressed)
            with self.assertRaisesRegex(ValueError, "decoded-byte safety limit"):
                self._parse(path, maximum_bytes=1024)

    def test_exact_decoded_budget_preserves_utf8_bom_crlf_and_security_provenance(self) -> None:
        row = _row(read_name="SYN_READ_é")
        decoded = (
            "\ufeff"
            + "\t".join(_HEADER)
            + "\r\n"
            + "\t".join(row[column] for column in _HEADER)
            + "\r\n"
        ).encode()
        with tempfile.TemporaryDirectory() as directory:
            for suffix in (".tsv", ".tsv.gz"):
                with self.subTest(suffix=suffix):
                    path = Path(directory) / ("exact" + suffix)
                    path.write_bytes(
                        gzip.compress(decoded, mtime=0) if suffix == ".tsv.gz" else decoded
                    )
                    parsed = self._parse(path, maximum_bytes=len(decoded))
                    self.assertEqual(parsed.summary.decoded_size_bytes, len(decoded))
                    self.assertEqual(parsed.summary.adapter_version, NANOPOLISH_ADAPTER_VERSION)
                    self.assertEqual(parsed.summary.input_security_profile, "bounded-tsv-gzip-v1")
                    self.assertEqual(tuple(parsed.calls_by_read), ("SYN_READ_é",))

    def test_header_size_width_and_multiline_header_are_bounded_before_csv(self) -> None:
        headers = {
            "oversized": (b"x" * 16_385, "header exceeds"),
            "wide": (
                ("\t".join((*_HEADER, *(f"extra{i}" for i in range(64)))) + "\n").encode(),
                "column safety limit",
            ),
            "multiline": (b'"chromosome\n', "invalid single-line header"),
        }
        with tempfile.TemporaryDirectory() as directory:
            for name, (contents, error) in headers.items():
                with self.subTest(name=name):
                    path = Path(directory) / (name + ".tsv.gz")
                    path.write_bytes(gzip.compress(contents, mtime=0))
                    with self.assertRaisesRegex(ValueError, error):
                        self._parse(path)

    def test_large_physical_data_line_is_rejected_with_bounded_read(self) -> None:
        # A reduced test limit exercises the production branch without a large fixture.
        decoded = ("\t".join(_HEADER) + "\n").encode() + b"x" * 300
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "line.tsv.gz"
            path.write_bytes(gzip.compress(decoded, mtime=0))
            with (
                unittest.mock.patch(
                    "ontseq_platform.methylation_mixture._NANOPOLISH_MAXIMUM_PHYSICAL_LINE_BYTES",
                    256,
                ),
                self.assertRaisesRegex(ValueError, "physical line exceeds the 256-byte"),
            ):
                self._parse(path)

    def test_quoted_read_identifiers_cannot_contain_ascii_or_c1_controls(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for control in ("\x00", "\t", "\n", "\r", "\x1f", "\x7f", "\x85", "\x9f"):
                with self.subTest(control=repr(control)):
                    row = _row(read_name="SYN" + control + "READ")
                    contents = io.StringIO(newline="")
                    writer = csv.DictWriter(contents, fieldnames=_HEADER, delimiter="\t")
                    writer.writeheader()
                    writer.writerow(row)
                    path = Path(directory) / "controls.tsv.gz"
                    path.write_bytes(gzip.compress(contents.getvalue().encode(), mtime=0))
                    with self.assertRaisesRegex(ValueError, "forbidden control character"):
                        self._parse(path)

    def test_ambiguous_digest_pools_are_rejected_and_valid_digest_bytes_are_unchanged(self) -> None:
        for names in (("a\nA\x00b", "c"), ("a", "b\nA\x00c")):
            with self.subTest(names=repr(names)):
                with self.assertRaisesRegex(ValueError, "forbidden control character"):
                    _selection_digest(names, (), context="SYNTHETIC_SAME_CONTEXT")
                with self.assertRaisesRegex(ValueError, "forbidden control character"):
                    _selection_digest((), names, context="SYNTHETIC_SAME_CONTEXT")
        legacy_bytes = b"context=synthetic\nA\0READ_A\nA\0READ_B\nB\0READ_C\n"
        self.assertEqual(
            _selection_digest(("READ_B", "READ_A"), ("READ_C",), context="synthetic"),
            hashlib.sha256(legacy_bytes).hexdigest(),
        )


class PairedSourceMixtureTests(unittest.TestCase):
    def test_known_proportions_are_recovered_with_a_constant_read_group_budget(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            report = _run(source_a, source_b)

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(
            report.source_fraction_definition,
            "fraction_of_held_out_read_groups_from_source_a",
        )
        self.assertEqual(len(report.reference_markers), 8)
        self.assertEqual(report.completed_levels, 10)
        self.assertEqual(report.no_call_levels, 0)
        self.assertEqual({level.total_read_groups for level in report.levels}, {20})
        self.assertEqual(
            [level.target_source_a_fraction for level in report.levels[::2]],
            [0.0, 0.25, 0.5, 0.75, 1.0],
        )
        for level in report.levels:
            self.assertEqual(level.status, ModuleRunStatus.COMPLETED)
            self.assertAlmostEqual(
                level.estimated_source_a_fraction or 0.0,
                level.realized_source_a_read_group_fraction,
            )
            self.assertEqual(level.markers_used, 8)
            self.assertEqual(level.valid_mixture_calls, 160)
            self.assertEqual(level.successful_uncertainty_draws, 100)
            self.assertLessEqual(
                level.confidence_interval_lower or 0.0,
                level.estimated_source_a_fraction or 0.0,
            )
            self.assertGreaterEqual(
                level.confidence_interval_upper or 0.0,
                level.estimated_source_a_fraction or 0.0,
            )
        self.assertAlmostEqual(report.mean_bias or 0.0, 0.0)
        self.assertAlmostEqual(report.mean_absolute_error or 0.0, 0.0)
        self.assertAlmostEqual(report.root_mean_squared_error or 0.0, 0.0)

    def test_same_seed_reproduces_entire_report_and_new_seed_changes_selections(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            first = _run(source_a, source_b)
            repeated = _run(source_a, source_b)
            changed = _run(source_a, source_b, policy=_policy(seed=42))

        self.assertEqual(first.model_dump(), repeated.model_dump())
        self.assertEqual(
            [level.selection_sha256 for level in first.levels],
            [level.selection_sha256 for level in repeated.levels],
        )
        self.assertNotEqual(
            first.calibration_selection_sha256,
            changed.calibration_selection_sha256,
        )
        self.assertNotEqual(
            [level.selection_sha256 for level in first.levels],
            [level.selection_sha256 for level in changed.levels],
        )

    def test_absent_reference_contrast_returns_no_call_without_numeric_estimates(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(
                Path(directory),
                source_b_ratio=5.0,
                source_b_sequence="ACGTA",
            )
            report = _run(source_a, source_b)

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.reference_markers, [])
        self.assertEqual(report.completed_levels, 0)
        self.assertEqual(report.no_call_levels, 10)
        self.assertIsNone(report.mean_bias)
        self.assertIsNone(report.mean_absolute_error)
        self.assertIsNone(report.root_mean_squared_error)
        for level in report.levels:
            self.assertEqual(level.status, ModuleRunStatus.NO_CALL)
            self.assertIsNone(level.estimated_source_a_fraction)
            self.assertIsNone(level.confidence_interval_lower)
            self.assertIsNone(level.confidence_interval_upper)
            self.assertIsNone(level.model_fit_standardized_rmse)
            self.assertIn("reference marker", level.no_call_reason or "")

    def test_insufficient_mixture_coverage_returns_no_call_and_zero_used_markers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            report = _run(
                source_a,
                source_b,
                policy=_policy(minimum_mixture_valid_calls=21),
            )

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(len(report.reference_markers), 8)
        for level in report.levels:
            self.assertEqual(level.status, ModuleRunStatus.NO_CALL)
            self.assertEqual(level.markers_used, 0)
            self.assertIsNone(level.estimated_source_a_fraction)
            self.assertIn("valid-call floor", level.no_call_reason or "")

    def test_strong_extrapolation_returns_no_call_without_a_bounded_estimate(self) -> None:
        policy = _policy(
            source_a_fractions=[0.0, 0.5, 1.0],
            replicates=1,
            calibration_read_group_fraction=0.5,
            total_mixture_read_groups=4,
            minimum_reference_valid_calls=4,
            minimum_mixture_valid_calls=4,
            minimum_absolute_delta_beta=0.2,
            minimum_markers=6,
            maximum_markers=6,
            maximum_standardized_model_fit_rmse=3.0,
            maximum_extrapolation=0.1,
        )
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_extrapolating_sources(Path(directory), policy)
            report = _run(source_a, source_b, policy=policy)

        by_target = {level.target_source_a_fraction: level for level in report.levels}
        self.assertEqual(by_target[0.0].status, ModuleRunStatus.NO_CALL)
        self.assertEqual(by_target[1.0].status, ModuleRunStatus.NO_CALL)
        self.assertEqual(by_target[0.5].status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.status, MethylationMixtureReportStatus.PARTIAL)
        self.assertEqual(report.completed_levels, 1)
        self.assertEqual(report.no_call_levels, 2)
        self.assertAlmostEqual(report.quantifiable_fraction, 1 / 3)
        self.assertFalse(report.grid_fully_quantified)
        self.assertEqual(report.accuracy_metrics_scope, "completed_levels_only")
        self.assertEqual(
            report.technical_recovery_assessment,
            TechnicalRecoveryAssessment.NOT_EVALUABLE,
        )
        self.assertEqual(len(report.technical_recovery_reasons), 1)
        self.assertIn("quantifiable share", report.technical_recovery_reasons[0])
        completed_error = (by_target[0.5].estimated_source_a_fraction or 0.0) - by_target[
            0.5
        ].realized_source_a_read_group_fraction
        self.assertAlmostEqual(report.mean_bias or 0.0, completed_error)
        self.assertAlmostEqual(report.mean_absolute_error or 0.0, abs(completed_error))
        self.assertAlmostEqual(report.root_mean_squared_error or 0.0, abs(completed_error))
        for target in (0.0, 1.0):
            level = by_target[target]
            self.assertIsNone(level.estimated_source_a_fraction)
            self.assertIsNone(level.unconstrained_source_a_fraction)
            self.assertIn("allowed extrapolation", level.no_call_reason or "")

    def test_differential_callability_does_not_bias_the_half_mixture_toward_source_a(
        self,
    ) -> None:
        policy = _policy(
            source_a_fractions=[0.0, 0.5, 1.0],
            replicates=1,
            minimum_reference_valid_calls=2,
            minimum_mixture_valid_calls=1,
            minimum_markers=6,
            maximum_markers=6,
            maximum_standardized_model_fit_rmse=3.0,
        )
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_differential_callability_sources(Path(directory), policy)
            report = _run(source_a, source_b, policy=policy)

        half = next(level for level in report.levels if level.target_source_a_fraction == 0.5)
        self.assertEqual(half.status, ModuleRunStatus.COMPLETED)
        self.assertAlmostEqual(half.realized_source_a_read_group_fraction, 0.5)
        self.assertIsNotNone(half.estimated_source_a_fraction)
        assert half.estimated_source_a_fraction is not None
        self.assertAlmostEqual(half.estimated_source_a_fraction, 0.5, delta=0.05)
        # A beta-only estimator would be 10 / (10 + sparse B calls), well above 0.5.
        self.assertLess(half.estimated_source_a_fraction, 0.6)

    def test_sparse_marker_disagreement_fails_the_standardized_fit_gate(self) -> None:
        policy = _policy(
            source_a_fractions=[0.0, 0.5, 1.0],
            replicates=1,
            calibration_read_group_fraction=0.25,
            total_mixture_read_groups=200,
            minimum_reference_valid_calls=5,
            minimum_mixture_valid_calls=5,
            minimum_markers=6,
            maximum_markers=6,
            maximum_standardized_model_fit_rmse=0.5,
            maximum_extrapolation=1.0,
        )
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_sparse_marker_disagreement_sources(Path(directory), policy)
            report = _run(source_a, source_b, policy=policy)

        half = next(level for level in report.levels if level.target_source_a_fraction == 0.5)
        self.assertEqual(half.realized_source_a_read_group_fraction, 0.5)
        self.assertEqual(half.status, ModuleRunStatus.NO_CALL)
        self.assertEqual(half.markers_used, 6)
        self.assertEqual(half.valid_mixture_calls, 60)
        self.assertIsNone(half.estimated_source_a_fraction)
        self.assertIsNone(half.model_fit_standardized_rmse)
        self.assertIn("standardized", half.no_call_reason or "")
        self.assertIn("exceeds 0.5000", half.no_call_reason or "")


class PersistedContractValidationTests(unittest.TestCase):
    def test_reference_marker_rejects_tampered_id_and_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            marker = _run(source_a, source_b).reference_markers[0]

        tampered_id = marker.model_dump(mode="json")
        tampered_id["marker_id"] = f"{tampered_id['marker_id']}.altered"
        with self.assertRaisesRegex(ValueError, "inconsistent with its exact coordinates"):
            MethylationReferenceMarker.model_validate(tampered_id)

        reversed_coordinate = marker.model_dump(mode="json")
        reversed_coordinate["start"] = reversed_coordinate["end"] + 1
        reversed_coordinate["marker_id"] = (
            f"{reversed_coordinate['chromosome']}:{reversed_coordinate['start']}-"
            f"{reversed_coordinate['end']}:{reversed_coordinate['strand']}:"
            f"motifs={reversed_coordinate['motif_count']}:"
            f"sequence_sha256={reversed_coordinate['sequence_sha256']}"
        )
        with self.assertRaisesRegex(ValueError, "end cannot be smaller than start"):
            MethylationReferenceMarker.model_validate(reversed_coordinate)

    def test_completed_level_must_still_satisfy_every_locked_policy_gate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            report = _run(source_a, source_b)

        too_few_markers = report.model_dump(mode="json")
        level = too_few_markers["levels"][0]
        for observation in level["marker_observations"][:3]:
            observation["modified_calls"] = 0
            observation["unmethylated_calls"] = 0
            observation["valid_calls"] = 0
            observation["beta"] = None
        level["markers_used"] = 5
        level["valid_mixture_calls"] = sum(
            observation["valid_calls"] for observation in level["marker_observations"]
        )
        with self.assertRaisesRegex(ValueError, "minimum-marker policy"):
            MethylationMixtureReport.model_validate(too_few_markers)

        wide_interval = report.model_dump(mode="json")
        wide_interval["policy"]["maximum_confidence_interval_width"] = 0.99
        wide_interval["levels"][0]["confidence_interval_lower"] = 0.0
        wide_interval["levels"][0]["confidence_interval_upper"] = 1.0
        with self.assertRaisesRegex(ValueError, "interval-width policy"):
            MethylationMixtureReport.model_validate(wide_interval)

        poor_fit = report.model_dump(mode="json")
        poor_fit["levels"][0]["model_fit_standardized_rmse"] = 3.01
        with self.assertRaisesRegex(ValueError, "fit is inconsistent with marker counts"):
            MethylationMixtureReport.model_validate(poor_fit)

        extrapolated = report.model_dump(mode="json")
        extrapolated["levels"][0]["unconstrained_source_a_fraction"] = 1.11
        extrapolated["levels"][0]["estimated_source_a_fraction"] = 1.0
        with self.assertRaisesRegex(ValueError, "raw estimate is inconsistent with marker counts"):
            MethylationMixtureReport.model_validate(extrapolated)

        too_few_draws = report.model_dump(mode="json")
        too_few_draws["levels"][0]["successful_uncertainty_draws"] = 79
        with self.assertRaisesRegex(ValueError, "uncertainty-draw policy"):
            MethylationMixtureReport.model_validate(too_few_draws)

        inconsistent_points = report.model_dump(mode="json")
        inconsistent_points["levels"][0]["unconstrained_source_a_fraction"] = 0.4
        inconsistent_points["levels"][0]["estimated_source_a_fraction"] = 0.5
        with self.assertRaisesRegex(ValueError, "raw estimate is inconsistent with marker counts"):
            MethylationMixtureReport.model_validate(inconsistent_points)

    def test_report_rejects_tampered_hash_counts_metrics_fraction_and_status(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            report = _run(source_a, source_b)

        marker_hash = report.model_dump(mode="json")
        marker_hash["reference_marker_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "marker lock checksum"):
            MethylationMixtureReport.model_validate(marker_hash)

        calibration_counts = report.model_dump(mode="json")
        calibration_counts["source_a_calibration_read_groups"] += 1
        calibration_counts["source_a_mixture_pool_read_groups"] -= 1
        with self.assertRaisesRegex(ValueError, "pool sizes do not match"):
            MethylationMixtureReport.model_validate(calibration_counts)

        aggregate_counts = report.model_dump(mode="json")
        aggregate_counts["completed_levels"] = 9
        with self.assertRaisesRegex(ValueError, "level counts are inconsistent"):
            MethylationMixtureReport.model_validate(aggregate_counts)

        metrics = report.model_dump(mode="json")
        metrics["mean_bias"] = 0.25
        with self.assertRaisesRegex(ValueError, "metrics are inconsistent"):
            MethylationMixtureReport.model_validate(metrics)

        observations = report.model_dump(mode="json")
        for observation in observations["levels"][0]["marker_observations"]:
            observation["modified_calls"], observation["unmethylated_calls"] = (
                observation["unmethylated_calls"],
                observation["modified_calls"],
            )
            observation["beta"] = (
                observation["modified_calls"] / observation["valid_calls"]
                if observation["valid_calls"]
                else None
            )
        with self.assertRaisesRegex(ValueError, "estimate is inconsistent with marker counts"):
            MethylationMixtureReport.model_validate(observations)

        quantifiable = report.model_dump(mode="json")
        quantifiable["quantifiable_fraction"] = 0.9
        with self.assertRaisesRegex(ValueError, "quantifiable_fraction is inconsistent"):
            MethylationMixtureReport.model_validate(quantifiable)

        status = report.model_dump(mode="json")
        status["status"] = "PARTIAL"
        with self.assertRaisesRegex(ValueError, "status is inconsistent"):
            MethylationMixtureReport.model_validate(status)

    def test_report_binds_calibration_split_and_level_seeds_to_policy(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            report = _run(source_a, source_b)

        changed_split = report.model_dump(mode="json")
        changed_split["policy"]["calibration_read_group_fraction"] = 0.5
        with self.assertRaisesRegex(ValueError, "pool sizes do not match"):
            MethylationMixtureReport.model_validate(changed_split)

        changed_seed = report.model_dump(mode="json")
        changed_seed["policy"]["seed"] += 1
        with self.assertRaisesRegex(ValueError, "level seed does not match"):
            MethylationMixtureReport.model_validate(changed_seed)

    def test_no_call_report_still_binds_split_and_level_seed_to_policy(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory), source_b_ratio=5.0)
            report = _run(source_a, source_b)

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        changed_seed = report.model_dump(mode="json")
        changed_seed["levels"][0]["seed"] += 1
        with self.assertRaisesRegex(ValueError, "level seed does not match"):
            MethylationMixtureReport.model_validate(changed_seed)

        changed_split = report.model_dump(mode="json")
        changed_split["source_a_calibration_read_groups"] += 1
        changed_split["source_a_mixture_pool_read_groups"] -= 1
        with self.assertRaisesRegex(ValueError, "pool sizes do not match"):
            MethylationMixtureReport.model_validate(changed_split)


class TechnicalRecoveryAssessmentTests(unittest.TestCase):
    def test_exact_delta_beta_threshold_is_stable_under_binary_rounding(self) -> None:
        sequence_sha256 = "a" * 64
        keys = [
            _MarkerKey("chr1", "+", 8_000 + index, 8_000 + index, 1, sequence_sha256)
            for index in range(20)
        ]
        selected = _reference_markers(
            {key: _Counts(modified=6, unmethylated=4) for key in keys},
            {key: _Counts(modified=4, unmethylated=6) for key in keys},
            source_a_total_read_groups=10,
            source_b_total_read_groups=10,
            policy=_policy(
                minimum_absolute_delta_beta=0.2,
                minimum_markers=20,
                maximum_markers=20,
            ),
        )

        self.assertEqual(len(selected), 20)

    def test_conditional_interval_contains_perfect_low_coverage_endpoints(self) -> None:
        markers = [
            MethylationReferenceMarker(
                marker_id=(
                    f"chr1:{10_000 + index}-{10_000 + index}:+:motifs=1:sequence_sha256={'a' * 64}"
                ),
                chromosome="chr1",
                strand="+",
                start=10_000 + index,
                end=10_000 + index,
                motif_count=1,
                sequence_sha256="a" * 64,
                source_a_modified_calls=7,
                source_a_unmethylated_calls=3,
                source_a_ambiguous_calls=0,
                source_a_valid_calls=10,
                source_a_total_read_groups=10,
                source_b_modified_calls=3,
                source_b_unmethylated_calls=7,
                source_b_ambiguous_calls=0,
                source_b_valid_calls=10,
                source_b_total_read_groups=10,
                source_a_beta=0.7,
                source_b_beta=0.3,
                delta_beta=0.4,
            )
            for index in range(200)
        ]
        cases = ((0.0, 6, 14), (1.0, 14, 6))
        for expected, modified, unmethylated in cases:
            observations = {
                marker.marker_id: _Counts(modified=modified, unmethylated=unmethylated)
                for marker in markers
            }
            raw, point, fit, _serialized = _point_estimate(
                markers,
                observations,
                minimum_mixture_valid_calls=3,
                mixture_total_read_groups=20,
            )
            self.assertAlmostEqual(raw or 0.0, expected)
            self.assertAlmostEqual(point or 0.0, expected)
            self.assertAlmostEqual(fit or 0.0, 0.0)
            rng = random.Random(20260904)
            draws = sorted(
                value
                for value in (
                    _draw_estimate(
                        markers,
                        observations,
                        minimum_mixture_valid_calls=3,
                        mixture_total_read_groups=20,
                        fitted_source_a_fraction=expected,
                        rng=rng,
                    )
                    for _ in range(500)
                )
                if value is not None
            )
            self.assertLessEqual(_quantile(draws, 0.025), expected)
            self.assertGreaterEqual(_quantile(draws, 0.975), expected)

    def test_exact_synthetic_recovery_passes_all_locked_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            report = _run(source_a, source_b)

        self.assertEqual(report.technical_recovery_assessment, TechnicalRecoveryAssessment.PASS)
        self.assertEqual(report.technical_recovery_reasons, [])
        self.assertEqual(report.quantifiable_fraction, 1.0)
        self.assertEqual(report.mean_bias, 0.0)
        self.assertEqual(report.mean_absolute_error, 0.0)
        self.assertEqual(report.root_mean_squared_error, 0.0)
        self.assertEqual(report.confidence_interval_coverage_fraction, 1.0)

    def test_observed_bad_recovery_fails_all_four_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            policy = _policy()
            source_a, source_b = _write_biased_recovery_sources(Path(directory), policy)
            failed = _run(source_a, source_b, policy=policy)

        self.assertEqual(failed.technical_recovery_assessment, TechnicalRecoveryAssessment.FAIL)
        self.assertGreater(abs(failed.mean_bias or 0.0), failed.policy.maximum_absolute_mean_bias)
        self.assertGreater(
            failed.mean_absolute_error or 0.0,
            failed.policy.maximum_mean_absolute_error,
        )
        self.assertGreater(
            failed.root_mean_squared_error or 0.0,
            failed.policy.maximum_root_mean_squared_error,
        )
        self.assertLess(
            failed.confidence_interval_coverage_fraction or 0.0,
            failed.policy.minimum_confidence_interval_coverage_fraction,
        )
        self.assertEqual(len(failed.technical_recovery_reasons), 4)

    def test_recovery_coverage_assessment_and_reasons_are_tamper_evident(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            report = _run(source_a, source_b)

        coverage = report.model_dump(mode="json")
        coverage["confidence_interval_coverage_fraction"] = 0.0
        with self.assertRaisesRegex(ValueError, "coverage is inconsistent"):
            MethylationMixtureReport.model_validate(coverage)

        assessment = report.model_dump(mode="json")
        assessment["technical_recovery_assessment"] = "FAIL"
        with self.assertRaisesRegex(ValueError, "assessment is inconsistent"):
            MethylationMixtureReport.model_validate(assessment)

        reasons = report.model_dump(mode="json")
        reasons["technical_recovery_reasons"] = ["Synthetic tampering"]
        with self.assertRaisesRegex(ValueError, "assessment is inconsistent"):
            MethylationMixtureReport.model_validate(reasons)


class InputContractAndArtifactTests(unittest.TestCase):
    def test_operator_must_explicitly_declare_distinct_biological_sources(self) -> None:
        missing = Path("SYNTHETIC_FILE_THAT_DOES_NOT_EXIST.tsv")
        with self.assertRaisesRegex(ValueError, "must explicitly declare"):
            _run(
                missing,
                missing,
                sources_declared_biologically_distinct=False,
            )

    def test_source_ids_builds_and_byte_identity_must_be_compatible(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_a, source_b = _write_distinct_sources(root)
            with self.assertRaisesRegex(ValueError, "different source IDs"):
                _run(source_a, source_b, source_b_id="SOURCE_A")
            with self.assertRaisesRegex(ValueError, "different genome builds"):
                _run(source_a, source_b, source_b_build=GenomeBuild.GRCH37)

            identical = root / "identical.tsv"
            identical.write_bytes(source_a.read_bytes())
            with self.assertRaisesRegex(ValueError, "byte-identical"):
                _run(source_a, identical)

    def test_sequence_context_is_part_of_exact_marker_identity(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory), source_b_sequence="TCGTA")
            report = _run(source_a, source_b)

        self.assertEqual(report.reference_markers, [])
        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)

    def test_adapter_and_reference_build_provenance_cannot_contradict_input(self) -> None:
        invalid = _known_metadata().model_dump(mode="json")
        invalid["caller_name"] = "modkit"
        with self.assertRaises(ValueError):
            NanopolishSourceMetadata.model_validate(invalid)

        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            mismatched = _known_metadata(reference_genome_build="GRCh37")
            with self.assertRaisesRegex(ValueError, "reference genome build contradicts"):
                _run(source_a, source_b, source_a_metadata=mismatched)

    def test_known_metadata_mismatch_is_refused_before_source_files_are_read(self) -> None:
        metadata_a = _known_metadata(caller_version="0.13.3")
        metadata_b = _known_metadata(caller_version="0.14.0")
        missing = Path("SYNTHETIC_FILE_THAT_DOES_NOT_EXIST.tsv")
        with self.assertRaisesRegex(ValueError, "incompatible.*caller_version"):
            _run(
                missing,
                missing,
                source_a_metadata=metadata_a,
                source_b_metadata=metadata_b,
            )

    def test_unknown_metadata_warns_and_software_git_provenance_round_trips(self) -> None:
        known = _known_metadata()
        git_commit = "c" * 40
        with tempfile.TemporaryDirectory() as directory:
            source_a, source_b = _write_distinct_sources(Path(directory))
            report = _run(
                source_a,
                source_b,
                source_a_metadata=known,
                source_b_metadata=NanopolishSourceMetadata(),
                software_version="9.8.7-test",
                git_commit=git_commit,
            )

        incomplete = [
            warning for warning in report.warnings if "Upstream provenance is incomplete" in warning
        ]
        self.assertEqual(len(incomplete), 1)
        self.assertIn("source B", incomplete[0])
        self.assertEqual(report.source_a.upstream, known)
        self.assertEqual(report.source_b.upstream.caller_name, "unknown")
        self.assertEqual(report.software_version, "9.8.7-test")
        self.assertEqual(report.git_commit, git_commit)
        serialized = report.model_dump_json()
        self.assertIn('"software_version":"9.8.7-test"', serialized)
        self.assertIn(f'"git_commit":"{git_commit}"', serialized)

    def test_analysis_id_accepts_70_characters_and_rejects_71_before_parsing(self) -> None:
        accepted = "A" * 70
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_a, source_b = _write_distinct_sources(root)
            report = _run(source_a, source_b, analysis_id=accepted)
            self.assertEqual(report.analysis_id, accepted)
            self.assertTrue(all(len(level.level_id) == 80 for level in report.levels))

            missing = root / "missing.tsv"
            with self.assertRaisesRegex(ValueError, "safe pseudonymous identifier"):
                _run(missing, missing, analysis_id="A" * 71)

    def test_reports_and_exports_omit_input_paths_and_read_names(self) -> None:
        sensitive_read_token = "SYNTHETIC_RAW_READ_IDENTIFIER"
        with tempfile.TemporaryDirectory(prefix="synthetic-private-path-") as directory:
            root = Path(directory)
            source_a, source_b = _write_distinct_sources(
                root,
                read_name_token=sensitive_read_token,
            )
            report = _run(source_a, source_b)
            csv_path = render_methylation_mixture_csv(report, root / "report.csv")
            html_path = render_methylation_mixture_html(report, root / "report.html")
            serialized = report.model_dump_json()
            csv_text = csv_path.read_text(encoding="utf-8")
            html_text = html_path.read_text(encoding="utf-8")

            for artifact in (serialized, csv_text, html_text):
                self.assertNotIn(directory, artifact)
                self.assertNotIn(source_a.name, artifact)
                self.assertNotIn(source_b.name, artifact)
                self.assertNotIn(sensitive_read_token, artifact)

        self.assertIn("target_source_a_fraction", csv_text.splitlines()[0])
        self.assertEqual(len(csv_text.splitlines()), len(report.levels) + 1)
        csv_rows = list(csv.DictReader(io.StringIO(csv_text)))
        self.assertTrue(csv_rows)
        self.assertTrue(all(row["technical_recovery_assessment"] == "PASS" for row in csv_rows))
        self.assertTrue(all(row["quantifiable_fraction"] == "1.0" for row in csv_rows))
        self.assertIn("RESEARCH USE ONLY", html_text)
        self.assertIn("fraction of held-out read groups", html_text)
        self.assertIn("Technical recovery<strong>PASS</strong>", html_text)
        self.assertIn("Quantifiable<strong>100.0%</strong>", html_text)
        self.assertIn(report.algorithm_version, html_text)
        self.assertIn(report.source_a.source_id, html_text)
        self.assertIn(report.source_b.source_id, html_text)

    def test_cli_generates_three_reloadable_artifacts_from_synthetic_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_a, source_b = _write_distinct_sources(root)
            policy_path = root / "policy.yaml"
            policy_path.write_text(
                yaml.safe_dump(_policy(replicates=1).model_dump(mode="json"), sort_keys=False),
                encoding="utf-8",
                newline="\n",
            )
            output_dir = root / "outputs"
            argv = [
                "ontseq",
                "methylation-mixture",
                "--source-a-calls",
                str(source_a),
                "--source-b-calls",
                str(source_b),
                "--source-a-id",
                "CLI_SOURCE_A",
                "--source-b-id",
                "CLI_SOURCE_B",
                "--confirm-biologically-distinct-sources",
                "--source-a-genome-build",
                "GRCh38",
                "--source-b-genome-build",
                "GRCh38",
                "--analysis-id",
                "SYNTHETIC_CLI",
                "--policy",
                str(policy_path),
                "--output-dir",
                str(output_dir),
            ]
            captured = io.StringIO()
            with unittest.mock.patch("sys.argv", argv), contextlib.redirect_stdout(captured):
                cli.main()

            stem = "SYNTHETIC_CLI.methylation-mixture"
            json_path = output_dir / f"{stem}.json"
            csv_path = output_dir / f"{stem}.csv"
            html_path = output_dir / f"{stem}.html"
            for path in (json_path, csv_path, html_path):
                self.assertTrue(path.is_file(), path)
                self.assertIn(str(path), captured.getvalue())
            reloaded = MethylationMixtureReport.model_validate_json(
                json_path.read_text(encoding="utf-8")
            )

        self.assertEqual(reloaded.analysis_id, "SYNTHETIC_CLI")
        self.assertEqual(reloaded.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(reloaded.completed_levels, 5)


if __name__ == "__main__":
    unittest.main()
