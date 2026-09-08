from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from ontseq_platform.methylation_mixture import MethylationMixturePolicy
from ontseq_platform.modbam import (
    ModbamAdapterPolicy,
    ModbamAlignment,
    ModbamSourceMetadata,
    ModbamSourceSummary,
    build_modkit_extract_command,
    decode_mm_ml,
    parse_modbam_source,
    parse_modkit_source,
)
from ontseq_platform.models import GenomeBuild
from ontseq_platform.reference import sha256_file


def _sam(
    name: str = "SYN_READ_001",
    *,
    sequence: str = "ACGCGT",
    flag: int = 0,
    position: int = 1,
    cigar: str = "6M",
    mm: str = "C+m?,0,0;",
    ml: str = "255,0",
    extra: str = "",
    chromosome: str = "chr1",
) -> str:
    return (
        f"{name}\t{flag}\t{chromosome}\t{position}\t60\t{cigar}\t*\t0\t0\t{sequence}\t*"
        f"\tMM:Z:{mm}\tML:B:C{',' + ml if ml else ''}{extra}\n"
    )


def _row(
    name: str = "SYN_READ_001",
    *,
    position: int = 1,
    query: int = 1,
    probability: str = "0.998046875",
    code: str = "m",
    flag: int = 0,
    inferred: str = "false",
    chromosome: str = "chr1",
) -> dict[str, str]:
    strand = "-" if flag & 16 else "+"
    return {
        "read_id": name,
        "forward_read_position": str(query),
        "ref_position": str(position),
        "chrom": chromosome,
        "mod_strand": "+",
        "ref_strand": strand,
        "ref_mod_strand": strand,
        "read_length": "6",
        "mod_qual": probability,
        "mod_code": code,
        "canonical_base": "C",
        "modified_primary_base": "C",
        "inferred": inferred,
        "flag": str(flag),
        "alignment_start": "0",
        "alignment_end": "6",
    }


class ModbamAdapterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.reference = self.directory / "synthetic.fa"
        self.reference.write_bytes(b">chr1 SYNTHETIC ONLY\nACGCGT\nACGCGT\n")
        self.policy = MethylationMixturePolicy(profile_id="synthetic-only", note="Synthetic test")
        self.upstream = ModbamSourceMetadata(
            caller_version="synthetic-1",
            reference_id="synthetic-reference-v1",
            reference_genome_build=GenomeBuild.GRCH38,
            reference_sha256=sha256_file(self.reference),
            sequencing_platform="ONT-synthetic",
            flow_cell_product_code="synthetic",
            library_kit="synthetic",
            basecaller_version="synthetic-1",
            basecaller_model="synthetic-cpg-v1",
            modification_model="synthetic-5mC-v1",
            cytosine_modifications="5mC",
            modkit_version="0.6.1",
            source_modbam_sha256="a" * 64,
        )

    def sam(self, rows: str, *, header: str | None = None, **kwargs: object):  # type: ignore[no-untyped-def]
        path = self.directory / "synthetic.sam"
        header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:12\n" if header is None else header
        path.write_text(header + rows, encoding="ascii")
        parameters = dict(
            reference_path=self.reference,
            source_id="SYN_SOURCE",
            genome_build=GenomeBuild.GRCH38,
            policy=self.policy,
            upstream=self.upstream,
        )
        parameters.update(kwargs)
        return parse_modbam_source(path, **parameters)  # type: ignore[arg-type]

    def tsv(self, rows: list[dict[str, str]], **kwargs: object):  # type: ignore[no-untyped-def]
        path = self.directory / "synthetic.tsv"
        with path.open("w", newline="", encoding="ascii") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        parameters = dict(
            reference_path=self.reference,
            source_id="SYN_SOURCE",
            genome_build=GenomeBuild.GRCH38,
            policy=self.policy,
            upstream=self.upstream,
        )
        parameters.update(kwargs)
        return parse_modkit_source(path, **parameters)  # type: ignore[arg-type]

    def test_sam_coordinates_and_probability_classes(self) -> None:
        source = self.sam(_sam())
        calls = source.calls_by_read["SYN_READ_001"]
        self.assertEqual([call.marker.start for call in calls], [1, 3])
        self.assertEqual([call.marker.end for call in calls], [1, 3])
        self.assertEqual([call.state for call in calls], ["methylated", "unmethylated"])
        self.assertEqual(calls[0].marker.sequence_sha256, hashlib.sha256(b"CG").hexdigest())
        self.assertEqual(source.summary.input_format, "sam_mm_ml")
        self.assertEqual(source.summary.parser_version, "0.1.1")
        self.assertEqual(source.summary.informative_call_rows, 2)
        self.assertEqual(
            source.summary.reference_fingerprint.sha256, self.upstream.reference_sha256
        )
        self.assertNotIn("SYN_READ_001", source.summary.model_dump_json())

    def test_reverse_sam_mm_deltas_use_original_sequence(self) -> None:
        # SAM SEQ is still ACGCGT (palindrome), but original C positions map to G.
        source = self.sam(_sam(flag=16))
        calls = source.calls_by_read["SYN_READ_001"]
        self.assertEqual(
            [(c.marker.start, c.marker.strand, c.state) for c in calls],
            [
                (4, "-", "methylated"),
                (2, "-", "unmethylated"),
            ],
        )

    def test_reverse_nonpalindromic_sequence_and_cigar(self) -> None:
        record = ModbamAlignment("SYN", "chr1", 16, 10, "GACGT", "1S2M1I1M", "C+m?,0;", (255,))
        forward, mapping, events = decode_mm_ml(record, expected_codes=frozenset({"m"}))
        self.assertEqual(forward, "ACGTC")
        self.assertEqual(mapping, [12, None, 11, 10, None])
        self.assertEqual(set(events), {1})

    def test_cigar_insertion_and_deletion_map_correctly(self) -> None:
        record = ModbamAlignment("SYN", "chr1", 0, 20, "ACGCGT", "1S2M1I1D2M", "C+m,0,0;", (255, 0))
        _, mapping, _ = decode_mm_ml(record, expected_codes=frozenset({"m"}))
        self.assertEqual(mapping, [None, 20, 21, None, 23, 24])

    def test_reference_mismatch_and_non_cpg_are_not_counted(self) -> None:
        source = self.sam(_sam(sequence="CCGCGT", mm="C+m?,0,0,0;", ml="255,255,0"))
        self.assertEqual(source.summary.canonical_call_rows, 2)
        self.assertEqual(source.summary.excluded_counts["non_cpg_or_reference_mismatch"], 1)

    def test_implicit_skips_are_ambiguous_even_dot(self) -> None:
        for mm in ("C+m.,0;", "C+m?,0;", "C+m,0;"):
            with self.subTest(mm=mm):
                source = self.sam(_sam(mm=mm, ml="255"))
                self.assertEqual(
                    [c.state for c in source.calls_by_read["SYN_READ_001"]],
                    [
                        "methylated",
                        "ambiguous",
                    ],
                )

    def test_empty_mm_events_are_unknown_not_canonical(self) -> None:
        source = self.sam(_sam(mm="C+m;", ml=""))
        self.assertEqual(source.summary.ambiguous_call_rows, 2)
        self.assertEqual(source.summary.informative_call_rows, 0)

    def test_high_hydroxymethylation_is_never_canonical(self) -> None:
        upstream = self.upstream.model_copy(update={"cytosine_modifications": "5mC_5hmC"})
        source = self.sam(_sam(mm="C+mh?,0,0;", ml="0,255,0,0"), upstream=upstream)
        self.assertEqual(
            [c.state for c in source.calls_by_read["SYN_READ_001"]],
            [
                "ambiguous",
                "unmethylated",
            ],
        )

    def test_separate_m_h_groups_have_same_result_as_interleaved(self) -> None:
        upstream = self.upstream.model_copy(update={"cytosine_modifications": "5mC_5hmC"})
        source = self.sam(_sam(mm="C+m?,0,0;C+h?,0,0;", ml="0,0,255,0"), upstream=upstream)
        self.assertEqual(
            [c.state for c in source.calls_by_read["SYN_READ_001"]],
            [
                "ambiguous",
                "unmethylated",
            ],
        )

    def test_declared_model_must_match_all_cytosine_tags(self) -> None:
        for mm, ml in (("C+mh,0;", "0,255"), ("C+h,0;", "255"), ("C+76792,0;", "255")):
            with self.subTest(mm=mm), self.assertRaisesRegex(ValueError, "declared model"):
                self.sam(_sam(mm=mm, ml=ml))

    def test_other_fundamental_base_probabilities_are_consumed(self) -> None:
        source = self.sam(_sam(mm="A+a,0;C+m,0,0;", ml="255,255,0"))
        self.assertEqual(source.summary.informative_call_rows, 2)

    def test_bad_mm_ml_and_opposite_mod_strand_fail_closed(self) -> None:
        cases = [
            ("C+m,4;", "255", "beyond"),
            ("C+m,0,0;", "255", "count"),
            ("C+m,0;", "255,0", "count"),
            ("C+m,0;", "256", "unsigned byte"),
            ("C-m,0;", "255", "Opposite"),
            ("C+m,0;C+m,0;", "255,0", "repeats"),
            ("C+m,-1;", "255", "grammar"),
            ("C+m,0", "255", "grammar"),
        ]
        for mm, ml, error in cases:
            with self.subTest(mm=mm), self.assertRaisesRegex(ValueError, error):
                self.sam(_sam(mm=mm, ml=ml))

    def test_missing_or_wrong_type_tags_fail_closed(self) -> None:
        for line in [
            _sam().replace("\tMM:Z:C+m?,0,0;", ""),
            _sam().replace("\tML:B:C,255,0", ""),
            _sam().replace("MM:Z:", "MM:i:"),
            _sam().replace("ML:B:C", "ML:B:i"),
            _sam().replace("MM:Z:", "Mm:Z:"),
        ]:
            with self.subTest(line=line), self.assertRaises(ValueError):
                self.sam(line)

    def test_stale_mn_and_unsafe_clipping_fail(self) -> None:
        for cigar, extra in (("6M", "\tMN:i:7"), ("2H6M", ""), ("2M1S3M", "\tMN:i:6")):
            with self.subTest(cigar=cigar), self.assertRaises(ValueError):
                self.sam(_sam(cigar=cigar, extra=extra))
        self.assertEqual(
            self.sam(_sam(cigar="2H6M", extra="\tMN:i:6")).summary.canonical_call_rows, 2
        )

    def test_flagged_alignments_excluded_before_missing_tags(self) -> None:
        excluded = "".join(
            _sam(f"SYN_EXCLUDED_{flag}", flag=flag).split("\tMM:")[0] + "\n"
            for flag in (4, 256, 512, 1024, 2048)
        )
        source = self.sam(_sam() + excluded)
        self.assertEqual(source.summary.read_group_count, 1)
        self.assertEqual(source.summary.excluded_counts["alignment_flags"], 5)

    def test_repeated_primary_alignments_and_paired_reads_fail(self) -> None:
        with self.assertRaisesRegex(ValueError, "duplicate primary"):
            self.sam(_sam() + _sam())
        with self.assertRaisesRegex(ValueError, "Paired"):
            self.sam(_sam(flag=1))

    def test_dorado_duplex_tag_rejected_without_program_header(self) -> None:
        # Dorado duplex offspring must not enter the simplex read denominator,
        # even after an upstream tool removes or renames the program header.
        for flag in (0, 16):
            with self.subTest(flag=flag), self.assertRaisesRegex(ValueError, "Duplex reads"):
                self.sam(_sam(flag=flag, extra="\tdx:i:1"))

    def test_dorado_simplex_tags_remain_accepted(self) -> None:
        baseline = self.sam(_sam())
        for dx in (-1, 0):
            with self.subTest(dx=dx):
                source = self.sam(_sam(extra=f"\tdx:i:{dx}"))
                self.assertEqual(source.calls_by_read, baseline.calls_by_read)
        # A duplicate-flagged duplex record remains an exclusion.
        source = self.sam(_sam() + _sam("SYN_DUPLEX", flag=1024, extra="\tdx:i:1"))
        self.assertEqual(source.summary.read_group_count, 1)
        self.assertEqual(source.summary.excluded_counts, {"alignment_flags": 1})

    def test_malformed_or_repeated_dorado_duplex_tag_rejected(self) -> None:
        for extra in ("\tdx:Z:1", "\tdx:i:2", "\tdx:i:-2", "\tdx:i:0\tdx:i:1"):
            with self.subTest(extra=extra), self.assertRaisesRegex(ValueError, "dx"):
                self.sam(_sam(extra=extra))

    def test_sam_and_modkit_explicit_probabilities_agree(self) -> None:
        sam = self.sam(_sam())
        tsv = self.tsv([_row(), _row(position=3, query=3, probability="0.001953125")])
        self.assertEqual(sam.calls_by_read, tsv.calls_by_read)
        self.assertEqual(tsv.summary.input_format, "modkit_extract_full_tsv")
        self.assertEqual(tsv.summary.parser_version, "0.1.1")

    def test_reverse_sam_and_modkit_agree(self) -> None:
        sam = self.sam(_sam(flag=16))
        tsv = self.tsv(
            [
                _row(position=4, query=1, flag=16),
                _row(position=2, query=3, flag=16, probability="0.001953125"),
            ]
        )
        self.assertEqual(sam.calls_by_read, tsv.calls_by_read)

    def test_modkit_hm_canonical_requires_both_classes(self) -> None:
        upstream = self.upstream.model_copy(update={"cytosine_modifications": "5mC_5hmC"})
        source = self.tsv(
            [
                _row(probability="0.001953125"),
                _row(code="h"),
                _row(position=3, query=3, probability="0.001953125"),
            ],
            upstream=upstream,
        )
        self.assertEqual(
            [c.state for c in source.calls_by_read["SYN_READ_001"]],
            [
                "ambiguous",
                "ambiguous",
            ],
        )

    def test_modkit_implicit_probabilities_are_ambiguous(self) -> None:
        source = self.tsv([_row(probability="0", inferred="true")])
        self.assertEqual(source.summary.ambiguous_call_rows, 1)

    def test_modkit_probability_domain_and_joint_sum(self) -> None:
        for value in ("255", "nan", "inf", "-0.1", "text"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "probability"):
                self.tsv([_row(probability=value)])
        upstream = self.upstream.model_copy(update={"cytosine_modifications": "5mC_5hmC"})
        with self.assertRaisesRegex(ValueError, "total probability"):
            self.tsv([_row(), _row(code="h")], upstream=upstream)

    def test_modkit_probability_must_be_an_untransformed_ml_midpoint(self) -> None:
        for value in ("0", "1", "0.9", "0.5"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "ML-bin midpoint"):
                self.tsv([_row(probability=value)])
        # Actual Rust f32 decimal formatting rounds the binary midpoint.
        source = self.tsv([_row(probability="0.9980469")])
        self.assertEqual(source.summary.informative_call_rows, 1)

    def test_sam_tsv_agree_at_exact_probability_threshold(self) -> None:
        adapter_policy = ModbamAdapterPolicy(minimum_class_probability=0.75)
        sam = self.sam(_sam(ml="192,63"), adapter_policy=adapter_policy)
        tsv = self.tsv(
            [
                _row(probability="0.7519531"),
                _row(position=3, query=3, probability="0.24804688"),
            ],
            adapter_policy=adapter_policy,
        )
        self.assertEqual(sam.calls_by_read, tsv.calls_by_read)
        self.assertEqual(tsv.summary.informative_call_rows, 2)

    def test_modkit_alignment_end_is_exclusive(self) -> None:
        with self.assertRaisesRegex(ValueError, "outside its alignment"):
            self.tsv([{**_row(), "alignment_end": "1"}])

    def test_modkit_alignment_span_must_fit_locked_reference(self) -> None:
        # The reported CpG remains valid; the contradictory alignment extent
        # must still invalidate the source rather than passing context checks.
        with self.assertRaisesRegex(ValueError, "beyond the locked reference"):
            self.tsv([{**_row(), "alignment_end": "13"}])
        with self.assertRaisesRegex(ValueError, "absent from the locked reference"):
            self.tsv([_row(chromosome="chr2")])
        source = self.tsv(
            [
                {
                    **_row(position=9, query=3),
                    "alignment_start": "6",
                    "alignment_end": "12",
                }
            ]
        )
        self.assertEqual(source.summary.informative_call_rows, 1)

    def test_modkit_release_is_locked_in_metadata_and_parser(self) -> None:
        for version in ("0.4.1", "0.5.0", "0.6.0", "0.6.4", "unknown"):
            with self.subTest(version=version):
                metadata = self.upstream.model_dump()
                metadata["modkit_version"] = version
                with self.assertRaises(ValueError):
                    ModbamSourceMetadata.model_validate(metadata)
                with self.assertRaisesRegex(ValueError, "supports only"):
                    self.tsv(
                        [_row()],
                        upstream=self.upstream.model_copy(
                            update={
                                "modkit_version": version,
                            }
                        ),
                    )

    def test_probability_gate_uses_conservative_ml_bin_bounds(self) -> None:
        # N=204 gives [0.796875, 0.80078125]; its midpoint must not pass 0.8.
        source = self.sam(_sam(ml="204,51"))
        self.assertEqual(source.summary.ambiguous_call_rows, 2)
        source = self.tsv([_row(probability="0.798828125")])
        self.assertEqual(source.summary.ambiguous_call_rows, 1)

    def test_modkit_duplicate_and_contradictory_sites_rejected(self) -> None:
        for rows in [
            [_row(), _row()],
            [_row(), _row(query=2)],
            [_row(), _row(position=3)],
            [_row(), {**_row(code="h"), "alignment_end": "5"}],
        ]:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                self.tsv(rows)

    def test_modkit_strand_and_reference_guards(self) -> None:
        for field, value in (
            ("ref_strand", "-"),
            ("ref_mod_strand", "-"),
            ("mod_strand", "-"),
            ("alignment_start", "2"),
            ("forward_read_position", "6"),
            ("inferred", "maybe"),
        ):
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.tsv([{**_row(), field: value}])

    def test_modkit_flagged_and_unmapped_rows_excluded(self) -> None:
        source = self.tsv(
            [_row(), _row("SYN_DUPLICATE", flag=1024), _row("SYN_UNALIGNED", position=-1)]
        )
        self.assertEqual(source.summary.read_group_count, 1)
        self.assertEqual(
            source.summary.excluded_counts, {"alignment_flags": 1, "unaligned_base": 1}
        )

    def test_modkit_requires_full_schema_and_source_provenance(self) -> None:
        row = _row()
        del row["mod_qual"]
        row["call_prob"] = "0.99"
        with self.assertRaisesRegex(ValueError, "extract full"):
            self.tsv([row])
        for updates in ({"modkit_version": None}, {"source_modbam_sha256": None}):
            with self.subTest(updates=updates), self.assertRaisesRegex(ValueError, "version.*SHA"):
                self.tsv([_row()], upstream=self.upstream.model_copy(update=updates))

    def test_reviewed_full_header_with_and_without_optional_motifs(self) -> None:
        # The names/order are from the reviewed 0.6.1 writer; all values are synthetic.
        full_columns = (
            "read_id",
            "forward_read_position",
            "ref_position",
            "chrom",
            "mod_strand",
            "ref_strand",
            "ref_mod_strand",
            "fw_soft_clipped_start",
            "fw_soft_clipped_end",
            "alignment_start",
            "alignment_end",
            "read_length",
            "mod_qual",
            "mod_code",
            "base_qual",
            "ref_kmer",
            "query_kmer",
            "canonical_base",
            "modified_primary_base",
            "inferred",
            "flag",
        )
        values = {
            **_row(),
            "fw_soft_clipped_start": "0",
            "fw_soft_clipped_end": "0",
            "base_qual": "30",
            "ref_kmer": "-ACGC",
            "query_kmer": "-ACGC",
        }
        full_row = {column: values[column] for column in full_columns}
        for row in (full_row, {**full_row, "motifs": "."}):
            with self.subTest(has_motifs="motifs" in row):
                self.assertEqual(self.tsv([row]).summary.informative_call_rows, 1)
        legacy = dict(full_row)
        del legacy["alignment_start"]
        del legacy["alignment_end"]
        with self.assertRaisesRegex(ValueError, "extract full"):
            self.tsv([legacy])

    def test_input_read_and_decompressed_limits(self) -> None:
        with self.assertRaisesRegex(ValueError, "row safety"):
            self.sam(_sam(), policy=self.policy.model_copy(update={"maximum_rows_per_source": 1}))
        with self.assertRaisesRegex(ValueError, "per-read"):
            self.sam(_sam(), adapter_policy=ModbamAdapterPolicy(maximum_bases_per_read=5))
        with self.assertRaisesRegex(ValueError, "byte safety"):
            self.sam(
                _sam(), policy=self.policy.model_copy(update={"maximum_input_bytes_per_source": 10})
            )
        path = self.directory / "synthetic.sam.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("@CO\t" + "A" * 1000 + "\n" + _sam())
        with self.assertRaisesRegex(ValueError, "decompressed"):
            parse_modbam_source(
                path,
                reference_path=self.reference,
                source_id="SYN_SOURCE",
                genome_build=GenomeBuild.GRCH38,
                upstream=self.upstream,
                policy=self.policy.model_copy(update={"maximum_input_bytes_per_source": 500}),
            )

    def test_source_and_reference_hash_drift_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Reference SHA"):
            self.sam(
                _sam(), upstream=self.upstream.model_copy(update={"reference_sha256": "0" * 64})
            )
        with (
            patch(
                "ontseq_platform.modbam.sha256_file",
                side_effect=["a" * 64, self.upstream.reference_sha256, "b" * 64],
            ),
            self.assertRaisesRegex(ValueError, "changed"),
        ):
            self.sam(_sam())

    def test_reference_index_handles_multiline_and_terminal_short_line(self) -> None:
        source = self.sam(_sam(position=7))
        self.assertEqual([c.marker.start for c in source.calls_by_read["SYN_READ_001"]], [7, 9])
        self.reference.write_bytes(b">chr1\r\nACGC\r\nGTAC\r\nGC\r\n")
        upstream = self.upstream.model_copy(
            update={"reference_sha256": sha256_file(self.reference)}
        )
        source = self.sam(_sam(), header="@SQ\tSN:chr1\tLN:10\n", upstream=upstream)
        self.assertEqual(source.summary.canonical_call_rows, 2)

    def test_reference_duplicate_aliases_rejected(self) -> None:
        self.reference.write_bytes(b">chr1\nACGCGT\n>1\nACGCGT\n")
        upstream = self.upstream.model_copy(
            update={"reference_sha256": sha256_file(self.reference)}
        )
        with self.assertRaisesRegex(ValueError, "ambiguous chromosome aliases"):
            self.sam(_sam(), upstream=upstream)

    def test_alignment_dictionary_must_match_reference_and_be_present(self) -> None:
        for header in (
            "",
            "@SQ\tSN:chr1\tLN:11\n",
            "@SQ\tSN:chr2\tLN:12\n",
            "@SQ\tSN:chr1\tLN:12\n@SQ\tSN:chr1\tLN:12\n",
        ):
            with self.subTest(header=header), self.assertRaisesRegex(ValueError, "@SQ"):
                self.sam(_sam(), header=header)

    def test_alignment_cigar_must_stay_inside_reference(self) -> None:
        with self.assertRaisesRegex(ValueError, "beyond the locked reference"):
            self.sam(_sam(position=10))

    def test_alignment_contig_must_be_declared_in_its_header(self) -> None:
        self.reference.write_bytes(b">chr1\nACGCGTACGCGT\n>chr2\nACGCGTACGCGT\n")
        upstream = self.upstream.model_copy(
            update={"reference_sha256": sha256_file(self.reference)}
        )
        with self.assertRaisesRegex(ValueError, "absent from its @SQ"):
            self.sam(_sam(chromosome="chr2"), upstream=upstream)

    def test_dorado_pg_version_and_rg_models_cross_checked_when_present(self) -> None:
        header = (
            "@SQ\tSN:chr1\tLN:12\n"
            "@PG\tID:basecaller\tPN:dorado\tVN:synthetic-1\n"
            "@RG\tID:SYN_GROUP\tPL:ONT\tDS:basecall_model=synthetic-cpg-v1 "
            "modbase_models=synthetic-5mC-v1\n"
        )
        source = self.sam(_sam(), header=header)
        self.assertEqual(
            set(source.summary.header_verified_fields),
            {
                "reference_sequence_dictionary",
                "caller_version",
                "basecaller_version",
                "basecaller_model",
                "modification_model",
            },
        )
        for altered in (
            header.replace("VN:synthetic-1", "VN:synthetic-2"),
            header.replace("PN:dorado", "PN:guppy"),
            header.replace("basecall_model=synthetic-cpg-v1", "basecall_model=other-model"),
            header.replace("modbase_models=synthetic-5mC-v1", "modbase_models=other-model"),
            header.replace("PL:ONT", "PL:ILLUMINA"),
        ):
            with self.subTest(header=altered), self.assertRaises(ValueError):
                self.sam(_sam(), header=altered)

    def test_missing_dorado_headers_remain_declared_not_verified(self) -> None:
        source = self.sam(_sam())
        self.assertEqual(source.summary.header_verified_fields, ["reference_sequence_dictionary"])
        self.assertEqual(
            source.summary.provenance_verification, "declared_with_available_header_cross_checks"
        )

    def test_actual_dorado_model_names_with_at_version_are_accepted(self) -> None:
        metadata = self.upstream.model_dump()
        metadata["basecaller_model"] = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
        metadata["modification_model"] = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v1"
        parsed = ModbamSourceMetadata.model_validate(metadata)
        self.assertIn("@v", parsed.basecaller_model)

    def test_summary_counts_and_fingerprints_validated(self) -> None:
        source = self.sam(_sam())
        for field, value in (
            ("ambiguous_call_rows", 2),
            ("reference_fingerprint", {"size_bytes": 1, "sha256": "0" * 64}),
        ):
            values = source.summary.model_dump()
            values[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                ModbamSourceSummary.model_validate(values)

    def test_bam_optional_dependency_message_is_explicit(self) -> None:
        path = self.directory / "synthetic.bam"
        path.write_bytes(b"synthetic placeholder, never a genome")
        with (
            patch("ontseq_platform.modbam.importlib.import_module", side_effect=ImportError),
            self.assertRaisesRegex(ValueError, "optional pysam"),
        ):
            parse_modbam_source(
                path,
                reference_path=self.reference,
                source_id="SYN_SOURCE",
                genome_build=GenomeBuild.GRCH38,
                upstream=self.upstream,
                policy=self.policy,
            )

    @unittest.skipUnless(importlib.util.find_spec("pysam"), "optional pysam is not installed")
    def test_real_bam_container_matches_sam_decoder(self) -> None:
        import pysam

        sam_source = self.sam(_sam() + _sam("SYN_READ_REVERSE", flag=16))
        bam_path = self.directory / "synthetic.bam"
        with (
            pysam.AlignmentFile(str(self.directory / "synthetic.sam"), "r") as sam,
            pysam.AlignmentFile(str(bam_path), "wb", template=sam) as bam,
        ):
            for record in sam:
                bam.write(record)
        bam_source = parse_modbam_source(
            bam_path,
            reference_path=self.reference,
            source_id="SYN_SOURCE",
            genome_build=GenomeBuild.GRCH38,
            upstream=self.upstream,
            policy=self.policy,
        )
        self.assertEqual(bam_source.calls_by_read, sam_source.calls_by_read)
        self.assertEqual(bam_source.summary.input_format, "modbam_mm_ml")
        self.assertEqual(bam_source.summary.parser_version, f"{pysam.__version__}/0.1.1")

    def test_modkit_command_is_version_bound_and_keeps_path_spaces(self) -> None:
        path = self.directory / "synthetic input.bam"
        path.write_bytes(b"synthetic command fixture")
        output = self.directory / "synthetic output.tsv"
        plan = build_modkit_extract_command(
            modbam_path=path,
            reference_path=self.reference,
            output_path=output,
            expected_version="0.6.1",
            threads=3,
        )
        self.assertEqual(plan.expected_version, "0.6.1")
        self.assertEqual(plan.version_check_argv, ("modkit", "--version"))
        self.assertEqual(plan.argv[:3], ("modkit", "extract", "full"))
        self.assertIn(str(path.resolve()), plan.argv)
        self.assertIn(str(output.resolve()), plan.argv)
        self.assertIn("--mapped-only", plan.argv)
        self.assertIn("--cpg", plan.argv)
        self.assertNotIn("--ignore", plan.argv)
        self.assertNotIn("--force", plan.argv)
        self.assertFalse(output.exists())

    def test_modkit_command_refuses_existing_output_or_unknown_version(self) -> None:
        path = self.directory / "synthetic.bam"
        path.write_bytes(b"synthetic command fixture")
        with self.assertRaisesRegex(ValueError, "overwrite"):
            build_modkit_extract_command(
                modbam_path=path,
                reference_path=self.reference,
                output_path=path,
                expected_version="0.6.1",
            )
        with self.assertRaises(ValueError):
            build_modkit_extract_command(
                modbam_path=path,
                reference_path=self.reference,
                output_path=self.directory / "new.tsv",
                expected_version="unknown",
            )
        with self.assertRaisesRegex(ValueError, "supports only"):
            build_modkit_extract_command(
                modbam_path=path,
                reference_path=self.reference,
                output_path=self.directory / "legacy.tsv",
                expected_version="0.4.1",
            )


class DecodeSafetyTests(unittest.TestCase):
    def test_standalone_decoder_enforces_simplex_dorado_status(self) -> None:
        base = ModbamAlignment("SYN", "chr1", 0, 0, "ACG", "3M", "C+m,0;", (255,))
        for dx in (-2, 1, 2, True):
            with self.subTest(dx=dx), self.assertRaises(ValueError):
                decode_mm_ml(replace(base, dx=dx), expected_codes=frozenset({"m"}))

    def test_standalone_decoder_validates_lengths(self) -> None:
        base = ModbamAlignment("SYN", "chr1", 0, 0, "ACG", "3M", "C+m,0;", (255,))
        for record in [
            replace(base, cigar="2M"),
            replace(base, cigar="0M3M"),
            replace(base, cigar="2M1P1M"),
            replace(base, cigar="999999999M"),
            replace(base, sequence=""),
            replace(base, reference_start=-1),
        ]:
            with self.subTest(record=record), self.assertRaises(ValueError):
                decode_mm_ml(record, expected_codes=frozenset({"m"}))

    def test_unknown_ml_value_cannot_be_a_canonical_call(self) -> None:
        record = ModbamAlignment("SYN", "chr1", 0, 0, "ACG", "3M", "C+m,0;", (-1,))
        with self.assertRaisesRegex(ValueError, "unsigned byte"):
            decode_mm_ml(record, expected_codes=frozenset({"m"}))


if __name__ == "__main__":
    unittest.main()
