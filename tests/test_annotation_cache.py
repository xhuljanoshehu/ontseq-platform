from __future__ import annotations

import gzip
import sqlite3
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.annotation_cache import (
    AnnotationCache,
    compile_annotation_cache,
    validate_annotation_cache,
)

FIXTURES = Path(__file__).parent / "fixtures" / "reference_cache"
GENCODE = FIXTURES / "gencode.v50.fragment.gtf"
MANE = FIXTURES / "MANE.GRCh38.v1.5.fragment.gff3"
CYTOBANDS = FIXTURES / "cytoBand.hg38.fragment.tsv"


def _compile(output: Path, *, gencode: Path = GENCODE) -> None:
    compile_annotation_cache(
        gencode,
        MANE,
        CYTOBANDS,
        output,
        metadata={
            "bundle_id": "GRCh38_GENCODE50_MANE1.5_v1_fixture",
            "bundle_version": "fixture-v1",
            "genome_build": "GRCh38",
            "gencode_release": "50-fixture",
            "mane_release": "1.5-fixture",
            "cytoband_release": "hg38-fixture",
        },
    )


class AnnotationCacheCompilerTests(unittest.TestCase):
    def test_gtf_gff3_and_cytobands_compile_to_zero_based_sqlite(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            database = Path(raw) / "annotation.sqlite"
            summary = compile_annotation_cache(
                GENCODE,
                MANE,
                CYTOBANDS,
                database,
                metadata={"genome_build": "GRCh38", "bundle_id": "fixture"},
            )
            self.assertEqual(summary.genes, 2)
            self.assertEqual(summary.transcripts, 3)
            self.assertEqual(summary.exons, 5)
            self.assertEqual(summary.cds, 3)
            self.assertEqual(summary.cytobands, 4)
            self.assertEqual(summary.mane_matched_transcripts, 2)

            connection = sqlite3.connect(database)
            try:
                bcr = connection.execute(
                    "SELECT chrom, start, end FROM genes WHERE gene_name = 'BCR'"
                ).fetchone()
                minus_cds = connection.execute(
                    """SELECT strand, start, end, phase FROM cds
                       WHERE transcript_id = 'ENST00000318560.7'"""
                ).fetchone()
                band = connection.execute(
                    "SELECT start, end, name, gie_stain FROM cytobands WHERE name = 'q11.23'"
                ).fetchone()
            finally:
                connection.close()

            self.assertEqual(bcr, ("chr22", 23280183, 23318037))
            self.assertEqual(minus_cds, ("-", 130854063, 130854351, 2))
            self.assertEqual(band, (22600000, 23500000, "q11.23", "gneg"))

    def test_mane_then_canonical_transcript_ranking_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            database = Path(raw) / "annotation.sqlite"
            _compile(database)
            transcripts = AnnotationCache(database).ranked_transcripts("BCR")
            by_gene_id = AnnotationCache(database).ranked_transcripts("ENSG00000186716.22")

        self.assertEqual(
            [item.transcript_id for item in transcripts],
            ["ENST00000305877.13", "ENST00000483320.6"],
        )
        self.assertEqual(transcripts[0].mane_status, "MANE Select")
        self.assertEqual(transcripts[0].mane_refseq_id, "NM_004327.4")
        self.assertEqual(transcripts[1].appris, "appris_principal_1")
        self.assertTrue(transcripts[1].is_canonical)
        self.assertEqual(by_gene_id, transcripts)

    def test_unversioned_mane_id_matches_one_versioned_gencode_transcript(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            database = Path(raw) / "annotation.sqlite"
            _compile(database)
            abl1 = AnnotationCache(database).ranked_transcripts("ABL1")

        self.assertEqual(len(abl1), 1)
        self.assertEqual(abl1[0].mane_status, "MANE Plus Clinical")
        self.assertEqual(abl1[0].mane_refseq_id, "NM_005157.6")

    def test_compilation_is_byte_deterministic_and_accepts_gzip_gtf(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            compressed = root / "fragment.gtf.gz"
            with gzip.open(compressed, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write(GENCODE.read_text(encoding="utf-8"))
            first = root / "first.sqlite"
            second = root / "second.sqlite"
            _compile(first, gencode=compressed)
            _compile(second, gencode=compressed)
            first_summary = validate_annotation_cache(first)
            second_summary = validate_annotation_cache(second)

        self.assertEqual(first_summary.sha256, second_summary.sha256)

    def test_unnamed_ucsc_chromosome_placeholder_is_not_a_cytoband(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            cytobands = root / "cytobands.tsv"
            cytobands.write_text(
                CYTOBANDS.read_text(encoding="utf-8") + "chrM\t0\t16569\t\tgneg\n",
                encoding="utf-8",
            )
            database = root / "annotation.sqlite"
            summary = compile_annotation_cache(
                GENCODE,
                MANE,
                cytobands,
                database,
                metadata={"genome_build": "GRCh38", "bundle_id": "fixture"},
            )

        self.assertEqual(summary.cytobands, 4)

    def test_wrong_build_is_refused_before_cache_activation(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            output = Path(raw) / "annotation.sqlite"
            with self.assertRaisesRegex(ValueError, "only genome_build=GRCh38"):
                compile_annotation_cache(
                    GENCODE,
                    MANE,
                    CYTOBANDS,
                    output,
                    metadata={"genome_build": "GRCh37"},
                )
            self.assertFalse(output.exists())

    def test_malformed_gtf_does_not_replace_an_existing_cache(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            database = root / "annotation.sqlite"
            _compile(database)
            before = database.read_bytes()
            malformed = root / "malformed.gtf"
            malformed.write_text("chr22\tbad\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "expected nine GTF columns"):
                _compile(database, gencode=malformed)

            self.assertEqual(database.read_bytes(), before)
            self.assertFalse(any(root.glob(".annotation.sqlite.*.tmp")))

    def test_mane_locus_disagreement_with_gencode_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            mismatched_mane = root / "mane.gff3"
            mismatched_mane.write_text(
                MANE.read_text(encoding="utf-8").replace(
                    "23280184\t23318037", "23280185\t23318037", 1
                ),
                encoding="utf-8",
            )
            output = root / "annotation.sqlite"
            with self.assertRaisesRegex(ValueError, "locus does not match GENCODE"):
                compile_annotation_cache(
                    GENCODE,
                    mismatched_mane,
                    CYTOBANDS,
                    output,
                    metadata={"genome_build": "GRCh38"},
                )

            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
