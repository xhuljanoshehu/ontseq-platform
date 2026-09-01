from __future__ import annotations

import sqlite3
import unittest
from contextlib import closing
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml

from ontseq_platform.models import PanelBundle
from ontseq_platform.panel_bundle import sha256_file
from ontseq_platform.panel_compiler import (
    PanelCompilerError,
    compile_panel_derivatives,
    materialize_and_pin_panel_derivatives,
)


def _annotation_cache(path: Path, *, genome_build: str = "GRCh38") -> None:
    with closing(sqlite3.connect(path)) as connection:
        connection.executescript(
            """
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE genes (
              gene_id TEXT PRIMARY KEY, gene_name TEXT NOT NULL, chrom TEXT NOT NULL,
              start INTEGER NOT NULL, end INTEGER NOT NULL, strand TEXT NOT NULL,
              gene_type TEXT
            );
            CREATE TABLE transcripts (
              transcript_id TEXT PRIMARY KEY, gene_id TEXT NOT NULL, transcript_name TEXT,
              chrom TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL,
              strand TEXT NOT NULL, transcript_type TEXT, tags_json TEXT, mane_status TEXT,
              mane_refseq_id TEXT, appris TEXT, is_canonical INTEGER NOT NULL,
              is_basic INTEGER NOT NULL, cds_length INTEGER NOT NULL,
              transcript_length INTEGER NOT NULL
            );
            """
        )
        connection.execute("INSERT INTO metadata VALUES ('genome_build', ?)", (genome_build,))
        connection.executemany(
            "INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?)",
            [
                ("ENSG_BCR", "BCR", "chr22", 23_180_000, 23_318_000, "+", "protein_coding"),
                ("ENSG_DUP1", "DUP", "chr1", 100, 200, "+", "protein_coding"),
                ("ENSG_DUP2", "DUP", "chr1", 300, 400, "+", "protein_coding"),
                ("ENSG_P2RY8_X", "P2RY8", "chrX", 1_462_571, 1_537_207, "-", "protein_coding"),
                ("ENSG_P2RY8_Y", "P2RY8", "chrY", 1_462_571, 1_537_207, "-", "protein_coding"),
            ],
        )
        connection.executemany(
            """
            INSERT INTO transcripts VALUES
            (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    "ENST_OTHER",
                    "ENSG_BCR",
                    "BCR-other",
                    "chr22",
                    23_180_000,
                    23_318_000,
                    "+",
                    "lncRNA",
                    "[]",
                    None,
                    None,
                    None,
                    0,
                    0,
                    0,
                    138_000,
                ),
                (
                    "ENST_CANONICAL",
                    "ENSG_BCR",
                    "BCR-canonical",
                    "chr22",
                    23_180_000,
                    23_318_000,
                    "+",
                    "protein_coding",
                    "[]",
                    None,
                    None,
                    "appris_principal_1",
                    0,
                    1,
                    2_000,
                    138_000,
                ),
                (
                    "ENST_PLUS",
                    "ENSG_BCR",
                    "BCR-plus",
                    "chr22",
                    23_180_000,
                    23_318_000,
                    "+",
                    "protein_coding",
                    "[]",
                    "MANE Plus Clinical",
                    "NM_PLUS",
                    None,
                    0,
                    1,
                    2_100,
                    138_000,
                ),
                (
                    "ENST_SELECT",
                    "ENSG_BCR",
                    "BCR-select",
                    "chr22",
                    23_180_000,
                    23_318_000,
                    "+",
                    "protein_coding",
                    "[]",
                    "MANE Select",
                    "NM_SELECT",
                    None,
                    0,
                    1,
                    1_500,
                    137_000,
                ),
                (
                    "ENST_P2RY8_X",
                    "ENSG_P2RY8_X",
                    "P2RY8-X",
                    "chrX",
                    1_462_571,
                    1_537_207,
                    "-",
                    "protein_coding",
                    "[]",
                    "MANE Select",
                    "NM_P2RY8",
                    None,
                    0,
                    1,
                    1_000,
                    74_636,
                ),
                (
                    "ENST_P2RY8_Y",
                    "ENSG_P2RY8_Y",
                    "P2RY8-Y",
                    "chrY",
                    1_462_571,
                    1_537_207,
                    "-",
                    "protein_coding",
                    "[]",
                    None,
                    None,
                    None,
                    0,
                    1,
                    1_000,
                    74_636,
                ),
            ],
        )
        connection.commit()


def _panel_bundle(path: Path) -> Path:
    derived = path / "derived"
    derived.mkdir(parents=True)
    selection = derived / "selection.bed"
    selection.write_text("chr22\t0\t30000000\tBCR\n", encoding="utf-8")
    manifest = {
        "schema_version": "1.0.0",
        "bundle_type": "panel",
        "bundle_id": "SYNTHETIC_GRCh38_v1",
        "version": "v1",
        "genome_build": "GRCh38",
        "assay_mode": "adaptive_sampling",
        "resources": [
            {
                "resource_id": "selection",
                "role": "selection_panel_buffered",
                "path": "derived/selection.bed",
                "sha256": sha256_file(selection),
                "coordinate_system": "zero_based_half_open",
            },
            {
                "resource_id": "roi",
                "role": "analysis_roi_unbuffered",
                "path": "derived/analysis_roi.bed",
                "generated": True,
                "coordinate_system": "zero_based_half_open",
                "derived_from": [
                    "selection",
                    "GRCh38_GENCODE50_MANE1.5_v1:annotation_cache",
                ],
            },
            {
                "resource_id": "transcripts",
                "role": "transcript_cache",
                "path": "derived/transcripts.tsv",
                "generated": True,
                "derived_from": [
                    "selection",
                    "GRCh38_GENCODE50_MANE1.5_v1:annotation_cache",
                ],
            },
        ],
        "selection_panel_resource_id": "selection",
        "analysis_roi_resource_id": "roi",
        "transcript_cache_resource_id": "transcripts",
        "unresolved_targets": [],
    }
    manifest_path = path / "bundle.yaml"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return manifest_path


class PanelDerivativeCompilerTests(unittest.TestCase):
    def test_roi_uses_exact_gene_body_and_excludes_unresolved_targets(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            cache = root / "annotations.sqlite"
            _annotation_cache(cache)
            selection = root / "selection.bed"
            selection.write_text(
                "chr22\t23170508\t23328037\tBCR\n"
                "chr5\t143396958\t143417420\tIGH_REVIEW_REQUIRED\n"
                "chr1\t0\t500\tDUP\n"
                "chr2\t0\t500\tMISSING\n",
                encoding="utf-8",
            )
            roi = root / "analysis_roi.bed"
            transcripts = root / "transcripts.tsv"
            summary = compile_panel_derivatives(selection, cache, roi, transcripts)
            self.assertEqual(summary.target_count, 4)
            self.assertEqual(summary.resolved_target_count, 1)
            self.assertEqual(
                summary.unresolved_targets,
                ("IGH_REVIEW_REQUIRED", "DUP", "MISSING"),
            )
            self.assertEqual(
                roi.read_text(encoding="utf-8"),
                "chr22\t23180000\t23318000\tBCR\tENSG_BCR\n",
            )

    def test_transcript_rank_is_mane_then_canonical_then_other(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            cache = root / "annotations.sqlite"
            _annotation_cache(cache)
            selection = root / "selection.bed"
            selection.write_text("chr22\t0\t30_000_000\tBCR\n".replace("_", ""))
            roi = root / "roi.bed"
            transcripts = root / "transcripts.tsv"
            summary = compile_panel_derivatives(selection, cache, roi, transcripts)
            rows = [line.split("\t") for line in transcripts.read_text().splitlines()]
            self.assertEqual(
                [row[4] for row in rows[1:]],
                ["ENST_SELECT", "ENST_PLUS", "ENST_CANONICAL", "ENST_OTHER"],
            )
            self.assertEqual([row[2] for row in rows[1:]], ["1", "2", "3", "5"])
            self.assertEqual([row[1] for row in rows[1:]], ["true", "false", "false", "false"])
            self.assertEqual(summary.transcript_count, 4)
            self.assertEqual(summary.preferred_transcript_count, 1)

    def test_par_gene_symbol_is_disambiguated_by_declared_target_chromosome(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            cache = root / "annotations.sqlite"
            _annotation_cache(cache)
            selection = root / "selection.bed"
            selection.write_text("chrX\t1452580\t1547186\tP2RY8\n", encoding="utf-8")
            roi = root / "roi.bed"
            transcripts = root / "transcripts.tsv"

            summary = compile_panel_derivatives(selection, cache, roi, transcripts)
            roi_text = roi.read_text(encoding="utf-8")

        self.assertEqual(summary.resolved_target_count, 1)
        self.assertEqual(summary.unresolved_targets, ())
        self.assertIn("ENSG_P2RY8_X", roi_text)
        self.assertNotIn("ENSG_P2RY8_Y", roi_text)

    def test_non_grch38_annotation_cache_fails_closed_without_outputs(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            cache = root / "annotations.sqlite"
            _annotation_cache(cache, genome_build="GRCh37")
            selection = root / "selection.bed"
            selection.write_text("chr22\t0\t30000000\tBCR\n")
            roi = root / "roi.bed"
            transcripts = root / "transcripts.tsv"
            with self.assertRaisesRegex(PanelCompilerError, "require a GRCh38"):
                compile_panel_derivatives(selection, cache, roi, transcripts)
            self.assertFalse(roi.exists())
            self.assertFalse(transcripts.exists())


class PanelMaterializationTests(unittest.TestCase):
    def test_derivatives_are_compiled_pinned_and_manifest_validated(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            bundle_root = root / "SYNTHETIC_GRCh38_v1"
            bundle_root.mkdir()
            manifest_path = _panel_bundle(bundle_root)
            annotation_cache = root / "annotations.sqlite"
            _annotation_cache(annotation_cache)
            summary = materialize_and_pin_panel_derivatives(bundle_root, annotation_cache)

            document = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
            bundle = PanelBundle.model_validate(document)
            roi = bundle.resource(bundle.analysis_roi_resource_id)
            transcripts = bundle.resource(bundle.transcript_cache_resource_id)
            roi_path = bundle_root / roi.path
            transcript_path = bundle_root / transcripts.path
            self.assertEqual(roi.sha256, sha256_file(roi_path))
            self.assertEqual(transcripts.sha256, sha256_file(transcript_path))
            self.assertEqual(roi.size_bytes, roi_path.stat().st_size)
            self.assertEqual(transcripts.size_bytes, transcript_path.stat().st_size)
            self.assertEqual(summary.analysis_roi_sha256, roi.sha256)
            self.assertEqual(summary.transcript_cache_sha256, transcripts.sha256)
            self.assertEqual(summary.manifest_sha256, sha256_file(manifest_path))
            self.assertEqual(summary.compilation.resolved_target_count, 1)

    def test_failed_compilation_leaves_manifest_and_outputs_untouched(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            bundle_root = root / "SYNTHETIC_GRCh38_v1"
            bundle_root.mkdir()
            manifest_path = _panel_bundle(bundle_root)
            original_manifest = manifest_path.read_bytes()
            annotation_cache = root / "annotations.sqlite"
            _annotation_cache(annotation_cache, genome_build="GRCh37")

            with self.assertRaisesRegex(PanelCompilerError, "require a GRCh38"):
                materialize_and_pin_panel_derivatives(bundle_root, annotation_cache)
            self.assertEqual(manifest_path.read_bytes(), original_manifest)
            self.assertFalse((bundle_root / "derived" / "analysis_roi.bed").exists())
            self.assertFalse((bundle_root / "derived" / "transcripts.tsv").exists())


if __name__ == "__main__":
    unittest.main()
