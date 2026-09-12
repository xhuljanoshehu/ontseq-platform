from __future__ import annotations

import gzip
import sqlite3
import unittest
from bisect import bisect_right
from contextlib import closing
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from ontseq_platform.breakpoint_annotation import (
    Breakpoint,
    BreakpointAnnotationError,
    ContextInterval,
    ContextIntervalIndex,
    PathBackedContextIntervalIndex,
    annotate_breakpoint_pair,
    annotate_events_from_cache,
)
from ontseq_platform.models import EventType, Evidence, GenomicEvent, Locus


def _cache(path: Path, *, build: str = "GRCh38") -> None:
    with closing(sqlite3.connect(path)) as connection:
        connection.executescript(
            """
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE genes (
              gene_id TEXT PRIMARY KEY, gene_name TEXT, chrom TEXT, start INTEGER,
              end INTEGER, strand TEXT, gene_type TEXT
            );
            CREATE TABLE transcripts (
              transcript_id TEXT PRIMARY KEY, gene_id TEXT, transcript_name TEXT,
              chrom TEXT, start INTEGER, end INTEGER, strand TEXT, transcript_type TEXT,
              tags_json TEXT, mane_status TEXT, mane_refseq_id TEXT, appris TEXT,
              is_canonical INTEGER, is_basic INTEGER, cds_length INTEGER,
              transcript_length INTEGER
            );
            CREATE TABLE exons (
              transcript_id TEXT, exon_id TEXT, exon_number INTEGER, chrom TEXT,
              start INTEGER, end INTEGER, strand TEXT
            );
            CREATE TABLE cds (
              transcript_id TEXT, exon_id TEXT, chrom TEXT, start INTEGER,
              end INTEGER, strand TEXT, phase INTEGER
            );
            CREATE TABLE cytobands (
              chrom TEXT, start INTEGER, end INTEGER, name TEXT, gie_stain TEXT
            );
            """
        )
        connection.execute("INSERT INTO metadata VALUES ('genome_build', ?)", (build,))
        connection.executemany(
            "INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?)",
            [
                ("G_BCR", "BCR", "chr22", 100, 400, "+", "protein_coding"),
                ("G_ABL1", "ABL1", "chr9", 1_000, 1_500, "-", "protein_coding"),
            ],
        )
        connection.executemany(
            "INSERT INTO transcripts VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                (
                    "T_BCR_SELECT",
                    "G_BCR",
                    "BCR-201",
                    "chr22",
                    100,
                    400,
                    "+",
                    "protein_coding",
                    "[]",
                    "MANE Select",
                    "NM_BCR",
                    None,
                    0,
                    1,
                    180,
                    300,
                ),
                (
                    "T_BCR_OTHER",
                    "G_BCR",
                    "BCR-202",
                    "chr22",
                    100,
                    400,
                    "+",
                    "protein_coding",
                    "[]",
                    None,
                    None,
                    None,
                    0,
                    1,
                    180,
                    300,
                ),
                (
                    "T_ABL1_SELECT",
                    "G_ABL1",
                    "ABL1-201",
                    "chr9",
                    1_000,
                    1_500,
                    "-",
                    "protein_coding",
                    "[]",
                    "MANE Select",
                    "NM_ABL1",
                    None,
                    0,
                    1,
                    200,
                    500,
                ),
            ],
        )
        connection.executemany(
            "INSERT INTO exons VALUES (?, ?, ?, ?, ?, ?, ?)",
            [
                ("T_BCR_SELECT", "BCR_E1", 1, "chr22", 100, 200, "+"),
                ("T_BCR_SELECT", "BCR_E2", 2, "chr22", 300, 400, "+"),
                ("T_BCR_OTHER", "BCR_E1", 1, "chr22", 100, 200, "+"),
                ("T_BCR_OTHER", "BCR_E2", 2, "chr22", 300, 400, "+"),
                ("T_ABL1_SELECT", "ABL_E1", 1, "chr9", 1_400, 1_500, "-"),
                ("T_ABL1_SELECT", "ABL_E2", 2, "chr9", 1_200, 1_300, "-"),
            ],
        )
        connection.executemany(
            "INSERT INTO cds VALUES (?, ?, ?, ?, ?, ?, ?)",
            [
                ("T_BCR_SELECT", "BCR_E1", "chr22", 120, 190, "+", 0),
                ("T_ABL1_SELECT", "ABL_E2", "chr9", 1_220, 1_280, "-", 2),
            ],
        )
        connection.executemany(
            "INSERT INTO cytobands VALUES (?, ?, ?, ?, ?)",
            [
                ("chr22", 0, 500, "q11.23", "gneg"),
                ("chr9", 900, 1_600, "q34.12", "gpos"),
            ],
        )
        connection.commit()


class BreakpointAnnotationTests(unittest.TestCase):
    def test_context_index_matches_linear_half_open_overlap_and_bisects_one_contig(self) -> None:
        intervals = [
            ContextInterval("chr2", start, start + 5, f"other-{start}")
            for start in range(0, 2_000, 5)
        ]
        intervals.extend(
            [
                ContextInterval("chr1", 50, 1_500, "nested-long"),
                ContextInterval("1", 1_000, 1_100, "first", 1.0),
                ContextInterval("chr1", 1_000, 1_050, "same-start"),
                ContextInterval("chr1", 1_099, 1_101, "edge"),
            ]
        )
        # Deliberately append expired rows out of order to exercise sorting and block pruning.
        intervals.extend(
            ContextInterval("1", start, start + 2, f"expired-{start}")
            for start in range(900, 0, -3)
        )
        index = ContextIntervalIndex(intervals)

        def payload(items: tuple[ContextInterval, ...]) -> list[tuple[int, int, str, float | None]]:
            return sorted((item.start, item.end, item.label, item.value) for item in items)

        for position in (49, 50, 999, 1_000, 1_049, 1_099, 1_100, 1_101, 1_499, 1_500):
            expected = tuple(
                interval
                for interval in intervals
                if interval.chromosome.removeprefix("chr") == "1"
                and interval.start <= position < interval.end
            )
            self.assertEqual(payload(index.overlaps("chr1", position)), payload(expected))

        with patch(
            "ontseq_platform.breakpoint_annotation.bisect_right",
            wraps=bisect_right,
        ) as mocked_bisect:
            index.overlaps("1", 1_000)

        mocked_bisect.assert_called_once()
        searched_starts = mocked_bisect.call_args.args[0]
        self.assertEqual(len(searched_starts), index.contig_interval_counts["chr1"])
        self.assertLess(len(searched_starts), len(intervals))

    def test_bundle_bed_index_uses_seekable_uncompressed_offsets_for_gzip(self) -> None:
        with TemporaryDirectory() as raw:
            bed = Path(raw) / "mappability.bed.gz"
            with gzip.open(bed, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write("chr1\t0\t10\tfirst\nchr2\t100\t200\tother\nchr1\t20\t30\tsecond\n")
            index = PathBackedContextIntervalIndex(bed, resource_type="mappability")

            self.assertEqual(index.contig_interval_counts, {"chr1": 2, "chr2": 1})
            self.assertEqual(index.overlaps("1", 9)[0].label, "first")
            self.assertEqual(index.overlaps("chr1", 25)[0].label, "second")
            self.assertEqual(index.cached_contigs, ("chr1",))

    def test_event_batch_indexes_a_plain_sequence_only_once(self) -> None:
        class SinglePassIntervals(list[ContextInterval]):
            iterations = 0

            def __iter__(self):  # type: ignore[no-untyped-def]
                self.iterations += 1
                if self.iterations > 1:
                    raise AssertionError("context intervals were re-scanned")
                return super().__iter__()

        with TemporaryDirectory() as raw:
            cache = Path(raw) / "annotations.sqlite"
            _cache(cache)
            intervals = SinglePassIntervals([ContextInterval("chr22", 100, 200, "context")])
            events = [
                GenomicEvent(
                    event_id=f"INS_{number}",
                    event_type=EventType.INSERTION,
                    primary=Locus(chromosome="chr22", start=150, end=151),
                )
                for number in range(2)
            ]

            annotated = annotate_events_from_cache(
                events,
                cache,
                context_resources={"repeatmasker": intervals},
            )

            self.assertEqual(intervals.iterations, 1)
            self.assertTrue(
                all(
                    event.breakpoint_annotations[0].contexts[0].label == "context"
                    for event in annotated
                )
            )

    def test_both_breakpoints_keep_transcript_and_context_evidence(self) -> None:
        with TemporaryDirectory() as raw:
            cache = Path(raw) / "annotations.sqlite"
            _cache(cache)
            pair = annotate_breakpoint_pair(
                cache,
                Breakpoint("22", 150),
                Breakpoint("chr9", 1_250),
                context_resources={
                    "repeatmasker": [ContextInterval("chr22", 140, 160, "LINE/L1")],
                    "blacklist": [ContextInterval("chr9", 1_200, 1_300, "ENCODE")],
                    "mappability": [ContextInterval("22", 100, 200, "unique", 1.0)],
                },
            )
            self.assertEqual(pair.primary.cytoband, "q11.23")
            self.assertEqual(pair.primary.genes, ("BCR",))
            self.assertEqual(pair.primary.transcripts[0].transcript_id, "T_BCR_SELECT")
            self.assertTrue(pair.primary.transcripts[0].preferred)
            self.assertEqual(pair.primary.transcripts[0].region, "exon")
            self.assertEqual(pair.primary.transcripts[0].exon_number, 1)
            self.assertEqual(pair.primary.transcripts[0].cds_phase, 0)
            self.assertEqual(
                {hit.resource_type for hit in pair.primary.contexts},
                {"repeatmasker", "mappability"},
            )
            self.assertIsNotNone(pair.secondary)
            assert pair.secondary is not None
            self.assertEqual(pair.secondary.cytoband, "q34.12")
            self.assertEqual(pair.secondary.genes, ("ABL1",))
            self.assertEqual(pair.secondary.transcripts[0].exon_number, 2)
            self.assertEqual(pair.secondary.transcripts[0].cds_phase, 2)
            self.assertEqual(pair.secondary.contexts[0].resource_type, "blacklist")

    def test_minus_strand_intron_number_uses_transcript_order(self) -> None:
        with TemporaryDirectory() as raw:
            cache = Path(raw) / "annotations.sqlite"
            _cache(cache)
            pair = annotate_breakpoint_pair(cache, Breakpoint("chr9", 1_350))
            hit = pair.primary.transcripts[0]
            self.assertEqual(hit.strand, "-")
            self.assertEqual(hit.region, "intron")
            self.assertEqual(hit.intron_number, 1)
            self.assertIsNone(hit.cds_phase)

    def test_grch37_cache_is_refused(self) -> None:
        with TemporaryDirectory() as raw:
            cache = Path(raw) / "annotations.sqlite"
            _cache(cache, build="GRCh37")
            with self.assertRaisesRegex(BreakpointAnnotationError, "requires GRCh38"):
                annotate_breakpoint_pair(cache, Breakpoint("chr22", 150))

    def test_pipeline_event_contract_retains_both_breakpoints_and_unknown_frame(self) -> None:
        with TemporaryDirectory() as raw:
            cache = Path(raw) / "annotations.sqlite"
            _cache(cache)
            event = GenomicEvent(
                event_id="BND_1",
                event_type=EventType.TRANSLOCATION,
                primary=Locus(chromosome="chr22", start=150, end=151),
                secondary=Locus(chromosome="chr9", start=1_250, end=1_251),
                evidence=[
                    Evidence(
                        caller="synthetic",
                        caller_version="1",
                        supporting_read_strands="+-",
                    )
                ],
            )

            annotated = annotate_events_from_cache([event], cache)[0]

            self.assertEqual(
                [item.label for item in annotated.breakpoint_annotations],
                ["primary", "secondary"],
            )
            self.assertEqual(annotated.genes, ["ABL1", "BCR"])
            self.assertIsNotNone(annotated.fusion_evidence)
            assert annotated.fusion_evidence is not None
            self.assertEqual(annotated.fusion_evidence.gene_a.gene, "BCR")
            self.assertEqual(annotated.fusion_evidence.gene_b.gene, "ABL1")
            self.assertEqual(annotated.fusion_evidence.orientation, "+-")
            self.assertEqual(annotated.fusion_evidence.frame_status, "unknown")


if __name__ == "__main__":
    unittest.main()
