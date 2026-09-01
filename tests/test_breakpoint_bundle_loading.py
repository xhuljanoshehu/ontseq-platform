from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ontseq_platform.breakpoint_annotation import PathBackedContextIntervalIndex
from ontseq_platform.pipeline.runner import _bundle_breakpoint_context_resources


class BreakpointBundleLoadingTests(unittest.TestCase):
    def test_bundle_bed_index_is_lazy_and_handles_discontiguous_contig_ranges(self) -> None:
        with TemporaryDirectory() as raw:
            bed = Path(raw) / "context.bed"
            bed.write_text(
                "# synthetic technical context\n"
                "chr2\t20\t40\tlate\n"
                "chr1\t100\t200\tfirst\n"
                "track name=context\n"
                "2\t0\t30\tearly\n"
                "chr1\t150\t250\tsecond\n"
                "chr3\t5\t10\n",
                encoding="utf-8",
            )

            resources = _bundle_breakpoint_context_resources({"repeatmasker": bed})
            index = resources["repeatmasker"]
            self.assertIsInstance(index, PathBackedContextIntervalIndex)
            assert isinstance(index, PathBackedContextIntervalIndex)
            self.assertEqual(index.cached_contigs, ())
            self.assertEqual(index.contig_interval_counts, {"chr2": 2, "chr1": 2, "chr3": 1})

            self.assertEqual(
                [item.label for item in index.overlaps("chr2", 25)],
                ["early", "late"],
            )
            self.assertEqual(index.cached_contigs, ("chr2",))
            self.assertEqual(
                [item.label for item in index.overlaps("1", 175)],
                ["first", "second"],
            )
            self.assertEqual(index.cached_contigs, ("chr2", "chr1"))
            self.assertEqual(index.overlaps("chr3", 10), ())
            self.assertEqual(index.overlaps("chr3", 9)[0].label, "repeatmasker")
            self.assertEqual(index.cached_contigs, ("chr1", "chr3"))


if __name__ == "__main__":
    unittest.main()
