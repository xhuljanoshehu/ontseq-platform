from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.annotation_cache import compile_annotation_cache
from ontseq_platform.breakpoint_annotation import ContextInterval, annotate_events_from_cache
from ontseq_platform.models import (
    EventType,
    Evidence,
    GenomeBuild,
    GenomicEvent,
    IntervalResourceLock,
    Locus,
)
from ontseq_platform.reference import sha256_file
from ontseq_platform.sv_annotation import annotate_sv_events
from ontseq_platform.sv_evidence import classify_sv_event


def _event() -> GenomicEvent:
    return GenomicEvent(
        event_id="SYNTHETIC-MAPPABILITY-TEST",
        event_type=EventType.DELETION,
        primary=Locus(chromosome="chr7", start=120, end=180),
        evidence=[Evidence(caller="synthetic", caller_version="1", support_reads=20)],
    )


def _cache(root: Path, build: str) -> Path:
    gtf = root / "synthetic.gtf"
    gtf.write_text(
        'chr7\tSYNTHETIC\tgene\t101\t300\t.\t+\t.\tgene_id "SYN_G"; '
        'gene_name "SYNTHETIC"; gene_type "protein_coding";\n'
        'chr7\tSYNTHETIC\ttranscript\t101\t300\t.\t+\t.\tgene_id "SYN_G"; '
        'transcript_id "SYN_T.1"; transcript_type "protein_coding"; tag "basic";\n'
        'chr7\tSYNTHETIC\texon\t101\t300\t.\t+\t.\tgene_id "SYN_G"; '
        'transcript_id "SYN_T.1"; exon_number "1";\n',
        encoding="utf-8",
    )
    bands = root / "synthetic.bands"
    bands.write_text("chr7\t0\t500\tp11\tgneg\n", encoding="utf-8")
    mane = root / "synthetic.mane.gff3" if build == "GRCh38" else None
    if mane is not None:
        mane.write_text(
            "##gff-version 3\n"
            "chr7\tSYNTHETIC\tmRNA\t101\t300\t.\t+\t.\tID=transcript:SYN_T.1;"
            "gene_id=SYN_G;tag=MANE%20Select;Dbxref=RefSeq:NM_SYN.1\n",
            encoding="utf-8",
        )
    output = root / "annotation.sqlite"
    compile_annotation_cache(
        gtf,
        mane,
        bands,
        output,
        metadata={
            "genome_build": build,
            "gencode_release": "GENCODE 19" if build == "GRCh37" else "GENCODE 50",
        },
    )
    return output


class MappabilityContextSemanticsTests(unittest.TestCase):
    def test_positive_overlap_keeps_context_without_penalty_for_both_builds(self) -> None:
        for build in ("GRCh37", "GRCh38"):
            with self.subTest(build=build), tempfile.TemporaryDirectory() as raw:
                cache = _cache(Path(raw), build)
                annotated = annotate_events_from_cache(
                    [_event()],
                    annotation_cache=cache,
                    expected_build=build,
                    context_resources={
                        "mappability": [ContextInterval("chr7", 100, 200, "Umap S100")]
                    },
                )[0]
                self.assertEqual(annotated.technical_flags, [])
                self.assertTrue(annotated.breakpoint_annotations)
                self.assertTrue(
                    all(annotation.contexts for annotation in annotated.breakpoint_annotations)
                )
                self.assertEqual(
                    annotated.breakpoint_annotations[0].contexts[0].resource_type, "mappability"
                )
                scored = classify_sv_event(annotated)
                self.assertEqual(scored.confidence, classify_sv_event(_event()).confidence)
                self.assertFalse(any("artifact_context_penalty" in note for note in scored.notes))
                self.assertFalse(scored.reportable)

    def test_absent_umap_row_does_not_create_low_mappability_evidence(self) -> None:
        for build in ("GRCh37", "GRCh38"):
            with self.subTest(build=build), tempfile.TemporaryDirectory() as raw:
                annotated = annotate_events_from_cache(
                    [_event()],
                    annotation_cache=_cache(Path(raw), build),
                    expected_build=build,
                    context_resources={
                        "mappability": [ContextInterval("chr7", 300, 400, "Umap S100")]
                    },
                )[0]
                self.assertEqual(annotated.technical_flags, [])
                self.assertTrue(
                    all(not annotation.contexts for annotation in annotated.breakpoint_annotations)
                )

    def test_legacy_interval_path_keeps_mappability_note_without_flag_for_both_builds(self) -> None:
        for build in (GenomeBuild.GRCH37, GenomeBuild.GRCH38):
            with self.subTest(build=build), tempfile.TemporaryDirectory() as raw:
                path = Path(raw) / "synthetic.umap.bed"
                path.write_text("chr7\t100\t200\tUmap S100\n", encoding="utf-8")
                lock = IntervalResourceLock(
                    resource_id="synthetic-umap",
                    resource_type="mappability",
                    source="synthetic-only",
                    release="synthetic-only",
                    genome_build=build,
                    sha256=sha256_file(path),
                )
                annotated = annotate_sv_events(
                    [_event()], genome_build=build, context_resources=[(path, lock)]
                )[0]
                self.assertEqual(annotated.technical_flags, [])
                self.assertTrue(
                    any("mappability reference context" in note for note in annotated.notes)
                )

    def test_old_mappability_flags_are_not_penalized_but_real_artifacts_still_are(self) -> None:
        event = _event().model_copy(
            update={"technical_flags": ["primary:mappability", "secondary:mappability"]}
        )
        self.assertFalse(
            any("artifact_context_penalty" in note for note in classify_sv_event(event).notes)
        )
        event = event.model_copy(
            update={"technical_flags": [*event.technical_flags, "primary:blacklist"]}
        )
        scored = classify_sv_event(event)
        self.assertTrue(any("artifact_context_penalty=1" in note for note in scored.notes))
        self.assertFalse(scored.reportable)


if __name__ == "__main__":
    unittest.main()
