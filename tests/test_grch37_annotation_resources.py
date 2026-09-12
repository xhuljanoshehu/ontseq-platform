from __future__ import annotations

import gzip
import hashlib
import io
import json
import sqlite3
import tempfile
import unittest
from contextlib import closing
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import yaml

from ontseq_platform.annotation_cache import (
    AnnotationCache,
    compile_annotation_cache,
    validate_annotation_cache,
)
from ontseq_platform.breakpoint_annotation import (
    Breakpoint,
    BreakpointAnnotationError,
    annotate_breakpoint_pair,
    annotate_events_from_cache,
)
from ontseq_platform.cnv.extension import (
    _annotate_cnv_cytobands,
    _load_cytobands,
    _verified_iscn_resource_provenance,
)
from ontseq_platform.cnv.qdnaseq import CnvFit, QDNAseqCallReport
from ontseq_platform.models import (
    CoordinateSystem,
    EventType,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    ReferenceBundle,
    ReferenceContig,
    ReferenceDictionaryContract,
    ResolvedResourceContext,
    ResourceFile,
)
from ontseq_platform.reference import sha256_file
from ontseq_platform.reference_catalog import (
    ReferenceBundleInstaller,
    ReferenceCatalog,
    ResourceValidationState,
    _normalize_grch37_blacklist,
    _normalize_ucsc_context_table,
    validate_reference_bundle_directory,
)

# Synthetic, version-19-shaped GTF records; these are not human reference coordinates.
GTF = (
    'chr7\tSYNTHETIC\tgene\t101\t400\t.\t+\t.\tgene_id "SYN_G1"; '
    'gene_name "SYNTHETIC_GENE"; gene_type "protein_coding";\n'
    'chr7\tSYNTHETIC\ttranscript\t101\t400\t.\t+\t.\tgene_id "SYN_G1"; '
    'transcript_id "SYN_T1"; transcript_type "protein_coding"; tag "basic";\n'
    'chr7\tSYNTHETIC\texon\t101\t200\t.\t+\t.\tgene_id "SYN_G1"; '
    'transcript_id "SYN_T1"; exon_number "1";\n'
    'chr7\tSYNTHETIC\tCDS\t111\t200\t.\t+\t0\tgene_id "SYN_G1"; '
    'transcript_id "SYN_T1";\n'
    'chr7\tSYNTHETIC\texon\t301\t400\t.\t+\t.\tgene_id "SYN_G1"; '
    'transcript_id "SYN_T1"; exon_number "2";\n'
)
CYTOBANDS = "chr7\t0\t500\tp11\tgneg\nchr7\t500\t1000\tq11\tgpos50\n"


def _cache(root: Path, *, name: str = "annotation.sqlite") -> Path:
    gtf = root / "synthetic.gencode19.gtf"
    bands = root / "synthetic.hg19.cytobands.tsv"
    gtf.write_text(GTF, encoding="utf-8")
    bands.write_text(CYTOBANDS, encoding="utf-8")
    output = root / name
    compile_annotation_cache(
        gtf,
        None,
        bands,
        output,
        metadata={"genome_build": "GRCh37", "gencode_release": "GENCODE 19"},
    )
    return output


def _cnv_report(events: list[GenomicEvent]) -> QDNAseqCallReport:
    fit = CnvFit(
        bin_size_kbp=500,
        cellularity=0.5,
        ploidy=2.0,
        fit_error=0.1,
        candidate_count=1,
        segment_count=len(events),
        segment_file="synthetic.segments.tsv",
        chromosome_file="synthetic.chromosomes.tsv",
        fit_plot="synthetic.fit.png",
        copy_number_plot="synthetic.copy-number.png",
        rds_file="synthetic.rds",
    )
    return QDNAseqCallReport(
        sample_id="SYNTHETIC",
        genome_build=GenomeBuild.GRCH37,
        status=ModuleRunStatus.COMPLETED if events else ModuleRunStatus.NO_CALL,
        primary_fit=fit,
        fits=[fit],
        chromosome_consensus=[],
        events=events,
        tools=[],
        output_files=[],
    )


def _recipe() -> tuple[ReferenceBundle, dict[str, bytes]]:
    contents = {
        "genome": b">chr7\n" + b"ACGT" * 250 + b"\n",
        "fai": b"chr7\t1000\t6\t1000\t1001\n",
        "gencode": GTF.encode(),
        "cytobands": CYTOBANDS.encode(),
    }
    roles = {
        "genome": ("genome_fasta", "genome.fa", None),
        "fai": ("fasta_index", "genome.fa.fai", None),
        "gencode": ("gencode_gtf", "sources/gencode.gtf", CoordinateSystem.ONE_BASED_INCLUSIVE),
        "cytobands": (
            "cytobands",
            "sources/cytobands.tsv",
            CoordinateSystem.ZERO_BASED_HALF_OPEN,
        ),
    }
    resources = [
        ResourceFile(
            resource_id=key,
            role=role,
            path=path,
            sha256=hashlib.sha256(contents[key]).hexdigest(),
            size_bytes=len(contents[key]),
            source_url=f"https://fixtures.invalid/{key}",
            release="GENCODE 19" if key == "gencode" else "synthetic-only",
            coordinate_system=coordinates,
        )
        for key, (role, path, coordinates) in roles.items()
    ]
    resources.extend(
        [
            ResourceFile(
                resource_id="reference_lock",
                role="reference_lock",
                path="reference.lock.json",
                generated=True,
                derived_from=["fai"],
            ),
            ResourceFile(
                resource_id="annotation_cache",
                role="annotation_cache",
                path="annotation.sqlite",
                generated=True,
                derived_from=["gencode", "cytobands"],
                coordinate_system=CoordinateSystem.ZERO_BASED_HALF_OPEN,
            ),
        ]
    )
    return (
        ReferenceBundle(
            bundle_id="GRCh37_SYNTHETIC_v1",
            version="synthetic-only",
            genome_build=GenomeBuild.GRCH37,
            resources=resources,
            reference_lock_resource_id="reference_lock",
            fasta_resource_id="genome",
            fai_resource_id="fai",
            annotation_cache_resource_id="annotation_cache",
        ),
        contents,
    )


class NativeGrch37AnnotationTests(unittest.TestCase):
    def test_native_cache_is_deterministic_explicit_and_has_no_mane(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            first = _cache(root)
            second = _cache(root, name="second.sqlite")
            self.assertEqual(first.read_bytes(), second.read_bytes())
            summary = validate_annotation_cache(first, expected_build="GRCh37")
            self.assertEqual(summary.mane_matched_transcripts, 0)
            self.assertEqual(summary.genes, 1)
            with closing(sqlite3.connect(first)) as connection:
                metadata = dict(connection.execute("SELECT key, value FROM metadata"))
            self.assertEqual(metadata["assembly_name"], "GRCh37.p13")
            self.assertEqual(metadata["mane_status"], "not_available_native_grch37")
            self.assertNotIn("mane_sha256", metadata)
            transcripts = AnnotationCache(first, expected_build="GRCh37").ranked_transcripts(
                "SYNTHETIC_GENE"
            )
            self.assertEqual(transcripts[0].transcript_id, "SYN_T1")
            self.assertIsNone(transcripts[0].mane_status)
            with self.assertRaisesRegex(ValueError, "requires GRCh38"):
                validate_annotation_cache(first)

    def test_native_compiler_refuses_mane_wrong_release_and_wrong_assembly(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            _cache(root)
            gtf = root / "synthetic.gencode19.gtf"
            bands = root / "synthetic.hg19.cytobands.tsv"
            for mane, metadata in (
                (gtf, {"genome_build": "GRCh37", "gencode_release": "GENCODE 19"}),
                (None, {"genome_build": "GRCh38"}),
                (None, {"genome_build": "GRCh37", "gencode_release": "GENCODE 50"}),
                (
                    None,
                    {
                        "genome_build": "GRCh37",
                        "gencode_release": "GENCODE 19",
                        "assembly_name": "GRCh38.p14",
                    },
                ),
            ):
                with self.subTest(metadata=metadata, mane=mane):
                    output = root / "must-not-publish.sqlite"
                    with self.assertRaises(ValueError):
                        compile_annotation_cache(gtf, mane, bands, output, metadata=metadata)
                    self.assertFalse(output.exists())

    def test_breakpoints_and_cnv_require_the_consuming_build(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            cache = _cache(Path(raw))
            pair = annotate_breakpoint_pair(cache, Breakpoint("chr7", 150), expected_build="GRCh37")
            self.assertEqual(pair.primary.genes, ("SYNTHETIC_GENE",))
            self.assertEqual(pair.primary.cytoband, "p11")
            self.assertEqual(len(_load_cytobands(cache, expected_build="GRCh37")), 2)
            with self.assertRaises(BreakpointAnnotationError):
                annotate_breakpoint_pair(cache, Breakpoint("chr7", 150))
            with self.assertRaises(ValueError):
                _load_cytobands(cache, expected_build="GRCh38")
            event = GenomicEvent(
                event_id="SYNTHETIC_EVENT",
                event_type=EventType.DELETION,
                primary=Locus(chromosome="chr7", start=100, end=400),
            )
            annotated = annotate_events_from_cache([event], cache, expected_build="GRCh37")
            self.assertIn("SYNTHETIC_GENE", annotated[0].genes)
            self.assertFalse(annotated[0].reportable)

    def test_cnv_artifact_reports_the_actual_build(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            cache = _cache(Path(raw))
            event = GenomicEvent(
                event_id="SYNTHETIC_CNV",
                event_type=EventType.DELETION,
                primary=Locus(chromosome="chr7", start=0, end=500),
            )
            report = _cnv_report([event])
            envelope = Mock()
            ctx = SimpleNamespace(
                config=SimpleNamespace(
                    annotation_cache=cache,
                    manifest=SimpleNamespace(
                        assay=SimpleNamespace(genome_build=GenomeBuild.GRCH37)
                    ),
                    reference_lock=SimpleNamespace(
                        contigs=[ReferenceContig(name="chr7", length=1000)]
                    ),
                ),
                envelope=envelope,
                path=lambda value: value.format(sample="SYNTHETIC"),
            )
            settings = SimpleNamespace(
                policy=SimpleNamespace(
                    cytoband_affected_fraction=0.66,
                    schema_version="0.1.0",
                    profile_id="synthetic-only",
                )
            )
            with patch("ontseq_platform.cnv.extension._settings", return_value=settings):
                _annotate_cnv_cytobands(ctx, report)
            payload = json.loads(envelope.atomic_write_text.call_args.args[1])
            self.assertEqual(payload["genome_build"], "GRCh37")

    def test_multi_group_cnv_projection_does_not_assign_an_iscn_band_range(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            cache = _cache(Path(raw))
            event = GenomicEvent(
                event_id="SYNTHETIC_MULTI_ARM_CNV",
                event_type=EventType.DELETION,
                primary=Locus(
                    chromosome="chr7",
                    start=0,
                    end=900,
                    cytoband_start="p11",
                    cytoband_end="q11",
                ),
                copy_number=1.0,
            )
            report = _cnv_report([event])
            envelope = Mock()
            ctx = SimpleNamespace(
                config=SimpleNamespace(
                    annotation_cache=cache,
                    manifest=SimpleNamespace(
                        assay=SimpleNamespace(genome_build=GenomeBuild.GRCH37)
                    ),
                    reference_lock=SimpleNamespace(
                        contigs=[ReferenceContig(name="chr7", length=1000)]
                    ),
                ),
                envelope=envelope,
                path=lambda value: value.format(sample="SYNTHETIC"),
            )
            settings = SimpleNamespace(
                policy=SimpleNamespace(
                    cytoband_affected_fraction=0.66,
                    schema_version="0.1.0",
                    profile_id="synthetic-only",
                )
            )

            with patch("ontseq_platform.cnv.extension._settings", return_value=settings):
                annotated, _artifact = _annotate_cnv_cytobands(ctx, report)

            projected = annotated.events[0]
            self.assertIsNone(projected.primary.cytoband_start)
            self.assertIsNone(projected.primary.cytoband_end)
            self.assertTrue(any("2 group(s)" in note for note in projected.notes))
            payload = json.loads(envelope.atomic_write_text.call_args.args[1])
            self.assertEqual(len(payload["affected_groups"]), 2)

    def test_multi_source_cytoband_group_does_not_assign_a_band_to_either_cnv(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            cache = _cache(Path(raw))
            events = [
                GenomicEvent(
                    event_id=event_id,
                    event_type=EventType.DELETION,
                    primary=Locus(
                        chromosome="chr7",
                        start=0,
                        end=400,
                        cytoband_start="p11",
                        cytoband_end="p11",
                    ),
                    copy_number=1.0,
                )
                for event_id in ("SYNTHETIC_CNV_A", "SYNTHETIC_CNV_B")
            ]
            report = _cnv_report(events)
            envelope = Mock()
            ctx = SimpleNamespace(
                config=SimpleNamespace(
                    annotation_cache=cache,
                    manifest=SimpleNamespace(
                        assay=SimpleNamespace(genome_build=GenomeBuild.GRCH37)
                    ),
                    reference_lock=SimpleNamespace(
                        contigs=[ReferenceContig(name="chr7", length=1000)]
                    ),
                ),
                envelope=envelope,
                path=lambda value: value.format(sample="SYNTHETIC"),
            )
            settings = SimpleNamespace(
                policy=SimpleNamespace(
                    cytoband_affected_fraction=0.66,
                    schema_version="0.1.0",
                    profile_id="synthetic-only",
                )
            )

            with patch("ontseq_platform.cnv.extension._settings", return_value=settings):
                annotated, _artifact = _annotate_cnv_cytobands(ctx, report)

            for projected in annotated.events:
                self.assertIsNone(projected.primary.cytoband_start)
                self.assertIsNone(projected.primary.cytoband_end)
                self.assertTrue(
                    any("1 group(s) spanning 2 source event(s)" in note for note in projected.notes)
                )
            payload = json.loads(envelope.atomic_write_text.call_args.args[1])
            self.assertEqual(len(payload["affected_groups"]), 1)
            self.assertEqual(
                payload["affected_groups"][0]["source_event_ids"],
                ["SYNTHETIC_CNV_A", "SYNTHETIC_CNV_B"],
            )

    def test_iscn_resource_verification_rejects_changed_source_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            gtf = root / "synthetic.gencode19.gtf"
            cytobands = root / "synthetic.hg19.cytobands.tsv"
            reference_lock = root / "reference.lock.json"
            cache = root / "annotation.sqlite"
            gtf.write_text(GTF, encoding="utf-8")
            cytobands.write_text(CYTOBANDS, encoding="utf-8")
            reference_lock.write_text("synthetic-reference-lock", encoding="utf-8")
            cytoband_sha256 = sha256_file(cytobands)
            compile_annotation_cache(
                gtf,
                None,
                cytobands,
                cache,
                metadata={
                    "genome_build": "GRCh37",
                    "gencode_release": "GENCODE 19",
                    "bundle_id": "GRCh37_SYNTHETIC_v1",
                    "bundle_version": "synthetic-v1",
                    "cytoband_release": "synthetic-cytobands-v1",
                    "cytoband_sha256": cytoband_sha256,
                },
            )
            resource_context = ResolvedResourceContext(
                profile_id="SYNTHETIC_GRCH37",
                profile_version="synthetic-v1",
                genome_build=GenomeBuild.GRCH37,
                reference_dictionary_contract=ReferenceDictionaryContract.EXACT_FULL,
                reference_bundle_id="GRCh37_SYNTHETIC_v1",
                reference_bundle_version="synthetic-v1",
                knowledge_bundle_id="SYNTHETIC_KNOWLEDGE",
                knowledge_bundle_version="synthetic-v1",
                resource_root=str(root),
                resource_paths={
                    "reference.reference_lock": str(reference_lock),
                    "reference.annotation_cache": str(cache),
                    "reference.cytobands": str(cytobands),
                },
                resource_checksums={
                    "reference.reference_lock": sha256_file(reference_lock),
                    "reference.annotation_cache": sha256_file(cache),
                    "reference.cytobands": cytoband_sha256,
                },
                resource_releases={
                    "reference.cytobands": "synthetic-cytobands-v1",
                },
            )
            ctx = SimpleNamespace(
                config=SimpleNamespace(
                    annotation_cache=cache,
                    resource_context=resource_context,
                ),
                manifest=SimpleNamespace(assay=SimpleNamespace(genome_build=GenomeBuild.GRCH37)),
            )

            provenance, blockers = _verified_iscn_resource_provenance(ctx)
            self.assertIsNotNone(provenance)
            self.assertEqual(blockers, [])

            cytobands.write_text(CYTOBANDS + "chr7\t1000\t1100\tq12\tgneg\n", encoding="utf-8")
            provenance, blockers = _verified_iscn_resource_provenance(ctx)
            self.assertIsNone(provenance)
            self.assertEqual(blockers[0].reason_code, "RESOURCE_CACHE_PROVENANCE_MISMATCH")
            self.assertIn("cytoband_sha256", blockers[0].detail)

    def test_relabelled_or_mane_claiming_cache_is_rejected(self) -> None:
        for change in (
            "UPDATE metadata SET value='GRCh38.p14' WHERE key='assembly_name'",
            "UPDATE transcripts SET mane_status='MANE Select'",
        ):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as raw:
                cache = _cache(Path(raw))
                with closing(sqlite3.connect(cache)) as connection:
                    connection.execute(change)
                    connection.commit()
                with self.assertRaises(ValueError):
                    validate_annotation_cache(cache, expected_build="GRCh37")
                if "assembly_name" in change:
                    with self.assertRaises(BreakpointAnnotationError):
                        annotate_breakpoint_pair(
                            cache, Breakpoint("chr7", 150), expected_build="GRCh37"
                        )


class NativeGrch37ReferenceTests(unittest.TestCase):
    def test_hg19_context_excludes_mitochondrial_and_alternate_sequences_only_for_native37(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            for role in ("repeatmasker", "simple_repeats", "segmental_duplication"):
                with self.subTest(role=role):
                    source = root / f"{role}.gz"
                    rows = []
                    for chrom in ("chr1", "chr22", "chrX", "chrY", "chrM", "chr6_apd_hap1"):
                        fields = (
                            ["0"] * 5 + [chrom, "10", "20", ".", ".", "repeat", "class", "family"]
                            if role == "repeatmasker"
                            else ["0", chrom, "10", "20", "synthetic_context"]
                        )
                        rows.append("\t".join(fields))
                    source.write_bytes(gzip.compress(("\n".join(rows) + "\n").encode(), mtime=0))
                    native = root / f"{role}.native.bed"
                    original = root / f"{role}.original.bed"
                    _normalize_ucsc_context_table(
                        source, native, role=role, genome_build=GenomeBuild.GRCH37
                    )
                    _normalize_ucsc_context_table(source, original, role=role)
                    lines = native.read_text().splitlines()
                    self.assertIn("excludes hg19 chrM/alternate contigs", lines[0])
                    self.assertEqual(
                        [line.split("\t")[0] for line in lines[1:]],
                        ["chr1", "chr22", "chrX", "chrY"],
                    )
                    self.assertEqual(len(original.read_text().splitlines()), 6)

    def test_native_blacklist_scope_and_bed3_validation(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            source, output = root / "blacklist.bed.gz", root / "blacklist.bed"
            source.write_bytes(
                gzip.compress(
                    b"track name=synthetic\nchr1\t10\t20\nchrY\t30\t40\tnamed\n"
                    b"chrM\t1\t10\nchr6_apd_hap1\t1\t10\n",
                    mtime=0,
                )
            )
            _normalize_grch37_blacklist(source, output)
            self.assertEqual(
                output.read_text().splitlines()[1:],
                [
                    "chr1\t10\t20\ttechnical_blacklist",
                    "chrY\t30\t40\tnamed",
                ],
            )
            source.write_bytes(gzip.compress(b"chr1\t20\t10\n", mtime=0))
            with self.assertRaisesRegex(ValueError, "invalid interval"):
                _normalize_grch37_blacklist(source, output)

    def test_pinned_native_recipe_installs_and_cache_build_is_validated(self) -> None:
        recipe, contents = _recipe()
        with (
            tempfile.TemporaryDirectory() as raw,
            patch("ontseq_platform.reference.validate_canonical_reference"),
            patch("ontseq_platform.reference_catalog.validate_canonical_reference"),
        ):
            root = Path(raw)
            manifest = root / "recipe.yaml"
            manifest.write_text(yaml.safe_dump(recipe.model_dump(mode="json")), encoding="utf-8")
            catalog = ReferenceCatalog.from_manifests([manifest])
            self.assertEqual(catalog.get(recipe.bundle_id).genome_build, GenomeBuild.GRCH37)
            installer = ReferenceBundleInstaller(
                root / "resources", opener=lambda url: io.BytesIO(contents[url.rsplit("/", 1)[-1]])
            )
            installed = installer.install(recipe)
            self.assertTrue(validate_reference_bundle_directory(installed.path).valid)
            cache = installed.path / "annotation.sqlite"
            validate_annotation_cache(cache, expected_build="GRCh37")

            # Re-pinning bytes does not legitimize a cache for the wrong assembly.
            with closing(sqlite3.connect(cache)) as connection:
                connection.execute("UPDATE metadata SET value='GRCh38' WHERE key='genome_build'")
                connection.commit()
            active_manifest = installed.path / "bundle.yaml"
            payload = yaml.safe_load(active_manifest.read_text(encoding="utf-8"))
            for resource in payload["resources"]:
                if resource["role"] == "annotation_cache":
                    resource["sha256"] = sha256_file(cache)
                    resource["size_bytes"] = cache.stat().st_size
            active_manifest.write_text(yaml.safe_dump(payload), encoding="utf-8")
            report = validate_reference_bundle_directory(installed.path)
            self.assertFalse(report.valid)
            cache_report = next(
                item for item in report.resources if item.role == "annotation_cache"
            )
            self.assertEqual(cache_report.state, ResourceValidationState.INVALID_DERIVED_ARTIFACT)

    def test_tiny_synthetic_recipe_cannot_bypass_production_reference_gate(self) -> None:
        recipe, contents = _recipe()
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            installer = ReferenceBundleInstaller(
                root, opener=lambda url: io.BytesIO(contents[url.rsplit("/", 1)[-1]])
            )
            with self.assertRaises(ValueError):
                installer.install(recipe)
            self.assertFalse((root / "references" / recipe.bundle_id / "bundle.yaml").exists())


if __name__ == "__main__":
    unittest.main()
