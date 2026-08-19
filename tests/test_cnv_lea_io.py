from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cnv.cytobands import Cytoband, CytobandTable
from ontseq_platform.cnv.lea_audit import audit_lea_ground_truth
from ontseq_platform.cnv.lea_io import load_lea_ace_call_set
from ontseq_platform.cnv.lea_truth_tables import LeaGroundTruthRow
from ontseq_platform.cnv.models import CnvDataBasis
from ontseq_platform.models import GenomeBuild, ReferenceContig, ReferenceLock


def _cn_text() -> str:
    lines = ["Chromosome,Copies,Ploidy,CNA"]
    for index in range(1, 23):
        lines.append(f"{index},2,2,0")
    lines.extend(["X,1,1,0", "Y,1,1,0"])
    return "\n".join(lines) + "\n"


def _lock() -> ReferenceLock:
    return ReferenceLock(
        reference_id="SYNTHETIC_GRCH37",
        genome_build=GenomeBuild.GRCH37,
        contigs=[ReferenceContig(name=str(index), length=1000) for index in range(1, 23)]
        + [ReferenceContig(name="X", length=1000), ReferenceContig(name="Y", length=1000)],
        source_fai_sha256="0" * 64,
    )


class LeaFileBoundaryTests(unittest.TestCase):
    def test_local_import_records_hashes_but_not_source_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cn = root / "CN.csv"
            dels_dups = root / "dels_dups.csv"
            bands = root / "hg19-cytoBand.txt"
            cn.write_text(_cn_text(), encoding="utf-8")
            dels_dups.write_text(
                "chromosome,name,event,frac_abr\nchr5,q13.2,del,0.75\n",
                encoding="utf-8",
            )
            bands.write_text(
                "chr5\t0\t400\tp15\tgneg\n"
                "chr5\t400\t600\tq13.1\tgpos25\n"
                "chr5\t600\t800\tq13.2\tgneg\n"
                "chr5\t800\t1000\tq13.3\tgpos50\n",
                encoding="utf-8",
            )
            band_sha = hashlib.sha256(bands.read_bytes()).hexdigest()
            result = load_lea_ace_call_set(
                cn_csv=cn,
                dels_dups_csv=dels_dups,
                cytoband_file=bands,
                expected_cytoband_sha256=band_sha,
                cytoband_resource_id="SYNTHETIC_HG19_BANDS",
                reference_lock=_lock(),
                call_set_id="LEA_SYNTHETIC_IO_001",
                sample_id="SYNTHETIC_AML_001",
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
            )
            self.assertIsNotNone(result.tool)
            assert result.tool is not None
            self.assertEqual(result.tool.parameters["cytoband_sha256"], band_sha)
            self.assertEqual(
                result.tool.parameters["cn_csv_sha256"],
                hashlib.sha256(cn.read_bytes()).hexdigest(),
            )
            serialized = result.model_dump_json()
            self.assertNotIn(str(root), serialized)

    def test_wrong_cytoband_checksum_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cn = root / "CN.csv"
            dels_dups = root / "dels_dups.csv"
            bands = root / "hg19-cytoBand.txt"
            cn.write_text(_cn_text(), encoding="utf-8")
            dels_dups.write_text("chromosome,name,event,frac_abr\n", encoding="utf-8")
            bands.write_text("chr5\t0\t1000\tq13\tgneg\n", encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                load_lea_ace_call_set(
                    cn_csv=cn,
                    dels_dups_csv=dels_dups,
                    cytoband_file=bands,
                    expected_cytoband_sha256="f" * 64,
                    cytoband_resource_id="SYNTHETIC_HG19_BANDS",
                    reference_lock=_lock(),
                    call_set_id="LEA_SYNTHETIC_IO_002",
                    sample_id="SYNTHETIC_AML_002",
                    data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
                )
            self.assertIn("SHA-256", str(context.exception))


class LeaTruthAuditTests(unittest.TestCase):
    def test_audit_is_aggregate_only_and_reconciles(self) -> None:
        table = CytobandTable(
            resource_id="SYNTHETIC_HG19_BANDS",
            genome_build=GenomeBuild.GRCH37,
            source="synthetic",
            bands=[Cytoband(contig="chr5", start=0, end=1000, name="q13", stain="gneg")],
        )
        summary = audit_lea_ground_truth(
            [
                LeaGroundTruthRow(sample_id="PRIVATE_A", karyotype="46,XX"),
                LeaGroundTruthRow(sample_id="PRIVATE_B", karyotype="46,XX,add(5)(q13)"),
            ],
            table,
        )
        self.assertEqual(summary.total_rows, 2)
        self.assertEqual(summary.fully_convertible_rows, 1)
        self.assertEqual(summary.incomplete_rows, 1)
        self.assertEqual(summary.unsupported_construct_count, 1)
        serialized = summary.model_dump_json()
        self.assertNotIn("PRIVATE_A", serialized)
        self.assertNotIn("PRIVATE_B", serialized)
        self.assertNotIn("add(5)(q13)", serialized)
        self.assertFalse(summary.contains_sample_identifiers)

    def test_audit_counts_invisible_unicode_without_normalizing_it(self) -> None:
        table = CytobandTable(
            resource_id="SYNTHETIC_HG19_BANDS",
            genome_build=GenomeBuild.GRCH37,
            source="synthetic",
            bands=[Cytoband(contig="chr5", start=0, end=1000, name="q13", stain="gneg")],
        )
        summary = audit_lea_ground_truth(
            [LeaGroundTruthRow(sample_id="PRIVATE_C", karyotype="\u200b46,XX")],
            table,
        )
        self.assertEqual(summary.formatting_issue_rows, 1)
        self.assertEqual(summary.formatting_issue_counts["zero_width_space_u200b"], 1)
        self.assertEqual(summary.incomplete_rows, 1)


if __name__ == "__main__":
    unittest.main()
