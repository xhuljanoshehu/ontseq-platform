#!/usr/bin/env python3
"""Verify that an ONTSeq wheel carries the immutable supported-build authority assets."""

from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZipFile

REQUIRED_SUFFIXES = (
    "ontseq_platform/iscn.py",
    "ontseq_platform/iscn_syntax.py",
    "ontseq_platform/methylation.py",
    "ontseq_platform/methylation_mixture.py",
    "ontseq_platform/modbam.py",
    "ontseq_platform/methylation_holdout_runner.py",
    "ontseq_platform/methylation_validation_runner.py",
    "share/ontseq/configs/methylation/modkit.technical.yaml",
    "share/ontseq/configs/methylation/paired_source_nanopolish.technical.yaml",
    "share/ontseq/configs/reference_bundles/GRCh38_GENCODE50_MANE1.5_v1/bundle.recipe.yaml",
    "share/ontseq/configs/reference_bundles/GRCh37_GENCODE19_v1/bundle.recipe.yaml",
    "share/ontseq/configs/reference_bundles/GRCh37_GENCODE19_HG19_v2/bundle.recipe.yaml",
    "share/ontseq/configs/knowledge_bundles/HEMATOLOGY_v1/bundle.yaml",
    "share/ontseq/configs/knowledge_bundles/HEMATOLOGY_v2/bundle.yaml",
    "share/ontseq/configs/knowledge_bundles/HEMATOLOGY_v3/bundle.yaml",
    "share/ontseq/configs/knowledge_bundles/HEMATOLOGY_GRCh37_v1/bundle.yaml",
    "share/ontseq/configs/panels/AML_AS_111_GRCh38_v1/bundle.yaml",
    "share/ontseq/configs/panels/AML_AS_111_GRCh38_v1/source/250611_fusion_panel_with_buffer.bed",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/bundle.yaml",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/mapping.lock.yaml",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/source/selection_panel.grch38.normalized.bed",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/source/target_gene_map.grch37.tsv",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/derived/liftover.minmatch0.99.grch37.raw.bed",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/derived/liftover.minmatch0.99.unmapped.txt",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/derived/roundtrip.hg38.bed",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/derived/roundtrip.unmapped.txt",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/derived/selection_panel.normalized.bed",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/derived/analysis_roi.bed",
    "share/ontseq/configs/panels/AML_AS_111_GRCh37_v1/derived/transcripts.tsv",
    "share/ontseq/configs/profiles/AML_LCWGS_GRCh38.yaml",
    "share/ontseq/configs/profiles/AML_AS_111_GRCh38.yaml",
    "share/ontseq/configs/profiles/AML_LCWGS_GRCh38_CANONICAL25.yaml",
    "share/ontseq/configs/profiles/AML_AS_111_GRCh38_CANONICAL25.yaml",
    "share/ontseq/configs/profiles/AML_LCWGS_GRCh37.yaml",
    "share/ontseq/configs/profiles/AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25.yaml",
    "share/ontseq/configs/profiles/AML_AS_111_GRCh37.yaml",
    "share/ontseq/configs/profiles/AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25.yaml",
    "share/ontseq/configs/sv/sniffles2.conservative.technical.yaml",
    "share/ontseq/configs/sv/cutesv.conservative.technical.yaml",
    "share/ontseq/configs/sv/sniffles2_cutesv.consensus.technical.yaml",
    "share/ontseq/configs/sv/evidence-priority.technical.yaml",
    "share/ontseq/scripts/run_qdnaseq_ace.R",
)

FORBIDDEN_PARTS = ("share/ontseq/configs/configs/",)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel", type=Path)
    args = parser.parse_args()
    with ZipFile(args.wheel) as archive:
        members = tuple(archive.namelist())
    missing = [
        suffix for suffix in REQUIRED_SUFFIXES if not any(x.endswith(suffix) for x in members)
    ]
    if missing:
        raise SystemExit(
            "wheel is missing required supported-build authority assets: " + ", ".join(missing)
        )
    forbidden = [part for part in FORBIDDEN_PARTS if any(part in item for item in members)]
    if forbidden:
        raise SystemExit(
            "wheel contains a duplicated nested configuration tree: " + ", ".join(forbidden)
        )
    print(f"Wheel resource check passed: {args.wheel}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
