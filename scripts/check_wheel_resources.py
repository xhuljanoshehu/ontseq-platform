#!/usr/bin/env python3
"""Verify that an ONTSeq wheel carries the immutable GRCh38 authority assets."""

from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZipFile

REQUIRED_SUFFIXES = (
    "share/ontseq/configs/reference_bundles/GRCh38_GENCODE50_MANE1.5_v1/bundle.recipe.yaml",
    "share/ontseq/configs/knowledge_bundles/HEMATOLOGY_v1/bundle.yaml",
    "share/ontseq/configs/panels/AML_AS_111_GRCh38_v1/bundle.yaml",
    "share/ontseq/configs/panels/AML_AS_111_GRCh38_v1/source/250611_fusion_panel_with_buffer.bed",
    "share/ontseq/configs/profiles/AML_LCWGS_GRCh38.yaml",
    "share/ontseq/configs/profiles/AML_AS_111_GRCh38.yaml",
)


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
        raise SystemExit("wheel is missing required GRCh38 authority assets: " + ", ".join(missing))
    print(f"Wheel resource check passed: {args.wheel}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
