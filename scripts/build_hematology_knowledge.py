#!/usr/bin/env python3
"""Build the checked-in hematology review knowledge from pinned public-source snapshots.

The output is intentionally panel-scoped and review-only.  A CIViC fusion is retained when
at least one partner is represented by the supplied panel and at least one associated CIViC
disease is a Disease Ontology descendant of ``hematologic cancer`` (DOID:2531).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

HEMATOLOGIC_CANCER = "DOID:2531"
CIVIC_SOURCE_ID = "CIVIC-2026-08-29"
DO_SOURCE_ID = "DO-2026-07-31"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_obo(path: Path) -> tuple[dict[str, str], dict[str, tuple[str, ...]]]:
    names: dict[str, str] = {}
    parents: dict[str, tuple[str, ...]] = {}
    current: dict[str, Any] | None = None

    def commit() -> None:
        nonlocal current
        if current is not None and "id" in current:
            term_id = str(current["id"])
            names[term_id] = str(current.get("name", ""))
            parents[term_id] = tuple(current.get("parents", ()))
        current = None

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line == "[Term]":
            commit()
            current = {}
            continue
        if line.startswith("["):
            commit()
            continue
        if current is None:
            continue
        if line.startswith("id: DOID:"):
            current["id"] = line[4:].strip()
        elif line.startswith("name: "):
            current["name"] = line[6:].strip()
        elif line.startswith("is_a: DOID:"):
            current.setdefault("parents", []).append(line[6:].split()[0])
        elif not line:
            commit()
    commit()
    return names, parents


def _is_descendant(
    term_id: str,
    parents: dict[str, tuple[str, ...]],
    *,
    root: str = HEMATOLOGIC_CANCER,
    visiting: frozenset[str] = frozenset(),
) -> bool:
    if term_id == root:
        return True
    if term_id in visiting:
        return False
    return any(
        _is_descendant(parent, parents, root=root, visiting=visiting | {term_id})
        for parent in parents.get(term_id, ())
    )


def _panel_genes(path: Path) -> set[str]:
    genes: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) >= 4:
            genes.add(columns[3].removesuffix("_REVIEW_REQUIRED").upper())
    return genes


def _fusion_nodes(paths: Iterable[Path]) -> list[dict[str, Any]]:
    nodes: list[dict[str, Any]] = []
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        nodes.extend(payload["data"]["fusions"]["nodes"])
    ids = [int(node["id"]) for node in nodes]
    if len(ids) != len(set(ids)):
        raise ValueError("CIViC fusion snapshots contain duplicate feature IDs")
    return nodes


def _pathology(
    disease: dict[str, Any],
    *,
    fusion: dict[str, Any],
) -> dict[str, Any]:
    disease_id = f"DOID:{disease['doid']}"
    return {
        "disease_id": disease_id,
        "name": disease["name"],
        "source_id": CIVIC_SOURCE_ID,
        "source_record_id": f"CIVIC-FUSION-{fusion['id']}/DISEASE-{disease['id']}",
        "source_url": f"https://civicdb.org{fusion['link']}",
        "evidence_item_count": int(disease["evidenceItemCount"]),
        "assertion_count": int(disease["assertionCount"]),
    }


def _civic_records(
    *,
    fusions: list[dict[str, Any]],
    disease_payload: dict[str, Any],
    panel_genes: set[str],
    parents: dict[str, tuple[str, ...]],
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for fusion in fusions:
        gene_a = str((fusion.get("fivePrimeGene") or {}).get("name", "")).upper()
        gene_b = str((fusion.get("threePrimeGene") or {}).get("name", "")).upper()
        if not gene_a or not gene_b or not ({gene_a, gene_b} & panel_genes):
            continue
        if fusion.get("deprecated"):
            continue
        if fusion.get("fivePrimePartnerStatus") != "KNOWN":
            continue
        if fusion.get("threePrimePartnerStatus") != "KNOWN":
            continue
        pathologies = []
        for disease in disease_payload.get(f"d{fusion['id']}", {}).get("nodes", []):
            if not disease.get("doid"):
                continue
            disease_id = f"DOID:{disease['doid']}"
            if _is_descendant(disease_id, parents):
                pathologies.append(_pathology(disease, fusion=fusion))
        if not pathologies:
            continue
        pathologies.sort(key=lambda item: (item["disease_id"], item["name"]))
        records.append(
            {
                "record_id": f"CIVIC-FUSION-{fusion['id']}",
                "pattern_type": "exact_pair",
                "genes": [gene_a, gene_b],
                "display_name": f"{gene_a}::{gene_b}",
                "relevance": "hematology_relevant_pattern",
                "source_ids": [CIVIC_SOURCE_ID, DO_SOURCE_ID],
                "pathologies": pathologies,
                "caveat": (
                    "CIViC disease association plus gene overlap is a review-prioritization "
                    "pattern only; it does not establish a productive fusion, lineage, "
                    "diagnosis, prognosis, or reportability."
                ),
            }
        )
    return sorted(records, key=lambda item: (item["display_name"], item["record_id"]))


def _merge_with_curated(
    public_records: list[dict[str, Any]], base_payload: dict[str, Any]
) -> list[dict[str, Any]]:
    by_pair = {
        frozenset(record["genes"]): record
        for record in public_records
        if record["pattern_type"] == "exact_pair"
    }
    merged: list[dict[str, Any]] = []
    consumed: set[frozenset[str]] = set()
    for record in base_payload["records"]:
        candidate = dict(record)
        candidate.setdefault("pathologies", [])
        if candidate["display_name"] == "PICALM::MLLT10":
            candidate["source_ids"] = [*candidate["source_ids"], DO_SOURCE_ID]
            candidate["pathologies"] = [
                {
                    "disease_id": "DOID:9119",
                    "name": "Acute Myeloid Leukemia",
                    "source_id": "BOREL-2012",
                    "source_record_id": "PMID:22871473",
                    "source_url": "https://pubmed.ncbi.nlm.nih.gov/22871473/",
                }
            ]
        pair = frozenset(candidate["genes"])
        public = by_pair.get(pair) if candidate["pattern_type"] == "exact_pair" else None
        if public is not None:
            consumed.add(pair)
            candidate["source_ids"] = list(
                dict.fromkeys([*candidate["source_ids"], *public["source_ids"]])
            )
            candidate["pathologies"] = public["pathologies"]
        merged.append(candidate)
    merged.extend(record for pair, record in by_pair.items() if pair not in consumed)
    return sorted(
        merged,
        key=lambda item: (
            item["pattern_type"] != "exact_pair",
            item["display_name"],
            item["record_id"],
        ),
    )


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def build(arguments: argparse.Namespace) -> None:
    _, parents = _parse_obo(arguments.disease_ontology)
    disease_payload = json.loads(arguments.civic_diseases.read_text(encoding="utf-8"))["data"]
    base_payload = json.loads(arguments.base_resource.read_text(encoding="utf-8"))
    public_records = _civic_records(
        fusions=_fusion_nodes(arguments.civic_fusions),
        disease_payload=disease_payload,
        panel_genes=_panel_genes(arguments.panel_bed),
        parents=parents,
    )
    records = _merge_with_curated(public_records, base_payload)
    pathology_count = sum(len(record.get("pathologies", [])) for record in records)
    source_fusion_sha256 = [_sha256(path) for path in arguments.civic_fusions]
    payload = {
        "schema_version": "0.2.0",
        "scope": {
            "panel": "AML_AS_111_GRCh38_v1",
            "ontology_root": HEMATOLOGIC_CANCER,
            "selection": (
                "Exact CIViC fusion with at least one panel gene and at least one associated "
                "disease descending from DOID:2531; plus locally curated AML patterns."
            ),
            "record_count": len(records),
            "pathology_association_count": pathology_count,
        },
        "records": records,
        "sources": {
            **base_payload["sources"],
            CIVIC_SOURCE_ID: {
                "title": "Clinical Interpretation of Variants in Cancer (CIViC)",
                "url": "https://civicdb.org/",
                "license": "CC0-1.0",
                "retrieved": "2026-08-29",
                "fusion_snapshot_sha256": source_fusion_sha256,
                "disease_snapshot_sha256": _sha256(arguments.civic_diseases),
            },
            DO_SOURCE_ID: {
                "title": "Human Disease Ontology",
                "release": "v2026-07-31",
                "url": (
                    "https://github.com/DiseaseOntology/HumanDiseaseOntology/releases/tag/"
                    "v2026-07-31"
                ),
                "license": "CC0-1.0",
                "sha256": _sha256(arguments.disease_ontology),
                "ontology_root": HEMATOLOGIC_CANCER,
            },
        },
    }
    _write_json(arguments.output_resource, payload)
    lock = {
        "schema_version": "0.1.0",
        "resource_id": "ontseq-hematology-rearrangements",
        "release": "0.3.0",
        "sha256": _sha256(arguments.output_resource),
        "source_ids": list(payload["sources"]),
        "note": (
            f"Panel-scoped hematology review knowledge: {len(records)} patterns and "
            f"{pathology_count} source-attributed pathology associations. Matching is gene-pair "
            "order independent. No record is a diagnostic, prognostic, or reportability rule."
        ),
    }
    _write_json(arguments.output_lock, lock)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-bed", type=Path, required=True)
    parser.add_argument("--civic-fusions", type=Path, action="append", required=True)
    parser.add_argument("--civic-diseases", type=Path, required=True)
    parser.add_argument("--disease-ontology", type=Path, required=True)
    parser.add_argument("--base-resource", type=Path, required=True)
    parser.add_argument("--output-resource", type=Path, required=True)
    parser.add_argument("--output-lock", type=Path, required=True)
    return parser


if __name__ == "__main__":
    build(_parser().parse_args())
