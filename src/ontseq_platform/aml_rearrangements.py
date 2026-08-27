from __future__ import annotations

import json
from pathlib import Path

from pydantic import TypeAdapter, ValidationError

from .models import (
    AmlKnowledgeLock,
    AmlRearrangementRecord,
    FusionSupportStatus,
    GenomicEvent,
    SvValidationStatus,
)
from .reference import sha256_file


def load_aml_knowledge(resource_path: Path, lock: AmlKnowledgeLock) -> list[AmlRearrangementRecord]:
    if not resource_path.is_file():
        raise ValueError(f"AML knowledge resource is missing: {resource_path}")
    observed = sha256_file(resource_path)
    if observed != lock.sha256:
        raise ValueError(
            f"AML knowledge checksum mismatch: expected {lock.sha256}, observed {observed}"
        )
    try:
        payload = json.loads(resource_path.read_text(encoding="utf-8"))
        records = TypeAdapter(list[AmlRearrangementRecord]).validate_python(payload["records"])
    except (json.JSONDecodeError, KeyError, ValidationError) as exc:
        raise ValueError("AML knowledge resource does not satisfy its typed contract") from exc
    source_ids = {source for record in records for source in record.source_ids}
    if not source_ids.issubset(set(lock.source_ids)):
        raise ValueError("AML knowledge record cites a source absent from its lock")
    return records


def prioritize_aml_rearrangements(
    events: list[GenomicEvent],
    *,
    resource_path: Path,
    lock: AmlKnowledgeLock,
) -> list[GenomicEvent]:
    """Attach known-pattern evidence without asserting or validating a gene fusion."""
    records = load_aml_knowledge(resource_path, lock)
    result: list[GenomicEvent] = []
    for event in events:
        genes = {gene.upper() for gene in event.genes}
        matches = []
        for record in records:
            required = {gene.upper() for gene in record.genes}
            exact_match = record.pattern_type == "exact_pair" and required.issubset(genes)
            open_partner_match = (
                record.pattern_type == "gene_any_partner"
                and required.issubset(genes)
                and len(genes) >= 2
            )
            if exact_match or open_partner_match:
                matches.append(record)
        if not matches:
            result.append(event)
            continue
        matches.sort(key=lambda record: (record.pattern_type != "exact_pair", record.record_id))
        best = matches[0]
        notes = list(event.notes)
        notes.extend(
            f"Known AML rearrangement pattern candidate {record.display_name}: {record.caveat}"
            for record in matches
        )
        result.append(
            event.model_copy(
                update={
                    "aml_relevance": best.relevance,
                    "known_rearrangement": best.display_name,
                    "fusion_status": FusionSupportStatus.CANDIDATE,
                    "validation_status": SvValidationStatus.BIOLOGICALLY_PRIORITIZED,
                    "notes": notes,
                    "reportable": False,
                }
            )
        )
    return result
