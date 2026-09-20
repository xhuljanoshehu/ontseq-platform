from __future__ import annotations

import hashlib
import json
from enum import StrEnum
from typing import Literal

from pydantic import Field

from .models import StrictModel

SHA256 = r"^[0-9a-f]{64}$"
SOURCE_ID = r"^[A-Za-z0-9][A-Za-z0-9._-]{1,127}$"
COLLECTION_ID = r"^mce-[0-9a-f]{40}$"


class CallerEvidenceSourceKind(StrEnum):
    NATIVE_ARTIFACT = "native_artifact"
    FULL_EVIDENCE = "full_evidence"
    NORMALIZED_EVENT = "normalized_event"


class MultiCallerCollectionAddress(StrictModel):
    lane_sha256: str = Field(pattern=SHA256)
    source_kind: CallerEvidenceSourceKind
    source_id: str = Field(pattern=SOURCE_ID)
    collection_id: str = Field(pattern=COLLECTION_ID)
    research_only: Literal[True] = True


def build_multicaller_collection_address(
    *,
    lane_sha256: str,
    source_kind: CallerEvidenceSourceKind,
    source_id: str,
) -> MultiCallerCollectionAddress:
    payload = {
        "lane_sha256": lane_sha256,
        "source_kind": source_kind.value,
        "source_id": source_id,
    }
    rendered = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )
    digest = hashlib.sha256(rendered.encode("utf-8")).hexdigest()
    return MultiCallerCollectionAddress(
        lane_sha256=lane_sha256,
        source_kind=source_kind,
        source_id=source_id,
        collection_id=f"mce-{digest[:40]}",
    )
