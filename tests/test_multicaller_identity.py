from __future__ import annotations

import hashlib
import unittest

from ontseq_platform.multicaller_identity import (
    CallerEvidenceSourceKind,
    build_multicaller_collection_address,
)
from pydantic import ValidationError


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


class MultiCallerCollectionIdentityTests(unittest.TestCase):
    def test_same_source_id_in_different_lanes_has_different_collection_identity(self) -> None:
        first = build_multicaller_collection_address(
            lane_sha256=_sha("lane-one"),
            source_kind=CallerEvidenceSourceKind.FULL_EVIDENCE,
            source_id="qdnaseq-fit-500-selected",
        )
        second = build_multicaller_collection_address(
            lane_sha256=_sha("lane-two"),
            source_kind=CallerEvidenceSourceKind.FULL_EVIDENCE,
            source_id="qdnaseq-fit-500-selected",
        )

        self.assertNotEqual(first.collection_id, second.collection_id)
        self.assertEqual(first.source_id, "qdnaseq-fit-500-selected")
        self.assertEqual(second.source_id, "qdnaseq-fit-500-selected")

    def test_collection_identity_is_deterministic(self) -> None:
        kwargs = {
            "lane_sha256": _sha("lane-one"),
            "source_kind": CallerEvidenceSourceKind.NATIVE_ARTIFACT,
            "source_id": "segments-500",
        }

        first = build_multicaller_collection_address(**kwargs)
        second = build_multicaller_collection_address(**kwargs)

        self.assertEqual(first.collection_id, second.collection_id)
        self.assertEqual(first, second)

    def test_source_kind_participates_in_collection_identity(self) -> None:
        native = build_multicaller_collection_address(
            lane_sha256=_sha("lane-one"),
            source_kind=CallerEvidenceSourceKind.NATIVE_ARTIFACT,
            source_id="same-source-id",
        )
        normalized = build_multicaller_collection_address(
            lane_sha256=_sha("lane-one"),
            source_kind=CallerEvidenceSourceKind.NORMALIZED_EVENT,
            source_id="same-source-id",
        )

        self.assertNotEqual(native.collection_id, normalized.collection_id)

    def test_existing_source_id_is_preserved_without_rewriting(self) -> None:
        original = "full-segment-500-1"
        address = build_multicaller_collection_address(
            lane_sha256=_sha("lane-one"),
            source_kind=CallerEvidenceSourceKind.FULL_EVIDENCE,
            source_id=original,
        )

        self.assertEqual(address.source_id, original)
        self.assertNotEqual(address.collection_id, original)

    def test_invalid_lane_hash_and_source_id_fail_closed(self) -> None:
        with self.assertRaises(ValidationError):
            build_multicaller_collection_address(
                lane_sha256="not-a-sha",
                source_kind=CallerEvidenceSourceKind.FULL_EVIDENCE,
                source_id="valid-source",
            )

        with self.assertRaises(ValidationError):
            build_multicaller_collection_address(
                lane_sha256=_sha("lane-one"),
                source_kind=CallerEvidenceSourceKind.FULL_EVIDENCE,
                source_id="bad source with spaces",
            )


if __name__ == "__main__":
    unittest.main()
