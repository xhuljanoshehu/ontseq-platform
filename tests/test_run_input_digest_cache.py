from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from ontseq_platform.pipeline.envelope import sha256_file
from ontseq_platform.pipeline.input_digest import RunInputDigestCache


def test_same_run_reads_an_unchanged_input_once(tmp_path: Path) -> None:
    source = tmp_path / "large.bam"
    source.write_bytes(b"A" * 1024)
    cache = RunInputDigestCache()

    with patch(
        "ontseq_platform.pipeline.input_digest.sha256_file", wraps=sha256_file
    ) as digest_file:
        first = cache.digest(source)
        second = cache.digest(source)

    assert first == second
    assert first[1] is True
    assert digest_file.call_count == 1


def test_mutation_invalidates_the_run_scoped_cache(tmp_path: Path) -> None:
    source = tmp_path / "sample.bam"
    source.write_bytes(b"before")
    cache = RunInputDigestCache()
    first, stable = cache.digest(source)
    assert stable

    source.write_bytes(b"after-with-a-different-size")
    second, stable = cache.digest(source)

    assert stable
    assert first != second


def test_an_unstable_read_is_never_cached(tmp_path: Path) -> None:
    source = tmp_path / "sample.bam"
    source.write_bytes(b"initial")
    cache = RunInputDigestCache()
    calls = 0

    def changing_digest(path: Path) -> str:
        nonlocal calls
        calls += 1
        digest = sha256_file(path)
        if calls == 1:
            source.write_bytes(b"changed-after-the-read")
        return digest

    with patch("ontseq_platform.pipeline.input_digest.sha256_file", side_effect=changing_digest):
        first, first_stable = cache.digest(source)
        second, second_stable = cache.digest(source)

    assert first_stable is False
    assert second_stable is True
    assert first != second
    assert calls == 2


def test_digest_state_is_not_shared_between_run_contexts(tmp_path: Path) -> None:
    source = tmp_path / "sample.bam"
    source.write_bytes(b"same bytes")
    first_run = RunInputDigestCache()
    second_run = RunInputDigestCache()

    with patch(
        "ontseq_platform.pipeline.input_digest.sha256_file", wraps=sha256_file
    ) as digest_file:
        first_run.digest(source)
        second_run.digest(source)

    assert digest_file.call_count == 2
