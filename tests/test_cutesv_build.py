import hashlib
from pathlib import Path

import pytest

from ontseq_platform import cutesv_build


def test_unknown_script_is_identified_but_never_rewritten(tmp_path: Path) -> None:
    script = tmp_path / "cuteSV"
    data = b"#!/bin/sh\nexit 0\n"
    script.write_bytes(data)
    script.chmod(0o700)
    identity = cutesv_build.executable_identity(str(script))
    assert identity["cutesv_source_sha256"] == hashlib.sha256(data).hexdigest()
    assert identity["cutesv_build_id"] == "unqualified"
    assert cutesv_build.prepare_executable(str(script), tmp_path, identity) == str(script)
    assert script.read_bytes() == data


def test_replaced_source_is_refused_before_staging(tmp_path: Path) -> None:
    script = tmp_path / "cuteSV"
    script.write_text("#!/bin/sh\nexit 0\n")
    script.chmod(0o700)
    identity = cutesv_build.executable_identity(str(script))
    script.write_text("#!/bin/sh\nexit 1\n")
    with pytest.raises(ValueError, match="changed before execution"):
        cutesv_build.prepare_executable(str(script), tmp_path, identity)


def test_correction_produces_integer_chunks_without_touching_source(
    tmp_path: Path, monkeypatch
) -> None:
    # Synthetic qualification fixture: exercise the transform and semantic result
    # without depending on an installed biological executable in ordinary CI.
    script = tmp_path / "cuteSV"
    body = (
        b"def chunk(local_ref_len, i, mapped_unit):\n"
        b"    batch_size=local_ref_len/(int(i[1]/mapped_unit)+1)\n"
        b"    return batch_size\n"
    )
    fixed = (
        b"def chunk(local_ref_len, i, mapped_unit):\n"
        b"    batch_size=max(1, local_ref_len//(int(i[1]/mapped_unit)+1))\n"
        b"    return batch_size\n"
    )
    # This small fixture qualifies only the partition transform; real worker
    # failure propagation is exercised in test_cutesv_worker_failures.py.
    monkeypatch.setattr(cutesv_build, "_WORKER_PATCHES", ())
    monkeypatch.setattr(cutesv_build, "STOCK_BODY_SHA256", hashlib.sha256(body).hexdigest())
    monkeypatch.setattr(cutesv_build, "FIXED_BODY_SHA256", hashlib.sha256(fixed).hexdigest())
    source = b"#!/usr/bin/python3\n" + body
    script.write_bytes(source)
    script.chmod(0o700)
    identity = cutesv_build.executable_identity(str(script))
    prepared = Path(cutesv_build.prepare_executable(str(script), tmp_path, identity))
    assert prepared != script
    namespace = {}
    exec(compile(prepared.read_bytes(), "synthetic-chunk", "exec"), namespace)
    assert namespace["chunk"](200_000, ["chr1", 80], 2) == 4878
    assert isinstance(namespace["chunk"](200_000, ["chr1", 80], 2), int)
    assert namespace["chunk"](1, ["chr1", 80], 2) == 1
    assert script.read_bytes() == source
    # Changing any body byte prevents the correction from being silently applied.
    script.write_bytes(source + b"# changed\n")
    assert cutesv_build.executable_identity(str(script))["cutesv_build_id"] == "unqualified"
