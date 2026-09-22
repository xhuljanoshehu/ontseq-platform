"""Identify and correct one pinned cuteSV entry script without altering its installation.

The original 2.1.3 script divides genomic work into floating-point intervals. pysam
truncates fetch bounds while the read-start check retains the float, losing reads
at boundaries. Only the exact known script body receives this one-line correction.
"""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

STOCK_BODY_SHA256 = "e34d5f414a9d5c44fb7212cb033260c22a74e636c79d340def27230d8aa6adb3"
FIXED_BODY_SHA256 = "2ff25693aefc586002b51582a574dcce5cf09c9f61b917189e793f8329b405a3"
BUILD_ID = "ontseq-cutesv-2.1.3-integer-chunks-v1"
_OLD = b"batch_size=local_ref_len/(int(i[1]/mapped_unit)+1)"
_NEW = b"batch_size=max(1, local_ref_len//(int(i[1]/mapped_unit)+1))"


def executable_identity(executable: str) -> dict[str, str]:
    resolved = shutil.which(executable)
    if resolved is None:
        return {}
    data = Path(resolved).read_bytes()
    body = data.partition(b"\n")[2]
    body_sha256 = hashlib.sha256(body).hexdigest()
    return {
        "cutesv_source_sha256": hashlib.sha256(data).hexdigest(),
        "cutesv_source_body_sha256": body_sha256,
        "cutesv_build_id": (
            BUILD_ID if body_sha256 in {STOCK_BODY_SHA256, FIXED_BODY_SHA256} else "unqualified"
        ),
        "cutesv_execution_body_sha256": (
            FIXED_BODY_SHA256 if body_sha256 == STOCK_BODY_SHA256 else body_sha256
        ),
    }


def prepare_executable(executable: str, directory: Path, expected: dict[str, str]) -> str:
    """Stage a verified copy, preserving the installed interpreter and all imports."""
    if executable_identity(executable) != expected:
        raise ValueError("cuteSV executable changed before execution")
    resolved = shutil.which(executable)
    if resolved is None or expected.get("cutesv_source_body_sha256") != STOCK_BODY_SHA256:
        return resolved or executable
    data = Path(resolved).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected["cutesv_source_sha256"]:
        raise ValueError("cuteSV executable changed while staging")
    if data.count(_OLD) != 1:
        raise ValueError("Pinned cuteSV partition correction does not apply exactly once")
    corrected = data.replace(_OLD, _NEW)
    if hashlib.sha256(corrected.partition(b"\n")[2]).hexdigest() != FIXED_BODY_SHA256:
        raise ValueError("Corrected cuteSV script does not match its qualified body")
    # Compile as a syntax check only; do not import or execute the tool in Core.
    compile(corrected, "cuteSV.integer-chunks", "exec")
    target = directory / "cuteSV.integer-chunks"
    target.write_bytes(corrected)
    target.chmod(0o700)
    return str(target)
