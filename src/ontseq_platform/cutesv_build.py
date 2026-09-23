"""Qualify a pinned cuteSV script and stage bounded correctness fixes.

The known 2.1.3 script loses reads at fractional task boundaries and can swallow
worker exceptions during extraction and clustering. Only exact reviewed bodies
receive the integer-boundary and exception-propagation fixes; installation files
are never changed. This is a technical qualification, not biological validation.
"""

from __future__ import annotations

import hashlib
import re
import shlex
import shutil
from pathlib import Path

STOCK_BODY_SHA256 = "e34d5f414a9d5c44fb7212cb033260c22a74e636c79d340def27230d8aa6adb3"
INTEGER_BODY_SHA256 = "2ff25693aefc586002b51582a574dcce5cf09c9f61b917189e793f8329b405a3"
FIXED_BODY_SHA256 = "205ad43af1b9eb9cd2d6a18a6e6c84211d519de7e4e84c7a22c4f94d4d0ffe89"
BUILD_ID = "ontseq-cutesv-2.1.3-worker-errors-v2"
LAUNCH_CONTRACT = "python-shebang-v1"
_PYTHON_INTERPRETER = re.compile(r"^python(?:3(?:\.\d+)?)?$")
_OLD = b"batch_size=local_ref_len/(int(i[1]/mapped_unit)+1)"
_NEW = b"batch_size=max(1, local_ref_len//(int(i[1]/mapped_unit)+1))"


_WORKER_PATCHES = (
    (
        b"        analysis_pools.map_async(multi_run_wrapper, paras)",
        b"        results.append(analysis_pools.map_async(multi_run_wrapper, paras))",
    ),
    (
        b"    analysis_pools.join()\n"
        b'    logging.info("Rebuilding signatures of structural variants.")',
        b"    analysis_pools.join()\n    for signature_result in results:\n"
        b"        signature_result.get()\n"
        b'    logging.info("Rebuilding signatures of structural variants.")',
    ),
    (
        b'            except:\n                logging.info("LocalError: {e}")',
        b"            except Exception:\n"
        b'                logging.exception("cuteSV clustering worker failed")\n'
        b"                raise",
    ),
)


def _launch_identity(first_line: bytes) -> tuple[str, str]:
    """Classify a portable Python shebang without pinning an installation-specific path."""
    try:
        text = first_line.decode("utf-8")
    except UnicodeDecodeError:
        return "unqualified", "unqualified"
    if not text.startswith("#!"):
        return "unqualified", "unqualified"
    try:
        parts = shlex.split(text[2:].strip(), posix=True)
    except ValueError:
        return "unqualified", "unqualified"
    if not parts:
        return "unqualified", "unqualified"
    command = Path(parts[0]).name
    if command == "env":
        if len(parts) == 2 and _PYTHON_INTERPRETER.fullmatch(parts[1]):
            return LAUNCH_CONTRACT, parts[1]
        return "unqualified", "unqualified"
    if len(parts) == 1 and _PYTHON_INTERPRETER.fullmatch(command):
        return LAUNCH_CONTRACT, parts[0]
    return "unqualified", "unqualified"


def executable_identity(executable: str) -> dict[str, str]:
    resolved = shutil.which(executable)
    if resolved is None:
        return {}
    data = Path(resolved).read_bytes()
    first_line, _separator, body = data.partition(b"\n")
    body_sha256 = hashlib.sha256(body).hexdigest()
    launch_contract, launch_interpreter = _launch_identity(first_line)
    return {
        "cutesv_source_sha256": hashlib.sha256(data).hexdigest(),
        "cutesv_source_body_sha256": body_sha256,
        "cutesv_source_shebang_sha256": hashlib.sha256(first_line).hexdigest(),
        "cutesv_launch_contract": launch_contract,
        "cutesv_launch_interpreter": launch_interpreter,
        "cutesv_build_id": (
            BUILD_ID
            if body_sha256 in {STOCK_BODY_SHA256, INTEGER_BODY_SHA256, FIXED_BODY_SHA256}
            else "unqualified"
        ),
        "cutesv_execution_body_sha256": (
            FIXED_BODY_SHA256
            if body_sha256 in {STOCK_BODY_SHA256, INTEGER_BODY_SHA256}
            else body_sha256
        ),
    }


def require_standard_lane_qualification(identity: dict[str, str]) -> None:
    """Reject a resolved cuteSV executable outside the reviewed body/launch contract."""
    if identity and (
        identity.get("cutesv_build_id") != BUILD_ID
        or identity.get("cutesv_execution_body_sha256") != FIXED_BODY_SHA256
        or identity.get("cutesv_launch_contract") != LAUNCH_CONTRACT
    ):
        raise ValueError("cuteSV executable body/launch chain is not qualified for the standard SV lane")


def prepare_executable(executable: str, directory: Path, expected: dict[str, str]) -> str:
    """Stage a verified and qualified copy for the standard SV lane."""
    if executable_identity(executable) != expected:
        raise ValueError("cuteSV executable changed before execution")
    require_standard_lane_qualification(expected)
    resolved = shutil.which(executable)
    if resolved is None or expected.get("cutesv_source_body_sha256") not in {
        STOCK_BODY_SHA256,
        INTEGER_BODY_SHA256,
    }:
        return resolved or executable
    data = Path(resolved).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected["cutesv_source_sha256"]:
        raise ValueError("cuteSV executable changed while staging")
    corrected = data
    patches: tuple[tuple[bytes, bytes], ...] = _WORKER_PATCHES
    if expected["cutesv_source_body_sha256"] == STOCK_BODY_SHA256:
        patches = ((_OLD, _NEW), *patches)
    for old, new in patches:
        if corrected.count(old) != 1:
            raise ValueError("Pinned cuteSV correction does not apply exactly once")
        corrected = corrected.replace(old, new)
    if hashlib.sha256(corrected.partition(b"\n")[2]).hexdigest() != FIXED_BODY_SHA256:
        raise ValueError("Corrected cuteSV script does not match its qualified body")
    # Compile as a syntax check only; do not import or execute the tool in Core.
    compile(corrected, "cuteSV.qualified", "exec")
    target = directory / "cuteSV.qualified"
    target.write_bytes(corrected)
    target.chmod(0o700)
    return str(target)
