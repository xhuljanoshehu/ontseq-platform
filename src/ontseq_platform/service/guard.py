"""The parts of the local service where a wrong answer is dangerous.

A service on ``127.0.0.1`` that accepts a filesystem path and starts subprocesses is a
local attack surface, not a convenience layer. Any page open in the same browser can issue
requests to it, and a path that escapes its allowed roots turns a viewer into a file
reader. So the four decisions that carry that weight live here, apart from the transport,
with no dependency beyond the standard library — which is also what lets them be tested in
a development environment where pydantic cannot be installed.

None of this makes the service safe to expose on a network. It is a single-user tool bound
to the loopback interface, and the checks below assume that.
"""

from __future__ import annotations

import re
import secrets
from collections.abc import Sequence
from pathlib import Path, PurePosixPath, PureWindowsPath

#: Header the page presents on every request. Named rather than reused from Authorization
#: so that a proxy or browser extension forwarding credentials cannot satisfy it by
#: accident.
TOKEN_HEADER = "X-ONTSeq-Token"

#: Loopback names a browser may legitimately use to reach this service.
LOOPBACK_HOSTS = frozenset({"127.0.0.1", "localhost", "[::1]", "::1"})

_DRIVE = re.compile(r"^([A-Za-z]):[\\/](.*)$", re.DOTALL)


class GuardError(Exception):
    """A request was refused. The message is safe to show the operator."""


def new_token() -> str:
    """A fresh per-process token. Never derived from anything predictable."""
    return secrets.token_urlsafe(32)


def token_matches(presented: str | None, expected: str) -> bool:
    """Constant-time comparison, so a wrong token cannot be found by timing it."""
    if not presented:
        return False
    return secrets.compare_digest(presented, expected)


def host_is_loopback(host_header: str | None, *, port: int) -> bool:
    """Whether the ``Host`` header names this service on the loopback interface.

    This is the DNS-rebinding check. An attacker who controls a domain can point it at
    127.0.0.1 and have the victim's browser make same-origin requests here; what the
    browser then sends as ``Host`` is the attacker's domain, not a loopback name. Checking
    the header — rather than trusting that binding to 127.0.0.1 is enough — is what closes
    that, because the socket really is local in that attack.
    """
    if not host_header:
        return False
    name, _, stated_port = host_header.rpartition(":")
    if not name:
        name, stated_port = host_header, ""
    if stated_port and stated_port != str(port):
        return False
    return name.lower() in LOOPBACK_HOSTS


def origin_is_loopback(origin_header: str | None, *, port: int) -> bool:
    """Whether ``Origin`` is absent or names this service.

    Absent is allowed: browsers omit it on same-origin GETs. Present and foreign is
    refused, which is what stops another open tab from driving this one.
    """
    if origin_header is None or origin_header == "null":
        return origin_header is None
    scheme, _, rest = origin_header.partition("://")
    if scheme not in {"http", "https"} or not rest:
        return False
    return host_is_loopback(rest, port=port)


def resolve_within(candidate: str | Path, roots: Sequence[Path]) -> Path:
    """Resolve a requested path and refuse anything outside the allowed roots.

    Resolution happens *before* the comparison and follows symbolic links, because a link
    inside an allowed root pointing out of it is the obvious way past a prefix check. An
    empty root list refuses everything: a service started without stating what it may read
    should read nothing, not everything.
    """
    if not roots:
        raise GuardError("no directory is allowed; start the service with --allow-root")
    resolved = Path(candidate).expanduser().resolve()
    for root in roots:
        base = Path(root).expanduser().resolve()
        if resolved == base or base in resolved.parents:
            return resolved
    allowed = ", ".join(str(Path(root).expanduser().resolve()) for root in roots)
    raise GuardError(f"path is outside the allowed directories ({allowed}): {resolved}")


def resolve_bam_index(bam_path: str | Path) -> Path:
    """Return the supported BAM index beside *bam_path*, with deterministic precedence.

    Both naming conventions are common: ``sample.bam.bai`` and ``sample.bai``. When both
    exist, the explicit ``sample.bam.bai`` form wins. No recursive search or unrelated
    index file is accepted.
    """
    bam = Path(bam_path)
    if bam.suffix.lower() != ".bam":
        raise GuardError(f"not a BAM path: {bam}")
    preferred = Path(f"{bam}.bai")
    alternative = bam.with_suffix(".bai")
    for candidate in (preferred, alternative):
        if candidate.is_file():
            return candidate
    raise GuardError(
        f"no BAM index found; expected {preferred.name} or {alternative.name} next to the BAM"
    )


def windows_to_wsl(path: str) -> str:
    """``P:\\Lab\\run.bam`` to ``/mnt/p/Lab/run.bam``.

    Translation is done here rather than by shelling out to ``wslpath`` so it is testable
    and so a missing ``wslpath`` is not a runtime failure. UNC paths are refused with the
    reason: WSL does not map them automatically, and silently inventing a mount point
    would produce a path that does not exist.
    """
    text = path.strip().strip('"')
    if text.startswith("\\\\") or text.startswith("//"):
        raise GuardError(
            f"network path {text!r} is not reachable from WSL on its own; mount the share "
            "first (mount -t drvfs) and give the mounted path instead"
        )
    match = _DRIVE.match(text)
    if not match:
        return text
    drive, remainder = match.groups()
    parts = [part for part in PureWindowsPath(remainder).parts if part not in ("\\", "/")]
    return str(PurePosixPath(f"/mnt/{drive.lower()}", *parts))


def wsl_to_windows(path: str) -> str:
    """``/mnt/p/Lab/run.bam`` back to ``P:\\Lab\\run.bam``, for display only.

    Only for showing the operator a path they recognise. What the pipeline records is the
    path it actually opened, which is the POSIX one — a provenance record naming a path
    that was never opened is worse than one that looks unfamiliar.
    """
    parts = PurePosixPath(path).parts
    if len(parts) < 3 or parts[0] != "/" or parts[1] != "mnt" or len(parts[2]) != 1:
        return path
    return str(PureWindowsPath(f"{parts[2].upper()}:\\", *parts[3:]))
