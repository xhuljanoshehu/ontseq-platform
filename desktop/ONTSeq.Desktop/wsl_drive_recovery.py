"""Embedded stdlib-only WSL helper. Repair only a proven stale drive mount.

Windows checks local drive availability first. Linux refuses other errors, other
sources and the system drive. Ordinary umount enforces the kernel's busy check;
there is deliberately no forced/lazy unmount or global WSL shutdown.
"""

import errno
import json
import os
import re
import subprocess
import sys
from pathlib import Path


def read_mounts():
    return Path("/proc/self/mountinfo").read_text()


def stat_mount(path):
    return os.stat(path)


def _ensure_reconstructable_drvfs_semantics(options, uid, gid):
    tokens = [token.strip() for token in re.split(r"[,;]", options) if token.strip()]

    semantic_overrides = set()
    if "metadata" in tokens:
        semantic_overrides.add("metadata")

    for key, default in (("umask", 0o22), ("fmask", 0), ("dmask", 0)):
        for token in tokens:
            if not token.startswith(key + "="):
                continue
            value = token.split("=", 1)[1]
            try:
                equivalent = int(value, 8) == default
            except ValueError:
                equivalent = False
            if not equivalent:
                semantic_overrides.add(token)

    for token in tokens:
        if token.startswith("case=") and token != "case=off":
            semantic_overrides.add(token)

    if semantic_overrides:
        rendered = ", ".join(sorted(semantic_overrides))
        raise ValueError(
            "Nichtstandard-DrvFs-Semantik kann nicht sicher rekonstruiert werden: " + rendered
        )

    for key, expected in (("uid", uid), ("gid", gid)):
        values = [token.split("=", 1)[1] for token in tokens if token.startswith(key + "=")]
        if values and any(value != str(expected) for value in values):
            raise ValueError(
                f"Bestehende DrvFs-{key}-Semantik stimmt nicht mit dem WSL-Benutzer überein."
            )


def recover(drive, uid, gid, *, repair=False):
    if not re.fullmatch("[D-Z]", drive) or uid <= 0 or gid <= 0:
        raise ValueError("Nur lokale Datenlaufwerke D: bis Z: mit normalem WSL-Benutzer.")
    target = "/mnt/" + drive.lower()
    stale = False
    try:
        stat_mount(target)
    except OSError as error:
        if error.errno != errno.ENODEV:
            raise
        stale = True
    matches = []
    for line in read_mounts().splitlines():
        left, right = line.split(" - ", 1)
        fields, filesystem = left.split(), right.split()
        mountpoint = fields[4]
        if mountpoint.startswith(target + "/"):
            raise ValueError("Untergeordnete Einbindung vorhanden; manuelle Prüfung erforderlich.")
        if mountpoint == target:
            matches.append(filesystem)
    if len(matches) != 1:
        raise ValueError("Keine eindeutige bestehende Laufwerkseinbindung gefunden.")
    kind, source, options = matches[0][:3]
    source = source.replace("\\134", "\\")
    if source not in {drive + ":", drive + ":\\"} or kind not in {"9p", "drvfs"}:
        raise ValueError("Die Einbindung gehört nicht zum angeforderten Windows-Laufwerk.")
    if kind == "9p" and "aname=drvfs;" not in options:
        raise ValueError("Keine bestätigte DrvFs-Einbindung.")
    if not stale:
        return {"state": "healthy", "drive": drive}
    if not repair:
        return {"state": "stale", "drive": drive}
    _ensure_reconstructable_drvfs_semantics(options, uid, gid)
    subprocess.run(["/bin/umount", target], check=True, capture_output=True, text=True, timeout=15)
    subprocess.run(
        ["/bin/mount", "-t", "drvfs", drive + ":", target, "-o", f"uid={uid},gid={gid}"],
        check=True,
        capture_output=True,
        text=True,
        timeout=15,
    )
    stat_mount(target)
    return {"state": "repaired", "drive": drive}


if __name__ == "__main__":
    try:
        print(
            json.dumps(
                recover(
                    sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), repair=sys.argv[4] == "repair"
                )
            )
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(str(error), file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError):
            print(error.stderr, file=sys.stderr)
        sys.exit(1)
