"""Fail-closed filesystem checks for resource-management mutation boundaries."""

from __future__ import annotations

import os
import stat
from pathlib import Path


def path_lexists(path: Path) -> bool:
    """Return whether a directory entry exists, including a dangling symbolic link."""

    return os.path.lexists(path)


def is_link_like(path: Path) -> bool:
    """Recognize symbolic links and Windows reparse points without following them.

    ``Path.is_junction`` was added in Python 3.12, while ONTSeq also supports Python 3.11.
    The ``lstat`` fallback therefore rejects junctions (and other reparse-point directory
    entries) on every supported Windows runtime.
    """

    if not path_lexists(path):
        return False
    junction_check = getattr(path, "is_junction", None)
    if path.is_symlink() or (junction_check is not None and junction_check()):
        return True
    if os.name != "nt":
        return False
    try:
        attributes = getattr(path.lstat(), "st_file_attributes", 0)
    except OSError:
        # The entry changed between lexists and lstat. Fail closed at a mutation boundary.
        return True
    reparse_point = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
    return bool(attributes & reparse_point)


def assert_safe_descendant(root: Path, target: Path, *, label: str = "resource path") -> None:
    """Require a lexically contained path with no existing link-like components.

    ``root`` is expected to be the caller's already-canonical resource root.  Both lexical and
    resolved containment are checked.  The component walk uses ``lexists`` so dangling links are
    rejected instead of being misclassified as absent destinations.
    """

    root = Path(os.path.abspath(root))
    target = Path(os.path.abspath(target))
    try:
        relative = target.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"{label} is outside the resource root: {target}") from exc

    if path_lexists(root) and is_link_like(root):
        raise ValueError(f"{label} has a link-like resource root: {root}")
    current = root
    for component in relative.parts:
        current /= component
        if path_lexists(current) and is_link_like(current):
            raise ValueError(f"{label} contains a symbolic link or junction: {current}")

    resolved_root = root.resolve(strict=False)
    resolved_target = target.resolve(strict=False)
    try:
        resolved_target.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"{label} resolves outside the resource root: {target}") from exc


def assert_plain_tree(root: Path, *, label: str = "resource tree") -> None:
    """Reject a tree root or descendant that is a link, junction, or special entry."""

    if not path_lexists(root):
        raise FileNotFoundError(f"{label} does not exist: {root}")
    if is_link_like(root) or not root.is_dir():
        raise ValueError(f"{label} root is not a plain directory: {root}")

    pending = [root]
    while pending:
        directory = pending.pop()
        with os.scandir(directory) as entries:
            for entry in entries:
                path = Path(entry.path)
                if entry.is_symlink() or is_link_like(path):
                    raise ValueError(f"{label} contains a symbolic link or junction: {path}")
                if entry.is_dir(follow_symlinks=False):
                    pending.append(path)
                elif not entry.is_file(follow_symlinks=False):
                    raise ValueError(f"{label} contains a special filesystem entry: {path}")


__all__ = [
    "assert_plain_tree",
    "assert_safe_descendant",
    "is_link_like",
    "path_lexists",
]
