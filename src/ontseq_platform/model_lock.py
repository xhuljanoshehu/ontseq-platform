"""Producing the checksum the basecalling policy locks a Dorado model to.

``BasecallPolicy.model_sha256`` exists so that a partially downloaded or altered model is
refused before hours of GPU time are spent on it. Until now nothing could *produce* that
value: preflight names it in a remedy, but only once a full manifest and policy already
exist and only for the model that run happens to point at. Locking a model is a setup task
that happens before any of that, usually once per site.

The fingerprint is deliberately identical to the one :func:`ontseq_platform.basecall.
model_signature` already computes — sorted relative file names interleaved with their
individual checksums — because changing it would silently invalidate every ``model_sha256``
a site has already recorded, in a way no error message would explain.

What this adds beyond the bare digest is the material to judge whether the directory is
worth locking at all. Any directory produces a valid-looking 64-character digest, including
one holding three zero-byte files from an interrupted download. The file count, the total
size and the listed concerns are what let a human notice that before the digest goes into a
policy and starts being enforced.

**Not verified against a real Dorado model.** This repository has never had one: no GPU, no
downloaded weights. The digest is computed by the same code path preflight compares against,
so the two agree by construction, but no claim is made here about what a correct Dorado
model directory contains — hence concerns are limited to defects that are defects under any
layout, and the per-file listing is offered instead of a structural check.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from .pipeline.envelope import sha256_file


class ModelLockError(RuntimeError):
    """The path given cannot be fingerprinted at all."""


@dataclass(frozen=True)
class ModelFile:
    """One file inside the model directory, as it contributes to the digest."""

    relative_path: str
    size_bytes: int


@dataclass(frozen=True)
class ModelFingerprint:
    """A model directory's checksum together with what went into it."""

    path: Path
    signature: str
    files: tuple[ModelFile, ...]
    #: Entries that exist as symbolic links but resolve to nothing.
    broken_links: tuple[str, ...]

    @property
    def file_count(self) -> int:
        return len(self.files)

    @property
    def total_bytes(self) -> int:
        return sum(item.size_bytes for item in self.files)

    @property
    def empty_files(self) -> tuple[str, ...]:
        return tuple(item.relative_path for item in self.files if item.size_bytes == 0)

    @property
    def concerns(self) -> tuple[str, ...]:
        """Reasons not to lock this directory, in the order a reader should meet them.

        Only defects that are defects under any layout. A model with an unusual file count
        is not reported, because this repository does not know what a usual one is and a
        threshold invented here would be a guess presented as a check.
        """
        found: list[str] = []
        if not self.files:
            found.append("the directory contains no files, so the checksum describes nothing")
        if self.broken_links:
            found.append(
                f"{len(self.broken_links)} broken symbolic link(s), which contribute "
                f"nothing to the checksum and will not be noticed later: "
                + ", ".join(self.broken_links)
            )
        empty = self.empty_files
        if empty:
            found.append(
                f"{len(empty)} empty file(s), the usual shape of an interrupted download: "
                + ", ".join(empty)
            )
        return tuple(found)


def fingerprint(directory: Path) -> ModelFingerprint:
    """Fingerprint a downloaded model directory.

    A bare model name is refused rather than resolved. Dorado keeps named models in its own
    cache, and guessing where that is would produce a checksum for a directory the run may
    not use — worse than having none, because it would look like a lock.
    """
    if not directory.exists():
        raise ModelLockError(f"no such path: {directory}")
    if not directory.is_dir():
        raise ModelLockError(
            f"{directory} is not a directory; a model can only be locked from a downloaded "
            "model directory, not from a model name"
        )

    digest = hashlib.sha256()
    files: list[ModelFile] = []
    broken: list[str] = []
    # Sorted exactly as ``basecall.model_signature`` sorts, over Path objects rather than
    # over strings: the two orderings differ once a name contains a separator, and the
    # digests would then disagree for the same directory.
    for item in sorted(directory.rglob("*")):
        relative = item.relative_to(directory).as_posix()
        if item.is_symlink() and not item.exists():
            broken.append(relative)
            continue
        if not item.is_file():
            continue
        digest.update(relative.encode("utf-8"))
        digest.update(sha256_file(item).encode("ascii"))
        files.append(ModelFile(relative_path=relative, size_bytes=item.stat().st_size))
    return ModelFingerprint(
        path=directory,
        signature=digest.hexdigest(),
        files=tuple(files),
        broken_links=tuple(broken),
    )


def human_size(size_bytes: int) -> str:
    """Bytes as a human reads them, with the exact figure kept beside it.

    Both, because the rounded figure is what makes "1.8 GiB" recognisable as a real model
    and the exact one is what makes two runs comparable.
    """
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    value = float(size_bytes)
    index = 0
    while value >= 1024 and index < len(units) - 1:
        value /= 1024
        index += 1
    rounded = f"{value:.0f} {units[index]}" if index == 0 else f"{value:.1f} {units[index]}"
    return f"{rounded} ({size_bytes} bytes)"


def render(result: ModelFingerprint, *, list_files: bool = False) -> str:
    """Render the fingerprint for a person setting up a site."""
    lines = [
        f"model directory: {result.path}",
        f"files hashed:    {result.file_count}",
        f"total size:      {human_size(result.total_bytes)}",
        f"sha256:          {result.signature}",
    ]
    if list_files:
        lines.append("")
        lines.append("files, in the order they enter the checksum:")
        for item in result.files:
            lines.append(f"    {item.relative_path:<60} {item.size_bytes}")
    if result.concerns:
        lines.append("")
        lines.append("DO NOT LOCK THIS DIRECTORY YET:")
        lines.extend(f"    {item}" for item in result.concerns)
        return "\n".join(lines)
    lines.append("")
    lines.append("Record both lines in the basecalling policy, together:")
    lines.append(f"    model: {result.path}")
    lines.append(f"    model_sha256: {result.signature}")
    lines.append("")
    lines.append(
        "The checksum belongs to this directory. Moving or re-downloading the model "
        "changes nothing about the checksum but everything about what it refers to, so "
        "re-run this command rather than copying the value across."
    )
    return "\n".join(lines)


def exit_code(result: ModelFingerprint) -> int:
    """``2`` when the directory should not be locked, matching the run failure code."""
    return 2 if result.concerns else 0
