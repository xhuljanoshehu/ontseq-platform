"""Refuse a tree whose declared versions disagree with each other.

Five files state what this release is: the package metadata, the package itself, the
citation record, the Windows desktop project and the README's status section. Nothing tied
them together, and they drifted — the README described the desktop on ``main`` as the
v0.2.1 engineering path while the project file, the desktop README and the changelog had
all moved to 0.3.4. The same class of drift recurred when the package and Desktop advanced
to 0.3.5 without a matching root README status or changelog release heading; this guard now
catches that mismatch in CI.

That matters more here than in most projects. A reader deciding whether a bundle is worth
testing against real data starts at the README, and the repository's own rule is that the
executable code at a commit is what counts. A stale version number in the one document
people read first sends them to check the wrong thing, and a run's provenance record names
a version that has to mean something specific.

Checked here rather than in a test because it is a property of the tree, not of the code,
and because ``make lint`` and CI already run this class of check. Deliberately not
automatic repair: which file is wrong is a decision, not something to guess at.
"""

from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

PYPROJECT = ROOT / "pyproject.toml"
PACKAGE_INIT = ROOT / "src" / "ontseq_platform" / "__init__.py"
CITATION = ROOT / "CITATION.cff"
DESKTOP_PROJECT = ROOT / "desktop" / "ONTSeq.Desktop" / "ONTSeq.Desktop.csproj"
CHANGELOG = ROOT / "CHANGELOG.md"
README = ROOT / "README.md"

#: The README section that states what the current build is. Version tokens elsewhere in
#: the README are historical context — a superseded milestone, a pinned tool — and are not
#: claims about this release.
STATUS_HEADING = "## 14. Entwicklungsstatus"

#: Matches both ``0.3.4`` and the ``v0.3.4`` form the prose uses. The optional ``v`` has to
#: be part of the pattern rather than stripped afterwards: the lookbehind would otherwise
#: treat the ``v`` as a word character and skip the match entirely, which is precisely how
#: the stale "v0.2.1" desktop claim would have slipped past this check.
SEMVER = re.compile(r"(?<![\w.])v?(\d+\.\d+\.\d+)(?![\w.])")


class Mismatch(Exception):
    """A declared version disagrees with the package metadata."""


def _declared_version() -> str:
    """The version in ``pyproject.toml``, which every other declaration must match."""
    with PYPROJECT.open("rb") as handle:
        payload = tomllib.load(handle)
    version = payload.get("project", {}).get("version")
    if not isinstance(version, str) or not version:
        raise Mismatch("pyproject.toml declares no [project] version")
    return version


def _package_version() -> str:
    """Read ``__version__`` as text rather than importing, so this runs uninstalled."""
    match = re.search(
        r'^__version__\s*=\s*["\']([^"\']+)["\']', PACKAGE_INIT.read_text(encoding="utf-8"), re.M
    )
    if match is None:
        raise Mismatch(f"{PACKAGE_INIT.relative_to(ROOT)} declares no __version__")
    return match.group(1)


def _citation_version() -> str:
    match = re.search(r"^version:\s*(\S+)\s*$", CITATION.read_text(encoding="utf-8"), re.M)
    if match is None:
        raise Mismatch("CITATION.cff declares no version")
    return match.group(1).strip("\"'")


def _desktop_version() -> str:
    match = re.search(
        r"<Version>\s*([^<\s]+)\s*</Version>", DESKTOP_PROJECT.read_text(encoding="utf-8")
    )
    if match is None:
        raise Mismatch(f"{DESKTOP_PROJECT.relative_to(ROOT)} declares no <Version>")
    return match.group(1)


def _readme_status_section() -> str:
    text = README.read_text(encoding="utf-8")
    start = text.find(STATUS_HEADING)
    if start < 0:
        raise Mismatch(f"README.md has no {STATUS_HEADING!r} section to check")
    end = text.find("\n## ", start + len(STATUS_HEADING))
    return text[start:] if end < 0 else text[start:end]


def _problems() -> list[str]:
    declared = _declared_version()
    found: list[str] = []

    for label, actual in (
        (f"{PACKAGE_INIT.relative_to(ROOT)} __version__", _package_version()),
        ("CITATION.cff version", _citation_version()),
        (f"{DESKTOP_PROJECT.relative_to(ROOT)} <Version>", _desktop_version()),
    ):
        if actual != declared:
            found.append(f"{label} is {actual}, but pyproject.toml declares {declared}")

    # The changelog may lead with an Unreleased section; what must exist is a heading for
    # the version being declared, so a released version is never undocumented.
    changelog = CHANGELOG.read_text(encoding="utf-8")
    if not re.search(rf"^## +{re.escape(declared)}\b", changelog, re.M):
        found.append(f"CHANGELOG.md has no '## {declared}' section for the declared version")

    # Every version the status section names is a claim about this build.
    stale = sorted({token for token in SEMVER.findall(_readme_status_section())} - {declared})
    if stale:
        found.append(
            f"README.md {STATUS_HEADING!r} names version(s) {', '.join(stale)} "
            f"alongside the declared {declared}"
        )
    return found


def main() -> int:
    try:
        problems = _problems()
    except (Mismatch, OSError, ValueError) as error:
        print(f"Version consistency check could not run: {error}", file=sys.stderr)
        return 2
    if problems:
        print("Version consistency check failed:", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        print(
            "\nUpdate the file that is behind. A run's provenance record names this "
            "version, so the declarations have to agree.",
            file=sys.stderr,
        )
        return 1
    print("Version consistency check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
