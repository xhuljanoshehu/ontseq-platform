"""Refuse a tree whose current release identities disagree with each other.

The package metadata, dependency lock, package, citation record, Windows project, executable
Desktop label, WSL installer, Desktop workflow and current operator documentation all state
what a release is. They previously drifted independently: a new build could keep an old runtime
directory or artifact name even after the package version advanced. This guard ties those
current-release surfaces together while leaving historical changelog entries untouched.

That matters more here than in most projects. A reader deciding whether a bundle is worth
testing against real data starts at the README, and the repository's own rule is that the
executable code at a commit is what counts. A stale version number in the one document
people read first sends them to check the wrong thing, and a run's provenance record names
a version that has to mean something specific.

Checked here because it is a property of the tree and ``make versions`` and CI run this
guard directly. Synthetic tests exercise its drift detection without modifying the tree.
Deliberately not automatic repair: which file is wrong is a decision, not something to guess at.
"""

from __future__ import annotations

import json
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

PYPROJECT = ROOT / "pyproject.toml"
UV_LOCK = ROOT / "uv.lock"
PACKAGE_INIT = ROOT / "src" / "ontseq_platform" / "__init__.py"
CITATION = ROOT / "CITATION.cff"
DESKTOP_PROJECT = ROOT / "desktop" / "ONTSeq.Desktop" / "ONTSeq.Desktop.csproj"
DESKTOP_VERSION = ROOT / "desktop" / "ONTSeq.Desktop" / "Version.cs"
DESKTOP_LAUNCHER = ROOT / "desktop" / "ONTSeq.Desktop" / "WslServiceLauncher.cs"
DESKTOP_MODELS = ROOT / "desktop" / "ONTSeq.Desktop" / "Models.cs"
DESKTOP_SETTINGS_EXAMPLE = ROOT / "desktop" / "desktop.settings.example.json"
DESKTOP_CHANGELOG = ROOT / "desktop" / "ONTSeq.Desktop" / "CHANGELOG.md"
DESKTOP_README = ROOT / "desktop" / "README.md"
DESKTOP_PACKAGED_README = ROOT / "desktop" / "ONTSeq.Desktop" / "README.txt"
DESKTOP_FIRST_RUN = ROOT / "desktop" / "README-FIRST-RUN.md"
DESKTOP_ISOLATED_TESTING = ROOT / "desktop" / "ISOLATED-TESTING.md"
OPERATOR_FIRST_RUN = ROOT / "docs" / "DESKTOP_FIRST_RUN.md"
DESKTOP_API_CONTRACT = ROOT / "docs" / "DESKTOP_API_CONTRACT.md"
REFERENCE_SYSTEM = ROOT / "docs" / "REFERENCE_SYSTEM.md"
DESKTOP_WORKFLOW = ROOT / ".github" / "workflows" / "desktop-ci.yml"
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


def _uv_lock_version() -> str:
    """Validate the TOML lock and read this editable project's locked version."""
    with UV_LOCK.open("rb") as handle:
        payload = tomllib.load(handle)
    matches = [
        package
        for package in payload.get("package", [])
        if package.get("name") == "ontseq-platform" and package.get("source") == {"editable": "."}
    ]
    if len(matches) != 1:
        raise Mismatch("uv.lock must contain exactly one editable ontseq-platform package")
    version = matches[0].get("version")
    if not isinstance(version, str) or not version:
        raise Mismatch("uv.lock declares no version for editable ontseq-platform")
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


def _desktop_core_version() -> str:
    match = re.search(
        r'\bCoreValue\s*=\s*"(\d+\.\d+\.\d+)"\s*;',
        DESKTOP_VERSION.read_text(encoding="utf-8"),
    )
    if match is None:
        raise Mismatch(f"{DESKTOP_VERSION.relative_to(ROOT)} declares no literal CoreValue")
    return match.group(1)


def _require_derived_declaration(path: Path, pattern: str, label: str) -> None:
    if not re.search(pattern, path.read_text(encoding="utf-8")):
        raise Mismatch(
            f"{path.relative_to(ROOT)} {label} must derive from DesktopVersion.CoreValue"
        )


def _desktop_label_version() -> str:
    _require_derived_declaration(
        DESKTOP_VERSION, r'\bValue\s*=\s*CoreValue\s*\+\s*"-engineering"\s*;', "Value"
    )
    return _desktop_core_version()


def _desktop_runtime_version() -> str:
    _require_derived_declaration(
        DESKTOP_LAUNCHER,
        r"\bReleaseVersion\s*=\s*DesktopVersion\.CoreValue\s*;",
        "ReleaseVersion",
    )
    return _desktop_core_version()


def _desktop_resource_root() -> str:
    _require_derived_declaration(
        DESKTOP_MODELS,
        r'\bDefaultResourceRootWsl\s*=\s*"~/.local/share/ontseq/resources-v"'
        r"\s*\+\s*DesktopVersion\.CoreValue\s*;",
        "DefaultResourceRootWsl",
    )
    return f"~/.local/share/ontseq/resources-v{_desktop_core_version()}"


def _desktop_workflow_version() -> str:
    match = re.search(
        r'^\s*ONTSEQ_VERSION:\s*["\']?(\d+\.\d+\.\d+)["\']?\s*$',
        DESKTOP_WORKFLOW.read_text(encoding="utf-8"),
        re.M,
    )
    if match is None:
        raise Mismatch(f"{DESKTOP_WORKFLOW.relative_to(ROOT)} declares no ONTSEQ_VERSION")
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
        ("uv.lock editable package version", _uv_lock_version()),
        (f"{PACKAGE_INIT.relative_to(ROOT)} __version__", _package_version()),
        ("CITATION.cff version", _citation_version()),
        (f"{DESKTOP_PROJECT.relative_to(ROOT)} <Version>", _desktop_version()),
        (f"{DESKTOP_VERSION.relative_to(ROOT)} CoreValue", _desktop_core_version()),
        (
            f"{DESKTOP_VERSION.relative_to(ROOT)} engineering version",
            _desktop_label_version(),
        ),
        (
            f"{DESKTOP_LAUNCHER.relative_to(ROOT)} ReleaseVersion",
            _desktop_runtime_version(),
        ),
        (
            f"{DESKTOP_WORKFLOW.relative_to(ROOT)} ONTSEQ_VERSION",
            _desktop_workflow_version(),
        ),
    ):
        if actual != declared:
            found.append(f"{label} is {actual}, but pyproject.toml declares {declared}")

    expected_resource_root = f"~/.local/share/ontseq/resources-v{declared}"
    resource_root = _desktop_resource_root()
    if resource_root != expected_resource_root:
        found.append(
            f"Desktop resource default is {resource_root}, expected {expected_resource_root}"
        )
    with DESKTOP_SETTINGS_EXAMPLE.open(encoding="utf-8") as handle:
        example = json.load(handle)
    if not isinstance(example, dict) or example.get("resourceRootWsl") != expected_resource_root:
        found.append(
            f"{DESKTOP_SETTINGS_EXAMPLE.relative_to(ROOT)} resourceRootWsl must be "
            f"{expected_resource_root}"
        )

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

    required_markers = (
        (DESKTOP_CHANGELOG, f"## {declared}-engineering"),
        (DESKTOP_README, f"Current engineering build: Desktop/Core v{declared}."),
        (DESKTOP_README, f"## v{declared} user path"),
        (DESKTOP_README, f"runtime-v{declared}"),
        (DESKTOP_README, expected_resource_root),
        (DESKTOP_PACKAGED_README, expected_resource_root),
        (DESKTOP_FIRST_RUN, expected_resource_root),
        (DESKTOP_ISOLATED_TESTING, expected_resource_root),
        (OPERATOR_FIRST_RUN, expected_resource_root),
        (DESKTOP_API_CONTRACT, expected_resource_root),
        (REFERENCE_SYSTEM, expected_resource_root),
        (DESKTOP_FIRST_RUN, f"# ONTSeq Desktop v{declared}"),
        (
            DESKTOP_FIRST_RUN,
            f"ontseq-desktop-v{declared}-win-x64-setup-engineering",
        ),
        (DESKTOP_WORKFLOW, "name: ontseq-linux-runtime-v${{ env.ONTSEQ_VERSION }}"),
        (
            DESKTOP_WORKFLOW,
            "name: ontseq-desktop-v${{ env.ONTSEQ_VERSION }}-win-x64-setup-engineering",
        ),
    )
    for path, marker in required_markers:
        if marker not in path.read_text(encoding="utf-8"):
            found.append(f"{path.relative_to(ROOT)} is missing current-release marker {marker!r}")
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
