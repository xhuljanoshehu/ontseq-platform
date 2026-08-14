from __future__ import annotations

import subprocess
from pathlib import Path

BANNED_SUFFIXES = (
    ".bam",
    ".bai",
    ".cram",
    ".crai",
    ".sam",
    ".pod5",
    ".fast5",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".vcf",
    ".vcf.gz",
    ".tbi",
    ".csi",
    ".bcf",
    ".bedmethyl",
    ".bigwig",
    ".bw",
)
MAX_BYTES = 5 * 1024 * 1024
IGNORED_PARTS = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    ".snakemake",
    ".venv",
    "__pycache__",
    "build",
    "dist",
    "logs",
    "results",
    "work",
}


def _candidate_files() -> list[Path]:
    completed = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode == 0:
        return [Path(line) for line in completed.stdout.splitlines() if line]
    return [
        path
        for path in Path(".").rglob("*")
        if path.is_file() and not IGNORED_PARTS.intersection(path.parts)
    ]


def main() -> None:
    failures: list[str] = []
    for path in _candidate_files():
        lowered = path.name.lower()
        if any(lowered.endswith(suffix) for suffix in BANNED_SUFFIXES):
            failures.append(f"prohibited genomic-data extension: {path}")
        if path.exists() and path.stat().st_size > MAX_BYTES:
            failures.append(f"file exceeds 5 MiB repository limit: {path}")
    if failures:
        raise SystemExit("Repository safety check failed:\n- " + "\n- ".join(failures))
    print("Repository safety check passed")


if __name__ == "__main__":
    main()
