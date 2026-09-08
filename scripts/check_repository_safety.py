from __future__ import annotations

import csv
import gzip
import json
import subprocess
from fnmatch import fnmatchcase
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
SENSITIVE_PATTERNS = (
    "*nanopolish*.tsv",
    "*nanopolish*.tsv.gz",
    "*methylation_calls*.tsv",
    "*methylation_calls*.tsv.gz",
    "*.methylation-mixture.json",
    "*.methylation-mixture.csv",
    "*.methylation-mixture.html",
    "*.methylation-mixture.json.gz",
    "*.methylation-mixture.csv.gz",
    "*.methylation-mixture.html.gz",
    "*modkit*.tsv",
    "*modkit*.tsv.gz",
    "*.methylation-validation.json",
    "*.methylation-validation.csv",
    "*.methylation-validation.html",
    "*.methylation-validation.json.gz",
    "*.methylation-validation.csv.gz",
    "*.methylation-validation.html.gz",
    "*.methylation-holdout.json",
    "*.methylation-holdout.csv",
    "*.methylation-holdout.html",
    "*.methylation-holdout.json.gz",
    "*.methylation-holdout.csv.gz",
    "*.methylation-holdout.html.gz",
)
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
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode == 0:
        return [Path(name) for name in completed.stdout.split("\0") if name]
    return [
        path
        for path in Path(".").rglob("*")
        if path.is_file() and not IGNORED_PARTS.intersection(path.parts)
    ]


def sensitive_content_reason(path: Path) -> str | None:
    """Inspect structured data, including renamed reports, without reading unlimited gzip data."""
    name = path.name.lower().removesuffix(".gz")
    suffix = Path(name).suffix
    if suffix not in {".tsv", ".csv", ".json", ".html", ".htm"}:
        return None
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    try:
        with opener(path, "rb") as handle:
            raw = handle.read(MAX_BYTES + 1)
    except (OSError, EOFError):
        return "cannot inspect structured data file"
    if len(raw) > MAX_BYTES:
        return "structured data exceeds 5 MiB inspection limit"
    try:
        content = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        return "cannot decode structured data file as UTF-8"
    if suffix in {".tsv", ".csv"}:
        try:
            header = next(
                csv.reader(content.splitlines(), delimiter="\t" if suffix == ".tsv" else ","), []
            )
        except csv.Error:
            return "cannot inspect delimited data header"
        columns = {column.strip().lower() for column in header}
        if {"read_name", "log_lik_ratio", "log_lik_methylated", "log_lik_unmethylated"} <= columns:
            return "Nanopolish read-level methylation data"
        if {"read_id", "ref_position", "mod_qual", "mod_code"} <= columns:
            return "modkit read-level methylation data"
        if {"read_group_budget", "replicate_seed", "estimated_source_a_fraction"} <= columns:
            return "sensitive methylation validation CSV report"
        if {"source_a_id", "source_b_id", "estimated_source_a_fraction"} <= columns:
            return "sensitive methylation mixture CSV report"
    elif suffix == ".json":
        try:
            payload = json.loads(content)
        except json.JSONDecodeError:
            return "cannot inspect invalid JSON data"
        if isinstance(payload, dict) and (
            payload.get("sensitive_output") is True
            or payload.get("estimand") == "methylation_signal_mixture_coefficient_for_source_a"
            or payload.get("estimand") == "fraction_of_retained_whole_read_groups_from_source_a"
        ):
            return "sensitive derived genomic JSON report"
    elif (
        "<title>Methylation source-fraction experiment</title>" in content
        or "<h1>Paired-source methylation mixture</h1>" in content
        or '<meta name="ontseq-sensitive-output" content="true">' in content
    ):
        return "sensitive methylation mixture HTML report"
    return None


def main() -> None:
    failures: list[str] = []
    for path in _candidate_files():
        lowered = path.name.lower()
        if any(lowered.endswith(suffix) for suffix in BANNED_SUFFIXES):
            failures.append(f"prohibited genomic-data extension: {path}")
        if any(fnmatchcase(lowered, pattern) for pattern in SENSITIVE_PATTERNS):
            failures.append(f"prohibited methylation artifact filename: {path}")
        if path.exists() and path.stat().st_size > MAX_BYTES:
            failures.append(f"file exceeds 5 MiB repository limit: {path}")
        elif path.is_file() and (reason := sensitive_content_reason(path)):
            failures.append(f"{reason}: {path}")
    if failures:
        raise SystemExit("Repository safety check failed:\n- " + "\n- ".join(failures))
    print("Repository safety check passed")


if __name__ == "__main__":
    main()
