"""Opt-in, local GATK Mutect2 research adapter; not an ONT-validated caller.

Invoke with python -m ontseq_platform.gatk_mutect2. The default is a command plan.
This module deliberately does not change the existing CNV/SV planner, clinical
report contract, or germline small-variant policy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .gatk_vcf import read_candidates

GATK_VERSION = "4.6.2.0"
SCHEMA_VERSION = "gatk-mutect2-evidence/0.1.0"
Runner = Callable[[list[str], Path, int], str]


class GATKError(ValueError):
    """A blocked or failed research run, never a negative biological finding."""


@dataclass(frozen=True)
class LockedFile:
    path: str
    sha256: str

    def local_path(self) -> str:
        if "://" in self.path or not Path(self.path).is_absolute():
            raise GATKError("Only absolute local input paths are accepted")
        return str(Path(self.path).resolve())

    def verify(self) -> None:
        if not re.fullmatch(r"[0-9a-f]{64}", self.sha256):
            raise GATKError("Invalid SHA-256 lock")
        path = Path(self.local_path())
        if not path.is_file() or not path.stat().st_size or digest(path) != self.sha256:
            raise GATKError(f"Missing, empty or checksum-mismatched input: {path}")


@dataclass(frozen=True)
class Sample:
    name: str
    bam: LockedFile
    index: LockedFile
    reference_sha256: str


@dataclass(frozen=True)
class VCFResource:
    vcf: LockedFile
    index: LockedFile
    reference_sha256: str
    technology: str = "not_applicable"
    assay_id: str = ""


@dataclass(frozen=True)
class GATKConfig:
    profile_id: str
    reference: LockedFile
    reference_fai: LockedFile
    reference_dict: LockedFile
    tumor: Sample
    intervals: LockedFile
    normal: Sample | None = None
    germline: VCFResource | None = None
    panel_of_normals: VCFResource | None = None
    population_sites: VCFResource | None = None
    threads: int = 2
    java_memory_gb: int = 4
    timeout_seconds: int = 86400
    gatk: str = "gatk"
    samtools: str = "samtools"

    @property
    def mode(self) -> str:
        return "tumor_normal" if self.normal else "tumor_only"


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _samples(config: GATKConfig) -> tuple[Sample, ...]:
    return (config.tumor, config.normal) if config.normal else (config.tumor,)


def _resources(config: GATKConfig) -> tuple[VCFResource, ...]:
    return tuple(
        x for x in (config.germline, config.panel_of_normals, config.population_sites) if x
    )


def _files(config: GATKConfig) -> list[LockedFile]:
    files = [config.reference, config.reference_fai, config.reference_dict, config.intervals]
    for sample in _samples(config):
        files.extend((sample.bam, sample.index))
    for resource in _resources(config):
        files.extend((resource.vcf, resource.index))
    return files


def _validate_config(config: GATKConfig) -> None:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", config.profile_id):
        raise GATKError("Invalid profile ID")
    for value in (config.threads, config.java_memory_gb, config.timeout_seconds):
        if type(value) is not int or value < 1:
            raise GATKError("Threads, memory and timeout must be positive integers")
    ref = config.reference.local_path()
    if config.reference_fai.local_path() != ref + ".fai":
        raise GATKError("FASTA index must be beside the reference")
    if config.reference_dict.local_path() != str(Path(ref).with_suffix(".dict")):
        raise GATKError("Sequence dictionary must be beside the reference")
    for sample in _samples(config):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", sample.name):
            raise GATKError("Sample names must be pseudonymous stable IDs")
        bam = sample.bam.local_path()
        if not bam.endswith(".bam") or sample.index.local_path() != bam + ".bai":
            raise GATKError("A BAM and its adjacent .bam.bai are required")
        if sample.reference_sha256 != config.reference.sha256:
            raise GATKError("Declared alignment reference differs from the reference lock")
    if config.normal and (
        config.normal.name == config.tumor.name
        or config.normal.bam.local_path() == config.tumor.bam.local_path()
        or config.normal.bam.sha256 == config.tumor.bam.sha256
    ):
        raise GATKError("Tumor and normal must be distinct samples and BAMs")
    for resource in _resources(config):
        vcf = resource.vcf.local_path()
        if not vcf.endswith(".vcf.gz") or resource.index.local_path() != vcf + ".tbi":
            raise GATKError("Resources require a .vcf.gz and adjacent .tbi")
        if resource.reference_sha256 != config.reference.sha256:
            raise GATKError("Resource reference differs from the reference lock")
    pon = config.panel_of_normals
    if pon and (pon.technology != "ONT" or pon.assay_id != config.profile_id):
        raise GATKError("PoN must be explicitly declared ONT and assay-matched")


def build_commands(config: GATKConfig, output_dir: Path) -> list[list[str]]:
    """Build shell-free commands without executing anything or acquiring resources."""
    _validate_config(config)
    out = output_dir.resolve()
    prefix = [config.gatk, "--java-options", f"-Xmx{config.java_memory_gb}g"]
    raw = str(out / "unfiltered.vcf.gz")
    mutect = prefix + ["Mutect2", "-R", config.reference.local_path()]
    for sample in _samples(config):
        mutect += ["-I", sample.bam.local_path()]
    if config.normal:
        mutect += ["-normal", config.normal.name]
    mutect += [
        "-L",
        config.intervals.local_path(),
        "--native-pair-hmm-threads",
        str(config.threads),
    ]
    for flag, resource in (
        ("--germline-resource", config.germline),
        ("--panel-of-normals", config.panel_of_normals),
    ):
        if resource:
            mutect += [flag, resource.vcf.local_path()]
    mutect += ["-O", raw]
    commands = [mutect]
    filtered = prefix + [
        "FilterMutectCalls",
        "-R",
        config.reference.local_path(),
        "-V",
        raw,
        "--stats",
        raw + ".stats",
    ]
    if config.population_sites:
        for role, pileup_sample in (("tumor", config.tumor), ("normal", config.normal)):
            if pileup_sample:
                commands.append(
                    prefix
                    + [
                        "GetPileupSummaries",
                        "-R",
                        config.reference.local_path(),
                        "-I",
                        pileup_sample.bam.local_path(),
                        "-V",
                        config.population_sites.vcf.local_path(),
                        "-L",
                        config.intervals.local_path(),
                        "-O",
                        str(out / f"{role}.pileups.table"),
                    ]
                )
        contamination = prefix + ["CalculateContamination", "-I", str(out / "tumor.pileups.table")]
        if config.normal:
            contamination += ["--matched-normal", str(out / "normal.pileups.table")]
        contamination += ["-O", str(out / "contamination.table")]
        commands.append(contamination)
        filtered += ["--contamination-table", str(out / "contamination.table")]
    commands.append(filtered + ["-O", str(out / "filtered.vcf.gz")])
    return commands


def subprocess_runner(argv: list[str], cwd: Path, timeout: int) -> str:
    completed = subprocess.run(
        argv,
        cwd=cwd,
        timeout=timeout,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=False,
    )
    return completed.stdout


def _sam_header(text: str) -> tuple[dict[str, int], set[str], bool]:
    contigs: dict[str, int] = {}
    names: set[str] = set()
    sorted_by_coordinate = False
    for line in text.splitlines():
        fields = line.split("\t")
        tags = dict(x.split(":", 1) for x in fields[1:] if ":" in x)
        if fields[0] == "@SQ":
            name, length = tags["SN"], int(tags["LN"])
            if name in contigs or length < 1:
                raise GATKError("Duplicate or invalid sequence dictionary entry")
            contigs[name] = length
        if fields[0] == "@RG":
            if not tags.get("SM"):
                raise GATKError("Read group lacks SM")
            names.add(tags["SM"])
        if fields[0] == "@HD":
            sorted_by_coordinate = tags.get("SO") == "coordinate"
    return contigs, names, sorted_by_coordinate


def _check_reference_inputs(config: GATKConfig) -> dict[str, int]:
    contigs, _, _ = _sam_header(Path(config.reference_dict.local_path()).read_text())
    fai = {}
    for line in Path(config.reference_fai.local_path()).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 5 or fields[0] in fai:
            raise GATKError("Malformed or duplicate FAI entry")
        fai[fields[0]] = int(fields[1])
    if not contigs or contigs != fai:
        raise GATKError("FAI and sequence dictionary are incompatible")
    count = 0
    for line in Path(config.intervals.local_path()).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise GATKError("Intervals must be a BED file")
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if chrom not in contigs or not 0 <= start < end <= contigs[chrom]:
            raise GATKError("BED interval incompatible with reference")
        count += 1
    if not count:
        raise GATKError("Empty BED scope")
    return contigs


def _verify_ref_alleles(config: GATKConfig, variants: list[dict[str, Any]]) -> None:
    offsets = {}
    for line in Path(config.reference_fai.local_path()).read_text().splitlines():
        fields = line.split("\t")
        offsets[fields[0]] = tuple(int(x) for x in fields[1:5])
    with Path(config.reference.local_path()).open("rb") as fasta:
        for variant in variants:
            length, offset, bases_per_line, bytes_per_line = offsets[variant["chromosome"]]
            start = variant["position"] - 1
            reference = variant["reference"].upper()
            if not 0 <= start < start + len(reference) <= length:
                raise GATKError("REF allele outside reference sequence")
            if bases_per_line < 1 or bytes_per_line < bases_per_line:
                raise GATKError("Invalid FAI line dimensions")
            observed = bytearray()
            while len(observed) < len(reference):
                position = start + len(observed)
                line_index, column = divmod(position, bases_per_line)
                take = min(len(reference) - len(observed), bases_per_line - column)
                fasta.seek(offset + line_index * bytes_per_line + column)
                chunk = fasta.read(take)
                if len(chunk) != take:
                    raise GATKError("REF allele cannot be read from indexed FASTA")
                observed.extend(chunk)
            if observed.decode("ascii").upper() != reference:
                raise GATKError("VCF REF allele differs from locked FASTA")


def _write_json(path: Path, content: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(content, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def run_gatk(
    config: GATKConfig,
    output_dir: Path,
    *,
    allow_experimental_ont: bool = False,
    runner: Runner = subprocess_runner,
) -> dict[str, Any]:
    """Run a bounded local research branch; require an unused output directory."""
    if allow_experimental_ont is not True:
        raise GATKError("ONT execution is experimental; explicit opt-in is required")
    out = output_dir.resolve()
    try:
        out.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise GATKError(
            "Output directory already exists; stale results must not be reused"
        ) from exc
    calls: list[list[str]] = []

    def call(argv: list[str]) -> str:
        calls.append(argv)
        log = out / f"command-{len(calls):02d}.log"
        try:
            text = runner(argv, out, config.timeout_seconds)
        except subprocess.CalledProcessError as exc:
            log.write_text(str(exc.output or "Command failed"))
            raise
        log.write_text(text)
        return text

    try:
        commands = build_commands(config, out)
        files = _files(config)
        for file in files:
            file.verify()
        contigs = _check_reference_inputs(config)
        _write_json(out / "config.lock.json", asdict(config))
        _write_json(out / "command-plan.json", {"commands": commands, "executed": False})
        version = call([config.gatk, "--version"])
        matches = re.findall(r"\bGATK\)?\s+v?(\d+\.\d+\.\d+\.\d+)(?![\w.+-])", version)
        if matches != [GATK_VERSION]:
            raise GATKError(f"GATK version must be exactly {GATK_VERSION}")
        samtools_version = call([config.samtools, "version"])
        for sample in _samples(config):
            call([config.samtools, "quickcheck", "-v", sample.bam.local_path()])
            header = call([config.samtools, "view", "-H", sample.bam.local_path()])
            observed, names, coordinate_sorted = _sam_header(header)
            if observed != contigs or names != {sample.name} or not coordinate_sorted:
                raise GATKError(
                    "BAM sample, reference dictionary or coordinate-sort header mismatch"
                )
        for argv in commands:
            if "FilterMutectCalls" in argv:
                stats = out / "unfiltered.vcf.gz.stats"
                if not stats.is_file() or not stats.stat().st_size:
                    raise GATKError("Mutect2 stats missing; refusing unfiltered fallback")
            call(argv)
            product = Path(argv[argv.index("-O") + 1])
            if not product.is_file() or not product.stat().st_size:
                raise GATKError(f"Missing output from successful process: {product.name}")
        variants = read_candidates(
            out / "filtered.vcf.gz", config.tumor.name, tuple(s.name for s in _samples(config))
        )
        for variant in variants:
            if variant["chromosome"] not in contigs or (
                variant["position"] + len(variant["reference"]) - 1 > contigs[variant["chromosome"]]
            ):
                raise GATKError("Output variant outside the reference dictionary")
        _verify_ref_alleles(config, variants)
        # Recheck every lock: changing inputs during a run invalidates publication.
        for file in files:
            file.verify()
        warnings = [
            "Research Use Only. Human Review Required. ONT assay is NOT_VALIDATED.",
            "Caller agreement is correlated evidence, not orthogonal confirmation.",
            "BAM reference identity is declared and checksum-bound; header checks compare SN/LN.",
            "Native alleles are not left-normalized; no consensus or clinical report integration.",
        ]
        if not config.germline:
            warnings.append("Germline population resource not supplied.")
        if not config.panel_of_normals:
            warnings.append("Assay-matched ONT panel of normals not supplied.")
        if not config.population_sites:
            warnings.append("Contamination was not estimated; it must not be assumed zero.")
        result = {
            "schema_version": SCHEMA_VERSION,
            "caller": "GATK Mutect2",
            "caller_version": GATK_VERSION,
            "samtools_version": samtools_version.strip(),
            "profile_id": config.profile_id,
            "mode": config.mode,
            "technology": "ONT",
            "status": "COMPLETED" if any(v["technical_pass"] for v in variants) else "NO_CALL",
            "research_only": True,
            "reportable": False,
            "analytical_validation": "NOT_VALIDATED",
            "reference_sha256": config.reference.sha256,
            "config_sha256": digest(out / "config.lock.json"),
            "input_sha256": {f.local_path(): f.sha256 for f in files},
            "artifact_sha256": {p.name: digest(p) for p in out.iterdir() if p.is_file()},
            "commands_executed": calls,
            "warnings": warnings,
            "variants": variants,
        }
        _write_json(out / "evidence.json", result)
        return result
    except (ValueError, OSError, KeyError, TypeError, subprocess.SubprocessError) as exc:
        _write_json(
            out / "failure.json",
            {
                "schema_version": SCHEMA_VERSION,
                "status": "FAILED",
                "reportable": False,
                "reason": str(exc),
                "commands_executed": calls,
            },
        )
        raise GATKError(str(exc)) from exc


def load_config(path: Path) -> GATKConfig:
    """Load a locked JSON manifest; unknown fields are rejected by constructors."""
    raw = json.loads(path.read_text())
    if not isinstance(raw, dict):
        raise GATKError("Config must be a JSON object")
    for key in ("reference", "reference_fai", "reference_dict", "intervals"):
        raw[key] = LockedFile(**raw[key])
    for key in ("tumor", "normal"):
        if key == "tumor" or raw.get(key) is not None:
            item = raw[key]
            if not isinstance(item, dict) or not item:
                raise GATKError(f"{key} must be a populated JSON object or null")
            raw[key] = Sample(
                **{**item, "bam": LockedFile(**item["bam"]), "index": LockedFile(**item["index"])}
            )
    for key in ("germline", "panel_of_normals", "population_sites"):
        if raw.get(key) is not None:
            item = raw[key]
            if not isinstance(item, dict) or not item:
                raise GATKError(f"{key} must be a populated JSON object or null")
            raw[key] = VCFResource(
                **{**item, "vcf": LockedFile(**item["vcf"]), "index": LockedFile(**item["index"])}
            )
    return GATKConfig(**raw)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--allow-experimental-ont", action="store_true")
    args = parser.parse_args()
    try:
        config = load_config(args.config)
        if args.execute:
            result = run_gatk(
                config, args.output_dir, allow_experimental_ont=args.allow_experimental_ont
            )
            print(json.dumps({"status": result["status"], "reportable": False}))
        else:
            print(
                json.dumps(
                    {
                        "status": "NOT_RUN",
                        "execution_enabled": False,
                        "commands": build_commands(config, args.output_dir),
                    },
                    indent=2,
                )
            )
    except (ValueError, OSError, KeyError, TypeError) as exc:
        parser.exit(2, f"GATK adapter blocked/failed: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
