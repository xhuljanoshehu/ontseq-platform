"""Shell-free, fail-closed, opt-in Mutect2 research adapter.

No automatic connection to ONTSeq's canonical runner, CNV/SV registry, or release path.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from .contracts import LockedFile, Mutect2Config, Mutect2Result, SmallVariantCandidate
from .vcf import extract_candidates, sha256_file


@dataclass(frozen=True)
class Step:
    name: str
    argv: tuple[str, ...]
    outputs: tuple[Path, ...] = ()


Executor = Callable[[Step, Path, int], int]


def config_sha256(config: Mutect2Config) -> str:
    encoded = json.dumps(config.model_dump(mode="json"), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode()).hexdigest()


def locked_inputs(config: Mutect2Config) -> dict[str, LockedFile]:
    reference = config.reference
    files = {
        "reference.fasta": reference.fasta,
        "reference.fai": reference.fai,
        "reference.dictionary": reference.dictionary,
        "gatk.jar": config.runtime.gatk_jar,
        "tumor.bam": config.tumor.bam,
        "tumor.index": config.tumor.bai,
    }
    if config.normal:
        files.update({"normal.bam": config.normal.bam, "normal.index": config.normal.bai})
    for label, resource in (
        ("germline", config.germline_resource),
        ("pon", config.panel_of_normals),
        ("contamination_sites", config.contamination_sites),
    ):
        if resource:
            files[label + ".vcf"] = resource.vcf
            files[label + ".index"] = resource.index
    if config.intervals:
        files["intervals"] = config.intervals
    return files


def _nonempty_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Required nonempty local artifact is unavailable: {path.name}")


def _reference_contigs(config: Mutect2Config) -> dict[str, int]:
    fai_contigs = {}
    for line in config.reference.fai.path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 5 or fields[0] in fai_contigs:
            raise ValueError("Invalid or duplicated FAI dictionary entry")
        length = int(fields[1])
        if length < 1:
            raise ValueError("Invalid reference length")
        fai_contigs[fields[0]] = length
    dictionary = {}
    for line in config.reference.dictionary.path.read_text().splitlines():
        if line.startswith("@SQ\t"):
            items = dict(item.split(":", 1) for item in line.split("\t")[1:])
            name, length = items["SN"], int(items["LN"])
            if name in dictionary:
                raise ValueError("Duplicated sequence dictionary entry")
            dictionary[name] = length
    if not fai_contigs or list(fai_contigs.items()) != list(dictionary.items()):
        raise ValueError("Reference FAI and sequence dictionary must agree in order and length")
    return fai_contigs


def preflight(config: Mutect2Config) -> dict[str, str]:
    """Hash every declared local input. Native GATK additionally validates BAMs and indices.

    Reference-binding metadata is an upstream attestation, not proof of correct historical
    alignment or of the absence of tumor/CHIP in a declared normal specimen.
    """
    files = locked_inputs(config)
    fingerprints = {}
    for label, artifact in files.items():
        _nonempty_file(artifact.path)
        digest = sha256_file(artifact.path)
        if digest != artifact.sha256:
            raise ValueError(f"Input checksum mismatch: {label}")
        fingerprints[label] = digest
    fasta = config.reference.fasta.path
    if config.reference.fai.path != Path(str(fasta) + ".fai"):
        raise ValueError("Reference index path must be <fasta>.fai")
    if config.reference.dictionary.path != fasta.with_suffix(".dict"):
        raise ValueError("Sequence dictionary must use the reference basename and .dict")
    if fasta.suffix.lower() not in {".fa", ".fasta", ".fna"}:
        raise ValueError("Only uncompressed FASTA reference files are supported by this adapter")
    _reference_contigs(config)
    samples = (config.tumor,) + ((config.normal,) if config.normal else ())
    for sample in samples:
        if sample.bam.path.suffix != ".bam":
            raise ValueError("This adapter accepts BAM, not SAM, CRAM, POD5 or FASTQ")
        valid = {Path(str(sample.bam.path) + ".bai"), sample.bam.path.with_suffix(".bai")}
        if sample.bai.path not in valid:
            raise ValueError("BAM index path does not match its BAM")
        # Do not let the tool discover a different, unregistered sidecar.
        alternatives = valid | {
            Path(str(sample.bam.path) + ".csi"),
            sample.bam.path.with_suffix(".csi"),
        }
        if any(path != sample.bai.path and path.exists() for path in alternatives):
            raise ValueError("Ambiguous or unregistered BAM index sidecar")
    for resource in (config.germline_resource, config.panel_of_normals, config.contamination_sites):
        if resource:
            if not str(resource.vcf.path).endswith(".vcf.gz"):
                raise ValueError("Resources must be bgzip VCF files with registered .tbi indices")
            if resource.index.path != Path(str(resource.vcf.path) + ".tbi"):
                raise ValueError("VCF index path does not match its VCF")
            if Path(str(resource.vcf.path) + ".csi").exists():
                raise ValueError("Unregistered VCF index sidecar")
    return fingerprints


def build_plan(config: Mutect2Config, workdir: Path) -> tuple[Step, ...]:
    """Plan only; creating a plan never executes GATK or authorizes a lane."""
    workdir = workdir.absolute()
    runtime = config.runtime
    prefix = (
        runtime.java_executable,
        f"-Xmx{runtime.heap_gb}g",
        "-jar",
        str(runtime.gatk_jar.path),
    )
    steps = [
        Step("java_version", (runtime.java_executable, "-version")),
        Step("version", (*prefix, "--version")),
    ]
    reference = str(config.reference.fasta.path)
    for role, sample in (("tumor", config.tumor), ("normal", config.normal)):
        if sample is None:
            continue
        validation = workdir / f"{role}.validation.txt"
        sample_file = workdir / f"{role}.sample.txt"
        steps += [
            Step(
                f"{role}_validate",
                (
                    *prefix,
                    "ValidateSamFile",
                    "--INPUT",
                    str(sample.bam.path),
                    "--REFERENCE_SEQUENCE",
                    reference,
                    "--MODE",
                    "SUMMARY",
                    "--INDEX_VALIDATION_STRINGENCY",
                    "EXHAUSTIVE",
                    "--OUTPUT",
                    str(validation),
                ),
                (validation,),
            ),
            Step(
                f"{role}_sample",
                (*prefix, "GetSampleName", "-I", str(sample.bam.path), "-O", str(sample_file)),
                (sample_file,),
            ),
        ]
    raw = workdir / "unfiltered.vcf.gz"
    stats = Path(str(raw) + ".stats")
    args = [*prefix, "Mutect2", "-R", reference, "-I", str(config.tumor.bam.path)]
    if config.normal:
        args += ["-I", str(config.normal.bam.path), "-normal", config.normal.sample_id]
    args += ["--germline-resource", str(config.germline_resource.vcf.path)]
    if config.panel_of_normals:
        args += ["--panel-of-normals", str(config.panel_of_normals.vcf.path)]
    if config.intervals:
        args += ["-L", str(config.intervals.path)]
    args += ["--native-pair-hmm-threads", str(runtime.pair_hmm_threads), "-O", str(raw)]
    steps.append(Step("mutect2", tuple(args), (raw, Path(str(raw) + ".tbi"), stats)))
    contamination = workdir / "contamination.table"
    segmentation = workdir / "segments.table"
    if config.estimate_contamination:
        if config.contamination_sites is None:
            raise ValueError("Contamination estimation requires a sites resource")
        sites = str(config.contamination_sites.vcf.path)
        for role, sample in (("tumor", config.tumor), ("normal", config.normal)):
            if sample is None:
                continue
            pileup = workdir / f"{role}.pileups.table"
            args = [
                *prefix,
                "GetPileupSummaries",
                "-R",
                reference,
                "-I",
                str(sample.bam.path),
                "-V",
                sites,
                "-L",
                sites,
            ]
            if config.intervals:
                args += ["-L", str(config.intervals.path), "--interval-set-rule", "INTERSECTION"]
            args += ["-O", str(pileup)]
            steps.append(Step(f"{role}_pileup", tuple(args), (pileup,)))
        args = [
            *prefix,
            "CalculateContamination",
            "-I",
            str(workdir / "tumor.pileups.table"),
            "-O",
            str(contamination),
            "--tumor-segmentation",
            str(segmentation),
        ]
        if config.normal:
            args += ["--matched-normal", str(workdir / "normal.pileups.table")]
        steps.append(Step("contamination", tuple(args), (contamination, segmentation)))
    filtered = workdir / "filtered.vcf.gz"
    filter_stats = workdir / "filtering.stats"
    args = [
        *prefix,
        "FilterMutectCalls",
        "-R",
        reference,
        "-V",
        str(raw),
        "--stats",
        str(stats),
        "--filtering-stats",
        str(filter_stats),
        "-O",
        str(filtered),
    ]
    if config.estimate_contamination:
        args += [
            "--contamination-table",
            str(contamination),
            "--tumor-segmentation",
            str(segmentation),
        ]
    steps.append(
        Step("filter", tuple(args), (filtered, Path(str(filtered) + ".tbi"), filter_stats))
    )
    return tuple(steps)


def _private_json(path: Path, data: object) -> None:
    # All files are local and potentially identifying. No external upload is implemented.
    temporary = path.with_name(path.name + ".tmp")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(temporary, path)


def _native_execute(step: Step, workdir: Path, timeout: int) -> int:
    environment = os.environ.copy()
    for key in ("JAVA_TOOL_OPTIONS", "JDK_JAVA_OPTIONS", "_JAVA_OPTIONS", "CLASSPATH"):
        environment.pop(key, None)
    stdout_fd = os.open(
        workdir / f"{step.name}.stdout.log", os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600
    )
    stderr_fd = os.open(
        workdir / f"{step.name}.stderr.log", os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600
    )
    with os.fdopen(stdout_fd, "wb") as stdout, os.fdopen(stderr_fd, "wb") as stderr:
        completed = subprocess.run(
            list(step.argv),
            cwd=workdir,
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            shell=False,
            timeout=timeout,
            check=False,
            umask=0o077 if os.name == "posix" else -1,
        )
    return completed.returncode


def _bounded_text(path: Path, limit: int = 1_048_576) -> str:
    with path.open("r", encoding="utf-8") as handle:
        text = handle.read(limit + 1)
    if len(text) > limit:
        raise ValueError(f"Unexpectedly large control artifact: {path.name}")
    return text


def _check_contamination(path: Path, sample_id: str) -> None:
    lines = [line for line in _bounded_text(path).splitlines() if line and not line.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if len(rows) != 1 or rows[0].get("sample") != sample_id:
        raise ValueError("Contamination estimate is absent or belongs to another sample")
    for column in ("contamination", "error"):
        number = float(rows[0][column])
        if not math.isfinite(number) or number < 0 or number > 1:
            raise ValueError("Contamination estimate or uncertainty is non-finite/out of range")


def run_mutect2(
    config: Mutect2Config,
    output_dir: Path,
    *,
    executor: Executor | None = None,
) -> Mutect2Result:
    """Execute one fresh local lane, or return NOT_RUN when either opt-in is absent.

    Injected executors are always labelled TEST_DOUBLE. They cannot establish real-tool
    qualification. A native run still does not establish analytical/clinical validation.
    """
    base = {
        "run_id": config.run_id,
        "sample_id": config.tumor.sample_id,
        "mode": config.mode,
        "config_sha256": config_sha256(config),
    }
    if not config.enabled or not config.allow_experimental_ont:
        return Mutect2Result(
            **base,
            status="NOT_RUN",
            messages=("Explicit enabled=true and allow_experimental_ont=true are both required.",),
        )
    output_dir = output_dir.absolute()
    # A fresh directory prevents stale successful results being mistaken for this run.
    output_dir.mkdir(mode=0o700, parents=True, exist_ok=False)
    execution_evidence: Literal["NONE", "TEST_DOUBLE", "NATIVE_EXECUTION"] = "NONE"
    execute = _native_execute if executor is None else executor
    inputs: dict[str, str] = {}
    exit_codes: dict[str, int] = {}
    candidates: tuple[SmallVariantCandidate, ...] = ()
    messages = [
        "Experimental ONT branch. Research Use Only; human review required.",
        "Caller agreement is not independent confirmation or biological truth.",
        "Zero extracted candidates does not establish a biological negative.",
        "Alleles retain native VCF representation; left-normalization has not been performed.",
    ]
    if config.mode == "tumor_only":
        messages.append("Tumor-only: germline versus somatic origin is not established.")
    if config.panel_of_normals is None:
        messages.append("No assay-matched ONT panel of normals was supplied.")
    if not config.estimate_contamination:
        messages.append("Cross-sample contamination was not estimated; it is not assumed zero.")
    status: Literal["FAILED", "COMPLETED"] = "FAILED"
    current_step = "preflight"
    try:
        inputs = preflight(config)
        plan = build_plan(config, output_dir)
        _private_json(output_dir / "config.lock.json", config.model_dump(mode="json"))
        _private_json(
            output_dir / "commands.json",
            [
                {"step": s.name, "argv": list(s.argv), "outputs": [str(p) for p in s.outputs]}
                for s in plan
            ],
        )
        for step in plan:
            current_step = step.name
            execution_evidence = "NATIVE_EXECUTION" if executor is None else "TEST_DOUBLE"
            return_code = execute(step, output_dir, config.step_timeout_seconds)
            exit_codes[step.name] = return_code
            if return_code != 0:
                raise RuntimeError(f"Native step exited with code {return_code}")
            for path in step.outputs:
                _nonempty_file(path)
            if step.name in {"version", "java_version"}:
                text = (
                    _bounded_text(output_dir / f"{step.name}.stdout.log")
                    + "\n"
                    + _bounded_text(output_dir / f"{step.name}.stderr.log")
                )
                if step.name == "version":
                    versions = re.findall(
                        r"(?:GATK\)?\s*v?|Version[:\s]+)(\d+\.\d+\.\d+\.\d+)(?![\w.+-])",
                        text,
                        flags=re.IGNORECASE,
                    )
                    if not versions or set(versions) != {config.runtime.expected_version}:
                        raise ValueError("GATK runtime version does not match the pinned baseline")
                else:
                    found = re.search(r'(?:openjdk|java) version "(\d+)', text)
                    if not found or int(found.group(1)) != config.runtime.java_major:
                        raise ValueError("Java runtime major version does not match its lock")
            if step.name.endswith("_sample"):
                sample = config.normal if step.name.startswith("normal") else config.tumor
                if sample is None:
                    raise ValueError("Normal-sample step requires an explicit normal")
                if _bounded_text(step.outputs[0]).strip() != sample.sample_id:
                    raise ValueError("BAM read-group sample name does not match its manifest")
            if step.name == "contamination":
                _check_contamination(output_dir / "contamination.table", config.tumor.sample_id)
        current_step = "input_reverification"
        if preflight(config) != inputs:
            raise ValueError("Input identity changed during execution")
        current_step = "vcf_extraction"
        candidates = extract_candidates(
            output_dir / "filtered.vcf.gz",
            run_id=config.run_id,
            sample_id=config.tumor.sample_id,
            normal_sample_id=config.normal.sample_id if config.normal else None,
        )
        contigs = _reference_contigs(config)
        if any(c.chromosome not in contigs or c.end0 > contigs[c.chromosome] for c in candidates):
            raise ValueError("Variant lies outside the locked reference dictionary")
        status = "COMPLETED"
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as exc:
        candidates = ()
        messages.append(f"{current_step}: {type(exc).__name__}: {exc}")
    outputs = {
        path.name: sha256_file(path)
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path.name != "result.json"
    }
    result = Mutect2Result(
        **base,
        status=status,
        messages=tuple(messages),
        candidates=candidates,
        input_sha256s=inputs,
        output_sha256s=outputs,
        step_exit_codes=exit_codes,
        execution_evidence=execution_evidence,
    )
    _private_json(output_dir / "result.json", result.model_dump(mode="json"))
    return result


class GATKAdapter:
    """Explicit API, deliberately not registered as an ONTSeq production provider."""

    def __init__(self, config: Mutect2Config):
        self.config = config

    def plan(self, output_dir: Path) -> tuple[Step, ...]:
        return build_plan(self.config, output_dir)

    def run(self, output_dir: Path) -> Mutect2Result:
        return run_mutect2(self.config, output_dir)


Mutect2Caller = GATKAdapter
