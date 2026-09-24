"""Build-aware MARLIN research adapter, independent of validated legacy bridges."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
import struct
import tempfile
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path

from .execution import CommandRunner, SubprocessRunner, ToolExecutionError
from .marlin_classification import MarlinClassAnnotation, load_marlin_class_annotations
from .marlin_contracts import (
    MarlinClassificationDecision,
    MarlinFeatureSummary,
    MarlinGroupedScore,
    MarlinModelUnitScore,
)
from .marlin_features import marlin_feature_vector_sha256
from .marlin_native_contracts import (
    ADAPTER_VERSION,
    MODEL_SHA256,
    NativeMarlinInstallation,
    NativeMarlinReadiness,
    NativeMarlinReport,
)
from .methylation import (
    _modkit_failed_processing_count,
    count_reads_with_independent_cytosine_mod_groups,
    count_reads_with_modified_base_tags,
    modkit_version,
)
from .models import (
    AlignedBamIntakeReport,
    FileFingerprint,
    GenomeBuild,
    InputKind,
    ModuleRunStatus,
    SampleManifest,
    ToolRecord,
    Verdict,
)
from .modkit_build import PR709_BUILD_ID, identify_modkit_binary
from .reference import sha256_file

_WORKER = Path(__file__).with_name("_marlin_native_worker.py")
_MAX_BYTES = 128 * 1024 * 1024
_MAX_ROWS = 2_000_000
_FEATURE_COUNT = 357340
_APPROVED_RUNTIME_ARCHIVE_SHA256 = (
    "b812927398645d3abaa3a56622bbf1b3279b70ff25067f42040d450083c822c6"
)
# Independently byte-verified restoration/relocation receipts, not installation self-attestation.
# Each exact manifest binds the complete file set, installed bytes and absolute runtime prefix.
_QUALIFIED_RUNTIME_PROFILES = {
    "1f7d1d996b06b851ab4116162bfff86ef9e86ec13ed1cf9e86de7ebfcb09250b": (
        "marlin-native-durable-linux-cpu-v1",
        "3.10.21",
    ),
    "b97f0b2468923415584b756306ed95ec5941d3b26ecb423095a2ac9135b7a461": (
        "marlin-native-source-linux-cpu-v1",
        "3.10.21",
    ),
}
# Relative directory links from the approved conda-pack archive; regular/file-link bytes
# are already covered by the pinned inventory. These links must not bypass enumeration.
_QUALIFIED_DIRECTORY_LINKS = {
    "lib/python3.1": "python3.10",
    "lib/terminfo": "../share/terminfo",
    "lib/icu/current": "78.3",
}
_FEATURES_SHA256 = "9c271460d790d207ce91abaa28c6834d52dde8fe21dc4529033b14c7bb3cc4f5"
_CLASSES_SHA256 = "976ad5fdf347ac4513a39a3e565dee318990b0400c7f35b7cb13e32eeb64afc3"
_PROBE_SHA256 = {
    GenomeBuild.GRCH37: "5fceeda5e76412e2d8ac116d9da320dd8e95645147aad9a4e7be1bef3029db0b",
    GenomeBuild.GRCH38: "9b6d1ce805be0456de71c84b9f2cf8278a1e4d5d510d7c45d136032398664147",
}


def _lines(path: Path) -> Iterator[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    total = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            total += len(raw.encode("utf-8"))
            if number > _MAX_ROWS or total > _MAX_BYTES:
                raise ValueError("MARLIN input exceeds bounded parser size")
            line = raw.strip()
            if line and not line.startswith("#"):
                yield line


def _probe_rows(path: Path) -> Iterator[tuple[str, int, int, str]]:
    seen = set()
    for line in _lines(path):
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError("MARLIN probe map requires BED4")
        chrom, start, end, probe = fields[0], int(fields[1]), int(fields[2]), fields[3]
        if not chrom or not probe or start < 0 or end != start + 1:
            raise ValueError("MARLIN probe map requires single-base positions")
        key = (chrom, start, probe)
        if key in seen:
            raise ValueError("duplicate MARLIN probe mapping")
        seen.add(key)
        yield chrom, start, end, probe


def combined_probe_fractions(path: Path, probe_map: Path) -> dict[str, float]:
    """Pool exact combined 5mC+5hmC call counts over mapped sites and strands."""
    mapping: dict[tuple[str, int], list[str]] = defaultdict(list)
    probe_chromosomes: dict[str, str] = {}
    for chrom, start, _end, probe in _probe_rows(probe_map):
        if probe in probe_chromosomes and probe_chromosomes[probe] != chrom:
            raise ValueError("MARLIN probe mapping spans chromosomes")
        probe_chromosomes[probe] = chrom
        mapping[chrom, start].append(probe)
    counts: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    seen: dict[tuple[str, int], set[str]] = defaultdict(set)
    for line in _lines(path):
        f = line.split()
        if len(f) != 18 or f[3] != "C":
            raise ValueError("MARLIN requires 18-column combined cytosine code C bedMethyl")
        start, end, valid, modified, canonical, other = (int(f[i]) for i in (1, 2, 9, 11, 12, 13))
        extra = [int(f[i]) for i in (14, 15, 16, 17)]
        fraction_percent = float(f[10])
        if (
            min(start, valid, modified, canonical, other, *extra) < 0
            or end != start + 1
            or not f[0]
        ):
            raise ValueError("invalid MARLIN bedMethyl coordinates/counts")
        if valid != modified + canonical + other or other:
            raise ValueError("combined-modification count identity mismatch")
        if (
            not math.isfinite(fraction_percent)
            or not 0 <= fraction_percent <= 100
            or (valid and abs(fraction_percent - 100 * modified / valid) > 0.011)
        ):
            raise ValueError("bedMethyl percent inconsistent with counts")
        strand = f[5]
        if strand not in {"+", "-", "."}:
            raise ValueError("invalid bedMethyl strand")
        key = (f[0], start)
        if strand in seen[key] or (seen[key] and (strand == "." or "." in seen[key])):
            raise ValueError("duplicate or overlapping combined bedMethyl site/strand")
        seen[key].add(strand)
        if key not in mapping:
            raise ValueError("bedMethyl position is outside the locked probe map")
        for probe in mapping[key]:
            counts[probe][0] += modified
            counts[probe][1] += valid
    return {probe: modified / valid for probe, (modified, valid) in counts.items() if valid}


def encode_native_features(
    features: Sequence[str], fractions: Mapping[str, float]
) -> tuple[int, ...]:
    if not features or len(set(features)) != len(features):
        raise ValueError("MARLIN requires nonempty unique ordered features")
    if any(not f or f != f.strip() for f in features):
        raise ValueError("invalid MARLIN feature identifier")
    if any(not math.isfinite(v) or not 0 <= v <= 1 for v in fractions.values()):
        raise ValueError("MARLIN fractions must be finite values from zero to one")
    return tuple(0 if f not in fractions else 1 if fractions[f] >= 0.5 else -1 for f in features)


def selected_marlin_bam_index(manifest: SampleManifest) -> Path:
    """Bind the manifest-selected index and reject ambiguous HTSlib auto-discovery."""
    if manifest.input.kind != InputKind.ALIGNED_BAM or not manifest.input.index_path:
        raise ValueError("MARLIN requires an explicitly selected aligned-BAM index")
    bam = Path(manifest.input.path).absolute()
    selected = Path(manifest.input.index_path).resolve(strict=True)
    if not selected.is_file():
        raise ValueError("MARLIN selected BAM index must be a readable file")
    candidates: set[Path] = set()
    for bam_name in {bam, bam.resolve(strict=True)}:
        for suffix in (".bai", ".csi"):
            candidates.update((Path(str(bam_name) + suffix), bam_name.with_suffix(suffix)))
    available: set[Path] = set()
    for candidate in candidates:
        if candidate.exists() or candidate.is_symlink():
            if not candidate.is_file():
                raise ValueError("MARLIN found an unreadable adjacent BAM index")
            available.add(candidate.resolve(strict=True))
    if selected not in available:
        raise ValueError(
            "MARLIN selected index must use a supported adjacent BAM .bai/.csi filename"
        )
    if available != {selected}:
        raise ValueError(
            "MARLIN found competing adjacent BAM indexes; retain only the selected "
            "BAI/CSI (or links to that same file) before running modkit"
        )
    return selected


def _asset_path(path: str, installation_path: Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else installation_path.parent / p


def _load_installation(path: Path, build: GenomeBuild) -> NativeMarlinInstallation:
    if not path.is_file():
        raise FileNotFoundError(
            "MARLIN installation.json is missing; install the MARLIN resource kit"
        )
    installation = NativeMarlinInstallation.model_validate_json(path.read_text(encoding="utf-8"))
    if installation.runtime_archive_sha256 != _APPROVED_RUNTIME_ARCHIVE_SHA256:
        raise ValueError("MARLIN installation does not pin the approved runtime archive")
    qualified = _QUALIFIED_RUNTIME_PROFILES.get(installation.runtime_manifest.sha256)
    if qualified is None or installation.python_version != qualified[1]:
        raise ValueError(
            "MARLIN requires an independently qualified runtime inventory and prefix; "
            "a rehashed installation manifest cannot qualify a new runtime or relocation"
        )
    if installation.model.sha256 != MODEL_SHA256:
        raise ValueError("MARLIN installation does not pin the approved original model")
    if build not in installation.probe_maps:
        raise ValueError(f"MARLIN installation has no probe map for {build}")
    if installation.probe_maps[build].sha256 != _PROBE_SHA256[build]:
        raise ValueError("MARLIN approved probe map identity differs from selected build")
    if (
        installation.features.sha256 != _FEATURES_SHA256
        or installation.classes.sha256 != _CLASSES_SHA256
    ):
        raise ValueError("MARLIN feature/class annotation identity differs from approved original")
    assets = (
        installation.model,
        installation.features,
        installation.classes,
        installation.probe_maps[build],
        installation.runtime_manifest,
    )
    for asset in assets:
        if not _asset_path(asset.path, path).is_file():
            raise FileNotFoundError(f"MARLIN installed asset is missing: {asset.path}")
    executable = Path(installation.python_executable)
    if executable.name != "python" or executable.parent.name != "bin":
        raise ValueError("MARLIN requires the bin/python executable of a qualified runtime")
    if not executable.is_absolute():
        raise ValueError("MARLIN Python executable must be an absolute path")
    if not Path(installation.python_executable).is_file():
        raise FileNotFoundError("MARLIN isolated Python runtime is missing")
    return installation


def _preflight(installation: NativeMarlinInstallation) -> None:
    with tempfile.TemporaryDirectory(prefix="ontseq-marlin-preflight-") as directory:
        config = Path(directory) / "preflight.json"
        config.write_text(
            json.dumps(
                {"threads": 1, "preflight": True, "python_version": installation.python_version}
            ),
            encoding="utf-8",
        )
        result = SubprocessRunner().run(
            [installation.python_executable, "-I", "-B", str(_WORKER), str(config)],
            timeout_seconds=30,
        )
        if result.returncode or json.loads(result.stdout).get("ready") is not True:
            raise ValueError(
                "MARLIN isolated runtime/network confinement preflight failed: "
                + result.stderr[-500:]
            )


def check_native_marlin_readiness(
    installation_path: Path | None,
    genome_build: GenomeBuild,
) -> NativeMarlinReadiness:
    if installation_path is None:
        return NativeMarlinReadiness(ready=False, reason="MARLIN installation is not configured")
    try:
        _load_installation(installation_path, genome_build)
    except (OSError, ValueError, ToolExecutionError) as exc:
        return NativeMarlinReadiness(ready=False, reason=str(exc))
    return NativeMarlinReadiness(
        ready=True,
        reason=(
            "MARLIN research runtime configured; installed bytes and execution are not yet "
            "verified. Full verification and confinement preflight run before analysis."
        ),
    )


def native_marlin_signature(installation_path: Path, genome_build: GenomeBuild) -> dict[str, str]:
    """Verify actual bytes on every signature request, including resume (no stat-only cache)."""
    installation = _load_installation(installation_path, genome_build)
    signature = {
        "adapter_version": ADAPTER_VERSION,
        "genome_build": genome_build.value,
        "installation_sha256": sha256_file(installation_path),
        "worker_sha256": sha256_file(_WORKER),
        "adapter_sha256": sha256_file(Path(__file__)),
        "contracts_sha256": sha256_file(Path(__file__).with_name("marlin_native_contracts.py")),
        "runtime_archive_sha256": installation.runtime_archive_sha256,
        "python_version": installation.python_version,
        "tensorflow_version": installation.tensorflow_version,
        "keras_version": installation.keras_version,
        "confinement": installation.confinement,
        "qualified_runtime_profile": _QUALIFIED_RUNTIME_PROFILES[
            installation.runtime_manifest.sha256
        ][0],
    }
    for name, asset in (
        ("model", installation.model),
        ("features", installation.features),
        ("classes", installation.classes),
        ("probes", installation.probe_maps[genome_build]),
        ("runtime_manifest", installation.runtime_manifest),
    ):
        digest = sha256_file(_asset_path(asset.path, installation_path))
        if digest != asset.sha256:
            raise ValueError(f"MARLIN {name} checksum differs from installation lock")
        signature[f"{name}_sha256"] = digest
    runtime_files = json.loads(
        _asset_path(installation.runtime_manifest.path, installation_path).read_text(
            encoding="utf-8"
        )
    )
    if (
        not isinstance(runtime_files, dict)
        or not runtime_files
        or installation.python_executable not in runtime_files
    ):
        raise ValueError("MARLIN runtime manifest must include the selected Python executable")
    runtime_root = Path(installation.python_executable).parent.parent
    runtime_entries = tuple(runtime_root.rglob("*"))
    if any(p.is_symlink() and not p.exists() for p in runtime_entries):
        raise ValueError("MARLIN qualified runtime contains a dangling link")
    directory_links = {
        str(p.relative_to(runtime_root)): str(p.readlink())
        for p in runtime_entries
        if p.is_symlink() and p.is_dir()
    }
    if directory_links != _QUALIFIED_DIRECTORY_LINKS:
        raise ValueError("MARLIN runtime directory link map differs from the approved archive")
    inventory = {str(p) for p in runtime_entries if p.is_file()}
    if inventory != set(runtime_files):
        raise ValueError("MARLIN runtime file inventory differs from locked manifest")
    for name, digest in runtime_files.items():
        if (
            not isinstance(name, str)
            or not Path(name).is_absolute()
            or not isinstance(digest, str)
            or not re.fullmatch("[0-9a-f]{64}", digest)
        ):
            raise ValueError("invalid MARLIN runtime file manifest")
        if sha256_file(Path(name)) != digest:
            raise ValueError(f"MARLIN runtime file checksum mismatch: {name}")
    signature["python_sha256"] = runtime_files[installation.python_executable]
    return signature


def _fingerprint(path: Path) -> FileFingerprint:
    return FileFingerprint(size_bytes=path.stat().st_size, sha256=sha256_file(path))


def _verify_selected_index(
    manifest: SampleManifest, index: Path, expected: FileFingerprint
) -> None:
    if selected_marlin_bam_index(manifest) != index or _fingerprint(index) != expected:
        raise ValueError("MARLIN selected BAM index changed during execution")


def _verify_reference_fai(reference_fasta: Path, expected: FileFingerprint) -> None:
    if _fingerprint(Path(str(reference_fasta) + ".fai")) != expected:
        raise ValueError("MARLIN reference FAI changed during execution")


def _write_json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, allow_nan=False, indent=2)


def _group_scores(
    scores: Sequence[MarlinModelUnitScore], annotations: Sequence[MarlinClassAnnotation], field: str
) -> list[MarlinGroupedScore]:
    members: dict[str, list[float]] = defaultdict(list)
    by_id = {a.model_id: a for a in annotations}
    if sorted(by_id) != list(range(1, 43)) or len(annotations) != 42:
        raise ValueError("MARLIN grouping requires exactly annotations 1..42")
    for unit in scores:
        members[getattr(by_id[unit.model_id], field)].append(unit.score)
    return sorted(
        (
            MarlinGroupedScore(label=label, score=math.fsum(values))
            for label, values in members.items()
        ),
        key=lambda s: (-s.score, s.label),
    )


def _worker_scores(
    payload: dict[str, object],
    *,
    tensor_sha256: str,
    installation: NativeMarlinInstallation,
    threads: int,
) -> list[MarlinModelUnitScore]:
    for key, expected in {
        "tensor_sha256": tensor_sha256,
        "model_sha256": MODEL_SHA256,
        "tensorflow": "2.13.1",
        "keras": "2.13.1",
        "python": installation.python_version,
        "threads": threads,
        "confinement": installation.confinement,
    }.items():
        if payload.get(key) != expected:
            raise ValueError(f"MARLIN worker output identity mismatch: {key}")
    values = payload.get("scores")
    if not isinstance(values, list) or len(values) != 42:
        raise ValueError("MARLIN worker must return exactly 42 scores")
    scores = [
        MarlinModelUnitScore(model_id=i, label=f"model_{i}", score=value)
        for i, value in enumerate(values, 1)
    ]
    if not math.isclose(math.fsum(s.score for s in scores), 1, abs_tol=1e-5):
        raise ValueError("MARLIN worker scores do not sum to one")
    digest = hashlib.sha256(struct.pack("<42f", *(s.score for s in scores))).hexdigest()
    if digest != payload.get("score_vector_sha256"):
        raise ValueError("MARLIN worker score checksum mismatch")
    return scores


def run_native_marlin(
    *,
    run_id: str,
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    reference_fasta: Path,
    installation_path: Path,
    output_dir: Path,
    runner: CommandRunner,
    modkit: str,
    samtools: str,
    threads: int,
) -> NativeMarlinReport:
    """Execute one immutable native research run, returning explicit nonfatal failures."""
    base = {
        "run_id": run_id,
        "sample_id": manifest.sample_id,
        "genome_build": manifest.assay.genome_build,
    }
    if output_dir.exists() and any(output_dir.iterdir()):
        return NativeMarlinReport(
            **base,
            status=ModuleRunStatus.FAILED,
            reason="Refusing existing MARLIN outputs; use a fresh stage directory",
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance: dict[str, object] = {}
    tools: list[ToolRecord] = []
    fingerprints: dict[str, FileFingerprint] = {}
    signature: dict[str, str] = {}
    try:
        installation = _load_installation(installation_path, manifest.assay.genome_build)
    except FileNotFoundError as exc:
        report = NativeMarlinReport(**base, status=ModuleRunStatus.NOT_RUN, reason=str(exc))
        _write_json(output_dir / "native-marlin.json", report.model_dump(mode="json"))
        return report
    except (OSError, ValueError) as exc:
        report = NativeMarlinReport(**base, status=ModuleRunStatus.FAILED, reason=str(exc))
        _write_json(output_dir / "native-marlin.json", report.model_dump(mode="json"))
        return report
    try:
        if (
            manifest.input.kind != InputKind.ALIGNED_BAM
            or manifest.sample_id != intake.sample_id
            or manifest.assay.genome_build != intake.genome_build
            or manifest.assay.reference_id != intake.reference_id
            or intake.verdict == Verdict.FAIL
        ):
            raise ValueError("MARLIN manifest/intake identity or aligned-BAM gate mismatch")
        if threads < 1:
            raise ValueError("MARLIN threads must be positive")
        model_threads = min(threads, 8)
        bam = Path(manifest.input.path)
        with reference_fasta.open("rb") as reference_handle:
            compressed = reference_handle.read(2) == b"\x1f\x8b"
        if compressed or reference_fasta.suffix.lower() in {".gz", ".bgz", ".bgzf"}:
            raise ValueError(
                "MARLIN requires uncompressed reference FASTA; "
                "compressed FASTA/GZI is not qualified"
            )
        fingerprints["reference_fai"] = _fingerprint(Path(str(reference_fasta) + ".fai"))
        fingerprints["bam"] = _fingerprint(bam)
        fingerprints["reference"] = _fingerprint(reference_fasta)
        index = selected_marlin_bam_index(manifest)
        fingerprints["bam_index"] = _fingerprint(index)
        if intake.index_fingerprint is None or intake.index_fingerprint.sha256 is None:
            raise ValueError("MARLIN requires an intake BAM-index fingerprint with SHA256")
        if fingerprints["bam_index"] != intake.index_fingerprint:
            raise ValueError("MARLIN BAM index differs from accepted intake fingerprint")
        provenance["bam_index_path"] = str(index)
        if intake.input_fingerprint is not None and (
            intake.input_fingerprint.size_bytes != fingerprints["bam"].size_bytes
            or (
                intake.input_fingerprint.sha256 is not None
                and intake.input_fingerprint.sha256 != fingerprints["bam"].sha256
            )
        ):
            raise ValueError("MARLIN BAM differs from accepted intake fingerprint")
        signature = native_marlin_signature(installation_path, manifest.assay.genome_build)
        _preflight(installation)
        features = tuple(_lines(_asset_path(installation.features.path, installation_path)))
        if len(features) != _FEATURE_COUNT:
            raise ValueError("MARLIN requires exactly 357340 ordered features")
        annotations = load_marlin_class_annotations(
            _asset_path(installation.classes.path, installation_path)
        )
        probes = _asset_path(
            installation.probe_maps[manifest.assay.genome_build].path, installation_path
        )
        include_bed = output_dir / "model-probes.bed"
        seen_positions = set()
        feature_set = set(features)
        with include_bed.open("x", encoding="utf-8") as handle:
            for chrom, start, end, probe in _probe_rows(probes):
                if probe not in feature_set:
                    continue
                key = (chrom, start, end)
                if key not in seen_positions:
                    handle.write(f"{chrom}\t{start}\t{end}\n")
                    seen_positions.add(key)
        if not seen_positions:
            raise ValueError("MARLIN selected build probe map has no model features")
        fingerprints["include_bed"] = _fingerprint(include_bed)
        binary = identify_modkit_binary(modkit)
        modkit = binary.executable
        version_result = runner.run([modkit, "--version"], timeout_seconds=30)
        version = modkit_version(version_result.stdout + "\n" + version_result.stderr)
        if version_result.returncode or version != "0.6.4":
            raise ValueError("MARLIN requires pinned modkit 0.6.4")
        tagged = count_reads_with_modified_base_tags(
            bam, runner=runner, samtools=samtools, threads=threads
        )
        if tagged is None or tagged == 0:
            raise ValueError("MARLIN requires verified MM/ML modified-base tags")
        independent = count_reads_with_independent_cytosine_mod_groups(
            bam,
            runner=runner,
            samtools=samtools,
            threads=threads,
        )
        if independent is None or (independent > 0 and binary.build_id != PR709_BUILD_ID):
            raise ValueError(
                "Independent cytosine MM groups require the exact qualified modkit PR709 binary"
            )
        pileup = output_dir / "combined.bedmethyl"
        log = output_dir / "modkit.log"
        argv = [
            modkit,
            "pileup",
            str(bam),
            str(pileup),
            "--threads",
            str(threads),
            "--modified-bases",
            "5mC",
            "5hmC",
            "--combine-mods",
            "--cpg",
            "--ref",
            str(reference_fasta),
            "--include-bed",
            str(include_bed),
            "--filter-threshold",
            "0.8",
            "--suppress-progress",
            "--log-filepath",
            str(log),
        ]
        provenance.update(
            {
                "pileup_argv": argv,
                "filter_threshold": 0.8,
                "combine_mods": True,
                "combine_strands": False,
                "observed_beta_threshold": 0.5,
                "model_score_threshold": 0.8,
                "assay_assessability": "NOT_ESTABLISHED",
                "weighted_pooling": "sum(N_mod)/sum(N_valid)",
                "threads": threads,
                "model_threads": model_threads,
                "model_interop_threads": 1,
                "independent_cytosine_groups": independent,
            }
        )
        tools.append(
            ToolRecord(
                name="modkit", version=version, parameters={**binary.parameters(), **provenance}
            )
        )
        if identify_modkit_binary(modkit).sha256 != binary.sha256:
            raise ValueError("modkit binary changed during MARLIN preflight")
        _verify_selected_index(manifest, index, fingerprints["bam_index"])
        _verify_reference_fai(reference_fasta, fingerprints["reference_fai"])
        result = runner.run(argv, timeout_seconds=14400)
        _verify_selected_index(manifest, index, fingerprints["bam_index"])
        _verify_reference_fai(reference_fasta, fingerprints["reference_fai"])
        if identify_modkit_binary(modkit).sha256 != binary.sha256:
            raise ValueError("modkit binary changed during MARLIN pileup")
        if result.returncode or _modkit_failed_processing_count(log, result.stdout, result.stderr):
            raise ValueError(
                "MARLIN modkit pileup failed or reported failed records: " + result.stderr[-500:]
            )
        if not pileup.is_file():
            raise ValueError("MARLIN modkit pileup produced no output")
        fingerprints["bedmethyl"] = _fingerprint(pileup)
        fractions = combined_probe_fractions(pileup, probes)
        if _fingerprint(pileup) != fingerprints["bedmethyl"]:
            raise ValueError("MARLIN bedMethyl changed during parsing")
        values = encode_native_features(features, fractions)
        observed = sum(f in fractions for f in features)
        summary = MarlinFeatureSummary(
            observed_model_feature_count=observed,
            explicit_na_feature_count=0,
            absent_feature_count=_FEATURE_COUNT - observed,
            non_model_probe_count=sum(f not in feature_set for f in fractions),
            observed_fraction=observed / _FEATURE_COUNT,
            feature_vector_sha256=marlin_feature_vector_sha256(values),
            feature_artifact_sha256=installation.features.sha256,
        )
        tensor = output_dir / "features.float32"
        with tensor.open("xb") as handle:
            handle.write(struct.pack(f"<{len(values)}f", *values))
        fingerprints["tensor"] = _fingerprint(tensor)
        tensor_sha256 = fingerprints["tensor"].sha256
        assert tensor_sha256 is not None
        shared = {
            **base,
            "feature_summary": summary,
            "tools": tools,
            "parameters": provenance,
            "input_fingerprints": fingerprints,
            "installation_signature": signature,
        }
        if not observed:
            report = NativeMarlinReport(
                **shared,
                status=ModuleRunStatus.NO_CALL,
                reason="No MARLIN model feature was observed; inference was not run",
                decision=MarlinClassificationDecision.UNKNOWN,
            )
        else:
            worker_output = output_dir / "worker-output.json"
            worker_config = output_dir / "worker-config.json"
            _write_json(
                worker_config,
                {
                    "threads": model_threads,
                    "python_version": installation.python_version,
                    "model_path": str(
                        _asset_path(installation.model.path, installation_path).resolve()
                    ),
                    "model_sha256": installation.model.sha256,
                    "tensor_path": str(tensor.resolve()),
                    "tensor_sha256": fingerprints["tensor"].sha256,
                    "output_path": str(worker_output.resolve()),
                },
            )
            worker_argv = [
                installation.python_executable,
                "-I",
                "-B",
                str(_WORKER.resolve()),
                str(worker_config.resolve()),
            ]
            result = runner.run(worker_argv, timeout_seconds=1800)
            (output_dir / "worker.log").write_text(
                result.stdout + "\n" + result.stderr, encoding="utf-8"
            )
            if result.returncode or not worker_output.is_file():
                raise ValueError("MARLIN isolated model worker failed: " + result.stderr[-500:])
            payload = json.loads(worker_output.read_text(encoding="utf-8"))
            scores = _worker_scores(
                payload,
                tensor_sha256=tensor_sha256,
                installation=installation,
                threads=model_threads,
            )
            fingerprints["worker_output"] = _fingerprint(worker_output)
            fingerprints["worker_config"] = _fingerprint(worker_config)
            tools.append(
                ToolRecord(
                    name="MARLIN-native",
                    version=ADAPTER_VERSION,
                    parameters={
                        "worker_argv": worker_argv,
                        **payload,
                        "scores": "see raw_model_scores",
                    },
                )
            )
            classes = _group_scores(scores, annotations, "class_name_current")
            families = _group_scores(scores, annotations, "methylation_class_family")
            lineages = _group_scores(scores, annotations, "lineage")
            report = NativeMarlinReport(
                **{**shared, "tools": tools, "input_fingerprints": fingerprints},
                status=ModuleRunStatus.COMPLETED,
                reason=(
                    "Original MARLIN model executed; assay assessability is not established. "
                    "The 0.8 model-score threshold does not establish specimen confidence; "
                    "the native classification remains UNKNOWN."
                ),
                raw_model_scores=scores,
                class_scores=classes,
                family_scores=families,
                lineage_scores=lineages,
                top_class=classes[0].label,
                top_class_score=classes[0].score,
                decision=MarlinClassificationDecision.UNKNOWN,
                model_score_threshold_met=classes[0].score >= 0.8,
            )
        # Bind the final result to unchanged inputs and installed executable evidence.
        if (
            _fingerprint(bam) != fingerprints["bam"]
            or _fingerprint(reference_fasta) != fingerprints["reference"]
        ):
            raise ValueError("MARLIN input changed during execution")
        if native_marlin_signature(installation_path, manifest.assay.genome_build) != signature:
            raise ValueError("MARLIN installation changed during execution")
        _verify_selected_index(manifest, index, fingerprints["bam_index"])
        _verify_reference_fai(reference_fasta, fingerprints["reference_fai"])
    except (OSError, ValueError, KeyError, TypeError, ToolExecutionError) as exc:
        report = NativeMarlinReport(
            **base,
            status=ModuleRunStatus.FAILED,
            reason=str(exc),
            tools=tools,
            parameters=provenance,
            input_fingerprints=fingerprints,
            installation_signature=signature,
        )
    _write_json(output_dir / "native-marlin.json", report.model_dump(mode="json"))
    return report
