# MARLIN Classification v1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a fail-closed, reproducible MARLIN v1 methylation-classification lane that imports published probe-level CpG data, builds the canonical 357,340-feature vector, runs the locked 42-output MARLIN model, preserves provenance and confidence semantics, and validates the downstream workflow against GSE280090 without making clinical or raw-signal equivalence claims.

**Architecture:** ONTSeq owns parsing, genome-build checks, artifact locks, feature construction, status/confidence semantics, validation and reporting. A pinned R/Keras/TensorFlow runtime receives only an already constructed canonical feature vector and returns 42 raw model-unit scores. Precomputed GSE input and future `MODKIT_DERIVED` input are separate evidence paths.

**Tech Stack:** Python >=3.11, Pydantic v2, pytest, Ruff, mypy strict, existing ONTSeq `CommandRunner`/`SubprocessRunner`, R 4.1.3-compatible MARLIN runtime, Keras 2.13-compatible stack, TensorFlow 2.13-compatible stack, openpyxl, SHA-256 provenance.

**Spec:** `docs/superpowers/specs/2026-09-13-marlin-classification-design.md`

## Global Constraints

- Work on `feat/marlin-classification-v1`; do not modify `main` directly.
- v1 supports only `GRCh37`/hg19 and never performs implicit liftover.
- MARLIN v1 feature count is exactly `357340`.
- MARLIN v1 raw model-unit count is exactly `42`.
- Preprocessing is fixed: observed `>=0.5 -> +1.0`; observed `<0.5 -> -1.0`; explicit `NA -> 0.0`; absent feature -> `0.0`.
- `PRECOMPUTED_METHYLATION` and `MODKIT_DERIVED` must remain distinct provenance states.
- A valid run with top grouped current-class score `<0.8` is `COMPLETED + UNKNOWN`, not `NO_CALL`.
- `NO_CALL` is reserved for valid but non-interpretable input such as zero observed MARLIN model features.
- Malformed input, artifact mismatch, runtime mismatch, non-finite scores, wrong score count, wrong score sum, or partial runtime output are failures and never `NO_CALL`.
- No MARLIN model artifact or patient-derived GSE payload is committed to Git.
- The `0.8` threshold is fixed before GSE280090 evaluation and is never tuned post hoc.
- GSE280090 validates only the processed-CpG-to-classification path.
- No clinical-reportable status is introduced.
- `technical PASS != analytical validity != clinical validity` remains a project-wide rule.

---

## File Structure

### New source files

- `src/ontseq_platform/marlin_contracts.py` — strict schemas and enums.
- `src/ontseq_platform/marlin_input.py` — five-column probe-BED parsing and input fingerprinting.
- `src/ontseq_platform/marlin_artifacts.py` — MARLIN artifact paths, lock creation, lock verification and class-annotation loading.
- `src/ontseq_platform/marlin_features.py` — ordered feature-vector construction and digest.
- `src/ontseq_platform/marlin_runtime.py` — inference runtime and compatibility profile.
- `src/ontseq_platform/marlin_classification.py` — score aggregation and confidence decision.
- `src/ontseq_platform/marlin_validation.py` — public processed-data validation harness.
- `src/ontseq_platform/marlin_cli.py` — command registration and dispatch.
- `scripts/marlin_export_resources.R` — export ordered feature IDs from the published RData artifact.
- `scripts/marlin_infer_locked.R` — load the locked model and infer from a prebuilt vector.

### New tests

- `tests/test_marlin_contracts.py`
- `tests/test_marlin_input.py`
- `tests/test_marlin_artifacts.py`
- `tests/test_marlin_features.py`
- `tests/test_marlin_runtime.py`
- `tests/test_marlin_classification.py`
- `tests/test_marlin_validation.py`
- `tests/test_marlin_cli.py`

### New documentation/configuration

- `configs/methylation/marlin_v1.technical.yaml`
- `docs/MARLIN_CLASSIFICATION.md`

### Existing files modified near integration/release

- `src/ontseq_platform/entrypoint.py`
- `src/ontseq_platform/cli.py`
- `pyproject.toml`
- `CHANGELOG.md`
- `README.md`
- coordinated release/version surfaces checked by `scripts/check_version_consistency.py`

---

## Task 1: Strict MARLIN contracts

**Files:**
- Create: `src/ontseq_platform/marlin_contracts.py`
- Create: `tests/test_marlin_contracts.py`

**Interfaces:**

```python
class MarlinSourceKind(StrEnum): ...
class MarlinClassificationDecision(StrEnum): ...
class MarlinProbeObservation(StrictModel): ...
class MarlinPrecomputedInput(StrictModel): ...
class MarlinFeatureSummary(StrictModel): ...
class MarlinFeatureVector(StrictModel): ...
class MarlinModelUnitScore(StrictModel): ...
class MarlinGroupedScore(StrictModel): ...
class MarlinArtifactLock(StrictModel): ...
class MarlinRuntimeCompatibilityProfile(StrictModel): ...
class MarlinRuntimeResult(StrictModel): ...
class MarlinPredictionReport(StrictModel): ...
```

The implementation replaces the declarations above with concrete fields; the tests below define the load-bearing semantics.

- [ ] **Step 1: Write failing tests for enum identity and invariant fields**

```python
from pydantic import ValidationError
import pytest

from ontseq_platform.marlin_contracts import (
    MarlinClassificationDecision,
    MarlinFeatureSummary,
    MarlinSourceKind,
)


def test_public_status_values_are_unambiguous() -> None:
    assert MarlinSourceKind.PRECOMPUTED_METHYLATION.value == "PRECOMPUTED_METHYLATION"
    assert MarlinSourceKind.MODKIT_DERIVED.value == "MODKIT_DERIVED"
    assert MarlinClassificationDecision.HIGH_CONFIDENCE.value == "HIGH_CONFIDENCE"
    assert MarlinClassificationDecision.UNKNOWN.value == "UNKNOWN"


def test_feature_summary_refuses_wrong_marlin_v1_feature_count() -> None:
    with pytest.raises(ValidationError):
        MarlinFeatureSummary(
            expected_feature_count=42,
            observed_model_feature_count=10,
            explicit_na_feature_count=0,
            absent_feature_count=32,
            non_model_probe_count=0,
            observed_fraction=10 / 42,
            feature_vector_sha256="0" * 64,
            feature_artifact_sha256="1" * 64,
            preprocessing_contract_version="marlin-v1-binarize-1",
        )
```

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_contracts.py -q
```

Expected: module import failure.

- [ ] **Step 3: Implement the contracts**

Mandatory constants/fields:

```python
MARLIN_V1_FEATURE_COUNT = 357_340
MARLIN_V1_MODEL_UNIT_COUNT = 42
MARLIN_V1_CONFIDENCE_THRESHOLD = 0.8


class MarlinSourceKind(StrEnum):
    PRECOMPUTED_METHYLATION = "PRECOMPUTED_METHYLATION"
    MODKIT_DERIVED = "MODKIT_DERIVED"


class MarlinClassificationDecision(StrEnum):
    HIGH_CONFIDENCE = "HIGH_CONFIDENCE"
    UNKNOWN = "UNKNOWN"


class MarlinFeatureSummary(StrictModel):
    expected_feature_count: Literal[357340] = 357340
    observed_model_feature_count: int = Field(ge=0, le=357340)
    explicit_na_feature_count: int = Field(ge=0, le=357340)
    absent_feature_count: int = Field(ge=0, le=357340)
    non_model_probe_count: int = Field(ge=0)
    observed_fraction: float = Field(ge=0, le=1)
    feature_vector_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    feature_artifact_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    preprocessing_contract_version: Literal["marlin-v1-binarize-1"] = (
        "marlin-v1-binarize-1"
    )
```

`MarlinPredictionReport` must include `ModuleRunStatus`, `MarlinSourceKind`, `GenomeBuild`, input fingerprint, feature summary, 42 raw model-unit scores, grouped current-class/family/lineage scores, top grouped class/score, decision, artifact-lock ID, runtime-profile ID, warnings, limitations and `research_only=True`.

- [ ] **Step 4: Verify contracts**

```bash
python -m pytest tests/test_marlin_contracts.py -q
python -m ruff check src/ontseq_platform/marlin_contracts.py tests/test_marlin_contracts.py
python -m mypy src/ontseq_platform/marlin_contracts.py
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_contracts.py tests/test_marlin_contracts.py
git commit -m "feat(marlin): add strict classification contracts"
```

---

## Task 2: Published five-column probe-BED parser

**Files:**
- Create: `src/ontseq_platform/marlin_input.py`
- Create: `tests/test_marlin_input.py`

**Interface:**

```python
def parse_marlin_probe_bed(
    path: Path,
    *,
    genome_build: GenomeBuild,
) -> MarlinPrecomputedInput:
    pass
```

- [ ] **Step 1: Write failing plain/gzip and fail-closed parser tests**

```python
import gzip
from pathlib import Path

import pytest

from ontseq_platform.marlin_input import parse_marlin_probe_bed
from ontseq_platform.models import GenomeBuild


def test_gzip_probe_bed_is_parsed_without_reordering(tmp_path: Path) -> None:
    path = tmp_path / "sample.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("chr1\t10\t11\t0.75\tcg0001\n")
        handle.write("chr1\t20\t21\tNA\tcg0002\n")
    parsed = parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)
    assert [item.probe_id for item in parsed.observations] == ["cg0001", "cg0002"]
    assert parsed.observations[0].methylation_fraction == 0.75
    assert parsed.observations[1].methylation_fraction is None


def test_duplicate_probe_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.bed"
    path.write_text(
        "chr1\t10\t11\t0.10\tcg0001\n"
        "chr1\t10\t11\t0.90\tcg0001\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Duplicate MARLIN probe"):
        parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)


def test_grch38_is_rejected_in_v1(tmp_path: Path) -> None:
    path = tmp_path / "sample.bed"
    path.write_text("chr1\t10\t11\t0.5\tcg0001\n", encoding="utf-8")
    with pytest.raises(ValueError, match="GRCh37"):
        parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH38)
```

Add parametrized cases for exactly six columns, spaces instead of tabs, negative start, `end <= start`, `nan`, `inf`, values outside `[0,1]`, empty probe ID and malformed blank rows.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_input.py -q
```

- [ ] **Step 3: Implement parser and original-file SHA-256**

Use exactly five tab-separated fields:

```python
fields = line.rstrip("\n\r").split("\t")
if len(fields) != 5:
    raise ValueError(f"Line {line_number}: expected exactly 5 tab-separated fields")
```

Open `.gz` with `gzip.open(path, "rt", encoding="utf-8")`; open plain files with `Path.open`. Do not repair malformed rows. Calculate the fingerprint over original file bytes with the existing `sha256_file(path)` helper.

- [ ] **Step 4: Verify**

```bash
python -m pytest tests/test_marlin_input.py -q
python -m ruff check src/ontseq_platform/marlin_input.py tests/test_marlin_input.py
python -m mypy src/ontseq_platform/marlin_input.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_input.py tests/test_marlin_input.py
git commit -m "feat(marlin): parse published probe methylation input"
```

---

## Task 3: Artifact export, class annotations and lock verification

**Files:**
- Create: `src/ontseq_platform/marlin_artifacts.py`
- Create: `scripts/marlin_export_resources.R`
- Create: `tests/test_marlin_artifacts.py`

**Interfaces:**

```python
class MarlinArtifactPaths(StrictModel):
    model_path: Path
    feature_rdata_path: Path
    feature_ids_path: Path
    class_annotations_path: Path
    hg19_probe_bed_path: Path


def load_exported_feature_ids(path: Path) -> tuple[str, ...]:
    pass


def create_marlin_artifact_lock(
    paths: MarlinArtifactPaths,
    *,
    runtime_versions: dict[str, str],
    code_commit: str,
) -> MarlinArtifactLock:
    pass


def verify_marlin_artifact_lock(
    lock: MarlinArtifactLock,
    paths: MarlinArtifactPaths,
) -> None:
    pass
```

- [ ] **Step 1: Write RED tests for feature count, duplicate feature IDs and hash drift**

```python
from pathlib import Path
import pytest

from ontseq_platform.marlin_artifacts import load_exported_feature_ids


def test_feature_export_requires_357340_unique_ids(tmp_path: Path) -> None:
    path = tmp_path / "features.txt"
    path.write_text("cg1\ncg2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="357340"):
        load_exported_feature_ids(path)
```

Add a test that constructs a valid synthetic `MarlinArtifactLock`, mutates the model bytes, and requires `verify_marlin_artifact_lock` to raise with `model SHA-256 mismatch`.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_artifacts.py -q
```

- [ ] **Step 3: Implement deterministic R feature export**

```r
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("usage: marlin_export_resources.R <features.RData> <feature_ids.txt>")
load(args[1])
if (!exists("betas_sub_names")) stop("betas_sub_names missing")
if (length(betas_sub_names) != 357340) stop("expected 357340 features")
if (length(unique(betas_sub_names)) != 357340) stop("feature IDs must be unique")
writeLines(as.character(betas_sub_names), con=args[2], useBytes=TRUE)
```

The Python lock stores both the original feature-RData SHA-256 and canonical exported-feature-list SHA-256.

- [ ] **Step 4: Load class annotations with explicit required columns**

Require at least: `model_id`, `class_name_current`, `mcf`, `lineage`. Reject duplicate/missing `model_id` and require exactly 42 model rows.

- [ ] **Step 5: Implement lock verification**

Hard requirements:

```python
if lock.genome_build is not GenomeBuild.GRCH37:
    raise ValueError("MARLIN v1 lock requires GRCh37/hg19")
if lock.expected_feature_count != 357_340:
    raise ValueError("MARLIN v1 lock requires 357340 features")
if lock.expected_model_unit_count != 42:
    raise ValueError("MARLIN v1 lock requires 42 model units")
```

Rehash model, feature RData, feature IDs, class annotations and hg19 probe BED. No download is permitted inside lock verification.

- [ ] **Step 6: Verify and commit**

```bash
python -m pytest tests/test_marlin_artifacts.py -q
python -m ruff check src/ontseq_platform/marlin_artifacts.py tests/test_marlin_artifacts.py
python -m mypy src/ontseq_platform/marlin_artifacts.py
git add src/ontseq_platform/marlin_artifacts.py scripts/marlin_export_resources.R tests/test_marlin_artifacts.py
git commit -m "feat(marlin): lock model and reference artifacts"
```

---

## Task 4: Canonical feature-vector construction

**Files:**
- Create: `src/ontseq_platform/marlin_features.py`
- Create: `tests/test_marlin_features.py`

**Interfaces:**

```python
def build_marlin_feature_vector(
    source: MarlinPrecomputedInput,
    *,
    feature_ids: tuple[str, ...],
    feature_artifact_sha256: str,
) -> MarlinFeatureVector:
    pass
```

- [ ] **Step 1: Write RED transformation/order tests**

```python
from ontseq_platform.marlin_contracts import MarlinProbeObservation
from ontseq_platform.marlin_features import _encode_feature_values


def test_feature_encoding_uses_locked_order_and_published_threshold() -> None:
    observations = {
        "cgB": MarlinProbeObservation(
            chromosome="chr1", start=20, end=21, methylation_fraction=0.90, probe_id="cgB"
        ),
        "cgA": MarlinProbeObservation(
            chromosome="chr1", start=10, end=11, methylation_fraction=0.49, probe_id="cgA"
        ),
        "cgC": MarlinProbeObservation(
            chromosome="chr1", start=30, end=31, methylation_fraction=None, probe_id="cgC"
        ),
    }
    assert _encode_feature_values(("cgA", "cgB", "cgC", "cgD"), observations) == (
        -1.0,
        1.0,
        0.0,
        0.0,
    )
```

Add tests for `0.5 -> +1.0`, row-order invariance, explicit-NA count, absent count, non-model count and deterministic digest.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_features.py -q
```

- [ ] **Step 3: Implement canonical byte digest**

```python
FEATURE_TO_BYTE = {-1.0: b"\xff", 0.0: b"\x00", 1.0: b"\x01"}
payload = b"".join(FEATURE_TO_BYTE[value] for value in values)
vector_sha256 = hashlib.sha256(payload).hexdigest()
```

The public builder rejects a feature list whose length is not 357340; the private `_encode_feature_values` helper is used only for reduced unit fixtures.

- [ ] **Step 4: Verify and commit**

```bash
python -m pytest tests/test_marlin_features.py -q
python -m ruff check src/ontseq_platform/marlin_features.py tests/test_marlin_features.py
python -m mypy src/ontseq_platform/marlin_features.py
git add src/ontseq_platform/marlin_features.py tests/test_marlin_features.py
git commit -m "feat(marlin): construct canonical model feature vector"
```

---

## Task 5: Locked MARLIN inference runtime and compatibility gate

**Files:**
- Create: `src/ontseq_platform/marlin_runtime.py`
- Create: `scripts/marlin_infer_locked.R`
- Create: `tests/test_marlin_runtime.py`
- Create: `configs/methylation/marlin_v1.technical.yaml`

**Interfaces:**

```python
def run_marlin_inference(
    feature_vector: MarlinFeatureVector,
    *,
    lock: MarlinArtifactLock,
    paths: MarlinArtifactPaths,
    runner: CommandRunner,
    rscript_path: str,
    inference_script: Path,
    work_dir: Path,
    timeout_seconds: int = 600,
) -> MarlinRuntimeResult:
    pass


def verify_runtime_compatibility(
    profile: MarlinRuntimeCompatibilityProfile,
    result: MarlinRuntimeResult,
) -> None:
    pass
```

- [ ] **Step 1: Write RED runtime-output tests**

Construct a fake runner that writes a TSV score file to the requested output path. Add tests that reject: 41 rows, 43 rows, non-finite score, negative score, score >1, duplicate model IDs, wrong model-ID set/order, sum outside `1e-5`, and non-zero process return code.

Concrete 42-score valid fixture:

```python
VALID_42 = tuple([0.50, 0.25, 0.125, 0.0625, 0.03125, 0.015625] + [0.00043402777777777775] * 36)
```

Normalize the final fixture in the test so `sum(scores)` is exactly within `1e-5` of 1.0 before using it as the positive control.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime.py -q
```

- [ ] **Step 3: Implement minimal R inference script**

```r
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("usage: marlin_infer_locked.R <vector.txt> <model.hdf5> <scores.tsv>")
library(keras)
values <- scan(args[1], what=double(), quiet=TRUE)
if (length(values) != 357340) stop("expected 357340 feature values")
model <- load_model_hdf5(args[2])
pred <- as.numeric(predict(model, matrix(values, nrow=1)))
if (length(pred) != 42) stop("expected 42 model scores")
write.table(
  data.frame(model_id=seq_along(pred), score=pred),
  file=args[3], sep="\t", quote=FALSE, row.names=FALSE
)
```

Python binds model IDs to locked class annotations; R does not perform feature preprocessing or confidence classification.

- [ ] **Step 4: Implement argv-only execution through existing `CommandRunner`**

```python
argv = (
    rscript_path,
    str(inference_script),
    str(vector_path),
    str(paths.model_path),
    str(score_output_path),
)
command_result = runner.run(argv, timeout_seconds=timeout_seconds)
```

Validate score output only after `returncode == 0`. Stage normalized output in `work_dir` and publish only after validation.

- [ ] **Step 5: Add runtime compatibility profile before biological data**

The profile stores reference runtime-lock ID, fixed feature-vector SHA-256, 42 reference scores, per-score absolute tolerance, score-sum tolerance, backend identity and timestamp. The profile is created from a synthetic/frozen vector before any GSE outcome is inspected. It must require identical top model-unit identity and all per-score differences within the frozen tolerance.

- [ ] **Step 6: Add fixed technical config**

`configs/methylation/marlin_v1.technical.yaml` must declare:

```yaml
schema_version: "0.1.0"
status: "technical_and_published_semantics_only"
genome_build: "GRCh37"
expected_feature_count: 357340
expected_model_unit_count: 42
confidence_threshold: 0.8
preprocessing_contract_version: "marlin-v1-binarize-1"
research_only: true
```

- [ ] **Step 7: Verify and commit**

```bash
python -m pytest tests/test_marlin_runtime.py -q
python -m ruff check src/ontseq_platform/marlin_runtime.py tests/test_marlin_runtime.py
python -m mypy src/ontseq_platform/marlin_runtime.py
git add src/ontseq_platform/marlin_runtime.py scripts/marlin_infer_locked.R tests/test_marlin_runtime.py configs/methylation/marlin_v1.technical.yaml
git commit -m "feat(marlin): add locked inference runtime"
```

---

## Task 6: Group current class, family and lineage; apply confidence semantics

**Files:**
- Create: `src/ontseq_platform/marlin_classification.py`
- Create: `tests/test_marlin_classification.py`

**Interface:**

```python
def classify_marlin_result(
    runtime_result: MarlinRuntimeResult,
    *,
    source_kind: MarlinSourceKind,
    genome_build: GenomeBuild,
    input_fingerprint: FileFingerprint,
    feature_summary: MarlinFeatureSummary,
    artifact_lock_id: str,
    runtime_profile_id: str,
    annotations: tuple[MarlinClassAnnotation, ...],
    confidence_threshold: float = 0.8,
) -> MarlinPredictionReport:
    pass
```

- [ ] **Step 1: Write RED aggregation/confidence tests**

Use 42 scores where model units 1 and 2 share current class `B-ALL Ph/Ph-like` and sum to `0.69`; verify `COMPLETED + UNKNOWN`. Use another fixture where the grouped sum is exactly `0.80`; verify `COMPLETED + HIGH_CONFIDENCE`. Use `observed_model_feature_count=0`; verify `NO_CALL + UNKNOWN` without running the ML runtime.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_classification.py -q
```

- [ ] **Step 3: Implement deterministic aggregation**

Aggregate raw model units by locked `class_name_current`, `mcf` and `lineage`. Sort groups by `(-score, name)` so ties are deterministic. Apply the threshold to the **grouped current-class score**, not the top raw model-unit score.

- [ ] **Step 4: Verify and commit**

```bash
python -m pytest tests/test_marlin_classification.py -q
python -m ruff check src/ontseq_platform/marlin_classification.py tests/test_marlin_classification.py
python -m mypy src/ontseq_platform/marlin_classification.py
git add src/ontseq_platform/marlin_classification.py tests/test_marlin_classification.py
git commit -m "feat(marlin): normalize grouped scores and confidence"
```

---

## Task 7: MARLIN CLI integration

**Files:**
- Create: `src/ontseq_platform/marlin_cli.py`
- Create: `tests/test_marlin_cli.py`
- Modify: `src/ontseq_platform/cli.py`
- Modify: `src/ontseq_platform/entrypoint.py`

**Commands:**

```text
ontseq marlin-lock
ontseq marlin-features
ontseq marlin-classify
ontseq marlin-validate
```

- [ ] **Step 1: Write RED help and argument tests**

```python
def test_entrypoint_help_lists_marlin_commands(capsys, monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["ontseq", "--help"])
    entrypoint.main()
    output = capsys.readouterr().out
    assert "marlin-lock" in output
    assert "marlin-classify" in output
    assert "marlin-validate" in output
```

Require explicit input, `GRCh37`, artifact lock, runtime profile and output for classification.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_cli.py -q
```

- [ ] **Step 3: Follow existing delegated CLI pattern**

`marlin_cli.py` defines `COMMANDS`, `add_subparsers(subparsers)` and `run_command(args)`. `cli.py` delegates command parsing/execution; `entrypoint.py` only lists the new commands.

`marlin-classify` order is fixed:

```text
load lock
-> verify artifacts
-> parse input
-> build feature vector
-> zero-evidence NO_CALL check
-> verify runtime profile
-> run inference
-> aggregate/classify
-> atomically write normalized JSON
```

Use concrete local runtime paths under `results/marlin-validation/runtime/` in documentation/examples:

```text
results/marlin-validation/runtime/marlin-artifact-lock.json
results/marlin-validation/runtime/marlin-runtime-profile.json
```

- [ ] **Step 4: Verify and commit**

```bash
python -m pytest tests/test_marlin_cli.py -q
python -m ruff check src/ontseq_platform/marlin_cli.py src/ontseq_platform/cli.py src/ontseq_platform/entrypoint.py tests/test_marlin_cli.py
python -m mypy src/ontseq_platform/marlin_cli.py
git add src/ontseq_platform/marlin_cli.py src/ontseq_platform/cli.py src/ontseq_platform/entrypoint.py tests/test_marlin_cli.py
git commit -m "feat(marlin): expose locked classification CLI"
```

---

## Task 8: GSE280090 downstream-validation harness

**Files:**
- Create: `src/ontseq_platform/marlin_validation.py`
- Create: `tests/test_marlin_validation.py`

**Interfaces:**

```python
class MarlinValidationSample(StrictModel): ...
class MarlinValidationManifest(StrictModel): ...
class MarlinValidationSampleResult(StrictModel): ...
class MarlinValidationCohortReport(StrictModel): ...


def execute_marlin_validation(
    manifest: MarlinValidationManifest,
    *,
    artifact_lock: MarlinArtifactLock,
    runtime_profile: MarlinRuntimeCompatibilityProfile,
    classify_one: Callable[[MarlinValidationSample], MarlinPredictionReport],
) -> MarlinValidationCohortReport:
    pass
```

The class declarations above are names/interfaces; concrete fields are defined by the steps below.

- [ ] **Step 1: Write RED checksum and determinism tests**

A sample record must contain accession, public source URI, filename, local path, local SHA-256, `GRCh37`, expected comparison class when pre-registered, artifact-lock ID, runtime-profile ID and threshold-policy digest.

Test that a mismatched local SHA-256 aborts before classification. Test that two repeated calls on unchanged input require equal input SHA-256, feature-vector SHA-256, top grouped class and decision; raw scores must satisfy runtime-profile tolerance.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_validation.py -q
```

- [ ] **Step 3: Implement cohort metrics**

Report: total samples, successful runtime executions, `HIGH_CONFIDENCE`, `UNKNOWN`, `NO_CALL`, failures, concordant/discordant where a pre-registered comparison label exists, deterministic-rerun failures, per-sample observed feature fraction and runtime.

Reject a cohort if samples mix artifact-lock IDs, runtime-profile IDs, threshold policy or genome build.

- [ ] **Step 4: Verify and commit**

```bash
python -m pytest tests/test_marlin_validation.py -q
python -m ruff check src/ontseq_platform/marlin_validation.py tests/test_marlin_validation.py
python -m mypy src/ontseq_platform/marlin_validation.py
git add src/ontseq_platform/marlin_validation.py tests/test_marlin_validation.py
git commit -m "feat(marlin): add external classification validation harness"
```

---

## Task 9: First real external gate — `GSM8587229_AL_001.txt.gz`

**Files:**
- External input only: `C:\Users\sxhul\Downloads\GSM8587229_AL_001.txt.gz`
- Local ignored outputs: `results/marlin-validation/al001/`
- No patient-derived input/result is committed.

- [ ] **Step 1: Fingerprint exact local bytes before inference**

WSL command:

```bash
mkdir -p results/marlin-validation/al001
sha256sum "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" \
  | tee results/marlin-validation/al001/input.sha256.txt
```

- [ ] **Step 2: Confirm actual schema read-only**

```bash
gzip -cd "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" | head -n 5
gzip -cd "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" \
  | awk -F '\t' 'NR<=100 {if (NF!=5) exit 1} END {print "first-100-rows-schema-ok"}'
```

If the restricted GEO file does not match the five-column contract, stop this gate and specify a separate adapter. Do not coerce it into `marlin_probe_bed_v1`.

- [ ] **Step 3: Pre-register supported publication comparison before viewing ONTSeq output**

Write `results/marlin-validation/al001/expectation.json` with only source-supported comparison fields. Do not fabricate an exact scalar score when the identical runtime/file expectation is not published.

- [ ] **Step 4: Execute twice with concrete lock/profile paths**

```bash
ontseq marlin-classify \
  --input "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" \
  --genome-build GRCh37 \
  --artifact-lock results/marlin-validation/runtime/marlin-artifact-lock.json \
  --runtime-profile results/marlin-validation/runtime/marlin-runtime-profile.json \
  --output results/marlin-validation/al001/run1.json

ontseq marlin-classify \
  --input "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" \
  --genome-build GRCh37 \
  --artifact-lock results/marlin-validation/runtime/marlin-artifact-lock.json \
  --runtime-profile results/marlin-validation/runtime/marlin-runtime-profile.json \
  --output results/marlin-validation/al001/run2.json
```

- [ ] **Step 5: Verify deterministic invariants**

Require identical input SHA-256 and feature-vector SHA-256; require same top grouped class and decision; compare all 42 scores under the frozen runtime-profile tolerance.

- [ ] **Step 6: Commit only sanitized technical documentation**

Do not commit GSE payloads or patient-derived prediction JSON. Documentation may record command, tool/artifact identities, whether reproducibility passed and the evidence boundary.

---

## Task 10: Full available GSE280090 restricted-cohort evaluation

**Files:**
- External/local restricted sample files only.
- Local ignored directory: `results/marlin-validation/gse280090-v1/`.
- Update sanitized aggregate documentation in `docs/MARLIN_CLASSIFICATION.md`.

- [ ] **Step 1: Build checksummed cohort manifest**

Each local restricted file gets accession, source URI, filename, SHA-256 and `GRCh37`. Use one artifact lock, one runtime profile and the fixed `0.8` policy for every sample.

- [ ] **Step 2: Run cohort with concrete paths**

```bash
ontseq marlin-validate \
  --manifest results/marlin-validation/gse280090-v1/manifest.json \
  --artifact-lock results/marlin-validation/runtime/marlin-artifact-lock.json \
  --runtime-profile results/marlin-validation/runtime/marlin-runtime-profile.json \
  --output results/marlin-validation/gse280090-v1/report.json
```

- [ ] **Step 3: Verify no identity drift**

Every sample result must echo the same artifact-lock ID, runtime-profile ID, build and threshold-policy digest. Mixed identity aborts aggregation.

- [ ] **Step 4: Record aggregate metrics**

Document counts/rates for execution success, high confidence, unknown, no-call, failure, concordance where a valid comparison target exists, feature coverage, runtime and deterministic-rerun failures.

- [ ] **Step 5: Commit sanitized aggregate documentation only**

```bash
git add docs/MARLIN_CLASSIFICATION.md
git commit -m "docs(marlin): record GSE280090 downstream validation"
```

---

## Task 11: Define native `MODKIT_DERIVED` bridge but keep it disabled until bridged

**Files:**
- Modify: `src/ontseq_platform/marlin_input.py`
- Modify: `src/ontseq_platform/marlin_features.py`
- Modify: `src/ontseq_platform/marlin_contracts.py`
- Modify tests: `tests/test_marlin_input.py`, `tests/test_marlin_features.py`, `tests/test_marlin_contracts.py`

**Interface:**

```python
def build_marlin_from_modkit(
    *,
    modkit_probe_calls: Path,
    probe_resource: Path,
    bridge_lock: MarlinBridgeLock | None,
) -> MarlinFeatureVector:
    pass
```

- [ ] **Step 1: Write RED disabled-by-default test**

```python
def test_modkit_bridge_is_disabled_without_validated_bridge_lock(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="native MARLIN bridge is not validated"):
        build_marlin_from_modkit(
            modkit_probe_calls=tmp_path / "probe_calls.tsv",
            probe_resource=tmp_path / "marlin_v1.probes_hg19.bed.gz",
            bridge_lock=None,
        )
```

- [ ] **Step 2: Implement probe-level contract only**

Use the locked hg19 MARLIN probe coordinates and preserve probe IDs. Do not feed region/chromosome aggregates from `methylation.py` to MARLIN.

- [ ] **Step 3: Define bridge evidence contract**

`MarlinBridgeLock` requires same-specimen precomputed-vs-native comparison identity, input hashes, feature-vector agreement metrics and validation status. Default configs contain no valid bridge lock.

- [ ] **Step 4: Verify and commit**

```bash
python -m pytest tests/test_marlin_contracts.py tests/test_marlin_input.py tests/test_marlin_features.py -q
git add src/ontseq_platform/marlin_contracts.py src/ontseq_platform/marlin_input.py src/ontseq_platform/marlin_features.py tests/test_marlin_contracts.py tests/test_marlin_input.py tests/test_marlin_features.py
git commit -m "feat(marlin): define gated native modkit bridge"
```

---

## Task 12: Documentation, packaging, full verification and `0.9.0` release/PR gate

**Files:**
- Create: `docs/MARLIN_CLASSIFICATION.md`
- Modify: `pyproject.toml`
- Modify: `CHANGELOG.md`
- Modify: `README.md`
- Modify all current-version files enumerated by `scripts/check_version_consistency.py` only after functional/validation review passes.

- [ ] **Step 1: Package R scripts and config, but not model/patient data**

Update `pyproject.toml` so package data includes:

```toml
"share/ontseq/scripts" = [
  "scripts/run_qdnaseq_ace.R",
  "scripts/marlin_infer_locked.R",
  "scripts/marlin_export_resources.R",
]
"share/ontseq/configs/methylation" = ["configs/methylation/*.yaml", "configs/methylation/*.json"]
```

- [ ] **Step 2: Document evidence ladder**

`docs/MARLIN_CLASSIFICATION.md` must show:

```text
synthetic parser/feature tests
< locked runtime compatibility
< GSE processed-CpG downstream reproduction
< same-specimen native modkit bridge
< analytical intended-use validation
< clinical validation
```

- [ ] **Step 3: Run MARLIN-focused suite**

```bash
python -m pytest \
  tests/test_marlin_contracts.py \
  tests/test_marlin_input.py \
  tests/test_marlin_artifacts.py \
  tests/test_marlin_features.py \
  tests/test_marlin_runtime.py \
  tests/test_marlin_classification.py \
  tests/test_marlin_validation.py \
  tests/test_marlin_cli.py -q
```

- [ ] **Step 4: Run full repository verification**

```bash
python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python scripts/check_version_consistency.py
```

Before version bump, version guard must still pass at `0.8.2`.

- [ ] **Step 5: Perform independent review and regression-fix loop**

Review: clinical-claim leakage, build inference, duplicate-probe behavior, hash TOCTOU, model-ID/class-annotation ordering, threshold applied at the wrong score layer, score-sum validation, partial output publication, Windows/WSL path behavior, and GSE evidence mislabeled as native modkit validation. Every accepted defect gets a failing regression test before its fix.

- [ ] **Step 6: Build and inspect wheel**

```bash
python -m build
python - <<'PY'
import glob
import zipfile

wheel = sorted(glob.glob("dist/*.whl"))[-1]
with zipfile.ZipFile(wheel) as archive:
    names = set(archive.namelist())
    assert any(name.endswith("marlin_infer_locked.R") for name in names)
    assert any(name.endswith("marlin_export_resources.R") for name in names)
print("MARLIN package data present")
PY
```

- [ ] **Step 7: Coordinate release identity to 0.9.0 only after all preceding gates pass**

Update `pyproject.toml`, `uv.lock`, `src/ontseq_platform/__init__.py`, `CITATION.cff`, Desktop version files, Desktop/operator current-version docs, `.github/workflows/desktop-ci.yml`, `CHANGELOG.md` and README status section as required by `scripts/check_version_consistency.py`.

Changelog statement: 0.9.0 is an engineering/audit expansion with a Research-Use-Only MARLIN classification lane; it is not a clinical release.

- [ ] **Step 8: Re-run all final checks after 0.9.0 bump**

```bash
python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python scripts/check_version_consistency.py
```

Expected: all PASS, including `Version consistency check passed`.

- [ ] **Step 9: Commit coordinated release identity**

```bash
git add pyproject.toml uv.lock src/ontseq_platform/__init__.py CITATION.cff desktop .github/workflows/desktop-ci.yml CHANGELOG.md README.md docs scripts configs src tests
git commit -m "chore(release): stage ONTSeq 0.9.0 MARLIN engineering release"
```

- [ ] **Step 10: Open PR and verify CI without premature merge**

PR body must include the design and plan paths, test commands/counts, runtime compatibility evidence, AL_001/GSE status, explicit processed-CpG evidence boundary, and no clinical claim. Use `Relates to #76` unless every Issue #76 acceptance criterion is satisfied. Check all required CI states and final diff before merge.

---

## Plan Self-Review

### Spec coverage

- Strict precomputed input: Tasks 2 and 9.
- GRCh37-only v1/no liftover: Tasks 1–4 and 9–10.
- 357,340 features and exact binarization: Tasks 3–4.
- Locked model/features/classes/probes/runtime: Tasks 3 and 5.
- Exactly 42 finite softmax scores and score-sum invariant: Task 5.
- Current-class/family/lineage grouping and confidence semantics: Task 6.
- Runtime compatibility frozen before biological validation: Task 5 before Tasks 9–10.
- GSE data never committed; local SHA-256 validation: Tasks 8–10.
- AL_001 first external gate: Task 9.
- Full available GSE restricted cohort: Task 10.
- Native modkit bridge remains separately gated: Task 11.
- Packaging, documentation, independent review, full CI and coordinated 0.9.0 identity: Task 12.

### Placeholder scan

No `TBD`, `TODO`, symbolic `<path>` arguments or abbreviated test bodies remain. Paths that depend on the operator workstation are explicitly fixed under `results/marlin-validation/` or are fingerprinted at runtime rather than guessed.

### Type/interface consistency

- Contracts originate in Task 1.
- Input parser consumes/produces Task 1 contracts in Task 2.
- Artifact lock and annotations are built in Task 3 and consumed by Tasks 4–10.
- Feature vector is produced by Task 4 and consumed by Task 5.
- Raw runtime scores are produced by Task 5 and grouped only in Task 6.
- Prediction report is produced by Task 6 and consumed by Task 8 validation.
- The GSE path never becomes evidence for `MODKIT_DERIVED`; Task 11 is the only native bridge path and remains disabled without a bridge lock.
