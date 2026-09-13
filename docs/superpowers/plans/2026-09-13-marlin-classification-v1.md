# MARLIN Classification v1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a fail-closed, reproducible MARLIN v1 methylation-classification lane that imports published probe-level CpG data, constructs the canonical 357,340-feature vector, runs the locked 42-class MARLIN model, preserves provenance/confidence semantics, and validates the downstream workflow against GSE280090 without making clinical or raw-signal equivalence claims.

**Architecture:** ONTSeq owns input parsing, build checks, artifact locks, feature construction, status/confidence semantics, validation and reporting. A pinned R/Keras/TensorFlow runtime receives only an already-built canonical feature vector and returns exactly 42 model-unit scores. Precomputed GSE inputs remain a separate evidence path from future `MODKIT_DERIVED` inputs.

**Tech Stack:** Python 3.11+, Pydantic v2, pytest, Ruff, mypy strict, existing ONTSeq `CommandRunner`/`SubprocessRunner`, R 4.1.3-compatible MARLIN runtime, Keras 2.13-compatible stack, TensorFlow 2.13-compatible stack, openpyxl for class annotations, SHA-256 provenance.

**Spec:** `docs/superpowers/specs/2026-09-13-marlin-classification-design.md`

## Global Constraints

- Start from `feat/marlin-classification-v1`; do not modify `main` directly.
- First supported build is `GRCh37`/hg19 only; no implicit liftover.
- MARLIN v1 expected feature count is exactly `357340`.
- MARLIN v1 expected raw model-unit count is exactly `42`.
- Preprocessing is exactly: observed `>=0.5 -> +1.0`, observed `<0.5 -> -1.0`, explicit `NA -> 0.0`, absent feature -> `0.0`.
- `PRECOMPUTED_METHYLATION` and `MODKIT_DERIVED` provenance must never be conflated.
- A valid run with top grouped-class score `<0.8` is `COMPLETED + UNKNOWN`, not `NO_CALL`.
- `NO_CALL` is reserved for valid but non-interpretable input, including zero observed model features.
- Malformed input, lock mismatch, runtime mismatch, non-finite scores, wrong class count, wrong score sum, or partial runtime output are `FAILED`/exceptions, never `NO_CALL`.
- No MARLIN artifact or patient-derived GSE payload is committed to Git.
- No post-hoc threshold tuning against GSE280090.
- GSE280090 validates only the downstream processed-CpG-to-classification path.
- No clinical-reportable status is introduced.
- Existing `technical PASS != analytical validity != clinical validity` project rule remains unchanged.

---

## File Structure Map

### New source files

- `src/ontseq_platform/marlin_contracts.py` — enums and strict Pydantic contracts only.
- `src/ontseq_platform/marlin_input.py` — strict five-column probe-BED parser, gzip support, input fingerprinting.
- `src/ontseq_platform/marlin_artifacts.py` — resource export/inspection, artifact lock creation and verification.
- `src/ontseq_platform/marlin_features.py` — ordered 357,340-feature construction and vector hashing.
- `src/ontseq_platform/marlin_runtime.py` — runtime preflight, score execution/parsing, compatibility profile.
- `src/ontseq_platform/marlin_classification.py` — grouped class/family/lineage aggregation and `HIGH_CONFIDENCE`/`UNKNOWN` semantics.
- `src/ontseq_platform/marlin_validation.py` — GSE validation manifests, deterministic reruns, cohort metrics.
- `src/ontseq_platform/marlin_cli.py` — CLI registration and command handlers.
- `scripts/marlin_infer_locked.R` — minimal runtime script that consumes an ONTSeq-built feature vector and emits 42 scores; no feature preprocessing.
- `scripts/marlin_export_resources.R` — deterministic one-time export of ordered feature IDs from `marlin_v1.features.RData` into canonical text for hashing/use by Python.

### New tests

- `tests/test_marlin_contracts.py`
- `tests/test_marlin_input.py`
- `tests/test_marlin_artifacts.py`
- `tests/test_marlin_features.py`
- `tests/test_marlin_runtime.py`
- `tests/test_marlin_classification.py`
- `tests/test_marlin_validation.py`
- `tests/test_marlin_cli.py`

### New docs/configs

- `docs/MARLIN_CLASSIFICATION.md`
- `configs/methylation/marlin_v1.technical.yaml`
- optional generated local validation outputs remain under ignored `results/` paths.

### Existing files modified late in the plan

- `src/ontseq_platform/entrypoint.py` — discover MARLIN commands.
- `src/ontseq_platform/cli.py` — delegate MARLIN command family or import the MARLIN CLI registrar.
- `pyproject.toml` — package the MARLIN R scripts/config and eventually bump release version.
- `CHANGELOG.md`, `README.md`, `CITATION.cff`, `uv.lock`, Desktop version surfaces — coordinated `0.9.0` release step only after code/validation gates pass.

---

### Task 1: Define MARLIN contracts and status semantics

**Files:**
- Create: `src/ontseq_platform/marlin_contracts.py`
- Create: `tests/test_marlin_contracts.py`

**Interfaces:**
- Produces: `MarlinSourceKind`, `MarlinClassificationDecision`, `MarlinProbeObservation`, `MarlinFeatureSummary`, `MarlinModelUnitScore`, `MarlinGroupedScore`, `MarlinArtifactLock`, `MarlinRuntimeCompatibilityProfile`, `MarlinPredictionReport`.
- Consumes: `GenomeBuild`, `ModuleRunStatus`, `StrictModel`, `FileFingerprint` from `ontseq_platform.models`.

- [ ] **Step 1: Write failing contract tests**

```python
from pydantic import ValidationError
import pytest

from ontseq_platform.marlin_contracts import (
    MarlinClassificationDecision,
    MarlinFeatureSummary,
    MarlinSourceKind,
)


def test_feature_summary_requires_marlin_v1_feature_count() -> None:
    with pytest.raises(ValidationError, match="357340"):
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


def test_decision_values_are_unambiguous() -> None:
    assert MarlinClassificationDecision.HIGH_CONFIDENCE.value == "HIGH_CONFIDENCE"
    assert MarlinClassificationDecision.UNKNOWN.value == "UNKNOWN"
    assert MarlinSourceKind.PRECOMPUTED_METHYLATION.value == "PRECOMPUTED_METHYLATION"
```

- [ ] **Step 2: Run tests and verify RED**

Run:

```bash
python -m pytest tests/test_marlin_contracts.py -q
```

Expected: import failure because `marlin_contracts.py` does not exist.

- [ ] **Step 3: Implement strict contracts**

Minimum definitions must include these invariants:

```python
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

`MarlinPredictionReport` must carry `ModuleRunStatus`, source kind, GRCh37 build, input fingerprint, feature summary, all 42 raw scores, grouped class/family/lineage scores, decision, top class/score, artifact lock ID, warnings/limitations, and `research_only=True`.

- [ ] **Step 4: Run contract tests and full type/lint checks for the new file**

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

### Task 2: Implement strict published five-column probe-BED import

**Files:**
- Create: `src/ontseq_platform/marlin_input.py`
- Create: `tests/test_marlin_input.py`

**Interfaces:**
- Consumes: `MarlinProbeObservation`, `MarlinSourceKind`.
- Produces: `parse_marlin_probe_bed(path: Path, *, genome_build: GenomeBuild) -> MarlinPrecomputedInput`.

- [ ] **Step 1: Write parser tests for plain/gzip, schema, coordinates and duplicates**

```python
def test_parse_probe_bed_gzip(tmp_path: Path) -> None:
    path = tmp_path / "sample.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("chr1\t10\t11\t0.75\tcg0001\n")
        handle.write("chr1\t20\t21\tNA\tcg0002\n")
    parsed = parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)
    assert parsed.observations[0].methylation_fraction == 0.75
    assert parsed.observations[1].methylation_fraction is None
    assert parsed.source_kind is MarlinSourceKind.PRECOMPUTED_METHYLATION


def test_duplicate_probe_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "dup.bed"
    path.write_text(
        "chr1\t10\t11\t0.1\tcg0001\nchr1\t10\t11\t0.9\tcg0001\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Duplicate MARLIN probe"):
        parse_marlin_probe_bed(path, genome_build=GenomeBuild.GRCH37)
```

Also add explicit failing cases for: six columns, whitespace delimiter instead of tabs, negative start, `end <= start`, `nan`, `inf`, value `<0`, value `>1`, blank probe ID, blank line policy, and GRCh38 rejection in v1.

- [ ] **Step 2: Run tests and verify RED**

```bash
python -m pytest tests/test_marlin_input.py -q
```

Expected: import/function failure.

- [ ] **Step 3: Implement deterministic parser and SHA-256 fingerprint**

Core parsing rule:

```python
fields = line.rstrip("\n\r").split("\t")
if len(fields) != 5:
    raise ValueError(f"Line {line_number}: expected exactly 5 tab-separated fields")
```

Use `gzip.open(..., "rt", encoding="utf-8")` only for `.gz`; otherwise `Path.open`. Do not accept or repair malformed rows. Compute SHA-256 over the original compressed file bytes using `sha256_file(path)`.

- [ ] **Step 4: Verify tests/lint/type**

```bash
python -m pytest tests/test_marlin_input.py -q
python -m ruff check src/ontseq_platform/marlin_input.py tests/test_marlin_input.py
python -m mypy src/ontseq_platform/marlin_input.py
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_input.py tests/test_marlin_input.py
git commit -m "feat(marlin): parse published probe methylation input"
```

---

### Task 3: Export and lock MARLIN reference resources

**Files:**
- Create: `src/ontseq_platform/marlin_artifacts.py`
- Create: `scripts/marlin_export_resources.R`
- Create: `tests/test_marlin_artifacts.py`
- Modify: `pyproject.toml` package data for `scripts/marlin_export_resources.R` when packaging is tested.

**Interfaces:**
- Produces: `create_marlin_artifact_lock(...) -> MarlinArtifactLock`, `verify_marlin_artifact_lock(lock, paths) -> None`, `load_exported_feature_ids(path) -> tuple[str, ...]`.
- Consumes: MARLIN upstream `marlin_v1.features.RData`, `marlin_v1.class_annotations.xlsx`, hg19 probe BED, model HDF5, explicit runtime versions.

- [ ] **Step 1: Write RED tests for hash mismatch, feature count, duplicate feature IDs and build mismatch**

```python
def test_lock_rejects_wrong_feature_count(tmp_path: Path) -> None:
    features = tmp_path / "features.txt"
    features.write_text("cg1\ncg2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="357340"):
        load_exported_feature_ids(features)


def test_verify_lock_detects_changed_model(valid_lock, artifact_paths) -> None:
    artifact_paths.model.write_bytes(b"changed")
    with pytest.raises(ValueError, match="model SHA-256"):
        verify_marlin_artifact_lock(valid_lock, artifact_paths)
```

- [ ] **Step 2: Run tests and verify RED**

```bash
python -m pytest tests/test_marlin_artifacts.py -q
```

Expected: missing module/functions.

- [ ] **Step 3: Implement deterministic feature-resource export script**

`scripts/marlin_export_resources.R` must do only resource conversion, not prediction:

```r
args <- commandArgs(trailingOnly=TRUE)
load(args[1])
stopifnot(exists("betas_sub_names"))
stopifnot(length(betas_sub_names) == 357340)
stopifnot(length(unique(betas_sub_names)) == 357340)
writeLines(as.character(betas_sub_names), con=args[2], useBytes=TRUE)
```

Python lock creation records both the original `.RData` SHA-256 and exported canonical feature-list SHA-256. The exported list is derived/local and can be regenerated; it is not silently substituted for the upstream artifact identity.

- [ ] **Step 4: Implement lock verification**

Require:

```python
if lock.genome_build is not GenomeBuild.GRCH37:
    raise ValueError("MARLIN v1 ONTSeq lock currently requires GRCh37/hg19")
if lock.expected_feature_count != 357340:
    raise ValueError("MARLIN v1 lock must declare 357340 model features")
if lock.expected_model_unit_count != 42:
    raise ValueError("MARLIN v1 lock must declare 42 model units")
```

Hash model, feature RData, canonical feature list, class annotations, probe BED and the pinned MARLIN code manifest/commit identity. Never download in this function.

- [ ] **Step 5: Run tests/lint/type**

```bash
python -m pytest tests/test_marlin_artifacts.py -q
python -m ruff check src/ontseq_platform/marlin_artifacts.py tests/test_marlin_artifacts.py
python -m mypy src/ontseq_platform/marlin_artifacts.py
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/ontseq_platform/marlin_artifacts.py scripts/marlin_export_resources.R tests/test_marlin_artifacts.py
git commit -m "feat(marlin): lock model and reference artifacts"
```

---

### Task 4: Build the canonical 357,340-value feature vector in Python

**Files:**
- Create: `src/ontseq_platform/marlin_features.py`
- Create: `tests/test_marlin_features.py`

**Interfaces:**
- Consumes: `MarlinPrecomputedInput`, ordered feature IDs, feature artifact SHA-256.
- Produces: `MarlinFeatureVector(values: tuple[float, ...], summary: MarlinFeatureSummary)` via `build_marlin_feature_vector(...)`.

- [ ] **Step 1: Write exact transformation and ordering tests**

Use a small internal helper test that allows a reduced test feature list while production wrapper enforces `357340`:

```python
def test_vector_is_ordered_by_locked_features_not_input_rows() -> None:
    observations = {
        "cgB": MarlinProbeObservation(chromosome="chr1", start=2, end=3, methylation_fraction=0.9, probe_id="cgB"),
        "cgA": MarlinProbeObservation(chromosome="chr1", start=1, end=2, methylation_fraction=0.49, probe_id="cgA"),
        "cgC": MarlinProbeObservation(chromosome="chr1", start=3, end=4, methylation_fraction=None, probe_id="cgC"),
    }
    vector = _build_feature_vector_for_ids(("cgA", "cgB", "cgC", "cgD"), observations)
    assert vector.values == (-1.0, 1.0, 0.0, 0.0)
```

Add boundary test `0.5 -> +1`, row-order invariance, explicit-NA versus absent counters, non-model probe counter, deterministic SHA-256.

- [ ] **Step 2: Run tests and verify RED**

```bash
python -m pytest tests/test_marlin_features.py -q
```

- [ ] **Step 3: Implement canonical vector encoding and digest**

Use a platform-independent byte contract. Recommended representation:

```python
FEATURE_TO_BYTE = {-1.0: b"\xff", 0.0: b"\x00", 1.0: b"\x01"}
payload = b"".join(FEATURE_TO_BYTE[value] for value in values)
digest = hashlib.sha256(payload).hexdigest()
```

The production `build_marlin_feature_vector` must reject a locked feature list whose length is not 357340.

- [ ] **Step 4: Run tests/lint/type**

```bash
python -m pytest tests/test_marlin_features.py -q
python -m ruff check src/ontseq_platform/marlin_features.py tests/test_marlin_features.py
python -m mypy src/ontseq_platform/marlin_features.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_features.py tests/test_marlin_features.py
git commit -m "feat(marlin): construct canonical model feature vector"
```

---

### Task 5: Add the minimal locked MARLIN inference runtime

**Files:**
- Create: `src/ontseq_platform/marlin_runtime.py`
- Create: `scripts/marlin_infer_locked.R`
- Create: `tests/test_marlin_runtime.py`

**Interfaces:**
- Consumes: `MarlinFeatureVector`, `MarlinArtifactLock`, verified paths, `CommandRunner`.
- Produces: `run_marlin_inference(...) -> MarlinRuntimeResult` with exactly 42 named scores and runtime provenance.

- [ ] **Step 1: Write RED tests using a fake `CommandRunner`**

```python
def test_runtime_rejects_41_scores(fake_lock, feature_vector, tmp_path: Path) -> None:
    runner = FakeRunner(stdout="class\tscore\n" + "\n".join(f"c{i}\t0.01" for i in range(41)))
    with pytest.raises(ValueError, match="exactly 42"):
        run_marlin_inference(feature_vector, fake_lock, runner=runner, work_dir=tmp_path)


def test_runtime_rejects_nonfinite_score(...):
    ...


def test_runtime_rejects_softmax_sum_outside_tolerance(...):
    ...
```

Also test duplicate labels, wrong label set/order, non-zero return code, missing output, and `1e-5` score-sum tolerance.

- [ ] **Step 2: Run tests and verify RED**

```bash
python -m pytest tests/test_marlin_runtime.py -q
```

- [ ] **Step 3: Implement R runtime script with no preprocessing**

Input file contains the already-built vector, one numeric value per line, in locked order. R script must refuse wrong vector length:

```r
values <- scan(args[1], what=double(), quiet=TRUE)
stopifnot(length(values) == 357340)
model <- load_model_hdf5(args[2])
pred <- as.numeric(predict(model, matrix(values, nrow=1)))
stopifnot(length(pred) == 42)
write.table(
  data.frame(model_id=seq_along(pred), score=pred),
  file=args[3], sep="\t", quote=FALSE, row.names=FALSE
)
```

Class names are bound in Python from the locked annotation resource; the runtime cannot rename classes.

- [ ] **Step 4: Implement Python invocation via existing `CommandRunner`**

Do not shell-expand. Build an argv sequence such as:

```python
argv = (
    rscript_path,
    str(inference_script),
    str(vector_path),
    str(model_path),
    str(score_output_path),
)
result = runner.run(argv, timeout_seconds=runtime_timeout_seconds)
```

Stage files in the supplied work directory; publish normalized JSON only after score validation.

- [ ] **Step 5: Run runtime unit tests/lint/type**

```bash
python -m pytest tests/test_marlin_runtime.py -q
python -m ruff check src/ontseq_platform/marlin_runtime.py tests/test_marlin_runtime.py
python -m mypy src/ontseq_platform/marlin_runtime.py
```

- [ ] **Step 6: Commit**

```bash
git add src/ontseq_platform/marlin_runtime.py scripts/marlin_infer_locked.R tests/test_marlin_runtime.py
git commit -m "feat(marlin): add locked 42-score inference runtime"
```

---

### Task 6: Freeze a runtime compatibility profile before biological validation

**Files:**
- Modify: `src/ontseq_platform/marlin_runtime.py`
- Modify: `tests/test_marlin_runtime.py`
- Create: `configs/methylation/marlin_v1.technical.yaml`

**Interfaces:**
- Produces: `create_runtime_compatibility_profile(...)`, `verify_runtime_compatibility(...)`.
- No GSE280090 outcome is consumed by this task.

- [ ] **Step 1: Write RED tests for tolerance lock**

```python
def test_compatibility_requires_same_top_model_unit() -> None:
    profile = compatibility_profile(reference_scores=(...))
    candidate = list(profile.reference_scores)
    candidate[0], candidate[1] = candidate[1], candidate[0]
    with pytest.raises(ValueError, match="top model unit"):
        verify_runtime_compatibility(profile, tuple(candidate))
```

Also test per-score difference above predeclared tolerance and incompatible runtime lock ID.

- [ ] **Step 2: Run RED tests**

```bash
python -m pytest tests/test_marlin_runtime.py -q
```

- [ ] **Step 3: Implement compatibility profile**

The profile stores: fixed synthetic/reference vector SHA-256, all 42 reference scores, absolute tolerance, score-sum tolerance, reference runtime lock ID, backend identity, and creation timestamp.

Do **not** set tolerance based on GSE data. Use the prior synthetic/fixed-vector runtime comparison to select and document a conservative tolerance; if no reproducible empirical tolerance is available in the checkout, initial profile creation remains a local explicit operator action and is not fabricated in source.

- [ ] **Step 4: Add technical config with only published/fixed semantics**

`configs/methylation/marlin_v1.technical.yaml` records `confidence_threshold: 0.8`, feature count `357340`, model-unit count `42`, build `GRCh37`, source/provenance notes, and `status: technical_defaults_and_published_semantics_only`.

- [ ] **Step 5: Verify**

```bash
python -m pytest tests/test_marlin_runtime.py -q
python -m ruff check src/ontseq_platform/marlin_runtime.py tests/test_marlin_runtime.py
python -m mypy src/ontseq_platform/marlin_runtime.py
```

- [ ] **Step 6: Commit**

```bash
git add src/ontseq_platform/marlin_runtime.py tests/test_marlin_runtime.py configs/methylation/marlin_v1.technical.yaml
git commit -m "feat(marlin): gate runtimes with frozen compatibility profile"
```

---

### Task 7: Implement class/family/lineage aggregation and confidence semantics

**Files:**
- Create: `src/ontseq_platform/marlin_classification.py`
- Create: `tests/test_marlin_classification.py`

**Interfaces:**
- Consumes: 42 raw model-unit scores plus locked class annotation rows.
- Produces: grouped class/family/lineage scores and final `MarlinPredictionReport` decision.

- [ ] **Step 1: Write RED tests for aggregation and `UNKNOWN` semantics**

```python
def test_valid_low_score_is_completed_unknown() -> None:
    report = classify_marlin_scores(
        runtime_result=runtime_result_with_grouped_top_score(0.69),
        feature_summary=feature_summary(observed=1000),
        threshold=0.8,
        ...,
    )
    assert report.status is ModuleRunStatus.COMPLETED
    assert report.decision is MarlinClassificationDecision.UNKNOWN
    assert report.top_grouped_score == pytest.approx(0.69)


def test_zero_observed_features_is_no_call() -> None:
    report = classify_marlin_scores(... feature_summary=feature_summary(observed=0), ...)
    assert report.status is ModuleRunStatus.NO_CALL
    assert report.decision is MarlinClassificationDecision.UNKNOWN
```

Also test `0.8` exactly is `HIGH_CONFIDENCE`, aggregation sums model units belonging to the same current class, and raw scores are retained unchanged.

- [ ] **Step 2: Run RED tests**

```bash
python -m pytest tests/test_marlin_classification.py -q
```

- [ ] **Step 3: Implement grouped classification**

Use annotation rows keyed by `model_id`, retaining current class, methylation-family and lineage mappings. Aggregate with deterministic iteration order. Classification uses the top **grouped current-class score**, not the top raw model unit.

- [ ] **Step 4: Verify**

```bash
python -m pytest tests/test_marlin_classification.py -q
python -m ruff check src/ontseq_platform/marlin_classification.py tests/test_marlin_classification.py
python -m mypy src/ontseq_platform/marlin_classification.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_classification.py tests/test_marlin_classification.py
git commit -m "feat(marlin): normalize grouped scores and confidence"
```

---

### Task 8: Add focused MARLIN CLI commands

**Files:**
- Create: `src/ontseq_platform/marlin_cli.py`
- Create: `tests/test_marlin_cli.py`
- Modify: `src/ontseq_platform/cli.py`
- Modify: `src/ontseq_platform/entrypoint.py`

**Interfaces:**
- Commands: `marlin-lock`, `marlin-features`, `marlin-classify`, `marlin-validate`.

- [ ] **Step 1: Write CLI discovery tests**

```python
def test_help_lists_marlin_commands(capsys, monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["ontseq", "--help"])
    entrypoint.main()
    out = capsys.readouterr().out
    assert "marlin-lock" in out
    assert "marlin-classify" in out
```

Add argument-validation tests requiring explicit `--genome-build GRCh37`, lock paths, input paths and output path.

- [ ] **Step 2: Run RED tests**

```bash
python -m pytest tests/test_marlin_cli.py -q
```

- [ ] **Step 3: Implement MARLIN CLI registration/dispatch**

Follow the existing `methylation_validation_cli.py` pattern: define `COMMANDS`, `add_subparsers`, and `run_command`. `cli.py` delegates instead of embedding MARLIN logic.

`marlin-classify` pipeline:

```text
load lock -> verify artifacts -> parse input -> build features -> verify runtime profile
-> run inference -> classify grouped scores -> write normalized JSON atomically
```

- [ ] **Step 4: Verify focused CLI and entrypoint tests**

```bash
python -m pytest tests/test_marlin_cli.py tests/test_entrypoint.py -q
python -m ruff check src/ontseq_platform/marlin_cli.py src/ontseq_platform/cli.py src/ontseq_platform/entrypoint.py tests/test_marlin_cli.py
python -m mypy src/ontseq_platform/marlin_cli.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_cli.py src/ontseq_platform/cli.py src/ontseq_platform/entrypoint.py tests/test_marlin_cli.py
git commit -m "feat(marlin): expose locked classification CLI"
```

---

### Task 9: Add external GSE280090 validation contracts and deterministic rerun checks

**Files:**
- Create: `src/ontseq_platform/marlin_validation.py`
- Create: `tests/test_marlin_validation.py`

**Interfaces:**
- Produces: `MarlinValidationManifest`, `MarlinValidationSample`, `MarlinValidationResult`, `MarlinValidationCohortReport`, `execute_marlin_validation(...)`.

- [ ] **Step 1: Write RED validation tests using synthetic local files**

```python
def test_validation_refuses_input_checksum_mismatch(tmp_path: Path, manifest) -> None:
    sample = tmp_path / "sample.txt.gz"
    write_probe_fixture(sample)
    with pytest.raises(ValueError, match="SHA-256"):
        execute_marlin_validation(manifest_with_wrong_sha(manifest), ...)


def test_repeated_identical_run_requires_same_feature_digest_and_decision(...):
    first = execute_marlin_validation(...)
    second = execute_marlin_validation(...)
    assert first.feature_vector_sha256 == second.feature_vector_sha256
    assert first.raw_score_digest == second.raw_score_digest
    assert first.decision == second.decision
```

- [ ] **Step 2: Run RED tests**

```bash
python -m pytest tests/test_marlin_validation.py -q
```

- [ ] **Step 3: Implement manifest and cohort metrics**

Each sample manifest includes accession, public source URI, filename, local SHA-256, build, expected comparison label if independently pre-registered, artifact-lock ID, runtime-profile ID and threshold-policy ID.

Cohort report computes: sample count, completed/high-confidence count, completed/unknown count, no-call count, failed count, concordant/discordant counts where expected labels exist, deterministic-rerun failures, and per-sample feature coverage.

No threshold modification is accepted by the validation runner after a manifest is locked.

- [ ] **Step 4: Verify**

```bash
python -m pytest tests/test_marlin_validation.py -q
python -m ruff check src/ontseq_platform/marlin_validation.py tests/test_marlin_validation.py
python -m mypy src/ontseq_platform/marlin_validation.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_validation.py tests/test_marlin_validation.py
git commit -m "feat(marlin): add external classification validation harness"
```

---

### Task 10: Run the first real external gate with `GSM8587229_AL_001.txt.gz`

**Files:**
- No patient-derived data committed.
- Generated local manifest/result under ignored `results/marlin-validation/`.
- Update after successful execution: `docs/MARLIN_CLASSIFICATION.md` with the exact local validation command and evidence boundary, but never embed patient-derived rows.

**Interfaces:**
- Consumes local file `C:\Users\sxhul\Downloads\GSM8587229_AL_001.txt.gz` through the WSL-visible path when executed on the user's workstation.
- Produces a local checksummed manifest and normalized JSON result.

- [ ] **Step 1: Fingerprint the exact local input before prediction**

Example on WSL:

```bash
sha256sum "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz"
```

Record this hash in the generated validation manifest; do not hardcode a guessed hash into source.

- [ ] **Step 2: Inspect schema without changing the file**

```bash
gzip -cd "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" | head -n 5
gzip -cd "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" | awk -F '\t' 'NR<=100 {if (NF!=5) exit 1} END {print "first-100-rows-schema-ok"}'
```

If the actual GEO restricted file does not match `marlin_probe_bed_v1`, stop and add a separately specified adapter; do not coerce it through the parser.

- [ ] **Step 3: Run locked classification twice**

```bash
ontseq marlin-classify \
  --input "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" \
  --genome-build GRCh37 \
  --artifact-lock <local-marlin-lock.json> \
  --runtime-profile <local-runtime-profile.json> \
  --output results/marlin-validation/AL_001.run1.json

ontseq marlin-classify \
  --input "/mnt/c/Users/sxhul/Downloads/GSM8587229_AL_001.txt.gz" \
  --genome-build GRCh37 \
  --artifact-lock <local-marlin-lock.json> \
  --runtime-profile <local-runtime-profile.json> \
  --output results/marlin-validation/AL_001.run2.json
```

- [ ] **Step 4: Verify deterministic invariants**

Require identical input SHA-256, feature-vector SHA-256, top grouped class, decision and scores within the frozen runtime compatibility tolerance. Do not require byte-identical TensorFlow output unless the frozen profile itself requires it.

- [ ] **Step 5: Compare to independently registered publication expectation**

The comparison target must be registered before reading the ONTSeq result. For AL_001, use only a publication-supported label/score expectation; if the publication does not support an exact scalar expectation under the identical file/runtime, record only the supported class/confidence-level expectation.

- [ ] **Step 6: Commit only code/docs resulting from the gate**

Do not commit the GSE file or patient-derived result payload. Commit sanitized technical findings and commands only.

---

### Task 11: Run cohort validation over available GSE280090 restricted files

**Files:**
- Local external inputs only.
- Generated local results under ignored `results/marlin-validation/gse280090-v1/`.
- Update: `docs/MARLIN_CLASSIFICATION.md` with aggregate, non-identifying technical metrics and evidence boundary.

**Interfaces:**
- Consumes: locked cohort manifest generated from locally downloaded restricted files.
- Produces: `MarlinValidationCohortReport`.

- [ ] **Step 1: Build the cohort manifest from local files and public accession metadata**

No automatic downloader in normal ONTSeq runtime. Each local file gets SHA-256 and explicit `GRCh37` declaration.

- [ ] **Step 2: Execute with the same artifact lock, compatibility profile and fixed `0.8` confidence policy**

```bash
ontseq marlin-validate \
  --manifest results/marlin-validation/gse280090-v1/manifest.json \
  --artifact-lock <local-marlin-lock.json> \
  --runtime-profile <local-runtime-profile.json> \
  --output results/marlin-validation/gse280090-v1/report.json
```

- [ ] **Step 3: Verify no threshold/preprocessing drift**

The report must echo the same lock/profile/policy digests for all samples. A mixed identity aborts cohort aggregation rather than producing pooled metrics.

- [ ] **Step 4: Record aggregate metrics**

Required metrics: total samples, successful runtime executions, `HIGH_CONFIDENCE`, `UNKNOWN`, `NO_CALL`, failures, concordant/discordant where truth/comparison target exists, feature-coverage distribution, runtime distribution, deterministic-rerun failures.

- [ ] **Step 5: Commit sanitized validation documentation only**

```bash
git add docs/MARLIN_CLASSIFICATION.md
git commit -m "docs(marlin): record GSE280090 downstream validation boundary"
```

---

### Task 12: Add the future native `MODKIT_DERIVED` bridge behind a disabled gate

**Files:**
- Modify: `src/ontseq_platform/marlin_input.py`
- Modify: `src/ontseq_platform/marlin_features.py`
- Create/modify tests in `tests/test_marlin_input.py`, `tests/test_marlin_features.py`

**Interfaces:**
- Produces a probe-level adapter contract only; does not claim equivalence until same-specimen bridge data exist.

- [ ] **Step 1: Write RED test that native bridge is disabled without explicit validated bridge identity**

```python
def test_modkit_derived_source_is_refused_without_bridge_lock() -> None:
    with pytest.raises(ValueError, match="native MARLIN bridge is not validated"):
        build_marlin_from_modkit(..., bridge_lock=None)
```

- [ ] **Step 2: Implement probe-coordinate mapping contract without enabling production use**

Use the locked `marlin_v1.probes_hg19.bed.gz` coordinate resource, preserving probe IDs and per-probe CpG fractions. Do not use chromosome/target aggregates from `methylation.py`.

- [ ] **Step 3: Add bridge-lock schema that requires same-specimen evidence before enablement**

The bridge lock stores precomputed-vs-native feature-vector comparison metrics and must remain absent in default configs until such data exist.

- [ ] **Step 4: Verify tests**

```bash
python -m pytest tests/test_marlin_input.py tests/test_marlin_features.py -q
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_input.py src/ontseq_platform/marlin_features.py tests/test_marlin_input.py tests/test_marlin_features.py
git commit -m "feat(marlin): define gated native modkit bridge contract"
```

---

### Task 13: Documentation, packaging and command-surface integration

**Files:**
- Create: `docs/MARLIN_CLASSIFICATION.md`
- Modify: `pyproject.toml`
- Modify: `CHANGELOG.md`
- Modify: relevant README command documentation.

**Interfaces:**
- Packages R scripts/configs; documents install-time external artifacts and evidence boundaries.

- [ ] **Step 1: Add documentation tests or static assertions where existing repository conventions support them**

Ensure `ontseq --help` lists commands and package-data tests confirm both MARLIN R scripts are included in built wheel/runtime bundle.

- [ ] **Step 2: Package scripts/config**

Add data-file entries for:

```toml
"share/ontseq/scripts" = [
  "scripts/run_qdnaseq_ace.R",
  "scripts/marlin_infer_locked.R",
  "scripts/marlin_export_resources.R",
]
"share/ontseq/configs/methylation" = ["configs/methylation/*.yaml", "configs/methylation/*.json"]
```

Do not package the MARLIN model or patient-derived validation data into the Python wheel.

- [ ] **Step 3: Document exact evidence ladder**

`docs/MARLIN_CLASSIFICATION.md` must explicitly state:

```text
Synthetic parser/feature tests
    < locked runtime compatibility
    < GSE processed-CpG downstream reproduction
    < same-specimen modkit bridge
    < analytical intended-use validation
    < clinical validation
```

- [ ] **Step 4: Build wheel and inspect contents**

```bash
python -m build
python - <<'PY'
import zipfile, glob
wheel = sorted(glob.glob('dist/*.whl'))[-1]
with zipfile.ZipFile(wheel) as z:
    names = set(z.namelist())
    assert any(name.endswith('marlin_infer_locked.R') for name in names)
    assert any(name.endswith('marlin_export_resources.R') for name in names)
print('MARLIN package data present')
PY
```

- [ ] **Step 5: Commit**

```bash
git add docs/MARLIN_CLASSIFICATION.md pyproject.toml CHANGELOG.md README.md
git commit -m "docs(marlin): document and package classification lane"
```

---

### Task 14: Full repository verification and independent review

**Files:**
- No new functional scope.
- Fix only defects found by verification/review, each with a regression test.

**Interfaces:**
- Produces a review-ready branch; no merge yet.

- [ ] **Step 1: Run MARLIN-focused suite**

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

- [ ] **Step 2: Run complete Python test suite**

```bash
python -m pytest -q
```

- [ ] **Step 3: Run lint/format/type/schema/version guards**

```bash
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python scripts/check_version_consistency.py
```

Expected at this stage: version consistency remains `0.8.2` until Task 15; all other checks PASS.

- [ ] **Step 4: Run build/install smoke in a clean environment**

Build wheel, install into a fresh virtual environment, run `ontseq --help`, parse a synthetic `.gz` probe fixture and execute fake-runner MARLIN classification tests from the installed package.

- [ ] **Step 5: Independent review**

Review specifically for: accidental clinical claims, build inference, duplicate-probe behavior, hash TOCTOU gaps, score-label ordering, threshold application to wrong score layer, runtime output partial-write hazards, Windows/WSL path behavior, and whether GSE evidence is mislabeled as native modkit validation.

- [ ] **Step 6: Fix every accepted review finding with a failing regression test first**

Each defect follows RED -> minimal fix -> focused PASS -> full relevant suite.

- [ ] **Step 7: Commit review fixes**

```bash
git add <reviewed-files>
git commit -m "fix(marlin): address independent integration review"
```

---

### Task 15: Coordinate the `0.9.0` release identity only after all gates pass

**Files:**
- Modify: `pyproject.toml`
- Modify: `uv.lock`
- Modify: `src/ontseq_platform/__init__.py`
- Modify: `CITATION.cff`
- Modify: `desktop/ONTSeq.Desktop/ONTSeq.Desktop.csproj`
- Modify: `desktop/ONTSeq.Desktop/Version.cs`
- Modify: `desktop/desktop.settings.example.json`
- Modify: current-version Desktop/operator docs required by `scripts/check_version_consistency.py`
- Modify: `.github/workflows/desktop-ci.yml`
- Modify: `CHANGELOG.md`
- Modify: README current-status section.

**Interfaces:**
- Produces a repository-consistent `0.9.0` engineering release identity; still Research Use Only.

- [ ] **Step 1: Read `scripts/check_version_consistency.py` and enumerate every required current-release surface**

Do not perform a partial version bump.

- [ ] **Step 2: Change all current release declarations from `0.8.2` to `0.9.0` together**

Changelog wording must explicitly say that 0.9.0 adds an engineering/audit MARLIN classification lane and does not represent clinical validation.

- [ ] **Step 3: Regenerate/update dependency lock using the repository's existing lock workflow**

Verify the editable `ontseq-platform` entry is `0.9.0`.

- [ ] **Step 4: Run version guard**

```bash
python scripts/check_version_consistency.py
```

Expected: `Version consistency check passed`.

- [ ] **Step 5: Re-run full repository verification**

```bash
python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python scripts/check_version_consistency.py
```

- [ ] **Step 6: Commit release identity**

```bash
git add pyproject.toml uv.lock src/ontseq_platform/__init__.py CITATION.cff desktop .github/workflows/desktop-ci.yml CHANGELOG.md README.md docs
git commit -m "chore(release): stage ONTSeq 0.9.0 engineering release"
```

---

### Task 16: Open PR, verify CI, and update Issue #76 without merging prematurely

**Files:**
- GitHub metadata only.

**Interfaces:**
- Produces a reviewable PR linked to Issue #76.

- [ ] **Step 1: Compare branch against `main` and confirm only intended files changed**

```bash
git diff --stat main...feat/marlin-classification-v1
git diff --check main...feat/marlin-classification-v1
```

- [ ] **Step 2: Open PR with explicit evidence boundary**

PR body must include:

- design/spec path;
- implementation-plan path;
- test counts/commands;
- runtime compatibility evidence;
- AL_001/GSE status, if executed;
- explicit statement that processed-CpG validation does not validate POD5/Dorado/MM-ML/modkit;
- no clinical claim;
- `Closes #76` only if every acceptance criterion in Issue #76 is actually satisfied; otherwise use `Relates to #76`.

- [ ] **Step 3: Wait for/check all required CI statuses**

No merge on partial/unknown CI.

- [ ] **Step 4: Perform final PR diff review against the approved spec**

Confirm all new contracts retain RUO semantics and no hidden download, threshold tuning or build inference was introduced.

- [ ] **Step 5: Update Issue #76 with evidence and remaining gaps**

If cohort validation or native bridge remains incomplete, leave the issue open and state the exact unmet acceptance criteria.

- [ ] **Step 6: Merge only after explicit final review/approval**

Preferred merge behavior should match repository convention; do not bypass CI or force-update `main`.

---

## Plan Self-Review

### Spec coverage

- Strict five-column precomputed input: Tasks 2, 10.
- GRCh37-only v1 and no liftover: Tasks 1–4, 9–10.
- 357,340 ordered features and exact binarization: Tasks 3–4.
- Locked model/features/classes/probe resources/runtime: Tasks 3, 5–6.
- Exactly 42 finite softmax scores and score-sum invariant: Task 5.
- Grouped current-class/family/lineage aggregation: Task 7.
- `HIGH_CONFIDENCE` vs `UNKNOWN`; zero-evidence `NO_CALL`: Task 7.
- Runtime compatibility frozen before biological validation: Task 6 before Tasks 10–11.
- No GSE payload in Git; checksummed local validation: Tasks 9–11.
- AL_001 first external gate: Task 10.
- Full GSE280090 restricted cohort: Task 11.
- Native modkit bridge kept separate/gated: Task 12.
- Packaging/docs: Task 13.
- Full verification and independent review: Task 14.
- Consistent 0.9.0 version surface: Task 15.
- PR/CI/evidence boundary: Task 16.

### Placeholder scan

The plan contains no `TBD`/`TODO` implementation placeholders. Local artifact paths/hashes that cannot be truthfully known in source are obtained by explicit lock/fingerprint commands rather than guessed.

### Type/interface consistency

- `MarlinFeatureSummary` is defined in Task 1 and consumed in Tasks 4 and 7.
- `MarlinArtifactLock` is defined in Task 1 and created/verified in Task 3.
- `MarlinFeatureVector` is produced in Task 4 and consumed in Task 5.
- `MarlinRuntimeCompatibilityProfile` is defined in Task 1 and created/verified in Task 6.
- Raw 42-score runtime output is produced in Task 5 and grouped only in Task 7.
- `MarlinPredictionReport` is assembled in Task 7 and consumed by validation in Task 9.
- The GSE path never becomes evidence for `MODKIT_DERIVED`; the native bridge remains Task 12 and separately gated.
