# MARLIN Dual-Runtime Compatibility Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a fail-closed RUO dual-runtime MARLIN qualification path that proves whether a live candidate runtime reproduces a frozen reference runtime's 42 model scores for the same immutable classifier artifact set and deterministic 357,340-feature fixture.

**Architecture:** Preserve the existing single-runtime `MarlinArtifactLock` and `marlin-freeze-runtime` behavior. Add a canonical immutable artifact-set identity derived from two execution locks, explicit live-probed runtime identities, a pure numerical comparison layer, and a separate `ontseq marlin-compare-runtimes` CLI command that writes one atomic machine-readable PASS/FAIL report. No biological data or model payloads enter Git.

**Tech Stack:** Python 3.11+, Pydantic strict models, existing `CommandRunner`/`SubprocessRunner`, existing MARLIN R/Keras/TensorFlow inference bridge, pytest, Ruff, mypy, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-09-15-marlin-dual-runtime-compatibility-design.md`

## Global Constraints

- Research Use Only; no analytical or clinical validity claim.
- Keep ONTSeq version at `0.8.2` for this feature.
- `main` is not modified or merged without explicit user approval.
- PR #77 remains Draft.
- Existing `ontseq marlin-freeze-runtime` semantics remain backward compatible.
- Reference and candidate executions use separate `MarlinArtifactLock` objects.
- Two execution locks may differ in `lock_id`, source URI formatting, runtime versions, backend and timestamps, but must have identical canonical immutable artifact-set identity before candidate inference.
- Canonical artifact-set identity includes model/code/features/class-annotation/probe hashes, model/code version identity, build, expected feature/model-unit counts and preprocessing/runtime contract versions. It excludes runtime identity, source URI formatting, local paths and timestamps.
- Candidate runtime identity is created only from a live runtime probe; manually claimed versions are not evidence.
- First engineering tolerances remain `absolute_score_tolerance = 1e-7` and `score_sum_tolerance = 1e-5`; do not tune after inspecting AL_001/GSE280090.
- Exactly 357,340 fixture values and exactly 42 model scores are required.
- No patient data, GSE payloads, trained model bytes, credentials or large runtime artifacts are committed.
- Artifact-set or fixture mismatch must fail before candidate model execution.
- Numerical incompatibility is PASS/FAIL engineering evidence, never `UNKNOWN`, `NO_CALL` or a biological diagnosis.

---

### Task 1: Add strict dual-runtime provenance contracts

**Files:**
- Modify: `src/ontseq_platform/marlin_contracts.py`
- Modify: `tests/test_marlin_contracts.py`

**Interfaces:**
- Consumes: existing `MarlinArtifactLock`, `MarlinRuntimeCompatibilityProfile`, `StrictModel`.
- Produces:

```python
class MarlinDualRuntimeVerdict(StrEnum):
    PASS = "PASS"
    FAIL = "FAIL"

class MarlinArtifactSetIdentity(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    artifact_set_sha256: str = Field(pattern=_SHA256)
    model_version: str = Field(min_length=1)
    model_sha256: str = Field(pattern=_SHA256)
    code_version_or_commit: str = Field(min_length=1)
    code_manifest_sha256: str = Field(pattern=_SHA256)
    feature_sha256: str = Field(pattern=_SHA256)
    canonical_feature_list_sha256: str = Field(pattern=_SHA256)
    class_annotation_sha256: str = Field(pattern=_SHA256)
    probe_resource_sha256: str = Field(pattern=_SHA256)
    genome_build: GenomeBuild
    expected_feature_count: Literal[357340] = 357340
    expected_model_unit_count: Literal[42] = 42
    preprocessing_contract_version: Literal["marlin-v1-binarize-1"] = "marlin-v1-binarize-1"
    runtime_contract_version: Literal["marlin-runtime-v1"] = "marlin-runtime-v1"
    research_only: Literal[True] = True

class MarlinRuntimeIdentity(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    runtime_id: str = Field(min_length=3)
    R_version: str = Field(min_length=1)
    keras_version: str = Field(min_length=1)
    tensorflow_version: str = Field(min_length=1)
    python_version: str = Field(pattern=r"^\d+\.\d+\.\d+$")
    execution_backend: Literal["cpu", "gpu"]
    runtime_contract_version: Literal["marlin-runtime-v1"] = "marlin-runtime-v1"
    created_at: datetime
    research_only: Literal[True] = True

class MarlinDualRuntimeCompatibilityReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    comparison_id: str = Field(min_length=3)
    artifact_set_identity: MarlinArtifactSetIdentity
    reference_artifact_lock_id: str = Field(min_length=3)
    candidate_artifact_lock_id: str = Field(min_length=3)
    feature_vector_sha256: str = Field(pattern=_SHA256)
    reference_profile_id: str = Field(min_length=1)
    reference_runtime_identity: MarlinRuntimeIdentity
    candidate_runtime_identity: MarlinRuntimeIdentity
    absolute_score_tolerance: float = Field(gt=0, le=1)
    score_sum_tolerance: float = Field(gt=0, le=1)
    reference_scores: list[float] = Field(min_length=42, max_length=42)
    candidate_scores: list[float] = Field(min_length=42, max_length=42)
    absolute_differences: list[float] = Field(min_length=42, max_length=42)
    max_absolute_difference: float = Field(ge=0)
    reference_top_model_unit_index: int = Field(ge=0, le=41)
    candidate_top_model_unit_index: int = Field(ge=0, le=41)
    all_scores_within_tolerance: bool
    top_model_unit_matches: bool
    softmax_invariants_pass: bool
    verdict: MarlinDualRuntimeVerdict
    created_at: datetime
    research_only: Literal[True] = True
```

- [ ] **Step 1: Write RED contract tests**

Add tests for timezone-aware timestamps, exact 42-element score/difference vectors, finite `[0,1]` scores, non-negative finite differences, consistent `max_absolute_difference`, and verdict consistency.

```python
def test_dual_runtime_report_rejects_pass_when_score_gate_failed() -> None:
    payload = _valid_dual_runtime_report_payload()
    payload["verdict"] = "PASS"
    payload["all_scores_within_tolerance"] = False
    with pytest.raises(ValueError, match="PASS"):
        MarlinDualRuntimeCompatibilityReport.model_validate(payload)
```

Add a test that a naive `MarlinRuntimeIdentity.created_at` is rejected.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_contracts.py -q
```

Expected: FAIL because the new contracts do not exist.

- [ ] **Step 3: Implement minimal models and validators**

The report validator independently checks every score/difference is finite, both score vectors are probability vectors, both timestamps are timezone-aware, and:

```python
expected_max = max(self.absolute_differences)
if not math.isclose(self.max_absolute_difference, expected_max, rel_tol=0, abs_tol=1e-15):
    raise ValueError("MARLIN dual-runtime max difference is inconsistent")
expected_pass = self.all_scores_within_tolerance and self.top_model_unit_matches and self.softmax_invariants_pass
if (self.verdict is MarlinDualRuntimeVerdict.PASS) != expected_pass:
    raise ValueError("MARLIN dual-runtime PASS verdict is inconsistent with gate results")
```

Do not embed artifact-lock comparison logic in the Pydantic model.

- [ ] **Step 4: Verify GREEN**

```bash
python -m pytest tests/test_marlin_contracts.py -q
python -m ruff check src/ontseq_platform/marlin_contracts.py tests/test_marlin_contracts.py
python -m mypy src/ontseq_platform/marlin_contracts.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_contracts.py tests/test_marlin_contracts.py
git commit -m "feat(marlin): add dual-runtime evidence contracts"
```

---

### Task 2: Derive canonical artifact-set identity and live runtime identity

**Files:**
- Create: `src/ontseq_platform/marlin_runtime_compare.py`
- Create: `tests/test_marlin_runtime_compare.py`

**Interfaces:**
- Consumes: `MarlinArtifactLock`, `MarlinRuntimeProbeReport`, Task 1 contracts.
- Produces:

```python
def derive_marlin_artifact_set_identity(lock: MarlinArtifactLock) -> MarlinArtifactSetIdentity: ...

def runtime_identity_from_probe(
    lock: MarlinArtifactLock,
    probe: MarlinRuntimeProbeReport,
    *,
    runtime_id: str,
    created_at: datetime,
) -> MarlinRuntimeIdentity: ...

def require_same_marlin_artifact_set(
    reference_lock: MarlinArtifactLock,
    candidate_lock: MarlinArtifactLock,
) -> MarlinArtifactSetIdentity: ...
```

- [ ] **Step 1: Write RED identity tests**

Construct two valid locks with identical immutable fields but different `lock_id`, source URIs, R/Python/Keras/TensorFlow versions, backend and timestamps. Assert equal `artifact_set_sha256`. Mutating any immutable SHA, `model_version`, `code_version_or_commit`, build, expected counts or contract version must change the digest or be rejected by the existing strict lock model.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime_compare.py -q
```

Expected: FAIL because `marlin_runtime_compare.py` does not exist.

- [ ] **Step 3: Implement canonical digest**

Use exactly this explicit payload and canonical encoding:

```python
payload = {
    "model_version": lock.model_version,
    "model_sha256": lock.model_sha256,
    "code_version_or_commit": lock.code_version_or_commit,
    "code_manifest_sha256": lock.code_manifest_sha256,
    "feature_sha256": lock.feature_sha256,
    "canonical_feature_list_sha256": lock.canonical_feature_list_sha256,
    "class_annotation_sha256": lock.class_annotation_sha256,
    "probe_resource_sha256": lock.probe_resource_sha256,
    "genome_build": lock.genome_build.value,
    "expected_feature_count": lock.expected_feature_count,
    "expected_model_unit_count": lock.expected_model_unit_count,
    "preprocessing_contract_version": lock.preprocessing_contract_version,
    "runtime_contract_version": lock.runtime_contract_version,
}
encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
digest = hashlib.sha256(encoded).hexdigest()
```

Never hash the whole lock serialization. Source URIs are excluded.

- [ ] **Step 4: Implement runtime identity binding**

Require exact equality between live probe and corresponding lock fields. Use the same field semantics as existing `verify_runtime_probe_matches_lock()` and create the identity only after all values match.

- [ ] **Step 5: Verify GREEN**

```bash
python -m pytest tests/test_marlin_runtime_compare.py -q
python -m ruff check src/ontseq_platform/marlin_runtime_compare.py tests/test_marlin_runtime_compare.py
python -m mypy src/ontseq_platform/marlin_runtime_compare.py
```

- [ ] **Step 6: Commit**

```bash
git add src/ontseq_platform/marlin_runtime_compare.py tests/test_marlin_runtime_compare.py
git commit -m "feat(marlin): separate artifact-set and runtime identities"
```

---

### Task 3: Implement pure 42-score comparison

**Files:**
- Modify: `src/ontseq_platform/marlin_runtime_compare.py`
- Modify: `tests/test_marlin_runtime_compare.py`

**Interfaces:**
- Produces:

```python
def compare_marlin_runtime_results(
    *,
    comparison_id: str,
    artifact_set_identity: MarlinArtifactSetIdentity,
    reference_lock: MarlinArtifactLock,
    candidate_lock: MarlinArtifactLock,
    reference_runtime_identity: MarlinRuntimeIdentity,
    candidate_runtime_identity: MarlinRuntimeIdentity,
    reference_profile: MarlinRuntimeCompatibilityProfile,
    candidate_result: MarlinRuntimeResult,
    created_at: datetime,
) -> MarlinDualRuntimeCompatibilityReport: ...
```

- [ ] **Step 1: Write RED numerical policy tests**

```python
def test_dual_runtime_exact_tolerance_boundary_passes() -> None:
    report = _compare_with_compensated_delta(1e-7)
    assert report.verdict is MarlinDualRuntimeVerdict.PASS


def test_dual_runtime_one_score_above_tolerance_fails() -> None:
    report = _compare_with_compensated_delta(1.0000001e-7)
    assert report.verdict is MarlinDualRuntimeVerdict.FAIL
```

Add top-unit mismatch -> FAIL and valid 42-score softmax drift -> FAIL tests.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime_compare.py -q
```

- [ ] **Step 3: Implement pure comparison**

Malformed/inconsistent evidence raises. Structurally valid numerical incompatibility returns a FAIL report.

```python
reference_scores = tuple(reference_profile.reference_scores)
candidate_scores = tuple(item.score for item in candidate_result.model_scores)
differences = tuple(abs(a - b) for a, b in zip(reference_scores, candidate_scores, strict=True))
all_within = all(delta <= reference_profile.absolute_score_tolerance for delta in differences)
reference_top = reference_profile.top_model_unit_index
candidate_top = max(range(42), key=candidate_scores.__getitem__)
softmax_ok = (
    abs(sum(reference_scores) - 1.0) <= reference_profile.score_sum_tolerance
    and abs(sum(candidate_scores) - 1.0) <= reference_profile.score_sum_tolerance
)
verdict = MarlinDualRuntimeVerdict.PASS if all_within and reference_top == candidate_top and softmax_ok else MarlinDualRuntimeVerdict.FAIL
```

Before comparison require the reference profile belongs to the reference lock, candidate result belongs to the candidate lock, candidate fixture digest matches the reference profile, the supplied artifact-set identity matches both locks, and both runtime identities match their execution locks.

- [ ] **Step 4: Verify GREEN**

```bash
python -m pytest tests/test_marlin_runtime_compare.py -q
python -m ruff check src/ontseq_platform/marlin_runtime_compare.py tests/test_marlin_runtime_compare.py
python -m mypy src/ontseq_platform/marlin_runtime_compare.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_runtime_compare.py tests/test_marlin_runtime_compare.py
git commit -m "feat(marlin): compare frozen and candidate runtime scores"
```

---

### Task 4: Add fail-closed candidate execution orchestration

**Files:**
- Modify: `src/ontseq_platform/marlin_runtime.py`
- Modify: `src/ontseq_platform/marlin_runtime_compare.py`
- Modify: `tests/test_marlin_runtime.py`
- Modify: `tests/test_marlin_runtime_compare.py`

**Interfaces:**
- Rename existing private `_load_frozen_runtime_fixture()` to public `load_frozen_runtime_fixture()` without changing parsing/hash semantics. Update `verify_frozen_runtime_fixture()` to call the public function.
- Produce:

```python
def execute_marlin_dual_runtime_comparison(
    *,
    comparison_id: str,
    reference_lock: MarlinArtifactLock,
    candidate_lock: MarlinArtifactLock,
    reference_runtime_identity: MarlinRuntimeIdentity,
    reference_profile: MarlinRuntimeCompatibilityProfile,
    runtime_fixture_path: Path,
    model_path: Path,
    inference_script: Path,
    candidate_probe_script: Path,
    runner: CommandRunner,
    work_dir: Path,
    candidate_rscript_path: str,
    created_at: datetime,
) -> MarlinDualRuntimeCompatibilityReport: ...
```

- [ ] **Step 1: Write RED ordering tests**

Use a recording fake runner. Artifact-set mismatch and fixture digest mismatch must raise before any candidate inference invocation. Candidate probe must run before candidate inference.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime_compare.py tests/test_marlin_runtime.py -q
```

- [ ] **Step 3: Publicize the single canonical fixture loader**

Rename `_load_frozen_runtime_fixture` to `load_frozen_runtime_fixture`; update only its callers and tests. Do not duplicate fixture parsing.

- [ ] **Step 4: Implement exact orchestration order**

```text
same artifact-set preflight
-> reference identity/profile checks
-> load fixture and require profile digest
-> candidate live probe
-> candidate identity binding
-> candidate inference once
-> pure comparison report
```

The reference runtime is not re-executed here; its frozen profile plus persisted live-probed runtime identity are the oracle.

- [ ] **Step 5: Verify GREEN**

```bash
python -m pytest tests/test_marlin_runtime_compare.py tests/test_marlin_runtime.py -q
python -m ruff check src/ontseq_platform/marlin_runtime_compare.py src/ontseq_platform/marlin_runtime.py tests/test_marlin_runtime_compare.py tests/test_marlin_runtime.py
python -m mypy src/ontseq_platform/marlin_runtime_compare.py src/ontseq_platform/marlin_runtime.py
```

- [ ] **Step 6: Commit**

```bash
git add src/ontseq_platform/marlin_runtime_compare.py src/ontseq_platform/marlin_runtime.py tests/test_marlin_runtime_compare.py tests/test_marlin_runtime.py
git commit -m "feat(marlin): execute fail-closed dual-runtime qualification"
```

---

### Task 5: Add `marlin-compare-runtimes` CLI and atomic report persistence

**Files:**
- Modify: `src/ontseq_platform/marlin_cli.py`
- Modify: `src/ontseq_platform/entrypoint.py`
- Create: `tests/test_marlin_cli_compare.py`

**Interfaces:**
- Add command `marlin-compare-runtimes` to `COMMANDS` with required arguments:

```text
--reference-artifact-lock
--candidate-artifact-lock
--reference-runtime-identity
--reference-profile
--runtime-fixture
--model
--inference-script
--candidate-runtime-probe-script
--candidate-rscript
--comparison-id
--output
```

- [ ] **Step 1: Write RED help/parser tests**

Global help must list `marlin-compare-runtimes`; missing required args must raise `SystemExit`; existing output must raise `FileExistsError` before execution.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_cli_compare.py -q
```

- [ ] **Step 3: Implement parser and handler**

Load both execution locks, reference identity and reference profile via `load_model()`, call `execute_marlin_dual_runtime_comparison()` with `SubprocessRunner()`, and persist with existing `_write_model_atomic()`.

Add:

```python
def _require_new_output(path: Path, *, label: str) -> Path:
    candidate = Path(path)
    if candidate.exists():
        raise FileExistsError(f"{label} output already exists: {candidate}")
    return candidate
```

Do not change `_require_new_freeze_outputs()`.

- [ ] **Step 4: Verify GREEN**

```bash
python -m pytest tests/test_marlin_cli_compare.py tests/test_marlin_cli.py tests/test_marlin_cli_freeze.py -q
python -m ruff check src/ontseq_platform/marlin_cli.py src/ontseq_platform/entrypoint.py tests/test_marlin_cli_compare.py
python -m mypy src/ontseq_platform/marlin_cli.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_cli.py src/ontseq_platform/entrypoint.py tests/test_marlin_cli_compare.py
git commit -m "feat(marlin): expose dual-runtime comparison CLI"
```

---

### Task 6: Persist the exact reference runtime identity during freeze

**Files:**
- Modify: `src/ontseq_platform/marlin_runtime_freeze.py`
- Modify: `src/ontseq_platform/marlin_cli.py`
- Modify: `tests/test_marlin_runtime_freeze.py`
- Modify: `tests/test_marlin_cli_freeze.py`

**Interfaces:**
- Add optional `marlin-freeze-runtime --output-runtime-identity <path>`.
- Use the already executed live probe from `_run_freeze()`; do not run a second probe.

- [ ] **Step 1: Write RED backward-compatibility tests**

Existing freeze invocation without the new option still requires and writes only fixture/profile. With the option, a third distinct non-existing output path is required and contains a `MarlinRuntimeIdentity` matching the reference execution lock and live probe.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime_freeze.py tests/test_marlin_cli_freeze.py -q
```

- [ ] **Step 3: Implement reference identity persistence**

Create identity with `runtime_identity_from_probe(reference_lock, live_runtime, runtime_id=f"{reference_lock.lock_id}:runtime", created_at=<same freeze timestamp>)`. Use one timezone-aware timestamp for profile and identity in that invocation. Atomically write every requested output; if any requested write fails, remove newly written siblings from that invocation.

- [ ] **Step 4: Verify GREEN**

```bash
python -m pytest tests/test_marlin_runtime_freeze.py tests/test_marlin_cli_freeze.py -q
python -m ruff check src/ontseq_platform/marlin_runtime_freeze.py src/ontseq_platform/marlin_cli.py tests/test_marlin_runtime_freeze.py tests/test_marlin_cli_freeze.py
python -m mypy src/ontseq_platform/marlin_runtime_freeze.py src/ontseq_platform/marlin_cli.py
```

- [ ] **Step 5: Commit**

```bash
git add src/ontseq_platform/marlin_runtime_freeze.py src/ontseq_platform/marlin_cli.py tests/test_marlin_runtime_freeze.py tests/test_marlin_cli_freeze.py
git commit -m "feat(marlin): persist frozen reference runtime identity"
```

---

### Task 7: Documentation, technical policy and wheel contract

**Files:**
- Modify: `docs/MARLIN_RUNTIME_FREEZE.md`
- Modify: `docs/MARLIN_CLASSIFICATION.md`
- Modify: `configs/methylation/marlin_v1.technical.yaml`
- Modify: `scripts/check_wheel_resources.py`
- Modify: `tests/test_repository_safety.py`

**Interfaces:**
- Add `ontseq_platform/marlin_runtime_compare.py` to `REQUIRED_SUFFIXES` in `scripts/check_wheel_resources.py`.
- Add technical config:

```yaml
dual_runtime_compatibility:
  status: engineering_gate
  absolute_score_tolerance: 1.0e-7
  score_sum_tolerance: 1.0e-5
  fixture_generator_version: sha256-index-mod3-v1
  research_only: true
```

- [ ] **Step 1: Document the operator sequence**

Document exactly:

```text
A. reconstruct/reference MARLIN environment and live-probe it
B. create reference execution lock
C. freeze fixture/profile/reference runtime identity
D. create candidate execution lock from its live probe using identical artifact bytes
E. run marlin-compare-runtimes
F. require report verdict PASS before AL_001/GSE280090
```

If the reference environment resolves to versions other than upstream's documented stack, record actual observed versions and never relabel them as R 4.1.3.

- [ ] **Step 2: Add wheel and safety assertions**

The wheel contract must require the new module. Repository-safety test must explicitly reject committed files named `marlin-runtime-profile.json`, `runtime-comparison.json`, `GSM8587229_AL_001.txt.gz`, and common trained-model extensions under validation result paths, while preserving existing fixture/code files.

- [ ] **Step 3: Verify docs/config/package gates**

```bash
python -m pytest tests/test_repository_safety.py tests/test_marlin_runtime_compare.py tests/test_marlin_cli_compare.py -q
python -m ruff check src tests scripts
python -m ruff format --check src tests scripts
python -m mypy src/ontseq_platform
```

- [ ] **Step 4: Commit**

```bash
git add docs/MARLIN_RUNTIME_FREEZE.md docs/MARLIN_CLASSIFICATION.md configs/methylation/marlin_v1.technical.yaml scripts/check_wheel_resources.py tests/test_repository_safety.py
git commit -m "docs(marlin): define dual-runtime qualification gate"
```

---

### Task 8: Full automated verification and real-model qualification gate

**Files:**
- No planned production-code edits. Verification defects are fixed in the owning Task 1-7 file and re-run through its targeted tests before this final gate.

- [ ] **Step 1: Run targeted MARLIN suite**

```bash
python -m pytest tests/test_marlin_contracts.py tests/test_marlin_runtime.py tests/test_marlin_runtime_freeze.py tests/test_marlin_runtime_compare.py tests/test_marlin_cli.py tests/test_marlin_cli_freeze.py tests/test_marlin_cli_compare.py tests/test_marlin_validation.py -q
```

- [ ] **Step 2: Run complete Python quality gates**

```bash
python -m pytest -q
python -m ruff check src tests scripts
python -m ruff format --check src tests scripts
python -m mypy src/ontseq_platform
python scripts/check_repository_safety.py
```

- [ ] **Step 3: Require successful GitHub Actions for final head**

```text
CI
Desktop CI
Desktop bundle contract
MARLIN runtime smoke
CodeQL when triggered
```

- [ ] **Step 4: Final code review**

Review the Task 1-7 diff against the spec. Mandatory findings to check: runtime fields absent from artifact-set digest; source URIs absent from digest; artifact-set/fixture mismatch before inference; valid numerical incompatibility returns report FAIL; malformed evidence raises; exact `1e-7` boundary passes; reference identity came from the same live probe used for freeze; no biological/clinical claim leakage.

- [ ] **Step 5: Execute controlled real-model qualification outside Git**

```text
1. Reference live probe -> reference execution lock.
2. Freeze deterministic fixture + reference profile + reference runtime identity.
3. Candidate live probe -> candidate execution lock over identical artifact bytes.
4. Run marlin-compare-runtimes.
5. Archive the machine-readable report locally.
6. Proceed to AL_001 only when verdict == PASS.
```

The trained model, generated scores/profile, GSE data and patient-derived data remain outside Git.
