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
- Two execution locks may differ in `lock_id`, runtime versions, backend, source URI formatting and timestamps, but must have identical canonical immutable artifact-set identity before candidate inference.
- Canonical artifact-set identity includes model/code/features/class-annotation/probe hashes, build, expected feature/model-unit counts and preprocessing/runtime contract versions; it excludes runtime identity and timestamps.
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
    artifact_set_sha256: str
    model_version: str
    model_sha256: str
    code_version_or_commit: str
    code_manifest_sha256: str
    feature_sha256: str
    canonical_feature_list_sha256: str
    class_annotation_sha256: str
    probe_resource_sha256: str
    genome_build: GenomeBuild
    expected_feature_count: Literal[357340] = 357340
    expected_model_unit_count: Literal[42] = 42
    preprocessing_contract_version: Literal["marlin-v1-binarize-1"] = "marlin-v1-binarize-1"
    runtime_contract_version: Literal["marlin-runtime-v1"] = "marlin-runtime-v1"
    research_only: Literal[True] = True

class MarlinRuntimeIdentity(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    runtime_id: str
    R_version: str
    keras_version: str
    tensorflow_version: str
    python_version: str
    execution_backend: Literal["cpu", "gpu"]
    runtime_contract_version: Literal["marlin-runtime-v1"] = "marlin-runtime-v1"
    created_at: datetime
    research_only: Literal[True] = True

class MarlinDualRuntimeCompatibilityReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    comparison_id: str
    artifact_set_identity: MarlinArtifactSetIdentity
    reference_artifact_lock_id: str
    candidate_artifact_lock_id: str
    feature_vector_sha256: str
    reference_profile_id: str
    reference_runtime_identity: MarlinRuntimeIdentity
    candidate_runtime_identity: MarlinRuntimeIdentity
    absolute_score_tolerance: float
    score_sum_tolerance: float
    reference_scores: list[float]
    candidate_scores: list[float]
    absolute_differences: list[float]
    max_absolute_difference: float
    reference_top_model_unit_index: int
    candidate_top_model_unit_index: int
    all_scores_within_tolerance: bool
    top_model_unit_matches: bool
    softmax_invariants_pass: bool
    verdict: MarlinDualRuntimeVerdict
    created_at: datetime
    research_only: Literal[True] = True
```

- [ ] **Step 1: Write RED contract tests**

Add tests proving timezone-aware timestamps, exact 42-element score/difference vectors, finite `[0,1]` scores, non-negative finite differences, `max_absolute_difference == max(absolute_differences)`, and report verdict consistency.

```python
def test_dual_runtime_report_rejects_pass_when_score_gate_failed() -> None:
    payload = _valid_dual_runtime_report_payload()
    payload["verdict"] = "PASS"
    payload["all_scores_within_tolerance"] = False
    with pytest.raises(ValueError, match="PASS"):
        MarlinDualRuntimeCompatibilityReport.model_validate(payload)
```

Also test that a `MarlinRuntimeIdentity` with a naive `created_at` is rejected.

- [ ] **Step 2: Run RED tests**

```bash
python -m pytest tests/test_marlin_contracts.py -q
```

Expected: FAIL because the new contract classes do not exist.

- [ ] **Step 3: Implement the minimal strict models and validators**

`MarlinDualRuntimeCompatibilityReport` validator must independently recompute:

```python
expected_max = max(self.absolute_differences)
if not math.isclose(self.max_absolute_difference, expected_max, rel_tol=0, abs_tol=1e-15):
    raise ValueError("MARLIN dual-runtime max difference is inconsistent")

expected_pass = (
    self.all_scores_within_tolerance
    and self.top_model_unit_matches
    and self.softmax_invariants_pass
)
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

- [ ] **Step 1: Write RED artifact identity tests**

Construct two valid locks with identical immutable fields but different `lock_id`, source URIs, R/Python/Keras/TensorFlow versions, backend and timestamps. Assert equal `artifact_set_sha256`.

Then mutate each immutable field one at a time (`model_sha256`, `code_manifest_sha256`, `feature_sha256`, `canonical_feature_list_sha256`, `class_annotation_sha256`, `probe_resource_sha256`, `genome_build` where model validation permits a synthetic dict-level test, `preprocessing_contract_version`, `runtime_contract_version`) and require a different digest or model rejection.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime_compare.py -q
```

Expected: FAIL because `marlin_runtime_compare.py` does not exist.

- [ ] **Step 3: Implement canonical digest**

Use an explicit JSON-compatible canonical payload with sorted keys and compact separators; never hash `model_dump_json()` wholesale because future unrelated fields could silently alter identity.

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

Source URIs are deliberately excluded.

- [ ] **Step 4: Implement live-probe identity binding**

`runtime_identity_from_probe()` must require exact equality between probe values and the corresponding execution-lock runtime fields. Reuse the same semantics as `verify_runtime_probe_matches_lock`; do not accept typed version overrides.

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

### Task 3: Implement pure 42-score dual-runtime comparison

**Files:**
- Modify: `src/ontseq_platform/marlin_runtime_compare.py`
- Modify: `tests/test_marlin_runtime_compare.py`

**Interfaces:**
- Consumes: Task 2 identity functions, existing `MarlinRuntimeResult`, existing `MarlinRuntimeCompatibilityProfile`.
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

Add independent tests for:

```python
def test_dual_runtime_exact_tolerance_boundary_passes() -> None:
    # reference unit 1 = 0.5000000; candidate unit 1 = 0.5000001
    # compensate another score by -1e-7 so both sums remain 1.0
    report = _compare_with_delta(1e-7)
    assert report.verdict is MarlinDualRuntimeVerdict.PASS


def test_dual_runtime_one_score_above_tolerance_fails() -> None:
    report = _compare_with_delta(1.0000001e-7)
    assert report.verdict is MarlinDualRuntimeVerdict.FAIL
```

Also test top-unit mismatch yields FAIL and a candidate vector violating score-sum tolerance cannot produce PASS.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime_compare.py -q
```

- [ ] **Step 3: Implement comparison without exceptions for valid incompatibility**

Precondition/evidence corruption remains an exception. A structurally valid candidate inference that differs numerically returns a report with `verdict=FAIL`.

```python
reference_scores = tuple(reference_profile.reference_scores)
candidate_scores = tuple(item.score for item in candidate_result.model_scores)
differences = tuple(abs(left - right) for left, right in zip(reference_scores, candidate_scores, strict=True))
all_within = all(delta <= reference_profile.absolute_score_tolerance for delta in differences)
reference_top = reference_profile.top_model_unit_index
candidate_top = max(range(42), key=candidate_scores.__getitem__)
softmax_ok = (
    abs(sum(reference_scores) - 1.0) <= reference_profile.score_sum_tolerance
    and abs(sum(candidate_scores) - 1.0) <= reference_profile.score_sum_tolerance
)
verdict = MarlinDualRuntimeVerdict.PASS if all_within and reference_top == candidate_top and softmax_ok else MarlinDualRuntimeVerdict.FAIL
```

Require:
- `reference_profile.reference_runtime_lock_id == reference_lock.lock_id`;
- `candidate_result.artifact_lock_id == candidate_lock.lock_id`;
- `candidate_result.feature_vector_sha256 == reference_profile.feature_vector_sha256`;
- artifact-set identity matches both locks;
- runtime identities match their respective locks.

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
- Modify: `src/ontseq_platform/marlin_runtime_compare.py`
- Modify: `tests/test_marlin_runtime_compare.py`

**Interfaces:**
- Consumes: existing `_load_frozen_runtime_fixture` behavior via a new public wrapper in `marlin_runtime.py` if needed, `probe_marlin_runtime()`, `run_marlin_inference()`, artifact verification helpers.
- Produces:

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

- [ ] **Step 1: Write RED ordering tests with a recording fake runner**

The fake runner records probe/inference calls. Assert artifact-set mismatch or fixture digest mismatch raises before any candidate inference command is issued.

```python
with pytest.raises(ValueError, match="artifact set"):
    execute_marlin_dual_runtime_comparison(...)
assert recording_runner.inference_calls == 0
```

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime_compare.py -q
```

- [ ] **Step 3: Expose a narrow public fixture loader if required**

If orchestration cannot reuse the existing private `_load_frozen_runtime_fixture()` cleanly, rename it to `load_frozen_runtime_fixture()` in `marlin_runtime.py`, retain behavior exactly, and update its internal callers/tests. Do not create a second fixture parser.

- [ ] **Step 4: Implement exact execution order**

Order is mandatory:

```text
same artifact-set preflight
-> reference identity/profile checks
-> load fixture + require profile digest
-> candidate live probe
-> candidate identity binding
-> one candidate inference
-> pure comparison
```

Do not execute the reference runtime again in this command; the reference oracle is the previously frozen profile plus its separately persisted live-probed reference runtime identity.

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
- Consumes: Task 4 orchestration, existing `_write_model_atomic()`, `load_model()` and CLI delegation pattern.
- Produces command:

```text
ontseq marlin-compare-runtimes \
  --reference-artifact-lock <reference-lock.json> \
  --candidate-artifact-lock <candidate-lock.json> \
  --reference-runtime-identity <reference-runtime-identity.json> \
  --reference-profile <reference-profile.json> \
  --runtime-fixture <runtime-fixture.txt> \
  --model <marlin-model.hdf5> \
  --inference-script scripts/marlin_infer_locked.R \
  --candidate-runtime-probe-script scripts/marlin_runtime_probe.R \
  --candidate-rscript <candidate-env/bin/Rscript> \
  --comparison-id MARLIN_V1_R413_VS_R423_CPU \
  --output results/marlin-validation/runtime-comparison.json
```

- [ ] **Step 1: Write RED parser/help tests**

```python
def test_global_help_lists_marlin_compare_runtimes(capsys, monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["ontseq", "--help"])
    entrypoint.main()
    assert "marlin-compare-runtimes" in capsys.readouterr().out
```

Require every input above and reject an existing output path to prevent silent evidence overwrite.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_cli_compare.py -q
```

- [ ] **Step 3: Implement parser and command handler**

Add `"marlin-compare-runtimes"` to `COMMANDS`. The handler loads both locks, reference identity and profile, invokes Task 4 with `SubprocessRunner()`, and writes the report using existing atomic JSON persistence.

Use a new helper:

```python
def _require_new_output(path: Path, *, label: str) -> Path:
    candidate = Path(path)
    if candidate.exists():
        raise FileExistsError(f"{label} output already exists: {candidate}")
    return candidate
```

Do not alter `_require_new_freeze_outputs()` semantics.

- [ ] **Step 4: Add fake-runner integration test at the Python API boundary**

The CLI parser test stays subprocess-free. End-to-end model execution remains covered in `test_marlin_runtime_compare.py` with fake probe/inference outputs.

- [ ] **Step 5: Verify GREEN**

```bash
python -m pytest tests/test_marlin_cli_compare.py tests/test_marlin_cli.py tests/test_marlin_cli_freeze.py -q
python -m ruff check src/ontseq_platform/marlin_cli.py src/ontseq_platform/entrypoint.py tests/test_marlin_cli_compare.py
python -m mypy src/ontseq_platform/marlin_cli.py
```

- [ ] **Step 6: Commit**

```bash
git add src/ontseq_platform/marlin_cli.py src/ontseq_platform/entrypoint.py tests/test_marlin_cli_compare.py
git commit -m "feat(marlin): expose dual-runtime comparison CLI"
```

---

### Task 6: Persist reference runtime identity during freeze without breaking existing outputs

**Files:**
- Modify: `src/ontseq_platform/marlin_runtime_freeze.py`
- Modify: `src/ontseq_platform/marlin_cli.py`
- Modify: `tests/test_marlin_runtime_freeze.py`
- Modify: `tests/test_marlin_cli_freeze.py`

**Interfaces:**
- Consumes: Task 1 `MarlinRuntimeIdentity`, Task 2 `runtime_identity_from_probe()`.
- Produces optional CLI argument:

```text
--output-runtime-identity <reference-runtime-identity.json>
```

- [ ] **Step 1: Write RED backward-compatibility tests**

Existing `marlin-freeze-runtime` invocation without `--output-runtime-identity` must still parse and retain its two original required outputs. A new invocation with the optional argument writes a third evidence artifact representing the exact live reference probe used before freeze.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_marlin_runtime_freeze.py tests/test_marlin_cli_freeze.py -q
```

- [ ] **Step 3: Extend freeze result without duplicating probe execution**

The live probe already occurs in `_run_freeze()`. Create `MarlinRuntimeIdentity` from that same probe and reference execution lock. Do not launch a second probe solely for persistence.

If `--output-runtime-identity` is supplied, require it to be a third distinct non-existing path and atomically write it. Failure writing any requested freeze evidence must remove newly written siblings from that invocation.

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

### Task 7: Documentation, packaging and evidence-boundary regression tests

**Files:**
- Modify: `docs/MARLIN_RUNTIME_FREEZE.md`
- Modify: `docs/MARLIN_CLASSIFICATION.md`
- Modify: `configs/methylation/marlin_v1.technical.yaml`
- Modify: `scripts/check_wheel_resources.py` only if the new Python module is not automatically packaged by the existing package discovery.
- Modify: `tests/test_repository_safety.py` only if needed to make the no-biological-payload contract explicit.

**Interfaces:**
- Consumes: completed CLI/contracts.
- Produces: operator sequence and precise claim language.

- [ ] **Step 1: Add the reference qualification workflow**

Document:

```text
A. build/reconstruct reference MARLIN runtime
B. create reference execution lock with live probe
C. marlin-freeze-runtime --output-runtime-identity ...
D. create candidate execution lock with candidate live probe
E. marlin-compare-runtimes ...
F. require verdict PASS before AL_001/GSE280090
```

State explicitly that if the reference runtime resolves to anything other than the upstream documented environment, the observed versions are reported as-is and must not be relabeled as R 4.1.3.

- [ ] **Step 2: Freeze policy identifiers in technical config**

Add non-biological policy metadata without changing existing classifier semantics:

```yaml
dual_runtime_compatibility:
  status: engineering_gate
  absolute_score_tolerance: 1.0e-7
  score_sum_tolerance: 1.0e-5
  fixture_generator_version: sha256-index-mod3-v1
  research_only: true
```

- [ ] **Step 3: Verify docs/config/package tests**

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

Only stage paths that actually changed.

---

### Task 8: Full automated verification and external real-model qualification instructions

**Files:**
- No production-code change unless verification finds a defect.
- Update PR #77 description only after all automated gates are green, summarizing evidence level accurately.

**Interfaces:**
- Consumes: Tasks 1-7.
- Produces: verified feature branch plus a reproducible operator procedure for the model-containing environment.

- [ ] **Step 1: Run targeted MARLIN suite**

```bash
python -m pytest \
  tests/test_marlin_contracts.py \
  tests/test_marlin_runtime.py \
  tests/test_marlin_runtime_freeze.py \
  tests/test_marlin_runtime_compare.py \
  tests/test_marlin_cli.py \
  tests/test_marlin_cli_freeze.py \
  tests/test_marlin_cli_compare.py \
  tests/test_marlin_validation.py -q
```

Expected: PASS.

- [ ] **Step 2: Run complete Python quality gates**

```bash
python -m pytest -q
python -m ruff check src tests scripts
python -m ruff format --check src tests scripts
python -m mypy src/ontseq_platform
python scripts/check_repository_safety.py
```

Expected: PASS.

- [ ] **Step 3: Verify GitHub Actions for the final head**

Require successful:

```text
CI
Desktop CI
Desktop bundle contract
MARLIN runtime smoke
CodeQL (when triggered)
```

Do not claim completion from local tests alone.

- [ ] **Step 4: Perform final code review before declaring engineering-complete**

Review the complete Task 1-7 diff against the spec, with special attention to:
- no runtime fields inside canonical artifact-set digest;
- no source URI influence on digest;
- artifact-set/fixture mismatch occurs before inference;
- valid numerical mismatch yields report FAIL rather than command failure;
- malformed evidence fails hard;
- exact-boundary `1e-7` behavior;
- reference identity is the actual live-probed environment used for freeze;
- no clinical/biological claim leakage.

- [ ] **Step 5: Real-model qualification remains external and must use controlled local resources**

On the workstation/runtime that has the trained MARLIN model, execute:

```text
1. Reference runtime: live probe -> reference execution lock.
2. Freeze deterministic fixture + reference profile + reference runtime identity.
3. Candidate runtime: live probe -> candidate execution lock using the identical model/artifact bytes.
4. Run marlin-compare-runtimes.
5. Archive the machine-readable report locally.
6. Proceed to AL_001 only if verdict == PASS.
```

The trained model, generated score profile, GSE data and patient-derived data remain outside Git.

- [ ] **Step 6: Commit only verification-derived fixes if necessary**

If no code changes are required, do not manufacture a verification commit.
