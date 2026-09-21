"""One-shot PR85 repair control; excluded from the resulting feature commit."""
from pathlib import Path
import subprocess
import os


def replace(path, old, new):
    p = Path(path)
    text = p.read_text()
    if text.count(old) != 1:
        raise RuntimeError(f"Expected one exact replacement in {path}: {old[:70]!r}")
    p.write_text(text.replace(old, new), encoding="utf-8")


REGRESSION_TEXT = '''

@pytest.mark.parametrize("version", ["4.6.2.0.1", "4.6.2.0-SNAPSHOT", "4.6.2.0+local"])
def test_review_rejects_nonexact_gatk_build(config_dict, tmp_path, version):
    from ontseq_platform.gatk_adapter.core import run_mutect2

    fake = SyntheticExecutor()

    def execute(step, workdir, timeout):
        code = fake(step, workdir, timeout)
        if step.name == "version":
            (workdir / "version.stdout.log").write_text(
                "The Genome Analysis Toolkit (GATK) v" + version
            )
        return code

    result = run_mutect2(enabled(config_dict), tmp_path / "out", executor=execute)
    assert result.status == "FAILED"
    assert "mutect2" not in fake.calls
'''
tests = Path("tests/test_gatk_adapter.py")
if "test_review_rejects_nonexact_gatk_build" in tests.read_text():
    raise RuntimeError("Review regressions already present; refusing duplicate repair")
tests.write_text(tests.read_text() + REGRESSION_TEXT, encoding="utf-8")
replace("tests/test_gatk_mutect2.py",
        '["4.6.1.0", "4.6.2.01", "4.7.0.0", "garbage"]',
        '["4.6.1.0", "4.6.2.01", "4.7.0.0", "garbage", "4.6.2.0-SNAPSHOT", "4.6.2.0+local"]')
red = subprocess.run(
    ["python", "-m", "pytest", "-q", "tests/test_gatk_adapter.py", "tests/test_gatk_mutect2.py",
     "-k", "review_rejects_nonexact or version_mismatch_blocks_calling"],
    env=dict(os.environ, PYTHONPATH="src"), text=True, capture_output=True,
)
print(red.stdout)
print(red.stderr)
if red.returncode != 1 or "5 failed" not in red.stdout:
    raise RuntimeError("Expected five exact-version regressions to fail before repair")

BASE = Path("src/ontseq_platform/gatk_adapter")
replace(BASE / "contracts.py", "GATK_VERSION = '4.6.2.0'  # Explicit compatibility baseline, not a claim of latest release.", "# Explicit compatibility baseline, not a claim of latest release.\nGATK_VERSION: Literal['4.6.2.0'] = '4.6.2.0'")
replace(BASE / "contracts.py", "subjects = [self.tumor, self.germline_resource]", "subjects: list[BamSample | VcfResource] = [self.tumor, self.germline_resource]")
replace(BASE / "core.py", "from typing import Callable", "from collections.abc import Callable\nfrom typing import Literal")
replace(BASE / "core.py", "from .contracts import LockedFile, Mutect2Config, Mutect2Result", "from .contracts import LockedFile, Mutect2Config, Mutect2Result, SmallVariantCandidate")
replace(BASE / "core.py", "        sites = str(config.contamination_sites.vcf.path)", "        if config.contamination_sites is None:\n            raise ValueError('Contamination estimation requires a sites resource')\n        sites = str(config.contamination_sites.vcf.path)")
replace(BASE / "core.py", "        options = {'umask': 0o077} if os.name == 'posix' else {}\n", "")
replace(BASE / "core.py", "check=False, **options,", "check=False, umask=0o077 if os.name == 'posix' else -1,")
replace(BASE / "core.py", "    execution_evidence = 'NONE'", "    execution_evidence: Literal['NONE', 'TEST_DOUBLE', 'NATIVE_EXECUTION'] = 'NONE'")
replace(BASE / "core.py", "    candidates = ()\n    messages = [", "    candidates: tuple[SmallVariantCandidate, ...] = ()\n    messages = [")
replace(BASE / "core.py", "    status = 'FAILED'", "    status: Literal['FAILED', 'COMPLETED'] = 'FAILED'")
replace(BASE / "core.py", "                if _bounded_text(step.outputs[0]).strip() != sample.sample_id:", "                if sample is None:\n                    raise ValueError('Normal-sample step requires an explicit normal')\n                if _bounded_text(step.outputs[0]).strip() != sample.sample_id:")
replace(BASE / "core.py", r"(\d+\.\d+\.\d+\.\d+)'", r"(\d+\.\d+\.\d+\.\d+)(?![\w.+-])'")
replace(BASE / "cli.py", "    config = Mutect2Config(\n", "    normal = None\n    if normal_bam is not None:\n        if normal_sample is None:\n            raise ValueError('A normal BAM requires its sample name')\n        normal = _sample(normal_bam, normal_sample, fasta.sha256)\n    config = Mutect2Config(\n")
replace(BASE / "cli.py", "normal=_sample(normal_bam, normal_sample, fasta.sha256) if normal_bam else None,", "normal=normal,")
replace(BASE / "cli.py", "description='Experimental local ONT/Mutect2 adapter. Research Use Only. No clinical release.'", "description=(\n            'Experimental local ONT/Mutect2 adapter. '\n            'Research Use Only. No clinical release.'\n        )")
replace(BASE / "cli.py", "'notice': 'Planning only; no native tools, integrity checks or biological analysis ran.',", "'notice': (\n                    'Planning only; no native tools, integrity checks '\n                    'or biological analysis ran.'\n                ),")
replace(BASE / "vcf.py", "from typing import TextIO", "from typing import Literal, TextIO, TypedDict, overload")
replace(BASE / "vcf.py", "def _numbers(value: str | None, length: int, *, integer: bool = False) -> tuple:", '''@overload
def _numbers(value: str | None, length: int, *, integer: Literal[True]) -> tuple[int | None, ...]: ...


@overload
def _numbers(
    value: str | None, length: int, *, integer: Literal[False] = False
) -> tuple[float | None, ...]: ...


def _numbers(
    value: str | None, length: int, *, integer: bool = False
) -> tuple[int | float | None, ...]:''')
replace(BASE / "vcf.py", "    parsed = []", "    parsed: list[int | float | None] = []")
replace(BASE / "vcf.py", "def _sample_fields(keys: list[str], value: str, n_alt: int) -> dict:", '''class SampleFields(TypedDict):
    af: tuple[float | None, ...]
    ad: tuple[int | None, ...]
    dp: int | None


def _sample_fields(keys: list[str], value: str, n_alt: int) -> SampleFields:''')
replace(BASE / "vcf.py", "def _variant_type(ref: str, alt: str) -> str:", "def _variant_type(ref: str, alt: str) -> Literal['SNV', 'MNV', 'INS', 'DEL', 'COMPLEX']:")
replace(BASE / "vcf.py", "            info = {}", "            info: dict[str, str | Literal[True]] = {}")
replace("src/ontseq_platform/gatk_mutect2.py", r"(?![\d.])", r"(?![\w.+-])")

replace("docs/gatk_adapter/README.md", '''It has **not** been committed, pushed, merged, or tested in the full upstream repository.
A complete repository checkout was not available in this execution environment.''', '''This package is tracked in PR #85. Merge eligibility requires full-repository checks on the
exact reviewed revision; the original standalone test results are historical, not a current
integration qualification. See the PR checks and REVIEW.md for engineering scope.''')
replace("examples/gatk/README.md", "**after applying the patch**", "with the reviewed adapter package present")
replace("examples/gatk/README.md", "Those upstream gates have not been executed in this delivery environment.", "The exact reviewed revision must pass all of those upstream gates before merge.")
replace("docs/gatk_adapter/CHANGELOG_FRAGMENT.md", "# Proposed changelog fragment (not yet inserted into upstream CHANGELOG.md)", "# GATK adapter change summary (also recorded in CHANGELOG.md)")
replace("docs/gatk_adapter/CHANGELOG_FRAGMENT.md", "qualification is claimed. Not committed or merged into the upstream repository.", "qualification is claimed. See PR #85 for the current integration status.")
replace("docs/gatk_adapter/VALIDATION_IMPACT.md", "## Not executed or not delivered", "## Original standalone-delivery boundary\n\nThis table describes the original standalone delivery. GitHub upload and software-gate status\nare superseded by PR #85 checks on the exact reviewed revision. Native GATK and analytical\nvalidation remain unqualified; see REVIEW.md.\n\n### Historical status")
replace("CHANGELOG.md", "## Unreleased\n", '''## Unreleased

- Add opt-in GATK 4.6.2.0 research adapters, separate candidate/config/result schemas,
  tumor-only and paired modes, native artifact provenance and synthetic contract tests.
- Repair strict typing and exact-version checks; reject suffixed or extended GATK builds.
  Include both adapter entry points in focused CI without weakening full-project gates.
- Validation impact: adapters remain separately invoked, non-reportable and not
  real-GATK-qualified or analytically validated. No canonical runner, CNV/SV policy,
  reference bundle, package version or automatic clinical release changes.
''')
Path("docs/gatk_adapter/REVIEW.md").write_text('''# PR #85 engineering review

Research Use Only. Human review required. No automatic clinical release.

## Reviewed scope

Review base: af9a1ba35f12fb9fbbeb3774789b259cd4794c2d.
Main integration base: 2c7aa2bf1a2d01f3f51aa75232c608407f8329b9.

Both explicitly invoked adapters are retained. `gatk_mutect2` is the legacy,
interval-bounded CLI; `gatk_adapter` is the separate versioned Pydantic API/CLI.
Their manifests, outputs and missing-resource policies are not interchangeable.
Neither is registered for automatic production analysis or clinical reporting.

## Findings and repairs

- Uploaded standalone files failed strict typing and repository lint/format gates.
  Add explicit literal, union, tuple, numeric and optional-input types; retain the
  existing Pydantic schemas and runtime opt-in defaults.
- Both version gates accepted some suffixed build identifiers as the pinned version.
  Add failing regressions before correction and reject extended/local/snapshot IDs
  before Mutect2 is called. Do not relax the fixed GATK 4.6.2.0 baseline.
- Focused CI omitted the newly added package. Include both implementations, tests
  and versioned schemas, with formatting enforced rather than merely displayed.
- Correct stale standalone-upload statements and add the root changelog entry.

## Verification required for this revision

The repair workflow runs the newly added regressions in RED before changes and
then the full project gates: `make safety`, `make versions`, `make lint`,
`make test`, schema freshness and wheel-resource verification. The final PR also
runs the repository's usual CI; only exact-head results are evidence for merging.
No ignore, skip, exclusion or reportability threshold is added by this repair.

## Unchanged qualification limits

Software contract tests use synthetic stand-ins. They do not execute GATK and
cannot establish native BAM/JAR interoperability, ONT sensitivity, specificity,
VAF bias, LoD, somatic origin or clinical validity. The versioned extractor
preserves native anchored alleles without left normalization or locus-by-locus
REF validation. Legacy and new outputs must not be joined by naive string
identity or treated as independent confirmation. All original native artifacts
remain the review authority. Real-tool and assay qualification remain separate.
''', encoding="utf-8")
