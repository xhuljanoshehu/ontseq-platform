from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def replace_once(path: Path, old: str, new: str, *, label: str) -> None:
    text = path.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{label}: expected one anchor, found {count}")
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


def write_input_digest_helper() -> None:
    path = ROOT / "src/ontseq_platform/pipeline/input_digest.py"
    path.write_text(
        '''"""Run-scoped content digests for large immutable analysis inputs.

Stage plans are rebuilt during every invocation so content-addressed resume can decide
whether a stage is reusable. Several plans fingerprint the same multi-gigabyte BAM. This
cache removes duplicate reads *within one RunContext* while preserving the deliberate
re-verification reads used by intake and release verification.
"""

from __future__ import annotations

from pathlib import Path

from .envelope import sha256_file


FileIdentity = tuple[str, int, int, int, int, int]


class RunInputDigestCache:
    """Memoize stable input digests for one pipeline run invocation only."""

    def __init__(self) -> None:
        self._digests: dict[FileIdentity, str] = {}

    @staticmethod
    def _identity(path: Path) -> FileIdentity:
        resolved = path.resolve(strict=True)
        stat = resolved.stat()
        return (
            str(resolved),
            int(stat.st_dev),
            int(stat.st_ino),
            int(stat.st_size),
            int(stat.st_mtime_ns),
            int(stat.st_ctime_ns),
        )

    def digest(self, path: Path) -> tuple[str, bool]:
        """Return SHA-256 plus whether the file stayed unchanged while it was read.

        The cache key includes the resolved path and filesystem identity/state. A changed
        file therefore misses the cache. A digest is inserted only after a second stat
        proves that the same target, size, mtime and ctime survived the full read.
        """

        before = self._identity(path)
        cached = self._digests.get(before)
        if cached is not None:
            return cached, True

        digest = sha256_file(Path(before[0]))
        try:
            after = self._identity(path)
        except OSError:
            return digest, False
        stable = before == after
        if stable:
            self._digests[before] = digest
        return digest, stable
''',
        encoding="utf-8",
    )


def write_input_digest_tests() -> None:
    path = ROOT / "tests/test_run_input_digest_cache.py"
    path.write_text(
        '''from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from ontseq_platform.pipeline.envelope import sha256_file
from ontseq_platform.pipeline.input_digest import RunInputDigestCache


def test_same_run_reads_an_unchanged_input_once(tmp_path: Path) -> None:
    source = tmp_path / "large.bam"
    source.write_bytes(b"A" * 1024)
    cache = RunInputDigestCache()

    with patch(
        "ontseq_platform.pipeline.input_digest.sha256_file", wraps=sha256_file
    ) as digest_file:
        first = cache.digest(source)
        second = cache.digest(source)

    assert first == second
    assert first[1] is True
    assert digest_file.call_count == 1


def test_mutation_invalidates_the_run_scoped_cache(tmp_path: Path) -> None:
    source = tmp_path / "sample.bam"
    source.write_bytes(b"before")
    cache = RunInputDigestCache()
    first, stable = cache.digest(source)
    assert stable

    source.write_bytes(b"after-with-a-different-size")
    second, stable = cache.digest(source)

    assert stable
    assert first != second


def test_an_unstable_read_is_never_cached(tmp_path: Path) -> None:
    source = tmp_path / "sample.bam"
    source.write_bytes(b"initial")
    cache = RunInputDigestCache()
    calls = 0

    def changing_digest(path: Path) -> str:
        nonlocal calls
        calls += 1
        digest = sha256_file(path)
        if calls == 1:
            source.write_bytes(b"changed-after-the-read")
        return digest

    with patch("ontseq_platform.pipeline.input_digest.sha256_file", side_effect=changing_digest):
        first, first_stable = cache.digest(source)
        second, second_stable = cache.digest(source)

    assert first_stable is False
    assert second_stable is True
    assert first != second
    assert calls == 2


def test_digest_state_is_not_shared_between_run_contexts(tmp_path: Path) -> None:
    source = tmp_path / "sample.bam"
    source.write_bytes(b"same bytes")
    first_run = RunInputDigestCache()
    second_run = RunInputDigestCache()

    with patch(
        "ontseq_platform.pipeline.input_digest.sha256_file", wraps=sha256_file
    ) as digest_file:
        first_run.digest(source)
        second_run.digest(source)

    assert digest_file.call_count == 2
''',
        encoding="utf-8",
    )


def patch_runner() -> None:
    path = ROOT / "src/ontseq_platform/pipeline/runner.py"
    replace_once(
        path,
        "from .envelope import Artifact, RunEnvelope, sha256_file, stage_signature\n",
        "from .envelope import Artifact, RunEnvelope, sha256_file, stage_signature\n"
        "from .input_digest import RunInputDigestCache\n",
        label="runner import",
    )
    replace_once(
        path,
        "    artifacts: dict[StageId, list[Artifact]] = field(default_factory=dict)\n\n    @property\n",
        "    artifacts: dict[StageId, list[Artifact]] = field(default_factory=dict)\n"
        "    input_digests: RunInputDigestCache = field(\n"
        "        default_factory=RunInputDigestCache, repr=False\n"
        "    )\n\n"
        "    def fingerprint_external_input(\n"
        "        self, path: Path, *, label: str | None = None\n"
        "    ) -> tuple[str, str]:\n"
        "        if not path.is_file():\n"
        "            raise StageFailure(\"required external input is missing\")\n"
        "        try:\n"
        "            digest, stable = self.input_digests.digest(path)\n"
        "        except OSError as exc:\n"
        "            raise StageFailure(\"required external input could not be fingerprinted\") from exc\n"
        "        if not stable:\n"
        "            raise StageFailure(\n"
        "                f\"{label or 'required external input'} changed while it was being fingerprinted\"\n"
        "            )\n"
        "        return (label or path.name, digest)\n\n"
        "    @property\n",
        label="RunContext cache",
    )
    replace_once(
        path,
        "def _external_fingerprint(path: Path, *, label: str | None = None) -> tuple[str, str]:\n"
        "    \"\"\"Fingerprint an input from outside the envelope by name and content.\"\"\"\n"
        "    digest, stable = _stable_digest(path)\n"
        "    if not stable:\n"
        "        raise StageFailure(\n"
        "            f\"{label or 'required external input'} changed while it was being fingerprinted\"\n"
        "        )\n"
        "    return (label or path.name, digest)\n",
        "def _external_fingerprint(\n"
        "    ctx: RunContext, path: Path, *, label: str | None = None\n"
        ") -> tuple[str, str]:\n"
        "    \"\"\"Fingerprint a plan input through this run's stable-digest cache.\"\"\"\n"
        "    return ctx.fingerprint_external_input(path, label=label)\n",
        label="external fingerprint helper",
    )
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"(?<!def )_external_fingerprint\(", "_external_fingerprint(ctx, ", text)
    path.write_text(text, encoding="utf-8")

    replace_once(
        path,
        "    if policy.cpg_only and ctx.config.reference_fasta is None:\n"
        "        raise StageFailure(\n"
        "            \"CpG-restricted methylation requires the locked reference FASTA\"\n"
        "        )\n",
        "    if ctx.config.reference_fasta is None:\n"
        "        raise StageFailure(\n"
        "            \"modkit --modified-bases requires the locked reference FASTA for every \"\n"
        "            \"methylation run\"\n"
        "        )\n",
        label="methylation plan reference gate",
    )


def patch_cnv_extension() -> None:
    path = ROOT / "src/ontseq_platform/cnv/extension.py"
    replace_once(
        path,
        "    external_inputs = [\n"
        "        (bam.name, sha256_file(bam)),\n"
        "        (settings.script.name, sha256_file(settings.script)),\n"
        "    ]\n"
        "    if ctx.config.annotation_cache is not None:\n"
        "        external_inputs.append((\"annotation_cache\", sha256_file(ctx.config.annotation_cache)))\n",
        "    external_inputs = [\n"
        "        ctx.fingerprint_external_input(bam),\n"
        "        ctx.fingerprint_external_input(settings.script),\n"
        "    ]\n"
        "    if ctx.config.annotation_cache is not None:\n"
        "        external_inputs.append(\n"
        "            ctx.fingerprint_external_input(\n"
        "                ctx.config.annotation_cache, label=\"annotation_cache\"\n"
        "            )\n"
        "        )\n",
        label="CNV plan cache",
    )


def patch_target_coverage_extension() -> None:
    path = ROOT / "src/ontseq_platform/target_coverage_extension.py"
    replace_once(
        path,
        "from .pipeline.envelope import sha256_file\n",
        "",
        label="target coverage digest import",
    )
    replace_once(
        path,
        "        external_inputs=(\n"
        "            (target_bed.name, sha256_file(target_bed)),\n"
        "            (bam.name, sha256_file(bam)),\n"
        "        ),\n",
        "        external_inputs=(\n"
        "            ctx.fingerprint_external_input(target_bed),\n"
        "            ctx.fingerprint_external_input(bam),\n"
        "        ),\n",
        label="target coverage plan cache",
    )


def patch_methylation() -> None:
    path = ROOT / "src/ontseq_platform/methylation.py"
    replace_once(
        path,
        '_VERSION = re.compile(r"(?<!\\d)(\\d+\\.\\d+\\.\\d+)(?!\\d)")\n',
        '_VERSION = re.compile(r"(?<!\\d)(\\d+\\.\\d+\\.\\d+)(?!\\d)")\n'
        '_MODKIT_FAILED_PROCESSING = re.compile(r"(?i)(\\d+)\\s+failed processing\\b")\n'
        '_INDEPENDENT_CYTOSINE_GROUP_EXPR = (\n'
        '    \'[MM] =~ "C[+-][^;]*;.*C[+-][^;]*;"\'\n'
        ')\n',
        label="methylation safety constants",
    )

    anchor = "    return int(text)\n\n\ndef _modkit_include_bed"
    insertion = '''    return int(text)\n\n\ndef count_reads_with_independent_cytosine_mod_groups(\n    bam: Path,\n    *,\n    runner: CommandRunner,\n    samtools: str = "samtools",\n    threads: int = 4,\n) -> int | None:\n    """Count reads with two independent cytosine MM groups.\n\n    modkit 0.6.4 can silently advance same-base MM groups together and lose or\n    misassign calls. A single multi-code group is not matched; only two separate C\n    groups are treated as the known-risk representation.\n    """\n\n    result = runner.run(\n        [\n            samtools,\n            "view",\n            "-c",\n            "-@",\n            str(threads),\n            "-e",\n            _INDEPENDENT_CYTOSINE_GROUP_EXPR,\n            str(bam),\n        ],\n        timeout_seconds=7200,\n    )\n    if result.returncode != 0:\n        return None\n    text = result.stdout.strip()\n    if not text.isdigit():\n        return None\n    return int(text)\n\n\ndef _modkit_failed_processing_count(\n    log_path: Path, *diagnostics: str\n) -> int:\n    """Return the largest record-failure count reported by a nominally successful pileup."""\n\n    texts = [item for item in diagnostics if item]\n    if log_path.is_file():\n        texts.append(log_path.read_text(encoding="utf-8", errors="replace"))\n    counts = [\n        int(match)\n        for text in texts\n        for match in _MODKIT_FAILED_PROCESSING.findall(text)\n    ]\n    return max(counts, default=0)\n\n\ndef _modkit_include_bed'''
    replace_once(path, anchor, insertion, label="methylation safety helpers")

    replace_once(
        path,
        "    warnings: list[str] = []\n    tagged_reads: int | None = None\n",
        "    warnings: list[str] = []\n"
        "    tagged_reads: int | None = None\n"
        "    independent_cytosine_groups: int | None = None\n",
        label="methylation risk counter",
    )

    anchor = '''        if tagged_reads is None:\n            warnings.append(\n                "The installed samtools could not evaluate a tag filter expression, so the "\n                "presence of MM/ML tags was not verified before the pileup ran."\n            )\n\n    include_bed = ('''
    insertion = '''        if tagged_reads is None:\n            warnings.append(\n                "The installed samtools could not evaluate a tag filter expression, so the "\n                "presence of MM/ML tags was not verified before the pileup ran."\n            )\n\n    if version == "0.6.4":\n        independent_cytosine_groups = count_reads_with_independent_cytosine_mod_groups(\n            Path(manifest.input.path),\n            runner=command_runner,\n            samtools=samtools,\n            threads=threads,\n        )\n        if independent_cytosine_groups is None:\n            raise ValueError(\n                "modkit 0.6.4 has a known silent-error path for independent same-base MM "\n                "groups, but samtools could not verify whether this BAM contains that "\n                "representation. Refusing an unverified methylation pileup"\n            )\n        if independent_cytosine_groups > 0:\n            raise ValueError(\n                f"The aligned BAM contains {independent_cytosine_groups} read(s) with "\n                "independent cytosine MM groups. modkit 0.6.4 can silently lose or "\n                "misassign calls for this valid representation; refusing the pileup until "\n                "a corrected pinned modkit release is validated"\n            )\n\n    include_bed = ('''
    replace_once(path, anchor, insertion, label="methylation same-base guard")

    anchor = '''    result = command_runner.run(argv, timeout_seconds=14400)\n    if result.returncode != 0:\n        detail = (result.stderr or result.stdout or "").strip().splitlines()\n        tail = detail[-1] if detail else "no diagnostic output"\n        raise ValueError(f"modkit pileup failed with exit code {result.returncode}: {tail}")\n    if not bedmethyl_path.is_file():\n        raise ValueError("modkit pileup reported success but produced no bedMethyl output")\n'''
    insertion = '''    result = command_runner.run(argv, timeout_seconds=14400)\n    if result.returncode != 0:\n        detail = (result.stderr or result.stdout or "").strip().splitlines()\n        tail = detail[-1] if detail else "no diagnostic output"\n        raise ValueError(f"modkit pileup failed with exit code {result.returncode}: {tail}")\n    failed_processing = _modkit_failed_processing_count(\n        log_path, result.stdout, result.stderr\n    )\n    if failed_processing:\n        bedmethyl_path.unlink(missing_ok=True)\n        raise ValueError(\n            f"modkit pileup exited successfully but reported {failed_processing} failed "\n            "record(s); refusing partial bedMethyl output"\n        )\n    if not bedmethyl_path.is_file():\n        raise ValueError("modkit pileup reported success but produced no bedMethyl output")\n'''
    replace_once(path, anchor, insertion, label="modkit failed-processing guard")

    replace_once(
        path,
        '        "expected_version": policy.expected_version,\n',
        '        "expected_version": policy.expected_version,\n'
        '        "modkit_064_independent_cytosine_group_guard": independent_cytosine_groups,\n',
        label="methylation safety provenance",
    )


def patch_methylation_tests() -> None:
    path = ROOT / "tests/test_methylation.py"
    replace_once(
        path,
        '    def __init__(self, *, version: str = "0.6.4", tagged_reads: str = "1200") -> None:\n'
        '        self.version = version\n'
        '        self.tagged_reads = tagged_reads\n'
        '        self.commands: list[list[str]] = []\n',
        '    def __init__(\n'
        '        self,\n'
        '        *,\n'
        '        version: str = "0.6.4",\n'
        '        tagged_reads: str = "1200",\n'
        '        independent_cytosine_groups: str = "0",\n'
        '        failed_processing: int = 0,\n'
        '    ) -> None:\n'
        '        self.version = version\n'
        '        self.tagged_reads = tagged_reads\n'
        '        self.independent_cytosine_groups = independent_cytosine_groups\n'
        '        self.failed_processing = failed_processing\n'
        '        self.commands: list[list[str]] = []\n',
        label="fake runner state",
    )
    replace_once(
        path,
        '        if "view" in argv:\n'
        '            return CommandResult(tuple(argv), 0, self.tagged_reads, "")\n'
        '        if "pileup" in argv:\n'
        '            Path(argv[3]).write_text("\\n".join(self.rows) + "\\n", encoding="utf-8")\n'
        '            return CommandResult(tuple(argv), 0, "", "")\n',
        '        if "view" in argv:\n'
        '            count = (\n'
        '                self.independent_cytosine_groups\n'
        '                if any("C[+-]" in item for item in argv)\n'
        '                else self.tagged_reads\n'
        '            )\n'
        '            return CommandResult(tuple(argv), 0, count, "")\n'
        '        if "pileup" in argv:\n'
        '            Path(argv[3]).write_text("\\n".join(self.rows) + "\\n", encoding="utf-8")\n'
        '            if self.failed_processing:\n'
        '                log_path = Path(argv[argv.index("--log-filepath") + 1])\n'
        '                log_path.write_text(\n'
        '                    f"[ERROR] ~{self.failed_processing} failed processing\\n",\n'
        '                    encoding="utf-8",\n'
        '                )\n'
        '            return CommandResult(tuple(argv), 0, "", "")\n',
        label="fake runner behavior",
    )

    anchor = '''    def test_a_completed_pileup_is_normalized_and_records_the_tag_count(self) -> None:\n'''
    tests = '''    def test_modkit_064_refuses_independent_same_base_cytosine_groups(self) -> None:\n        with tempfile.TemporaryDirectory() as directory:\n            root = Path(directory)\n            manifest, intake = self._fixture(root)\n            runner = _FakeRunner(independent_cytosine_groups="2")\n            with self.assertRaises(ValueError) as raised:\n                run_methylation(\n                    manifest,\n                    intake,\n                    _policy(cpg_only=False, combine_strands=False),\n                    output_dir=root / "out",\n                    reference_fasta=root / "reference.fa",\n                    runner=runner,\n                )\n        self.assertIn("independent cytosine MM groups", str(raised.exception))\n\n    def test_success_exit_with_failed_processing_is_refused_and_output_removed(self) -> None:\n        with tempfile.TemporaryDirectory() as directory:\n            root = Path(directory)\n            manifest, intake = self._fixture(root)\n            runner = _FakeRunner(failed_processing=3)\n            with self.assertRaises(ValueError) as raised:\n                run_methylation(\n                    manifest,\n                    intake,\n                    _policy(cpg_only=False, combine_strands=False),\n                    output_dir=root / "out",\n                    reference_fasta=root / "reference.fa",\n                    runner=runner,\n                )\n            output = root / "out" / "SYNTHETIC_001.modkit.bedmethyl"\n            self.assertFalse(output.exists())\n        self.assertIn("reported 3 failed record", str(raised.exception))\n\n''' + anchor
    replace_once(path, anchor, tests, label="methylation safety unit tests")


def patch_methylation_integration_test() -> None:
    path = ROOT / "tests/test_methylation_pipeline_integration.py"
    anchor = '''def test_optional_methylation_reaches_schema_030_and_preserves_resource_and_iscn_contracts(\n'''
    test = '''def test_non_cpg_methylation_plan_still_requires_the_reference(tmp_path: Path) -> None:\n    fixture = _Fixture(tmp_path)\n    policy = fixture.policy.model_copy(\n        update={"cpg_only": False, "combine_strands": False}\n    )\n    config = replace(\n        fixture.config, methylation_policy=policy, reference_fasta=None\n    )\n    context = RunContext(\n        config=config,\n        envelope=fixture.envelope(),\n        runner=_ModkitVersionOnly(),\n        manifest=config.manifest,\n    )\n\n    with pytest.raises(StageFailure, match="reference FASTA"):\n        IMPLEMENTATIONS[StageId.METHYLATION].plan(context)\n\n\n''' + anchor
    replace_once(path, anchor, test, label="methylation planner regression")


def patch_real_modkit_test() -> None:
    path = ROOT / "tests/test_modkit_real_tool.py"
    anchor = '''    def test_real_bam_without_modification_tags_never_becomes_unmethylated(self) -> None:\n'''
    test = '''    def test_real_samtools_detects_independent_same_base_groups_before_modkit(self) -> None:\n        reference = self.root / "independent.fa"\n        reference.write_text(">chr1\\nCC\\n", encoding="utf-8")\n        self.pysam.faidx(str(reference))\n        bam = self.root / "independent.bam"\n        header = {\n            "HD": {"VN": "1.6", "SO": "coordinate"},\n            "SQ": [{"SN": "chr1", "LN": 2}],\n        }\n        with self.pysam.AlignmentFile(str(bam), "wb", header=header) as target:\n            record = self.pysam.AlignedSegment()\n            record.query_name = "INDEPENDENT_GROUPS"\n            record.query_sequence = "CC"\n            record.query_qualities = array.array("B", [40, 40])\n            record.flag = 0\n            record.reference_id = 0\n            record.reference_start = 0\n            record.mapping_quality = 60\n            record.cigarstring = "2M"\n            record.set_tag("MN", 2)\n            record.set_tag("MM", "C+m?,0;C+h?,1;")\n            record.set_tag("ML", array.array("B", [200, 250]))\n            target.write(record)\n        self.pysam.index(str(bam))\n        manifest = SampleManifest(\n            sample_id="INDEPENDENT_MODKIT",\n            run_id="SYNTHETIC_RUN",\n            input=InputSpec(\n                kind=InputKind.ALIGNED_BAM,\n                path=str(bam),\n                index_path=str(bam) + ".bai",\n            ),\n            assay=AssaySpec(\n                mode=AssayMode.LOW_COVERAGE_WGS,\n                genome_build=GenomeBuild.GRCH38,\n                reference_id="synthetic-independent-reference",\n            ),\n            analysis=AnalysisSpec(\n                profile="synthetic", modules=[AnalysisModule.METHYLATION]\n            ),\n        )\n        intake = AlignedBamIntakeReport(\n            sample_id=manifest.sample_id,\n            reference_id=manifest.assay.reference_id,\n            genome_build=GenomeBuild.GRCH38,\n            checks=[],\n            verdict=Verdict.PASS,\n        )\n        policy = MethylationPolicy(\n            profile_id="SYNTHETIC_INDEPENDENT_GROUPS",\n            expected_version="0.6.4",\n            status="technical_defaults_only",\n            modification_codes=["m", "h"],\n            cpg_only=False,\n            combine_strands=False,\n            region_source="chromosome",\n            note="Regression for the modkit 0.6.4 independent-group safety gate",\n        )\n\n        with self.assertRaisesRegex(ValueError, "independent cytosine MM groups"):\n            run_methylation(\n                manifest,\n                intake,\n                policy,\n                output_dir=self.root / "independent-output",\n                reference_fasta=reference,\n                threads=1,\n            )\n        self.assertFalse(\n            (\n                self.root\n                / "independent-output"\n                / "INDEPENDENT_MODKIT.modkit.bedmethyl"\n            ).exists()\n        )\n\n''' + anchor
    replace_once(path, anchor, test, label="real modkit risk regression")


def patch_docs() -> None:
    path = ROOT / "docs/METHYLATION_LANE.md"
    replace_once(
        path,
        "## Three refusals\n",
        "## Five refusals\n",
        label="methylation refusal heading",
    )
    anchor = '''**An empty pileup is never reported as "unmethylated".**'''
    addition = '''**Known modkit 0.6.4 same-base MM-group risk fails closed.** The pinned\nrelease can silently lose or misassign calls when one read carries independent modification\ngroups for the same canonical cytosine (for example separate `C+m` and `C+h` groups). ONTSeq\nchecks the complete BAM with the pinned samtools filter-expression engine before pileup and\nrefuses this representation. A single multi-code group is not rejected. This is an upstream\ntool limitation, not a methylation threshold.\n\n**A nominally successful pileup with failed records is not accepted.** modkit 0.6.4 can exit\nzero while its log reports records that failed processing. ONTSeq inspects that diagnostic;\nany non-zero `failed processing` count deletes the partial bedMethyl result and fails the\nstage.\n\n''' + anchor
    replace_once(path, anchor, addition, label="methylation upstream safety docs")

    replace_once(
        path,
        "- The real-binary CI lane proves the adapter meets the pinned modkit 0.6.4 on synthetic\n"
        "  fixtures. It proves nothing about analytical recovery on biological data.\n",
        "- The real-binary CI lane proves selected interoperability cases against pinned modkit\n"
        "  0.6.4. It does not establish correctness for every valid MM/ML representation, and\n"
        "  current upstream 0.6.4 limitations require the explicit same-base-group and failed-\n"
        "  processing refusals above. It proves nothing about analytical recovery on biological\n"
        "  data.\n",
        label="methylation limitations docs",
    )

    changelog = ROOT / "CHANGELOG.md"
    anchor = "## Unreleased\n\n### Changed\n"
    addition = '''## Unreleased\n\n### Fixed\n\n- Reuse stable external-input SHA-256 values within one run invocation instead of reading\n  the same multi-gigabyte BAM once per stage plan. The cache is run-scoped and keyed by\n  resolved path plus filesystem identity/state; unstable reads are never cached. Intake\n  mutation detection and release/artifact verification remain uncached and still re-read\n  bytes deliberately. CNV and Adaptive-Sampling coverage extension plans share the same\n  run-scoped cache.\n- Make methylation planning match execution after the modkit 0.6.4 migration: every\n  `--modified-bases` run requires the locked reference FASTA, including non-CpG policies.\n- Fail closed around two modkit 0.6.4 scientific-correctness hazards: independent same-base\n  cytosine MM groups are rejected before pileup, and a zero-exit pileup that reports any\n  `failed processing` records has its partial bedMethyl output discarded. These guards do\n  not validate methylation accuracy; they prevent known silent-tool failure modes from\n  becoming plausible normalized fractions.\n\n### Validation impact\n\n- No CNV/SV/methylation threshold, caller score, reportability boundary or biological\n  classifier changes. Hashing changes only duplicate I/O. Methylation acceptance becomes\n  stricter for inputs exposed to known modkit 0.6.4 failure modes.\n\n### Changed\n'''
    replace_once(changelog, anchor, addition, label="audit changelog")


def main() -> None:
    write_input_digest_helper()
    write_input_digest_tests()
    patch_runner()
    patch_cnv_extension()
    patch_target_coverage_extension()
    patch_methylation()
    patch_methylation_tests()
    patch_methylation_integration_test()
    patch_real_modkit_test()
    patch_docs()


if __name__ == "__main__":
    main()
