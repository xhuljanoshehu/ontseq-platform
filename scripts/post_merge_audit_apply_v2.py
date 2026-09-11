from __future__ import annotations

import re
from pathlib import Path

import post_merge_audit_apply as audit


ROOT = Path(__file__).resolve().parents[1]


def patch_runner() -> None:
    path = ROOT / "src/ontseq_platform/pipeline/runner.py"
    audit.replace_once(
        path,
        "from .envelope import Artifact, RunEnvelope, sha256_file, stage_signature\n",
        "from .envelope import Artifact, RunEnvelope, sha256_file, stage_signature\n"
        "from .input_digest import RunInputDigestCache\n",
        label="runner import",
    )
    audit.replace_once(
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
    audit.replace_once(
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

    audit.replace_once(
        path,
        "    if policy.cpg_only and ctx.config.reference_fasta is None:\n"
        "        raise StageFailure(\n"
        "            \"the methylation policy restricts the pileup to CpG sites, which is a property \"\n"
        "            \"of the reference; pass --reference-fasta\"\n"
        "        )\n",
        "    if ctx.config.reference_fasta is None:\n"
        "        raise StageFailure(\n"
        "            \"modkit --modified-bases requires the locked reference FASTA for every \"\n"
        "            \"methylation run\"\n"
        "        )\n",
        label="methylation plan reference gate",
    )


def main() -> None:
    audit.write_input_digest_helper()
    audit.write_input_digest_tests()
    patch_runner()
    audit.patch_cnv_extension()
    audit.patch_target_coverage_extension()
    audit.patch_methylation()
    audit.patch_methylation_tests()
    audit.patch_methylation_integration_test()
    audit.patch_real_modkit_test()
    audit.patch_docs()


if __name__ == "__main__":
    main()
