from __future__ import annotations

import sys
from pathlib import Path

from ..io import load_model
from ..runtime_cli import main as runtime_main
from .extension import QDNAseqExtensionSettings, register_qdnaseq_extension
from .qdnaseq import QDNAseqPolicy


def _take_option(name: str, default: str) -> str:
    try:
        index = sys.argv.index(name)
    except ValueError:
        return default
    if index + 1 >= len(sys.argv):
        raise SystemExit(f"ERROR: {name} requires a value")
    value = sys.argv[index + 1]
    del sys.argv[index : index + 2]
    return value


def main() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    policy_path = Path(
        _take_option(
            "--cnv-policy",
            str(repo_root / "configs/cnv/qdnaseq_ace.technical.yaml"),
        )
    )
    rscript = _take_option("--qdnaseq-rscript", "Rscript")
    script = Path(
        _take_option(
            "--qdnaseq-script",
            str(repo_root / "scripts/run_qdnaseq_ace.R"),
        )
    )
    if policy_path.is_file():
        policy = load_model(policy_path, QDNAseqPolicy)
    else:
        policy = QDNAseqPolicy(
            profile_id="qdnaseq-ace-multibin-v1",
            note="Built-in fallback matching configs/cnv/qdnaseq_ace.technical.yaml",
        )
    register_qdnaseq_extension(
        QDNAseqExtensionSettings(policy=policy, rscript=rscript, script=script)
    )
    runtime_main()


if __name__ == "__main__":
    main()
