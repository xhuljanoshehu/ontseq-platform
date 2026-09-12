"""Record exact local source content without inventing a release git commit."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "build-metadata" / "SOURCE-TREE.json"


def main() -> None:
    raw = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=ROOT
    )
    names = sorted(set(raw.decode("utf-8").strip("\0").split("\0")))
    inventory = []
    for name in names:
        path = ROOT / name
        if path == OUTPUT or not path.is_file():
            continue
        inventory.append({"path": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    canonical = "".join(f"{row['sha256']}  {row['path']}\n" for row in inventory)
    base_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    output = {
        "schema_version": "1.0.0",
        "kind": "local_engineering_source_inventory_not_a_release_commit",
        "generated_at": datetime.now(UTC).isoformat(),
        "base_commit": base_commit,
        "working_tree_sha256": hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        "hash_format": "sha256 over sorted '<file_sha256>  <relative_path>\\n' entries",
        "excluded": ["build-metadata/SOURCE-TREE.json", "git-ignored build outputs"],
        "files": inventory,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"Source inventory: {len(inventory)} files, SHA256 {output['working_tree_sha256']}")


if __name__ == "__main__":
    main()
