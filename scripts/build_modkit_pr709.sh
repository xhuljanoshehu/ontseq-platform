#!/usr/bin/env bash
# Build a reviewable candidate. This script does not modify the runtime allowlist.
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
if [[ $# != 1 ]]; then
  echo 'Usage: scripts/build_modkit_pr709.sh /absolute/build-directory' >&2
  exit 2
fi
build_dir="$1"
[[ "$build_dir" = /* ]] || { echo 'Use an absolute build directory' >&2; exit 2; }
mkdir -p "$build_dir"
source_dir="$build_dir/source"
source_commit=9fb9aea763aa1b78ac1172fa6d6734e7955ccc70
[[ "$(rustc --version | cut -d ' ' -f 2)" = 1.93.1 ]] || {
  echo 'This build recipe requires Rust 1.93.1' >&2; exit 1;
}
if [[ ! -d "$source_dir" ]]; then
  git clone --filter=blob:none --no-checkout \
    https://github.com/SuhasSrinivasan/modkit.git "$source_dir"
  git -C "$source_dir" sparse-checkout init --cone
  git -C "$source_dir" sparse-checkout set \
    modkit/src modkit-core/src modkit-logging/src ochm/src safe-record/src
  git -C "$source_dir" checkout --detach "$source_commit"
fi
[[ "$(git -C "$source_dir" rev-parse HEAD)" = "$source_commit" ]]
git -C "$source_dir" diff --exit-code HEAD
python - "$source_dir" <<'PY_SOURCE'
import pathlib, subprocess, sys
root = pathlib.Path(sys.argv[1])
extras = set()
for mode in ([], ['--ignored']):
    result = subprocess.check_output(['git', '-C', str(root), 'ls-files', '--others',
                                      '--exclude-standard', '-z', *mode])
    extras.update(p.decode() for p in result.split(b'\0') if p)
extras.discard('Cargo.lock')  # Replaced with the separately checksummed lock below.
if extras:
    raise SystemExit('Unexpected source/build inputs: ' + ', '.join(sorted(extras)))
PY_SOURCE
python - "$repo_dir/build-support/modkit-pr709" <<'PY'
import hashlib, json, pathlib, sys
root = pathlib.Path(sys.argv[1])
record = json.loads((root / 'source.json').read_text())
assert hashlib.sha256((root / 'Cargo.lock').read_bytes()).hexdigest() == record['cargo_lock_sha256']
PY
cp "$repo_dir/build-support/modkit-pr709/Cargo.lock" "$source_dir/Cargo.lock"
export CARGO_TARGET_DIR="$build_dir/target"
(
  cd "$source_dir"
  cargo build --locked --release -p modkit
)
binary="$CARGO_TARGET_DIR/release/modkit"
cd "$repo_dir"
PATH="$(dirname "$binary"):$PATH" \
PYTHONPATH="$repo_dir/src${PYTHONPATH:+:$PYTHONPATH}" \
PYTEST_ADDOPTS= ONTSEQ_MODKIT_REAL_TOOL=1 ONTSEQ_MODKIT_PR709_REAL_TOOL=1 \
ONTSEQ_MODKIT_PR709_BINARY="$binary" \
python -m pytest -o addopts= -q tests/test_modkit_real_tool.py \
  tests/test_modkit_pr709_real.py --junitxml="$build_dir/qualification.xml"
python - "$build_dir/qualification.xml" <<'PY_TESTS'
import sys, xml.etree.ElementTree as ET
cases = list(ET.parse(sys.argv[1]).getroot().iter('testcase'))
assert len(cases) == 76, f'Expected all 76 qualification cases, observed {len(cases)}'
assert not any(case.find(tag) is not None for case in cases
               for tag in ('failure', 'error', 'skipped')), 'Qualification must contain only passing cases'
PY_TESTS
sha256sum "$binary" "$source_dir/Cargo.lock" > "$build_dir/SHA256SUMS"
cp "$source_dir/LICENCE.txt" "$build_dir/LICENCE.modkit.txt"
printf '%s\n' 'Candidate built and tested. Register its checksum through source review before installation.'
