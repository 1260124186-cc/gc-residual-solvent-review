#!/bin/bash
set -euo pipefail

task_dir="$(cd "$(dirname "$0")/.." && pwd)"
evidence_dir="$task_dir/solution/evidence"
image="harbor-gc-review-quality:$$"
mkdir -p "$evidence_dir"

docker build -q -t "$image" "$task_dir/environment" >/dev/null

prepare_workspace() {
  local name="$1"
  shift
  local directory
  directory="$(mktemp -d "${TMPDIR:-/tmp}/gc-review-${name}.XXXXXX")"
  cp -R "$task_dir/environment/workspace/." "$directory"/
  for module in "$@"; do
    cp "$task_dir/solution/$module.py" "$directory/gc_review/$module.py"
  done
  printf '%s\n' "$directory"
}

run_case() {
  local name="$1"
  local workspace="$2"
  local logs
  logs="$(mktemp -d "${TMPDIR:-/tmp}/gc-review-logs-${name}.XXXXXX")"
  chmod 0777 "$logs"
  {
    printf '$ docker run --rm -v %s:/workspace -v %s/tests:/tests:ro -v %s:/logs %s bash /tests/test.sh\n' \
      "$workspace" "$task_dir" "$logs" "$image"
    docker run --rm \
      -v "$workspace:/workspace" \
      -v "$task_dir/tests:/tests:ro" \
      -v "$logs:/logs" \
      "$image" \
      bash /tests/test.sh
    printf 'reward=%s\n' "$(cat "$logs/verifier/reward.txt")"
    printf 'summary='
    python3 - "$logs/verifier/ctrf.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    summary = json.load(handle)["results"]["summary"]
print(
    f"tests={summary['tests']} passed={summary['passed']} "
    f"failed={summary['failed']}"
)
PY
  } >"$evidence_dir/$name.log" 2>&1
}

baseline_workspace="$(mktemp -d "${TMPDIR:-/tmp}/gc-review-baseline.XXXXXX")"
cp -R "$task_dir/environment/workspace/." "$baseline_workspace"/
run_case baseline "$baseline_workspace"

oracle_workspace="$(prepare_workspace oracle chromatogram calibration reporting cli)"
run_case oracle "$oracle_workspace"

partial_c1_workspace="$(prepare_workspace c1 chromatogram)"
run_case partial-c1 "$partial_c1_workspace"

partial_c2_workspace="$(prepare_workspace c2 calibration)"
run_case partial-c2 "$partial_c2_workspace"

partial_c3_workspace="$(prepare_workspace c3 reporting)"
run_case partial-c3 "$partial_c3_workspace"

empty_workspace="$(prepare_workspace empty chromatogram calibration reporting)"
cp "$task_dir/solution/bypasses/empty/cli.py" "$empty_workspace/gc_review/cli.py"
run_case bypass-empty "$empty_workspace"

hardcoded_workspace="$(prepare_workspace hardcoded chromatogram calibration reporting)"
cp "$task_dir/solution/bypasses/hardcoded/cli.py" "$hardcoded_workspace/gc_review/cli.py"
run_case bypass-hardcoded "$hardcoded_workspace"

malformed_workspace="$(prepare_workspace malformed chromatogram calibration reporting)"
cp "$task_dir/solution/bypasses/malformed/cli.py" "$malformed_workspace/gc_review/cli.py"
run_case bypass-malformed "$malformed_workspace"

docker image rm -f "$image" >/dev/null 2>&1 || true
