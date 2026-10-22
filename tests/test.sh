#!/bin/bash
set -u

reward_file=/logs/verifier/reward.txt
ctrf_file=/logs/verifier/ctrf.json
mkdir -p /logs/verifier
echo "0" > "$reward_file"

write_reward_zero() {
  echo "0" > "$reward_file"
}

trap 'write_reward_zero; exit 130' INT TERM HUP

if ! python -m pip install --user --no-cache-dir --disable-pip-version-check \
  --require-hashes --only-binary=:all: -r /tests/requirements.txt; then
  write_reward_zero
  exit 0
fi

cd /workspace
python -m pytest --ctrf "$ctrf_file" /tests/test_outputs.py -q

python3 - "$ctrf_file" "$reward_file" <<'PY'
import json
import sys

ctrf_path, reward_path = sys.argv[1:]
groups = {
    "C1": (0.30, "test_c1_"),
    "C2": (0.30, "test_c2_"),
    "C3": (0.30, "test_c3_"),
    "E2E": (0.10, "test_e2e_"),
}
try:
    with open(ctrf_path, encoding="utf-8") as handle:
        tests = json.load(handle)["results"]["tests"]
except (OSError, KeyError, TypeError, json.JSONDecodeError):
    tests = []

statuses = {}
for test in tests:
    name = str(test.get("name", "")).split("::")[-1]
    statuses[name] = test.get("status") == "passed"

score = 0.0
for weight, prefix in groups.values():
    relevant = [passed for name, passed in statuses.items() if name.startswith(prefix)]
    if relevant and all(relevant):
        score += weight

score = min(1.0, max(0.0, score))
with open(reward_path, "w", encoding="utf-8") as handle:
    handle.write(f"{score:.6f}\n")
PY

trap - INT TERM HUP
exit 0
