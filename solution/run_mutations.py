#!/usr/bin/env python3
"""Run per-mutation counter-evidence for the GC review task."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


TASK = Path(__file__).resolve().parents[1]
EVIDENCE = TASK / "solution" / "evidence"
MODULES = ["chromatogram", "calibration", "reporting", "cli"]

MUTATIONS = [
    {
        "id": "M1_no_baseline",
        "module": "chromatogram",
        "old": "signal - _linear_baseline(time, start, end)",
        "new": "signal",
        "expected": [
            "test_c1_linear_baseline_is_removed",
            "test_c1_trapezoid_peak_area_matches_reference",
        ],
        "description": "跳过线性基线扣除，直接用原始信号积分",
    },
    {
        "id": "M2_unweighted_fit",
        "module": "calibration",
        "old": "1.0 / (x * x) for x in xs",
        "new": "1.0 for _ in xs",
        "expected": ["test_c2_weighted_least_squares_parameters"],
        "description": "把 1/x^2 权重退化成等权重",
    },
    {
        "id": "M3_wrong_loq",
        "module": "calibration",
        "old": "loq = 10.0 * sdyx / slope",
        "new": "loq = 3.3 * sdyx / slope",
        "expected": ["test_c2_lod_and_loq_use_weighted_residual_error"],
        "description": "把 LOQ 系数从 10 改成 3.3",
    },
    {
        "id": "M4_no_blank_correction",
        "module": "reporting",
        "old": "raw * dilution_factor - blank_correction_mg_l",
        "new": "raw * dilution_factor",
        "expected": ["test_c3_dilution_and_blank_correction_are_applied"],
        "description": "保留稀释倍数但跳过空白扣除",
    },
    {
        "id": "M5_wrong_rpd_denominator",
        "module": "reporting",
        "old": "mean = 0.5 * (first + second)",
        "new": "mean = first + second",
        "expected": ["test_c3_duplicate_rpd_boundary"],
        "description": "RPD 分母从平行样均值改成总和",
    },
]


def build_image() -> str:
    image = f"harbor-gc-review-mutation:{os.getpid()}"
    subprocess.run(
        ["docker", "build", "-q", "-t", image, str(TASK / "environment")],
        check=True,
        capture_output=True,
    )
    return image


def prepare(name: str) -> Path:
    workspace = Path(tempfile.mkdtemp(prefix=f"gc-review-mut-{name}-"))
    shutil.copytree(
        TASK / "environment" / "workspace", workspace, dirs_exist_ok=True
    )
    for module in MODULES:
        shutil.copy(
            TASK / "solution" / f"{module}.py",
            workspace / "gc_review" / f"{module}.py",
        )
    return workspace


def apply_mutation(workspace: Path, mutation: dict) -> None:
    path = workspace / "gc_review" / f"{mutation['module']}.py"
    text = path.read_text(encoding="utf-8")
    if mutation["old"] not in text:
        raise SystemExit(f"anchor not found for {mutation['id']}: {mutation['old']!r}")
    path.write_text(
        text.replace(mutation["old"], mutation["new"], 1), encoding="utf-8"
    )


def run_tests(image: str, workspace: Path, logs: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{workspace}:/workspace",
            "-v",
            f"{TASK / 'tests'}:/tests:ro",
            "-v",
            f"{logs}:/logs",
            image,
            "bash",
            "/tests/test.sh",
        ],
        capture_output=True,
        text=True,
    )


def main() -> int:
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    image = build_image()
    results = []
    try:
        for mutation in MUTATIONS:
            workspace = prepare(mutation["id"])
            apply_mutation(workspace, mutation)
            logs = Path(tempfile.mkdtemp(prefix=f"gc-review-mut-logs-{mutation['id']}-"))
            os.chmod(logs, 0o777)
            process = run_tests(image, workspace, logs)
            log_path = EVIDENCE / f"mutation-{mutation['id']}.log"
            log_path.write_text(
                process.stdout + "\n--- stderr ---\n" + process.stderr,
                encoding="utf-8",
            )
            ctrf = logs / "verifier" / "ctrf.json"
            tests = (
                json.loads(ctrf.read_text(encoding="utf-8"))["results"]["tests"]
                if ctrf.exists()
                else []
            )
            failed = [
                str(test.get("name", "")).split("::")[-1]
                for test in tests
                if test.get("status") != "passed"
            ]
            killed = all(name in failed for name in mutation["expected"])
            results.append(
                {
                    "id": mutation["id"],
                    "description": mutation["description"],
                    "expected_tests": mutation["expected"],
                    "failed_tests": failed,
                    "killed": killed,
                    "evidence": str(log_path.relative_to(TASK)),
                }
            )
            print(f"{mutation['id']}: killed={killed} failed={failed}")
    finally:
        subprocess.run(
            ["docker", "image", "rm", "-f", image], capture_output=True, check=False
        )
    (TASK / "solution" / "mutation-results.json").write_text(
        json.dumps({"mutations": results}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return 0 if all(item["killed"] for item in results) else 1


if __name__ == "__main__":
    sys.exit(main())
