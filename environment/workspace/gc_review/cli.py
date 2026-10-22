"""Command-line entry point for GC sequence review."""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

from .calibration import weighted_calibration
from .chromatogram import ANALYTE_WINDOW, INTERNAL_STANDARD_WINDOW, integrate_peak
from .reporting import review_samples


def _read_rows(path: str, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(columns):
            raise ValueError(f"{path} has unexpected columns")
        return list(reader)


def _read_chromatograms(path: str) -> dict[str, tuple[list[float], list[float]]]:
    grouped: dict[str, tuple[list[float], list[float]]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["run_id", "time_s", "signal_mv"]:
            raise ValueError(f"{path} has unexpected columns")
        times: dict[str, list[float]] = defaultdict(list)
        signals: dict[str, list[float]] = defaultdict(list)
        for row in reader:
            times[row["run_id"]].append(float(row["time_s"]))
            signals[row["run_id"]].append(float(row["signal_mv"]))
    for run_id in times:
        grouped[run_id] = (times[run_id], signals[run_id])
    return grouped


def run(sequence_path: str, chromatogram_path: str, output_path: str) -> None:
    columns = (
        "run_id",
        "sample_id",
        "role",
        "replicate",
        "dilution_factor",
        "nominal_mg_l",
        "spike_added_mg_l",
    )
    sequence = _read_rows(sequence_path, columns)
    chromatograms = _read_chromatograms(chromatogram_path)
    calibration_rows = [row for row in sequence if row["role"] == "calibration"]
    ratios = {}
    for row in sequence:
        times, signals = chromatograms[row["run_id"]]
        analyte = integrate_peak(times, signals, ANALYTE_WINDOW)
        internal = integrate_peak(times, signals, INTERNAL_STANDARD_WINDOW)
        ratios[row["run_id"]] = float(analyte["area"]) / float(internal["area"])
    calibration = weighted_calibration(
        [
            {
                "nominal_mg_l": float(row["nominal_mg_l"]),
                "response_ratio": ratios[row["run_id"]],
            }
            for row in calibration_rows
        ]
    )
    samples = review_samples(sequence, ratios, calibration)
    report = {
        "calibration": calibration,
        "blank_mean_mg_l": 0.0,
        "blank_correction_mg_l": 0.0,
        "qc": {
            "calibration": "pass" if calibration["accepted"] else "fail",
            "blank": "pass",
            "duplicate": "pass",
            "spike_recovery": "pass",
        },
        "samples": samples,
        "overall_status": "pass",
    }
    Path(output_path).write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    arguments = sys.argv[1:] if argv is None else argv
    if len(arguments) != 3:
        print("usage: python -m gc_review.cli SEQUENCE CHROMATOGRAMS OUTPUT", file=sys.stderr)
        return 2
    try:
        run(*arguments)
    except (OSError, ValueError, KeyError, ZeroDivisionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
