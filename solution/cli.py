"""Command-line entry point for GC sequence review."""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

from .calibration import weighted_calibration
from .chromatogram import ANALYTE_WINDOW, INTERNAL_STANDARD_WINDOW, integrate_peak
from .reporting import review_samples


SEQUENCE_COLUMNS = (
    "run_id",
    "sample_id",
    "role",
    "replicate",
    "dilution_factor",
    "nominal_mg_l",
    "spike_added_mg_l",
)
CHROMATOGRAM_COLUMNS = ("run_id", "time_s", "signal_mv")
VALID_ROLES = {"calibration", "blank", "sample", "spike"}


def _number(value: str, field: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} is not numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"{field} must be finite")
    return number


def _optional_number(value: str, field: str) -> float | None:
    if value.strip() == "":
        return None
    return _number(value, field)


def _read_sequence(path: str) -> list[dict[str, object]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != SEQUENCE_COLUMNS:
            raise ValueError("sequence columns do not match the required schema")
        rows: list[dict[str, object]] = []
        seen_run_ids: set[str] = set()
        for raw in reader:
            run_id = (raw["run_id"] or "").strip()
            sample_id = (raw["sample_id"] or "").strip()
            role = (raw["role"] or "").strip()
            if not run_id or not sample_id:
                raise ValueError("run_id and sample_id must be non-empty")
            if run_id in seen_run_ids:
                raise ValueError(f"duplicate run_id: {run_id}")
            seen_run_ids.add(run_id)
            if role not in VALID_ROLES:
                raise ValueError(f"unknown role: {role}")
            replicate = _number(raw["replicate"], "replicate")
            if not replicate.is_integer() or replicate < 1:
                raise ValueError("replicate must be a positive integer")
            dilution = _number(raw["dilution_factor"], "dilution_factor")
            if dilution <= 0:
                raise ValueError("dilution_factor must be positive")
            nominal = _optional_number(raw["nominal_mg_l"], "nominal_mg_l")
            spike_added = _optional_number(
                raw["spike_added_mg_l"], "spike_added_mg_l"
            )
            if role == "calibration":
                if nominal is None or nominal <= 0 or spike_added is not None:
                    raise ValueError("calibration rows require only positive nominal_mg_l")
            elif role == "spike":
                if spike_added is None or spike_added <= 0:
                    raise ValueError("spike rows require positive spike_added_mg_l")
            else:
                if nominal is not None or spike_added is not None:
                    raise ValueError("unused concentration fields must be empty")
            rows.append(
                {
                    "run_id": run_id,
                    "sample_id": sample_id,
                    "role": role,
                    "replicate": int(replicate),
                    "dilution_factor": dilution,
                    "nominal_mg_l": nominal,
                    "spike_added_mg_l": spike_added,
                }
            )
    if not rows:
        raise ValueError("sequence is empty")
    return rows


def _read_chromatograms(path: str) -> dict[str, tuple[list[float], list[float]]]:
    times: dict[str, list[float]] = defaultdict(list)
    signals: dict[str, list[float]] = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != CHROMATOGRAM_COLUMNS:
            raise ValueError("chromatogram columns do not match the required schema")
        for raw in reader:
            run_id = (raw["run_id"] or "").strip()
            if not run_id:
                raise ValueError("chromatogram run_id must be non-empty")
            times[run_id].append(_number(raw["time_s"], "time_s"))
            signals[run_id].append(_number(raw["signal_mv"], "signal_mv"))
    for run_id, run_times in times.items():
        if any(right <= left for left, right in zip(run_times, run_times[1:])):
            raise ValueError(f"time values are not strictly increasing for {run_id}")
        if len(run_times) != len(signals[run_id]):
            raise ValueError(f"incomplete chromatogram for {run_id}")
    return {run_id: (times[run_id], signals[run_id]) for run_id in times}


def _validate_sequence_groups(sequence: list[dict[str, object]]) -> None:
    calibration_ids = {
        str(row["sample_id"]) for row in sequence if row["role"] == "calibration"
    }
    if len(calibration_ids) < 5:
        raise ValueError("at least five calibration levels are required")
    sample_ids = {
        str(row["sample_id"]) for row in sequence if row["role"] == "sample"
    }
    if not sample_ids:
        raise ValueError("at least one sample is required")
    for sample_id in sample_ids:
        sample_rows = [
            row
            for row in sequence
            if row["role"] == "sample" and str(row["sample_id"]) == sample_id
        ]
        spike_rows = [
            row
            for row in sequence
            if row["role"] == "spike" and str(row["sample_id"]) == sample_id
        ]
        if len(sample_rows) != 2 or len(spike_rows) != 1:
            raise ValueError(
                f"{sample_id} requires exactly two sample rows and one spike row"
            )
    if not any(row["role"] == "blank" for row in sequence):
        raise ValueError("at least one method blank is required")


def build_report(
    sequence: list[dict[str, object]],
    chromatograms: dict[str, tuple[list[float], list[float]]],
) -> dict[str, object]:
    _validate_sequence_groups(sequence)
    sequence_ids = {str(row["run_id"]) for row in sequence}
    if sequence_ids != set(chromatograms):
        missing = sorted(sequence_ids - set(chromatograms))
        extra = sorted(set(chromatograms) - sequence_ids)
        raise ValueError(f"run_id mismatch; missing={missing}, extra={extra}")

    ratios: dict[str, float] = {}
    for row in sequence:
        run_id = str(row["run_id"])
        times, signals = chromatograms[run_id]
        analyte = integrate_peak(times, signals, ANALYTE_WINDOW)
        internal = integrate_peak(times, signals, INTERNAL_STANDARD_WINDOW)
        internal_area = float(internal["area"])
        if internal_area <= 0:
            raise ValueError(f"non-positive internal-standard area for {run_id}")
        ratios[run_id] = float(analyte["area"]) / internal_area

    calibration_rows = [
        row for row in sequence if row["role"] == "calibration"
    ]
    calibration = weighted_calibration(
        [
            {
                "nominal_mg_l": float(row["nominal_mg_l"]),
                "response_ratio": ratios[str(row["run_id"])],
            }
            for row in calibration_rows
        ]
    )
    review = review_samples(sequence, ratios, calibration)
    return {
        "calibration": calibration,
        "blank_mean_mg_l": review["blank_mean_mg_l"],
        "blank_correction_mg_l": review["blank_correction_mg_l"],
        "qc": review["qc"],
        "samples": review["samples"],
        "overall_status": review["overall_status"],
    }


def run(sequence_path: str, chromatogram_path: str, output_path: str) -> None:
    report = build_report(
        _read_sequence(sequence_path),
        _read_chromatograms(chromatogram_path),
    )
    Path(output_path).write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    arguments = sys.argv[1:] if argv is None else argv
    if len(arguments) != 3:
        print(
            "usage: python -m gc_review.cli SEQUENCE CHROMATOGRAMS OUTPUT",
            file=sys.stderr,
        )
        return 2
    try:
        run(*arguments)
    except (OSError, ValueError, KeyError, ZeroDivisionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
