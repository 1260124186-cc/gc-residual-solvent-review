#!/usr/bin/env python3
"""Generate the fixed-seed simulated GC sequence shipped with the task."""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path


SEED = 20261022
START_S = 36.0
STOP_S = 62.5
STEP_S = 0.5
ANALYTE_CENTER_S = 42.0
ANALYTE_SIGMA_S = 0.58
INTERNAL_CENTER_S = 57.0
INTERNAL_SIGMA_S = 0.66
INTERNAL_AMPLITUDE = 420.0


RUNS = [
    ("CAL-1", "CAL1", "calibration", 1, 1.0, 1.0, None, 0.73),
    ("CAL-2", "CAL2", "calibration", 1, 1.0, 2.0, None, 1.43),
    ("CAL-3", "CAL3", "calibration", 1, 1.0, 4.0, None, 2.83),
    ("CAL-4", "CAL4", "calibration", 1, 1.0, 8.0, None, 5.63),
    ("CAL-5", "CAL5", "calibration", 1, 1.0, 16.0, None, 11.23),
    ("BLK-1", "BLANK", "blank", 1, 1.0, None, None, 0.0),
    ("BLK-2", "BLANK", "blank", 2, 1.0, None, None, 0.0),
    ("S1-R1", "S1", "sample", 1, 5.0, None, None, 0.870),
    ("S1-R2", "S1", "sample", 2, 5.0, None, None, 0.887),
    ("S1-S", "S1", "spike", 1, 5.0, None, 8.0, 1.993),
    ("S2-R1", "S2", "sample", 1, 10.0, None, None, 0.0),
    ("S2-R2", "S2", "sample", 2, 10.0, None, None, 0.0),
    ("S2-S", "S2", "spike", 1, 10.0, None, 5.0, 0.380004622),
]


def gaussian(time_s: float, center_s: float, sigma_s: float, amplitude: float) -> float:
    distance = (time_s - center_s) / sigma_s
    return amplitude * math.exp(-0.5 * distance * distance)


def write_sequence(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "run_id",
                "sample_id",
                "role",
                "replicate",
                "dilution_factor",
                "nominal_mg_l",
                "spike_added_mg_l",
            ]
        )
        for row in RUNS:
            writer.writerow(["" if value is None else value for value in row[:7]])


def write_chromatograms(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["run_id", "time_s", "signal_mv"])
        for run in RUNS:
            run_id, target_ratio = run[0], float(run[7])
            is_amplitude = INTERNAL_AMPLITUDE + rng.uniform(-1.5, 1.5)
            analyte_amplitude = (
                target_ratio
                * is_amplitude
                * INTERNAL_SIGMA_S
                / ANALYTE_SIGMA_S
            )
            steps = int(round((STOP_S - START_S) / STEP_S))
            for index in range(steps + 1):
                time_s = START_S + index * STEP_S
                drift = 0.85 + 0.0022 * (time_s - START_S)
                noise = rng.gauss(0.0, 0.060)
                signal = (
                    drift
                    + noise
                    + gaussian(
                        time_s,
                        ANALYTE_CENTER_S,
                        ANALYTE_SIGMA_S,
                        analyte_amplitude,
                    )
                    + gaussian(
                        time_s,
                        INTERNAL_CENTER_S,
                        INTERNAL_SIGMA_S,
                        is_amplitude,
                    )
                )
                writer.writerow([run_id, f"{time_s:.1f}", f"{signal:.4f}"])


def main() -> None:
    task_root = Path(__file__).resolve().parents[1]
    data_dir = task_root / "environment" / "workspace" / "data"
    write_sequence(data_dir / "sequence.csv")
    write_chromatograms(data_dir / "chromatograms.csv")


if __name__ == "__main__":
    main()
