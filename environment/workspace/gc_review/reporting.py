"""Sample concentration and sequence review helpers."""

from __future__ import annotations

import math


def instrument_concentration(response_ratio: float, calibration: dict[str, object]) -> float:
    slope = float(calibration["slope"])
    intercept = float(calibration["intercept"])
    return (response_ratio - intercept) / slope


def corrected_concentration(
    response_ratio: float,
    dilution_factor: float,
    calibration: dict[str, object],
    blank_correction_mg_l: float,
) -> float:
    raw = instrument_concentration(response_ratio, calibration)
    return max(0.0, raw - blank_correction_mg_l)


def relative_percent_difference(values: list[float]) -> float:
    if len(values) != 2:
        raise ValueError("exactly two replicate values are required")
    first, second = values
    mean = 0.5 * (first + second)
    if mean == 0:
        return 0.0
    return 100.0 * abs(first - second) / mean


def spike_recovery_percent(
    spike_concentration_mg_l: float,
    native_concentration_mg_l: float,
    spike_added_mg_l: float,
) -> float:
    return 100.0 * (spike_concentration_mg_l - native_concentration_mg_l) / spike_added_mg_l


def review_samples(
    sequence: list[dict[str, object]],
    ratios: dict[str, float],
    calibration: dict[str, object],
) -> dict[str, object]:
    """Build a report with the legacy single-dilution behavior."""
    blanks = [
        instrument_concentration(ratios[row["run_id"]], calibration)
        for row in sequence
        if row["role"] == "blank"
    ]
    blank_mean = sum(blanks) / len(blanks) if blanks else 0.0
    samples: dict[str, dict[str, object]] = {}
    for sample_id in sorted({str(row["sample_id"]) for row in sequence if row["role"] == "sample"}):
        values = [
            instrument_concentration(ratios[row["run_id"]], calibration)
            for row in sequence
            if row["role"] == "sample" and str(row["sample_id"]) == sample_id
        ]
        samples[sample_id] = {
            "concentration_mg_l": sum(values) / len(values),
            "replicates_mg_l": values,
            "duplicate_rpd_percent": relative_percent_difference(values),
            "spike_recovery_percent": math.nan,
            "status": "pass",
        }
    return samples
