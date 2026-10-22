"""Sample concentration and sequence review helpers."""

from __future__ import annotations

import math


def instrument_concentration(response_ratio: float, calibration: dict[str, object]) -> float:
    slope = float(calibration["slope"])
    intercept = float(calibration["intercept"])
    if slope == 0 or not math.isfinite(slope):
        raise ValueError("calibration slope must be finite and non-zero")
    return (response_ratio - intercept) / slope


def corrected_concentration(
    response_ratio: float,
    dilution_factor: float,
    calibration: dict[str, object],
    blank_correction_mg_l: float,
) -> float:
    raw = instrument_concentration(response_ratio, calibration)
    return max(0.0, raw * dilution_factor - blank_correction_mg_l)


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
    if spike_added_mg_l <= 0:
        raise ValueError("spike addition must be positive")
    return 100.0 * (spike_concentration_mg_l - native_concentration_mg_l) / spike_added_mg_l


def sample_summary(
    sample_id: str,
    sample_rows: list[dict[str, object]],
    spike_row: dict[str, object],
    ratios: dict[str, float],
    calibration: dict[str, object],
    blank_correction_mg_l: float,
) -> dict[str, object]:
    replicates = [
        corrected_concentration(
            ratios[str(row["run_id"])],
            float(row["dilution_factor"]),
            calibration,
            blank_correction_mg_l,
        )
        for row in sorted(sample_rows, key=lambda item: int(item["replicate"]))
    ]
    native = sum(replicates) / len(replicates)
    spike = corrected_concentration(
        ratios[str(spike_row["run_id"])],
        float(spike_row["dilution_factor"]),
        calibration,
        blank_correction_mg_l,
    )
    rpd = relative_percent_difference(replicates)
    recovery = spike_recovery_percent(
        spike,
        native,
        float(spike_row["spike_added_mg_l"]),
    )
    return {
        "concentration_mg_l": native,
        "replicates_mg_l": replicates,
        "duplicate_rpd_percent": rpd,
        "spike_recovery_percent": recovery,
    }


def review_samples(
    sequence: list[dict[str, object]],
    ratios: dict[str, float],
    calibration: dict[str, object],
) -> dict[str, object]:
    blank_values = [
        max(
            0.0,
            instrument_concentration(ratios[str(row["run_id"])], calibration)
            * float(row["dilution_factor"]),
        )
        for row in sequence
        if row["role"] == "blank"
    ]
    blank_mean = sum(blank_values) / len(blank_values) if blank_values else 0.0
    blank_correction = max(0.0, blank_mean)
    calibration_ok = bool(calibration["accepted"])
    blank_ok = blank_correction <= 0.10
    loq = float(calibration["loq_mg_l"])

    sample_ids = sorted(
        {
            str(row["sample_id"])
            for row in sequence
            if row["role"] == "sample"
        }
    )
    samples: dict[str, dict[str, object]] = {}
    all_duplicate_ok = True
    all_spike_ok = True
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
        summary = sample_summary(
            sample_id,
            sample_rows,
            spike_rows[0],
            ratios,
            calibration,
            blank_correction,
        )
        rpd = float(summary["duplicate_rpd_percent"])
        recovery = float(summary["spike_recovery_percent"])
        duplicate_ok = math.isfinite(rpd) and rpd <= 10.0
        spike_ok = math.isfinite(recovery) and 80.0 <= recovery <= 120.0
        all_duplicate_ok = all_duplicate_ok and duplicate_ok
        all_spike_ok = all_spike_ok and spike_ok
        if not calibration_ok or not blank_ok:
            status = "qc_failed"
        elif float(summary["concentration_mg_l"]) < loq:
            status = "below_loq"
        elif not duplicate_ok or not spike_ok:
            status = "qc_failed"
        else:
            status = "pass"
        summary["status"] = status
        samples[sample_id] = summary

    overall_fail = any(item["status"] == "qc_failed" for item in samples.values())
    return {
        "blank_mean_mg_l": blank_mean,
        "blank_correction_mg_l": blank_correction,
        "qc": {
            "calibration": "pass" if calibration_ok else "fail",
            "blank": "pass" if blank_ok else "fail",
            "duplicate": "pass" if all_duplicate_ok else "fail",
            "spike_recovery": "pass" if all_spike_ok else "fail",
        },
        "samples": samples,
        "overall_status": "fail" if overall_fail else "pass",
    }
