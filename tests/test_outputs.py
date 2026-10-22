import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, "/workspace")

from gc_review.calibration import weighted_calibration
from gc_review.chromatogram import (
    ANALYTE_WINDOW,
    INTERNAL_STANDARD_WINDOW,
    PeakWindow,
    integrate_peak,
)
from gc_review.reporting import (
    corrected_concentration,
    relative_percent_difference,
    review_samples,
    spike_recovery_percent,
)


DATA = Path("/workspace/data")


def _reference_anchor(times, signals, bounds):
    selected = [
        signal
        for time, signal in zip(times, signals)
        if bounds[0] <= time <= bounds[1]
    ]
    if len(selected) < 2:
        raise ValueError("anchor needs two points")
    return 0.5 * (bounds[0] + bounds[1]), sum(selected) / len(selected)


def _reference_integrate(times, signals, window):
    start = _reference_anchor(times, signals, window.pre_anchor)
    end = _reference_anchor(times, signals, window.post_anchor)

    def baseline(time):
        fraction = (time - start[0]) / (end[0] - start[0])
        return start[1] + fraction * (end[1] - start[1])

    points = [
        (time, signal - baseline(time))
        for time, signal in zip(times, signals)
        if window.integration[0] <= time <= window.integration[1]
    ]
    area = sum(
        0.5 * (s0 + s1) * (t1 - t0)
        for (t0, s0), (t1, s1) in zip(points, points[1:])
    )
    return {
        "area": area,
        "baseline_start_mv": baseline(points[0][0]),
        "baseline_end_mv": baseline(points[-1][0]),
        "points": len(points),
    }


def _reference_calibration(points):
    ordered = sorted(points, key=lambda item: float(item["nominal_mg_l"]))
    xs = [float(item["nominal_mg_l"]) for item in ordered]
    ys = [float(item["response_ratio"]) for item in ordered]
    weights = [1.0 / (x * x) for x in xs]
    wsum = sum(weights)
    xmean = sum(w * x for w, x in zip(weights, xs)) / wsum
    ymean = sum(w * y for w, y in zip(weights, ys)) / wsum
    sxx = sum(w * (x - xmean) ** 2 for w, x in zip(weights, xs))
    sxy = sum(
        w * (x - xmean) * (y - ymean)
        for w, x, y in zip(weights, xs, ys)
    )
    slope = sxy / sxx
    intercept = ymean - slope * xmean
    residuals = [
        y - intercept - slope * x for x, y in zip(xs, ys)
    ]
    ss_res = sum(
        w * residual * residual
        for w, residual in zip(weights, residuals)
    )
    ss_tot = sum(w * (y - ymean) ** 2 for w, y in zip(weights, ys))
    r2 = 1.0 - ss_res / ss_tot
    sdyx = math.sqrt(ss_res / (len(points) - 2))
    loq = 10.0 * sdyx / slope
    recoveries = [
        100.0 * ((y - intercept) / slope) / x for x, y in zip(xs, ys)
    ]
    accepted = (
        len(ordered) >= 5
        and slope > 0
        and r2 >= 0.995
        and 80.0 <= recoveries[0] <= 120.0
        and all(85.0 <= value <= 115.0 for value in recoveries[1:])
    )
    return {
        "slope": slope,
        "intercept": intercept,
        "weighted_r2": r2,
        "sdyx": sdyx,
        "lod_mg_l": 3.3 * sdyx / slope,
        "loq_mg_l": loq,
        "accepted": accepted,
        "recoveries": recoveries,
    }


def _load_trace(run_id):
    times = []
    signals = []
    with (DATA / "chromatograms.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["run_id"] == run_id:
                times.append(float(row["time_s"]))
                signals.append(float(row["signal_mv"]))
    if not times:
        raise ValueError(f"missing trace {run_id}")
    return times, signals


def _load_sequence():
    with (DATA / "sequence.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _reference_report():
    rows = _load_sequence()
    ratios = {}
    for row in rows:
        times, signals = _load_trace(row["run_id"])
        analyte = _reference_integrate(times, signals, ANALYTE_WINDOW)
        internal = _reference_integrate(times, signals, INTERNAL_STANDARD_WINDOW)
        ratios[row["run_id"]] = analyte["area"] / internal["area"]
    calibration_points = [
        {
            "nominal_mg_l": float(row["nominal_mg_l"]),
            "response_ratio": ratios[row["run_id"]],
        }
        for row in rows
        if row["role"] == "calibration"
    ]
    calibration = _reference_calibration(calibration_points)
    blank_values = [
        max(
            0.0,
            (ratios[row["run_id"]] - calibration["intercept"])
            / calibration["slope"]
            * float(row["dilution_factor"]),
        )
        for row in rows
        if row["role"] == "blank"
    ]
    blank_mean = sum(blank_values) / len(blank_values)
    correction = max(0.0, blank_mean)

    samples = {}
    duplicate_ok = True
    spike_ok = True
    for sample_id in sorted(
        {row["sample_id"] for row in rows if row["role"] == "sample"}
    ):
        sample_rows = sorted(
            [
                row
                for row in rows
                if row["role"] == "sample" and row["sample_id"] == sample_id
            ],
            key=lambda item: int(item["replicate"]),
        )
        spike_row = next(
            row
            for row in rows
            if row["role"] == "spike" and row["sample_id"] == sample_id
        )
        values = [
            max(
                0.0,
                (ratios[row["run_id"]] - calibration["intercept"])
                / calibration["slope"]
                * float(row["dilution_factor"])
                - correction,
            )
            for row in sample_rows
        ]
        native = sum(values) / 2.0
        spike_value = max(
            0.0,
            (ratios[spike_row["run_id"]] - calibration["intercept"])
            / calibration["slope"]
            * float(spike_row["dilution_factor"])
            - correction,
        )
        rpd = (
            0.0
            if abs(values[0] + values[1]) < 1e-30
            else 100.0 * abs(values[0] - values[1]) / (0.5 * sum(values))
        )
        recovery = (
            100.0
            * (spike_value - native)
            / float(spike_row["spike_added_mg_l"])
        )
        duplicate_ok = duplicate_ok and rpd <= 10.0
        spike_ok = spike_ok and 80.0 <= recovery <= 120.0
        if not calibration["accepted"] or correction > 0.10:
            status = "qc_failed"
        elif native < calibration["loq_mg_l"]:
            status = "below_loq"
        elif rpd > 10.0 or not 80.0 <= recovery <= 120.0:
            status = "qc_failed"
        else:
            status = "pass"
        samples[sample_id] = {
            "concentration_mg_l": native,
            "replicates_mg_l": values,
            "duplicate_rpd_percent": rpd,
            "spike_recovery_percent": recovery,
            "status": status,
        }
    return {
        "calibration": calibration,
        "blank_mean_mg_l": blank_mean,
        "blank_correction_mg_l": correction,
        "qc": {
            "calibration": "pass" if calibration["accepted"] else "fail",
            "blank": "pass" if correction <= 0.10 else "fail",
            "duplicate": "pass" if duplicate_ok else "fail",
            "spike_recovery": "pass" if spike_ok else "fail",
        },
        "samples": samples,
        "overall_status": (
            "fail"
            if any(item["status"] == "qc_failed" for item in samples.values())
            else "pass"
        ),
    }


def _run_cli(sequence_path, chromatogram_path, output_path):
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "gc_review.cli",
            str(sequence_path),
            str(chromatogram_path),
            str(output_path),
        ],
        cwd="/workspace",
        capture_output=True,
        text=True,
        timeout=60,
    )


def _close(actual, expected, rel_tol=1e-9, abs_tol=1e-12):
    return math.isfinite(actual) and math.isclose(
        actual, expected, rel_tol=rel_tol, abs_tol=abs_tol
    )


def test_c1_linear_baseline_is_removed():
    times = [float(index) for index in range(11)]
    signals = [2.0 + 0.5 * time for time in times]
    signals[5] += 3.0
    window = PeakWindow((3.0, 7.0), (0.0, 2.0), (8.0, 10.0))
    actual = integrate_peak(times, signals, window)
    expected = _reference_integrate(times, signals, window)
    assert _close(float(actual["area"]), expected["area"])
    assert _close(float(actual["baseline_start_mv"]), expected["baseline_start_mv"])
    assert _close(float(actual["baseline_end_mv"]), expected["baseline_end_mv"])


def test_c1_trapezoid_peak_area_matches_reference():
    times = [36.0 + 0.1 * index for index in range(121)]
    signals = [
        1.1 + 0.013 * (time - 36.0)
        + 2.7 * math.exp(-0.5 * ((time - 42.0) / 0.48) ** 2)
        + 1.6 * math.exp(-0.5 * ((time - 57.0) / 0.55) ** 2)
        for time in times
    ]
    actual = integrate_peak(times, signals, ANALYTE_WINDOW)
    expected = _reference_integrate(times, signals, ANALYTE_WINDOW)
    assert _close(float(actual["area"]), expected["area"])
    assert actual["points"] == expected["points"]


def test_c1_analyte_and_internal_standard_use_distinct_windows():
    times, signals = _load_trace("CAL-3")
    analyte = _reference_integrate(times, signals, ANALYTE_WINDOW)
    internal = _reference_integrate(times, signals, INTERNAL_STANDARD_WINDOW)
    actual_analyte = integrate_peak(times, signals, ANALYTE_WINDOW)
    actual_internal = integrate_peak(times, signals, INTERNAL_STANDARD_WINDOW)
    assert _close(float(actual_analyte["area"]), analyte["area"], rel_tol=1e-10)
    assert _close(float(actual_internal["area"]), internal["area"], rel_tol=1e-10)
    assert not _close(float(actual_analyte["area"]), float(actual_internal["area"]))


def test_c1_rejects_nonmonotonic_or_incomplete_trace():
    with pytest.raises(ValueError):
        integrate_peak([0.0, 2.0, 1.0], [1.0, 1.0, 1.0], ANALYTE_WINDOW)
    with pytest.raises(ValueError):
        integrate_peak([39.5, 40.0, 41.0], [1.0, 1.0, 1.0], ANALYTE_WINDOW)
    with pytest.raises(ValueError):
        integrate_peak([39.2], [1.0], ANALYTE_WINDOW)


def test_c2_weighted_least_squares_parameters():
    points = [
        {"nominal_mg_l": 1.0, "response_ratio": 0.734},
        {"nominal_mg_l": 2.0, "response_ratio": 1.428},
        {"nominal_mg_l": 4.0, "response_ratio": 2.842},
        {"nominal_mg_l": 8.0, "response_ratio": 5.624},
        {"nominal_mg_l": 16.0, "response_ratio": 11.238},
    ]
    expected = _reference_calibration(points)
    actual = weighted_calibration(points)
    assert _close(float(actual["slope"]), expected["slope"])
    assert _close(float(actual["intercept"]), expected["intercept"])


def test_c2_weighted_r2_uses_weighted_mean():
    points = [
        {"nominal_mg_l": 1.0, "response_ratio": 0.741},
        {"nominal_mg_l": 2.0, "response_ratio": 1.415},
        {"nominal_mg_l": 4.0, "response_ratio": 2.836},
        {"nominal_mg_l": 8.0, "response_ratio": 5.647},
        {"nominal_mg_l": 16.0, "response_ratio": 11.221},
    ]
    expected = _reference_calibration(points)
    actual = weighted_calibration(points)
    assert _close(float(actual["weighted_r2"]), expected["weighted_r2"])


def test_c2_standard_recovery_acceptance_rules():
    points = [
        {"nominal_mg_l": 1.0, "response_ratio": 0.73},
        {"nominal_mg_l": 2.0, "response_ratio": 1.43},
        {"nominal_mg_l": 4.0, "response_ratio": 2.83},
        {"nominal_mg_l": 8.0, "response_ratio": 5.63},
        {"nominal_mg_l": 16.0, "response_ratio": 11.23},
    ]
    good = weighted_calibration(points)
    assert good["accepted"] is True
    bad_points = [dict(point) for point in points]
    bad_points[0]["response_ratio"] *= 1.30
    bad = weighted_calibration(bad_points)
    expected_bad = _reference_calibration(bad_points)
    assert bad["accepted"] is False
    low_recovery = next(
        item["recovery_percent"]
        for item in bad["standards"]
        if item["nominal_mg_l"] == 1.0
    )
    assert _close(low_recovery, expected_bad["recoveries"][0])


def test_c2_lod_and_loq_use_weighted_residual_error():
    points = [
        {"nominal_mg_l": 1.0, "response_ratio": 0.71},
        {"nominal_mg_l": 2.0, "response_ratio": 1.46},
        {"nominal_mg_l": 4.0, "response_ratio": 2.82},
        {"nominal_mg_l": 8.0, "response_ratio": 5.68},
        {"nominal_mg_l": 16.0, "response_ratio": 11.18},
    ]
    expected = _reference_calibration(points)
    actual = weighted_calibration(points)
    assert _close(float(actual["sdyx"]), expected["sdyx"])
    assert _close(float(actual["lod_mg_l"]), expected["lod_mg_l"])
    assert _close(float(actual["loq_mg_l"]), expected["loq_mg_l"])


def test_c3_dilution_and_blank_correction_are_applied():
    calibration = {"slope": 2.0, "intercept": 0.1}
    actual = corrected_concentration(1.1, 10.0, calibration, 0.2)
    assert _close(actual, 4.8)
    assert corrected_concentration(0.1, 10.0, calibration, 0.2) == 0.0


def test_c3_duplicate_rpd_boundary():
    values = [100.0, 110.5263157894737]
    assert _close(relative_percent_difference(values), 10.0, rel_tol=1e-12)
    assert relative_percent_difference([0.0, 0.0]) == 0.0
    assert relative_percent_difference([100.0, 112.0]) > 10.0


def test_c3_spike_recovery_formula_and_boundaries():
    assert _close(spike_recovery_percent(18.0, 10.0, 10.0), 80.0)
    assert _close(spike_recovery_percent(22.0, 10.0, 10.0), 120.0)
    assert spike_recovery_percent(17.9, 10.0, 10.0) < 80.0


def test_c3_sample_status_propagation():
    calibration = {
        "slope": 1.0,
        "intercept": 0.0,
        "loq_mg_l": 0.5,
        "accepted": True,
    }
    sequence = [
        {"run_id": "BLK", "sample_id": "BLANK", "role": "blank", "replicate": 1, "dilution_factor": 1.0},
        {"run_id": "A1", "sample_id": "A", "role": "sample", "replicate": 1, "dilution_factor": 1.0},
        {"run_id": "A2", "sample_id": "A", "role": "sample", "replicate": 2, "dilution_factor": 1.0},
        {"run_id": "AS", "sample_id": "A", "role": "spike", "replicate": 1, "dilution_factor": 1.0, "spike_added_mg_l": 8.0},
        {"run_id": "B1", "sample_id": "B", "role": "sample", "replicate": 1, "dilution_factor": 1.0},
        {"run_id": "B2", "sample_id": "B", "role": "sample", "replicate": 2, "dilution_factor": 1.0},
        {"run_id": "BS", "sample_id": "B", "role": "spike", "replicate": 1, "dilution_factor": 1.0, "spike_added_mg_l": 5.0},
        {"run_id": "C1", "sample_id": "C", "role": "sample", "replicate": 1, "dilution_factor": 1.0},
        {"run_id": "C2", "sample_id": "C", "role": "sample", "replicate": 2, "dilution_factor": 1.0},
        {"run_id": "CS", "sample_id": "C", "role": "spike", "replicate": 1, "dilution_factor": 1.0, "spike_added_mg_l": 8.0},
    ]
    ratios = {
        "BLK": 0.0,
        "A1": 10.0,
        "A2": 10.2,
        "AS": 18.1,
        "B1": 0.2,
        "B2": 0.21,
        "BS": 5.1,
        "C1": 10.0,
        "C2": 13.0,
        "CS": 18.0,
    }
    review = review_samples(sequence, ratios, calibration)
    assert review["samples"]["A"]["status"] == "pass"
    assert review["samples"]["B"]["status"] == "below_loq"
    assert review["samples"]["C"]["status"] == "qc_failed"
    assert review["qc"]["duplicate"] == "fail"
    assert review["overall_status"] == "fail"


def test_e2e_report_schema_and_summary_values():
    reference = _reference_report()
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "report.json"
        process = _run_cli(DATA / "sequence.csv", DATA / "chromatograms.csv", output)
        report = json.loads(output.read_text(encoding="utf-8"))
    assert process.returncode == 0
    assert set(report) == {
        "calibration",
        "blank_mean_mg_l",
        "blank_correction_mg_l",
        "qc",
        "samples",
        "overall_status",
    }
    assert set(report["calibration"]) == {
        "slope",
        "intercept",
        "weighted_r2",
        "sdyx",
        "lod_mg_l",
        "loq_mg_l",
        "accepted",
        "standards",
    }
    assert report["overall_status"] == reference["overall_status"]
    assert report["samples"]["S1"]["status"] == "pass"
    assert report["samples"]["S2"]["status"] == "below_loq"


def test_e2e_full_sequence_matches_independent_reference():
    expected = _reference_report()
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "report.json"
        process = _run_cli(DATA / "sequence.csv", DATA / "chromatograms.csv", output)
        actual = json.loads(output.read_text(encoding="utf-8"))
    assert process.returncode == 0
    for key in (
        "slope",
        "intercept",
        "weighted_r2",
        "sdyx",
        "lod_mg_l",
        "loq_mg_l",
    ):
        assert _close(float(actual["calibration"][key]), expected["calibration"][key], rel_tol=1e-8)
    assert actual["calibration"]["accepted"] == expected["calibration"]["accepted"]
    for index, standard in enumerate(actual["calibration"]["standards"]):
        nominal = 1.0 * (2 ** index)
        assert _close(float(standard["nominal_mg_l"]), nominal)
        assert _close(
            float(standard["recovery_percent"]),
            expected["calibration"]["recoveries"][index],
            rel_tol=1e-8,
        )
    for sample_id in ("S1", "S2"):
        actual_sample = actual["samples"][sample_id]
        expected_sample = expected["samples"][sample_id]
        for key in ("concentration_mg_l", "duplicate_rpd_percent", "spike_recovery_percent"):
            assert _close(float(actual_sample[key]), expected_sample[key], rel_tol=1e-8)
        assert actual_sample["status"] == expected_sample["status"]


def test_e2e_rejects_invalid_inputs_without_report():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source_sequence = (DATA / "sequence.csv").read_text(encoding="utf-8")
        source_chromatograms = (DATA / "chromatograms.csv").read_text(encoding="utf-8")

        missing_column = root / "missing.csv"
        lines = source_sequence.splitlines()
        missing_column.write_text(
            "\n".join(line.rsplit(",", 1)[0] for line in lines) + "\n",
            encoding="utf-8",
        )
        non_numeric = root / "bad_number.csv"
        non_numeric.write_text(
            source_sequence.replace(",1.0,1.0,", ",1.0,not-a-number,", 1),
            encoding="utf-8",
        )
        duplicate = root / "duplicate.csv"
        duplicate_lines = source_sequence.splitlines()
        duplicate_lines[2] = duplicate_lines[1]
        duplicate.write_text("\n".join(duplicate_lines) + "\n", encoding="utf-8")
        mismatch = root / "mismatch.csv"
        mismatch_chromatograms = root / "mismatch_chrom.csv"
        mismatch_chromatograms.write_text(
            source_chromatograms.replace("S1-R1", "UNKNOWN", 196),
            encoding="utf-8",
        )
        mismatch.write_text(source_sequence, encoding="utf-8")

        cases = [
            (missing_column, DATA / "chromatograms.csv"),
            (non_numeric, DATA / "chromatograms.csv"),
            (duplicate, DATA / "chromatograms.csv"),
            (mismatch, mismatch_chromatograms),
        ]
        for index, (sequence_path, chromatogram_path) in enumerate(cases):
            output = root / f"out-{index}.json"
            process = _run_cli(sequence_path, chromatogram_path, output)
            assert process.returncode != 0
            assert not output.exists()
