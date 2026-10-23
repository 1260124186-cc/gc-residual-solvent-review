"""Check GC integration, calibration, sample QC, and report validation."""

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
    """Return the anchor midpoint and mean signal."""
    selected = [
        signal
        for time, signal in zip(times, signals, strict=True)
        if bounds[0] <= time <= bounds[1]
    ]
    if len(selected) < 2:
        raise ValueError("anchor needs two points")
    return 0.5 * (bounds[0] + bounds[1]), sum(selected) / len(selected)


def _reference_integrate(times, signals, window):
    """Integrate a peak using an independent linear baseline."""
    start = _reference_anchor(times, signals, window.pre_anchor)
    end = _reference_anchor(times, signals, window.post_anchor)

    def baseline(time):
        """Interpolate the local baseline between the two anchors."""
        fraction = (time - start[0]) / (end[0] - start[0])
        return start[1] + fraction * (end[1] - start[1])

    points = [
        (time, signal - baseline(time))
        for time, signal in zip(times, signals, strict=True)
        if window.integration[0] <= time <= window.integration[1]
    ]
    area = sum(
        0.5 * (s0 + s1) * (t1 - t0)
        for (t0, s0), (t1, s1) in zip(points, points[1:], strict=False)
    )
    return {
        "area": area,
        "baseline_start_mv": baseline(points[0][0]),
        "baseline_end_mv": baseline(points[-1][0]),
        "points": len(points),
    }


def _reference_calibration(points):
    """Fit an independent weighted calibration line."""
    ordered = sorted(points, key=lambda item: float(item["nominal_mg_l"]))
    xs = [float(item["nominal_mg_l"]) for item in ordered]
    ys = [float(item["response_ratio"]) for item in ordered]
    weights = [1.0 / (x * x) for x in xs]
    wsum = sum(weights)
    xmean = sum(w * x for w, x in zip(weights, xs, strict=True)) / wsum
    ymean = sum(w * y for w, y in zip(weights, ys, strict=True)) / wsum
    sxx = sum(w * (x - xmean) ** 2 for w, x in zip(weights, xs, strict=True))
    sxy = sum(
        w * (x - xmean) * (y - ymean)
        for w, x, y in zip(weights, xs, ys, strict=True)
    )
    slope = sxy / sxx
    intercept = ymean - slope * xmean
    residuals = [
        y - intercept - slope * x for x, y in zip(xs, ys, strict=True)
    ]
    ss_res = sum(
        w * residual * residual
        for w, residual in zip(weights, residuals, strict=True)
    )
    ss_tot = sum(w * (y - ymean) ** 2 for w, y in zip(weights, ys, strict=True))
    r2 = 1.0 - ss_res / ss_tot
    sdyx = math.sqrt(ss_res / (len(points) - 2))
    loq = 10.0 * sdyx / slope
    recoveries = [
        100.0 * ((y - intercept) / slope) / x for x, y in zip(xs, ys, strict=True)
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
    """Load one run's time and signal arrays."""
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
    """Load the sequence rows."""
    with (DATA / "sequence.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _reference_report():
    """Build an independent reference report."""
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
    """Run the GC review CLI in a subprocess."""
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


def test_c1_linear_baseline_is_removed():
    """T01 题面「把前后两个基线点连成直线，作为整个积分窗口的基线…纵坐标是扣除该点线性基线后的信号」：线性基线必须被扣除，面积与两端基线值都要匹配独立实现。"""
    times = [float(index) for index in range(11)]
    signals = [2.0 + 0.5 * time for time in times]
    signals[5] += 3.0
    window = PeakWindow((3.0, 7.0), (0.0, 2.0), (8.0, 10.0))
    actual = integrate_peak(times, signals, window)
    expected = _reference_integrate(times, signals, window)
    assert float(actual["area"]) == pytest.approx(expected["area"], rel=1e-9, abs=1e-12)
    assert float(actual["baseline_start_mv"]) == pytest.approx(
        expected["baseline_start_mv"], rel=1e-9, abs=1e-12
    )
    assert float(actual["baseline_end_mv"]) == pytest.approx(
        expected["baseline_end_mv"], rel=1e-9, abs=1e-12
    )


def test_c1_trapezoid_peak_area_matches_reference():
    """T02 题面「面积用相邻采样点梯形法累加」：梯形积分面积和采样点数必须与独立实现一致。"""
    times = [36.0 + 0.1 * index for index in range(121)]
    signals = [
        1.1 + 0.013 * (time - 36.0)
        + 2.7 * math.exp(-0.5 * ((time - 42.0) / 0.48) ** 2)
        + 1.6 * math.exp(-0.5 * ((time - 57.0) / 0.55) ** 2)
        for time in times
    ]
    actual = integrate_peak(times, signals, ANALYTE_WINDOW)
    expected = _reference_integrate(times, signals, ANALYTE_WINDOW)
    assert float(actual["area"]) == pytest.approx(expected["area"], rel=1e-9, abs=1e-12)
    assert actual["points"] == expected["points"]


def test_c1_analyte_and_internal_standard_use_distinct_windows():
    """T03 题面「乙酸乙酯的积分窗口为39.2至44.8秒…内标峰使用54.4至59.6秒的积分窗口」：两个峰必须用各自窗口积分，面积不同。"""
    times, signals = _load_trace("CAL-3")
    analyte = _reference_integrate(times, signals, ANALYTE_WINDOW)
    internal = _reference_integrate(times, signals, INTERNAL_STANDARD_WINDOW)
    actual_analyte = integrate_peak(times, signals, ANALYTE_WINDOW)
    actual_internal = integrate_peak(times, signals, INTERNAL_STANDARD_WINDOW)
    assert float(actual_analyte["area"]) == pytest.approx(
        analyte["area"], rel=1e-10, abs=1e-12
    )
    assert float(actual_internal["area"]) == pytest.approx(
        internal["area"], rel=1e-10, abs=1e-12
    )
    assert float(actual_analyte["area"]) != pytest.approx(
        float(actual_internal["area"]), rel=1e-9, abs=1e-12
    )


def test_c1_rejects_nonmonotonic_trace():
    """T04 题面「同一进样的时间必须严格递增」：倒序时间轴必须被拒绝。"""
    with pytest.raises(ValueError, match="time|point|window|trace"):
        integrate_peak([0.0, 2.0, 1.0], [1.0, 1.0, 1.0], ANALYTE_WINDOW)


def test_c1_rejects_incomplete_window():
    """T05 题面「积分前锚点区间为36.4至38.0秒」：锚点区间采样不足时必须被拒绝。"""
    with pytest.raises(ValueError, match="time|point|window|trace"):
        integrate_peak([39.5, 40.0, 41.0], [1.0, 1.0, 1.0], ANALYTE_WINDOW)


def test_c1_rejects_single_point():
    """T06 题面「面积用相邻采样点梯形法累加」：单点轨迹无法积分，必须被拒绝。"""
    with pytest.raises(ValueError, match="time|point|window|trace"):
        integrate_peak([39.2], [1.0], ANALYTE_WINDOW)


def test_c2_weighted_least_squares_parameters():
    """T07 题面「每个校准点的权重为 1/x^2…求加权最小二乘解」：斜率、截距必须与独立加权拟合一致。"""
    points = [
        {"nominal_mg_l": 1.0, "response_ratio": 0.734},
        {"nominal_mg_l": 2.0, "response_ratio": 1.428},
        {"nominal_mg_l": 4.0, "response_ratio": 2.842},
        {"nominal_mg_l": 8.0, "response_ratio": 5.624},
        {"nominal_mg_l": 16.0, "response_ratio": 11.238},
    ]
    expected = _reference_calibration(points)
    actual = weighted_calibration(points)
    assert float(actual["slope"]) == pytest.approx(expected["slope"], rel=1e-9, abs=1e-12)
    assert float(actual["intercept"]) == pytest.approx(
        expected["intercept"], rel=1e-9, abs=1e-12
    )


def test_c2_weighted_r2_uses_weighted_mean():
    """T08 题面「用同一组权重计算加权决定系数 R2」：R2必须用加权均值计算。"""
    points = [
        {"nominal_mg_l": 1.0, "response_ratio": 0.741},
        {"nominal_mg_l": 2.0, "response_ratio": 1.415},
        {"nominal_mg_l": 4.0, "response_ratio": 2.836},
        {"nominal_mg_l": 8.0, "response_ratio": 5.647},
        {"nominal_mg_l": 16.0, "response_ratio": 11.221},
    ]
    expected = _reference_calibration(points)
    actual = weighted_calibration(points)
    assert float(actual["weighted_r2"]) == pytest.approx(
        expected["weighted_r2"], rel=1e-9, abs=1e-12
    )


def test_c2_standard_recovery_acceptance_rules():
    """T09 题面「最低浓度允许80%至120%，其余浓度允许85%至115%」：最低点与其余点用不同回收率区间判定。"""
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
    assert low_recovery == pytest.approx(
        expected_bad["recoveries"][0], rel=1e-9, abs=1e-12
    )


def test_c2_lod_and_loq_use_weighted_residual_error():
    """T10 题面「s_yx = sqrt(…)，LOD = 3.3*s_yx/b，LOQ = 10*s_yx/b」：sdyx、LOD、LOQ必须按加权残差标准误计算。"""
    points = [
        {"nominal_mg_l": 1.0, "response_ratio": 0.71},
        {"nominal_mg_l": 2.0, "response_ratio": 1.46},
        {"nominal_mg_l": 4.0, "response_ratio": 2.82},
        {"nominal_mg_l": 8.0, "response_ratio": 5.68},
        {"nominal_mg_l": 16.0, "response_ratio": 11.18},
    ]
    expected = _reference_calibration(points)
    actual = weighted_calibration(points)
    assert float(actual["sdyx"]) == pytest.approx(expected["sdyx"], rel=1e-9, abs=1e-12)
    assert float(actual["lod_mg_l"]) == pytest.approx(
        expected["lod_mg_l"], rel=1e-9, abs=1e-12
    )
    assert float(actual["loq_mg_l"]) == pytest.approx(
        expected["loq_mg_l"], rel=1e-9, abs=1e-12
    )


def test_c3_dilution_and_blank_correction_are_applied():
    """T11 题面「样品最终浓度等于 max(0, 仪器浓度*稀释倍数 - 空白校正值)」：稀释倍数和空白扣除必须同时生效。"""
    calibration = {"slope": 2.0, "intercept": 0.1}
    actual = corrected_concentration(1.1, 10.0, calibration, 0.2)
    assert actual == pytest.approx(4.8, rel=1e-9, abs=1e-12)
    assert corrected_concentration(0.1, 10.0, calibration, 0.2) == 0.0


def test_c3_duplicate_rpd_boundary():
    """T12 题面「RPD=100*abs(c1-c2)/((c1+c2)/2)；两个浓度均为零时RPD记为0」：RPD公式和零值分支必须正确。"""
    values = [100.0, 110.5263157894737]
    assert relative_percent_difference(values) == pytest.approx(
        10.0, rel=1e-12, abs=1e-12
    )
    assert relative_percent_difference([0.0, 0.0]) == 0.0
    assert relative_percent_difference([100.0, 112.0]) > 10.0


def test_c3_spike_recovery_formula_and_boundaries():
    """T13 题面「加标回收率为 100*(加标样最终浓度-平行样平均浓度)/spike_added_mg_l」：回收率公式和80/120边界必须正确。"""
    assert spike_recovery_percent(18.0, 10.0, 10.0) == pytest.approx(
        80.0, rel=1e-9, abs=1e-12
    )
    assert spike_recovery_percent(22.0, 10.0, 10.0) == pytest.approx(
        120.0, rel=1e-9, abs=1e-12
    )
    assert spike_recovery_percent(17.9, 10.0, 10.0) < 80.0


def test_c3_sample_status_propagation():
    """T14 题面「低于LOQ为below_loq，平行样或加标不合格为qc_failed；任一qc_failed整批fail」：状态传播必须正确。"""
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
    """T15 题面「输出的根对象固定包含 calibration、blank_mean_mg_l、blank_correction_mg_l、qc、samples、overall_status」：报告字段和样例结论必须齐全。"""
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
    """T16 题面「把原始色谱和序列表整理成可以直接交给QA复核的JSON报告」：整条序列的校准、样品和状态必须与独立实现一致。"""
    expected = _reference_report()
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "report.json"
        process = _run_cli(DATA / "sequence.csv", DATA / "chromatograms.csv", output)
        actual = json.loads(output.read_text(encoding="utf-8"))
    assert process.returncode == 0
    calibration_keys = (
        "slope",
        "intercept",
        "weighted_r2",
        "sdyx",
        "lod_mg_l",
        "loq_mg_l",
    )
    actual_calibration = {
        key: float(actual["calibration"][key]) for key in calibration_keys
    }
    expected_calibration = {
        key: expected["calibration"][key] for key in calibration_keys
    }
    assert actual_calibration == pytest.approx(expected_calibration, rel=1e-8, abs=1e-12)
    assert actual["calibration"]["accepted"] == expected["calibration"]["accepted"]
    standards = actual["calibration"]["standards"]
    assert len(standards) == len(expected["calibration"]["recoveries"]) == 5
    assert [float(item["nominal_mg_l"]) for item in standards] == pytest.approx(
        [1.0, 2.0, 4.0, 8.0, 16.0]
    )
    assert [float(item["recovery_percent"]) for item in standards] == pytest.approx(
        expected["calibration"]["recoveries"], rel=1e-8, abs=1e-12
    )
    assert set(actual["samples"]) == set(expected["samples"]) == {"S1", "S2"}
    numeric_sample_keys = (
        "concentration_mg_l", "duplicate_rpd_percent", "spike_recovery_percent"
    )
    actual_samples = {
        (sample_id, key): float(actual["samples"][sample_id][key])
        for sample_id in ("S1", "S2")
        for key in numeric_sample_keys
    }
    expected_samples = {
        (sample_id, key): expected["samples"][sample_id][key]
        for sample_id in ("S1", "S2")
        for key in numeric_sample_keys
    }
    assert actual_samples == pytest.approx(expected_samples, rel=1e-8, abs=1e-12)
    assert {
        sample_id: sample["status"] for sample_id, sample in actual["samples"].items()
    } == {
        sample_id: sample["status"] for sample_id, sample in expected["samples"].items()
    }


def _invalid_input_cases(root):
    """Build one file pair for each invalid-input scenario."""
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

    return {
        "missing-column": (missing_column, DATA / "chromatograms.csv"),
        "non-numeric": (non_numeric, DATA / "chromatograms.csv"),
        "duplicate": (duplicate, DATA / "chromatograms.csv"),
        "trace-mismatch": (mismatch, mismatch_chromatograms),
    }


def test_e2e_rejects_missing_column_without_report():
    """T17 题面「输入缺列…程序应以非零状态退出，不能生成报告」：缺列必须拒绝且不写半成品。"""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        sequence_path, chromatogram_path = _invalid_input_cases(root)["missing-column"]
        output = root / "report.json"
        process = _run_cli(sequence_path, chromatogram_path, output)
        assert process.returncode != 0
        assert not output.exists()


def test_e2e_rejects_non_numeric_field_without_report():
    """T18 题面「数值非有限…程序应以非零状态退出，不能生成报告」：非数值必须拒绝且不写半成品。"""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        sequence_path, chromatogram_path = _invalid_input_cases(root)["non-numeric"]
        output = root / "report.json"
        process = _run_cli(sequence_path, chromatogram_path, output)
        assert process.returncode != 0
        assert not output.exists()


def test_e2e_rejects_duplicate_run_without_report():
    """T19 题面「进样号在序列中必须唯一」：重复进样号必须拒绝且不写半成品。"""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        sequence_path, chromatogram_path = _invalid_input_cases(root)["duplicate"]
        output = root / "report.json"
        process = _run_cli(sequence_path, chromatogram_path, output)
        assert process.returncode != 0
        assert not output.exists()


def test_e2e_rejects_trace_mismatch_without_report():
    """T20 题面「色谱与序列的进样号集合必须完全一致」：进样号不一致必须拒绝且不写半成品。"""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        sequence_path, chromatogram_path = _invalid_input_cases(root)["trace-mismatch"]
        output = root / "report.json"
        process = _run_cli(sequence_path, chromatogram_path, output)
        assert process.returncode != 0
        assert not output.exists()
