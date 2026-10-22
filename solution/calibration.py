"""Internal-standard calibration helpers."""

from __future__ import annotations

import math


def weighted_calibration(points: list[dict[str, float]]) -> dict[str, object]:
    """Fit slope and intercept by 1/x^2 weighted least squares."""
    if len(points) < 2:
        raise ValueError("at least two calibration points are required")
    ordered = sorted(points, key=lambda point: float(point["nominal_mg_l"]))
    xs = [float(point["nominal_mg_l"]) for point in ordered]
    ys = [float(point["response_ratio"]) for point in ordered]
    if any(not math.isfinite(value) or value <= 0 for value in xs):
        raise ValueError("calibration concentrations must be positive and finite")
    if any(not math.isfinite(value) for value in ys):
        raise ValueError("calibration ratios must be finite")
    if len(set(xs)) != len(xs):
        raise ValueError("calibration concentrations must be distinct")

    weights = [1.0 / (x * x) for x in xs]
    weight_sum = sum(weights)
    x_mean = sum(weight * x for weight, x in zip(weights, xs)) / weight_sum
    y_mean = sum(weight * y for weight, y in zip(weights, ys)) / weight_sum
    sxx = sum(weight * (x - x_mean) ** 2 for weight, x in zip(weights, xs))
    sxy = sum(
        weight * (x - x_mean) * (y - y_mean)
        for weight, x, y in zip(weights, xs, ys)
    )
    if sxx == 0:
        raise ValueError("calibration concentrations must span more than one value")
    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    predicted = [intercept + slope * x for x in xs]
    residuals = [y - fit for y, fit in zip(ys, predicted)]
    ss_res = sum(
        weight * residual * residual
        for weight, residual in zip(weights, residuals)
    )
    ss_tot = sum(
        weight * (y - y_mean) ** 2 for weight, y in zip(weights, ys)
    )
    r2 = 1.0 - ss_res / ss_tot if ss_tot else 1.0
    sdyx = math.sqrt(ss_res / (len(points) - 2)) if len(points) > 2 else 0.0
    lod = 3.3 * sdyx / slope if slope > 0 else math.inf
    loq = 10.0 * sdyx / slope if slope > 0 else math.inf

    standards = []
    recoveries_ok = True
    for index, (x, y) in enumerate(zip(xs, ys)):
        recovery = 100.0 * ((y - intercept) / slope) / x if slope else math.nan
        lower, upper = (80.0, 120.0) if index == 0 else (85.0, 115.0)
        if not math.isfinite(recovery) or not lower <= recovery <= upper:
            recoveries_ok = False
        standards.append(
            {
                "nominal_mg_l": x,
                "response_ratio": y,
                "recovery_percent": recovery,
            }
        )

    accepted = (
        len(points) >= 5
        and slope > 0
        and math.isfinite(r2)
        and r2 >= 0.995
        and recoveries_ok
    )
    return {
        "slope": slope,
        "intercept": intercept,
        "weighted_r2": r2,
        "sdyx": sdyx,
        "lod_mg_l": lod,
        "loq_mg_l": loq,
        "accepted": accepted,
        "standards": standards,
    }
