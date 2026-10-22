"""Internal-standard calibration helpers."""

from __future__ import annotations

import math


def weighted_calibration(points: list[dict[str, float]]) -> dict[str, object]:
    """Fit the calibration points and return method figures."""
    if len(points) < 2:
        raise ValueError("at least two calibration points are required")
    xs = [float(point["nominal_mg_l"]) for point in points]
    ys = [float(point["response_ratio"]) for point in points]
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    sxx = sum((x - x_mean) ** 2 for x in xs)
    sxy = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    predicted = [intercept + slope * x for x in xs]
    residual = [y - fit for y, fit in zip(ys, predicted)]
    total = [y - y_mean for y in ys]
    ss_res = sum(value * value for value in residual)
    ss_tot = sum(value * value for value in total)
    r2 = 1.0 - ss_res / ss_tot if ss_tot else 1.0
    sdyx = math.sqrt(ss_res / max(1, len(points) - 2))
    lod = 3.3 * sdyx / slope if slope > 0 else math.inf
    loq = 10.0 * sdyx / slope if slope > 0 else math.inf
    standards = [
        {
            "nominal_mg_l": x,
            "response_ratio": y,
            "recovery_percent": 100.0 * ((y - intercept) / slope) / x if slope else math.nan,
        }
        for x, y in zip(xs, ys)
    ]
    return {
        "slope": slope,
        "intercept": intercept,
        "weighted_r2": r2,
        "sdyx": sdyx,
        "lod_mg_l": lod,
        "loq_mg_l": loq,
        "accepted": len(points) >= 5 and slope > 0 and r2 >= 0.995,
        "standards": standards,
    }
