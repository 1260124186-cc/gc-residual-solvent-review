"""Peak integration helpers for GC-FID traces."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PeakWindow:
    integration: tuple[float, float]
    pre_anchor: tuple[float, float]
    post_anchor: tuple[float, float]


ANALYTE_WINDOW = PeakWindow(
    integration=(39.2, 44.8),
    pre_anchor=(36.4, 38.0),
    post_anchor=(46.0, 47.6),
)
INTERNAL_STANDARD_WINDOW = PeakWindow(
    integration=(54.4, 59.6),
    pre_anchor=(52.0, 53.6),
    post_anchor=(60.4, 62.0),
)


def _anchor(
    times: list[float],
    signals: list[float],
    bounds: tuple[float, float],
) -> tuple[float, float]:
    selected = [
        signal
        for time, signal in zip(times, signals)
        if bounds[0] <= time <= bounds[1]
    ]
    if len(selected) < 2:
        raise ValueError("not enough points in baseline anchor")
    center = 0.5 * (bounds[0] + bounds[1])
    return center, sum(selected) / len(selected)


def _linear_baseline(
    time: float,
    start: tuple[float, float],
    end: tuple[float, float],
) -> float:
    if end[0] <= start[0]:
        raise ValueError("baseline anchors must be ordered")
    fraction = (time - start[0]) / (end[0] - start[0])
    return start[1] + fraction * (end[1] - start[1])


def integrate_peak(
    times: list[float],
    signals: list[float],
    window: PeakWindow,
) -> dict[str, float | int]:
    """Integrate baseline-corrected signal with the trapezoid rule."""
    if len(times) != len(signals):
        raise ValueError("time and signal lengths differ")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("time values must be strictly increasing")
    start = _anchor(times, signals, window.pre_anchor)
    end = _anchor(times, signals, window.post_anchor)
    points = [
        (time, signal)
        for time, signal in zip(times, signals)
        if window.integration[0] <= time <= window.integration[1]
    ]
    if len(points) < 3:
        raise ValueError("not enough points in integration window")
    corrected = [
        (time, signal - _linear_baseline(time, start, end))
        for time, signal in points
    ]
    area = 0.0
    for (t0, s0), (t1, s1) in zip(corrected, corrected[1:]):
        area += 0.5 * (s0 + s1) * (t1 - t0)
    return {
        "area": area,
        "baseline_start_mv": _linear_baseline(points[0][0], start, end),
        "baseline_end_mv": _linear_baseline(points[-1][0], start, end),
        "points": len(points),
    }
