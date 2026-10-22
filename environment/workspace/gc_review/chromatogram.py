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


def integrate_peak(
    times: list[float],
    signals: list[float],
    window: PeakWindow,
) -> dict[str, float | int]:
    """Return a placeholder integration over the requested retention window."""
    if len(times) != len(signals):
        raise ValueError("time and signal lengths differ")
    points = [
        (time, signal)
        for time, signal in zip(times, signals)
        if window.integration[0] <= time <= window.integration[1]
    ]
    if len(points) < 3:
        raise ValueError("not enough points in integration window")
    area = 0.0
    for (t0, s0), (t1, s1) in zip(points, points[1:]):
        area += 0.5 * (s0 + s1) * (t1 - t0)
    return {
        "area": area,
        "baseline_start_mv": 0.0,
        "baseline_end_mv": 0.0,
        "points": len(points),
    }
