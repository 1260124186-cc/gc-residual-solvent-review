"""Deliberately hard-coded response used to probe the scoring boundary."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    arguments = sys.argv[1:] if argv is None else argv
    report = {
        "calibration": {
            "slope": 0.70013853923367,
            "intercept": 0.02994759170186545,
            "weighted_r2": 0.9999999498511366,
            "sdyx": 0.000133717,
            "lod_mg_l": 0.000630,
            "loq_mg_l": 0.0019097618292766051,
            "accepted": True,
            "standards": [],
        },
        "blank_mean_mg_l": 0.0,
        "blank_correction_mg_l": 0.0,
        "qc": {
            "calibration": "pass",
            "blank": "pass",
            "duplicate": "pass",
            "spike_recovery": "pass",
        },
        "samples": {
            "S1": {
                "concentration_mg_l": 6.062404561166957,
                "replicates_mg_l": [6.003335797479239, 6.121473324854676],
                "duplicate_rpd_percent": 1.948690922611359,
                "spike_recovery_percent": 99.49300648716657,
                "status": "pass",
            },
            "S2": {
                "concentration_mg_l": 0.0,
                "replicates_mg_l": [0.0, 0.0],
                "duplicate_rpd_percent": 0.0,
                "spike_recovery_percent": 100.04914505935027,
                "status": "below_loq",
            },
        },
        "overall_status": "pass",
    }
    if len(arguments) == 3:
        Path(arguments[2]).write_text(json.dumps(report), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
