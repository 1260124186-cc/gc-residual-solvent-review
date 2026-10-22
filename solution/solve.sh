#!/bin/bash
set -euo pipefail

install -m 0644 /solution/chromatogram.py /workspace/gc_review/chromatogram.py
install -m 0644 /solution/calibration.py /workspace/gc_review/calibration.py
install -m 0644 /solution/reporting.py /workspace/gc_review/reporting.py
install -m 0644 /solution/cli.py /workspace/gc_review/cli.py
