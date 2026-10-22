"""Deliberately empty response used to probe the scoring boundary."""

from __future__ import annotations

import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    arguments = sys.argv[1:] if argv is None else argv
    if len(arguments) == 3:
        Path(arguments[2]).write_text("{}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
