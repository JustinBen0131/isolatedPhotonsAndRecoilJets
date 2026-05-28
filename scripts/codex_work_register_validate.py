#!/usr/bin/env python3
"""Validate the ThesisAnalysis Codex work register."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from codex_work_register_common import DEFAULT_REGISTER, RegisterError, load_register, validate_register


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("register", nargs="?", default=str(DEFAULT_REGISTER))
    args = parser.parse_args()

    try:
        data = load_register(Path(args.register))
        errors = validate_register(data)
    except (RegisterError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    workstream_count = len(data.get("workstreams") or [])
    print(f"OK: {args.register} validates with {workstream_count} workstreams")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
