#!/usr/bin/env python3
"""Retired ThesisAnalysis super-heartbeat entrypoint."""

from __future__ import annotations

import argparse
import sys


def retired_message() -> str:
    return (
        "ERROR: `python3 scripts/codex_os_nightly_heartbeat.py nightly` is retired. "
        "Use one of the lane commands instead, for example "
        "`python3 scripts/codex_os_dream.py lane --lane-id status_provenance`."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("nightly", help="retired")
    parser.parse_args()
    print(retired_message(), file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
