#!/usr/bin/env python3
"""Print the public tree names, branch counts, and entry counts."""

from __future__ import annotations

import argparse
from pathlib import Path

import uproot


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    with uproot.open(args.input) as root:
        for name in root.keys(cycle=False):
            obj = root[name]
            if hasattr(obj, "num_entries"):
                print(f"{name:16s} entries={obj.num_entries:8d} branches={len(obj.keys()):4d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
