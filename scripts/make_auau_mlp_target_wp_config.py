#!/usr/bin/env python3
"""Inject AuAu tight-MLP model paths and validation-derived WPs into a YAML."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--template", type=Path, required=True)
    ap.add_argument("--working-points", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--model-dir", type=Path, default=None)
    return ap.parse_args()


def replace_or_append(text: str, key: str, value: str) -> str:
    line = f"{key}: {value}"
    pat = re.compile(rf"^{re.escape(key)}\s*:.*$", re.MULTILINE)
    if pat.search(text):
        return pat.sub(line, text)
    if text and not text.endswith("\n"):
        text += "\n"
    return text + line + "\n"


def main() -> int:
    args = parse_args()
    text = args.template.read_text()
    wp = json.loads(args.working_points.read_text())
    entries = wp.get("runtime_entries", [])
    if not entries:
        raise SystemExit(f"No runtime_entries found in {args.working_points}")
    quoted = "[" + ", ".join(json.dumps(x) for x in entries) + "]"
    text = re.sub(r"^auau_tight_mlp_working_point_entries\s*:.*\n?", "", text, flags=re.MULTILINE)
    text = replace_or_append(text, "auau_tight_mlp_working_point_entries", quoted)
    if args.model_dir is not None:
        text = replace_or_append(text, "auau_tight_mlp_model_dir", str(args.model_dir))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(f"[OK] wrote MLP target-WP YAML: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
