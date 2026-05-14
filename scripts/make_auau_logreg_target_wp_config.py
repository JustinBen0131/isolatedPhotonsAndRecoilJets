#!/usr/bin/env python3
"""Inject AuAu tight-logistic-regression artifact and WP entries into YAML."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--template", type=Path, required=True)
    ap.add_argument("--working-points", type=Path, required=True)
    ap.add_argument("--artifact", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--product", default="", help="Product key to keep from the working-point JSON.")
    ap.add_argument("--apply-pt-min", default="15.0")
    ap.add_argument("--apply-pt-max", default="35.0")
    return ap.parse_args()


def replace_or_append(text: str, key: str, value: str) -> str:
    line = f"{key}: {value}"
    pat = re.compile(rf"^{re.escape(key)}\s*:.*$", re.MULTILINE)
    if pat.search(text):
        return pat.sub(line, text)
    if text and not text.endswith("\n"):
        text += "\n"
    return text + line + "\n"


def replace_photon_id_sets(text: str) -> str:
    replacement = "photon_id_sets:\n  - [reference, auauTightLogReg, auauLogRegComplement]\n"
    pat = re.compile(r"^photon_id_sets\s*:\s*\n(?:^[ \t]*-.*\n?)+", re.MULTILINE)
    if pat.search(text):
        return pat.sub(replacement, text)
    if text and not text.endswith("\n"):
        text += "\n"
    return text + replacement


def main() -> int:
    args = parse_args()
    text = args.template.read_text()
    wp = json.loads(args.working_points.read_text())
    products = wp.get("products", {})
    if not products:
        raise SystemExit(f"No products found in {args.working_points}")
    if args.product:
        if args.product not in products:
            raise SystemExit(f"Product {args.product!r} not found in {args.working_points}; available={sorted(products)}")
        selected_products = {args.product: products[args.product]}
    else:
        selected_key = next(iter(products))
        selected_products = {selected_key: products[selected_key]}
    entries = [item["entry"] for item in selected_products.values() if item.get("entry")]
    if len(entries) != 1:
        raise SystemExit(f"Expected exactly one runtime entry after product selection, got {len(entries)}")
    quoted = "[" + ", ".join(json.dumps(x) for x in entries) + "]"
    text = replace_photon_id_sets(text)
    text = replace_or_append(text, "auau_tight_logreg_model_file", str(args.artifact))
    text = replace_or_append(text, "auau_tight_logreg_apply_pt_min", args.apply_pt_min)
    text = replace_or_append(text, "auau_tight_logreg_apply_pt_max", args.apply_pt_max)
    text = replace_or_append(text, "auau_tight_logreg_min_intercept", "0.80")
    text = replace_or_append(text, "auau_tight_logreg_min_slope", "0.0")
    text = replace_or_append(text, "auau_tight_logreg_max", "1.0")
    text = replace_or_append(text, "auau_nontight_logreg_min_intercept", "0.20")
    text = replace_or_append(text, "auau_nontight_logreg_min_slope", "0.0")
    text = replace_or_append(text, "auau_nontight_logreg_max_intercept", "0.80")
    text = replace_or_append(text, "auau_nontight_logreg_max_slope", "0.0")
    text = replace_or_append(text, "auau_tight_logreg_working_point_entries", quoted)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(f"[OK] wrote logreg target-WP YAML: {args.out}")
    print(f"[OK] selected logreg product: {next(iter(selected_products))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
