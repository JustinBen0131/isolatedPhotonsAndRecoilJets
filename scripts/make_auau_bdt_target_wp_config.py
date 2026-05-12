#!/usr/bin/env python3
"""Generate a RecoilJets YAML with validation-derived AuAu BDT working points."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--template", type=Path, required=True)
    ap.add_argument("--working-points", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--model-dir", type=Path, default=None)
    ap.add_argument(
        "--product-map",
        action="append",
        default=[],
        help="Variant/product mapping, e.g. auauEtFineCent3BDT=ptFine_cent3. May be repeated or comma-separated.",
    )
    ap.add_argument("--allow-global-fallback", action="store_true")
    return ap.parse_args()


def runtime_entry(variant: str, wp: dict) -> str:
    if wp.get("status") != "ready":
        raise ValueError(f"product working point is not ready: {wp.get('reason', 'unknown reason')}")
    if wp.get("mode") == "grid2d":
        pt_edges = ";".join(f"{float(x):.8g}" for x in wp.get("pt_edges", []))
        cent_edges = ";".join(f"{float(x):.8g}" for x in wp.get("cent_edges", []))
        vals = []
        for row in wp.get("grid_thresholds", []):
            vals.extend(float(x) for x in row)
        thresholds = ";".join(f"{x:.8g}" for x in vals)
        return (
            f"{variant}|grid2d|{pt_edges}|{cent_edges}|{thresholds}|"
            f"{float(wp['pt_min']):.8g}|{float(wp['pt_max']):.8g}|{float(wp.get('max_score', 1.0)):.8g}"
        )
    if wp.get("mode") == "binned":
        edges = ";".join(f"{float(x):.8g}" for x in wp.get("binned_edges", []))
        thresholds = ";".join(f"{float(x):.8g}" for x in wp.get("binned_thresholds", []))
        return (
            f"{variant}|binned|{edges}|{thresholds}|"
            f"{float(wp['pt_min']):.8g}|{float(wp['pt_max']):.8g}|{float(wp.get('max_score', 1.0)):.8g}"
        )
    return (
        f"{variant}|linear|{float(wp['intercept']):.8g}|{float(wp.get('slope', 0.0)):.8g}|"
        f"{float(wp['pt_min']):.8g}|{float(wp['pt_max']):.8g}|{float(wp.get('max_score', 1.0)):.8g}"
    )


def parse_product_map(items: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for item in items:
        for token in item.split(","):
            token = token.strip()
            if not token:
                continue
            if "=" not in token:
                raise SystemExit(f"Bad --product-map token, expected variant=product: {token}")
            variant, product = token.split("=", 1)
            mapping[variant.strip()] = product.strip()
    if not mapping:
        raise SystemExit("At least one --product-map variant=product entry is required")
    return mapping


def replace_model_dir(text: str, model_dir: Path | None) -> str:
    if model_dir is None:
        return text
    # Campaign templates keep a single expanded-model directory key. Replace the
    # value only, preserving comments and any other YAML content.
    replacement = f"auau_tight_bdt_expanded_model_dir: {model_dir}"
    text = re.sub(r"^auau_tight_bdt_expanded_model_dir\s*:.*$", replacement, text, flags=re.MULTILINE)
    return text


def main() -> int:
    args = parse_args()
    if not args.template.is_file():
        raise SystemExit(f"Template YAML does not exist: {args.template}")
    if not args.working_points.is_file():
        raise SystemExit(f"Working-point JSON does not exist: {args.working_points}")

    manifest = json.loads(args.working_points.read_text())
    products = manifest.get("products", {})
    mapping = parse_product_map(args.product_map)
    entries = []
    resolved = {}
    for variant, product in mapping.items():
        if product not in products:
            raise SystemExit(f"Product '{product}' from map '{variant}={product}' not found in {args.working_points}")
        entry = runtime_entry(variant, products[product])
        entries.append(entry)
        resolved[variant] = {"product": product, "entry": entry, "working_point": products[product]}

    text = replace_model_dir(args.template.read_text(), args.model_dir)
    text = re.sub(r"^auau_tight_bdt_working_point_entries\s*:.*\n?", "", text, flags=re.MULTILINE)
    if "auau_tight_bdt_min_intercept" in text and not args.allow_global_fallback:
        text += (
            "\n# Generated from validation-derived per-model target-efficiency cuts.\n"
            "# If a listed AuAu BDT tight mode is active, this overrides the global flat cut.\n"
        )
    quoted = ", ".join(json.dumps(x) for x in entries)
    text += f"auau_tight_bdt_working_point_entries: [{quoted}]\n"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    args.out.with_suffix(args.out.suffix + ".working_points.json").write_text(
        json.dumps(
            {
                "schema": "RECOILJETS_AUAU_BDT_RUNTIME_WP_CONFIG_V1",
                "template": str(args.template),
                "working_points": str(args.working_points),
                "model_dir": str(args.model_dir) if args.model_dir else None,
                "target_signal_efficiency": manifest.get("target_signal_efficiency"),
                "resolved": resolved,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"generated_yaml={args.out}")
    print(f"runtime_entries={len(entries)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
