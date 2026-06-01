#!/usr/bin/env python3
"""Generate a RecoilJets YAML with validation-derived AuAu BDT working points."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable


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


def fmt_yaml_number(x: float) -> str:
    return f"{float(x):.8g}"


def fmt_int3(x: float) -> str:
    return f"{int(round(float(x))):03d}"


def pt_tag(lo: float, hi: float) -> str:
    return f"pt_{fmt_int3(lo)}_{fmt_int3(hi)}"


def cent_tag(lo: float, hi: float) -> str:
    return f"cent_{fmt_int3(lo)}_{fmt_int3(hi)}"


def yaml_list(items: Iterable[object]) -> str:
    return "[" + ", ".join(json.dumps(str(x)) for x in items) + "]"


def yaml_number_list(items: Iterable[float]) -> str:
    return "[" + ", ".join(fmt_yaml_number(float(x)) for x in items) + "]"


def model_path(model_dir: Path, product: str, suffix: str) -> Path:
    return model_dir / f"auau_tight_bdt_{product}_{suffix}_tmva.root"


def metadata_path(model_dir: Path, product: str, suffix: str) -> Path:
    return model_dir / f"auau_tight_bdt_{product}_{suffix}_tmva.metadata.json"


def pt_cent_suffixes(pt_edges: list[float], cent_edges: list[float]) -> list[str]:
    suffixes = []
    # Runtime pt x centrality model lookup is pT-major, centrality-minor.
    for ip in range(len(pt_edges) - 1):
        for ic in range(len(cent_edges) - 1):
            suffixes.append(
                f"{pt_tag(pt_edges[ip], pt_edges[ip + 1])}_{cent_tag(cent_edges[ic], cent_edges[ic + 1])}"
            )
    return suffixes


def fallback_cent_suffixes(pt_edges: list[float], cent_edges: list[float]) -> list[str]:
    # auauPtCent* modes require one fallback model per centrality bin even for
    # 15-35 GeV target-WP applications. Reuse the highest trained pT interval
    # as the configured fallback; the target-WP cut itself still rejects ET
    # outside the validation range.
    if len(pt_edges) < 2:
        return []
    ip = len(pt_edges) - 2
    return [
        f"{pt_tag(pt_edges[ip], pt_edges[ip + 1])}_{cent_tag(cent_edges[ic], cent_edges[ic + 1])}"
        for ic in range(len(cent_edges) - 1)
    ]


def read_feature_list(model_dir: Path, product: str, suffixes: list[str]) -> list[str]:
    for suffix in suffixes:
        meta = metadata_path(model_dir, product, suffix)
        if not meta.is_file():
            continue
        data = json.loads(meta.read_text())
        features = data.get("features") or data.get("feature_names") or []
        if features:
            return [str(x) for x in features]
    return []


def strip_generated_runtime_keys(text: str) -> str:
    keys = [
        "auau_tight_bdt_pt_bin_edges",
        "auau_tight_bdt_etfine_pt_bin_edges",
        "auau_tight_bdt_cent3_edges",
        "auau_tight_bdt_cent7_edges",
        "auau_tight_bdt_ptCent3_model_files",
        "auau_tight_bdt_ptCent7_model_files",
        "auau_tight_bdt_ptCent3_fallback_model_files",
        "auau_tight_bdt_ptCent7_fallback_model_files",
        "auau_tight_bdt_etFineCent3_model_files",
        "auau_tight_bdt_etFineCent7_model_files",
        "auau_tight_bdt_centDep_features",
    ]
    for key in keys:
        text = re.sub(rf"^{re.escape(key)}\s*:.*\n?", "", text, flags=re.MULTILINE)
    return text


def routed_runtime_yaml(variant: str, product: str, wp: dict, model_dir: Path | None) -> tuple[str, dict]:
    if model_dir is None or wp.get("mode") != "grid2d":
        return "", {}

    pt_edges = [float(x) for x in wp.get("pt_edges", [])]
    cent_edges = [float(x) for x in wp.get("cent_edges", [])]
    if len(pt_edges) < 2 or len(cent_edges) < 2:
        return "", {}

    suffixes = pt_cent_suffixes(pt_edges, cent_edges)
    model_files = [str(model_path(model_dir, product, suffix)) for suffix in suffixes]
    missing = [path for path in model_files if not Path(path).is_file()]
    if missing:
        preview = "\n  ".join(missing[:5])
        raise SystemExit(
            f"Model files needed for {variant}={product} are missing ({len(missing)} total):\n  {preview}"
        )

    cent_key = "auau_tight_bdt_cent7_edges" if len(cent_edges) == 8 else "auau_tight_bdt_cent3_edges"
    lines = [
        "",
        "# Generated routed-model metadata for validation-derived target-WP cuts.",
        f"auau_tight_bdt_pt_bin_edges: {yaml_number_list(pt_edges)}",
        f"{cent_key}: {yaml_number_list(cent_edges)}",
    ]
    resolved: dict[str, object] = {
        "pt_edges": pt_edges,
        "cent_edges": cent_edges,
        "model_files": model_files,
    }

    if variant == "auauPtCent3BDT":
        lines.append(f"auau_tight_bdt_ptCent3_model_files: {yaml_list(model_files)}")
        fallback_suffixes = fallback_cent_suffixes(pt_edges, cent_edges)
        fallback_files = [str(model_path(model_dir, product, suffix)) for suffix in fallback_suffixes]
        lines.append(f"auau_tight_bdt_ptCent3_fallback_model_files: {yaml_list(fallback_files)}")
        resolved["fallback_model_files"] = fallback_files
    elif variant == "auauPtCent7BDT":
        lines.append(f"auau_tight_bdt_ptCent7_model_files: {yaml_list(model_files)}")
        fallback_suffixes = fallback_cent_suffixes(pt_edges, cent_edges)
        fallback_files = [str(model_path(model_dir, product, suffix)) for suffix in fallback_suffixes]
        lines.append(f"auau_tight_bdt_ptCent7_fallback_model_files: {yaml_list(fallback_files)}")
        resolved["fallback_model_files"] = fallback_files
    elif variant == "auauEtFineCent3BDT":
        lines.append(f"auau_tight_bdt_etfine_pt_bin_edges: {yaml_number_list(pt_edges)}")
        lines.append(f"auau_tight_bdt_etFineCent3_model_files: {yaml_list(model_files)}")
    elif variant == "auauEtFineCent7BDT":
        lines.append(f"auau_tight_bdt_etfine_pt_bin_edges: {yaml_number_list(pt_edges)}")
        lines.append(f"auau_tight_bdt_etFineCent7_model_files: {yaml_list(model_files)}")
    else:
        return "", {}

    features = read_feature_list(model_dir, product, suffixes)
    if features:
        lines.append(f"auau_tight_bdt_centDep_features: {yaml_list(features)}")
        resolved["features"] = features
    return "\n".join(lines) + "\n", resolved


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
    # value when present, preserving comments and any other YAML content. Some
    # lean production templates omit the key, so append it explicitly for routed
    # BDT modes that discover model files from this directory.
    replacement = f"auau_tight_bdt_expanded_model_dir: {model_dir}"
    pattern = re.compile(r"^auau_tight_bdt_expanded_model_dir\s*:.*$", flags=re.MULTILINE)
    if pattern.search(text):
        return pattern.sub(replacement, text)
    return text.rstrip() + "\n" + replacement + "\n"


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
    text = strip_generated_runtime_keys(text)
    text = re.sub(r"^auau_tight_bdt_working_point_entries\s*:.*\n?", "", text, flags=re.MULTILINE)
    if "auau_tight_bdt_min_intercept" in text and not args.allow_global_fallback:
        text += (
            "\n# Generated from validation-derived per-model target-efficiency cuts.\n"
            "# If a listed AuAu BDT tight mode is active, this overrides the global flat cut.\n"
        )
    routed_blocks = []
    routed_resolved = {}
    for variant, product in mapping.items():
        block, info = routed_runtime_yaml(variant, product, products[product], args.model_dir)
        if block:
            routed_blocks.append(block)
            routed_resolved[variant] = info
    if routed_blocks:
        text += "\n".join(routed_blocks)
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
                "routed_runtime": routed_resolved,
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
