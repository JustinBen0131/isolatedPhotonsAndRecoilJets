#!/usr/bin/env python3
"""Prepare compact CSV tables for expanded AuAu tight-BDT validation plots."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


LABELS = {
    "centINDcontrol_allRange": "PPG12-like",
    "centAsFeat_allRange": "Centrality input",
    "centAsFeatMinOpt_pt5to40": "Minority optimized",
    "centDepBDTs_allRange": "3 centrality BDTs",
    "centDepFineBDTs_allRange": "7 centrality BDTs",
    "ptBinCentAsFeat": "pT-binned, cent input",
    "ptCentDep3": "pT x 3 centrality",
    "ptCentDepFine": "pT x fine centrality",
}

PLOT_PRODUCTS = list(LABELS)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def fval(value, default=math.nan) -> float:
    try:
        return float(value)
    except Exception:
        return default


def range_label(key: str, suffix: str = "") -> str:
    lo, hi = key.split("_", 1)
    return f"{lo}-{hi}{suffix}"


def centrality_sort_key(key: str) -> float:
    return fval(key.split("_", 1)[0])


def pt_sort_key(key: str) -> float:
    return fval(key.split("_", 1)[0])


def product_label(product: str) -> str:
    return LABELS.get(product, product.replace("_", " "))


def build_product_summary(report_dir: Path, out_dir: Path) -> None:
    rankings = read_csv(report_dir / "validation_model_rankings.csv")
    by_product = {row["product"]: row for row in rankings}
    rows = []
    for idx, product in enumerate(PLOT_PRODUCTS, start=1):
        row = by_product.get(product)
        if not row:
            continue
        sig = fval(row.get("signal_score_mean"))
        bkg = fval(row.get("background_score_mean"))
        rows.append(
            {
                "rank_order": idx,
                "product": product,
                "label": product_label(product),
                "auc_inclusive": fval(row.get("auc_inclusive")),
                "auc_pt_mean": fval(row.get("auc_pt_mean")),
                "auc_cent_mean": fval(row.get("auc_cent_mean")),
                "auc_pt_cent_mean": fval(row.get("auc_pt_cent_mean")),
                "finite_score_fraction": fval(row.get("finite_score_fraction")),
                "eligible_entries": int(fval(row.get("eligible_entries"), 0)),
                "signal_score_mean": sig,
                "background_score_mean": bkg,
                "score_separation": sig - bkg,
            }
        )
    write_csv(
        out_dir / "money_product_summary.csv",
        rows,
        [
            "rank_order",
            "product",
            "label",
            "auc_inclusive",
            "auc_pt_mean",
            "auc_cent_mean",
            "auc_pt_cent_mean",
            "finite_score_fraction",
            "eligible_entries",
            "signal_score_mean",
            "background_score_mean",
            "score_separation",
        ],
    )


def build_diagnostic_tables(report_dir: Path, out_dir: Path) -> None:
    diag = json.loads((report_dir / "validation_deep_diagnostics.json").read_text())
    centrality_rows = []
    pt_rows = []
    threshold_rows = []
    for product in PLOT_PRODUCTS:
        pdata = diag.get("products", {}).get(product)
        if not pdata:
            continue
        for key, cell in sorted(pdata.get("auc_by_centrality", {}).items(), key=lambda kv: centrality_sort_key(kv[0])):
            centrality_rows.append(
                {
                    "product": product,
                    "label": product_label(product),
                    "centrality_bin": key,
                    "centrality_label": range_label(key, "%"),
                    "auc": cell.get("auc", ""),
                    "entries": cell.get("entries", ""),
                    "signal_entries": cell.get("signal_entries", ""),
                    "background_entries": cell.get("background_entries", ""),
                }
            )
        for key, cell in sorted(pdata.get("auc_by_pt", {}).items(), key=lambda kv: pt_sort_key(kv[0])):
            if key == "35_40":
                continue
            pt_rows.append(
                {
                    "product": product,
                    "label": product_label(product),
                    "pt_bin": key,
                    "pt_label": range_label(key, " GeV"),
                    "auc": cell.get("auc", ""),
                    "entries": cell.get("entries", ""),
                    "signal_entries": cell.get("signal_entries", ""),
                    "background_entries": cell.get("background_entries", ""),
                }
            )
        for item in pdata.get("thresholds_inclusive", {}).get("by_signal_efficiency", []):
            target = fval(item.get("target_signal_efficiency"))
            if target not in (0.7, 0.8, 0.9):
                continue
            threshold_rows.append(
                {
                    "product": product,
                    "label": product_label(product),
                    "target_signal_efficiency": target,
                    "threshold": item.get("threshold", ""),
                    "background_fake_rate": item.get("background_fake_rate", ""),
                    "background_rejection": item.get("background_rejection", ""),
                }
            )

    write_csv(
        out_dir / "money_auc_by_centrality.csv",
        centrality_rows,
        ["product", "label", "centrality_bin", "centrality_label", "auc", "entries", "signal_entries", "background_entries"],
    )
    write_csv(
        out_dir / "money_auc_by_pt.csv",
        pt_rows,
        ["product", "label", "pt_bin", "pt_label", "auc", "entries", "signal_entries", "background_entries"],
    )
    write_csv(
        out_dir / "money_thresholds.csv",
        threshold_rows,
        ["product", "label", "target_signal_efficiency", "threshold", "background_fake_rate", "background_rejection"],
    )


def build_feature_separation(report_dir: Path, out_dir: Path) -> None:
    rows = read_csv(report_dir / "validation_feature_summary.csv")
    inclusive = {}
    for row in rows:
        if row.get("scope") != "inclusive":
            continue
        inclusive.setdefault(row["feature"], {})[row["class"]] = row
    out_rows = []
    for feature, by_class in inclusive.items():
        sig = by_class.get("signal")
        bkg = by_class.get("background")
        if not sig or not bkg:
            continue
        ms, mb = fval(sig["mean"]), fval(bkg["mean"])
        ss, sb = fval(sig["std"]), fval(bkg["std"])
        denom = math.sqrt(max(1.0e-12, 0.5 * (ss * ss + sb * sb)))
        out_rows.append(
            {
                "feature": feature,
                "signal_mean": ms,
                "background_mean": mb,
                "separation": abs(ms - mb) / denom,
                "signed_delta": ms - mb,
            }
        )
    out_rows.sort(key=lambda r: r["separation"], reverse=True)
    write_csv(
        out_dir / "money_feature_separation.csv",
        out_rows[:16],
        ["feature", "signal_mean", "background_mean", "separation", "signed_delta"],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report_dir", type=Path)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()
    report_dir = args.report_dir
    out_dir = args.out_dir or report_dir / "money_tables"
    required = [
        "validation_model_rankings.csv",
        "validation_deep_diagnostics.json",
        "validation_feature_summary.csv",
    ]
    for name in required:
        if not (report_dir / name).is_file():
            raise SystemExit(f"Missing required validation file: {report_dir / name}")
    build_product_summary(report_dir, out_dir)
    build_diagnostic_tables(report_dir, out_dir)
    build_feature_separation(report_dir, out_dir)
    print(f"[OK] wrote money plot tables to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
