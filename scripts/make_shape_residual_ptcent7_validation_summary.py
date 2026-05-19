#!/usr/bin/env python3
"""Make a slide-style validation summary for shape-residual BDT sidecars."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


BASE_PRODUCT = "globalEtCent1535_bdt_noIso_ptCent7"

PRODUCT_LABELS = {
    BASE_PRODUCT: "Reference 8 p_{T} x 7 centrality BDT",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeResiduals": "Full feature list + residual ratios",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateDiag": "Full feature list + template distance",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateAll": "Full feature list + residuals + template distance",
    "baseV3E_w33_cent_ptCent7_shapeTemplateDiag": "Compact base v3E+w33 + template distance",
    "baseV3E_w33_cent_ptCent7_shapeTemplateAll": "Compact base v3E+w33 + residuals + template distance",
}

PRODUCT_FAMILY = {
    BASE_PRODUCT: "reference",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeResiduals": "full-list variants",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateDiag": "full-list variants",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateAll": "full-list variants",
    "baseV3E_w33_cent_ptCent7_shapeTemplateDiag": "compact variants",
    "baseV3E_w33_cent_ptCent7_shapeTemplateAll": "compact variants",
}

FAMILY_COLORS = {
    "reference": "#6B7280",
    "full-list variants": "#0072B2",
    "compact variants": "#009E73",
}


def load_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def product_metrics(metrics: dict, product: str) -> dict:
    products = metrics.get("products", {})
    if product not in products:
        raise KeyError(f"Product {product!r} not found in validation metrics")
    return products[product]


def feature_counts_from_registry(path: Path | None) -> dict[str, int]:
    if path is None or not path.is_file():
        return {}
    registry = load_json(path)
    models = registry.get("models", [])
    out: dict[str, int] = {}
    for item in models:
        product = str(item.get("product", ""))
        features = item.get("features") or item.get("feature_names") or []
        if product and features:
            out.setdefault(product, len(features))
    return out


def collect_rows(baseline_metrics: Path, shape_metrics: Path, shape_registry: Path | None) -> list[dict]:
    baseline = load_json(baseline_metrics)
    shape = load_json(shape_metrics)
    feature_counts = feature_counts_from_registry(shape_registry)
    feature_counts.setdefault(BASE_PRODUCT, 32)

    rows: list[dict] = []
    base = product_metrics(baseline, BASE_PRODUCT)
    products = [BASE_PRODUCT] + list(shape.get("products", {}).keys())
    for product in products:
        metrics = base if product == BASE_PRODUCT else product_metrics(shape, product)
        auc_cent = metrics.get("auc_by_centrality", {})
        rows.append(
            {
                "product": product,
                "label": PRODUCT_LABELS.get(product, product),
                "family": PRODUCT_FAMILY.get(product, "shape variants"),
                "feature_count": int(feature_counts.get(product, 0)),
                "model_count": 56,
                "auc_inclusive": float(metrics["auc_inclusive"]),
                "auc_0_20": float(auc_cent.get("0_20", "nan")),
                "auc_20_50": float(auc_cent.get("20_50", "nan")),
                "auc_50_80": float(auc_cent.get("50_80", "nan")),
                "finite_score_fraction": float(metrics.get("finite_score_fraction", "nan")),
                "signal_score_mean": float(metrics.get("signal_score_mean", "nan")),
                "background_score_mean": float(metrics.get("background_score_mean", "nan")),
            }
        )
    rows.sort(key=lambda row: row["auc_0_20"], reverse=True)
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    fieldnames = [
        "product",
        "label",
        "family",
        "model_count",
        "feature_count",
        "auc_inclusive",
        "auc_0_20",
        "auc_20_50",
        "auc_50_80",
        "finite_score_fraction",
        "signal_score_mean",
        "background_score_mean",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def draw(rows: list[dict], out: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, ax = plt.subplots(figsize=(15.6, 8.75), dpi=220)
    fig.patch.set_facecolor("white")

    y = list(range(len(rows)))
    values = [row["auc_0_20"] for row in rows]
    colors = [FAMILY_COLORS.get(row["family"], "#BDBDBD") for row in rows]
    labels = [
        f"{row['label']}\n{row['model_count']} BDTs, {row['feature_count']} inputs"
        for row in rows
    ]
    ax.barh(y, values, color=colors, height=0.64)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=12.6)
    ax.invert_yaxis()
    ax.set_xlim(0.76, 0.89)
    ax.set_xlabel("Validation AUC in 0-20% centrality", fontsize=16.5)
    ax.tick_params(axis="x", labelsize=13.8)
    ax.grid(axis="x", color="#D1D5DB", linewidth=1.0, alpha=0.8)
    ax.set_axisbelow(True)

    for yy, val, row in zip(y, values, rows):
        ax.text(
            val + 0.002,
            yy,
            f"0-20% {val:.3f}   incl. {row['auc_inclusive']:.3f}",
            va="center",
            ha="left",
            fontsize=13.3,
            fontweight="bold" if yy == 0 else "normal",
        )

    handles = [
        Line2D([0], [0], color=FAMILY_COLORS[name], lw=9, label=name)
        for name in ("reference", "full-list variants", "compact variants")
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        frameon=True,
        facecolor="white",
        edgecolor="#E5E7EB",
        fontsize=12.5,
        title="Model family",
        title_fontsize=13.2,
    )

    fig.text(0.052, 0.962, "Shape-residual routed BDT validation", ha="left", va="top", fontsize=24, fontweight="bold")
    fig.text(
        0.052,
        0.919,
        r"Photon12+20 & Jet12+20+30, 15 < $p_{T}$ < 35 GeV; 8 $p_{T}$ x 7 centrality routing",
        ha="left",
        va="top",
        fontsize=15.0,
        color="#374151",
    )
    fig.text(
        0.052,
        0.885,
        "All variants avoid isolation inputs; ranked by central 0-20% AUC because that is the hardest AuAu region",
        ha="left",
        va="top",
        fontsize=13.4,
        color="#111827",
    )
    fig.text(0.835, 0.958, "sPHENIX", ha="left", va="top", fontsize=19, fontstyle="italic", fontweight="bold")
    fig.text(0.918, 0.958, " Internal", ha="left", va="top", fontsize=19)

    fig.tight_layout(rect=[0.045, 0.055, 0.985, 0.825])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-metrics", type=Path, required=True)
    parser.add_argument("--shape-metrics", type=Path, required=True)
    parser.add_argument("--shape-registry", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tag", default="shape_residual_ptcent7_validation_slide31_style")
    args = parser.parse_args()

    rows = collect_rows(args.baseline_metrics, args.shape_metrics, args.shape_registry)
    csv_path = args.outdir / f"{args.tag}.csv"
    png_path = args.outdir / f"{args.tag}.png"
    write_csv(rows, csv_path)
    draw(rows, png_path)
    print(png_path)
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
