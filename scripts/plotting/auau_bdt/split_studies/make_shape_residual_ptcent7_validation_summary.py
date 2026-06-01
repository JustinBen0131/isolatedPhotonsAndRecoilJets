#!/usr/bin/env python3
"""Make a slide-style validation summary for shape-residual BDT sidecars."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import math

import matplotlib.pyplot as plt


BASE_PRODUCT = "globalEtCent1535_bdt_noIso_ptCent7"

PRODUCT_LABELS = {
    BASE_PRODUCT: r"Reference 8 $p_{T}$ x 7 centrality BDT",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeResiduals": "Full feature list + residual ratios",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateDiag": "Full feature list + template distance",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateAll": "Full feature list + residuals + template distance",
    "baseV3E_w33_cent_ptCent7_shapeTemplateDiag": "Compact base v3E+w33 + template distance",
    "baseV3E_w33_cent_ptCent7_shapeTemplateAll": "Compact base v3E+w33 + residuals + template distance",
}

PLOT_LABELS = {
    BASE_PRODUCT: "Reference\n8 $p_{T}$ x 7 centrality",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeResiduals": "+ residual ratios\nfull feature list",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateDiag": "+ template distance\nfull feature list",
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateAll": "+ residuals\n+ template distance\nfull feature list",
    "baseV3E_w33_cent_ptCent7_shapeTemplateDiag": "compact base v3E+w33\n+ template distance",
    "baseV3E_w33_cent_ptCent7_shapeTemplateAll": "compact base v3E+w33\n+ residuals\n+ template distance",
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

PRODUCT_COLORS = {
    BASE_PRODUCT: ("#6B7280", "#4B5563", "#F9FAFB"),
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateAll": ("#16A34A", "#047857", "#F0FDF4"),
    "globalEtCent1535_bdt_noIso_ptCent7_shapeTemplateDiag": ("#2563EB", "#1D4ED8", "#EFF6FF"),
    "globalEtCent1535_bdt_noIso_ptCent7_shapeResiduals": ("#F97316", "#C2410C", "#FFF7ED"),
    "baseV3E_w33_cent_ptCent7_shapeTemplateAll": ("#14B8A6", "#0F766E", "#F0FDFA"),
    "baseV3E_w33_cent_ptCent7_shapeTemplateDiag": ("#64748B", "#475569", "#F8FAFC"),
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


def order_for_slide22(rows: list[dict]) -> list[dict]:
    ref = [row for row in rows if row["product"] == BASE_PRODUCT]
    variants = [row for row in rows if row["product"] != BASE_PRODUCT]
    variants.sort(key=lambda row: row["auc_0_20"], reverse=True)
    return ref + variants


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
    rows = order_for_slide22(rows)
    fig, ax = plt.subplots(figsize=(15.8, 8.9), dpi=220)
    fig.patch.set_facecolor("white")

    y = list(range(len(rows)))
    values = [row["auc_0_20"] for row in rows]
    finite_values = [value for value in values if math.isfinite(value)]
    ref_auc = next((row["auc_0_20"] for row in rows if row["product"] == BASE_PRODUCT), finite_values[0])
    xmin = max(0.0, ref_auc - 0.0032)
    xmax = min(1.0, max(finite_values) + 0.0048)

    for yy, row in zip(y, rows):
        color, edge, fill = PRODUCT_COLORS.get(row["product"], ("#94A3B8", "#64748B", "#F8FAFC"))
        if row["product"] != BASE_PRODUCT:
            ax.axhspan(yy - 0.43, yy + 0.43, color=fill, alpha=0.78, zorder=0)
            ax.hlines(yy, ref_auc, row["auc_0_20"], color=color, linewidth=11, alpha=0.33, zorder=2)
        ax.scatter(
            row["auc_0_20"],
            yy,
            s=300 if row["product"] != BASE_PRODUCT else 240,
            color=color,
            edgecolor=edge,
            linewidth=2.2,
            zorder=4,
        )

    labels = [PLOT_LABELS.get(row["product"], row["label"]) for row in rows]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=15.5)
    ax.invert_yaxis()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(len(rows) - 0.35, -0.65)
    ax.set_xlabel("0-20% centrality validation AUC", fontsize=20.5)
    ax.tick_params(axis="x", labelsize=15.8, pad=7)
    ax.tick_params(axis="y", length=0, pad=18)
    ax.grid(axis="x", color="#D1D5DB", linewidth=1.0, alpha=0.8)
    ax.set_axisbelow(True)
    ax.axvline(ref_auc, color="#4B5563", linewidth=2.0, linestyle=(0, (4, 6)), alpha=0.95, zorder=1)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)

    for yy, val, row in zip(y, values, rows):
        delta = val - ref_auc
        pct = 100.0 * delta / ref_auc if ref_auc else 0.0
        color, edge, _ = PRODUCT_COLORS.get(row["product"], ("#94A3B8", "#64748B", "#F8FAFC"))
        badge_x = ref_auc - 0.00075
        ax.text(
            badge_x,
            yy,
            f"{row['model_count']} BDTs\n{row['feature_count']} inputs",
            ha="right",
            va="center",
            fontsize=11.9,
            color="#4B5563",
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#CBD5E1",
                "linewidth": 1.0,
                "alpha": 0.97,
            },
            zorder=5,
        )
        text_x = min(val + 0.00055, xmax - 0.0013)
        ax.text(
            text_x,
            yy - 0.10,
            f"AUC {val:.4f}",
            va="center",
            ha="left",
            fontsize=18.0,
            color="#111827",
            fontweight="bold",
            zorder=5,
        )
        gain_text = "reference" if row["product"] == BASE_PRODUCT else f"{pct:+.2f}% gain in AUC"
        ax.text(
            text_x,
            yy + 0.22,
            gain_text,
            va="center",
            ha="left",
            fontsize=14.4,
            color="#4B5563" if row["product"] == BASE_PRODUCT else color,
            zorder=5,
        )

    fig.text(
        0.255,
        0.955,
        "AUC gain from shape-residual inputs (0-20%)",
        ha="left",
        va="top",
        fontsize=24.0,
        fontweight="bold",
    )
    fig.text(
        0.255,
        0.902,
        r"Photon12+20 signal & Jet12+20+30 inclusive background; 15 < $p_{T}$ < 35 GeV",
        ha="left",
        va="top",
        fontsize=15.3,
        color="#4B5563",
    )
    fig.text(
        0.255,
        0.870,
        r"All rows use 8 $p_{T}$ x 7 centrality routing; gains are relative to the current 56-BDT reference",
        ha="left",
        va="top",
        fontsize=13.9,
        color="#4B5563",
    )
    fig.text(0.817, 0.918, "sPHENIX", ha="left", va="top", fontsize=20.0, fontstyle="italic", fontweight="bold")
    fig.text(0.903, 0.918, " Internal", ha="left", va="top", fontsize=20.0)
    fig.text(
        0.255,
        0.055,
        "0-20% centrality is the ranking metric. Template-distance inputs drive the visible gain; residual ratios alone are nearly neutral.",
        ha="left",
        va="bottom",
        fontsize=11.8,
        color="#6B7280",
    )

    fig.tight_layout(rect=[0.065, 0.115, 0.975, 0.795])
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
    plot_rows = order_for_slide22(rows)
    write_csv(plot_rows, csv_path)
    draw(plot_rows, png_path)
    print(png_path)
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
