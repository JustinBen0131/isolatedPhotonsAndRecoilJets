#!/usr/bin/env python3
"""Regenerate the 3x3-vs-full-width BDT comparison plots."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_FULL_REPORT = Path(
    "dataOutput/auauTightBDTValidation/model_validation_condor_20260509_192942"
)
DEFAULT_3X3_REPORT = Path(
    "dataOutput/auauTightBDTValidation/model_validation_condor_cent3x3_20260510_221419"
)
DEFAULT_OUT_DIR = DEFAULT_3X3_REPORT / "money_plots_3x3"

FULL_PRODUCT = "centAsFeat_pt5to40"
FULL_PRODUCT_FALLBACK = "centAsFeat_allRange"
THREEX3_PRODUCT = "centAsFeat3x3_pt5to40"

SCOPES = [
    ("inclusive", "Inclusive"),
    ("cent_0_20", "0-20%"),
    ("cent_20_50", "20-50%"),
    ("cent_50_80", "50-80%"),
]

CENTRALITY_KEYS = [
    ("0_20", "0-20% central"),
    ("20_50", "20-50% mid-central"),
    ("50_80", "50-80% peripheral"),
]

SERIES = [
    ("full", "Full-cluster $w_\\eta/w_\\phi$", "#2C7FB8"),
    ("3x3", "3x3-window $w_\\eta/w_\\phi$", "#2CA02C"),
]


def sphenix_fig_label(fig, x=0.61, y=0.94, size=14):
    fig.text(x, y, "sPHENIX", ha="left", va="top", fontsize=size, fontstyle="italic", fontweight="bold")
    fig.text(x + 0.105, y, "Internal", ha="left", va="top", fontsize=size)
    fig.text(x, y - 0.048, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=10)


def style_axes(ax):
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=11)
    for spine in ax.spines.values():
        spine.set_color("black")


def product_from_threshold_table(report: Path, product: str, fallback: str | None = None):
    path = report / "validation_threshold_table.csv"
    rows = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row_product = row["product"]
            if row_product != product and (fallback is None or row_product != fallback):
                continue
            if row["target_type"] != "signal_efficiency" or float(row["target"]) != 0.8:
                continue
            if row_product == product or row["scope"] not in rows:
                rows[row["scope"]] = float(row["background_rejection"])
    missing = [scope for scope, _ in SCOPES if scope not in rows]
    if missing:
        raise RuntimeError(f"{path} missing 80% signal rows for {missing}")
    return rows


def plot_background_rejection(full_report: Path, threex3_report: Path, out_dir: Path):
    full = product_from_threshold_table(full_report, FULL_PRODUCT, FULL_PRODUCT_FALLBACK)
    three = product_from_threshold_table(threex3_report, THREEX3_PRODUCT)

    x = np.arange(len(SCOPES))
    width = 0.34
    fig, ax = plt.subplots(figsize=(10.6, 6.6))
    fig.subplots_adjust(left=0.16, right=0.985, top=0.78, bottom=0.16)

    vals = {
        "full": np.asarray([full[scope] for scope, _ in SCOPES]),
        "3x3": np.asarray([three[scope] for scope, _ in SCOPES]),
    }
    for offset, (key, label, color) in zip([-width / 2, width / 2], SERIES):
        bars = ax.bar(x + offset, vals[key], width=width, label=label, color=color)
        for bar, val in zip(bars, vals[key]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                val + 0.00045,
                f"{100.0 * val:.1f}%",
                ha="center",
                va="bottom",
                fontsize=10,
            )

    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in SCOPES], fontsize=12)
    ax.set_ylim(0.795, 0.872)
    ax.set_ylabel("Background rejection at 80% signal efficiency", fontsize=18)
    ax.set_xlabel("Selection scope", fontsize=18, loc="right")
    fig.suptitle("Operating-point check for the centrality-input BDT", fontsize=16, y=0.985)
    ax.legend(loc="upper left", bbox_to_anchor=(0.055, 0.975), frameon=False, fontsize=13)
    sphenix_fig_label(fig, x=0.60, y=0.90, size=13)
    style_axes(ax)

    out = out_dir / "bdt3x3_centfeat_background_rejection_80sig_clean.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def load_roc_product(report: Path, product: str, fallback: str | None = None):
    path = report / "roc_points.json"
    data = json.loads(path.read_text())
    products = data["products"]
    if product in products:
        return products[product]
    if fallback and fallback in products:
        return products[fallback]
    raise RuntimeError(f"{path} missing product {product}")


def plot_roc_by_centrality(full_report: Path, threex3_report: Path, out_dir: Path):
    full = load_roc_product(full_report, FULL_PRODUCT, FULL_PRODUCT_FALLBACK)
    three = load_roc_product(threex3_report, THREEX3_PRODUCT)

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.9), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.73, bottom=0.17, wspace=0.16)

    for ax, (cent_key, cent_label) in zip(axes, CENTRALITY_KEYS):
        for payload, key, label, color in (
            (full, "full", SERIES[0][1], SERIES[0][2]),
            (three, "3x3", SERIES[1][1], SERIES[1][2]),
        ):
            cdata = payload["by_centrality"][cent_key]
            fpr = np.asarray(cdata["background_fake_rate"], dtype=float)
            tpr = np.asarray(cdata["signal_efficiency"], dtype=float)
            ax.plot(fpr, tpr, color=color, linewidth=2.2, label=f"{label}  AUC={cdata['auc']:.3f}")
        ax.plot([0, 1], [0, 1], color="#888888", linewidth=1.0, linestyle="--", zorder=-1)
        ax.set_title(cent_label, fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(True, color="#DDDDDD", linewidth=0.6, alpha=0.65)
        style_axes(ax)

    axes[0].set_ylabel("Truth-photon signal efficiency", fontsize=12)
    for ax in axes:
        ax.set_xlabel("Inclusive-MC background fake rate", fontsize=11)
    for ax in axes:
        ax.legend(loc="lower right", frameon=False, fontsize=8.5)
    fig.suptitle(r"ROC comparison by centrality: full-cluster vs 3x3 $w_\eta/w_\phi$", fontsize=16, y=0.97)
    sphenix_fig_label(fig, x=0.075, y=0.90)

    out = out_dir / "bdt3x3_weta_roc_by_centrality_inclusive_mc.png"
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-report", type=Path, default=DEFAULT_FULL_REPORT)
    parser.add_argument("--threex3-report", type=Path, default=DEFAULT_3X3_REPORT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    bar = plot_background_rejection(args.full_report, args.threex3_report, args.out_dir)
    roc = plot_roc_by_centrality(args.full_report, args.threex3_report, args.out_dir)
    print(bar)
    print(roc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
