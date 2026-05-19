#!/usr/bin/env python3
"""Make 2x3 score separation plot for base-v3E centrality vs +3x3-width BDT."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


PRODUCTS = [
    ("baseBDT_v3E_withCentrality", "Base PPG12 baseV3E\n+ centrality"),
    ("baseBDT_v3E_withCentrality_w33", "Only change:\nadd weta33 and wphi33"),
]
CENT_BINS = ["0-20", "20-50", "50-80"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hist-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="basev3e_cent_vs_w33_score_split_2x3")
    return parser.parse_args()


def step(ax, edges: np.ndarray, density: np.ndarray, color: str, label: str) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", lw=2.0, color=color, label=label)
    ax.fill_between(edges, y, step="post", color=color, alpha=0.08, linewidth=0)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.hist_csv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    sig_color = "#1f77b4"
    bkg_color = "#ff7f0e"

    prepared: dict[tuple[str, str], dict[str, object]] = {}
    max_y = 0.0
    for product, _ in PRODUCTS:
        for cent in CENT_BINS:
            part = df[(df["product"] == product) & (df["centrality_bin"] == cent)].sort_values("score_low")
            if part.empty:
                raise SystemExit(f"Missing rows for {product} {cent}")
            edges = np.r_[part["score_low"].to_numpy(float), float(part["score_high"].iloc[-1])]
            widths = part["score_high"].to_numpy(float) - part["score_low"].to_numpy(float)
            sig = part["signal_count"].to_numpy(float)
            bkg = part["background_count"].to_numpy(float)
            sig_total = float(sig.sum())
            bkg_total = float(bkg.sum())
            sig_density = np.divide(sig, sig_total * widths, out=np.zeros_like(sig), where=(sig_total * widths) > 0)
            bkg_density = np.divide(bkg, bkg_total * widths, out=np.zeros_like(bkg), where=(bkg_total * widths) > 0)
            max_y = max(max_y, float(np.nanmax(sig_density)), float(np.nanmax(bkg_density)))
            prepared[(product, cent)] = {
                "edges": edges,
                "signal_density": sig_density,
                "background_density": bkg_density,
                "auc": float(part["auc"].iloc[0]),
            }

    fig, axes = plt.subplots(2, 3, figsize=(13.8, 6.75), dpi=240, sharex=True, sharey=True)
    for row_i, (product, row_label) in enumerate(PRODUCTS):
        for col_i, cent in enumerate(CENT_BINS):
            ax = axes[row_i, col_i]
            cell = prepared[(product, cent)]
            step(ax, cell["edges"], cell["signal_density"], sig_color, "Signal")
            step(ax, cell["edges"], cell["background_density"], bkg_color, "Background")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, max_y * 1.16)
            ax.grid(True, color="0.90", linewidth=0.55)
            ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=10.5)
            for spine in ax.spines.values():
                spine.set_linewidth(1.0)
            if row_i == 0:
                ax.set_title(f"{cent}%", fontsize=15, fontweight="bold", pad=8)
            ax.text(
                0.055,
                0.90,
                f"AUC {cell['auc']:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=12,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.82", alpha=0.92),
            )
            if col_i == 0:
                ax.text(
                    0.055,
                    0.60 if row_i == 0 else 0.69,
                    row_label,
                    transform=ax.transAxes,
                    ha="left",
                    va="center",
                    rotation=0,
                    fontsize=10.8,
                    fontweight="bold",
                    color="0.10",
                    bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="0.78", alpha=0.94),
                )

    axes[0, 2].legend(loc="upper right", frameon=False, fontsize=12)
    for ax in axes[-1, :]:
        ax.set_xlabel("BDT score", fontsize=13)
    fig.text(
        0.025,
        0.49,
        "Area-normalized density",
        rotation=90,
        ha="center",
        va="center",
        fontsize=12.8,
    )

    fig.suptitle("Base v3E BDT score separation", fontsize=17, fontweight="bold", y=0.978)
    fig.text(
        0.52,
        0.931,
        r"Photon12+20 signal vs Jet12+20+30 background, $15<E_T<35$ GeV",
        ha="center",
        va="top",
        fontsize=10.8,
        color="0.32",
    )
    fig.tight_layout(rect=[0.055, 0.06, 0.995, 0.895], w_pad=1.0, h_pad=0.95)

    png = args.outdir / f"{args.prefix}.png"
    csv_out = args.outdir / f"{args.prefix}.csv"
    fig.savefig(png)
    plt.close(fig)
    df.to_csv(csv_out, index=False)
    print(png)
    print(csv_out)


if __name__ == "__main__":
    main()
