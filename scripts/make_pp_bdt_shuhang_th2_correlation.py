#!/usr/bin/env python3
"""Make TH2-style same-row score correlation plots for pp baseV3E audit."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm


INK = "#111827"
GRID = "#D1D5DB"
RED = "#D62728"
BLUE = "#2563EB"


def draw_sphenix_label(ax: plt.Axes) -> None:
    ax.text(
        0.04,
        0.96,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontstyle="italic",
        fontweight="bold",
        color=INK,
    )
    ax.text(
        0.245,
        0.96,
        "Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        color=INK,
    )


def make_panel(
    ax: plt.Axes,
    frame: pd.DataFrame,
    y_col: str,
    y_label: str,
    *,
    bins: int,
    title: str,
    show_label: bool,
) -> None:
    x = frame["our_score"].to_numpy(dtype=float)
    y = frame[y_col].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & (x >= 0.0) & (x <= 1.0) & (y >= 0.0) & (y <= 1.0)
    x = x[valid]
    y = y[valid]

    counts, xedges, yedges, image = ax.hist2d(
        x,
        y,
        bins=bins,
        range=[[0.0, 1.0], [0.0, 1.0]],
        norm=LogNorm(vmin=1),
        cmap="viridis",
    )
    ax.plot([0, 1], [0, 1], color="white", lw=2.0, ls="--", alpha=0.9)
    ax.plot([0, 1], [0, 1], color=INK, lw=0.8, ls="--", alpha=0.45)

    # Median response profile makes score compression/calibration differences visible.
    centers = 0.5 * (xedges[:-1] + xedges[1:])
    med = np.full_like(centers, np.nan)
    q16 = np.full_like(centers, np.nan)
    q84 = np.full_like(centers, np.nan)
    bin_idx = np.digitize(x, xedges) - 1
    for i in range(len(centers)):
        vals = y[bin_idx == i]
        if vals.size >= 20:
            med[i] = np.median(vals)
            q16[i], q84[i] = np.quantile(vals, [0.16, 0.84])
    mask = np.isfinite(med)
    ax.fill_between(centers[mask], q16[mask], q84[mask], color="white", alpha=0.18, lw=0)
    ax.plot(centers[mask], med[mask], color="white", lw=2.8, label="median Shuhang score")
    ax.plot(centers[mask], med[mask], color=RED, lw=1.5)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("This analysis frozen XGBoost score", fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", color=INK, pad=10)
    ax.grid(color=GRID, lw=0.7, alpha=0.45)
    ax.tick_params(axis="both", labelsize=12, direction="in", top=True, right=True)
    if show_label:
        draw_sphenix_label(ax)
        ax.text(
            0.04,
            0.84,
            r"$p{+}p$ $\sqrt{s}=200$ GeV" "\n" r"$|\eta|<0.7$, $22<E_T<28$ GeV" "\n"
            f"same-row audit: {x.size:,} rows",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=12.5,
            color=INK,
            bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="#E5E7EB", alpha=0.88),
        )
    return image


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--bins", type=int, default=100)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.csv)
    needed = {"our_score", "shuhang_split_score", "shuhang_nosplit_score"}
    missing = needed.difference(frame.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {args.csv}: {sorted(missing)}")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": INK,
            "axes.linewidth": 1.2,
        }
    )

    # Split model: the primary out-of-the-box baseV3E comparison.
    fig, ax = plt.subplots(figsize=(9.6, 8.2), constrained_layout=True)
    image = make_panel(
        ax,
        frame,
        "shuhang_split_score",
        "Shuhang split TMVA baseV3E score",
        bins=args.bins,
        title="Same-row BDT score correlation",
        show_label=True,
    )
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("rows per TH2 bin", fontsize=13)
    cbar.ax.tick_params(labelsize=11)
    split_out = args.outdir / "pp_basev3e_same_row_shuhang_split_vs_this_analysis_TH2.png"
    fig.savefig(split_out, dpi=220, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    # Side-by-side split/nosplit model diagnostic.
    fig, axes = plt.subplots(1, 2, figsize=(16.5, 7.8), constrained_layout=True)
    image0 = make_panel(
        axes[0],
        frame,
        "shuhang_split_score",
        "Shuhang split TMVA score",
        bins=args.bins,
        title="Split TMVA model",
        show_label=True,
    )
    image1 = make_panel(
        axes[1],
        frame,
        "shuhang_nosplit_score",
        "Shuhang nosplit TMVA score",
        bins=args.bins,
        title="Nosplit TMVA cross-check",
        show_label=False,
    )
    cbar = fig.colorbar(image1, ax=axes.ravel().tolist(), pad=0.015)
    cbar.set_label("rows per TH2 bin", fontsize=13)
    cbar.ax.tick_params(labelsize=11)
    both_out = args.outdir / "pp_basev3e_same_row_shuhang_split_nosplit_vs_this_analysis_TH2.png"
    fig.savefig(both_out, dpi=220, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    print(split_out)
    print(both_out)


if __name__ == "__main__":
    main()
