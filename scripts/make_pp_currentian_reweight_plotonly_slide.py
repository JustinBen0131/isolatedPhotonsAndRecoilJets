#!/usr/bin/env python3
"""Make a plot-only 16:9 pp ET/eta reweighting QA slide."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_BASE = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_fullsim_20260521_1811"
    / "validation/ppg12_exact_reweight_closure"
)
DEFAULT_OUT = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_fullsim_20260521_1811"
    / "slide_assets/pp_currentIAN_et_eta_reweighting_plotonly_slide9_style.png"
)

SIGNAL_COLOR = "#009E73"
BACKGROUND_COLOR = "#6B7280"
INK = "#111827"
GRID = "#E5E7EB"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, default=DEFAULT_BASE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return ap.parse_args()


def density(rows: pd.DataFrame, value_col: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grouped = (
        rows.groupby(["bin_low", "bin_high"], as_index=False)[value_col]
        .sum()
        .sort_values("bin_low")
    )
    lows = grouped["bin_low"].to_numpy(float)
    highs = grouped["bin_high"].to_numpy(float)
    vals = grouped[value_col].to_numpy(float)
    widths = highs - lows
    total = float(vals.sum())
    dens = np.divide(vals, total * widths, out=np.zeros_like(vals), where=(total > 0) & (widths > 0))
    return lows, highs, dens


def axis_data(df: pd.DataFrame, axis: str, cls: int, weighted: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    value_col = "sum_ppg12_exact_weight" if weighted else "n_rows"
    rows = df[(df["axis"] == axis) & (df["class"] == cls)]
    return density(rows, value_col)


def step_hist(ax, lows, highs, vals, *, color: str, label: str) -> None:
    if len(lows) == 0:
        return
    x = np.r_[lows, highs[-1]]
    y = np.r_[vals, vals[-1]]
    ax.step(x, y, where="post", color=color, lw=2.4, label=label)


def draw_panel(ax, df: pd.DataFrame, axis: str, weighted: bool, title: str, xlabel: str) -> None:
    for cls, color, label in ((1, SIGNAL_COLOR, "Signal"), (0, BACKGROUND_COLOR, "Background")):
        lows, highs, vals = axis_data(df, axis, cls, weighted)
        step_hist(ax, lows, highs, vals, color=color, label=label)
    ax.set_title(title, fontsize=20, fontweight="bold", pad=10)
    ax.set_xlabel(xlabel, fontsize=17)
    ax.set_ylabel("Area-normalized density", fontsize=17)
    ax.grid(True, color=GRID, lw=0.8)
    ax.tick_params(direction="in", top=True, right=True, labelsize=14)
    for spine in ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(1.1)
    ax.legend(frameon=False, fontsize=15, loc="best")


def main() -> int:
    args = parse_args()
    csv_path = args.base / "ppg12_exact_sample_inventory_binned.csv"
    if not csv_path.is_file():
        raise SystemExit(f"Missing binned inventory: {csv_path}")
    df = pd.read_csv(csv_path)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(16, 9), dpi=180)
    fig.patch.set_facecolor("white")
    fig.suptitle(
        r"pp current-IAN training prior control: $E_T$ and $\eta$ before/after PPG12-style weights",
        fontsize=24,
        fontweight="bold",
        y=0.958,
    )
    draw_panel(axes[0, 0], df, "cluster_Et", False, r"Raw cluster $E_T$", r"cluster $E_T$ [GeV]")
    draw_panel(axes[0, 1], df, "cluster_Et", True, r"Weighted cluster $E_T$", r"cluster $E_T$ [GeV]")
    draw_panel(axes[1, 0], df, "cluster_Eta", False, r"Raw cluster $\eta$", r"cluster $\eta$")
    draw_panel(axes[1, 1], df, "cluster_Eta", True, r"Weighted cluster $\eta$", r"cluster $\eta$")
    fig.subplots_adjust(left=0.075, right=0.975, top=0.86, bottom=0.09, wspace=0.18, hspace=0.38)
    fig.savefig(args.out)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
