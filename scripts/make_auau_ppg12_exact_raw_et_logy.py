#!/usr/bin/env python3
"""Make a standalone log-y raw cluster-ET spectrum panel for Au+Au PPG12-exact QA."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_BASE = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439"
    / "slideReady/ppg12_exact_reweight_bdt_remote"
)
DEFAULT_OUT = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439"
    / "slideReady/ppg12_exact_reweight_bdt"
    / "auau_ppg12_exact_raw_clusterEt_signal_background_logy.png"
)

SIGNAL_COLOR = "#059669"
BACKGROUND_COLOR = "#EA580C"
INK = "#0F172A"
MUTED = "#475569"
GRID = "#E5E7EB"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, default=DEFAULT_BASE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return ap.parse_args()


def step_hist(ax, rows: pd.DataFrame, *, color: str, label: str) -> np.ndarray:
    rows = rows.sort_values("bin_low")
    lows = rows["bin_low"].to_numpy(float)
    highs = rows["bin_high"].to_numpy(float)
    vals = rows["raw_density"].to_numpy(float)
    if len(lows) == 0:
        return vals
    x = np.r_[lows, highs[-1]]
    y = np.r_[vals, vals[-1]]
    ax.step(x, y, where="post", color=color, lw=2.8, label=label)
    return vals


def main() -> int:
    args = parse_args()
    csv_path = args.base / "ppg12_exact_fine_reweight_histograms.csv"
    if not csv_path.is_file():
        raise SystemExit(f"Missing fine histogram CSV: {csv_path}")

    df = pd.read_csv(csv_path)
    et = df[df["axis"] == "cluster_Et"].copy()
    if et.empty:
        raise SystemExit(f"No cluster_Et rows in {csv_path}")

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

    fig, ax = plt.subplots(figsize=(8.8, 6.2), dpi=220)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    vals = []
    vals.extend(step_hist(ax, et[et["class"] == 1], color=SIGNAL_COLOR, label="Signal"))
    vals.extend(step_hist(ax, et[et["class"] == 0], color=BACKGROUND_COLOR, label="Background"))
    positive = np.asarray([v for v in vals if v > 0], dtype=float)

    ax.set_yscale("log")
    if positive.size:
        ax.set_ylim(max(positive.min() * 0.55, 1e-6), positive.max() * 2.0)
    ax.set_xlim(5, 35)
    ax.set_xlabel(r"cluster $E_T$ [GeV]", fontsize=19)
    ax.set_ylabel("Area-normalized density", fontsize=19)
    ax.set_title(r"Raw Au+Au training $E_T$ spectra before reweighting", fontsize=22, fontweight="bold", pad=12)
    ax.grid(True, which="major", color=GRID, lw=0.95)
    ax.grid(True, which="minor", color=GRID, lw=0.45, alpha=0.62)
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=16, length=6)
    ax.tick_params(which="minor", length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.15)
        spine.set_color(INK)

    ax.add_patch(
        Rectangle(
            (0.055, 0.045),
            0.39,
            0.13,
            transform=ax.transAxes,
            facecolor="white",
            edgecolor="none",
            alpha=0.96,
            zorder=5,
        )
    )
    ax.text(
        0.068,
        0.145,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=16,
        color=INK,
        va="top",
        zorder=6,
    )
    ax.text(
        0.068,
        0.09,
        "Au+Au PPG12-exact training rows",
        transform=ax.transAxes,
        fontsize=13.5,
        color=MUTED,
        va="top",
        zorder=6,
    )
    ax.legend(frameon=False, fontsize=16, loc="upper right", handlelength=2.6)

    fig.tight_layout(pad=1.0)
    fig.savefig(args.out)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
