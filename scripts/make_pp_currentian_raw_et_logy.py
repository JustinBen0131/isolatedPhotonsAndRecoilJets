#!/usr/bin/env python3
"""Make a standalone log-y raw cluster-ET spectrum panel for pp current-IAN QA."""

from __future__ import annotations

import argparse
import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_BASE = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_inSituStitch_20260522_0143"
    / "validation/ppg12_exact_reweight_closure"
)
DEFAULT_OUT = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_inSituStitch_20260522_0143"
    / "slide_assets/pp_currentIAN_raw_clusterEt_signal_background_logy.png"
)

SIGNAL_COLOR = "#009E73"
BACKGROUND_COLOR = "#4B5563"
INK = "#111827"
GRID = "#E5E7EB"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, default=DEFAULT_BASE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return ap.parse_args()


def density(rows: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grouped = (
        rows.groupby(["bin_low", "bin_high"], as_index=False)["n_rows"]
        .sum()
        .sort_values("bin_low")
    )
    lows = grouped["bin_low"].to_numpy(float)
    highs = grouped["bin_high"].to_numpy(float)
    vals = grouped["n_rows"].to_numpy(float)
    widths = highs - lows
    total = float(vals.sum())
    dens = np.divide(vals, total * widths, out=np.zeros_like(vals), where=(total > 0) & (widths > 0))
    return lows, highs, dens


def step_hist(ax, lows: np.ndarray, highs: np.ndarray, vals: np.ndarray, *, color: str, label: str) -> None:
    if len(lows) == 0:
        return
    x = np.r_[lows, highs[-1]]
    y = np.r_[vals, vals[-1]]
    ax.step(x, y, where="post", color=color, lw=3.0, label=label)
    centers = 0.5 * (lows + highs)
    mask = vals > 0
    ax.scatter(centers[mask], vals[mask], s=34, color=color, edgecolor="white", linewidth=0.55, zorder=4)


def read_inventory_csv(path: Path) -> pd.DataFrame:
    text = path.read_text()
    header = "source_sample,class,axis,bin_low,bin_high,n_rows,sum_ppg12_exact_weight"
    start = text.find(header)
    if start < 0:
        raise SystemExit(f"Could not find inventory CSV header in {path}")
    return pd.read_csv(io.StringIO(text[start:]))


def main() -> int:
    args = parse_args()
    csv_path = args.base / "ppg12_exact_sample_inventory_binned.csv"
    if not csv_path.is_file():
        raise SystemExit(f"Missing binned inventory: {csv_path}")

    df = read_inventory_csv(csv_path)
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

    fig, ax = plt.subplots(figsize=(8.6, 6.1), dpi=220)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    positive = []
    for cls, color, label in ((1, SIGNAL_COLOR, "Signal"), (0, BACKGROUND_COLOR, "Background")):
        lows, highs, vals = density(et[et["class"] == cls])
        step_hist(ax, lows, highs, vals, color=color, label=label)
        positive.extend(vals[vals > 0])

    ymin = max(min(positive) * 0.55, 1e-7) if positive else 1e-6
    ymax = max(positive) * 2.1 if positive else 1.0
    ax.set_yscale("log")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(15, 35)
    ax.set_xlabel(r"cluster $E_T$ [GeV]", fontsize=19)
    ax.set_ylabel("Area-normalized density", fontsize=19)
    ax.set_title(r"Raw pp training $E_T$ spectra before reweighting", fontsize=22, fontweight="bold", pad=12)
    ax.grid(True, which="major", color=GRID, lw=0.95)
    ax.grid(True, which="minor", color=GRID, lw=0.45, alpha=0.6)
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=16, length=6)
    ax.tick_params(which="minor", length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.15)
        spine.set_color(INK)

    ax.add_patch(
        Rectangle(
            (0.043, 0.835),
            0.42,
            0.13,
            transform=ax.transAxes,
            facecolor="white",
            edgecolor="none",
            alpha=0.96,
            zorder=5,
        )
    )
    ax.text(
        0.055,
        0.93,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=16,
        color=INK,
        va="top",
        zorder=6,
    )
    ax.text(
        0.055,
        0.875,
        "pp current-IAN base-v3E training rows",
        transform=ax.transAxes,
        fontsize=13.5,
        color="#374151",
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
