#!/usr/bin/env python3
"""Make a slide-ready pp current-IAN ET/eta reweighting verification PNG."""

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
import json
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
    / "slide_assets/pp_currentIAN_et_eta_reweighting_verification_slide.png"
)

SIGNAL_COLOR = "#009E73"
BACKGROUND_COLOR = "#6B7280"
INK = "#111827"
MUTED = "#475569"
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


def closure_stats(df: pd.DataFrame, axis: str) -> tuple[float, float]:
    _, _, sig = axis_data(df, axis, 1, True)
    _, _, bkg = axis_data(df, axis, 0, True)
    mask = np.isfinite(sig) & np.isfinite(bkg) & (sig > 0) & (bkg > 0)
    if not np.any(mask):
        return float("nan"), float("nan")
    ratio = sig[mask] / bkg[mask]
    return float(np.nanmax(np.abs(ratio - 1.0))), float(np.sqrt(np.nanmean((ratio - 1.0) ** 2)))


def step_hist(ax, lows, highs, vals, *, color: str, label: str) -> None:
    if len(lows) == 0:
        return
    x = np.r_[lows, highs[-1]]
    y = np.r_[vals, vals[-1]]
    ax.step(x, y, where="post", color=color, lw=2.5, label=label)


def draw_panel(ax, df: pd.DataFrame, axis: str, weighted: bool, title: str, xlabel: str) -> None:
    for cls, color, label in ((1, SIGNAL_COLOR, "Signal"), (0, BACKGROUND_COLOR, "Background")):
        lows, highs, vals = axis_data(df, axis, cls, weighted)
        step_hist(ax, lows, highs, vals, color=color, label=label)
    ax.set_title(title, fontsize=19.5, fontweight="bold", pad=9)
    ax.set_xlabel(xlabel, fontsize=15.8)
    ax.set_ylabel("Area-normalized density", fontsize=15.8)
    ax.grid(True, color=GRID, lw=0.8)
    ax.tick_params(direction="in", top=True, right=True, labelsize=12.8)
    for spine in ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(1.1)
    ax.legend(frameon=False, fontsize=13.8, loc="best")
    if weighted:
        max_dev, rms = closure_stats(df, axis)
        ax.text(
            0.03,
            0.92,
            f"weighted S/B ratio\nmax dev. {100*max_dev:.1f}%, RMS {100*rms:.1f}%",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11.2,
            color=MUTED,
            bbox=dict(boxstyle="round,pad=0.24", fc="white", ec="#CBD5E1", lw=0.8, alpha=0.9),
        )


def inventory_text(meta: dict) -> str:
    inv = meta.get("sample_validation", {}).get("inventory", [])
    sig = sum(int(row.get("n_signal", 0)) for row in inv)
    bkg = sum(int(row.get("n_background", 0)) for row in inv)
    weighting = meta.get("weighting", {})
    class_counts = weighting.get("class_counts", {})
    sig_cap = int(class_counts.get("1", 0))
    bkg_cap = int(class_counts.get("0", 0))
    return (
        f"Training rows before cap: {sig/1e6:.2f}M signal, {bkg/1e6:.2f}M background; "
        f"weighted training cap: {sig_cap/1e6:.1f}M + {bkg_cap/1e6:.1f}M."
    )


def main() -> int:
    args = parse_args()
    csv_path = args.base / "ppg12_exact_sample_inventory_binned.csv"
    meta_path = args.base / "ppg12_exact_reweighting_metadata.json"
    if not csv_path.is_file():
        raise SystemExit(f"Missing binned inventory: {csv_path}")
    if not meta_path.is_file():
        raise SystemExit(f"Missing metadata: {meta_path}")
    df = pd.read_csv(csv_path)
    meta = json.loads(meta_path.read_text())
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
    fig.text(
        0.055,
        0.965,
        r"pp training-prior check: PPG12-style weights flatten $E_T$ and $\eta$ before BDT training",
        fontsize=23.5,
        fontweight="bold",
        ha="left",
        va="top",
        color=INK,
    )
    fig.text(
        0.055,
        0.915,
        "Current-IAN pp baseV3E/no-centrality lane; PhotonJet5/10/20 signal vs Jet8/12/20/30 background.",
        fontsize=13.8,
        ha="left",
        va="top",
        color=MUTED,
    )
    draw_panel(axes[0, 0], df, "cluster_Et", False, r"Raw cluster $E_T$", r"cluster $E_T$ [GeV]")
    draw_panel(axes[0, 1], df, "cluster_Et", True, r"After PPG12-style $E_T$ weight", r"cluster $E_T$ [GeV]")
    draw_panel(axes[1, 0], df, "cluster_Eta", False, r"Raw cluster $\eta$", r"cluster $\eta$")
    draw_panel(axes[1, 1], df, "cluster_Eta", True, r"After PPG12-style $\eta$ weight", r"cluster $\eta$")
    fig.text(
        0.055,
        0.045,
        "Readout: the weighted signal/background priors are deliberately made comparable, so the pp baseV3E BDT is tested on shower-shape separation rather than sample-composition shortcuts.",
        fontsize=13.0,
        ha="left",
        va="bottom",
        color=INK,
    )
    fig.text(0.055, 0.018, inventory_text(meta), fontsize=10.8, ha="left", va="bottom", color=MUTED)
    fig.subplots_adjust(left=0.075, right=0.975, top=0.84, bottom=0.14, wspace=0.18, hspace=0.40)
    fig.savefig(args.out)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
