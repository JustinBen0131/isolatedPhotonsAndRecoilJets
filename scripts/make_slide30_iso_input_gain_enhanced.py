#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
OUTDIR = ROOT / "slideReady" / "bdt_iso_input_score_comparison"
HIST_CSV = OUTDIR / "binned_bdt_ptcent7_signal_background_score_split_iso_input_comparison_hist.csv"
OUT_PNG = OUTDIR / "binned_bdt_ptcent7_signal_background_score_split_iso_input_comparison_slide30_enhanced.png"

CENTRALITIES = ["0-20%", "20-50%", "50-80%"]
MODEL_ROWS = [
    ("baseline_inputs", "No isolation input"),
    ("baseline_plus_isolation_inputs", "Isolation input included"),
]
SIGNAL_COLOR = "#0072B2"
BACKGROUND_COLOR = "#D55E00"
GAIN_COLOR = "#2A9D8F"


def add_sphenix_label(fig: plt.Figure) -> None:
    fig.text(0.855, 0.965, "sPHENIX", ha="right", va="top", fontsize=15.4, fontstyle="italic", fontweight="bold")
    fig.text(0.945, 0.965, "Internal", ha="right", va="top", fontsize=15.4)


def style_axis(ax: plt.Axes) -> None:
    ax.grid(True, color="0.88", linewidth=0.85)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", direction="in", top=True, right=True, labelsize=9.4, length=4.2, width=1.0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.05)
        spine.set_color("0.24")


def step_fill(ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str, label: str, linestyle: str = "-") -> None:
    ax.fill_between(x, y, step="post", color=color, alpha=0.105, linewidth=0)
    ax.step(x, y, where="post", color=color, linewidth=2.15, linestyle=linestyle, label=label)


def plot_score_panel(
    ax: plt.Axes,
    sub: pd.DataFrame,
    auc: float,
    gain: float | None,
    is_bottom: bool,
    ylabel: str | None,
) -> None:
    x = np.r_[sub["bin_lo"].to_numpy(), sub["bin_hi"].to_numpy()[-1]]
    sig = np.r_[sub["signal_density"].to_numpy(), sub["signal_density"].to_numpy()[-1]]
    bkg = np.r_[sub["background_density"].to_numpy(), sub["background_density"].to_numpy()[-1]]
    step_fill(ax, x, sig, SIGNAL_COLOR, "Signal")
    step_fill(ax, x, bkg, BACKGROUND_COLOR, "Background", "--")
    ymax = max(float(np.nanmax(sig)), float(np.nanmax(bkg)), 1.0)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, ymax * 1.22)
    style_axis(ax)
    ax.text(
        0.045,
        0.90,
        f"AUC {auc:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.3,
        fontweight="bold",
        color="0.18",
        bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="0.78", linewidth=0.7),
    )
    if gain is not None:
        ax.text(
            0.045,
            0.73,
            f"+{100.0 * gain:.1f}% AUC",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9.1,
            fontweight="bold",
            color=GAIN_COLOR,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#B7DED7", linewidth=0.8),
        )
    if is_bottom:
        ax.set_xlabel("BDT score", fontsize=10.8, labelpad=5)
    else:
        ax.tick_params(labelbottom=False)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=10.0, labelpad=8)


def main() -> None:
    df = pd.read_csv(HIST_CSV)
    auc_by = (
        df.groupby(["model_slug", "centrality_bin"], as_index=False)
        .agg(auc=("auc", "first"), signal_entries=("signal_entries", "first"), background_entries=("background_entries", "first"))
        .set_index(["model_slug", "centrality_bin"])
    )
    gains = {
        cent: float(auc_by.loc[("baseline_plus_isolation_inputs", cent), "auc"] - auc_by.loc[("baseline_inputs", cent), "auc"])
        for cent in CENTRALITIES
    }

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.05,
        }
    )
    fig = plt.figure(figsize=(14.9, 8.38), dpi=250, facecolor="white")
    gs = fig.add_gridspec(
        nrows=2,
        ncols=3,
        width_ratios=[1.0, 1.0, 1.0],
        left=0.075,
        right=0.982,
        top=0.79,
        bottom=0.13,
        hspace=0.255,
        wspace=0.135,
    )
    axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]

    for c, cent in enumerate(CENTRALITIES):
        axes[0][c].set_title(f"{cent} centrality", fontsize=12.6, fontweight="bold", pad=7)
        for r, (slug, _label) in enumerate(MODEL_ROWS):
            sub = df[(df["model_slug"].eq(slug)) & (df["centrality_bin"].eq(cent))].sort_values("bin_lo")
            auc = float(auc_by.loc[(slug, cent), "auc"])
            gain = gains[cent] if slug == "baseline_plus_isolation_inputs" else None
            ylabel = None
            if c == 0:
                ylabel = ("No isolation input\nDensity" if r == 0 else "Isolation input included\nDensity")
            plot_score_panel(axes[r][c], sub, auc, gain, r == 1, ylabel)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.392, 0.862), ncol=2, frameon=False, fontsize=11.3, handlelength=2.4)
    fig.text(0.075, 0.962, "Isolation input gives a centrality-dependent separation gain", ha="left", va="top", fontsize=17.0, fontweight="bold", color="#111827")
    fig.text(
        0.075,
        0.925,
        r"8 $E_T$ x 7 centrality binned BDTs, Photon12+20 signal vs Jet12+20+30 background, 15 < cluster $E_T$ < 35 GeV",
        ha="left",
        va="top",
        fontsize=11.1,
        color="0.25",
    )
    add_sphenix_label(fig)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, facecolor="white")
    print(OUT_PNG)


if __name__ == "__main__":
    main()
