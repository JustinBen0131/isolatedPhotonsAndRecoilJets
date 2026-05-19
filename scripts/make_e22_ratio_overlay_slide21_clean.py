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
VALIDATION_DIR = ROOT / "validation" / "basev3e_w33_e22ratio_20260518_192630" / "validation"
OUTDIR = ROOT / "slideReady" / "basev3e_w33_e22ratio_diagnostics"

HIST_CSV = VALIDATION_DIR / "energy_ratio_histograms_E22E37_E22E53_fullstat.csv"
GAIN_CSV = OUTDIR / "basev3e_w33_E22E37_E22E53_split_gain.csv"
OUT_PNG = OUTDIR / "basev3e_w33_E22E37_E22E53_energy_ratio_overlays_by_cent_slide21_clean.png"

FEATURE_LABELS = {
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "e32_over_e35": r"$E_{32}/E_{35}$",
}

CENTRALITY_ORDER = ["0-20%", "20-50%", "50-80%"]
SIGNAL_COLOR = "#0072B2"
BACKGROUND_COLOR = "#D55E00"


def add_sphenix_label(fig: plt.Figure) -> None:
    fig.text(
        0.835,
        0.958,
        "sPHENIX",
        ha="right",
        va="top",
        fontsize=16.5,
        fontstyle="italic",
        fontweight="bold",
    )
    fig.text(0.945, 0.958, "Internal", ha="right", va="top", fontsize=16.5)


def style_axis(ax: plt.Axes) -> None:
    ax.grid(True, color="0.90", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", direction="in", top=True, right=True, length=4.5, width=1.0, labelsize=10.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.05)
        spine.set_color("0.25")


def draw_panel(ax: plt.Axes, data: pd.DataFrame, feature: str, cent: str) -> None:
    panel = data[data["feature"].eq(feature) & data["cent_bin"].eq(cent)]
    for cls, color, label in [
        ("signal", SIGNAL_COLOR, "signal"),
        ("background", BACKGROUND_COLOR, "background"),
    ]:
        sub = panel[panel["class"].eq(cls)].sort_values("bin_low")
        x = 0.5 * (sub["bin_low"].to_numpy() + sub["bin_high"].to_numpy())
        y = sub["density"].to_numpy()
        ax.fill_between(x, y, step="mid", color=color, alpha=0.115, linewidth=0)
        ax.step(x, y, where="mid", color=color, linewidth=2.15, label=label)
    ymax = max(float(panel["density"].max()) * 1.22, 0.1)
    ax.set_ylim(0.0, ymax)
    ax.set_xlim(0.0, 1.02)
    style_axis(ax)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    hist = pd.read_csv(HIST_CSV)
    gain = pd.read_csv(GAIN_CSV).set_index("feature")
    feature_order = [
        feature
        for feature in gain.sort_values("split_gain_fraction", ascending=False).index
        if feature in FEATURE_LABELS and feature in set(hist["feature"])
    ]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.titleweight": "bold",
        }
    )
    fig, axes = plt.subplots(
        len(feature_order),
        len(CENTRALITY_ORDER),
        figsize=(15.2, 8.55),
        dpi=240,
        sharex=True,
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.103, right=0.982, top=0.805, bottom=0.112, wspace=0.105, hspace=0.205)

    for r, feature in enumerate(feature_order):
        gain_pct = 100.0 * float(gain.loc[feature, "split_gain_fraction"])
        for c, cent in enumerate(CENTRALITY_ORDER):
            ax = axes[r, c]
            draw_panel(ax, hist, feature, cent)
            if r == 0:
                ax.set_title(cent, fontsize=15.0, pad=8)
            if c == 0:
                ax.set_ylabel(f"{FEATURE_LABELS[feature]}\nDensity", fontsize=12.0, labelpad=10)
                ax.text(
                    -0.015,
                    0.88,
                    f"gain {gain_pct:.1f}%",
                    transform=ax.transAxes,
                    ha="left",
                    va="center",
                    fontsize=10.8,
                    fontweight="bold",
                    color="0.20",
                    bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="0.76", linewidth=0.7),
                )
            else:
                ax.set_ylabel("")
            if r == len(feature_order) - 1:
                ax.set_xlabel("Feature value", fontsize=12.2, labelpad=7)
            else:
                ax.tick_params(labelbottom=False)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.52, 0.885),
        ncol=2,
        frameon=False,
        fontsize=12.6,
        handlelength=2.3,
        columnspacing=2.2,
    )
    fig.text(
        0.103,
        0.948,
        "Two-ratio BDT input-shape check",
        ha="left",
        va="top",
        fontsize=17.0,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.103,
        0.913,
        r"Normalized validation shapes, 15 < cluster $E_T$ < 35 GeV; rows ordered by split-gain usage",
        ha="left",
        va="top",
        fontsize=12.4,
        color="0.25",
    )
    add_sphenix_label(fig)
    fig.savefig(OUT_PNG, facecolor="white")
    print(OUT_PNG)


if __name__ == "__main__":
    main()
