#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
RUN_ROOT = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439"
)
HIST_CSV = (
    RUN_ROOT
    / "validation/basev3e_w33_e22ratio_20260518_192630/validation/"
    / "energy_ratio_histograms_E22E37_E22E53_fullstat.csv"
)
GAIN_CSV = (
    RUN_ROOT
    / "slideReady/basev3e_w33_e22ratio_diagnostics/basev3e_w33_E22E37_E22E53_split_gain.csv"
)
OUTDIR = REPO / "dataOutput/slide_assets/WP_GammaJets_5_20_26/slide19_e22_shape_separation"
OUT_PNG = OUTDIR / "slide19_e22_shape_separation_candidate.png"
ROUTED_GAIN_CSV = OUTDIR / "ptcent3_routed_split_gain_by_centrality.csv"


FEATURE_LABELS = {
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "e32_over_e35": r"$E_{32}/E_{35}$",
}


def add_box(
    fig: plt.Figure,
    xy: tuple[float, float],
    wh: tuple[float, float],
    face: str,
    *,
    edge: str = "none",
    lw: float = 0.0,
    radius: float = 0.018,
    zorder: int = -1,
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            xy,
            wh[0],
            wh[1],
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def add_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    weight: str = "normal",
    color: str = "#111827",
    ha: str = "left",
    va: str = "top",
    linespacing: float = 1.08,
) -> None:
    fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
    )


def ordered_features(hist: pd.DataFrame) -> list[str]:
    rows = (
        hist[["feature", "feature_order_by_gain", "split_gain"]]
        .drop_duplicates()
        .sort_values(["feature_order_by_gain", "split_gain"], ascending=[True, False])
    )
    return [f for f in rows["feature"].tolist() if f in FEATURE_LABELS]


def draw() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    hist = pd.read_csv(HIST_CSV)
    global_gain_lookup = pd.read_csv(GAIN_CSV).set_index("feature")["split_gain_fraction"].to_dict()
    routed_gain = pd.read_csv(ROUTED_GAIN_CSV)
    routed_gain_lookup = {
        (row["centrality"], row["feature"]): float(row["split_gain_percent"])
        for row in routed_gain.to_dict("records")
    }
    features = ordered_features(hist)
    cent_order = ["0-20%", "20-50%", "50-80%"]

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.unicode_minus": False,
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=180, facecolor="white")
    fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, facecolor="white", zorder=-5))

    navy = "#101828"
    gray = "#475467"
    soft_yellow = "#FFF2B8"
    signal_color = "#1F77B4"
    bkg_color = "#D95F02"
    gain_box = "#F2F4F7"

    add_text(
        fig,
        0.035,
        0.958,
        "E22 ratios give the clearest signal/background shape separation",
        size=27.5,
        weight="bold",
        color=navy,
    )
    add_text(
        fig,
        0.035,
        0.913,
        "Feature-shape validation for the two-ratio ablation, 15 < cluster $E_T$ < 35 GeV; rows follow global split-gain ordering.",
        size=14.2,
        color=gray,
    )

    add_box(fig, (0.045, 0.148), (0.910, 0.715), "white", edge="#D0D5DD", lw=1.0, radius=0.012)
    add_text(
        fig,
        0.500,
        0.858,
        "Centrality-independent split gain from the global two-ratio BDT",
        size=12.7,
        weight="bold",
        color=navy,
        ha="center",
    )

    gain_box_y = 0.789
    gain_box_h = 0.032
    gain_box_w = 0.145
    gain_box_gap = 0.018
    gain_box_left = 0.5 - (4 * gain_box_w + 3 * gain_box_gap) / 2
    for i, feature in enumerate(features):
        x = gain_box_left + i * (gain_box_w + gain_box_gap)
        add_box(fig, (x, gain_box_y), (gain_box_w, gain_box_h), gain_box, edge="#D0D5DD", lw=0.7, radius=0.008)
        add_text(
            fig,
            x + gain_box_w / 2,
            gain_box_y + 0.024,
            f"{FEATURE_LABELS[feature]}  {100.0 * float(global_gain_lookup.get(feature, 0.0)):.1f}%",
            size=10.4,
            weight="bold",
            color=navy,
            ha="center",
        )

    add_text(
        fig,
        0.500,
        0.725,
        r"Panel labels show the 32-input $8\,E_T \times 3$ routed-BDT gain, summed over fine $E_T$ routes inside each centrality column.",
        size=9.7,
        weight="bold",
        color=gray,
        ha="center",
    )

    axes = []
    left, right = 0.090, 0.940
    bottom, top = 0.225, 0.675
    wspace, hspace = 0.040, 0.038
    ncols, nrows = 3, len(features)
    ax_w = (right - left - wspace * (ncols - 1)) / ncols
    ax_h = (top - bottom - hspace * (nrows - 1)) / nrows

    for r, feature in enumerate(features):
        row_axes = []
        row = hist[hist["feature"].eq(feature)]
        ymax = max(float(row["density"].max()) * 1.18, 0.1)
        for c, cent in enumerate(cent_order):
            x0 = left + c * (ax_w + wspace)
            y0 = top - (r + 1) * ax_h - r * hspace
            ax = fig.add_axes([x0, y0, ax_w, ax_h])
            ax.set_facecolor("white")
            for cls, color, label in [
                ("signal", signal_color, "Signal"),
                ("background", bkg_color, "Background"),
            ]:
                sub = row[row["cent_bin"].eq(cent) & row["class"].eq(cls)].sort_values("bin_low")
                x = 0.5 * (sub["bin_low"].to_numpy() + sub["bin_high"].to_numpy())
                y = sub["density"].to_numpy()
                ax.step(x, y, where="mid", color=color, linewidth=2.0, label=label)
                ax.fill_between(x, y, step="mid", color=color, alpha=0.13)
            gain_percent = routed_gain_lookup.get((cent, feature))
            if gain_percent is not None:
                ax.text(
                    0.045,
                    0.840,
                    f"routed gain {gain_percent:.1f}%",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=9.2,
                    fontweight="bold",
                    color=navy,
                    bbox={
                        "boxstyle": "round,pad=0.20,rounding_size=0.02",
                        "facecolor": "white",
                        "edgecolor": "#D0D5DD",
                        "linewidth": 0.6,
                        "alpha": 0.88,
                    },
                )
            ax.set_xlim(0.0, 1.2)
            ax.set_ylim(0.0, ymax)
            ax.grid(color="#E4E7EC", linewidth=0.7)
            ax.set_axisbelow(True)
            ax.tick_params(axis="both", labelsize=8.8, direction="in", top=True, right=True, length=3.5)
            for spine in ax.spines.values():
                spine.set_linewidth(0.85)
                spine.set_color("#344054")
            if r == 0:
                ax.set_title(cent, fontsize=14.5, fontweight="bold", color=navy, pad=7)
            if c == 0:
                ax.set_ylabel(f"{FEATURE_LABELS[feature]}\nDensity", fontsize=12.8, color=navy, labelpad=8)
            else:
                ax.set_ylabel("")
            if r == nrows - 1:
                ax.set_xlabel("Feature value", fontsize=11.8, color=navy, labelpad=5)
            else:
                ax.set_xticklabels([])
            row_axes.append(ax)
        axes.append(row_axes)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.500, 0.777),
        ncol=2,
        frameon=False,
        fontsize=17.5,
        handlelength=3.0,
        columnspacing=4.0,
    )

    add_box(fig, (0.055, 0.045), (0.890, 0.075), soft_yellow, radius=0.018)
    add_text(fig, 0.075, 0.099, "Takeaway", size=16.2, weight="bold", color=navy)
    add_text(
        fig,
        0.075,
        0.073,
        r"$E_{22}/E_{37}$ and $E_{22}/E_{53}$ sharply separate photon-like showers from jet background across centrality; "
        r"the older $E_{11}/E_{33}$ and $E_{32}/E_{35}$ ratios are still useful but less dominant in the trained BDT.",
        size=11.9,
        color="#1F2937",
    )

    fig.savefig(OUT_PNG, dpi=180, facecolor="white")
    plt.close(fig)
    print(OUT_PNG)


if __name__ == "__main__":
    draw()
