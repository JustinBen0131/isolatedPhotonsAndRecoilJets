#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.patches import FancyBboxPatch, Rectangle


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
SOURCE_DIR = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439"
    / "slideReady/basev3e_w33_feature_diagnostics"
)
CSV_PATH = SOURCE_DIR / "basev3e_w33_split_gain.csv"
OUTDIR = REPO / "dataOutput/slide_assets/WP_GammaJets_5_20_26/slide17_split_gain"
OUT_PATH = OUTDIR / "slide17_split_gain_explainer_candidate.png"


LABELS = {
    "e32_over_e35": r"$E_{32}/E_{35}$",
    "cluster_et1": "cluster et1",
    "cluster_wphi33_cogx": r"$w_{\phi}^{3x3}$",
    "cluster_weta_cogx": r"$w_{\eta}$",
    "cluster_weta33_cogx": r"$w_{\eta}^{3x3}$",
    "cluster_wphi_cogx": r"$w_{\phi}$",
    "centrality": "centrality",
    "cluster_et3": "cluster et3",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "cluster_et2": "cluster et2",
    "cluster_Et": r"cluster $E_T$",
    "cluster_Eta": r"cluster $\eta$",
    "cluster_et4": "cluster et4",
    "vertexz": r"$z_{vtx}$",
}


def add_box(
    fig: plt.Figure,
    xy: tuple[float, float],
    wh: tuple[float, float],
    face: str,
    edge: str = "none",
    radius: float = 0.02,
    alpha: float = 1.0,
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
            linewidth=0.0 if edge == "none" else 1.0,
            alpha=alpha,
            zorder=-1,
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
    style: str = "normal",
    linespacing: float = 1.12,
) -> None:
    fig.text(
        x,
        y,
        text,
        ha=ha,
        va=va,
        fontsize=size,
        fontweight=weight,
        color=color,
        fontstyle=style,
        linespacing=linespacing,
    )


def draw() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CSV_PATH).copy()
    df["percent"] = 100.0 * df["split_gain_fraction"]
    df = df.sort_values("percent", ascending=True)

    top = df.sort_values("percent", ascending=False).iloc[0]
    top_label = LABELS.get(str(top["feature"]), str(top["feature"]))
    top_pct = float(top["percent"])

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
    soft_gray = "#F2F4F7"
    soft_blue = "#EAF2FB"
    soft_yellow = "#FFF2B8"
    orange = "#F26A21"

    add_text(
        fig,
        0.035,
        0.955,
        "Split gain shows which inputs the trained BDT actually uses",
        size=30.5,
        weight="bold",
        color=navy,
    )
    add_text(
        fig,
        0.035,
        0.905,
        "Base v3E + centrality + 3x3 widths, trained on Photon12+20 vs Jet12+20+30, 15 < cluster $E_T$ < 35 GeV",
        size=15.2,
        color=gray,
    )

    add_box(fig, (0.035, 0.230), (0.355, 0.630), soft_gray, radius=0.025)
    add_text(fig, 0.055, 0.839, "What the metric is measuring", size=21.0, weight="bold", color=navy)

    add_box(fig, (0.055, 0.653), (0.315, 0.132), "white", radius=0.018)
    add_text(fig, 0.072, 0.762, "XGBoost BDT", size=17.6, weight="bold", color=navy)
    xgb_text = (
        "XGBoost trains the BDT by adding\n"
        "many small decision trees. Each tree\n"
        "is a chain of yes/no cuts."
    )
    add_text(fig, 0.072, 0.727, xgb_text, size=13.5, color="#1F2937", linespacing=1.03)

    add_box(fig, (0.055, 0.507), (0.315, 0.122), "white", radius=0.018)
    add_text(fig, 0.072, 0.604, "Tree split and loss", size=17.4, weight="bold", color=navy)
    split_text = (
        "A split is one cut. Training loss is the\n"
        "penalty for wrong labels. Many splits can\n"
        "still have low gain if each cut is weak."
    )
    add_text(fig, 0.072, 0.570, split_text, size=12.9, color="#1F2937", linespacing=1.02)

    add_box(fig, (0.055, 0.361), (0.315, 0.122), "white", radius=0.018)
    add_text(fig, 0.072, 0.458, "Gain", size=17.4, weight="bold", color=navy)
    gain_text = (
        "Gain is how much a split reduces loss.\n"
        "Larger gain means that cut improved\n"
        "the training classification more."
    )
    add_text(fig, 0.072, 0.424, gain_text, size=13.2, color="#1F2937", linespacing=1.03)

    add_box(fig, (0.055, 0.254), (0.315, 0.083), "white", radius=0.018)
    add_text(fig, 0.072, 0.320, "How I calculate it", size=16.2, weight="bold", color=navy)
    calc_text = (
        "From the trained JSON, sum gain by\n"
        "feature; divide by total gain."
    )
    add_text(fig, 0.072, 0.289, calc_text, size=12.7, color="#1F2937", linespacing=1.00)

    add_box(fig, (0.415, 0.230), (0.555, 0.630), "white", edge="#D0D5DD", radius=0.012)
    ax = fig.add_axes([0.475, 0.345, 0.455, 0.400])

    y = np.arange(len(df))
    perc = df["percent"].to_numpy()
    colors = plt.cm.Blues(np.linspace(0.38, 0.88, len(df)))
    feature_order = df["feature"].tolist()
    added_features = {"centrality", "cluster_weta33_cogx", "cluster_wphi33_cogx"}
    for i, feat in enumerate(feature_order):
        if feat in added_features:
            colors[i] = to_rgba(orange)

    bars = ax.barh(y, perc, color=colors, edgecolor="none", height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS.get(v, v) for v in feature_order], fontsize=12.5)
    ax.set_xlabel("Fraction of total split gain [%]", fontsize=14.2, labelpad=7)
    ax.set_xlim(0.0, 34.0)
    ax.set_xticks(np.arange(0, 35, 5))
    ax.tick_params(axis="x", labelsize=12.0, direction="in", top=True, length=5)
    ax.tick_params(axis="y", labelsize=12.2, direction="in", right=True, length=4)
    ax.grid(axis="x", color="#E4E7EC", linewidth=0.9)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
        spine.set_color("#344054")

    for bar, (_, row) in zip(bars, df.iterrows()):
        value = float(row["percent"])
        count = int(row["split_count"])
        label = f"{value:.1f}%  ({count} splits)"
        x = value + 0.42
        ha = "left"
        color = "#111827"
        if value > 26.0:
            label = f"{value:.1f}%\n({count} splits)"
            x = value + 0.34
        ax.text(
            x,
            bar.get_y() + bar.get_height() / 2.0,
            label,
            va="center",
            ha=ha,
            fontsize=10.8,
            fontweight="bold",
            color=color,
            linespacing=0.88,
        )

    add_text(
        fig,
        0.435,
        0.832,
        "Feature usage in the trained global BDT",
        size=20.5,
        weight="bold",
        color=navy,
    )
    add_text(
        fig,
        0.435,
        0.800,
        "Bars sum XGBoost gain over all tree splits, then normalize by total gain across inputs.",
        size=13.8,
        color=gray,
    )
    fig.patches.append(
        Rectangle((0.695, 0.371), 0.013, 0.013, transform=fig.transFigure, facecolor=orange, edgecolor="none", zorder=3)
    )
    add_text(
        fig,
        0.713,
        0.385,
        "added beyond base v3E/PPG12 inputs",
        size=11.6,
        color=gray,
        va="center",
    )

    add_box(fig, (0.035, 0.065), (0.445, 0.126), soft_yellow, radius=0.020)
    add_text(fig, 0.055, 0.161, "Interpretation of Percentages", size=20.2, weight="bold", color=navy)
    add_text(
        fig,
        0.055,
        0.120,
        f"{top_label} at {top_pct:.1f}% means {top_pct:.1f}% of split-level loss reduction\ncame from cuts on that input. This is model usage, not exact AUC loss.",
        size=13.5,
        color="#1F2937",
        linespacing=1.06,
    )

    add_box(fig, (0.510, 0.065), (0.460, 0.126), soft_blue, radius=0.020)
    add_text(fig, 0.530, 0.161, "Physics interpretation", size=20.2, weight="bold", color=navy)
    add_text(
        fig,
        0.530,
        0.120,
        "After reweighting, ET gain means the BDT may use ET as context\n"
        "for which shower-shape cuts are optimal, not as a raw spectrum shortcut.",
        size=13.2,
        color="#1F2937",
        linespacing=1.06,
    )

    fig.savefig(OUT_PATH, dpi=180, facecolor="white")
    plt.close(fig)
    print(OUT_PATH)


if __name__ == "__main__":
    draw()
