#!/usr/bin/env python3
"""One-plot verdict slide for the overfitdiag split study."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


OUTDIR = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "ppg12_weighted_centinput_overfitdiag_20260527_1458/slide_pngs"
)
OUTPNG = OUTDIR / "overfitdiag_generalization_frontier.png"

SPLITS = np.array(["90/10", "50/50", "10/90"])
TRAIN_FRAC = np.array([0.90, 0.50, 0.10])
TRAIN_AUC = np.array([0.812006, 0.813699, 0.830674])
HOLDOUT_AUC = np.array([0.809383, 0.809211, 0.805470])
FULL_AUC = np.array([0.821053, 0.820764, 0.820117])
TRAIN_LOGLOSS = np.array([0.506162, 0.504753, 0.491257])
HOLDOUT_LOGLOSS = np.array([0.507955, 0.508777, 0.512486])

AUC_GAP = TRAIN_AUC - HOLDOUT_AUC
LOGLOSS_GAP = HOLDOUT_LOGLOSS - TRAIN_LOGLOSS

COLORS = {
    "ink": "#151a22",
    "muted": "#475569",
    "panel": "#f7f9fc",
    "edge": "#cbd6e4",
    "grid": "#e0e7f1",
    "good": "#2f855a",
    "mid": "#2563ad",
    "bad": "#c2410c",
    "gold": "#d97706",
}


def add_card(fig, xywh, face=COLORS["panel"], edge=COLORS["edge"]):
    x, y, w, h = xywh
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=1.2,
            zorder=-10,
        )
    )


def style_axis(ax):
    ax.set_axisbelow(True)
    ax.grid(color=COLORS["grid"], linewidth=1.0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#9aa7b8")
    ax.tick_params(labelsize=14.0, colors=COLORS["muted"])
    ax.xaxis.label.set_color(COLORS["ink"])
    ax.yaxis.label.set_color(COLORS["ink"])


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "white",
            "axes.titleweight": "bold",
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=160)

    fig.text(
        0.055,
        0.925,
        "90/10 is the clean split choice",
        fontsize=34,
        fontweight="bold",
        color=COLORS["ink"],
        ha="left",
        va="center",
    )
    fig.text(
        0.055,
        0.870,
        "A high training AUC alone is not success; the check is whether that gain survives on holdout rows the model never trained on.",
        fontsize=17.0,
        color=COLORS["muted"],
        ha="left",
        va="center",
    )

    plot_card = (0.055, 0.105, 0.650, 0.705)
    verdict_card = (0.735, 0.105, 0.210, 0.705)
    add_card(fig, plot_card)
    add_card(fig, verdict_card)

    ax = fig.add_axes([0.104, 0.180, 0.548, 0.550])
    ax.set_facecolor("white")

    ax.axvspan(0.0, 0.0065, color="#dcfce7", alpha=0.28, zorder=0)
    ax.axvspan(0.018, 0.030, color="#fee2e2", alpha=0.26, zorder=0)

    point_colors = [COLORS["good"], COLORS["mid"], COLORS["bad"]]
    marker_size = 420

    ax.plot(
        AUC_GAP,
        HOLDOUT_AUC,
        color="#7c8798",
        linewidth=2.2,
        linestyle=(0, (4, 3)),
        zorder=2,
    )
    ax.annotate(
        "less training data\n90/10 -> 50/50 -> 10/90",
        xy=(0.0189, 0.8061),
        xytext=(0.0106, 0.80705),
        arrowprops=dict(arrowstyle="-|>", color="#64748b", lw=2.0),
        color="#334155",
        fontsize=15.2,
        ha="center",
        va="center",
        fontweight="bold",
    )

    for i, split in enumerate(SPLITS):
        ax.scatter(
            AUC_GAP[i],
            HOLDOUT_AUC[i],
            s=marker_size,
            color=point_colors[i],
            edgecolor="white",
            linewidth=2.4,
            zorder=5,
        )
        if split == "90/10":
            label = "90/10"
            xytext = (AUC_GAP[i] - 0.00010, HOLDOUT_AUC[i] + 0.00055)
            ha = "left"
        elif split == "50/50":
            label = "50/50"
            xytext = (AUC_GAP[i] + 0.00035, HOLDOUT_AUC[i] - 0.00055)
            ha = "left"
        else:
            label = "10/90"
            xytext = (AUC_GAP[i] - 0.00120, HOLDOUT_AUC[i] + 0.00045)
            ha = "right"
        ax.text(
            xytext[0],
            xytext[1],
            label,
            ha=ha,
            va="center",
            fontsize=15.0,
            color=COLORS["ink"],
            fontweight="bold",
            zorder=6,
        )

    ax.set_title(
        "Upper-left is best: high holdout AUC, low gap",
        loc="left",
        fontsize=20.5,
        color=COLORS["ink"],
        pad=10,
    )
    ax.set_xlabel("Train - holdout AUC gap  (lower is better)", fontsize=16.0, labelpad=12)
    ax.set_ylabel("Holdout validation AUC  (higher is better)", fontsize=16.0, labelpad=12)
    ax.set_xlim(0.000, 0.029)
    ax.set_ylim(0.8038, 0.8103)
    style_axis(ax)

    ax_note = fig.add_axes([0.750, 0.145, 0.175, 0.635])
    ax_note.axis("off")
    ax_note.set_xlim(0, 1)
    ax_note.set_ylim(0, 1)

    ax_note.text(0.0, 0.96, "Verdict", fontsize=26, fontweight="bold", color=COLORS["ink"], ha="left", va="top")
    ax_note.text(
        0.0,
        0.835,
        "Choose 90/10.",
        fontsize=23.0,
        fontweight="bold",
        color=COLORS["good"],
        ha="left",
        va="top",
    )

    rows = [
        ("Best holdout AUC", "0.8094"),
        ("Smallest AUC gap", "0.0026"),
        ("Smallest logloss gap", "0.0018"),
    ]
    y = 0.660
    for metric, value in rows:
        ax_note.scatter([0.020], [y - 0.006], s=58, color=COLORS["good"])
        ax_note.text(0.090, y, metric, fontsize=15.0, fontweight="bold", color=COLORS["ink"], ha="left", va="top")
        ax_note.text(1.00, y, value, fontsize=15.0, color=COLORS["muted"], ha="right", va="top")
        y -= 0.135

    ax_note.plot([0.0, 1.0], [0.250, 0.250], color=COLORS["grid"], lw=1.4)
    ax_note.text(0.0, 0.205, "Overfitting call", fontsize=17.2, fontweight="bold", color=COLORS["ink"], ha="left", va="top")
    ax_note.text(
        0.0,
        0.135,
        "10/90 is the overfit case:\ntraining improves while holdout degrades.\n90/10 does not show that gap growth.",
        fontsize=14.0,
        color=COLORS["muted"],
        ha="left",
        va="top",
    )

    fig.savefig(OUTPNG, dpi=160)
    plt.close(fig)
    print(OUTPNG.resolve())


if __name__ == "__main__":
    main()
