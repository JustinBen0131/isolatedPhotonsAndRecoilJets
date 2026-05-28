#!/usr/bin/env python3
"""Build a slide-ready PNG summarizing which split overfits worst."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch
from textwrap import fill


OUTDIR = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "ppg12_weighted_centinput_overfitdiag_20260527_1458/slide_pngs"
)
OUTPNG = OUTDIR / "overfitdiag_which_split_overfits_worst.png"


SPLITS = ["90/10", "50/50", "10/90"]
TRAIN_FRACTION = np.array([0.90, 0.50, 0.10])
TRAIN_AUC = np.array([0.812006, 0.813699, 0.830674])
HOLDOUT_AUC = np.array([0.809383, 0.809211, 0.805470])
TRAIN_LOGLOSS = np.array([0.506162, 0.504753, 0.491257])
HOLDOUT_LOGLOSS = np.array([0.507955, 0.508777, 0.512486])
FULLSTAT_AUC = np.array([0.821053, 0.820764, 0.820117])

AUC_GAP = TRAIN_AUC - HOLDOUT_AUC
LOGLOSS_GAP = HOLDOUT_LOGLOSS - TRAIN_LOGLOSS


COLORS = {
    "ink": "#1b1f23",
    "muted": "#5b6472",
    "grid": "#d8dee9",
    "blue": "#2b6cb0",
    "orange": "#d97706",
    "red": "#c53030",
    "green": "#2f855a",
    "panel": "#f7f8fb",
    "panel_edge": "#c9d1dc",
    "cream": "#fff7ed",
}


def add_panel_background(fig, ax, pad=0.012, face="#f7f8fb", edge="#d6dde8"):
    bbox = ax.get_position()
    rect = FancyBboxPatch(
        (bbox.x0 - pad, bbox.y0 - pad),
        bbox.width + 2 * pad,
        bbox.height + 2 * pad,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        transform=fig.transFigure,
        facecolor=face,
        edgecolor=edge,
        linewidth=1.2,
        zorder=-10,
    )
    fig.patches.append(rect)


def style_axis(ax):
    ax.grid(axis="y", color=COLORS["grid"], linewidth=1.0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#9aa4b2")
    ax.spines["bottom"].set_color("#9aa4b2")
    ax.tick_params(colors=COLORS["muted"], labelsize=12)
    ax.yaxis.label.set_color(COLORS["ink"])
    ax.xaxis.label.set_color(COLORS["ink"])


def bar_labels(ax, bars, values, fmt, color=COLORS["ink"], dy=0.001):
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + dy,
            fmt.format(value),
            ha="center",
            va="bottom",
            fontsize=12,
            color=color,
            fontweight="bold",
        )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    x = np.arange(len(SPLITS))

    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    gs = fig.add_gridspec(
        3,
        6,
        height_ratios=[0.55, 2.1, 1.55],
        width_ratios=[1, 1, 1, 1, 1, 1],
        hspace=0.34,
        wspace=0.36,
    )

    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis("off")
    title_ax.text(
        0.0,
        0.74,
        "Which split overfits worst?",
        fontsize=33,
        fontweight="bold",
        color=COLORS["ink"],
        ha="left",
        va="center",
    )
    title_ax.text(
        0.0,
        0.24,
        "Instrumented BDT split replay: train-vs-holdout diagnostics + full-stat validation",
        fontsize=15,
        color=COLORS["muted"],
        ha="left",
        va="center",
    )
    title_ax.text(
        0.995,
        0.58,
        "Verdict: 10/90 is the overfit case",
        fontsize=17,
        color="white",
        fontweight="bold",
        ha="right",
        va="center",
        bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.16", fc=COLORS["red"], ec=COLORS["red"]),
    )

    ax_auc = fig.add_subplot(gs[1, 0:2])
    ax_loss = fig.add_subplot(gs[1, 2:4])
    ax_ind = fig.add_subplot(gs[1, 4:6])
    ax_gap = fig.add_subplot(gs[2, 0:3])
    ax_msg = fig.add_subplot(gs[2, 3:6])

    for ax in (ax_auc, ax_loss, ax_ind, ax_gap, ax_msg):
        add_panel_background(fig, ax)

    width = 0.34
    b1 = ax_auc.bar(x - width / 2, TRAIN_AUC, width, color=COLORS["orange"], label="Train")
    b2 = ax_auc.bar(x + width / 2, HOLDOUT_AUC, width, color=COLORS["blue"], label="Holdout")
    ax_auc.set_title("Train vs holdout AUC", loc="left", fontsize=17, fontweight="bold", color=COLORS["ink"])
    ax_auc.set_ylabel("AUC")
    ax_auc.set_xticks(x, SPLITS)
    ax_auc.set_ylim(0.798, 0.836)
    ax_auc.legend(frameon=False, fontsize=12, loc="upper left")
    style_axis(ax_auc)
    for i, gap in enumerate(AUC_GAP):
        y = max(TRAIN_AUC[i], HOLDOUT_AUC[i]) + 0.0010
        ax_auc.plot([i - width / 2, i + width / 2], [y, y], color=COLORS["red"], linewidth=1.6)
        ax_auc.text(i, y + 0.0011, f"gap {gap:.4f}", ha="center", fontsize=11, color=COLORS["red"], fontweight="bold")
    bar_labels(ax_auc, b1, TRAIN_AUC, "{:.3f}", color=COLORS["orange"], dy=0.00055)
    bar_labels(ax_auc, b2, HOLDOUT_AUC, "{:.3f}", color=COLORS["blue"], dy=0.00055)

    b3 = ax_loss.bar(x - width / 2, TRAIN_LOGLOSS, width, color=COLORS["orange"], label="Train")
    b4 = ax_loss.bar(x + width / 2, HOLDOUT_LOGLOSS, width, color=COLORS["blue"], label="Holdout")
    ax_loss.set_title("Train vs holdout logloss", loc="left", fontsize=17, fontweight="bold", color=COLORS["ink"])
    ax_loss.set_ylabel("Logloss")
    ax_loss.set_xticks(x, SPLITS)
    ax_loss.set_ylim(0.486, 0.516)
    style_axis(ax_loss)
    for i, gap in enumerate(LOGLOSS_GAP):
        y = max(TRAIN_LOGLOSS[i], HOLDOUT_LOGLOSS[i]) + 0.0007
        ax_loss.plot([i - width / 2, i + width / 2], [y, y], color=COLORS["red"], linewidth=1.6)
        ax_loss.text(i, y + 0.0009, f"gap {gap:.4f}", ha="center", fontsize=11, color=COLORS["red"], fontweight="bold")
    bar_labels(ax_loss, b3, TRAIN_LOGLOSS, "{:.3f}", color=COLORS["orange"], dy=0.00045)
    bar_labels(ax_loss, b4, HOLDOUT_LOGLOSS, "{:.3f}", color=COLORS["blue"], dy=0.00045)

    ax_ind.plot(TRAIN_FRACTION * 100, FULLSTAT_AUC, marker="o", markersize=8, color=COLORS["green"], linewidth=2.8)
    for frac, val, split in zip(TRAIN_FRACTION * 100, FULLSTAT_AUC, SPLITS):
        ax_ind.text(frac, val + 0.00017, f"{split}\n{val:.6f}", ha="center", fontsize=11, color=COLORS["ink"])
    ax_ind.set_title("Full-stat validation AUC", loc="left", fontsize=17, fontweight="bold", color=COLORS["ink"])
    ax_ind.set_xlabel("Training fraction (%)")
    ax_ind.set_ylabel("AUC")
    ax_ind.set_xlim(5, 95)
    ax_ind.set_xticks([10, 50, 90])
    ax_ind.set_ylim(0.8197, 0.82135)
    style_axis(ax_ind)
    ax_ind.text(
        0.04,
        0.08,
        "Independent full-stat score pass:\nless training is not better globally.",
        transform=ax_ind.transAxes,
        fontsize=12,
        color=COLORS["muted"],
        ha="left",
        va="bottom",
    )

    gap_x = np.arange(2)
    gap_width = 0.24
    gap_values = np.vstack([AUC_GAP, LOGLOSS_GAP]).T
    gap_colors = [COLORS["blue"], COLORS["orange"], COLORS["red"]]
    for i, split in enumerate(SPLITS):
        bars = ax_gap.bar(gap_x + (i - 1) * gap_width, gap_values[i], gap_width, label=split, color=gap_colors[i])
        for bar, val in zip(bars, gap_values[i]):
            ax_gap.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.0008,
                f"{val:.4f}",
                ha="center",
                fontsize=11,
                color=gap_colors[i],
                fontweight="bold",
            )
    ax_gap.set_title("Overfit gaps: 10/90 is the outlier", loc="left", fontsize=17, fontweight="bold", color=COLORS["ink"])
    ax_gap.set_xticks(gap_x, ["Train AUC - holdout AUC", "Holdout logloss - train logloss"])
    ax_gap.set_ylim(0, 0.030)
    ax_gap.set_ylabel("Gap")
    ax_gap.legend(frameon=False, fontsize=12, ncol=3, loc="upper left")
    style_axis(ax_gap)
    ax_gap.text(
        0.985,
        0.91,
        f"10/90 = {AUC_GAP[2] / AUC_GAP[0]:.1f}x the 90/10 AUC gap\n"
        f"and {LOGLOSS_GAP[2] / LOGLOSS_GAP[0]:.1f}x the 90/10 logloss gap",
        transform=ax_gap.transAxes,
        ha="right",
        va="top",
        fontsize=12,
        color=COLORS["red"],
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.28,rounding_size=0.08", fc="white", ec="#e5e7eb", alpha=0.94),
    )

    ax_msg.axis("off")
    ax_msg.text(0.03, 0.84, "What proves overfitting here?", fontsize=18, fontweight="bold", color=COLORS["ink"], ha="left")
    bullets = [
        ("Train AUC jumps only for 10/90", "0.812 -> 0.831 while holdout AUC drops to 0.805."),
        ("Train loss keeps improving", "10/90 train logloss is lowest, but holdout logloss is worst."),
        ("Full-stat validation does not reward 10/90", "Global AUC: 90/10 highest, 10/90 lowest."),
    ]
    y = 0.70
    for lead, detail in bullets:
        ax_msg.text(0.06, y, lead, fontsize=14.2, fontweight="bold", color=COLORS["ink"], ha="left", va="center")
        ax_msg.text(0.06, y - 0.088, detail, fontsize=12.6, color=COLORS["muted"], ha="left", va="center")
        y -= 0.205
    fig.subplots_adjust(left=0.055, right=0.975, top=0.94, bottom=0.07)
    fig.savefig(OUTPNG, dpi=160)
    plt.close(fig)
    print(OUTPNG.resolve())


if __name__ == "__main__":
    main()
