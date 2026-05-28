#!/usr/bin/env python3
"""Build a cleaner slide showing whether split-study gains transfer."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


OUTDIR = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "ppg12_weighted_centinput_overfitdiag_20260527_1458/slide_pngs"
)
OUTPNG = OUTDIR / "overfitdiag_transfer_story.png"

SPLITS = ["90/10", "50/50", "10/90"]
X = np.arange(len(SPLITS))
TRAIN_AUC = np.array([0.812006, 0.813699, 0.830674])
HOLDOUT_AUC = np.array([0.809383, 0.809211, 0.805470])
FULLSTAT_AUC = np.array([0.821053, 0.820764, 0.820117])
TRAIN_LOGLOSS = np.array([0.506162, 0.504753, 0.491257])
HOLDOUT_LOGLOSS = np.array([0.507955, 0.508777, 0.512486])
AUC_GAP = TRAIN_AUC - HOLDOUT_AUC
LOGLOSS_GAP = HOLDOUT_LOGLOSS - TRAIN_LOGLOSS

COLORS = {
    "ink": "#171b21",
    "muted": "#5e6878",
    "grid": "#d9e0ea",
    "blue": "#2b6cb0",
    "orange": "#d97706",
    "green": "#2f855a",
    "red": "#c53030",
    "panel": "#f7f8fb",
    "edge": "#d3dae6",
}


def add_card(fig, ax, pad=0.012):
    bbox = ax.get_position()
    fig.patches.append(
        FancyBboxPatch(
            (bbox.x0 - pad, bbox.y0 - pad),
            bbox.width + 2 * pad,
            bbox.height + 2 * pad,
            boxstyle="round,pad=0.008,rounding_size=0.012",
            transform=fig.transFigure,
            facecolor=COLORS["panel"],
            edgecolor=COLORS["edge"],
            linewidth=1.1,
            zorder=-10,
        )
    )


def style_ax(ax):
    ax.grid(axis="y", color=COLORS["grid"], linewidth=1.0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#a7b0be")
    ax.spines["bottom"].set_color("#a7b0be")
    ax.tick_params(labelsize=12, colors=COLORS["muted"])
    ax.yaxis.label.set_color(COLORS["ink"])


def line_label(ax, x, y, label, color, dy=0.0, dx=0.05, size=11):
    ax.text(
        x + dx,
        y + dy,
        label,
        color=color,
        fontsize=size,
        fontweight="bold",
        ha="left",
        va="center",
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    gs = fig.add_gridspec(
        3,
        6,
        height_ratios=[0.58, 2.1, 1.62],
        width_ratios=[1, 1, 1, 1, 1, 1],
        hspace=0.33,
        wspace=0.36,
    )

    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis("off")
    title_ax.text(
        0.0,
        0.72,
        "Does the 10/90 training gain transfer?",
        fontsize=32,
        fontweight="bold",
        color=COLORS["ink"],
        ha="left",
    )
    title_ax.text(
        0.0,
        0.24,
        "The decisive pattern is training-only improvement while holdout/full-stat behavior worsens or stays lower.",
        fontsize=15.5,
        color=COLORS["muted"],
        ha="left",
    )

    ax_auc = fig.add_subplot(gs[1, 0:3])
    ax_loss = fig.add_subplot(gs[1, 3:6])
    ax_gap = fig.add_subplot(gs[2, 0:3])
    ax_delta = fig.add_subplot(gs[2, 3:6])
    for ax in (ax_auc, ax_loss, ax_gap, ax_delta):
        add_card(fig, ax)

    ax_auc.plot(X, TRAIN_AUC, "-o", color=COLORS["orange"], linewidth=3.0, markersize=8)
    ax_auc.plot(X, HOLDOUT_AUC, "-o", color=COLORS["blue"], linewidth=3.0, markersize=8)
    ax_auc.plot(X, FULLSTAT_AUC, "-o", color=COLORS["green"], linewidth=3.0, markersize=8)
    ax_auc.set_title("AUC: training improves, transfer samples do not", loc="left", fontsize=17, fontweight="bold", color=COLORS["ink"])
    ax_auc.set_ylabel("AUC")
    ax_auc.set_xticks(X, SPLITS)
    ax_auc.set_ylim(0.802, 0.8335)
    style_ax(ax_auc)
    line_label(ax_auc, 2, TRAIN_AUC[-1], "train 0.831", COLORS["orange"], dy=0.0008, dx=-0.42)
    line_label(ax_auc, 2, HOLDOUT_AUC[-1], "holdout 0.805", COLORS["blue"], dy=-0.0012, dx=-0.42)
    line_label(ax_auc, 2, FULLSTAT_AUC[-1], "full-stat 0.820", COLORS["green"], dy=0.0012, dx=-0.45)
    ax_auc.text(0.02, 0.05, "less training data  ->", transform=ax_auc.transAxes, fontsize=12.5, color=COLORS["muted"])

    ax_loss.plot(X, TRAIN_LOGLOSS, "-o", color=COLORS["orange"], linewidth=3.0, markersize=8)
    ax_loss.plot(X, HOLDOUT_LOGLOSS, "-o", color=COLORS["blue"], linewidth=3.0, markersize=8)
    ax_loss.set_title("Logloss: train loss improves, holdout loss worsens", loc="left", fontsize=17, fontweight="bold", color=COLORS["ink"])
    ax_loss.set_ylabel("Logloss")
    ax_loss.set_xticks(X, SPLITS)
    ax_loss.set_ylim(0.489, 0.514)
    style_ax(ax_loss)
    line_label(ax_loss, 2, TRAIN_LOGLOSS[-1], "train 0.491", COLORS["orange"], dy=-0.0009, dx=-0.46)
    line_label(ax_loss, 2, HOLDOUT_LOGLOSS[-1], "holdout 0.512", COLORS["blue"], dy=0.0008, dx=-0.49)
    ax_loss.text(0.02, 0.05, "lower is better", transform=ax_loss.transAxes, fontsize=12.5, color=COLORS["muted"])

    width = 0.35
    ax_gap.bar(X - width / 2, AUC_GAP, width, color=COLORS["blue"], label="AUC gap")
    ax_gap.bar(X + width / 2, LOGLOSS_GAP, width, color=COLORS["red"], label="Logloss gap")
    ax_gap.set_title("Transfer gaps grow sharply at 10/90", loc="left", fontsize=17, fontweight="bold", color=COLORS["ink"])
    ax_gap.set_ylabel("Train-transfer gap")
    ax_gap.set_xticks(X, SPLITS)
    ax_gap.set_ylim(0, 0.029)
    ax_gap.legend(frameon=False, fontsize=12, ncol=2, loc="upper left")
    style_ax(ax_gap)
    for i in range(3):
        ax_gap.text(i - width / 2, AUC_GAP[i] + 0.0008, f"{AUC_GAP[i]:.4f}", ha="center", fontsize=11, color=COLORS["blue"], fontweight="bold")
        ax_gap.text(i + width / 2, LOGLOSS_GAP[i] + 0.0008, f"{LOGLOSS_GAP[i]:.4f}", ha="center", fontsize=11, color=COLORS["red"], fontweight="bold")
    ax_delta.axis("off")
    ax_delta.text(0.035, 0.86, "90/10 -> 10/90 numerical deltas", fontsize=17, fontweight="bold", color=COLORS["ink"])
    rows = [
        ("Train AUC", TRAIN_AUC[-1] - TRAIN_AUC[0], "+0.0187", COLORS["orange"], "improves"),
        ("Holdout AUC", HOLDOUT_AUC[-1] - HOLDOUT_AUC[0], "-0.0039", COLORS["blue"], "worsens"),
        ("Full-stat AUC", FULLSTAT_AUC[-1] - FULLSTAT_AUC[0], "-0.0009", COLORS["green"], "worsens"),
        ("Train logloss", TRAIN_LOGLOSS[-1] - TRAIN_LOGLOSS[0], "-0.0149", COLORS["orange"], "improves"),
        ("Holdout logloss", HOLDOUT_LOGLOSS[-1] - HOLDOUT_LOGLOSS[0], "+0.0045", COLORS["blue"], "worsens"),
    ]
    y = 0.70
    for name, _delta, value, color, verdict in rows:
        ax_delta.text(0.06, y, name, fontsize=13.3, color=COLORS["ink"], fontweight="bold", ha="left", va="center")
        ax_delta.text(0.48, y, value, fontsize=13.3, color=color, fontweight="bold", ha="left", va="center")
        ax_delta.text(0.68, y, verdict, fontsize=13.0, color=color, ha="left", va="center")
        y -= 0.125
    fig.subplots_adjust(left=0.055, right=0.975, top=0.94, bottom=0.07)
    fig.savefig(OUTPNG, dpi=160)
    plt.close(fig)
    print(OUTPNG.resolve())


if __name__ == "__main__":
    main()
