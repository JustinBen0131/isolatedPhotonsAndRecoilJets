#!/usr/bin/env python3
"""Audience-facing split diagnostic slide with compact term definitions."""

from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


OUTDIR = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "ppg12_weighted_centinput_overfitdiag_20260527_1458/slide_pngs"
)
OUTPNG = OUTDIR / "overfitdiag_transfer_story_clean.png"

SPLITS = ["90/10\n90% train", "50/50\n50% train", "10/90\n10% train"]
X = np.arange(len(SPLITS))

TRAIN_AUC = np.array([0.812006, 0.813699, 0.830674])
HOLDOUT_AUC = np.array([0.809383, 0.809211, 0.805470])
FULLSTAT_AUC = np.array([0.821053, 0.820764, 0.820117])

TRAIN_LOGLOSS = np.array([0.506162, 0.504753, 0.491257])
HOLDOUT_LOGLOSS = np.array([0.507955, 0.508777, 0.512486])

AUC_GAP = TRAIN_AUC - HOLDOUT_AUC
LOGLOSS_GAP = HOLDOUT_LOGLOSS - TRAIN_LOGLOSS

PALETTE = {
    "ink": "#161a20",
    "muted": "#5d6675",
    "light": "#f5f7fb",
    "line": "#d7dee9",
    "edge": "#cbd5e1",
    "train": "#d97706",
    "holdout": "#2563ad",
    "full": "#2f855a",
    "gap": "#b91c1c",
}


def card(fig, ax, pad=0.011):
    bbox = ax.get_position()
    fig.patches.append(
        FancyBboxPatch(
            (bbox.x0 - pad, bbox.y0 - pad),
            bbox.width + 2 * pad,
            bbox.height + 2 * pad,
            boxstyle="round,pad=0.008,rounding_size=0.014",
            transform=fig.transFigure,
            facecolor=PALETTE["light"],
            edgecolor=PALETTE["edge"],
            linewidth=1.0,
            zorder=-10,
        )
    )


def card_at(fig, xywh):
    x, y, w, h = xywh
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.010,rounding_size=0.018",
            transform=fig.transFigure,
            facecolor=PALETTE["light"],
            edgecolor=PALETTE["edge"],
            linewidth=1.1,
            zorder=-10,
        )
    )


def plot_axis(fig, xywh, title):
    x, y, w, h = xywh
    card_at(fig, xywh)
    fig.text(
        x + 0.022,
        y + h - 0.025,
        title,
        fontsize=17.2,
        fontweight="bold",
        color=PALETTE["ink"],
        ha="left",
        va="top",
    )
    ax = fig.add_axes([x + 0.035, y + 0.055, w - 0.060, h - 0.125])
    ax.set_facecolor("white")
    return ax


def style_axis(ax):
    ax.grid(axis="y", color=PALETTE["line"], linewidth=1.0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#9aa6b7")
    ax.spines["bottom"].set_color("#9aa6b7")
    ax.tick_params(colors=PALETTE["muted"], labelsize=11.5)
    ax.yaxis.label.set_color(PALETTE["ink"])
    ax.xaxis.label.set_color(PALETTE["ink"])


def draw_series(ax, y, color, marker, label):
    ax.plot(
        X,
        y,
        marker=marker,
        markersize=8.5,
        linewidth=3.0,
        color=color,
        label=label,
    )


def endpoint_label(ax, y, text, color, dy=0.0):
    ax.text(
        X[-1] + 0.06,
        y[-1] + dy,
        text,
        color=color,
        fontsize=11.6,
        fontweight="bold",
        ha="left",
        va="center",
    )


def delta_text(start, end, higher_is_better=True):
    delta = end - start
    sign = "+" if delta >= 0 else ""
    if higher_is_better:
        word = "better" if delta > 0 else "worse"
    else:
        word = "better" if delta < 0 else "worse"
    return f"{sign}{delta:.4f}", word


def wrapped(ax, x, y, text, width, **kwargs):
    ax.text(x, y, textwrap.fill(text, width=width), **kwargs)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "axes.titleweight": "bold",
            "figure.facecolor": "white",
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=160)

    fig.text(
        0.055,
        0.925,
        "Does the 10/90 training gain transfer?",
        fontsize=30,
        fontweight="bold",
        color=PALETTE["ink"],
        ha="left",
        va="center",
    )
    fig.text(
        0.055,
        0.870,
        "Each x-axis point is a separately trained model: 90/10, 50/50, or 10/90 train/holdout split.",
        fontsize=15.0,
        color=PALETTE["muted"],
        ha="left",
        va="center",
    )

    top_y, top_h = 0.485, 0.345
    bot_y, bot_h = 0.065, 0.345
    left = (0.055, top_y, 0.445, top_h)
    right = (0.535, top_y, 0.410, top_h)
    bottom_left = (0.055, bot_y, 0.430, bot_h)
    bottom_right = (0.515, bot_y, 0.430, bot_h)

    ax_auc = plot_axis(fig, left, "AUC: train rises; holdout/full sample do not")
    ax_loss = plot_axis(fig, right, "Logloss: train vs holdout for each model")
    ax_gap = plot_axis(fig, bottom_left, "Train-validation gaps")
    card_at(fig, bottom_right)
    ax_defs = fig.add_axes([bottom_right[0] + 0.035, bottom_right[1] + 0.040, bottom_right[2] - 0.060, bottom_right[3] - 0.070])

    draw_series(ax_auc, TRAIN_AUC, PALETTE["train"], "o", "Train")
    draw_series(ax_auc, HOLDOUT_AUC, PALETTE["holdout"], "o", "Holdout validation")
    draw_series(ax_auc, FULLSTAT_AUC, PALETTE["full"], "o", "Full scored sample")
    ax_auc.set_ylabel("AUC  (higher is better)")
    ax_auc.set_xticks(X, SPLITS)
    ax_auc.set_xlim(-0.10, 2.34)
    ax_auc.set_ylim(0.800, 0.838)
    style_axis(ax_auc)
    ax_auc.legend(frameon=False, loc="upper left", fontsize=11.8, ncol=3, handlelength=2.0)

    draw_series(ax_loss, TRAIN_LOGLOSS, PALETTE["train"], "o", "Train")
    draw_series(ax_loss, HOLDOUT_LOGLOSS, PALETTE["holdout"], "o", "Holdout validation")
    ax_loss.set_ylabel("Logloss  (lower is better)")
    ax_loss.set_xticks(X, SPLITS)
    ax_loss.set_xlim(-0.10, 2.34)
    ax_loss.set_ylim(0.489, 0.514)
    style_axis(ax_loss)
    ax_loss.legend(frameon=False, loc="upper left", fontsize=11.8, ncol=2, handlelength=2.0)

    width = 0.33
    ax_gap.bar(X - width / 2, AUC_GAP, width, color=PALETTE["holdout"], label="AUC gap")
    ax_gap.bar(X + width / 2, LOGLOSS_GAP, width, color=PALETTE["gap"], label="Logloss gap")
    ax_gap.set_ylabel("Gap size")
    ax_gap.set_xticks(X, SPLITS)
    ax_gap.set_ylim(0, 0.029)
    style_axis(ax_gap)
    ax_gap.legend(frameon=False, fontsize=11.8, ncol=2, loc="upper left")
    for i in range(len(X)):
        ax_gap.text(X[i] - width / 2, AUC_GAP[i] + 0.0007, f"{AUC_GAP[i]:.4f}", ha="center", fontsize=10.8, color=PALETTE["holdout"], fontweight="bold")
        ax_gap.text(X[i] + width / 2, LOGLOSS_GAP[i] + 0.0007, f"{LOGLOSS_GAP[i]:.4f}", ha="center", fontsize=10.8, color=PALETTE["gap"], fontweight="bold")
    ax_defs.axis("off")
    ax_defs.set_facecolor("none")
    ax_defs.set_xlim(0, 1)
    ax_defs.set_ylim(0, 1)
    ax_defs.text(0.00, 0.98, "How to read this slide", fontsize=19.2, fontweight="bold", color=PALETTE["ink"], ha="left", va="top")
    ax_defs.text(
        0.00,
        0.855,
        "Each column is one separately trained split model.",
        fontsize=12.8,
        color=PALETTE["muted"],
        ha="left",
        va="top",
    )

    guide_rows = [
        (PALETTE["train"], "Train", "rows used to train that model."),
        (PALETTE["holdout"], "Holdout validation", "reserved split sample; not seen during training."),
        (PALETTE["full"], "Full scored sample", "same 15.85M-row pool; includes train + holdout."),
        (PALETTE["gap"], "Gap bars", "size of train-holdout mismatch."),
    ]
    y = 0.710
    for color, term, desc in guide_rows:
        ax_defs.scatter([0.028], [y - 0.006], s=100, color=color, clip_on=False)
        ax_defs.text(0.078, y, term, fontsize=12.6, fontweight="bold", color=PALETTE["ink"], ha="left", va="top")
        wrapped(ax_defs, 0.480, y, desc, 34, fontsize=11.8, color=PALETTE["muted"], ha="left", va="top")
        y -= 0.146

    ax_defs.plot([0.00, 0.98], [0.145, 0.145], color=PALETTE["line"], linewidth=1.2)
    wrapped(
        ax_defs,
        0.00,
        0.095,
        "Bottom line: 10/90 looks better only on training rows; holdout AUC falls and holdout logloss worsens.",
        68,
        fontsize=13.2,
        color=PALETTE["ink"],
        ha="left",
        va="top",
        fontweight="bold",
    )

    fig.savefig(OUTPNG, dpi=160)
    plt.close(fig)
    print(OUTPNG.resolve())


if __name__ == "__main__":
    main()
