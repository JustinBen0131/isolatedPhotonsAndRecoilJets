#!/usr/bin/env python3
"""Build a cleaner slide showing whether split-study gains transfer."""

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


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--internal-metrics-json", type=Path, default=None)
    ap.add_argument("--compact-json", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=OUTPNG)
    ap.add_argument("--title", default="Does the 10/90 training gain transfer?")
    ap.add_argument(
        "--subtitle",
        default="The decisive pattern is training-only improvement while holdout/full-stat behavior worsens or stays lower.",
    )
    return ap.parse_args()


def load_metrics(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    train_auc = TRAIN_AUC.copy()
    holdout_auc = HOLDOUT_AUC.copy()
    fullstat_auc = FULLSTAT_AUC.copy()
    train_logloss = TRAIN_LOGLOSS.copy()
    holdout_logloss = HOLDOUT_LOGLOSS.copy()

    if args.internal_metrics_json is not None:
        payload = json.loads(args.internal_metrics_json.read_text())
        train_auc = np.asarray([payload[s]["train_auc"] for s in SPLITS], dtype=float)
        holdout_auc = np.asarray([payload[s]["auc"] for s in SPLITS], dtype=float)
        train_logloss = np.asarray([payload[s]["train_logloss"] for s in SPLITS], dtype=float)
        holdout_logloss = np.asarray([payload[s]["holdout_logloss"] for s in SPLITS], dtype=float)

    if args.compact_json is not None:
        compact = json.loads(args.compact_json.read_text())
        rows = {row["label"]: row for row in compact["rows"]}
        fullstat_auc = np.asarray([rows[s]["global_auc"] for s in SPLITS], dtype=float)

    return train_auc, holdout_auc, fullstat_auc, train_logloss, holdout_logloss


def padded_limits(values: np.ndarray, frac: float = 0.18, min_pad: float = 0.001) -> tuple[float, float]:
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    pad = max(min_pad, (hi - lo) * frac)
    return lo - pad, hi + pad


def padded_zero_limits(values: np.ndarray, frac: float = 0.22, min_pad: float = 0.0005) -> tuple[float, float]:
    lo = min(0.0, float(np.nanmin(values)))
    hi = max(0.0, float(np.nanmax(values)))
    pad = max(min_pad, (hi - lo) * frac)
    return lo - pad, hi + pad


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


def add_fig_card(fig, xywh, pad=0.0):
    x, y, w, h = xywh
    fig.patches.append(
        FancyBboxPatch(
            (x - pad, y - pad),
            w + 2 * pad,
            h + 2 * pad,
            boxstyle="round,pad=0.010,rounding_size=0.018",
            transform=fig.transFigure,
            facecolor=COLORS["panel"],
            edgecolor=COLORS["edge"],
            linewidth=1.2,
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
    args = parse_args()
    out = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    train_auc, holdout_auc, fullstat_auc, train_logloss, holdout_logloss = load_metrics(args)
    auc_gap = train_auc - holdout_auc
    logloss_gap = holdout_logloss - train_logloss

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "white",
            "axes.titleweight": "bold",
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(
        0.055,
        0.925,
        args.title,
        fontsize=34,
        fontweight="bold",
        color=COLORS["ink"],
        ha="left",
        va="center",
    )
    fig.text(
        0.055,
        0.870,
        "Each x-axis point is a separately trained model: 90/10, 50/50, or 10/90 train/holdout split.",
        fontsize=17.0,
        color=COLORS["muted"],
        ha="left",
        va="center",
    )

    card_auc = (0.055, 0.515, 0.450, 0.330)
    card_loss = (0.540, 0.515, 0.405, 0.330)
    card_gap = (0.055, 0.105, 0.440, 0.330)
    card_read = (0.530, 0.105, 0.415, 0.330)
    for card in (card_auc, card_loss, card_gap, card_read):
        add_fig_card(fig, card)

    ax_auc = fig.add_axes([0.110, 0.575, 0.360, 0.198])
    ax_auc.plot(X, train_auc, "-o", color=COLORS["orange"], linewidth=3.0, markersize=8)
    ax_auc.plot(X, holdout_auc, "-o", color=COLORS["blue"], linewidth=3.0, markersize=8)
    ax_auc.plot(X, fullstat_auc, "-o", color=COLORS["green"], linewidth=3.0, markersize=8)
    ax_auc.set_title("AUC: train rises; holdout/full sample do not", loc="left", fontsize=18.5, fontweight="bold", color=COLORS["ink"], y=1.12)
    ax_auc.set_ylabel("AUC  (higher is better)")
    ax_auc.set_xticks(X, [f"{s}\n{s.split('/')[0]}% train" for s in SPLITS])
    auc_lo = float(np.nanmin(np.concatenate([train_auc, holdout_auc, fullstat_auc])))
    auc_hi = float(np.nanmax(np.concatenate([train_auc, holdout_auc, fullstat_auc])))
    ax_auc.set_ylim(auc_lo - 0.0010, auc_hi + 0.0026)
    ax_auc.legend(["Train", "Holdout validation", "Full scored sample"], frameon=False, fontsize=11.5, loc="upper left", ncol=3)
    style_ax(ax_auc)

    ax_loss = fig.add_axes([0.590, 0.575, 0.330, 0.198])
    ax_loss.plot(X, train_logloss, "-o", color=COLORS["orange"], linewidth=3.0, markersize=8)
    ax_loss.plot(X, holdout_logloss, "-o", color=COLORS["blue"], linewidth=3.0, markersize=8)
    ax_loss.set_title("Logloss: train vs holdout for each model", loc="left", fontsize=18.5, fontweight="bold", color=COLORS["ink"], y=1.12)
    ax_loss.set_ylabel("Logloss  (lower is better)")
    ax_loss.set_xticks(X, [f"{s}\n{s.split('/')[0]}% train" for s in SPLITS])
    ax_loss.set_ylim(*padded_limits(np.concatenate([train_logloss, holdout_logloss]), frac=0.28, min_pad=0.0014))
    ax_loss.legend(["Train", "Holdout validation"], frameon=False, fontsize=11.5, loc="upper left", ncol=2)
    style_ax(ax_loss)

    ax_gap = fig.add_axes([0.105, 0.165, 0.355, 0.198])
    width = 0.35
    ax_gap.bar(X - width / 2, auc_gap, width, color=COLORS["blue"], label="AUC gap")
    ax_gap.bar(X + width / 2, logloss_gap, width, color=COLORS["red"], label="Logloss gap")
    ax_gap.set_title("Train-validation gaps", loc="left", fontsize=18.5, fontweight="bold", color=COLORS["ink"], y=1.12)
    ax_gap.set_ylabel("Gap size")
    ax_gap.set_xticks(X, [f"{s}\n{s.split('/')[0]}% train" for s in SPLITS])
    ax_gap.axhline(0, color="#9aa7b8", linewidth=1.0)
    ax_gap.set_ylim(*padded_zero_limits(np.concatenate([auc_gap, logloss_gap]), frac=0.22, min_pad=0.0005))
    ax_gap.legend(frameon=False, fontsize=12, ncol=2, loc="upper left")
    style_ax(ax_gap)
    for i in range(3):
        yspan = ax_gap.get_ylim()[1] - ax_gap.get_ylim()[0]
        yoff = 0.025 * yspan
        for xpos, value, color in [
            (i - width / 2, auc_gap[i], COLORS["blue"]),
            (i + width / 2, logloss_gap[i], COLORS["red"]),
        ]:
            va = "bottom" if value >= 0 else "top"
            text_y = value + yoff if value >= 0 else value - yoff
            ax_gap.text(xpos, text_y, f"{value:.4f}", ha="center", va=va, fontsize=11, color=color, fontweight="bold")

    ax_read = fig.add_axes([0.565, 0.145, 0.345, 0.250])
    ax_read.axis("off")
    ax_read.set_xlim(0, 1)
    ax_read.set_ylim(0, 1)
    ax_read.text(0.00, 0.94, "How to read this slide", fontsize=20.0, fontweight="bold", color=COLORS["ink"], ha="left", va="top")
    ax_read.text(0.00, 0.78, "Each column is one separately trained split model.", fontsize=12.8, color=COLORS["muted"], ha="left", va="top")
    guide_rows = [
        (COLORS["orange"], "Train", "rows used to train that model."),
        (COLORS["blue"], "Holdout validation", "reserved split sample; not seen\nduring training."),
        (COLORS["green"], "Full scored sample", "same 20.19M-row pool; includes\ntrain + holdout."),
        (COLORS["red"], "Gap bars", "size of train-holdout mismatch."),
    ]
    y = 0.63
    for color, label, text in guide_rows:
        ax_read.scatter([0.03], [y], s=92, color=color)
        ax_read.text(0.085, y + 0.02, label, fontsize=13.6, fontweight="bold", color=COLORS["ink"], ha="left", va="top")
        ax_read.text(0.47, y + 0.02, text, fontsize=11.6, color=COLORS["muted"], ha="left", va="top", linespacing=1.05)
        y -= 0.145
    ax_read.plot([0.0, 0.96], [0.115, 0.115], color=COLORS["grid"], lw=1.4)
    ax_read.text(
        0.00,
        0.055,
        "Bottom line: 90/10 remains the clean reference;\n"
        "the reduced-training rows do not buy a cleaner\n"
        "validation picture.",
        fontsize=12.0,
        fontweight="bold",
        color=COLORS["ink"],
        ha="left",
        va="top",
    )

    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(out.resolve())


if __name__ == "__main__":
    main()
