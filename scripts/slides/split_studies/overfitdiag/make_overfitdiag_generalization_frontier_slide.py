#!/usr/bin/env python3
"""One-plot verdict slide for the overfitdiag split study."""

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


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--internal-metrics-json", type=Path, default=None)
    ap.add_argument("--compact-json", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=OUTPNG)
    ap.add_argument("--title", default="90/10 is the clean split choice")
    ap.add_argument(
        "--subtitle",
        default=(
            "A high training AUC alone is not success; the check is whether that gain "
            "survives on holdout rows the model never trained on."
        ),
    )
    return ap.parse_args()


def load_metrics(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    train_auc = TRAIN_AUC.copy()
    holdout_auc = HOLDOUT_AUC.copy()
    full_auc = FULL_AUC.copy()
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
        full_auc = np.asarray([rows[s]["global_auc"] for s in SPLITS], dtype=float)

    return train_auc, holdout_auc, full_auc, train_logloss, holdout_logloss


def padded_limits(values: np.ndarray, frac: float = 0.22, min_pad: float = 0.0005) -> tuple[float, float]:
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    pad = max(min_pad, (hi - lo) * frac)
    return lo - pad, hi + pad


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
    args = parse_args()
    out = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    train_auc, holdout_auc, full_auc, train_logloss, holdout_logloss = load_metrics(args)
    auc_gap = train_auc - holdout_auc
    logloss_gap = holdout_logloss - train_logloss
    best_idx = int(np.nanargmax(holdout_auc))
    smallest_auc_gap_idx = int(np.nanargmin(auc_gap))
    smallest_logloss_gap_idx = int(np.nanargmin(logloss_gap))

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
        args.subtitle,
        fontsize=17.0,
        color=COLORS["muted"],
        ha="left",
        va="center",
    )

    plot_card = (0.030, 0.105, 0.675, 0.705)
    verdict_card = (0.735, 0.105, 0.210, 0.705)
    add_card(fig, plot_card)
    add_card(fig, verdict_card)

    ax = fig.add_axes([0.104, 0.180, 0.548, 0.550])
    ax.set_facecolor("white")

    point_colors = [COLORS["good"], COLORS["mid"], COLORS["bad"]]
    marker_size = 420

    ax.plot(
        auc_gap,
        holdout_auc,
        color="#7c8798",
        linewidth=2.2,
        linestyle=(0, (4, 3)),
        zorder=2,
    )
    xlim = padded_limits(auc_gap, frac=0.28, min_pad=0.00045)
    ylim = padded_limits(holdout_auc, frac=0.32, min_pad=0.00045)
    xspan = xlim[1] - xlim[0]
    yspan = ylim[1] - ylim[0]
    ax.text(
        0.975,
        0.965,
        "Less training data: 90/10 -> 50/50 -> 10/90",
        transform=ax.transAxes,
        color="#334155",
        fontsize=15.8,
        ha="right",
        va="top",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.82),
        zorder=7,
    )

    for i, split in enumerate(SPLITS):
        ax.scatter(
            auc_gap[i],
            holdout_auc[i],
            s=marker_size,
            color=point_colors[i],
            edgecolor="white",
            linewidth=2.4,
            zorder=5,
        )
        if split == "90/10":
            label = "90/10"
            xytext = (auc_gap[i], holdout_auc[i] + 0.060 * yspan)
            ha = "center"
            va = "bottom"
        elif split == "50/50":
            label = "50/50"
            xytext = (auc_gap[i] + 0.030 * xspan, holdout_auc[i] - 0.070 * yspan)
            ha = "left"
            va = "top"
        else:
            label = "10/90"
            xytext = (auc_gap[i] - 0.025 * xspan, holdout_auc[i] + 0.065 * yspan)
            ha = "right"
            va = "bottom"
        ax.text(
            xytext[0],
            xytext[1],
            label,
            ha=ha,
            va=va,
            fontsize=17.4,
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
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
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
        color=point_colors[best_idx],
        ha="left",
        va="top",
    )

    rows = [
        ("Best holdout AUC", f"{SPLITS[best_idx]}  {holdout_auc[best_idx]:.4f}"),
        ("Smallest AUC gap", f"{SPLITS[smallest_auc_gap_idx]}  {auc_gap[smallest_auc_gap_idx]:.4f}"),
        ("Smallest logloss gap", f"{SPLITS[smallest_logloss_gap_idx]}  {logloss_gap[smallest_logloss_gap_idx]:+.4f}"),
    ]
    y = 0.660
    for metric, value in rows:
        ax_note.scatter([0.020], [y - 0.006], s=58, color=point_colors[best_idx])
        ax_note.text(0.090, y, metric, fontsize=14.4, fontweight="bold", color=COLORS["ink"], ha="left", va="top")
        ax_note.text(0.090, y - 0.052, value, fontsize=14.2, color=COLORS["muted"], ha="left", va="top")
        y -= 0.155

    ax_note.plot([0.0, 1.0], [0.250, 0.250], color=COLORS["grid"], lw=1.4)
    ax_note.text(0.0, 0.205, "Overfitting call", fontsize=17.2, fontweight="bold", color=COLORS["ink"], ha="left", va="top")
    ax_note.text(
        0.0,
        0.135,
        "The reduced-training rows do not\nbuy a cleaner validation frontier.\nThe full-stat AUC changes only weakly.",
        fontsize=14.0,
        color=COLORS["muted"],
        ha="left",
        va="top",
    )

    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(out.resolve())


if __name__ == "__main__":
    main()
