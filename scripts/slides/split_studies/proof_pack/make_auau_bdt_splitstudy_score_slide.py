#!/usr/bin/env python3
"""Render a full-slide BDT train/test split stress-test comparison."""

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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


INK = "#111827"
MUTED = "#4B5563"
GRID = "#E5E7EB"
SIGNAL = "#0072B2"
BACKGROUND = "#CC334E"
PINK = "#FDE7F1"
PINK_EDGE = "#F4B9D4"
ORANGE = "#FFF0D8"
ORANGE_EDGE = "#FDBA74"
GREEN = "#EAF7EA"
GREEN_EDGE = "#86C98A"
BLUE = "#EAF2FF"
BLUE_EDGE = "#BFD7FF"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True, help="Compact JSON from score-cache reduction")
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


def step_xy(edges: list[float], density: list[float]) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    for i, val in enumerate(density):
        lo = edges[i]
        hi = edges[i + 1]
        if i == 0:
            xs.extend([lo, lo])
            ys.extend([0.0, val])
        else:
            xs.extend([lo, lo])
            ys.extend([ys[-1], val])
        xs.append(hi)
        ys.append(val)
    xs.append(edges[-1])
    ys.append(0.0)
    return xs, ys


def add_card(fig, xywh, title, body, face, edge):
    ax = fig.add_axes(xywh)
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=face, edgecolor=edge, linewidth=1.0))
    ax.text(0.035, 0.68, title, ha="left", va="center", fontsize=12.5, fontweight="bold", color=INK)
    ax.text(0.035, 0.29, body, ha="left", va="center", fontsize=10.6, color=MUTED, linespacing=1.12)


def main() -> None:
    args = parse_args()
    payload = json.loads(args.input.read_text())
    rows = payload["rows"]
    cols = payload["centrality_bins"]
    hist = payload["histograms"]
    hist_map = {(h["split_label"], h["cent_label"]): h for h in hist}

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 0.95,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(16, 9), dpi=220, facecolor="white")
    fig.text(0.045, 0.962, "BDT train/validation split stress test", ha="left", va="top", fontsize=25.5, fontweight="bold", color=INK)
    fig.text(
        0.046,
        0.922,
        r"Same global $\bf{centInput\_pt1535}$ feature family; rows change only how much data is used to train the BDT.",
        ha="left",
        va="top",
        fontsize=13.1,
        color=MUTED,
    )
    fig.text(0.805, 0.962, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=13.5, color=INK)
    fig.text(0.805, 0.934, "PYTHIA8 Au+Au embedded validation", ha="left", va="top", fontsize=10.8, color=INK)

    add_card(fig, [0.046, 0.835, 0.255, 0.060], "Model", r"global BDT, baseV3E + centrality + $w_{\eta,33}/w_{\phi,33}$", PINK, PINK_EDGE)
    add_card(fig, [0.315, 0.835, 0.315, 0.060], "Fixed evaluation sample", "same S/B counts shown in every row; only training fraction changes", BLUE, BLUE_EDGE)
    add_card(fig, [0.644, 0.835, 0.185, 0.060], "Weights", r"PPG12-style $E_T$ and $\eta$ weights only", GREEN, GREEN_EDGE)

    leg = fig.add_axes([0.842, 0.835, 0.145, 0.060])
    leg.axis("off")
    leg.set_xlim(0, 1)
    leg.set_ylim(0, 1)
    leg.plot([0.02, 0.22], [0.68, 0.68], color=SIGNAL, lw=4, solid_capstyle="butt")
    leg.text(0.28, 0.68, "Signal", fontsize=11.0, va="center", ha="left")
    leg.plot([0.02, 0.22], [0.30, 0.30], color=BACKGROUND, lw=4, solid_capstyle="butt")
    leg.text(0.28, 0.30, "Background", fontsize=11.0, va="center", ha="left")

    left = 0.155
    bottom = 0.175
    width = 0.805
    height = 0.610
    nrows = len(rows)
    ncols = len(cols)
    xgap = 0.030
    ygap = 0.052
    pw = (width - (ncols - 1) * xgap) / ncols
    ph = (height - (nrows - 1) * ygap) / nrows
    ymax = max(max(h["signal_density"] + h["background_density"]) for h in hist) * 1.16
    row_colors = ["#FCE7F3", "#FFEDD5", "#ECFDF5"]

    for ir, row in enumerate(rows):
        y0 = bottom + (nrows - 1 - ir) * (ph + ygap)
        band = fig.add_axes([0.046, y0, 0.060, ph])
        band.axis("off")
        band.add_patch(plt.Rectangle((0, 0), 1, 1, transform=band.transAxes, facecolor=row_colors[ir], edgecolor="#D1D5DB", linewidth=0.8))
        band.text(0.5, 0.60, row["label"], rotation=90, ha="center", va="center", fontsize=17.0, fontweight="bold", color=INK)
        band.text(0.5, 0.20, "train/test", rotation=90, ha="center", va="center", fontsize=8.6, color=MUTED)
        for ic, col in enumerate(cols):
            ax = fig.add_axes([left + ic * (pw + xgap), y0, pw, ph])
            h = hist_map[(row["label"], col["label"])]
            sx, sy = step_xy(h["bin_edges"], h["signal_density"])
            bx, by = step_xy(h["bin_edges"], h["background_density"])
            ax.plot(bx, by, color=BACKGROUND, lw=1.95)
            ax.plot(sx, sy, color=SIGNAL, lw=1.95)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, ymax)
            ax.grid(True, color=GRID, lw=0.50)
            ax.set_axisbelow(True)
            ax.tick_params(labelsize=9.7, length=4.0, width=0.8)
            if ir == 0:
                ax.set_title(col["label"], fontsize=14.0, fontweight="bold", pad=5)
            if ic == 0:
                ax.set_ylabel("Area-normalized density", fontsize=10.2, labelpad=8)
            else:
                ax.set_yticklabels([])
            if ir == nrows - 1:
                ax.set_xlabel("BDT score", fontsize=12.2, labelpad=4)
            else:
                ax.set_xticklabels([])
            ax.text(
                0.045,
                0.900,
                f"AUC {h['auc']:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=11.2,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#D1D5DB", lw=0.8, alpha=0.95),
            )
            ax.text(
                0.045,
                0.710,
                f"S {h['signal_entries']:,}\nB {h['background_entries']:,}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.9,
                color=MUTED,
                linespacing=1.10,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.76),
            )

    aucs = {r["label"]: r["global_auc"] for r in rows}
    delta_50 = aucs["50/50"] - aucs["90/10"]
    delta_10 = aucs["10/90"] - aucs["90/10"]
    cent_labels = [c["label"] for c in cols]
    centrality_deltas = []
    for label in cent_labels:
        base = hist_map[("90/10", label)]["auc"]
        centrality_deltas.append((hist_map[("50/50", label)]["auc"] - base, hist_map[("10/90", label)]["auc"] - base))
    min_cent_delta = min(min(pair) for pair in centrality_deltas)
    max_cent_delta = max(max(pair) for pair in centrality_deltas)
    band = fig.add_axes([0.046, 0.043, 0.914, 0.092])
    band.axis("off")
    band.add_patch(plt.Rectangle((0, 0), 1, 1, transform=band.transAxes, facecolor=ORANGE, edgecolor=ORANGE_EDGE, linewidth=1.0))
    band.text(0.020, 0.62, "Readout", fontsize=14.5, fontweight="bold", color=INK, ha="left", va="center")
    band.text(
        0.125,
        0.62,
        f"Global full-stat AUC is stable: 90/10 = {aucs['90/10']:.3f}, 50/50 = {aucs['50/50']:.3f}, 10/90 = {aucs['10/90']:.3f}.",
        fontsize=12.4,
        color=INK,
        ha="left",
        va="center",
    )
    band.text(
        0.125,
        0.25,
        f"Within centrality bins, AUC changes only {min_cent_delta:+.4f} to {max_cent_delta:+.4f}; this is a stress test, not proof that less training data is better.",
        fontsize=10.9,
        color=MUTED,
        ha="left",
        va="center",
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220)
    plt.close(fig)
    print(args.out)


if __name__ == "__main__":
    main()
