#!/usr/bin/env python3
"""Render compact/full stack algorithm score-separation comparison slide."""

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
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


def step_xy(edges, density):
    xs = []
    ys = []
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
    ax.text(0.035, 0.68, title, ha="left", va="center", fontsize=12.6, fontweight="bold", color=INK)
    ax.text(0.035, 0.29, body, ha="left", va="center", fontsize=10.7, color=MUTED, linespacing=1.12)


def main() -> None:
    args = parse_args()
    payload = json.loads(args.input.read_text())
    panels = payload["panels"]
    panel_map = {(p["family"], p["algorithm"]): p for p in panels}
    families = ["compact", "full features"]
    algorithms = ["NN", "Logistic", "GBM"]
    family_labels = {"compact": "compact\nscores + E_T + cent", "full features": "full features\nscores + baseV3E"}

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
    fig.text(0.045, 0.962, "Stack combiner check: simple models are competitive", ha="left", va="top", fontsize=25.0, fontweight="bold", color=INK)
    fig.text(
        0.046,
        0.922,
        r"Same upstream BDT and MLP scores; panels compare stack algorithm choices in 0-20% central Au+Au, $15<E_T<35$ GeV.",
        ha="left",
        va="top",
        fontsize=13.0,
        color=MUTED,
    )
    fig.text(0.805, 0.962, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=13.5, color=INK)
    fig.text(0.805, 0.934, "PYTHIA8 Au+Au embedded validation", ha="left", va="top", fontsize=10.8, color=INK)

    add_card(fig, [0.046, 0.835, 0.260, 0.060], "Inputs held fixed", "BDT score + MLP score are identical across rows", BLUE, BLUE_EDGE)
    add_card(fig, [0.322, 0.835, 0.260, 0.060], "Question", "does a more complex NN combiner beat logistic or GBM?", PINK, PINK_EDGE)
    add_card(fig, [0.598, 0.835, 0.230, 0.060], "Displayed split", "held-out test rows, centrality 0-20%", GREEN, GREEN_EDGE)

    leg = fig.add_axes([0.842, 0.835, 0.145, 0.060])
    leg.axis("off")
    leg.set_xlim(0, 1)
    leg.set_ylim(0, 1)
    leg.plot([0.02, 0.22], [0.68, 0.68], color=SIGNAL, lw=4, solid_capstyle="butt")
    leg.text(0.28, 0.68, "Signal", fontsize=11.0, va="center", ha="left")
    leg.plot([0.02, 0.22], [0.30, 0.30], color=BACKGROUND, lw=4, solid_capstyle="butt")
    leg.text(0.28, 0.30, "Background", fontsize=11.0, va="center", ha="left")

    left = 0.130
    bottom = 0.215
    width = 0.830
    height = 0.565
    xgap = 0.030
    ygap = 0.075
    pw = (width - 2 * xgap) / 3
    ph = (height - ygap) / 2
    ymax = max(max(p["signal_density"] + p["background_density"]) for p in panels) * 1.16
    row_colors = ["#FFEDD5", "#ECFDF5"]

    for ir, fam in enumerate(families):
        y0 = bottom + (len(families) - 1 - ir) * (ph + ygap)
        band = fig.add_axes([0.046, y0, 0.060, ph])
        band.axis("off")
        band.add_patch(plt.Rectangle((0, 0), 1, 1, transform=band.transAxes, facecolor=row_colors[ir], edgecolor="#D1D5DB", linewidth=0.8))
        band.text(0.5, 0.50, family_labels[fam], rotation=90, ha="center", va="center", fontsize=13.1, fontweight="bold", color=INK, linespacing=1.0)
        for ic, alg in enumerate(algorithms):
            ax = fig.add_axes([left + ic * (pw + xgap), y0, pw, ph])
            p = panel_map[(fam, alg)]
            sx, sy = step_xy(p["bin_edges"], p["signal_density"])
            bx, by = step_xy(p["bin_edges"], p["background_density"])
            ax.plot(bx, by, color=BACKGROUND, lw=2.0)
            ax.plot(sx, sy, color=SIGNAL, lw=2.0)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, ymax)
            ax.grid(True, color=GRID, lw=0.52)
            ax.set_axisbelow(True)
            ax.tick_params(labelsize=10.0, length=4.0, width=0.8)
            if ir == 0:
                ax.set_title(alg, fontsize=15.0, fontweight="bold", pad=6)
            if ic == 0:
                ax.set_ylabel("Area-normalized density", fontsize=11.1)
            else:
                ax.set_yticklabels([])
            if ir == len(families) - 1:
                ax.set_xlabel("Stack score", fontsize=12.5, labelpad=4)
            else:
                ax.set_xticklabels([])
            ax.text(
                0.045,
                0.900,
                f"AUC {p['auc']:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=12.0,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#D1D5DB", lw=0.8, alpha=0.95),
            )
            ax.text(
                0.045,
                0.725,
                f"S {p['signal_entries']:,}\nB {p['background_entries']:,}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.2,
                color=MUTED,
                linespacing=1.10,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.76),
            )

    auc = {(p["family"], p["algorithm"]): p["auc"] for p in panels}
    compact_best = max(algorithms, key=lambda alg: auc[("compact", alg)])
    full_best = max(algorithms, key=lambda alg: auc[("full features", alg)])
    band = fig.add_axes([0.046, 0.055, 0.914, 0.102])
    band.axis("off")
    band.add_patch(plt.Rectangle((0, 0), 1, 1, transform=band.transAxes, facecolor=ORANGE, edgecolor=ORANGE_EDGE, linewidth=1.0))
    band.text(0.020, 0.62, "Readout", fontsize=14.5, fontweight="bold", color=INK, ha="left", va="center")
    band.text(
        0.125,
        0.62,
        f"Compact inputs: NN {auc[('compact','NN')]:.3f}, logistic {auc[('compact','Logistic')]:.3f}, GBM {auc[('compact','GBM')]:.3f}; best = {compact_best}.",
        fontsize=12.5,
        color=INK,
        ha="left",
        va="center",
    )
    band.text(
        0.125,
        0.25,
        f"Full-feature inputs: NN {auc[('full features','NN')]:.3f}, logistic {auc[('full features','Logistic')]:.3f}, GBM {auc[('full features','GBM')]:.3f}; best = {full_best}.",
        fontsize=11.4,
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
