#!/usr/bin/env python3
"""Render an audience-facing slide for the recommended Au+Au WP80 cut map."""

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
PANEL = "#F8FAFC"
PANEL_EDGE = "#CBD5E1"
YELLOW_TINT = "#FFF4C7"
YELLOW_EDGE = "#F4D35E"
BLUE_TINT = "#EAF2FF"
BLUE_EDGE = "#BFD7FF"
GREEN_TINT = "#E8F6EF"
GREEN_EDGE = "#A7D8BF"
RED = "#CC334E"
BLUE = "#0072B2"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


def rows_by_centrality(cells: list[dict], labels: list[str]) -> list[list[dict]]:
    out = []
    for label in labels:
        rows = [r for r in cells if r["centrality_label"] == label]
        rows.sort(key=lambda r: float(r["pt_min"]))
        out.append(rows)
    return out


def add_header_card(fig, xywh, title, body, *, face=PANEL, edge=PANEL_EDGE, title_color=INK, body_color=MUTED):
    ax = fig.add_axes(xywh)
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=face, edgecolor=edge, linewidth=1.0))
    ax.text(0.035, 0.68, title, fontsize=13.5, fontweight="bold", color=title_color, ha="left", va="center")
    ax.text(0.035, 0.28, body, fontsize=11.2, color=body_color, ha="left", va="center", linespacing=1.15)


def draw(payload: dict, out: Path) -> None:
    meta = payload["metadata"]
    cells = payload["cells"]
    target = float(meta["target_signal_efficiency"])
    pt_edges = [float(x) for x in meta["pt_edges"]]
    cent_labels = ["0-20%", "20-50%", "50-80%"]
    grouped = rows_by_centrality(cells, cent_labels)

    cuts = np.array([[float(r["threshold"]) for r in rows] for rows in grouped], dtype=float)
    fakes = np.array([[float(r["background_fake_rate"]) for r in rows] for rows in grouped], dtype=float)
    effs = np.array([[float(r["signal_efficiency"]) for r in rows] for rows in grouped], dtype=float)
    sig_counts = np.array([[int(r["signal_entries"]) for r in rows] for rows in grouped], dtype=int)
    bkg_counts = np.array([[int(r["background_entries"]) for r in rows] for rows in grouped], dtype=int)

    et_labels = [f"{pt_edges[i]:.0f}-{pt_edges[i + 1]:.0f}" for i in range(len(pt_edges) - 1)]
    max_eff_dev = float(np.max(np.abs(effs - target)))
    fake_ranges = [(float(np.min(row)), float(np.max(row))) for row in fakes]
    total_sig = int(np.sum(sig_counts))
    total_bkg = int(np.sum(bkg_counts))

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.0,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "mathtext.fontset": "dejavuserif",
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=220, facecolor="white")
    fig.text(
        0.045,
        0.965,
        "Recommended Au+Au WP80 BDT cut: calibrate in the actual analysis bins",
        fontsize=23.2,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.046,
        0.918,
        rf"Use the PPG12 signal-efficiency idea, but make the Au+Au threshold a centrality-aware lookup table rather than a forced pp-style line.",
        fontsize=13.0,
        color=MUTED,
        ha="left",
        va="top",
    )
    fig.text(0.810, 0.955, r"$\bf{\it{sPHENIX}}$ Internal", fontsize=13.6, ha="left", va="top", color=INK)
    fig.text(0.810, 0.925, "PYTHIA8 Au+Au embedded validation", fontsize=10.8, ha="left", va="top", color=INK)

    add_header_card(fig, [0.045, 0.828, 0.278, 0.060], "Nominal prescription", "WP80 lookup in ET x centrality cells", face=GREEN_TINT, edge=GREEN_EDGE, title_color="#14532D", body_color="#14532D")
    add_header_card(fig, [0.338, 0.828, 0.292, 0.060], "Training sample", meta["training_sample"], face=PANEL, edge=PANEL_EDGE)
    add_header_card(fig, [0.645, 0.828, 0.230, 0.060], "Inputs", meta["training_inputs"], face=BLUE_TINT, edge=BLUE_EDGE, title_color="#1E3A8A", body_color="#1E3A8A")

    ax = fig.add_axes([0.065, 0.265, 0.605, 0.500])
    im = ax.imshow(cuts, cmap="viridis", vmin=0.55, vmax=0.69, aspect="auto")
    ax.set_xticks(np.arange(len(et_labels)))
    ax.set_xticklabels(et_labels, fontsize=10.8)
    ax.set_yticks(np.arange(len(cent_labels)))
    ax.set_yticklabels(cent_labels, fontsize=12.4, fontweight="bold")
    ax.set_xlabel(r"cluster $E_T$ bin [GeV]", fontsize=12.6, labelpad=8)
    ax.set_ylabel("centrality", fontsize=12.6, labelpad=10)
    ax.set_title("Recommended BDT score threshold per analysis bin", fontsize=15.5, fontweight="bold", pad=10)
    ax.set_xticks(np.arange(-0.5, len(et_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(cent_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i in range(cuts.shape[0]):
        for j in range(cuts.shape[1]):
            color = "white" if cuts[i, j] < 0.600 or cuts[i, j] > 0.650 else INK
            ax.text(j, i, f"{cuts[i, j]:.3f}", ha="center", va="center", fontsize=11.2, fontweight="bold", color=color)

    cax = fig.add_axes([0.680, 0.350, 0.011, 0.330])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=9.5)
    fig.text(0.675, 0.695, "BDT cut", fontsize=9.8, color=MUTED, ha="left", va="bottom")

    method = fig.add_axes([0.720, 0.502, 0.245, 0.263])
    method.axis("off")
    method.add_patch(plt.Rectangle((0, 0), 1, 1, transform=method.transAxes, facecolor=PANEL, edgecolor=PANEL_EDGE, linewidth=1.0))
    method.text(0.05, 0.84, "Why this is defensible", fontsize=13.4, fontweight="bold", color=INK, ha="left", va="center")
    bullets = [
        rf"Each cell keeps {100*target:.0f}% of signal by construction.",
        "No extrapolation across centrality.",
        "Smooth fits become cross-checks/systematics, not the definition.",
    ]
    y = 0.62
    for bullet in bullets:
        method.text(0.07, y, u"\u2022", fontsize=13.5, color=INK, ha="left", va="center")
        method.text(0.12, y, bullet, fontsize=10.9, color=INK, ha="left", va="center", wrap=True)
        y -= 0.20

    closure = fig.add_axes([0.720, 0.265, 0.245, 0.205])
    closure.axis("off")
    closure.add_patch(plt.Rectangle((0, 0), 1, 1, transform=closure.transAxes, facecolor="#FFFFFF", edgecolor=PANEL_EDGE, linewidth=1.0))
    closure.text(0.05, 0.80, "Closure from current score cache", fontsize=12.8, fontweight="bold", color=INK, ha="left", va="center")
    closure.text(0.05, 0.56, f"max |signal eff - 0.800| = {max_eff_dev:.1e}", fontsize=10.6, color=INK, ha="left", va="center")
    closure.text(0.05, 0.37, f"signal entries = {total_sig:,}", fontsize=10.6, color=INK, ha="left", va="center")
    closure.text(0.05, 0.19, f"background entries = {total_bkg:,}", fontsize=10.6, color=INK, ha="left", va="center")

    table = fig.add_axes([0.065, 0.045, 0.555, 0.145])
    table.axis("off")
    table.add_patch(plt.Rectangle((0, 0), 1, 1, transform=table.transAxes, facecolor=PANEL, edgecolor=PANEL_EDGE, linewidth=1.0))
    table.text(0.03, 0.78, "Background fake rate after the recommended cuts", fontsize=12.6, fontweight="bold", color=INK, ha="left", va="center")
    xs = [0.18, 0.43, 0.68]
    for x, label, (lo, hi) in zip(xs, cent_labels, fake_ranges):
        table.text(x, 0.48, label, fontsize=11.2, fontweight="bold", color=INK, ha="center", va="center")
        table.text(x, 0.22, f"{lo:.2f}-{hi:.2f}", fontsize=17.0, fontweight="bold", color=BLUE, ha="center", va="center")
    table.text(0.88, 0.45, "lower fake rate in\nmore peripheral cells", fontsize=10.7, color=MUTED, ha="center", va="center", linespacing=1.15)

    take = fig.add_axes([0.645, 0.045, 0.320, 0.145])
    take.axis("off")
    take.add_patch(plt.Rectangle((0, 0), 1, 1, transform=take.transAxes, facecolor=YELLOW_TINT, edgecolor=YELLOW_EDGE, linewidth=1.0))
    take.text(0.05, 0.72, "Audience-facing takeaway", fontsize=12.7, fontweight="bold", color=INK, ha="left", va="center")
    take.text(
        0.05,
        0.42,
        "The manageable Au+Au working point is a documented cut map:\n"
        "simple to apply, exact in the chosen bins, and honest about\n"
        "centrality dependence.",
        fontsize=10.8,
        color=INK,
        ha="left",
        va="center",
        linespacing=1.18,
    )
    take.text(0.05, 0.13, "A smooth 2D fit can be added later only after closure criteria are met.", fontsize=9.4, color=MUTED, ha="left", va="center")

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    draw(json.loads(args.input.read_text()), args.out)
    print(args.out)


if __name__ == "__main__":
    main()
