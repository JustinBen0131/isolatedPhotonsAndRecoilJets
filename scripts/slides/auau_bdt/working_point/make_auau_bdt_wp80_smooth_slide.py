#!/usr/bin/env python3
"""Render an audience-facing slide for a continuous Au+Au WP80 BDT cut surface."""

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
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import SmoothBivariateSpline


INK = "#111827"
MUTED = "#4B5563"
GRID = "#E5E7EB"
PANEL = "#F8FAFC"
PANEL_EDGE = "#CBD5E1"
YELLOW_TINT = "#FFF4C7"
YELLOW_EDGE = "#F4D35E"
GREEN_TINT = "#E8F6EF"
GREEN_EDGE = "#A7D8BF"
BLUE_TINT = "#EAF2FB"
BLUE_EDGE = "#9EC5E8"
BLUE = "#0072B2"
RED = "#CC334E"
PURPLE = "#7A3E9D"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


def fit_surface(cells: list[dict], n_cent_bins: int) -> tuple[SmoothBivariateSpline, dict]:
    et = np.array([float(r["pt_center"]) for r in cells], dtype=float)
    cent = np.array([0.5 * (float(r["centrality_min"]) + float(r["centrality_max"])) for r in cells], dtype=float)
    cut = np.array([float(r["threshold"]) for r in cells], dtype=float)
    weights = np.sqrt(np.array([max(1, int(r["signal_entries"])) for r in cells], dtype=float))
    weights = weights / weights.mean()

    ky = min(3, max(1, n_cent_bins - 1))
    # Keep the fit smooth enough to suppress bin noise while still following
    # the measured WP80 anchors at analysis-bin centers.
    spline = SmoothBivariateSpline(et, cent, cut, w=weights, kx=3, ky=ky, s=0.04)
    pred = spline.ev(et, cent)
    residual = cut - pred
    return spline, {
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rms_residual": float(np.sqrt(np.mean(residual * residual))),
        "residuals": residual,
        "et": et,
        "cent": cent,
        "cut": cut,
        "pred": pred,
        "ky": ky,
    }


def draw(payload: dict, out: Path) -> None:
    meta = payload["metadata"]
    cells = payload["cells"]
    target = float(meta["target_signal_efficiency"])
    pt_edges = [float(x) for x in meta["pt_edges"]]
    et_centers = np.array([(pt_edges[i] + pt_edges[i + 1]) / 2.0 for i in range(len(pt_edges) - 1)], dtype=float)
    cent_edges = [float(x) for x in meta["cent_edges"]]
    cent_labels = [f"{cent_edges[i]:.0f}-{cent_edges[i + 1]:.0f}%" for i in range(len(cent_edges) - 1)]
    cent_centers = np.array(
        [(cent_edges[i] + cent_edges[i + 1]) / 2.0 for i in range(len(cent_edges) - 1)],
        dtype=float,
    )
    et_labels = [f"{pt_edges[i]:.0f}-{pt_edges[i + 1]:.0f}" for i in range(len(pt_edges) - 1)]

    spline, fit = fit_surface(cells, len(cent_centers))
    smooth_cuts = np.array([[float(spline.ev(et, cent)) for et in et_centers] for cent in cent_centers], dtype=float)
    raw_cuts = np.array(
        [
            [float(r["threshold"]) for r in sorted([c for c in cells if c["centrality_label"] == label], key=lambda x: x["pt_min"])]
            for label in cent_labels
        ],
        dtype=float,
    )
    fakes = np.array(
        [
            [float(r["background_fake_rate"]) for r in sorted([c for c in cells if c["centrality_label"] == label], key=lambda x: x["pt_min"])]
            for label in cent_labels
        ],
        dtype=float,
    )
    effs = np.array([float(r["signal_efficiency"]) for r in cells], dtype=float)
    total_sig = sum(int(r["signal_entries"]) for r in cells)
    total_bkg = sum(int(r["background_entries"]) for r in cells)
    min_sig = min(int(r["signal_entries"]) for r in cells)
    min_bkg = min(int(r["background_entries"]) for r in cells)
    max_eff_dev = float(np.max(np.abs(effs - target)))
    fake_ranges = [(float(np.min(row)), float(np.max(row))) for row in fakes]
    overall_fake_range = (float(np.min(fakes)), float(np.max(fakes)))
    rep_et_idx = [0, 3, 6, 7] if len(et_labels) >= 8 else list(range(len(et_labels)))

    grid_et = np.linspace(pt_edges[0], pt_edges[-1], 240)
    grid_cent = np.linspace(0, 80, 180)
    grid_x, grid_y = np.meshgrid(grid_et, grid_cent)
    grid_z = spline.ev(grid_x.ravel(), grid_y.ravel()).reshape(grid_x.shape)

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
        "Recommended Au+Au WP80 BDT cut: smooth the measured 80% efficiency map",
        fontsize=21.0,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.046,
        0.918,
        f"Build one continuous threshold T($E_T$, centrality) from {len(cent_centers)} x {len(et_centers)} measured WP80 anchor points.",
        fontsize=13.0,
        color=MUTED,
        ha="left",
        va="top",
    )

    method = fig.add_axes([0.655, 0.640, 0.320, 0.265])
    method.axis("off")
    method.add_patch(plt.Rectangle((0, 0), 1, 1, transform=method.transAxes, facecolor=GREEN_TINT, edgecolor=GREEN_EDGE, linewidth=1.0))
    method.text(0.050, 0.84, "Fitting strategy", fontsize=15.2, fontweight="bold", color="#14532D", ha="left", va="center")
    method.text(
        0.060,
        0.52,
        "1. Measure the 80% signal-efficiency BDT cut\n"
        "   in every $E_T$ x centrality cell.\n"
        "2. Fit a weighted smooth surface through\n"
        "   those measured anchor cuts.\n"
        "3. Freeze the surface only after applying it\n"
        "   back to the cells and verifying closure.",
        fontsize=11.1,
        color=INK,
        ha="left",
        va="center",
        linespacing=1.24,
    )
    method.text(
        0.060,
        0.105,
        "Calibrated from measured cuts;\n"
        "smoothness only regularizes bin-to-bin noise.",
        fontsize=9.8,
        color=MUTED,
        ha="left",
        va="center",
        linespacing=1.10,
    )

    ax = fig.add_axes([0.060, 0.390, 0.555, 0.425])
    zmin = float(np.floor((min(np.min(grid_z), np.min(raw_cuts)) - 0.005) * 100) / 100)
    zmax = float(np.ceil((max(np.max(grid_z), np.max(raw_cuts)) + 0.005) * 100) / 100)
    levels = np.linspace(zmin, zmax, 15)
    cf = ax.contourf(grid_x, grid_y, grid_z, levels=levels, cmap="viridis", extend="both")
    cs = ax.contour(grid_x, grid_y, grid_z, levels=np.arange(zmin, zmax + 0.001, 0.02), colors="white", linewidths=0.8, alpha=0.72)
    ax.clabel(cs, inline=True, fontsize=9.2, fmt="%.2f")
    sc = ax.scatter(fit["et"], fit["cent"], c=fit["cut"], cmap="viridis", vmin=zmin, vmax=zmax, s=78, edgecolor="black", linewidth=0.8, zorder=5)
    ax.set_xlim(pt_edges[0], pt_edges[-1])
    ax.set_ylim(80, 0)
    ax.set_xlabel(r"cluster $E_T$ [GeV]", fontsize=13.4)
    ax.set_ylabel("centrality percentile", fontsize=13.4)
    ax.set_title("Measured WP80 anchors -> smoothed BDT threshold", fontsize=15.6, fontweight="bold", pad=9)
    ax.set_yticks(cent_centers)
    ax.set_yticklabels(cent_labels, fontsize=11.0, fontweight="bold")
    ax.tick_params(labelsize=11.3)
    ax.text(
        0.020,
        0.035,
        "black circles = measured 80% signal-efficiency cuts\nwhite contours = fitted BDT-score threshold",
        transform=ax.transAxes,
        fontsize=10.2,
        color=INK,
        ha="left",
        va="bottom",
        bbox=dict(facecolor="white", edgecolor=PANEL_EDGE, alpha=0.90, boxstyle="round,pad=0.32"),
    )
    cax = fig.add_axes([0.165, 0.292, 0.340, 0.020])
    cb = fig.colorbar(cf, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=9.8, length=3, pad=2)
    fig.text(0.335, 0.322, "fitted BDT-score threshold", fontsize=10.4, color=MUTED, ha="center", va="bottom")

    qual = fig.add_axes([0.655, 0.450, 0.320, 0.145])
    qual.axis("off")
    qual.add_patch(plt.Rectangle((0, 0), 1, 1, transform=qual.transAxes, facecolor=BLUE_TINT, edgecolor=BLUE_EDGE, linewidth=1.0))
    qual.text(0.050, 0.78, "Plot inputs", fontsize=13.8, fontweight="bold", color=INK, ha="left", va="center")
    qual.text(
        0.060,
        0.40,
        f"{len(cells)} WP80 anchor cuts from 109 score-cache shards\n"
        "15 <= cluster $E_T$ < 35 GeV, centrality 0-80%\n"
        f"Total entries: S={total_sig:,}, B={total_bkg:,}",
        fontsize=10.6,
        color=INK,
        ha="left",
        va="center",
        linespacing=1.20,
    )

    check = fig.add_axes([0.655, 0.305, 0.320, 0.105])
    check.axis("off")
    check.add_patch(plt.Rectangle((0, 0), 1, 1, transform=check.transAxes, facecolor=PANEL, edgecolor=PANEL_EDGE, linewidth=1.0))
    check.text(0.050, 0.74, "Current closure check", fontsize=13.0, fontweight="bold", color=INK, ha="left", va="center")
    check.text(
        0.060,
        0.33,
        f"Anchor residual: max {fit['max_abs_residual']:.3f}, RMS {fit['rms_residual']:.3f} score\n"
        f"Quantile closure: max |eff - 0.800| = {max_eff_dev:.1e}",
        fontsize=10.0,
        color=INK,
        ha="left",
        va="center",
        linespacing=1.18,
    )

    table = fig.add_axes([0.060, 0.018, 0.555, 0.245])
    table.axis("off")
    table.add_patch(plt.Rectangle((0, 0), 1, 1, transform=table.transAxes, facecolor=PANEL, edgecolor=PANEL_EDGE, linewidth=1.0))
    table.text(0.030, 0.865, "Selected fitted cut values", fontsize=14.0, fontweight="bold", color=INK, ha="left", va="center")
    table.text(0.030, 0.735, "Representative $E_T$ bins shown; full 7x8 table is written next to the PNG.", fontsize=10.3, color=MUTED, ha="left", va="center")
    x0, dx = 0.365, 0.140
    for out_j, j in enumerate(rep_et_idx):
        label = et_labels[j]
        table.text(x0 + out_j * dx, 0.615, label, fontsize=10.3, color=MUTED, ha="center", va="center")
    yrows = np.linspace(0.525, 0.085, len(cent_labels))
    for i, (label, y) in enumerate(zip(cent_labels, yrows)):
        table.text(0.060, y, label, fontsize=10.6, fontweight="bold", color=INK, ha="left", va="center")
        for out_j, j in enumerate(rep_et_idx):
            val = smooth_cuts[i, j]
            table.text(x0 + out_j * dx, y, f"{val:.3f}", fontsize=10.6, color=INK, ha="center", va="center")

    take = fig.add_axes([0.655, 0.018, 0.320, 0.245])
    take.axis("off")
    take.add_patch(plt.Rectangle((0, 0), 1, 1, transform=take.transAxes, facecolor=YELLOW_TINT, edgecolor=YELLOW_EDGE, linewidth=1.0))
    take.text(0.050, 0.78, "Audience-facing claim", fontsize=14.2, fontweight="bold", color=INK, ha="left", va="center")
    take.text(
        0.060,
        0.50,
        "A continuous, centrality-aware BDT cut\n"
        "is more defensible than separate noisy\n"
        "linear fits, provided the frozen surface\n"
        "passes cell-by-cell efficiency closure.",
        fontsize=11.0,
        color=INK,
        ha="left",
        va="center",
        linespacing=1.18,
    )
    take.text(
        0.060,
        0.15,
        f"WP80 anchor fake-rate range: {overall_fake_range[0]:.2f}-{overall_fake_range[1]:.2f}.",
        fontsize=10.0,
        color=MUTED,
        ha="left",
        va="center",
        linespacing=1.10,
        wrap=True,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.with_suffix(".csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["centrality_label", "centrality_center"] + [f"cut_{label}_GeV" for label in et_labels])
        for i, label in enumerate(cent_labels):
            writer.writerow([label, f"{cent_centers[i]:.3f}"] + [f"{x:.6f}" for x in smooth_cuts[i]])
    fig.savefig(out, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    draw(json.loads(args.input.read_text()), args.out)
    print(args.out)


if __name__ == "__main__":
    main()
