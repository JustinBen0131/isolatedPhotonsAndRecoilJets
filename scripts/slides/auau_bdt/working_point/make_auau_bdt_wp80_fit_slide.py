#!/usr/bin/env python3
"""Render a full-slide PNG for Au+Au BDT WP80 threshold fits."""

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
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


INK = "#111827"
MUTED = "#4B5563"
GRID = "#E5E7EB"
RED = "#CC334E"
BLUE = "#0072B2"
GREEN = "#009E73"
YELLOW_TINT = "#FFF4C7"
YELLOW_EDGE = "#F4D35E"
BLUE_TINT = "#EAF2FF"
BLUE_EDGE = "#BFD7FF"
PINK_TINT = "#FDE7F1"
PINK_EDGE = "#F4B9D4"
ORANGE = "#D55E00"
PURPLE = "#7A3E9D"


PPG12_NOTE_LINE = {
    "label": "PPG12 note Table 5",
    "intercept": 0.8156,
    "slope": -0.00156,
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


def add_card(fig, xywh, title, body, *, face, edge, title_color=INK, body_color=MUTED):
    ax = fig.add_axes(xywh)
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=face, edgecolor=edge, linewidth=1.0))
    ax.text(0.035, 0.69, title, fontsize=13.5, fontweight="bold", color=title_color, ha="left", va="center")
    ax.text(0.035, 0.28, body, fontsize=11.6, color=body_color, ha="left", va="center", linespacing=1.15)


def ppg12_fit(rows: list[dict]) -> dict:
    """PPG12-style fit: equal-weight linear fit to the per-ET WP thresholds."""
    good = [r for r in rows if r["source"] == "cell" and math.isfinite(float(r["threshold"]))]
    if len(good) < 2:
        return {"slope": math.nan, "intercept": math.nan, "max_abs_residual": math.nan, "rms_residual": math.nan}
    x = np.array([float(r["pt_center"]) for r in good], dtype=float)
    y = np.array([float(r["threshold"]) for r in good], dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    resid = y - pred
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "max_abs_residual": float(np.max(np.abs(resid))),
        "rms_residual": float(np.sqrt(np.mean(resid * resid))),
    }


def line_value(line: dict, et: float) -> float:
    return float(line["intercept"]) + float(line["slope"]) * et


def draw(payload: dict, out: Path) -> None:
    meta = payload["metadata"]
    cells = payload["cells"]
    target = float(meta["target_signal_efficiency"])
    cent_edges = [float(x) for x in meta["cent_edges"]]
    pt_edges = [float(x) for x in meta["pt_edges"]]
    ppg12_fits = {}
    for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
        key = f"{int(clo)}_{int(chi)}"
        cent_rows = [r for r in cells if float(r["centrality_min"]) == clo and float(r["centrality_max"]) == chi]
        ppg12_fits[key] = ppg12_fit(cent_rows)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "mathtext.fontset": "dejavuserif",
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=220, facecolor="white")
    fig.text(
        0.045,
        0.965,
        r"PPG12-style WP80 propagation exposes centrality-dependent Au+Au BDT cuts",
        fontsize=24.0,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.046,
        0.918,
        rf"Verbatim method: choose the BDT cutoff in each $E_T$ bin to retain {100*target:.0f}% signal, then fit cutoff = $a+bE_T$ separately in each centrality bin.",
        fontsize=13.2,
        color=MUTED,
        ha="left",
        va="top",
    )
    fig.text(0.805, 0.955, r"$\bf{\it{sPHENIX}}$ Internal", fontsize=13.6, ha="left", va="top", color=INK)
    fig.text(0.805, 0.925, "PYTHIA8 Au+Au embedded validation", fontsize=10.8, ha="left", va="top", color=INK)

    add_card(fig, [0.045, 0.828, 0.270, 0.060], "Model shown", meta["model_label"], face=PINK_TINT, edge=PINK_EDGE)
    add_card(fig, [0.328, 0.828, 0.300, 0.060], "Training sample", meta["training_sample"], face="#F8FAFC", edge="#CBD5E1")
    add_card(
        fig,
        [0.642, 0.828, 0.230, 0.060],
        "Training inputs",
        meta["training_inputs"],
        face=BLUE_TINT,
        edge=BLUE_EDGE,
        title_color="#1E3A8A",
        body_color="#1E3A8A",
    )

    legend_ax = fig.add_axes([0.875, 0.797, 0.120, 0.108])
    legend_ax.axis("off")
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.plot(
        [0.08],
        [0.78],
        color=INK,
        lw=0,
        marker="o",
        markerfacecolor="white",
        markeredgecolor=INK,
        markersize=6.5,
    )
    legend_ax.text(0.32, 0.78, "WP80 cells", fontsize=10.8, color=INK, ha="left", va="center", clip_on=False)
    legend_ax.plot([0.02, 0.26], [0.48, 0.48], color=INK, lw=2.3, ls="--")
    legend_ax.text(0.32, 0.48, "Au+Au linear fit", fontsize=10.3, color=INK, ha="left", va="center", clip_on=False)
    legend_ax.plot([0.02, 0.26], [0.19, 0.19], color=PURPLE, lw=2.0, ls=":")
    legend_ax.text(0.32, 0.19, "PPG12 note Table 5", fontsize=10.3, color=INK, ha="left", va="center", clip_on=False)

    left, bottom, width, height = 0.058, 0.285, 0.890, 0.475
    xgap = 0.045
    panel_w = (width - 2 * xgap) / 3
    colors = [RED, GREEN, BLUE]
    axes = [
        fig.add_axes([left + i * (panel_w + xgap), bottom, panel_w, height])
        for i in range(3)
    ]
    thresholds = [float(r["threshold"]) for r in cells if math.isfinite(float(r["threshold"]))]
    ref_vals = [line_value(PPG12_NOTE_LINE, x) for x in pt_edges]
    ylo = max(0.0, min(min(thresholds), min(ref_vals)) - 0.030)
    yhi = min(1.0, max(max(thresholds), max(ref_vals)) + 0.035)
    xfit = np.linspace(pt_edges[0], pt_edges[-1], 160)

    for idx, (ax, clo, chi) in enumerate(zip(axes, cent_edges[:-1], cent_edges[1:])):
        rows = [r for r in cells if float(r["centrality_min"]) == clo and float(r["centrality_max"]) == chi]
        x = np.array([float(r["pt_center"]) for r in rows], dtype=float)
        y = np.array([float(r["threshold"]) for r in rows], dtype=float)
        eff = np.array([float(r["signal_efficiency"]) for r in rows], dtype=float)
        sig_n = np.array([max(1, int(r["signal_entries"])) for r in rows], dtype=float)
        yerr = np.sqrt(np.maximum(eff * (1.0 - eff), 0.0) / sig_n)
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            color=colors[idx],
            markeredgecolor="black",
            markeredgewidth=0.5,
            markersize=6.5,
            capsize=2.5,
            linestyle="none",
        )
        key = f"{int(clo)}_{int(chi)}"
        fit = ppg12_fits[key]
        fit_y = float(fit["slope"]) * xfit + float(fit["intercept"])
        ax.plot(xfit, fit_y, color=INK, linewidth=2.0, linestyle="--")
        ax.plot(
            xfit,
            line_value(PPG12_NOTE_LINE, xfit),
            color=PURPLE,
            linewidth=1.8,
            linestyle=":",
            alpha=0.90,
        )
        fake_vals = [float(r["background_fake_rate"]) for r in rows]
        ax.text(
            0.575,
            0.705,
            f"Au+Au fit: {fit['intercept']:.3f} {fit['slope']:+.4f} $E_T$\nmax |resid| {fit['max_abs_residual']:.3f}\nbkg fake {min(fake_vals):.2f}-{max(fake_vals):.2f}",
            transform=ax.transAxes,
            fontsize=9.0,
            ha="left",
            va="top",
            color=INK,
            bbox={"facecolor": "white", "edgecolor": "#D1D5DB", "alpha": 0.95, "pad": 3.2},
        )
        ax.set_title(f"{int(clo)}-{int(chi)}% centrality", fontsize=15.0, fontweight="bold", pad=8)
        ax.set_xlim(pt_edges[0] - 0.5, pt_edges[-1] + 0.5)
        ax.set_ylim(ylo, yhi)
        ax.grid(True, color=GRID, linewidth=0.65)
        ax.tick_params(labelsize=11.2, length=4.5, width=0.9)
        ax.set_xlabel(r"cluster $E_T$ bin center [GeV]", fontsize=12.6)
        if idx == 0:
            ax.set_ylabel("BDT score threshold for WP80", fontsize=12.6)
        else:
            ax.set_yticklabels([])

    residuals = [float(ppg12_fits[f"{int(clo)}_{int(chi)}"]["max_abs_residual"]) for clo, chi in zip(cent_edges[:-1], cent_edges[1:])]
    inclusive = meta.get("inclusive", {})
    compare = fig.add_axes([0.055, 0.025, 0.540, 0.160])
    compare.axis("off")
    compare.add_patch(plt.Rectangle((0, 0), 1, 1, transform=compare.transAxes, facecolor="#F8FAFC", edgecolor="#CBD5E1", linewidth=1.0))
    compare.text(0.025, 0.86, "PPG12-style method output", fontsize=12.7, fontweight="bold", color=INK, ha="left", va="center")
    compare.text(0.435, 0.68, "linear cut a + b ET", fontsize=10.2, color=MUTED, ha="center", va="center")
    compare.text(0.675, 0.68, "@15", fontsize=10.2, color=MUTED, ha="center", va="center")
    compare.text(0.800, 0.68, "@25", fontsize=10.2, color=MUTED, ha="center", va="center")
    compare.text(0.925, 0.68, "@35", fontsize=10.2, color=MUTED, ha="center", va="center")
    table_rows = []
    for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
        table_rows.append((f"Au+Au {int(clo)}-{int(chi)}%", ppg12_fits[f"{int(clo)}_{int(chi)}"], INK, "bold"))
    table_rows.extend(
        [
            ("PPG12 note T5", PPG12_NOTE_LINE, PURPLE, "normal"),
        ]
    )
    yrows = [0.52, 0.385, 0.250, 0.095]
    for yrow, (label, fit, color, weight) in zip(yrows, table_rows):
        compare.text(0.035, yrow, label, fontsize=10.2, fontweight=weight, color=color, ha="left", va="center")
        compare.text(0.435, yrow, f"{fit['intercept']:.3f} {fit['slope']:+.4f} ET", fontsize=9.9, color=color, ha="center", va="center")
        for xpos, et in [(0.675, 15.0), (0.800, 25.0), (0.925, 35.0)]:
            compare.text(xpos, yrow, f"{line_value(fit, et):.3f}", fontsize=9.9, color=color, ha="center", va="center")

    refs = fig.add_axes([0.620, 0.025, 0.343, 0.160])
    refs.axis("off")
    refs.add_patch(plt.Rectangle((0, 0), 1, 1, transform=refs.transAxes, facecolor=YELLOW_TINT, edgecolor=YELLOW_EDGE, linewidth=1.0))
    refs.text(0.030, 0.82, "Problem with inheriting a pp line", fontsize=12.6, fontweight="bold", color=INK, ha="left", va="center")
    refs.text(
        0.030,
        0.56,
        (
            "The PPG12 note uses one pp ET-dependent BDT line;\n"
            "the Au+Au WP80 line moves with centrality\n"
            "and changes slope sign."
        ),
        fontsize=10.0,
        color=INK,
        ha="left",
        va="center",
        linespacing=1.12,
        wrap=True,
    )
    refs.text(
        0.030,
        0.31,
        (
            f"Max linear residual is {max(residuals):.3f}; use a cell map\n"
            "or validated smooth 2D calibration, not the pp cut."
        ),
        fontsize=10.0,
        color="#7A4A00",
        ha="left",
        va="center",
        linespacing=1.15,
        wrap=True,
    )
    refs.text(
        0.030,
        0.09,
        f"Rows: {meta['rows_loaded']:,}; inclusive WP80 cut {float(inclusive.get('threshold', float('nan'))):.3f}, bkg fake {float(inclusive.get('background_fake_rate', float('nan'))):.3f}.",
        fontsize=9.0,
        color=MUTED,
        ha="left",
        va="center",
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    draw(json.loads(args.input.read_text()), args.out)
    print(args.out)


if __name__ == "__main__":
    main()
