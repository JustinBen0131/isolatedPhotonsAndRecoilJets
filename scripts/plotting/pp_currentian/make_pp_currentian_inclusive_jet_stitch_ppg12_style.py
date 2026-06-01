#!/usr/bin/env python3
"""Render the pp current-IAN inclusive-jet stitching plot in PPG12 Fig. 6 style."""

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

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_fullsim_20260521_1811"
BASE = REPO / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN
IN_CSV = BASE / "validation/currentIAN_stitching/ppg12_currentian_truth_spectrum_root_histograms.csv"
OUT = BASE / "slide_assets/pp_currentIAN_inclusive_jet_truth_stitch_ppg12_identical_style.png"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]
COLORS = {
    "jet8": "#d62aa0",
    "jet12": "#2ca02c",
    "jet20": "#1296f3",
    "jet30": "#ff6f00",
    "jet40": "#d62aa0",
}


def load_data() -> dict[str, dict[str, np.ndarray]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    with IN_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["group"] == "jet":
                grouped.setdefault(row["sample"], []).append(row)
    out: dict[str, dict[str, np.ndarray]] = {}
    for sample, rows in grouped.items():
        out[sample] = {
            "x": np.array([float(r["bin_center"]) for r in rows]),
            "y": np.array([float(r["density_pb_per_gev"]) for r in rows]),
            "ey": np.array([float(r["density_err_pb_per_gev"]) for r in rows]),
            "used": np.array([int(r["used_in_stitch"]) for r in rows], dtype=bool),
        }
    return out


def combined(data: dict[str, dict[str, np.ndarray]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = data[SAMPLES[0]]["x"]
    y = np.zeros_like(x)
    e2 = np.zeros_like(x)
    for sample in SAMPLES:
        d = data[sample]
        used = d["used"]
        y[used] += d["y"][used]
        e2[used] += d["ey"][used] ** 2
    return x, y, np.sqrt(e2)


def fit_curve(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = (x >= 10) & (x <= 50) & (y > 0) & np.isfinite(y)
    coeff = np.polyfit(np.log(x[mask]), np.log(y[mask]), deg=4)
    grid = np.linspace(9.5, 50.0, 800)
    pred = np.exp(np.polyval(coeff, np.log(grid)))
    return grid, pred


def main() -> None:
    data = load_data()
    x, y, ey = combined(data)
    grid, pred = fit_curve(x, y)
    fit_y = np.interp(x, grid, pred, left=np.nan, right=np.nan)
    ratio = np.divide(y, fit_y, out=np.full_like(y, np.nan), where=(fit_y > 0) & (y > 0))
    ratio_err = np.divide(ey, fit_y, out=np.full_like(ey, np.nan), where=(fit_y > 0) & (y > 0))

    # PPG12 Fig. 6 labels the weighted histogram as counts. Scale our compact
    # pb/GeV spectrum onto the same visual decade range; MC/fit is unchanged.
    visual_scale = 1.0e12 / np.nanmax(y[(x >= 9) & (x <= 12)])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.4,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 8,
            "ytick.major.size": 8,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )

    fig = plt.figure(figsize=(7.2, 8.1), dpi=180)
    ax = fig.add_axes([0.14, 0.35, 0.79, 0.58])
    rax = fig.add_axes([0.14, 0.095, 0.79, 0.25], sharex=ax)

    ax.set_yscale("log")
    ax.set_xlim(8, 50)
    ax.set_ylim(1e4, 1e12)
    rax.set_ylim(0.85, 1.15)

    fit_top = pred * visual_scale
    ax.plot(grid, fit_top, color="red", lw=2.2)
    for sample in SAMPLES:
        d = data[sample]
        mask = d["used"] & (d["x"] >= 8) & (d["x"] <= 50) & (d["y"] > 0)
        ax.errorbar(
            d["x"][mask],
            d["y"][mask] * visual_scale,
            yerr=d["ey"][mask] * visual_scale,
            fmt="o",
            ms=4.4,
            lw=1.0,
            color=COLORS[sample],
            label=sample,
        )

    rmask = (x >= 8) & (x <= 50) & np.isfinite(ratio) & (y > 0)
    rax.axhline(1.0, color="black", lw=1.0, ls=(0, (6, 6)))
    rax.errorbar(x[rmask], ratio[rmask], yerr=ratio_err[rmask], fmt="o", ms=4.2, lw=1.0, color="black")

    ax.set_ylabel("counts", fontsize=20, fontweight="bold")
    rax.set_ylabel("MC / Fit", fontsize=20, fontweight="bold")
    rax.set_xlabel(r"Leading   $p_T^{\mathrm{jet}}$ [GeV]", fontsize=20, fontweight="bold")
    ax.tick_params(labelsize=17, which="both")
    rax.tick_params(labelsize=17, which="both")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.minorticks_on()
    rax.minorticks_on()

    ax.text(
        0.47,
        0.98,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p{+}p$ $\sqrt{s}=200$ GeV" + "\nPYTHIA8",
        transform=ax.transAxes,
        fontsize=15,
        va="top",
    )
    ax.legend(
        frameon=False,
        fontsize=13,
        loc="lower left",
        bbox_to_anchor=(0.18, 0.18),
        ncol=2,
        columnspacing=1.0,
        handlelength=1.6,
        handletextpad=0.45,
        borderaxespad=0.0,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=180)
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
