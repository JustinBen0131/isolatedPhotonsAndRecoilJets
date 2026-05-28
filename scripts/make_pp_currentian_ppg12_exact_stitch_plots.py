#!/usr/bin/env python3
"""Render PPG12-source-matched pp stitch plots for photon and jet spectra."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_fullsim_20260521_1811"
BASE = REPO / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN
STITCH_DIR = BASE / "validation/currentIAN_stitching"
ASSET_DIR = BASE / "slide_assets"

ROOTFIT_CSV = STITCH_DIR / "ppg12_currentian_rootfit_stitch_points.csv"

PHOTON_OUT = ASSET_DIR / "pp_currentIAN_photon_truth_stitch_ppg12_exact_style.png"
JET_OUT = ASSET_DIR / "pp_currentIAN_inclusive_jet_truth_stitch_ppg12_exact_style.png"

COLORS = {
    "photon5": "#d62aa0",
    "photon10": "#2ca02c",
    "photon20": "#1296f3",
    "jet8": "#d62aa0",
    "jet12": "#2ca02c",
    "jet20": "#1296f3",
    "jet30": "#ff6f00",
    "jet40": "#d62aa0",
}


def load_group(path: Path, group: str) -> dict[str, dict[str, np.ndarray]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            if row["group"] == group:
                grouped.setdefault(row["sample"], []).append(row)
    out: dict[str, dict[str, np.ndarray]] = {}
    for sample, rows in grouped.items():
        out[sample] = {
            "x": np.array([float(r["bin_center"]) for r in rows], dtype=float),
            "y": np.array([float(r["value"]) for r in rows], dtype=float),
            "ey": np.array([float(r["error"]) for r in rows], dtype=float),
            "used": np.array([int(r["used_in_stitch"]) for r in rows], dtype=bool),
            "fit": np.array([float(r["root_fit_value"] or "nan") for r in rows], dtype=float),
            "ratio": np.array([float(r["root_fit_ratio"] or "nan") for r in rows], dtype=float),
            "ratio_err": np.array([float(r["root_fit_ratio_error"] or "nan") for r in rows], dtype=float),
            "fit_params": rows[0].get("root_fit_params", ""),
        }
    return out


def combine(data: dict[str, dict[str, np.ndarray]], samples: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = data[samples[0]]["x"]
    y = np.zeros_like(x)
    e2 = np.zeros_like(x)
    for sample in samples:
        d = data[sample]
        used = d["used"]
        y[used] += d["y"][used]
        e2[used] += d["ey"][used] ** 2
    return x, y, np.sqrt(e2)


def eval_root_hagedorn(grid: np.ndarray, params: str) -> np.ndarray:
    """Evaluate the ROOT TF1 fitted by `extract_ppg12_rootfit_stitch_points.py`."""

    p = np.array([float(x) for x in params.split(";")], dtype=float)
    return p[0] * np.power(p[1] / grid, p[2] + p[3] * np.log(grid / p[1]) + p[4] * grid)


def draw(
    *,
    data: dict[str, dict[str, np.ndarray]],
    samples: list[str],
    out: Path,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    ratio_ylim: tuple[float, float],
    fit_range: tuple[float, float],
    ylabel: str,
    ratio_ylabel: str,
    xlabel: str,
    legend_anchor: tuple[float, float],
) -> None:
    stitched = data["stitched"]
    x = stitched["x"]
    y = stitched["y"]
    ratio = stitched["ratio"]
    ratio_err = stitched["ratio_err"]
    grid = np.linspace(fit_range[0], xlim[1], 900)
    pred = eval_root_hagedorn(grid, stitched["fit_params"])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
        }
    )

    fig = plt.figure(figsize=(6.4, 7.45), dpi=180)
    ax = fig.add_axes([0.14, 0.35, 0.79, 0.58])
    rax = fig.add_axes([0.14, 0.095, 0.79, 0.25], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    rax.set_ylim(*ratio_ylim)

    ax.plot(grid, pred, color="red", lw=1.7)
    for sample in samples:
        d = data[sample]
        mask = d["used"] & (d["x"] >= xlim[0]) & (d["x"] <= xlim[1]) & (d["y"] > 0)
        ax.errorbar(
            d["x"][mask],
            d["y"][mask],
            yerr=d["ey"][mask],
            fmt="o",
            ms=3.0,
            lw=0.9,
            color=COLORS[sample],
            label=sample,
        )

    rmask = (x >= xlim[0]) & (x <= xlim[1]) & np.isfinite(ratio) & (y > 0)
    rax.axhline(1.0, color="black", lw=0.9, ls=(0, (6, 6)))
    rax.errorbar(x[rmask], ratio[rmask], yerr=ratio_err[rmask], fmt="o", ms=3.0, lw=0.9, color="black")

    ax.set_ylabel(ylabel, fontsize=15)
    rax.set_ylabel(ratio_ylabel, fontsize=15)
    rax.set_xlabel(xlabel, fontsize=16)
    ax.tick_params(labelsize=13, which="both")
    rax.tick_params(labelsize=13, which="both")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.minorticks_on()
    rax.minorticks_on()

    ax.text(
        0.47,
        0.985,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p{+}p$ $\sqrt{s}=200$ GeV" + "\nPYTHIA8",
        transform=ax.transAxes,
        fontsize=13.5,
        va="top",
    )
    ax.legend(
        frameon=False,
        fontsize=12.5,
        loc="upper left",
        bbox_to_anchor=legend_anchor,
        handlelength=1.4,
        handletextpad=0.5,
        borderaxespad=0,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    photon = load_group(ROOTFIT_CSV, "photon")
    draw(
        data=photon,
        samples=["photon5", "photon10", "photon20"],
        out=PHOTON_OUT,
        xlim=(10, 40),
        ylim=(0.03, 3.5e4),
        ratio_ylim=(0.85, 1.15),
        fit_range=(10, 36),
        ylabel=r"$d\sigma / dE_T^\gamma$ [pb / GeV]",
        ratio_ylabel="Data / Fit",
        xlabel=r"Leading   $E_T^\gamma$ [GeV]",
        legend_anchor=(0.60, 0.79),
    )
    print(PHOTON_OUT)

    jet = load_group(ROOTFIT_CSV, "jet")
    draw(
        data=jet,
        samples=["jet8", "jet12", "jet20", "jet30", "jet40"],
        out=JET_OUT,
        xlim=(9, 50),
        ylim=(1e4, 1e12),
        ratio_ylim=(0.85, 1.15),
        fit_range=(10, 50),
        ylabel="counts",
        ratio_ylabel="MC / Fit",
        xlabel=r"Leading   $p_T^\mathrm{jet}$ [GeV]",
        legend_anchor=(0.60, 0.70),
    )
    print(JET_OUT)


if __name__ == "__main__":
    main()
