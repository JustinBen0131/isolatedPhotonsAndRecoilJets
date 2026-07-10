#!/usr/bin/env python3
"""Render the THE-94 inclusive-jet Fig.6 jet12-40 parity overlay.

This figure intentionally excludes jet8.  It compares the corrected THE-94
inclusive-jet output to the true period-combined PPG12 Fig.6 reference for the
jet12/20/30/40 slices only, with no global scale applied.  The red curve is a
smooth visual guide fitted to the PPG12 jet12-40 reference points used in this
figure, not a hidden normalization of the THE-94 points.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
THE94_BASE = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
)
INPUT_CSV = (
    THE94_BASE
    / "no_scale_overlay/inclusivejet_fig6_true_period_combined_no_scale_data_over_fit_points.csv"
)
OUT_DIR = THE94_BASE / "jet12to40_overlay"
OUT_PNG = OUT_DIR / "inclusivejet_fig6_jet12to40_ppg12_vs_the94_overlay.png"
OUT_POINTS = OUT_DIR / "inclusivejet_fig6_jet12to40_ppg12_vs_the94_overlay_points.csv"
OUT_MANIFEST = OUT_DIR / "inclusivejet_fig6_jet12to40_ppg12_vs_the94_overlay_manifest.json"

SAMPLES = ["jet12", "jet20", "jet30", "jet40"]
COLORS = {
    "jet12": "#239b35",
    "jet20": "#169ce8",
    "jet30": "#ff6f00",
    "jet40": "#c12ac7",
}


def read_input() -> list[dict[str, str]]:
    with INPUT_CSV.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["sample"] in SAMPLES]
    rows.sort(key=lambda row: (SAMPLES.index(row["sample"]), float(row["bin_center"])))
    return rows


def ff(row: dict[str, str], key: str) -> float:
    return float(row[key])


def fit_log_poly(rows: list[dict[str, str]], degree: int = 4) -> np.ndarray:
    x = np.array([ff(row, "bin_center") for row in rows], dtype=float)
    y = np.array([ff(row, "ppg12_true_period_value") for row in rows], dtype=float)
    mask = (x > 0.0) & (y > 0.0) & np.isfinite(x) & np.isfinite(y)
    return np.polyfit(np.log(x[mask]), np.log(y[mask]), degree)


def eval_fit(coeff: np.ndarray, x: np.ndarray | float) -> np.ndarray | float:
    return np.exp(np.polyval(coeff, np.log(x)))


def ratio_error(error: float, fit_value: float) -> float:
    if fit_value <= 0.0 or not math.isfinite(fit_value):
        return float("nan")
    return error / fit_value


def build_points(rows: list[dict[str, str]], coeff: np.ndarray) -> list[dict[str, object]]:
    points: list[dict[str, object]] = []
    for row in rows:
        x = ff(row, "bin_center")
        fit_value = float(eval_fit(coeff, x))
        ppg12 = ff(row, "ppg12_true_period_value")
        ppg12_err = ff(row, "ppg12_true_period_error")
        current = ff(row, "current_kept_value")
        current_err = ff(row, "current_kept_error")
        points.append(
            {
                "sample": row["sample"],
                "bin_low": ff(row, "bin_low"),
                "bin_high": ff(row, "bin_high"),
                "bin_center": x,
                "ppg12_true_period_value": ppg12,
                "ppg12_true_period_error": ppg12_err,
                "the94_current_value": current,
                "the94_current_error": current_err,
                "jet12to40_ppg12_fit_value": fit_value,
                "ppg12_over_fit": ppg12 / fit_value,
                "ppg12_over_fit_error": ratio_error(ppg12_err, fit_value),
                "the94_over_fit": current / fit_value,
                "the94_over_fit_error": ratio_error(current_err, fit_value),
                "the94_over_ppg12": ff(row, "current_over_ppg12_true_period"),
            }
        )
    return points


def write_points(points: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_POINTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(points[0].keys()))
        writer.writeheader()
        writer.writerows(points)


def sample_summary(points: list[dict[str, object]]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for sample in SAMPLES:
        ratios = [float(row["the94_over_ppg12"]) for row in points if row["sample"] == sample]
        fit_ratios = [float(row["the94_over_fit"]) for row in points if row["sample"] == sample]
        out[sample] = {
            "n_bins": len(ratios),
            "median_the94_over_ppg12": float(median(ratios)),
            "min_the94_over_ppg12": float(min(ratios)),
            "max_the94_over_ppg12": float(max(ratios)),
            "median_the94_over_fit": float(median(fit_ratios)),
        }
    return out


def render(points: list[dict[str, object]], coeff: np.ndarray) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.35,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 3.8,
            "ytick.minor.size": 3.8,
        }
    )

    fig = plt.figure(figsize=(8.1, 8.7), dpi=200)
    top = fig.add_axes([0.13, 0.395, 0.80, 0.565])
    bot = fig.add_axes([0.13, 0.095, 0.80, 0.29], sharex=top)

    top.set_yscale("log")
    top.set_xlim(13.5, 50.5)
    top.set_ylim(6.0e5, 4.0e11)
    bot.set_ylim(0.965, 1.035)

    fit_x = np.linspace(14.0, 50.0, 700)
    fit_y = eval_fit(coeff, fit_x)
    top.plot(fit_x, fit_y, color="red", lw=1.7, zorder=2)
    bot.axhline(1.0, color="0.45", lw=1.1, ls=(0, (5, 5)), zorder=1)

    for sample in SAMPLES:
        rows = [row for row in points if row["sample"] == sample]
        color = COLORS[sample]
        x = [float(row["bin_center"]) for row in rows]

        top.errorbar(
            x,
            [float(row["ppg12_true_period_value"]) for row in rows],
            yerr=[float(row["ppg12_true_period_error"]) for row in rows],
            fmt="o",
            ms=5.0,
            mfc="white",
            mec=color,
            mew=1.3,
            ecolor=color,
            elinewidth=0.65,
            linestyle="none",
            zorder=3,
        )
        top.errorbar(
            x,
            [float(row["the94_current_value"]) for row in rows],
            yerr=[float(row["the94_current_error"]) for row in rows],
            fmt="o",
            ms=4.7,
            mfc=color,
            mec=color,
            mew=0.7,
            ecolor=color,
            elinewidth=0.55,
            linestyle="none",
            zorder=4,
        )

        bot.errorbar(
            x,
            [float(row["ppg12_over_fit"]) for row in rows],
            yerr=[float(row["ppg12_over_fit_error"]) for row in rows],
            fmt="o",
            ms=4.5,
            mfc="white",
            mec="black",
            mew=1.05,
            ecolor="black",
            elinewidth=0.55,
            linestyle="none",
            zorder=3,
        )
        bot.errorbar(
            x,
            [float(row["the94_over_fit"]) for row in rows],
            yerr=[float(row["the94_over_fit_error"]) for row in rows],
            fmt="o",
            ms=4.1,
            mfc="black",
            mec="black",
            mew=0.65,
            ecolor="black",
            elinewidth=0.48,
            linestyle="none",
            zorder=4,
        )

    top.set_ylabel("counts", fontsize=25)
    bot.set_ylabel("Data / Fit", fontsize=22)
    bot.set_xlabel(r"Leading $p_T^\mathrm{jet}$ [GeV]", fontsize=26)
    top.tick_params(axis="x", labelbottom=False)
    for ax in (top, bot):
        ax.minorticks_on()
        ax.tick_params(labelsize=16)

    top.text(
        0.61,
        0.91,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=top.transAxes,
        ha="left",
        va="top",
        fontsize=20,
    )
    top.text(
        0.61,
        0.84,
        r"$p$+$p$ $\sqrt{s}$ = 200 GeV",
        transform=top.transAxes,
        ha="left",
        va="top",
        fontsize=18,
    )
    top.text(
        0.61,
        0.78,
        "PYTHIA8",
        transform=top.transAxes,
        ha="left",
        va="top",
        fontsize=18,
    )
    top.text(
        0.08,
        0.16,
        "inclusive+jet stitch comparison\njet12-40 only, no global scale",
        transform=top.transAxes,
        ha="left",
        va="bottom",
        fontsize=14,
    )

    sample_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=COLORS[sample],
            markeredgecolor=COLORS[sample],
            markersize=7,
            label=sample,
        )
        for sample in SAMPLES
    ]
    source_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="black",
            markerfacecolor="white",
            markeredgecolor="black",
            lw=0,
            markersize=7,
            label="PPG12",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="black",
            markerfacecolor="black",
            markeredgecolor="black",
            lw=0,
            markersize=7,
            label="THE-94 output",
        ),
        Line2D([0], [0], color="red", lw=1.8, label="PPG12 jet12-40 fit"),
    ]

    leg1 = top.legend(
        handles=sample_handles,
        title="sample",
        frameon=False,
        fontsize=13,
        title_fontsize=13,
        loc="lower left",
        bbox_to_anchor=(0.06, 0.33),
        handletextpad=0.6,
        borderaxespad=0,
    )
    top.add_artist(leg1)
    top.legend(
        handles=source_handles,
        title="source",
        frameon=False,
        fontsize=13,
        title_fontsize=13,
        loc="lower left",
        bbox_to_anchor=(0.30, 0.30),
        handlelength=2.0,
        handletextpad=0.7,
        borderaxespad=0,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)


def write_manifest(points: list[dict[str, object]], coeff: np.ndarray) -> None:
    manifest = {
        "artifact": "THE-94 inclusive-jet Fig.6 jet12-40 focused overlay",
        "interpretation": "nominal no-scale parity for jet12/20/30/40 only; jet8 intentionally excluded",
        "input_csv": str(INPUT_CSV),
        "output_png": str(OUT_PNG),
        "output_points_csv": str(OUT_POINTS),
        "samples": SAMPLES,
        "jet8_policy": "excluded from this figure; not rescaled or hidden",
        "fit_guide": {
            "type": "degree-4 log-log polynomial fitted to PPG12 true period-combined jet12-40 points",
            "coefficients_high_to_low": [float(v) for v in coeff],
            "note": "visual data/fit guide only; no scale is applied to THE-94 points",
        },
        "ratio_summary": sample_summary(points),
    }
    with OUT_MANIFEST.open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    rows = read_input()
    coeff = fit_log_poly(rows)
    points = build_points(rows, coeff)
    write_points(points)
    render(points, coeff)
    write_manifest(points, coeff)
    print(OUT_PNG)
    for sample, values in sample_summary(points).items():
        print(
            sample,
            "n=",
            values["n_bins"],
            "median THE94/PPG12=",
            f'{values["median_the94_over_ppg12"]:.9f}',
        )


if __name__ == "__main__":
    main()
