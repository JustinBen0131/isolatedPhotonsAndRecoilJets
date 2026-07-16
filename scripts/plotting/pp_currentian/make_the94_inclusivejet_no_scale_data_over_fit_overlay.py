#!/usr/bin/env python3
"""Render THE-94 Fig.6 inclusive-jet no-scale overlay with MC/Fit subpanel.

This is a presentation-only renderer for the current THE-94 diagnostic state:
the top panel compares the true period-combined PPG12 Fig.6 source to the fixed
THE-94 output with no global scale, while the bottom panel overlays both
sources divided by the PPG12 Fig.6 fit guide.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
THE94_BASE = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
)
NO_SCALE_DIR = THE94_BASE / "no_scale_overlay"
IN_RATIO_CSV = NO_SCALE_DIR / "inclusivejet_fig6_true_period_combined_no_scale_ratio.csv"
OUT_PNG = NO_SCALE_DIR / "inclusivejet_fig6_true_period_combined_no_scale_overlay.png"
OUT_POINTS = NO_SCALE_DIR / "inclusivejet_fig6_true_period_combined_no_scale_data_over_fit_points.csv"
OUT_MANIFEST = NO_SCALE_DIR / "inclusivejet_fig6_true_period_combined_no_scale_overlay_manifest.json"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]
COLORS = {
    "jet8": "#d62aa0",
    "jet12": "#239b35",
    "jet20": "#169ce8",
    "jet30": "#ff6f00",
    "jet40": "#d62aa0",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def ratio_error(value: float, err: float, fit: float) -> float:
    if value <= 0 or fit <= 0 or not all(math.isfinite(v) for v in (value, err, fit)):
        return float("nan")
    return err / fit


def fit_ppg12_reference(rows: list[dict[str, str]]) -> tuple[dict[tuple[str, float], float], np.ndarray, np.ndarray]:
    """Build the smooth Fig.6 guide from the plotted PPG12 period-combined points."""
    x = np.array([f(row, "bin_center") for row in rows], dtype=float)
    y = np.array([f(row, "ppg12_true_period_value") for row in rows], dtype=float)
    mask = (x >= 9.0) & (x <= 50.0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    coeff = np.polyfit(np.log(x[mask]), np.log(y[mask]), deg=4)
    grid = np.linspace(9.0, 50.0, 600)
    curve = np.exp(np.polyval(coeff, np.log(grid)))
    fit_by_bin: dict[tuple[str, float], float] = {}
    for row in rows:
        sample = row["sample"]
        center = round(f(row, "bin_center"), 6)
        fit_by_bin[(sample, center)] = float(np.exp(np.polyval(coeff, np.log(center))))
    return fit_by_bin, grid, curve


def build_points() -> tuple[list[dict[str, object]], np.ndarray, np.ndarray]:
    source_rows = read_csv(IN_RATIO_CSV)
    fit, fit_grid, fit_curve = fit_ppg12_reference(source_rows)
    rows: list[dict[str, object]] = []
    for row in source_rows:
        sample = row["sample"]
        center = round(f(row, "bin_center"), 6)
        fit_value = fit.get((sample, center), float("nan"))
        ppg12 = f(row, "ppg12_true_period_value")
        ppg12_err = f(row, "ppg12_true_period_error")
        current = f(row, "current_kept_value")
        current_err = f(row, "current_kept_error")
        rows.append(
            {
                "sample": sample,
                "bin_low": f(row, "bin_low"),
                "bin_high": f(row, "bin_high"),
                "bin_center": center,
                "ppg12_true_period_value": ppg12,
                "ppg12_true_period_error": ppg12_err,
                "current_kept_value": current,
                "current_kept_error": current_err,
                "ppg12_fit_value": fit_value,
                "ppg12_true_period_over_fit": ppg12 / fit_value if fit_value > 0 else float("nan"),
                "ppg12_true_period_over_fit_error": ratio_error(ppg12, ppg12_err, fit_value),
                "current_kept_over_fit": current / fit_value if fit_value > 0 else float("nan"),
                "current_kept_over_fit_error": ratio_error(current, current_err, fit_value),
                "current_over_ppg12_true_period": f(row, "current_over_ppg12_true_period"),
            }
        )
    rows.sort(key=lambda r: (SAMPLES.index(str(r["sample"])), float(r["bin_center"])))
    return rows, fit_grid, fit_curve


def write_points(rows: list[dict[str, object]]) -> None:
    fields = [
        "sample",
        "bin_low",
        "bin_high",
        "bin_center",
        "ppg12_true_period_value",
        "ppg12_true_period_error",
        "current_kept_value",
        "current_kept_error",
        "ppg12_fit_value",
        "ppg12_true_period_over_fit",
        "ppg12_true_period_over_fit_error",
        "current_kept_over_fit",
        "current_kept_over_fit_error",
        "current_over_ppg12_true_period",
    ]
    with OUT_POINTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def stats(rows: list[dict[str, object]]) -> dict[str, object]:
    by_sample: dict[str, dict[str, object]] = {}
    for sample in SAMPLES:
        sample_rows = [r for r in rows if r["sample"] == sample]
        ratios = [float(r["current_over_ppg12_true_period"]) for r in sample_rows]
        if ratios:
            by_sample[sample] = {
                "n_bins": len(ratios),
                "median_current_over_ppg12": float(np.median(ratios)),
                "min_current_over_ppg12": min(ratios),
                "max_current_over_ppg12": max(ratios),
            }
    return by_sample


def render(rows: list[dict[str, object]], fit_grid: np.ndarray, fit_curve: np.ndarray) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.45,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7.5,
            "ytick.major.size": 7.5,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )

    fig = plt.figure(figsize=(8.0, 8.89), dpi=200)
    top = fig.add_axes([0.13, 0.40, 0.79, 0.56])
    bot = fig.add_axes([0.13, 0.09, 0.79, 0.30], sharex=top)
    top.set_yscale("log")
    top.set_xlim(8, 50)
    top.set_ylim(1.0e6, 1.2e13)
    bot.set_ylim(0.60, 2.90)

    for sample in SAMPLES:
        pts = [r for r in rows if r["sample"] == sample]
        if not pts:
            continue
        color = COLORS[sample]
        x = [float(r["bin_center"]) for r in pts]
        top.errorbar(
            x,
            [float(r["ppg12_true_period_value"]) for r in pts],
            yerr=[float(r["ppg12_true_period_error"]) for r in pts],
            fmt="o",
            ms=5.2,
            mfc="white",
            mec=color,
            mew=1.25,
            ecolor=color,
            elinewidth=0.65,
            linestyle="none",
            zorder=3,
        )
        top.errorbar(
            x,
            [float(r["current_kept_value"]) for r in pts],
            yerr=[float(r["current_kept_error"]) for r in pts],
            fmt="o",
            ms=4.8,
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
            [float(r["ppg12_true_period_over_fit"]) for r in pts],
            yerr=[float(r["ppg12_true_period_over_fit_error"]) for r in pts],
            fmt="o",
            ms=4.6,
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
            [float(r["current_kept_over_fit"]) for r in pts],
            yerr=[float(r["current_kept_over_fit_error"]) for r in pts],
            fmt="o",
            ms=4.1,
            mfc="black",
            mec="black",
            mew=0.6,
            ecolor="black",
            elinewidth=0.48,
            linestyle="none",
            zorder=4,
        )

    top.plot(fit_grid, fit_curve, color="red", lw=1.6, zorder=2)
    bot.axhline(1.0, color="0.45", lw=1.0, ls=(0, (5, 5)))

    top.set_ylabel("counts", fontsize=24)
    bot.set_ylabel("MC / Fit", fontsize=22)
    bot.set_xlabel(r"Leading $p_T^\mathrm{jet}$ [GeV]", fontsize=25)
    top.tick_params(which="both", labelbottom=False, labelsize=18)
    bot.tick_params(which="both", labelsize=18)
    top.minorticks_on()
    bot.minorticks_on()

    top.text(
        0.56,
        0.97,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p$+$p$ $\sqrt{s}=200$ GeV" + "\nPYTHIA8",
        transform=top.transAxes,
        fontsize=21,
        va="top",
    )
    top.text(
        0.06,
        0.18,
        "no scale applied\ntrue period-combined PPG12 reference",
        transform=top.transAxes,
        fontsize=18,
        va="bottom",
    )

    sample_handles = [
        Line2D([0], [0], marker="o", color=COLORS[s], mfc=COLORS[s], lw=0, ms=8, label=s)
        for s in SAMPLES
    ]
    source_handles = [
        Line2D([0], [0], marker="o", color="black", mfc="white", lw=0, ms=8, label="PPG12 period-combined"),
        Line2D([0], [0], marker="o", color="black", mfc="black", lw=0, ms=8, label="THE-94 fixed output"),
        Line2D([0], [0], color="red", lw=1.8, label="PPG12 fit guide"),
    ]
    leg1 = top.legend(
        handles=sample_handles,
        title="sample",
        frameon=False,
        fontsize=15,
        title_fontsize=15,
        loc="lower left",
        bbox_to_anchor=(0.07, 0.31),
        handletextpad=0.6,
        labelspacing=0.52,
    )
    top.add_artist(leg1)
    top.legend(
        handles=source_handles,
        title="source",
        frameon=False,
        fontsize=13.5,
        title_fontsize=13.5,
        loc="upper left",
        bbox_to_anchor=(0.53, 0.61),
        handletextpad=0.7,
        labelspacing=0.55,
    )

    fig.savefig(OUT_PNG)
    plt.close(fig)


def write_manifest(rows: list[dict[str, object]]) -> None:
    manifest = {
        "current_roots": str(THE94_BASE / "final_roots/inclusivejet"),
        "decision": "not_parity_complete_jet8_fails_no_scale",
        "input_ratio_csv": str(IN_RATIO_CSV),
        "no_global_scale_applied": True,
        "png": str(OUT_PNG),
        "points_csv": str(OUT_POINTS),
        "ppg12_reference": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom.root::h_max_truth_jet_pT, extracted read-only and rebinned to 1 GeV",
        "fit_source": "degree-4 log-polynomial guide computed from the plotted true period-combined PPG12 Fig.6 points",
        "bottom_panel": "PPG12 true period-combined / same-reference PPG12 fit guide and THE-94 fixed output / same-reference PPG12 fit guide",
        "summary": stats(rows),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    rows, fit_grid, fit_curve = build_points()
    write_points(rows)
    render(rows, fit_grid, fit_curve)
    write_manifest(rows)
    print(OUT_PNG)
    print(OUT_POINTS)
    print(OUT_MANIFEST)


if __name__ == "__main__":
    main()
