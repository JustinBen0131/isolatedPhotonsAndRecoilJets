#!/usr/bin/env python3
"""Render a slide-6-style Jet12+20+30+40 embedded-inclusive stitch spectrum."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SAMPLES = ["Jet12", "Jet20", "Jet30", "Jet40"]
COLORS = {
    "Jet12": "#2ca02c",
    "Jet20": "#1296f3",
    "Jet30": "#ff6f00",
    "Jet40": "#7b2cbf",
}
MARKERS = {"Jet12": "o", "Jet20": "s", "Jet30": "^", "Jet40": "v"}
OWNERSHIP = {
    "Jet12": "12-21 GeV",
    "Jet20": "21-31 GeV",
    "Jet30": "31-41 GeV",
    "Jet40": ">=41 GeV",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--x-min", type=float, default=12.0)
    parser.add_argument("--x-max", type=float, default=50.0)
    return parser.parse_args()


def load_rows(path: Path) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, dict[str, float]]]:
    rows: dict[str, list[dict[str, str]]] = {sample: [] for sample in SAMPLES}
    meta: dict[str, dict[str, float]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            sample = row["sample"]
            if sample not in rows:
                continue
            rows[sample].append(row)
            meta[sample] = {
                "sigma_eff_pb": float(row["sigma_eff_pb"]),
                "nraw": float(row["nraw"]),
                "weight_rel": float(row["weight_rel"]),
            }

    data: dict[str, dict[str, np.ndarray]] = {}
    for sample, sample_rows in rows.items():
        if not sample_rows:
            raise RuntimeError(f"missing rows for {sample}")
        sample_rows.sort(key=lambda item: float(item["bin_low"]))
        data[sample] = {
            "x": np.array([float(r["bin_center"]) for r in sample_rows]),
            "xlo": np.array([float(r["bin_low"]) for r in sample_rows]),
            "xhi": np.array([float(r["bin_high"]) for r in sample_rows]),
            "y": np.array([float(r["weighted_entries"]) for r in sample_rows]),
            "ey": np.array([float(r["weighted_error"]) for r in sample_rows]),
            "raw": np.array([float(r["raw_entries"]) for r in sample_rows]),
        }
    return data, meta


def combine(data: dict[str, dict[str, np.ndarray]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = data[SAMPLES[0]]["x"]
    y = np.zeros_like(x)
    e2 = np.zeros_like(x)
    for sample in SAMPLES:
        y += data[sample]["y"]
        e2 += data[sample]["ey"] ** 2
    return x, y, np.sqrt(e2)


def fit_curve(x: np.ndarray, y: np.ndarray, x_min: float, x_max: float) -> tuple[np.ndarray, np.ndarray]:
    mask = (x >= x_min) & (x <= x_max) & (y > 0.0) & np.isfinite(y)
    if np.count_nonzero(mask) < 6:
        raise RuntimeError("not enough positive bins for stitch-spectrum fit")
    coeff = np.polyfit(np.log(x[mask]), np.log(y[mask]), deg=4)
    grid = np.linspace(x_min, x_max, 800)
    pred = np.exp(np.polyval(coeff, np.log(grid)))
    return grid, pred


def add_sphenix_label(ax: plt.Axes) -> None:
    ax.text(
        0.52,
        0.96,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\nAu+Au embedded PYTHIA8" + "\nInclusive Jet12+20+30+40",
        transform=ax.transAxes,
        fontsize=12.5,
        va="top",
    )


def main() -> None:
    args = parse_args()
    data, meta = load_rows(args.csv)
    x, y, ey = combine(data)
    grid, pred = fit_curve(x, y, args.x_min, args.x_max)
    fit_at_x = np.interp(x, grid, pred, left=np.nan, right=np.nan)
    ratio = np.divide(y, fit_at_x, out=np.full_like(y, np.nan), where=(fit_at_x > 0.0) & (y > 0.0))
    ratio_err = np.divide(ey, fit_at_x, out=np.full_like(ey, np.nan), where=(fit_at_x > 0.0) & (y > 0.0))

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.35,
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

    fig = plt.figure(figsize=(7.4, 8.2), dpi=190)
    ax = fig.add_axes([0.13, 0.35, 0.80, 0.58])
    rax = fig.add_axes([0.13, 0.095, 0.80, 0.25], sharex=ax)

    ax.set_yscale("log")
    ax.set_xlim(args.x_min, args.x_max)
    positive = y[(x >= args.x_min) & (x <= args.x_max) & (y > 0.0)]
    ymin = max(np.nanmin(positive) * 0.25, 1.0)
    ymax = np.nanmax(positive) * 35.0
    ax.set_ylim(ymin, ymax)
    rax.set_ylim(0.75, 1.25)

    ax.plot(grid, pred, color="red", lw=2.1, label="log-polynomial fit")
    for boundary in [21.0, 31.0, 41.0]:
        ax.axvline(boundary, color="0.65", lw=1.0, ls=(0, (4, 4)), zorder=0)
        rax.axvline(boundary, color="0.70", lw=1.0, ls=(0, (4, 4)), zorder=0)

    for sample in SAMPLES:
        d = data[sample]
        mask = (d["x"] >= args.x_min) & (d["x"] <= args.x_max) & (d["y"] > 0.0)
        ax.errorbar(
            d["x"][mask],
            d["y"][mask],
            yerr=d["ey"][mask],
            fmt=MARKERS[sample],
            ms=4.8,
            lw=1.0,
            color=COLORS[sample],
            label=sample,
        )

    rmask = (x >= args.x_min) & (x <= args.x_max) & np.isfinite(ratio) & (y > 0.0)
    rax.axhline(1.0, color="black", lw=1.0, ls=(0, (6, 6)))
    rax.errorbar(x[rmask], ratio[rmask], yerr=ratio_err[rmask], fmt="o", ms=4.0, lw=1.0, color="black")

    ax.set_ylabel("weighted counts", fontsize=17, fontweight="bold")
    rax.set_ylabel("MC / Fit", fontsize=17, fontweight="bold")
    rax.set_xlabel(r"max truth-jet filter $p_T$ [GeV]", fontsize=17, fontweight="bold")
    ax.tick_params(labelsize=13.5, which="both")
    rax.tick_params(labelsize=13.5, which="both")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.minorticks_on()
    rax.minorticks_on()
    add_sphenix_label(ax)

    ax.legend(
        frameon=False,
        fontsize=11.0,
        loc="lower left",
        bbox_to_anchor=(0.05, 0.08),
        ncol=2,
        columnspacing=0.9,
        handlelength=1.3,
        handletextpad=0.45,
    )

    lines = ["ownership / weights"]
    for sample in SAMPLES:
        m = meta[sample]
        lines.append(
            f"{sample}: {OWNERSHIP[sample]}\n"
            f"  sigma={m['sigma_eff_pb']:.4g} pb, w={m['weight_rel']:.4g}"
        )
    ax.text(
        0.985,
        0.59,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.3,
        bbox={"facecolor": "white", "edgecolor": "0.82", "boxstyle": "round,pad=0.35", "alpha": 0.92},
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=190)
    plt.close(fig)
    print(args.out)


if __name__ == "__main__":
    main()
