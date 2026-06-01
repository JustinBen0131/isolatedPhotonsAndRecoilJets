#!/usr/bin/env python3
"""Render slide-8/9 stitching candidates from Justin/RecoilJets contract CSV.

The input CSV must be produced by compare_pp_currentian_insitu_stitch_contract.py
from RecoilJets in-situ ROOT outputs.  This script deliberately refuses the
older PPG12/Shuhang-source CSV names so a reference plot cannot be mislabeled as
"This Analysis Output".
"""

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
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


LABELS = {
    "run28_photonjet5": "photon5",
    "run28_photonjet10": "photon10",
    "run28_photonjet20": "photon20",
    "run28_jet8": "jet8",
    "run28_jet12": "jet12",
    "run28_jet20": "jet20",
    "run28_jet30": "jet30",
    "run28_jet40": "jet40",
}

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


def load_group(csv_path: Path, group: str) -> dict[str, dict[str, np.ndarray]]:
    if "ppg12_currentian_truth_spectrum_root_histograms" in csv_path.name:
        raise RuntimeError("refusing PPG12/Shuhang-source CSV; use in-situ contract CSV")
    grouped: dict[str, list[dict[str, str]]] = {}
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            if row["group"] != group:
                continue
            grouped.setdefault(row["sample"], []).append(row)
    out: dict[str, dict[str, np.ndarray]] = {}
    for sample, rows in grouped.items():
        out[sample] = {
            "x": np.array([float(r["bin_center"]) for r in rows], dtype=float),
            "y": np.array([float(r["density_pb_per_gev"]) for r in rows], dtype=float),
            "ey": np.array([float(r["density_err_pb_per_gev"]) for r in rows], dtype=float),
            "count": np.array([float(r["raw_kept_events"]) for r in rows], dtype=float),
        }
    return out


def combine(data: dict[str, dict[str, np.ndarray]], samples: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = data[samples[0]]["x"]
    y = np.zeros_like(x)
    e2 = np.zeros_like(x)
    for sample in samples:
        d = data[sample]
        mask = d["count"] > 0
        y[mask] += d["y"][mask]
        e2[mask] += d["ey"][mask] ** 2
    return x, y, np.sqrt(e2)


def ppg12_modified_power_law(
    x: np.ndarray,
    p0: float,
    p1: float,
    p2: float,
    p3: float,
    p4: float,
) -> np.ndarray:
    return p0 * np.power(p1 / x, p2 + p3 * np.log(x / p1) + p4 * x)


def fit_curve(
    x: np.ndarray,
    y: np.ndarray,
    ey: np.ndarray,
    xmin: float,
    xmax: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = (x >= xmin) & (x <= xmax) & np.isfinite(y) & np.isfinite(ey) & (y > 0) & (ey > 0)
    if np.count_nonzero(mask) < 8:
        raise RuntimeError("not enough populated points to fit stitching spectrum")
    # Match the PPG12 ROOT macro fit:
    # [0]*pow([1]/x,[2]+[3]*log(x/[1])+[4]*x), seeded identically.
    p0 = np.array([2.09375e9, 1.0, 1.0, 2.0, 0.01], dtype=float)
    params, _ = curve_fit(
        ppg12_modified_power_law,
        x[mask],
        y[mask],
        p0=p0,
        sigma=ey[mask],
        absolute_sigma=False,
        method="lm",
        maxfev=200000,
    )
    grid = np.linspace(xmin, xmax, 800)
    pred = ppg12_modified_power_law(grid, *params)
    if not np.all(np.isfinite(pred)) or np.nanmax(pred) <= 0:
        raise RuntimeError("PPG12 modified-power-law fit produced non-finite prediction")
    return grid, pred, params


def render(csv_path: Path, out_dir: Path, group: str) -> Path:
    data = load_group(csv_path, group)
    if group == "photon":
        samples = ["run28_photonjet5", "run28_photonjet10", "run28_photonjet20"]
        out = out_dir / "pp_currentIAN_photon_truth_stitch_slide8_side_by_side_insitu_contract.png"
        xlim = (10, 40)
        ylim = (0.03, 4e4)
        fit_range = (10, 40)
        ylabel = r"$d\sigma / dE_T^\gamma$ [pb / GeV]"
        ratio_ylabel = "Data / Fit"
        xlabel = r"Leading   $E_T^\gamma$ [GeV]"
    else:
        samples = ["run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30", "run28_jet40"]
        out = out_dir / "pp_currentIAN_inclusive_jet_truth_stitch_slide9_side_by_side_insitu_contract.png"
        xlim = (9, 50)
        ylim = (1e4, 1e12)
        fit_range = (10, 50)
        ylabel = "counts"
        ratio_ylabel = "MC / Fit"
        xlabel = r"Leading   $p_T^\mathrm{jet}$ [GeV]"

    missing = [s for s in samples if s not in data]
    if missing:
        raise RuntimeError(f"missing samples in {csv_path}: {missing}")
    x, y, ey = combine(data, samples)
    grid, pred, fit_params = fit_curve(x, y, ey, *fit_range)
    fit_y = np.interp(x, grid, pred, left=np.nan, right=np.nan)
    ratio = np.divide(y, fit_y, out=np.full_like(y, np.nan), where=(fit_y > 0) & (y > 0))
    ratio_err = np.divide(ey, fit_y, out=np.full_like(ey, np.nan), where=(fit_y > 0) & (y > 0))

    visual_scale = 1.0
    if group == "jet":
        first = y[(x >= 9) & (x <= 12) & (y > 0)]
        visual_scale = 1.0e12 / float(np.nanmax(first)) if len(first) else 1.0

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.3,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(6.4, 7.45), dpi=180)
    ax = fig.add_axes([0.14, 0.35, 0.79, 0.58])
    rax = fig.add_axes([0.14, 0.095, 0.79, 0.25], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    rax.set_ylim(0.85, 1.15)

    ax.plot(grid, pred * visual_scale, color="red", lw=1.8)
    for sample in samples:
        d = data[sample]
        label = LABELS[sample]
        mask = (d["count"] > 0) & (d["x"] >= xlim[0]) & (d["x"] <= xlim[1])
        ax.errorbar(d["x"][mask], d["y"][mask] * visual_scale, yerr=d["ey"][mask] * visual_scale,
                    fmt="o", ms=3.2, lw=0.9, color=COLORS[label], label=label)

    rmask = (x >= xlim[0]) & (x <= xlim[1]) & np.isfinite(ratio) & (y > 0)
    rax.axhline(1.0, color="black", lw=0.9, ls=(0, (6, 6)))
    rax.errorbar(x[rmask], ratio[rmask], yerr=ratio_err[rmask], fmt="o", ms=3.0, lw=0.9, color="black")

    ax.set_ylabel(ylabel, fontsize=15)
    rax.set_ylabel(ratio_ylabel, fontsize=15)
    rax.set_xlabel(xlabel, fontsize=16)
    ax.tick_params(labelsize=13, which="both")
    rax.tick_params(labelsize=13, which="both")
    ax.minorticks_on()
    rax.minorticks_on()
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.text(0.47, 0.985, r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p{+}p$ $\sqrt{s}=200$ GeV" + "\nPYTHIA8",
            transform=ax.transAxes, fontsize=13.5, va="top")
    ax.text(0.025, 0.035, "Fit: PPG12 modified power law", transform=ax.transAxes,
            fontsize=8.7, va="bottom", ha="left")
    ax.legend(frameon=False, fontsize=12.5, loc="upper left", bbox_to_anchor=(0.58, 0.78),
              ncol=1 if group == "photon" else 2, columnspacing=0.9, handlelength=1.4)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    fit_out = out.with_suffix(".fit.txt")
    fit_out.write_text(
        "formula=[0]*pow([1]/x,[2]+[3]*log(x/[1])+[4]*x)\n"
        + "fit_range={:.6g},{:.6g}\n".format(*fit_range)
        + "params=" + ",".join(f"{p:.12g}" for p in fit_params) + "\n"
    )
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--groups", default="photon,jet", help="Comma-separated groups to render, e.g. photon or photon,jet")
    args = parser.parse_args()
    groups = [x.strip() for x in args.groups.split(",") if x.strip()]
    for group in groups:
        if group not in {"photon", "jet"}:
            raise RuntimeError(f"unknown group: {group}")
        print(render(args.csv, args.out_dir, group))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
