#!/usr/bin/env python3
"""Make slide-5/6-style pp current-IAN truth-stitching spectrum PNGs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_fullsim_20260521_1811"
BASE = REPO / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN
IN_CSV = BASE / "validation/currentIAN_stitching/ppg12_currentian_truth_spectrum_root_histograms.csv"
IN_JSON = BASE / "validation/currentIAN_stitching/ppg12_currentian_truth_spectrum_root_histograms_summary.json"
OUT_DIR = BASE / "slide_assets"


def load_rows(group: str) -> dict[str, dict[str, np.ndarray | str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    with IN_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["group"] != group:
                continue
            grouped.setdefault(row["sample"], []).append(row)
    out: dict[str, dict[str, np.ndarray | str]] = {}
    for sample, rows in grouped.items():
        out[sample] = {
            "x": np.array([float(r["bin_center"]) for r in rows], dtype=float),
            "y": np.array([float(r["density_pb_per_gev"]) for r in rows], dtype=float),
            "ey": np.array([float(r["density_err_pb_per_gev"]) for r in rows], dtype=float),
            "used": np.array([int(r["used_in_stitch"]) for r in rows], dtype=bool),
            "color": rows[0]["color"],
            "lo": float(rows[0]["stitch_window_low"]),
            "hi": float(rows[0]["stitch_window_high"]),
            "xsec": float(rows[0]["xsec_pb"]),
        }
    return out


def combined_xy(data: dict[str, dict[str, np.ndarray | str]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs = None
    ysum = None
    esum2 = None
    for d in data.values():
        x = d["x"]  # type: ignore[assignment]
        y = d["y"]  # type: ignore[assignment]
        ey = d["ey"]  # type: ignore[assignment]
        used = d["used"]  # type: ignore[assignment]
        if xs is None:
            xs = np.array(x, copy=True)
            ysum = np.zeros_like(xs)
            esum2 = np.zeros_like(xs)
        ysum[used] += y[used]  # type: ignore[index]
        esum2[used] += ey[used] ** 2  # type: ignore[index]
    assert xs is not None and ysum is not None and esum2 is not None
    return xs, ysum, np.sqrt(esum2)


def fit_curve(x: np.ndarray, y: np.ndarray, xmin: float, xmax: float) -> tuple[np.ndarray, np.ndarray]:
    mask = (x >= xmin) & (x <= xmax) & np.isfinite(y) & (y > 0)
    xx = x[mask]
    yy = y[mask]
    # Stable log-polynomial form of a modified power-law fit.
    coeff = np.polyfit(np.log(xx), np.log(yy), deg=4)
    grid = np.linspace(xmin, xmax, 600)
    pred = np.exp(np.polyval(coeff, np.log(grid)))
    return grid, pred


def eval_fit(x: np.ndarray, grid: np.ndarray, pred: np.ndarray) -> np.ndarray:
    return np.interp(x, grid, pred, left=np.nan, right=np.nan)


def fmt_window(lo: float, hi: float) -> str:
    if hi >= 100:
        return rf"$\geq {lo:g}$"
    return rf"${lo:g}$--${hi:g}$"


def draw_slide(group: str) -> Path:
    meta = json.loads(IN_JSON.read_text())
    data = load_rows(group)
    x, y, ey = combined_xy(data)
    if group == "photon":
        samples = ["photon5", "photon10", "photon20"]
        out = OUT_DIR / "pp_currentIAN_photon_truth_stitch_slide5_style.png"
        title = r"Photon+jet pp samples stitch smoothly in leading truth photon $E_T$"
        subtitle = r"Current PPG12 IAN windows: photon5 0--14 GeV, photon10 14--22 GeV, photon20 $\geq$22 GeV"
        xlim = (10, 40)
        fit_range = (10, 40)
        xlabel = r"Leading truth photon $E_T^\gamma$ [GeV]"
        ylabel = r"$d\sigma/dE_T^\gamma$ [pb / GeV]"
        fig_label = "Photon+jet signal-MC stitching"
        caption = (
            r"Leading truth photon $E_T^\gamma$ distributions from PYTHIA-8 photon 5, 10, and 20 GeV samples. "
            r"Each sample is shown only in its current-IAN ownership window; the combined spectrum is fit with a smooth modified power-law form."
        )
    else:
        samples = ["jet8", "jet12", "jet20", "jet30", "jet40"]
        out = OUT_DIR / "pp_currentIAN_inclusive_jet_truth_stitch_slide6_style.png"
        title = r"Inclusive pp jet samples stitch smoothly in leading truth jet $p_T$"
        subtitle = r"Current PPG12 IAN windows: jet8 9--14, jet12 14--21, jet20 21--32, jet30 32--42, jet40 $\geq$42 GeV"
        xlim = (9, 50)
        fit_range = (10, 50)
        xlabel = r"Leading truth jet $p_T^{jet}$ [GeV]"
        ylabel = r"$d\sigma/dp_T^{jet}$ [pb / GeV]"
        fig_label = "Inclusive-jet background-MC stitching"
        caption = (
            r"Leading truth jet $p_T$ distributions from PYTHIA-8 jet 8, 12, 20, 30, and 40 GeV samples. "
            r"The stitched spectrum uses non-overlapping current-IAN truth-jet ownership windows and the recorded PPG12 truth-spectrum weights."
        )
    grid, pred = fit_curve(x, y, *fit_range)
    fit_y = eval_fit(x, grid, pred)
    ratio = np.divide(y, fit_y, out=np.full_like(y, np.nan), where=(fit_y > 0) & (y > 0))
    ratio_err = np.divide(ey, fit_y, out=np.full_like(ey, np.nan), where=(fit_y > 0) & (y > 0))

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.3,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(0.055, 0.942, title, fontsize=30, fontweight="bold", color="#0f172a", va="top")
    fig.text(0.055, 0.895, subtitle, fontsize=16.5, color="#475569", va="top")

    left, width = 0.23, 0.45
    ax = fig.add_axes([left, 0.38, width, 0.42])
    rax = fig.add_axes([left, 0.22, width, 0.15], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    positive = y[(x >= xlim[0]) & (x <= xlim[1]) & (y > 0)]
    ax.set_ylim(max(positive.min() * 0.35, 1e-4), positive.max() * 5.0)
    rax.set_ylim(0.85, 1.15)

    for sample in samples:
        d = data[sample]
        sx = d["x"]  # type: ignore[assignment]
        sy = d["y"]  # type: ignore[assignment]
        sey = d["ey"]  # type: ignore[assignment]
        used = d["used"]  # type: ignore[assignment]
        color = str(d["color"])
        mask = used & (sx >= xlim[0]) & (sx <= xlim[1]) & (sy > 0)  # type: ignore[operator]
        ax.errorbar(sx[mask], sy[mask], yerr=sey[mask], fmt="o", ms=3.1, lw=0.8, color=color, label=sample)  # type: ignore[index]

    ax.plot(grid, pred, color="#e11d48", lw=1.6, label="smooth fit")
    rmask = (x >= xlim[0]) & (x <= xlim[1]) & np.isfinite(ratio) & (y > 0)
    rax.axhline(1.0, color="#6b7280", lw=1.0, ls=(0, (4, 4)))
    rax.errorbar(x[rmask], ratio[rmask], yerr=ratio_err[rmask], fmt="o", ms=2.8, lw=0.8, color="black")

    ax.grid(True, which="major", alpha=0.16)
    rax.grid(True, which="major", alpha=0.16)
    ax.set_ylabel(ylabel, fontsize=16)
    rax.set_ylabel("MC / fit", fontsize=15)
    rax.set_xlabel(xlabel, fontsize=17)
    ax.tick_params(labelsize=13)
    rax.tick_params(labelsize=13)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.text(
        0.48,
        0.96,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p{+}p$ $\sqrt{s}=200$ GeV" + "\nPYTHIA8",
        transform=ax.transAxes,
        fontsize=13,
        va="top",
    )
    ax.legend(frameon=False, fontsize=12, loc="upper right")

    # Reference labels and exact weight/window callout sit outside the plot, like a slide caption.
    fig.text(0.055, 0.815, fig_label, fontsize=18, fontweight="bold", color="#0f172a")
    window_lines = []
    for sample in samples:
        d = data[sample]
        window_lines.append(
            f"{sample}: {fmt_window(float(d['lo']), float(d['hi']))} GeV, "
            + rf"$\sigma$={float(d['xsec']):.4g} pb"
        )
    fig.text(
        0.70,
        0.78,
        "Weights and ownership",
        fontsize=17,
        fontweight="bold",
        color="#0f172a",
    )
    fig.text(
        0.70,
        0.745,
        "Histogram entries are weighted by the PPG12\ntruth-spectrum per-entry weight; plotted\nas pb/GeV using the ROOT metadata.",
        fontsize=13.5,
        color="#334155",
        linespacing=1.25,
        va="top",
    )
    fig.text(0.70, 0.63, "\n".join(window_lines), fontsize=12.5, color="#334155", linespacing=1.22, va="top")

    max_dev = np.nanmax(np.abs(ratio[rmask] - 1.0)) * 100.0
    fig.text(
        0.70,
        0.315,
        "Validation readout",
        fontsize=17,
        fontweight="bold",
        color="#0f172a",
    )
    fig.text(
        0.70,
        0.282,
        f"Max plotted |MC/fit - 1| = {max_dev:.1f}%.\nThe key check is smooth behavior across\nsample handoff boundaries, not the absolute\nnormalization of an individual slice.",
        fontsize=13.5,
        color="#334155",
        linespacing=1.25,
        va="top",
    )
    fig.text(0.055, 0.095, caption, fontsize=15.5, color="#111827", ha="left", va="bottom", wrap=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def main() -> None:
    for group in ("photon", "jet"):
        print(draw_slide(group))


if __name__ == "__main__":
    main()
