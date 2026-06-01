#!/usr/bin/env python3
"""Build a slide-6-style full-slide candidate with Jet12+20+30+40."""

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
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


SAMPLES = ["Jet12", "Jet20", "Jet30", "Jet40"]
COLORS = {
    "Jet12": "#1f4cff",
    "Jet20": "#ff7f0e",
    "Jet30": "#d62aa0",
    "Jet40": "#7b2cbf",
}
MARKERS = {"Jet12": "o", "Jet20": "s", "Jet30": "^", "Jet40": "v"}
CONFIGS = {
    "Jet12": "phpythia8_10GeV_JS_MDC2.cfg",
    "Jet20": "phpythia8_20GeV_JS_MDC2.cfg",
    "Jet30": "phpythia8_30GeV_JS_MDC2.cfg",
    "Jet40": "phpythia8_40GeV_JS_MDC2.cfg",
}
OWNERSHIP = {
    "Jet12": r"$12 \leq p_T < 21$",
    "Jet20": r"$21 \leq p_T < 31$",
    "Jet30": r"$31 \leq p_T < 41$",
    "Jet40": r"$p_T \geq 41$",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def load_summary(path: Path) -> dict:
    return json.loads(path.read_text())


def load_1gev(csv_path: Path) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, dict[str, float]]]:
    by_sample_bin: dict[tuple[str, int], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    meta: dict[str, dict[str, float]] = {}
    with csv_path.open() as handle:
        for row in csv.DictReader(handle):
            sample = row["sample"]
            if sample not in SAMPLES:
                continue
            lo = float(row["bin_low"])
            hi = float(row["bin_high"])
            cen = float(row["bin_center"])
            if hi <= 12.0 or lo >= 50.0:
                continue
            one_gev_lo = int(np.floor(cen))
            if one_gev_lo < 12 or one_gev_lo >= 50:
                continue
            raw = float(row["raw_entries"])
            raw_err = float(row["raw_error"])
            sigma = float(row["sigma_eff_pb"])
            nraw = float(row["nraw"])
            weight_rel = float(row["weight_rel"])
            key = (sample, one_gev_lo)
            by_sample_bin[key]["raw"] += raw
            by_sample_bin[key]["raw_err2"] += raw_err * raw_err
            by_sample_bin[key]["sigma"] = sigma
            by_sample_bin[key]["nraw"] = nraw
            meta[sample] = {"sigma_eff_pb": sigma, "nraw": nraw, "weight_rel": weight_rel}

    data: dict[str, dict[str, np.ndarray]] = {}
    for sample in SAMPLES:
        xs: list[float] = []
        ys: list[float] = []
        eys: list[float] = []
        raws: list[float] = []
        for lo in range(12, 50):
            item = by_sample_bin[(sample, lo)]
            raw = item.get("raw", 0.0)
            sigma = item.get("sigma", meta[sample]["sigma_eff_pb"])
            nraw = item.get("nraw", meta[sample]["nraw"])
            err = item.get("raw_err2", 0.0) ** 0.5
            xs.append(lo + 0.5)
            ys.append(raw * sigma / nraw if nraw > 0.0 else 0.0)
            eys.append(err * sigma / nraw if nraw > 0.0 else 0.0)
            raws.append(raw)
        data[sample] = {
            "x": np.array(xs),
            "y": np.array(ys),
            "ey": np.array(eys),
            "raw": np.array(raws),
        }
    return data, meta


def combined(data: dict[str, dict[str, np.ndarray]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    x = data[SAMPLES[0]]["x"]
    y = np.zeros_like(x)
    e2 = np.zeros_like(x)
    owner: list[str] = []
    for idx, xc in enumerate(x):
        vals = {sample: data[sample]["y"][idx] for sample in SAMPLES}
        sample = max(vals, key=vals.get)
        owner.append(sample)
        y[idx] = vals[sample]
        e2[idx] = data[sample]["ey"][idx] ** 2
    return x, y, np.sqrt(e2), owner


def fit_spectrum(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = (x >= 12.0) & (x <= 50.0) & (y > 0.0) & np.isfinite(y)
    xfit = x[mask]
    yfit = y[mask]
    popt = np.polyfit(np.log(xfit), np.log(yfit), deg=4)
    grid = np.linspace(12.0, 50.0, 800)
    return popt, grid, np.exp(np.polyval(popt, np.log(grid)))


def add_box(fig: plt.Figure, xywh: tuple[float, float, float, float]) -> plt.Axes:
    ax = fig.add_axes(xywh)
    ax.set_axis_off()
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor="#f8fafc", edgecolor="#d5dce6", linewidth=1.2))
    return ax


def draw_left_boxes(fig: plt.Figure, summary: dict, meta: dict[str, dict[str, float]]) -> None:
    def fig_box(x: float, y: float, w: float, h: float) -> None:
        fig.patches.append(
            Rectangle((x, y), w, h, transform=fig.transFigure, facecolor="#f8fafc", edgecolor="#d5dce6", linewidth=1.2, zorder=0)
        )

    def put(x: float, y: float, text: str, **kwargs) -> None:
        kwargs.setdefault("family", "serif")
        kwargs.setdefault("zorder", 5)
        fig.text(x, y, text, **kwargs)

    fig_box(0.030, 0.455, 0.460, 0.390)
    put(0.042, 0.805, "Generator reproduction + event ownership", fontsize=21, fontweight="bold")
    put(0.042, 0.740, "Base path:", fontsize=13.5)
    put(0.120, 0.740, "/sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/embed_2025", fontsize=10.6, color="0.35")
    put(0.042, 0.690, "Configs:", fontsize=14.5)
    y = 0.657
    for sample in SAMPLES:
        put(0.064, y, f"• {sample} → {CONFIGS[sample]}", fontsize=12.2)
        y -= 0.027
    put(0.042, 0.535, "Jet filter: generator-level inclusive-jet sample", fontsize=12.0)
    put(0.042, 0.505, r"Define: $p_T$(jet, filter) = highest-$p_T$ generator jet passing the filter", fontsize=12.0)
    put(0.042, 0.465, "Ownership:", fontsize=14.0)
    y = 0.442
    for sample in SAMPLES:
        put(0.064, y, f"• {sample}: {OWNERSHIP[sample]}", fontsize=10.9, color=COLORS[sample])
        y -= 0.019

    fig_box(0.030, 0.205, 0.460, 0.220)
    put(0.042, 0.387, "Derived stitched weights", fontsize=21, fontweight="bold")
    headers = ["Sample", "Ownership region", "Nraw", r"$\sigma_{\mathrm{eff}}$", "rel. w"]
    xs = [0.042, 0.126, 0.262, 0.338, 0.430]
    for xh, h in zip(xs, headers):
        put(xh, 0.348, h, fontsize=11.6)
    fig.lines.append(plt.Line2D([0.042, 0.470], [0.333, 0.333], transform=fig.transFigure, color="0.72", lw=1, zorder=5))
    y = 0.303
    for sample in SAMPLES:
        m = meta[sample]
        put(xs[0], y, sample, fontsize=10.8, color=COLORS[sample])
        put(xs[1], y, OWNERSHIP[sample], fontsize=10.5)
        put(xs[2], y, f"{m['nraw']/1e6:.3f}M", fontsize=10.5)
        put(xs[3], y, f"{m['sigma_eff_pb']:.4g} pb", fontsize=10.5)
        put(xs[4], y, f"{m['weight_rel']:.3g}", fontsize=10.5)
        y -= 0.027
    put(0.042, 0.222, r"$w_i = (\sigma_{\mathrm{eff},i}/N_{\mathrm{raw},i}) / (\sigma_{\mathrm{eff},40}/N_{\mathrm{raw},40})$", fontsize=10.0)


def draw_plot(fig: plt.Figure, data: dict[str, dict[str, np.ndarray]]) -> dict[str, float]:
    x, y, ey, owner = combined(data)
    coeff, grid, pred = fit_spectrum(x, y)
    fit_at_x = np.exp(np.polyval(coeff, np.log(x)))
    ratio = y / fit_at_x
    ratio_err = ey / fit_at_x

    ax = fig.add_axes([0.558, 0.405, 0.410, 0.390])
    rax = fig.add_axes([0.558, 0.305, 0.410, 0.090], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(12, 50)
    ax.set_ylim(2e0, 2e6)
    rax.set_ylim(0.88, 1.12)
    for b in (21, 31, 41):
        ax.axvline(b, color="0.55", lw=0.9, ls=":")
        rax.axvline(b, color="0.55", lw=0.9, ls=":")
    ax.plot(grid, pred, color="0.18", lw=1.7, ls=(0, (1.1, 1.1)), label="Fit: modified power law")

    for sample in SAMPLES:
        d = data[sample]
        mask = np.array([o == sample for o in owner])
        ax.errorbar(d["x"][mask], d["y"][mask], yerr=d["ey"][mask], color=COLORS[sample], marker=MARKERS[sample], ms=4.0, lw=0.9, label=f"{sample} stitched")
        rax.errorbar(d["x"][mask], ratio[mask], yerr=ratio_err[mask], color=COLORS[sample], marker=MARKERS[sample], ms=3.1, lw=0.8)

    ax.text(0.985, 1.03, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", fontsize=11.5, family="serif")
    ax.text(0.985, 0.93, r"PYTHIA8 embedded inclusive jets, $\sqrt{s_{NN}}=200$ GeV", transform=ax.transAxes, ha="right", fontsize=10.4, family="serif")
    ax.text(0.04, 0.31, "Embedded inclusive-jet stitch", transform=ax.transAxes, fontsize=11.5, family="serif")
    ax.text(0.04, 0.245, r"Jet12: $12 \leq p_T^{jet,filter}<21$ GeV", transform=ax.transAxes, fontsize=8.6, family="serif")
    ax.text(0.04, 0.195, r"Jet20: $21 \leq p_T^{jet,filter}<31$ GeV", transform=ax.transAxes, fontsize=8.6, family="serif")
    ax.text(0.04, 0.145, r"Jet30: $31 \leq p_T^{jet,filter}<41$ GeV", transform=ax.transAxes, fontsize=8.6, family="serif")
    ax.text(0.04, 0.095, r"Jet40: $p_T^{jet,filter}\geq41$ GeV", transform=ax.transAxes, fontsize=8.6, family="serif")

    handles, labels = ax.get_legend_handles_labels()
    order = [1, 2, 3, 4, 0]
    ax.legend([handles[i] for i in order], [labels[i] for i in order], frameon=False, fontsize=8.5, loc="upper right", bbox_to_anchor=(0.82, 0.76))
    ax.set_ylabel(r"$\sigma_{\mathrm{eff}}\times N$ scaled entries [pb / bin]", fontsize=13.5, family="serif")
    rax.set_ylabel("stitched / fit", fontsize=8.5, family="serif")
    rax.set_xlabel(r"$p_T^{jet,filter}$ [GeV]", fontsize=11.5, family="serif", loc="right")
    ax.tick_params(labelbottom=False, labelsize=9.5, which="both", direction="in", top=True, right=True)
    rax.tick_params(labelsize=8.0, which="both", direction="in", top=True, right=True)
    rax.axhline(1.0, color="0.45", lw=0.9, ls=":")
    rax.set_yticks([0.9, 1.0, 1.1])
    ax.minorticks_on()
    rax.minorticks_on()

    jumps: dict[str, float] = {}
    for boundary in (21, 31, 41):
        before_idx = np.where(np.isclose(x, boundary - 0.5))[0][0]
        after_idx = np.where(np.isclose(x, boundary + 0.5))[0][0]
        jumps[str(boundary)] = float(ratio[after_idx] / ratio[before_idx])
    return jumps


def draw_bottom(fig: plt.Figure, jumps: dict[str, float]) -> None:
    ax = add_box(fig, [0.030, 0.055, 0.940, 0.120])
    ax.text(
        0.5,
        0.62,
        r"The bounded 21/31/41 GeV ownership rule removes double counting across the Jet12 $\rightarrow$ Jet20 $\rightarrow$ Jet30 $\rightarrow$ Jet40 stitched spectrum.",
        ha="center",
        va="center",
        fontsize=18.5,
        family="serif",
    )
    ax.text(
        0.5,
        0.24,
        f"Boundary checks: 21 GeV jump = {jumps['21']:.3f}; 31 GeV jump = {jumps['31']:.3f}; 41 GeV jump = {jumps['41']:.3f}.",
        ha="center",
        va="center",
        fontsize=13.5,
        color="0.45",
        family="serif",
    )


def main() -> None:
    args = parse_args()
    summary = load_summary(args.summary)
    data, meta = load_1gev(args.csv)
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
    })
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    fig.text(0.025, 0.925, "Embedded Inclusive Jet 12+20+30+40 Stitching", fontsize=35, fontweight="bold", family="serif")
    draw_left_boxes(fig, summary, meta)
    jumps = draw_plot(fig, data)
    draw_bottom(fig, jumps)
    fig.text(0.985, 0.025, "6", ha="right", va="bottom", fontsize=22, fontweight="bold", family="serif")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=160)
    plt.close(fig)
    print(args.out)


if __name__ == "__main__":
    main()
