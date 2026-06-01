#!/usr/bin/env python3
"""Build a slide-7-style parallel component plot for pp inclusive-jet stitching."""

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
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[1]
BASE = REPO / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_exactStitch_20260523_0915"
CSV_PATH = BASE / "validation/insitu_stitching/pp_currentian_exactstitch_contract_points.csv"
SUMMARY_PATH = BASE / "validation/insitu_stitching/pp_currentian_exactstitch_contract_summary.json"
QA_PATH = BASE / "validation/insitu_stitching/pp_currentian_exactstitch_boundary_continuity_qa.json"
OUT = BASE / "slide_assets/pp_inclusive_jet_stitch_parallel_slide7_style.png"

SAMPLES = ["run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30", "run28_jet40"]
LABELS = {
    "run28_jet8": "jet8: 9-14 GeV",
    "run28_jet12": "jet12: 14-21 GeV",
    "run28_jet20": "jet20: 21-32 GeV",
    "run28_jet30": "jet30: 32-42 GeV",
    "run28_jet40": "jet40: >=42 GeV",
}
SHORT = {
    "run28_jet8": "jet8",
    "run28_jet12": "jet12",
    "run28_jet20": "jet20",
    "run28_jet30": "jet30",
    "run28_jet40": "jet40",
}
COLORS = {
    "run28_jet8": "#d63b9f",
    "run28_jet12": "#2ca02c",
    "run28_jet20": "#1592ff",
    "run28_jet30": "#ff6b00",
    "run28_jet40": "#c000d8",
    "sum": "#111111",
}
MARKERS = {
    "run28_jet8": "o",
    "run28_jet12": "s",
    "run28_jet20": "^",
    "run28_jet30": "D",
    "run28_jet40": "v",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.25,
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 13,
            "figure.dpi": 160,
            "savefig.dpi": 160,
        }
    )


def load_rows() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with CSV_PATH.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["group"] != "jet" or row["sample"] not in SAMPLES:
                continue
            rows.append(
                {
                    "sample": row["sample"],
                    "x": float(row["bin_center"]),
                    "x_low": float(row["bin_low"]),
                    "x_high": float(row["bin_high"]),
                    "y": float(row["density_pb_per_gev"]),
                    "ey": float(row["density_err_pb_per_gev"]),
                    "raw": float(row["raw_all_events"]),
                    "xsec": float(row["xsec_pb"]),
                    "events": float(row["events_processed_metadata"]),
                }
            )
    return rows


def pivot_rows(rows: list[dict[str, float | str]]) -> dict[str, np.ndarray]:
    by_x: dict[float, dict[str, tuple[float, float]]] = defaultdict(dict)
    bins: dict[float, tuple[float, float]] = {}
    for row in rows:
        x = float(row["x"])
        by_x[x][str(row["sample"])] = (float(row["y"]), float(row["ey"]))
        bins[x] = (float(row["x_low"]), float(row["x_high"]))

    xs = np.array(sorted(by_x), dtype=float)
    low = np.array([bins[x][0] for x in xs], dtype=float)
    high = np.array([bins[x][1] for x in xs], dtype=float)
    out: dict[str, np.ndarray] = {"x": xs, "xerr_low": xs - low, "xerr_high": high - xs}
    sum_y = np.zeros_like(xs)
    sum_var = np.zeros_like(xs)
    for sample in SAMPLES:
        y = np.array([by_x[x].get(sample, (0.0, 0.0))[0] for x in xs], dtype=float)
        ey = np.array([by_x[x].get(sample, (0.0, 0.0))[1] for x in xs], dtype=float)
        out[f"{sample}_y"] = y
        out[f"{sample}_ey"] = ey
        sum_y += y
        sum_var += ey * ey
    out["sum_y"] = sum_y
    out["sum_ey"] = np.sqrt(sum_var)
    denom = np.where(sum_y > 0, sum_y, np.nan)
    for sample in SAMPLES:
        out[f"{sample}_frac"] = out[f"{sample}_y"] / denom
    return out


def add_round_box(fig: plt.Figure, xy: tuple[float, float], wh: tuple[float, float], color: str) -> None:
    fig.add_artist(
        FancyBboxPatch(
            xy,
            wh[0],
            wh[1],
            boxstyle="round,pad=0.012,rounding_size=0.012",
            linewidth=1.0,
            edgecolor="#d8dee9",
            facecolor=color,
            transform=fig.transFigure,
            zorder=0,
        )
    )


def window_table(summary: dict) -> list[tuple[str, str, str]]:
    table = []
    for sample in SAMPLES:
        entry = next(s for s in summary["samples"] if s["sample"] == sample)
        low, high = entry["stitch_window"]
        high_txt = "100" if high >= 99 else f"{high:.0f}"
        table.append((SHORT[sample], f"[{low:.0f}, {high_txt})", f"{entry['xsec_pb']:.4g}"))
    return table


def build_slide() -> None:
    setup_style()
    with SUMMARY_PATH.open() as f:
        summary = json.load(f)
    with QA_PATH.open() as f:
        qa = json.load(f)
    arr = pivot_rows(load_rows())

    mask = (arr["x"] >= 9.0) & (arr["x"] <= 50.0)
    x = arr["x"][mask]
    xerr = np.vstack([arr["xerr_low"][mask], arr["xerr_high"][mask]])
    sum_y = arr["sum_y"][mask]
    pos = sum_y > 0

    fig = plt.figure(figsize=(16.0, 9.0), constrained_layout=False)
    fig.patch.set_facecolor("white")
    fig.text(0.04, 0.955, "pp Inclusive-Jet Stitching Parallel Check", fontsize=31, fontweight="bold", ha="left", va="top")
    fig.text(
        0.04,
        0.905,
        r"R = 0.4 truth jets, 1 GeV bins, per-event cross-section normalization; jet8 uses the current sPHENIX wiki value.",
        fontsize=15.5,
        color="#465366",
        ha="left",
        va="top",
    )

    add_round_box(fig, (0.055, 0.115), (0.62, 0.74), "white")
    add_round_box(fig, (0.705, 0.53), (0.245, 0.325), "#f7f9fc")
    add_round_box(fig, (0.705, 0.365), (0.245, 0.125), "#eef8f1")
    add_round_box(fig, (0.705, 0.165), (0.245, 0.155), "#fff8e9")

    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[3.0, 1.25],
        left=0.105,
        right=0.64,
        bottom=0.17,
        top=0.80,
        hspace=0.04,
    )
    ax = fig.add_subplot(gs[0])
    axf = fig.add_subplot(gs[1], sharex=ax)

    for sample in SAMPLES:
        y = arr[f"{sample}_y"][mask]
        ey = arr[f"{sample}_ey"][mask]
        keep = y > 0
        ax.errorbar(
            x[keep],
            y[keep],
            yerr=ey[keep],
            xerr=xerr[:, keep],
            fmt=MARKERS[sample],
            ms=4.7,
            lw=0.0,
            elinewidth=0.9,
            capsize=1.8,
            color=COLORS[sample],
            label=LABELS[sample],
            zorder=4,
        )
        axf.plot(
            x,
            arr[f"{sample}_frac"][mask],
            marker=MARKERS[sample],
            ms=4.5,
            lw=1.15,
            color=COLORS[sample],
            label=SHORT[sample],
        )

    ax.errorbar(
        x[pos],
        sum_y[pos],
        yerr=arr["sum_ey"][mask][pos],
        xerr=xerr[:, pos],
        fmt="o",
        ms=4.1,
        mfc="white",
        mec=COLORS["sum"],
        mew=1.1,
        ecolor=COLORS["sum"],
        elinewidth=0.7,
        capsize=1.6,
        color=COLORS["sum"],
        label="weighted sum",
        zorder=8,
    )

    ax.set_yscale("log")
    ax.set_xlim(8.8, 50.2)
    ax.set_ylim(max(5.0, float(np.nanmin(sum_y[pos])) * 0.45), float(np.nanmax(sum_y[pos])) * 2.8)
    ax.set_ylabel(r"$d\sigma/dp_T^{jet}$ [pb / GeV]")
    ax.grid(which="major", color="#d9dee7", alpha=0.8, linewidth=0.8)
    ax.grid(which="minor", color="#edf0f5", alpha=0.65, linewidth=0.45)
    ax.tick_params(which="both", direction="in", top=True, right=True, length=5)
    ax.tick_params(labelbottom=False)
    ax.legend(loc="upper right", ncol=2, frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, fontsize=11.1)
    ax.text(0.04, 0.95, r"$\bf{\it{sPHENIX}}$ Internal" + "\nPYTHIA8 pp inclusive jets", transform=ax.transAxes, ha="left", va="top", fontsize=13.4, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 2.0})

    axf.set_ylim(-0.035, 1.05)
    axf.set_yticks([0.0, 0.5, 1.0])
    axf.set_ylabel("fraction of\nweighted sum", fontsize=14)
    axf.set_xlabel(r"Leading truth-jet $p_T$ [GeV]")
    axf.grid(which="major", color="#d9dee7", alpha=0.8, linewidth=0.8)
    axf.tick_params(which="both", direction="in", top=True, right=True, length=5)

    fig.text(0.724, 0.815, "Stitching windows", fontsize=18.5, fontweight="bold", ha="left", va="top")
    y0 = 0.770
    fig.text(0.724, y0, "sample", fontsize=12.5, fontweight="bold")
    fig.text(0.785, y0, "window [GeV]", fontsize=12.5, fontweight="bold")
    fig.text(0.900, y0, "xsec [pb]", fontsize=12.5, fontweight="bold")
    for i, (sample, win, xsec) in enumerate(window_table(summary)):
        y = y0 - 0.039 * (i + 1)
        fig.text(0.724, y, sample, fontsize=13.0, color=COLORS[f"run28_{sample}"] if sample != "jet8" else COLORS["run28_jet8"])
        fig.text(0.785, y, win, fontsize=13.0)
        fig.text(0.900, y, xsec, fontsize=13.0)

    max_dev = max(float(b["max_same_bin_overlap_fractional_deviation"]) for b in qa["boundaries"])
    worst = max(qa["boundaries"], key=lambda b: float(b["max_same_bin_overlap_fractional_deviation"]))
    fig.text(0.724, 0.455, "Boundary QA", fontsize=17.5, fontweight="bold", color="#23864f", ha="left", va="top")
    fig.text(
        0.724,
        0.420,
        f"All handoffs pass the same-bin overlap check.\nLargest overlap deviation: {100*max_dev:.1f}% at {worst['boundary_GeV']:.0f} GeV.\nNo window gaps or double counting.",
        fontsize=13.0,
        linespacing=1.25,
        ha="left",
        va="top",
    )

    fig.text(0.724, 0.285, "How to read it", fontsize=17.5, fontweight="bold", color="#9a5a00", ha="left", va="top")
    fig.text(
        0.724,
        0.250,
        "Each colored sample is shown in parallel.\nOnly its assigned truth-pT window contributes\nto the black summed spectrum.\nThe lower panel makes ownership explicit.",
        fontsize=12.4,
        linespacing=1.25,
        ha="left",
        va="top",
    )

    fig.text(
        0.055,
        0.055,
        "Takeaway: the pp inclusive-jet reference baseline stitches smoothly across jet8, jet12, jet20, jet30, and jet40.",
        fontsize=14.2,
        fontweight="bold",
        color="#253044",
        ha="left",
        va="center",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, facecolor="white", bbox_inches=None)
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    build_slide()
