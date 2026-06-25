#!/usr/bin/env python3
"""Build THE-69 ABCD signal-leakage overlay slide from current AuAu signal MC."""

from __future__ import annotations

import csv
import json
import math
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


INPUT_ROOT = (
    REPO
    / "InputFiles/the69_default_auau_physicsqa/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
OUT_DIR = (
    REPO
    / "dataOutput/auauPhysicsQA/THE69_defaultAuAuBDT_physicsQA_20260620/abcd_leakage_overlay_slide"
)
OUT_PNG = OUT_DIR / "the69_abcd_signal_leakage_overlay_slide.png"
OUT_CSV = OUT_DIR / "the69_abcd_signal_leakage_points.csv"
OUT_JSON = OUT_DIR / "the69_abcd_signal_leakage_overlay_manifest.json"
OUT_SCRIPT = OUT_DIR / "the69_abcd_signal_leakage_overlay_speaker_script.md"

TOPDIR = "SIM"
ISO_TAG = "isoR40_isSliding"
PT_BINS = [(14, 16), (16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]
CENTRALITIES = [
    ("0_20", "0-20%", "#111827", "o", -0.10),
    ("20_50", "20-50%", "#2563EB", "s", 0.00),
    ("50_80", "50-80%", "#059669", "D", 0.10),
]
LEAKAGE_REGIONS = [
    ("B", "tight, non-isolated", "B leakage: tight but non-isolated", "#B45309"),
    ("C", "isolated, non-tight", "C leakage: isolated but non-tight", "#7C3AED"),
    ("D", "non-isolated, non-tight", "D leakage: non-isolated and non-tight", "#0F766E"),
]


@dataclass(frozen=True)
class Count:
    value: float
    error: float


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#1F2937",
            "axes.linewidth": 1.0,
            "axes.labelsize": 12.5,
            "axes.titlesize": 15.5,
            "xtick.labelsize": 10.5,
            "ytick.labelsize": 10.5,
            "legend.fontsize": 10.5,
            "mathtext.fontset": "dejavuserif",
        }
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return root_file


def hist_path(lo: int, hi: int, cent_key: str) -> str:
    return f"{TOPDIR}/h_sigABCD_MC_{ISO_TAG}_pT_{lo}_{hi}_cent_{cent_key}"


def bin_count(hist, bin_index: int) -> Count:
    value = float(hist.GetBinContent(bin_index))
    error = float(hist.GetBinError(bin_index))
    if error <= 0.0 and value > 0.0:
        error = math.sqrt(value)
    return Count(value, error)


def ratio_with_error(num: Count, den: Count) -> tuple[float, float]:
    if den.value <= 0.0 or not np.isfinite(den.value):
        return float("nan"), float("nan")
    value = num.value / den.value
    if num.value > 0.0:
        rel2 = (num.error / num.value) ** 2 + (den.error / den.value) ** 2
        error = abs(value) * math.sqrt(max(0.0, rel2))
    else:
        error = num.error / den.value
    return value, error


def collect_points() -> list[dict]:
    root_file = open_root(INPUT_ROOT)
    rows: list[dict] = []
    try:
        for cent_key, cent_label, _, _, _ in CENTRALITIES:
            for lo, hi in PT_BINS:
                path = hist_path(lo, hi, cent_key)
                hist = root_file.Get(path)
                if not hist:
                    raise KeyError(f"Missing histogram: {path}")
                counts = {
                    "A": bin_count(hist, 1),
                    "B": bin_count(hist, 2),
                    "C": bin_count(hist, 3),
                    "D": bin_count(hist, 4),
                }
                for region, _, _, _ in LEAKAGE_REGIONS:
                    leakage, leakage_err = ratio_with_error(counts[region], counts["A"])
                    rows.append(
                        {
                            "region": region,
                            "centrality": cent_label,
                            "cent_key": cent_key,
                            "pt_lo": lo,
                            "pt_hi": hi,
                            "pt_mid": 0.5 * (lo + hi),
                            "a_signal": counts["A"].value,
                            "a_signal_err": counts["A"].error,
                            "sideband_signal": counts[region].value,
                            "sideband_signal_err": counts[region].error,
                            "leakage_fraction": leakage,
                            "leakage_fraction_err": leakage_err,
                            "histogram": path,
                        }
                    )
    finally:
        root_file.Close()
    return rows


def write_points(rows: list[dict]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "region",
        "centrality",
        "cent_key",
        "pt_lo",
        "pt_hi",
        "pt_mid",
        "a_signal",
        "a_signal_err",
        "sideband_signal",
        "sideband_signal_err",
        "leakage_fraction",
        "leakage_fraction_err",
        "histogram",
    ]
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def add_definition_table(fig) -> None:
    ax = fig.add_axes([0.045, 0.705, 0.91, 0.182])
    ax.set_axis_off()
    box = FancyBboxPatch(
        (0.0, 0.0),
        1.0,
        1.0,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        facecolor="#F8FAFC",
        edgecolor="#CBD5E1",
        linewidth=1.2,
        transform=ax.transAxes,
    )
    ax.add_patch(box)
    ax.text(0.025, 0.83, "Leakage definition", fontsize=15.0, fontweight="bold", color="#111827", va="center")
    ax.text(
        0.025,
        0.64,
        r"$A_{\rm sig}$ is the truth-matched prompt-signal count reconstructed in the selected region "
        r"(tight + isolated).  Each sideband leakage is $f_X = N_{\rm sig}(X) / N_{\rm sig}(A)$.",
        fontsize=11.6,
        color="#374151",
        va="center",
    )

    x0 = 0.025
    y0 = 0.13
    col_w = [0.12, 0.32, 0.22, 0.22]
    headers = ["Region", "Reco class", "Leakage ratio", "ABCD use"]
    cells = [
        ["B", "tight, non-isolated", r"$f_B = B_{\rm sig}/A_{\rm sig}$", "signal contamination in B"],
        ["C", "isolated, non-tight", r"$f_C = C_{\rm sig}/A_{\rm sig}$", "signal contamination in C"],
        ["D", "non-isolated, non-tight", r"$f_D = D_{\rm sig}/A_{\rm sig}$", "signal contamination in D"],
    ]
    row_h = 0.115
    for i, header in enumerate(headers):
        ax.text(
            x0 + sum(col_w[:i]),
            y0 + 3.18 * row_h,
            header,
            fontsize=10.6,
            fontweight="bold",
            color="#111827",
        )
    for row_idx, row in enumerate(cells):
        y = y0 + (2 - row_idx) * row_h
        for col_idx, cell in enumerate(row):
            color = next((c for r, _, _, c in LEAKAGE_REGIONS if r == row[0]), "#111827") if col_idx == 0 else "#334155"
            weight = "bold" if col_idx == 0 else "normal"
            ax.text(x0 + sum(col_w[:col_idx]), y, cell, fontsize=10.2, color=color, fontweight=weight)

    ax.text(0.755, 0.83, "Centrality overlay", fontsize=11.0, fontweight="bold", color="#111827", va="center")
    for idx, (_, cent_label, color, marker, _) in enumerate(CENTRALITIES):
        y = 0.65 - idx * 0.18
        ax.plot([0.77], [y], marker=marker, markersize=7.2, color=color, markeredgecolor="white", markeredgewidth=0.7)
        ax.text(0.79, y, cent_label, fontsize=10.8, color="#111827", va="center")


def plot_panel(ax, rows: list[dict], region: str, title: str, accent: str) -> None:
    region_rows = [row for row in rows if row["region"] == region]
    for cent_key, cent_label, color, marker, offset in CENTRALITIES:
        subset = [row for row in region_rows if row["cent_key"] == cent_key]
        x = np.array([row["pt_mid"] + offset for row in subset], dtype=float)
        y = np.array([row["leakage_fraction"] for row in subset], dtype=float)
        yerr = np.array([row["leakage_fraction_err"] for row in subset], dtype=float)
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=marker,
            linestyle="none",
            markersize=6.0,
            capsize=2.6,
            elinewidth=1.2,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.6,
            color=color,
            label=cent_label,
            zorder=3,
        )
    ax.set_title(title, color=accent, fontweight="bold", pad=8)
    ax.set_xlabel(r"cluster $E_T$ bin center [GeV]")
    ax.set_ylabel(r"Leakage ratio")
    ax.grid(True, which="major", color="#E5E7EB", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.set_xlim(13.4, 31.1)
    ax.set_xticks([15, 17, 19, 21, 23, 25, 30.5])
    ax.set_xticklabels(["14-16", "16-18", "18-20", "20-22", "22-24", "24-26", "26-35"])
    yvals = np.array([row["leakage_fraction"] for row in region_rows], dtype=float)
    max_y = float(np.nanmax(yvals)) if yvals.size else 1.0
    ax.set_ylim(0.0, max(0.08, max_y * 1.17))
    ax.axhline(1.0, color="#94A3B8", linewidth=1.0, linestyle="--", zorder=1)
    if max_y > 1.0:
        ax.text(
            0.025,
            0.92,
            "ratio can exceed 1",
            transform=ax.transAxes,
            fontsize=9.5,
            color="#475569",
            ha="left",
            va="top",
        )


def render_slide(rows: list[dict]) -> None:
    setup_style()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.text(
        0.045,
        0.953,
        "Prompt-signal ABCD leakage from the current Au+Au baseline",
        ha="left",
        va="top",
        fontsize=26,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.045,
        0.912,
        "Current THE-69 signal MC, default Au+Au BDT WP80; centralities overlaid in each sideband leakage ratio.",
        ha="left",
        va="top",
        fontsize=13.8,
        color="#334155",
    )
    fig.text(
        0.955,
        0.915,
        "Photon12+20 signal MC\nsliding iso R=0.4",
        ha="right",
        va="top",
        fontsize=12.4,
        color="#475569",
    )
    add_definition_table(fig)

    axes = [
        fig.add_axes([0.060, 0.185, 0.275, 0.415]),
        fig.add_axes([0.372, 0.185, 0.275, 0.415]),
        fig.add_axes([0.684, 0.185, 0.275, 0.415]),
    ]
    for ax, (region, _, title, accent) in zip(axes, LEAKAGE_REGIONS):
        plot_panel(ax, rows, region, title, accent)

    takeaway = (
        "Read this as a signal-contamination correction template for ABCD sidebands, not as a purity. "
        "B leakage is the largest broad component; C and D are dominated by the lowest ET bin, then fall to smaller "
        "centrality-ordered values across the working range."
    )
    fig.text(
        0.060,
        0.073,
        textwrap.fill(takeaway, 135),
        ha="left",
        va="center",
        fontsize=13.2,
        color="#111827",
        bbox=dict(facecolor="#F8FAFC", edgecolor="#CBD5E1", linewidth=1.05, boxstyle="round,pad=0.5"),
    )

    fig.savefig(OUT_PNG, dpi=SLIDE_DPI, bbox_inches=None)
    plt.close(fig)


def write_manifest(rows: list[dict]) -> None:
    max_by_region = {}
    for region, _, _, _ in LEAKAGE_REGIONS:
        vals = [row["leakage_fraction"] for row in rows if row["region"] == region]
        max_by_region[region] = float(np.nanmax(vals)) if vals else None
    manifest = {
        "slide": str(OUT_PNG),
        "points_csv": str(OUT_CSV),
        "speaker_script": str(OUT_SCRIPT),
        "generator": str(THIS_FILE),
        "input_root": str(INPUT_ROOT),
        "topdir": TOPDIR,
        "histogram_template": f"{TOPDIR}/h_sigABCD_MC_{ISO_TAG}_pT_<lo>_<hi>_cent_<cent>",
        "definition": {
            "A": "truth-matched prompt signal reconstructed as tight + isolated",
            "B": "truth-matched prompt signal reconstructed as tight + non-isolated",
            "C": "truth-matched prompt signal reconstructed as isolated + non-tight",
            "D": "truth-matched prompt signal reconstructed as non-isolated + non-tight",
            "leakage": "f_X = N_sig(X) / N_sig(A), X in {B,C,D}",
        },
        "centralities": [cent[1] for cent in CENTRALITIES],
        "pt_bins": [list(p) for p in PT_BINS],
        "max_leakage_by_region": max_by_region,
        "notes": [
            "Ratios are not bounded probabilities; f_X can exceed 1 when the sideband true-signal count exceeds A-region true-signal count.",
            "This is candidate-level sigABCD leakage, not xJ-leading-photon leakage.",
            "Plots use markers only by slide policy; no connecting lines or fits are applied.",
        ],
    }
    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n")


def write_speaker_script() -> None:
    script = """# THE-69 ABCD Signal Leakage Overlay

This slide defines the prompt-signal leakage factors used to correct ABCD sidebands.
The denominator is the truth-matched prompt signal reconstructed in the selected A region: tight and isolated.
B leakage is tight but non-isolated signal relative to A. C leakage is isolated but non-tight signal relative to A. D leakage is non-isolated and non-tight signal relative to A.

The important readout is that these are sideband contamination ratios, not probabilities. They can exceed one in the lowest ET bin if more prompt-signal candidates land in a sideband than in selected A.
The current THE-69 signal MC shows B leakage as the dominant broad leakage component, while C and D are largest in the lowest ET bin and then settle lower across centrality.

Use this as the audience-facing definition and first-look shape check before applying or debating detailed ABCD sideband leakage corrections.
"""
    OUT_SCRIPT.write_text(script)


def main() -> None:
    rows = collect_points()
    write_points(rows)
    render_slide(rows)
    write_manifest(rows)
    write_speaker_script()
    print(json.dumps({"png": str(OUT_PNG), "csv": str(OUT_CSV), "manifest": str(OUT_JSON)}, indent=2))


if __name__ == "__main__":
    main()
