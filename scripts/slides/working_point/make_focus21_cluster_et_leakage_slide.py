#!/usr/bin/env python3
"""Build a full-slide PNG for the inclusive-jet reco-cluster ET leakage check."""

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

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[1]
OUTDIR = REPO / "dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_jet1234_20260526"
COARSE_CSV_PATH = OUTDIR / "inclusive_jet123_reco_cluster_et_weighted_components_cent_summed.csv"
COARSE_SUMMARY_PATH = OUTDIR / "inclusive_jet123_reco_cluster_et_weighted_components_summary.json"
FINE_CSV_PATH = OUTDIR / "inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50.csv"
FINE_SUMMARY_PATH = OUTDIR / "inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50_summary.json"
PNG_PATH = OUTDIR / "inclusive_jet1234_reco_cluster_et_leakage_slide7_followup_12to50.png"
COARSE_BACKUP_PATH = OUTDIR / "inclusive_jet1234_reco_cluster_et_leakage_slide7_followup_coarse_savedbins_backup.png"

SAMPLES = ["Jet12", "Jet20", "Jet30", "Jet40"]
COLORS = {
    "Jet12": "#1f5eff",
    "Jet20": "#ff6b00",
    "Jet30": "#159447",
    "Jet40": "#8f3ffc",
    "Sum": "#111111",
}
MARKERS = {"Jet12": "o", "Jet20": "s", "Jet30": "^", "Jet40": "v", "Sum": "D"}
LABELS = {
    "Jet12": r"Jet12: $12 \leq p_{T}^{truth\,jet} < 21$",
    "Jet20": r"Jet20: $21 \leq p_{T}^{truth\,jet} < 31$",
    "Jet30": r"Jet30: $31 \leq p_{T}^{truth\,jet} < 41$",
    "Jet40": r"Jet40: $p_{T}^{truth\,jet} \geq 41$",
    "Sum": "weighted sum",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.15,
            "axes.labelsize": 14,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 11,
            "figure.dpi": 160,
            "savefig.dpi": 220,
        }
    )


def weighted_pivot(df: pd.DataFrame) -> pd.DataFrame:
    pivot = (
        df.pivot_table(
            index=["x_low", "x_high", "x_center", "x_err_low", "x_err_high"],
            columns="sample",
            values=["weighted_entries", "weighted_error"],
            aggfunc="sum",
        )
        .sort_index()
        .reset_index()
    )
    pivot.columns = [
        "_".join(str(x) for x in col if str(x))
        if isinstance(col, tuple)
        else str(col)
        for col in pivot.columns
    ]
    for sample in SAMPLES:
        for base in ["weighted_entries", "weighted_error"]:
            col = f"{base}_{sample}"
            if col not in pivot:
                pivot[col] = 0.0
    pivot["weighted_entries_Sum"] = sum(
        pivot[f"weighted_entries_{sample}"] for sample in SAMPLES
    )
    pivot["weighted_error_Sum"] = np.sqrt(
        sum(pivot[f"weighted_error_{sample}"] ** 2 for sample in SAMPLES)
    )
    denom = pivot["weighted_entries_Sum"].replace(0, np.nan)
    for sample in SAMPLES:
        pivot[f"fraction_{sample}"] = pivot[f"weighted_entries_{sample}"] / denom
    return pivot


def add_round_box(fig: plt.Figure, xy: tuple[float, float], wh: tuple[float, float], color: str) -> None:
    fig.add_artist(
        FancyBboxPatch(
            xy,
            wh[0],
            wh[1],
            boxstyle="round,pad=0.012,rounding_size=0.018",
            linewidth=0,
            facecolor=color,
            transform=fig.transFigure,
            zorder=0,
        )
    )


def add_in_axes_internal_label(ax: plt.Axes) -> None:
    ax.text(
        0.985,
        0.965,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=13.2,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.8},
        zorder=30,
    )
    ax.text(
        0.985,
        0.910,
        "PYTHIA8 embedded inclusive jet, 0-80%",
        transform=ax.transAxes,
        fontsize=10.8,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
        zorder=30,
    )


def pct(value: float) -> str:
    return f"{100.0 * value:.0f}%" if value == 0 else f"{100.0 * value:.1f}%"


def input_paths() -> tuple[Path, Path, bool]:
    if FINE_CSV_PATH.exists() and FINE_SUMMARY_PATH.exists():
        return FINE_CSV_PATH, FINE_SUMMARY_PATH, True
    return COARSE_CSV_PATH, COARSE_SUMMARY_PATH, False


def range_fraction(summary: dict, label: str, sample: str) -> float:
    for entry in summary.get("bins_of_interest", []):
        if entry.get("x_range") == label:
            return float(entry["fractions"].get(sample, 0.0))
    return float("nan")


def backup_existing_png() -> None:
    if PNG_PATH.exists() and not COARSE_BACKUP_PATH.exists() and not FINE_CSV_PATH.exists():
        COARSE_BACKUP_PATH.write_bytes(PNG_PATH.read_bytes())


def build_slide() -> None:
    setup_style()
    csv_path, summary_path, is_fine = input_paths()
    backup_existing_png()
    df = pd.read_csv(csv_path)
    pivot = weighted_pivot(df)
    with summary_path.open() as f:
        summary = json.load(f)

    jet12_2224 = pct(range_fraction(summary, "22-24", "Jet12"))
    jet12_2426 = pct(range_fraction(summary, "24-26", "Jet12"))
    jet12_2635 = pct(range_fraction(summary, "26-35", "Jet12"))
    jet20_2224 = pct(range_fraction(summary, "22-24", "Jet20"))
    jet20_2426 = pct(range_fraction(summary, "24-26", "Jet20"))
    jet20_2635 = pct(range_fraction(summary, "26-35", "Jet20"))
    is_fine = is_fine and pivot["x_high"].max() >= 35.0
    pivot_plot = pivot[pivot["x_high"] <= 50.0].copy()

    x = pivot_plot["x_center"].to_numpy(float)
    xerr = np.vstack(
        [pivot_plot["x_err_low"].to_numpy(float), pivot_plot["x_err_high"].to_numpy(float)]
    )

    fig = plt.figure(figsize=(13.333, 7.5), constrained_layout=False)
    fig.patch.set_facecolor("white")

    fig.text(
        0.045,
        0.955,
        r"Reco-cluster $E_T$ leakage check",
        fontsize=27,
        fontweight="bold",
        ha="left",
        va="top",
    )
    add_round_box(fig, (0.045, 0.731), (0.43, 0.132), "#f3f5f8")
    add_round_box(fig, (0.505, 0.731), (0.45, 0.132), "#eef8f1")

    fig.text(0.062, 0.846, "What is plotted", fontsize=16.0, fontweight="bold", va="top")
    fig.text(
        0.062,
        0.810,
        (
            ("Reco photon-cluster ET in true 1 GeV bins from 12-50 GeV.\n" if is_fine
             else "Reco photon-cluster ET from saved bins, shown from 12-50 GeV.\n") +
            "Black points are the weighted sum.\n"
            "Bottom panel shows each sample fraction."
        ),
        fontsize=13.4,
        va="top",
    )
    fig.text(0.522, 0.846, "Main takeaway", fontsize=16.0, fontweight="bold", va="top")
    fig.text(
        0.522,
        0.810,
        (
            f"Jet12 is absent above 22 GeV; Jet20 carries the transition.\n"
            f"J20 fraction: {jet20_2224} (22-24), {jet20_2426} (24-26), {jet20_2635} (26-35).\n"
            r"Jet20 ends at $p_T^{truth\,jet}=31$ GeV; Jet30 owns 31-41 GeV; Jet40 then takes over."
        ),
        fontsize=12.6,
        va="top",
    )

    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[3.1, 1.25],
        left=0.08,
        right=0.96,
        bottom=0.145,
        top=0.685,
        hspace=0.06,
    )
    ax = fig.add_subplot(gs[0])
    ax_frac = fig.add_subplot(gs[1], sharex=ax)

    for sample in SAMPLES + ["Sum"]:
        y = pivot_plot[f"weighted_entries_{sample}"].to_numpy(float)
        yerr = pivot_plot[f"weighted_error_{sample}"].to_numpy(float)
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            xerr=xerr,
            fmt=MARKERS[sample],
            color=COLORS[sample],
            ecolor=COLORS[sample],
            elinewidth=1.0,
            capsize=2.2,
            markersize=5.3 if sample != "Sum" else 5.8,
            markeredgewidth=1.0,
            markerfacecolor=COLORS[sample] if sample != "Sum" else "white",
            linestyle="none",
            label=LABELS[sample],
            zorder=6 if sample == "Sum" else 4,
        )

    for sample in SAMPLES:
        ax_frac.plot(
            x,
            pivot_plot[f"fraction_{sample}"].to_numpy(float),
            marker=MARKERS[sample],
            markersize=5.3,
            linestyle="none",
            color=COLORS[sample],
            label=sample,
        )

    ax.set_yscale("log")
    ax.set_xlim(12.0, 50.0)
    y_positive = pivot_plot[[f"weighted_entries_{s}" for s in SAMPLES + ["Sum"]]].to_numpy(float)
    y_positive = y_positive[y_positive > 0]
    y_min = max(1.0, np.nanmin(y_positive) * 0.45) if len(y_positive) else 1.0
    y_max = np.nanmax(y_positive) * 2.4 if len(y_positive) else 1.0e7
    ax.set_ylim(y_min, y_max)
    ax.set_ylabel(
        "weighted entries / 1 GeV bin"
        if is_fine else
        "weighted entries / saved bin",
        labelpad=10,
    )
    ax.grid(which="major", color="#d7dce3", linewidth=0.8, alpha=0.75)
    ax.grid(which="minor", color="#edf0f4", linewidth=0.45, alpha=0.55)
    ax.tick_params(which="both", direction="in", top=True, right=True, length=5)
    ax.tick_params(which="minor", length=2.5)
    ax.tick_params(labelbottom=False)
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.755, 0.965),
        ncol=2,
        frameon=True,
        facecolor="white",
        edgecolor="white",
        framealpha=0.88,
        fontsize=10.3,
        handlelength=1.55,
        columnspacing=0.95,
        handletextpad=0.45,
        borderpad=0.45,
        labelspacing=0.30,
    )
    add_in_axes_internal_label(ax)

    ax_frac.set_ylim(-0.035, 1.05)
    ax_frac.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_frac.set_ylabel("fraction of\nweighted sum", labelpad=10)
    ax_frac.set_xlabel(r"reco photon-cluster $E_T$ ($p_T^\gamma$) [GeV]", labelpad=2)
    ax_frac.grid(which="major", color="#d7dce3", linewidth=0.8, alpha=0.75)
    ax_frac.tick_params(which="both", direction="in", top=True, right=True, length=5)
    fig.text(
        0.08,
        0.044,
        "1 GeV reco-cluster ET bins, 12-50 GeV",
        fontsize=13.0,
        color="#5e6675",
        ha="left",
        va="top",
    )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_PATH, bbox_inches=None, facecolor="white")
    plt.close(fig)
    print(PNG_PATH)


if __name__ == "__main__":
    build_slide()
