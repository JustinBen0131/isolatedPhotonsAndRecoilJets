#!/usr/bin/env python3
"""Make Blair-facing reco cluster ET source-composition plot for Jet12+20+30+40."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
OUTDIR = REPO / "dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_jet1234_20260526"
COARSE_CSV_PATH = OUTDIR / "inclusive_jet123_reco_cluster_et_weighted_components_cent_summed.csv"
COARSE_SUMMARY_PATH = OUTDIR / "inclusive_jet123_reco_cluster_et_weighted_components_summary.json"
FINE_CSV_PATH = OUTDIR / "inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50.csv"
FINE_SUMMARY_PATH = OUTDIR / "inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50_summary.json"
PNG_PATH = OUTDIR / "inclusive_jet1234_reco_cluster_et_weighted_components_blair_12to50.png"
COARSE_BACKUP_PATH = OUTDIR / "inclusive_jet1234_reco_cluster_et_weighted_components_blair_coarse_savedbins_backup.png"

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
            "font.family": "DejaVu Serif",
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.2,
            "axes.labelsize": 15,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 13,
            "figure.dpi": 160,
            "savefig.dpi": 220,
        }
    )


def weighted_pivot(df: pd.DataFrame) -> pd.DataFrame:
    pivot = (
        df.pivot_table(
            index=["x_low", "x_high", "x_center", "x_err_low", "x_err_high"],
            columns="sample",
            values=["weighted_entries", "weighted_error", "raw_entries"],
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
        for base in ["weighted_entries", "weighted_error", "raw_entries"]:
            col = f"{base}_{sample}"
            if col not in pivot:
                pivot[col] = 0.0
    pivot["weighted_entries_Sum"] = sum(
        pivot[f"weighted_entries_{sample}"] for sample in SAMPLES
    )
    pivot["weighted_error_Sum"] = np.sqrt(
        sum(pivot[f"weighted_error_{sample}"] ** 2 for sample in SAMPLES)
    )
    for sample in SAMPLES:
        denom = pivot["weighted_entries_Sum"].replace(0, np.nan)
        pivot[f"fraction_{sample}"] = pivot[f"weighted_entries_{sample}"] / denom
    return pivot


def make_summary(df: pd.DataFrame, pivot: pd.DataFrame, csv_path: Path) -> dict:
    metadata = (
        df[["sample", "truth_jet_gate", "sigma_eff_pb", "merge_scale", "root_path"]]
        .drop_duplicates()
        .to_dict(orient="records")
    )

    def row_for(lo: float, hi: float) -> dict:
        rows = pivot[(pivot["x_low"] >= lo) & (pivot["x_high"] <= hi)]
        weighted_entries = {
            sample: float(rows[f"weighted_entries_{sample}"].sum())
            for sample in SAMPLES
        }
        weighted_sum = sum(weighted_entries.values())
        return {
            "x_range": f"{lo:g}-{hi:g}",
            "weighted_sum": weighted_sum,
            "fractions": {
                sample: (
                    weighted_entries[sample] / weighted_sum
                    if weighted_sum > 0
                    else 0.0
                )
                for sample in SAMPLES
            },
            "weighted_entries": weighted_entries,
        }

    ranges = [row_for(22, 24), row_for(24, 26), row_for(26, 35), row_for(35, 40), row_for(40, 50)]

    return {
        "input_csv": str(csv_path),
        "output_png": str(PNG_PATH),
        "observable": "ABCD-summed reco photon-cluster ET, using pTgamma == cluster ET",
        "centrality_sum": "0-80%",
        "truth_jet_gates": {
            "Jet12": "12 <= pT_truth_jet < 21",
            "Jet20": "21 <= pT_truth_jet < 31",
            "Jet30": "31 <= pT_truth_jet < 41",
            "Jet40": "pT_truth_jet >= 41",
        },
        "sample_metadata": metadata,
        "bins_of_interest": ranges,
        "plot_caveat": "Rendered diagnostic view uses true 1 GeV bins from 12 to 50 GeV.",
    }


def input_paths() -> tuple[Path, Path, bool]:
    if FINE_CSV_PATH.exists():
        return FINE_CSV_PATH, FINE_SUMMARY_PATH, True
    return COARSE_CSV_PATH, COARSE_SUMMARY_PATH, False


def backup_existing_png() -> None:
    if PNG_PATH.exists() and not COARSE_BACKUP_PATH.exists() and not FINE_CSV_PATH.exists():
        COARSE_BACKUP_PATH.write_bytes(PNG_PATH.read_bytes())


def add_in_axes_internal_label(ax: plt.Axes) -> None:
    ax.text(
        0.985,
        0.965,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=14.0,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.8},
        zorder=30,
    )
    ax.text(
        0.985,
        0.912,
        "PYTHIA8 embedded inclusive jet, 0-80%",
        transform=ax.transAxes,
        fontsize=11.3,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
        zorder=30,
    )


def plot(df: pd.DataFrame, csv_path: Path) -> dict:
    pivot = weighted_pivot(df)
    summary = make_summary(df, pivot, csv_path)
    pivot_plot = pivot[pivot["x_high"] <= 50.0].copy()

    x = pivot_plot["x_center"].to_numpy(float)
    xerr = np.vstack(
        [pivot_plot["x_err_low"].to_numpy(float), pivot_plot["x_err_high"].to_numpy(float)]
    )

    fig = plt.figure(figsize=(14.6, 8.2), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[3.1, 1.35],
        left=0.082,
        right=0.985,
        bottom=0.13,
        top=0.755,
        hspace=0.06,
    )
    ax = fig.add_subplot(gs[0])
    ax_frac = fig.add_subplot(gs[1], sharex=ax)

    for sample in SAMPLES + ["Sum"]:
        y = pivot_plot[f"weighted_entries_{sample}"].to_numpy(float)
        yerr = pivot_plot[f"weighted_error_{sample}"].to_numpy(float)
        zorder = 6 if sample == "Sum" else 4
        lw = 2.6 if sample == "Sum" else 1.8
        ms = 6.4 if sample == "Sum" else 5.8
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            xerr=xerr,
            fmt=MARKERS[sample],
            color=COLORS[sample],
            ecolor=COLORS[sample],
            elinewidth=1.15,
            capsize=2.5,
            markersize=ms,
            markeredgewidth=1.0,
            markerfacecolor=COLORS[sample] if sample != "Sum" else "white",
            linewidth=lw,
            linestyle="none",
            label=LABELS[sample],
            zorder=zorder,
        )

    for sample in SAMPLES:
        frac = pivot_plot[f"fraction_{sample}"].to_numpy(float)
        ax_frac.plot(
            x,
            frac,
            marker=MARKERS[sample],
            markersize=6.0,
            linewidth=0.0,
            linestyle="none",
            color=COLORS[sample],
            label=sample,
        )

    ax.set_yscale("log")
    y_positive = pivot_plot[[f"weighted_entries_{s}" for s in SAMPLES + ["Sum"]]].to_numpy(float)
    y_positive = y_positive[y_positive > 0]
    y_min = max(1.0, np.nanmin(y_positive) * 0.45) if len(y_positive) else 1.0
    y_max = np.nanmax(y_positive) * 2.4 if len(y_positive) else 1.0e7
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(12.0, 50.0)
    ax.set_ylabel("weighted entries / 1 GeV bin", labelpad=12)
    ax.grid(which="major", color="#d4d8df", linewidth=0.8, alpha=0.75)
    ax.grid(which="minor", color="#edf0f4", linewidth=0.45, alpha=0.55)
    ax.tick_params(which="both", direction="in", top=True, right=True, length=6)
    ax.tick_params(which="minor", length=3)
    ax.tick_params(labelbottom=False)

    ax_frac.set_ylim(-0.035, 1.05)
    ax_frac.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_frac.set_ylabel("fraction of\nweighted sum", labelpad=10)
    ax_frac.set_xlabel(r"reco photon-cluster $E_T$  ($p_T^\gamma$) [GeV]")
    ax_frac.grid(which="major", color="#d4d8df", linewidth=0.8, alpha=0.75)
    ax_frac.tick_params(which="both", direction="in", top=True, right=True, length=6)

    meta = (
        "corrected truth-jet gates, weights applied\n"
        r"Jet12: $12 \leq p_T^{truth\,jet}<21$, "
        r"Jet20: $21 \leq p_T^{truth\,jet}<31$, "
        r"Jet30: $31 \leq p_T^{truth\,jet}<41$, "
        r"Jet40: $p_T^{truth\,jet}\geq41$"
    )
    fig.text(
        0.085,
        0.875,
        meta,
        fontsize=11.8,
        color="#222222",
        va="top",
        ha="left",
    )

    fig.text(
        0.085,
        0.965,
        r"Weighted reco cluster $E_T$ from inclusive-jet samples",
        fontsize=21,
        fontweight="bold",
        ha="left",
        va="top",
    )
    fig.text(
        0.085,
        0.922,
        (
            "ABCD A+B+C+D photon-candidate spectrum; all centralities summed."
        ),
        fontsize=13.0,
        ha="left",
        va="top",
    )
    legend = ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.99, 0.80),
        frameon=False,
        handlelength=2.8,
        borderaxespad=0.0,
    )
    for text in legend.get_texts():
        text.set_fontsize(13.5)
    add_in_axes_internal_label(ax)

    ax_frac.text(
        0.012,
        -0.36,
        "Dedicated diagnostic histogram: true 1 GeV reco-cluster ET bins, rendered from 12 to 50 GeV.",
        transform=ax_frac.transAxes,
        fontsize=10.3,
        color="#555b66",
        ha="left",
        va="top",
    )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_PATH, facecolor="white")
    plt.close(fig)
    return summary


def main() -> None:
    setup_style()
    csv_path, summary_path, is_fine = input_paths()
    backup_existing_png()
    df = pd.read_csv(csv_path)
    summary = plot(df, csv_path)
    print(f"Wrote {PNG_PATH}")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {summary_path}")
    print(json.dumps(summary["bins_of_interest"], indent=2))


if __name__ == "__main__":
    main()
