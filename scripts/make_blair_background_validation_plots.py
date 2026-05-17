#!/usr/bin/env python3
"""Make Blair-facing AuAu photon-ID background validation plots.

Inputs are count tables from ``scripts/audit_auau_single_bdt_sample_stats.py``.
The script can append a focused Jet30 count table to the existing Jet12+20
audit output, then recompute weighted spectra with the inclusive3 stitching
constants.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ET_BINS = ["15-17", "17-19", "19-21", "21-23", "23-25", "25-27", "27-30", "30-35"]
CENT7 = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-80"]
CENT3 = ["0-20", "20-50", "50-80"]
CENT3_MAP = {
    "0-10": "0-20",
    "10-20": "0-20",
    "20-30": "20-50",
    "30-40": "20-50",
    "40-50": "20-50",
    "50-60": "50-80",
    "60-80": "50-80",
}

SCALES_PB = {
    "embeddedPhoton12": 2598.12425 / 7685793.0,
    "embeddedPhoton20": 133.317866 / 8060542.0,
    "embeddedJet12": 1.21692467e6 / 7116378.0,
    "embeddedJet20": 5.44464934e4 / 7597337.0,
    "embeddedJet30": 2.40291630e3 / 8060542.0,
}
SIGNAL_SAMPLES = ["embeddedPhoton12", "embeddedPhoton20"]
BACKGROUND_SAMPLES = ["embeddedJet12", "embeddedJet20", "embeddedJet30"]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.1,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )


def read_counts(base_counts: Path, jet30_counts: Path | None) -> pd.DataFrame:
    frames = [pd.read_csv(base_counts)]
    if jet30_counts:
        frames.append(pd.read_csv(jet30_counts))
    df = pd.concat(frames, ignore_index=True)
    df = df[df["sample"].isin(SIGNAL_SAMPLES + BACKGROUND_SAMPLES)].copy()
    df["cent3"] = df["cent_bin"].map(CENT3_MAP)
    if df["cent3"].isna().any():
        bad = sorted(df.loc[df["cent3"].isna(), "cent_bin"].astype(str).unique())
        raise SystemExit(f"Unmapped centrality bins: {bad}")
    return df


def combine_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows = [df]
    for combo, parts, role in [
        ("embeddedPhoton12plus20", SIGNAL_SAMPLES, "signal"),
        ("embeddedJet12plus20plus30", BACKGROUND_SAMPLES, "background"),
    ]:
        sub = df[df["sample"].isin(parts)]
        grouped = (
            sub.groupby(["et_bin", "cent_bin", "cent3"], as_index=False)[["signal", "background", "entries"]]
            .sum()
        )
        grouped.insert(0, "role", role)
        grouped.insert(0, "sample", combo)
        grouped["background_to_signal"] = grouped["background"] / grouped["signal"].replace(0, np.nan)
        rows.append(grouped[df.columns])
    return pd.concat(rows, ignore_index=True)


def weighted_table(df: pd.DataFrame, samples: list[str], count_col: str, kind: str) -> pd.DataFrame:
    rows = []
    for cent in CENT3:
        for et in ET_BINS:
            local = df[(df["cent3"] == cent) & (df["et_bin"] == et) & (df["sample"].isin(samples))]
            y = 0.0
            e2 = 0.0
            raw = 0.0
            sample_counts = {}
            for sample in samples:
                n = float(local.loc[local["sample"] == sample, count_col].sum())
                w = SCALES_PB[sample]
                sample_counts[f"{sample}_raw"] = n
                raw += n
                y += n * w
                e2 += n * w * w
            rows.append(
                {
                    "kind": kind,
                    "cent3": cent,
                    "et_bin": et,
                    "raw_count": raw,
                    "weighted_yield_pb": y,
                    "weighted_yield_pb_err": math.sqrt(e2),
                    **sample_counts,
                }
            )
    return pd.DataFrame(rows)


def summed_table(cent3_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for kind in ["signal", "background"]:
        sub_kind = cent3_table[cent3_table["kind"] == kind]
        for et in ET_BINS:
            sub = sub_kind[sub_kind["et_bin"] == et]
            rows.append(
                {
                    "kind": kind,
                    "et_bin": et,
                    "raw_count": float(sub["raw_count"].sum()),
                    "weighted_yield_pb": float(sub["weighted_yield_pb"].sum()),
                    "weighted_yield_pb_err": math.sqrt(float((sub["weighted_yield_pb_err"] ** 2).sum())),
                }
            )
    return pd.DataFrame(rows)


def draw_header(fig, title: str, subtitle: str, y_title: float = 0.965) -> None:
    fig.suptitle(title, y=y_title, fontsize=19, fontweight="bold")
    fig.text(0.115, 0.925, "sPHENIX", ha="left", va="center", fontsize=13.5, fontweight="bold", fontstyle="italic")
    fig.text(0.205, 0.925, "Internal", ha="left", va="center", fontsize=13.5)
    fig.text(0.115, 0.882, subtitle, ha="left", va="center", fontsize=10.8)


def plot_background_panels(table: pd.DataFrame, out: Path) -> None:
    colors = {"0-20": "#B2182B", "20-50": "#D55E00", "50-80": "#2166AC"}
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 6.6), dpi=180, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.17, top=0.79, wspace=0.08)
    draw_header(
        fig,
        "Weighted Jet12+20+30 background spectrum by centrality",
        "Inclusive embedded jet stitched weights; marker-only spectra; errors = sqrt(sum w_i^2 N_i)",
    )
    for ax, cent in zip(axes, CENT3):
        sub = table[(table["kind"] == "background") & (table["cent3"] == cent)].set_index("et_bin").loc[ET_BINS]
        x = np.arange(len(ET_BINS))
        ax.errorbar(
            x,
            sub["weighted_yield_pb"],
            yerr=sub["weighted_yield_pb_err"],
            fmt="o",
            linestyle="none",
            color=colors[cent],
            markerfacecolor="white",
            markeredgewidth=2.2,
            markersize=7.0,
            elinewidth=1.6,
            capsize=4,
            capthick=1.6,
        )
        ax.set_title(f"{cent}% centrality", fontsize=16, fontweight="bold")
        ax.set_yscale("log")
        ax.set_xticks(x, ET_BINS, rotation=32, ha="right", fontsize=11.5)
        ax.set_xlabel(r"Candidate $E_T$ bin [GeV]", fontsize=13.5)
        ax.grid(axis="y", which="both", color="0.90", linewidth=0.8)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Weighted background candidates / bin [pb-equiv.]", fontsize=13.5)
    fig.savefig(out / "weighted_background_spectrum_cent3_markers_errorbars.png", bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)


def plot_background_overlay(table: pd.DataFrame, out: Path) -> None:
    colors = {"0-20": "#B2182B", "20-50": "#D55E00", "50-80": "#2166AC"}
    fig, ax = plt.subplots(figsize=(10.8, 6.7), dpi=180)
    fig.subplots_adjust(left=0.12, right=0.74, top=0.80, bottom=0.16)
    draw_header(
        fig,
        "Weighted Jet12+20+30 background spectra",
        "Three centrality bins overlaid; inclusive embedded jet stitched weights; no connected lines",
    )
    x = np.arange(len(ET_BINS))
    offsets = {"0-20": -0.09, "20-50": 0.0, "50-80": 0.09}
    for cent in CENT3:
        sub = table[(table["kind"] == "background") & (table["cent3"] == cent)].set_index("et_bin").loc[ET_BINS]
        ax.errorbar(
            x + offsets[cent],
            sub["weighted_yield_pb"],
            yerr=sub["weighted_yield_pb_err"],
            fmt="o",
            linestyle="none",
            color=colors[cent],
            markerfacecolor="white",
            markeredgewidth=2.2,
            markersize=7.2,
            elinewidth=1.6,
            capsize=4,
            capthick=1.6,
            label=f"{cent}% centrality",
        )
    ax.set_yscale("log")
    ax.set_xticks(x, ET_BINS, rotation=28, ha="right")
    ax.set_xlabel(r"Candidate $E_T$ bin [GeV]", fontsize=15)
    ax.set_ylabel("Weighted background candidates / bin [pb-equiv.]", fontsize=15)
    ax.grid(axis="y", which="both", color="0.90", linewidth=0.8)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), title="Centrality")
    fig.savefig(out / "weighted_background_spectrum_cent3_overlay_markers_errorbars.png", bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)


def plot_signal_background_panels(table: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 6.6), dpi=180, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.17, top=0.79, wspace=0.08)
    draw_header(
        fig,
        "Weighted signal and background spectra by centrality",
        "Signal = Photon12+20 labels; background = Jet12+20+30 labels; errors = sqrt(sum w_i^2 N_i)",
    )
    x = np.arange(len(ET_BINS))
    for ax, cent in zip(axes, CENT3):
        for kind, color, marker, label, dx in [
            ("signal", "#0072B2", "o", "Signal labels", -0.06),
            ("background", "#D55E00", "s", "Background labels", 0.06),
        ]:
            sub = table[(table["kind"] == kind) & (table["cent3"] == cent)].set_index("et_bin").loc[ET_BINS]
            ax.errorbar(
                x + dx,
                sub["weighted_yield_pb"],
                yerr=sub["weighted_yield_pb_err"],
                fmt=marker,
                linestyle="none",
                color=color,
                markerfacecolor="white",
                markeredgewidth=2.2,
                markersize=7.0,
                elinewidth=1.6,
                capsize=4,
                capthick=1.6,
                label=label,
            )
        ax.set_title(f"{cent}% centrality", fontsize=16, fontweight="bold")
        ax.set_yscale("log")
        ax.set_xticks(x, ET_BINS, rotation=32, ha="right", fontsize=11.5)
        ax.set_xlabel(r"Candidate $E_T$ bin [GeV]", fontsize=13.5)
        ax.grid(axis="y", which="both", color="0.90", linewidth=0.8)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Weighted candidates / bin [pb-equiv.]", fontsize=13.5)
    axes[-1].legend(frameon=False, loc="upper right")
    fig.savefig(out / "weighted_signal_background_spectrum_cent3_markers_errorbars.png", bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)


def plot_signal_background_summed(table: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(12.2, 7.0), dpi=180)
    fig.subplots_adjust(left=0.115, right=0.705, top=0.805, bottom=0.155)
    draw_header(
        fig,
        "Weighted signal and background spectra",
        "Signal Photon12+20; background Jet12+20+30 stitched weights; statistical errors from sqrt(sum w_i^2 N_i)",
    )
    x = np.arange(len(ET_BINS))
    for kind, color, marker, label, dx in [
        ("signal", "#0072B2", "o", "Signal labels, Photon12+20", -0.06),
        ("background", "#D55E00", "s", "Background labels, Jet12+20+30", 0.06),
    ]:
        sub = table[table["kind"] == kind].set_index("et_bin").loc[ET_BINS]
        ax.errorbar(
            x + dx,
            sub["weighted_yield_pb"],
            yerr=sub["weighted_yield_pb_err"],
            fmt=marker,
            linestyle="none",
            color=color,
            markerfacecolor="white",
            markeredgewidth=2.2,
            markersize=7.8,
            elinewidth=1.8,
            capsize=5,
            capthick=1.8,
            label=label,
        )
    ax.set_yscale("log")
    ax.set_xticks(x, ET_BINS)
    ax.set_xlabel(r"Candidate $E_T$ bin [GeV]", fontsize=17)
    ax.set_ylabel("Weighted candidates / bin [pb-equiv.]", fontsize=17)
    ax.grid(axis="y", which="both", color="0.90", linewidth=0.8)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.025, 0.995), title="Sample", fontsize=12.5, title_fontsize=12.5)
    fig.savefig(out / "weighted_signal_background_spectrum_summed_centrality_markers_errorbars.png", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def plot_heatmap(matrix_rows: pd.DataFrame, sample: str, value: str, outpath: Path, title: str, cmap: str) -> None:
    matrix = np.full((len(CENT7), len(ET_BINS)), np.nan)
    sub = matrix_rows[matrix_rows["sample"] == sample]
    for _, row in sub.iterrows():
        if row["cent_bin"] in CENT7 and row["et_bin"] in ET_BINS:
            matrix[CENT7.index(row["cent_bin"]), ET_BINS.index(row["et_bin"])] = float(row[value])
    fig, ax = plt.subplots(figsize=(10.6, 6.1), dpi=180)
    im = ax.imshow(matrix, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(ET_BINS)), ET_BINS, fontsize=11)
    ax.set_yticks(np.arange(len(CENT7)), CENT7, fontsize=11)
    ax.set_xlabel(r"Candidate $E_T$ bin [GeV]", fontsize=14)
    ax.set_ylabel("Centrality bin [%]", fontsize=14)
    ax.set_title(title, fontsize=17, fontweight="bold")
    finite_vals = matrix[np.isfinite(matrix)]
    threshold = float(np.nanmax(finite_vals) * 0.55) if finite_vals.size else 0.0
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            if np.isfinite(val):
                color = "white" if val > threshold else "black"
                ax.text(j, i, f"{val:,.0f}", ha="center", va="center", fontsize=8.6, color=color)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Label count", fontsize=12)
    fig.text(0.02, 0.97, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.12, 0.97, "Internal", ha="left", va="top", fontsize=13)
    fig.tight_layout(rect=(0.02, 0.02, 0.98, 0.91))
    fig.savefig(outpath, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-counts", required=True, type=Path)
    parser.add_argument("--jet30-counts", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    setup_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    centrality_dir = args.outdir / "centrality_count_views"
    centrality_dir.mkdir(parents=True, exist_ok=True)

    counts = combine_rows(read_counts(args.base_counts, args.jet30_counts))
    counts.to_csv(args.outdir / "counts_by_et_cent_single_bdt_centInput_pt1535_jet12plus20plus30.csv", index=False)

    signal = weighted_table(counts, SIGNAL_SAMPLES, "signal", "signal")
    background = weighted_table(counts, BACKGROUND_SAMPLES, "background", "background")
    weighted_cent3 = pd.concat([background, signal], ignore_index=True)
    weighted_cent3.to_csv(args.outdir / "weighted_signal_background_spectrum_cent3_counts.csv", index=False)
    background.to_csv(args.outdir / "weighted_background_spectrum_cent3_counts.csv", index=False)
    summed = summed_table(weighted_cent3)
    summed.to_csv(args.outdir / "weighted_signal_background_spectrum_summed_centrality_counts.csv", index=False)

    plot_background_panels(weighted_cent3, args.outdir)
    plot_background_overlay(weighted_cent3, args.outdir)
    plot_signal_background_panels(weighted_cent3, args.outdir)
    plot_signal_background_summed(summed, args.outdir)

    plot_heatmap(
        counts,
        "embeddedJet12plus20plus30",
        "background",
        centrality_dir / "background_count_tile_chart_et_centrality.png",
        r"Background label counts by $E_T$ and centrality",
        "Reds",
    )
    plot_heatmap(
        counts,
        "embeddedPhoton12plus20",
        "signal",
        centrality_dir / "signal_count_tile_chart_et_centrality.png",
        r"Signal label counts by $E_T$ and centrality",
        "Blues",
    )

    print(args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
