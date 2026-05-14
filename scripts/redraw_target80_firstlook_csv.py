#!/usr/bin/env python3
"""Redraw target80 first-look CSVs into cleaner slide-candidate PNGs."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import csv


REPO = Path(__file__).resolve().parents[1]
BASE = (
    REPO
    / "dataOutput/target80_all_available/bdt_target80_gated_20260512_001012/"
    / "analysis_config_expanded_5to40_target80/first_look_plots"
)
OUT = BASE / "slide_clean_redraws"

CENTS = [("0_20", "0-20% central"), ("20_50", "20-50% mid-central"), ("50_80", "50-80% peripheral")]

MODEL_STYLE = {
    "box": ("Box-cuts", "#111111", "o"),
    "flat050": ("Centrality-input BDT, score > 0.50", "#0072B2", "s"),
    "target80": ("Centrality-input BDT, target 80% signal efficiency", "#D55E00", "D"),
    "centInput": ("Centrality as input", "#0072B2", "o"),
    "minOpt": ("Minority-balanced", "#009E73", "P"),
    "cent3": ("3 centrality-bin BDTs", "#E69F00", "s"),
    "cent7": ("7 centrality-bin BDTs", "#6F3FB5", "^"),
    "ptCent3": (r"$E_T$ x 3 centrality-bin BDTs", "#D55E00", "D"),
}


def setup_axes(ax, ylabel: str | None, xlabel: str | None):
    ax.tick_params(direction="in", which="both", top=True, right=True, length=6)
    ax.minorticks_on()
    ax.grid(True, color="#D0D0D0", alpha=0.45, linewidth=0.7)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=15)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=15)
    ax.tick_params(labelsize=12)


def header(ax, y=0.90, line2=r"Pythia overlay, $\sqrt{s_{NN}} = 200$ GeV"):
    ax.text(0.055, y, "sPHENIX", transform=ax.transAxes, ha="left", va="top",
            fontsize=13, fontweight="bold", fontstyle="italic")
    ax.text(0.265, y, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=13)
    ax.text(0.055, y - 0.095, line2, transform=ax.transAxes, ha="left", va="top", fontsize=11)
    ax.text(0.055, y - 0.160, r"$R = 0.3$, sliding isolation", transform=ax.transAxes,
            ha="left", va="top", fontsize=11)


def redraw_points(csv_name: str, out_name: str, model_keys: list[str], ylabel: str,
                  ylim: tuple[float, float], x_label: str, xlim=(4.8, 35.6)):
    df = read_point_csv(BASE / csv_name)
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.2), sharey=True)
    shifts = {key: shift for key, shift in zip(model_keys, pd.Series(range(len(model_keys))).sub((len(model_keys) - 1) / 2) * 0.12)}
    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        for key in model_keys:
            label, color, marker = MODEL_STYLE[key]
            sub = df[(df["model_key"] == key) & (df["centrality"] == cent)].copy()
            if sub.empty:
                continue
            sub["x"] = sub["pt_mid"] + shifts[key]
            ax.errorbar(
                sub["x"],
                sub["value"],
                yerr=sub["error"],
                fmt=marker,
                linestyle="none",
                color=color,
                markerfacecolor=color,
                markeredgecolor=color,
                markersize=6.8,
                elinewidth=1.3,
                capsize=2.3,
                label=label,
                zorder=3,
            )
        ax.set_title(title, fontsize=17, pad=12)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        setup_axes(ax, ylabel if ic == 0 else None, x_label if ic == 1 else None)
        if ic == 0:
            header(ax)
        if ic == 2:
            leg = ax.legend(loc="upper right", frameon=True, framealpha=0.92, fontsize=10.5,
                            borderpad=0.55, handletextpad=0.45, ncol=1)
            leg.get_frame().set_linewidth(0.0)
    fig.tight_layout(w_pad=1.5)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / out_name
    fig.savefig(path, dpi=190)
    plt.close(fig)
    return path


def read_point_csv(path: Path) -> pd.DataFrame:
    rows = []
    with path.open(newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return pd.DataFrame()
        for raw in reader:
            if len(raw) == len(header):
                rows.append(dict(zip(header, raw)))
                continue
            # The legacy CSV writer did not quote labels containing commas.
            # Repair by treating the final six columns as fixed numeric fields.
            if len(raw) > len(header):
                repaired = [
                    raw[0],
                    raw[1],
                    ",".join(raw[2:-6]),
                    raw[-6],
                    raw[-5],
                    raw[-4],
                    raw[-3],
                    raw[-2],
                    raw[-1],
                ]
                rows.append(dict(zip(header, repaired)))
    df = pd.DataFrame(rows)
    for col in ["pt_mid", "numerator", "denominator", "value", "error"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def tradeoff(out_name: str):
    df = pd.read_csv(BASE / "target80_signal_recovery_vs_background_tight_fraction.csv")
    model_keys = ["box", "centInput", "minOpt", "cent3", "cent7", "ptCent3"]
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.2), sharey=True)
    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        for key in model_keys:
            label, color, marker = MODEL_STYLE[key]
            sub = df[(df["model_key"] == key) & (df["centrality"] == cent)]
            if sub.empty:
                continue
            row = sub.iloc[0]
            ax.errorbar(
                row["background_tight_fraction"],
                row["signal_recovery"],
                xerr=row["background_tight_fraction_err"],
                yerr=row["signal_recovery_err"],
                fmt=marker,
                linestyle="none",
                color=color,
                markersize=8.0,
                elinewidth=1.3,
                capsize=2.5,
                label=label,
                zorder=3,
            )
        ax.set_title(title, fontsize=17, pad=12)
        ax.set_xlim(0.0, 0.43)
        ax.set_ylim(0.0, 0.50)
        setup_axes(ax, "Truth-photon recovery" if ic == 0 else None,
                   "Inclusive-jet tight fraction" if ic == 1 else None)
        if ic == 0:
            header(ax, y=0.94, line2=r"Pythia overlay, $15 < E_T^\gamma < 35$ GeV")
        if ic == 2:
            leg = ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=10.0,
                            borderpad=0.55, handletextpad=0.45)
            leg.get_frame().set_linewidth(0.0)
    fig.tight_layout(w_pad=1.4)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / out_name
    fig.savefig(path, dpi=190)
    plt.close(fig)
    return path


def metric_table():
    df = pd.read_csv(BASE / "target80_integrated_metric_summary_15to35.csv")
    keep = ["box", "centInput", "minOpt", "cent3", "cent7", "ptCent3"]
    df = df[df["model_key"].isin(keep)].copy()
    df["signal_recovery_gain_vs_box"] = 0.0
    df["jet_tight_fraction_change_vs_box"] = 0.0
    for cent, _ in CENTS:
        box = df[(df["centrality"] == cent) & (df["model_key"] == "box")].iloc[0]
        mask = df["centrality"] == cent
        df.loc[mask, "signal_recovery_gain_vs_box"] = df.loc[mask, "signal_recovery"] - box["signal_recovery"]
        df.loc[mask, "jet_tight_fraction_change_vs_box"] = df.loc[mask, "background_tight_fraction"] - box["background_tight_fraction"]
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "target80_integrated_metrics_for_slide.csv"
    df.to_csv(out, index=False)
    return out


def main():
    paths = []
    paths.append(redraw_points(
        "centinput_reference_flat050_target80_truth_photon_recovery_1x3.csv",
        "slide_recovery_box_vs_flat050_vs_target80_centinput.png",
        ["box", "flat050", "target80"],
        "Truth-photon recovery",
        (0.0, 0.64),
        r"Truth photon $E_T$ [GeV]",
    ))
    paths.append(redraw_points(
        "centinput_reference_flat050_target80_inclusive_jet_tight_fraction_1x3.csv",
        "slide_jet_tight_fraction_box_vs_flat050_vs_target80_centinput.png",
        ["box", "flat050", "target80"],
        "Inclusive-jet tight fraction",
        (0.0, 0.58),
        r"Photon candidate $E_T$ [GeV]",
    ))
    paths.append(redraw_points(
        "target80_model_family_truth_photon_recovery_1x3.csv",
        "slide_recovery_target80_family_compact.png",
        ["box", "centInput", "minOpt", "cent3", "cent7", "ptCent3"],
        "Truth-photon recovery",
        (0.0, 0.64),
        r"Truth photon $E_T$ [GeV]",
    ))
    paths.append(redraw_points(
        "target80_model_family_truth_photon_recovery_1x3.csv",
        "slide_recovery_box_vs_centinput_vs_etcent3_target80.png",
        ["box", "centInput", "ptCent3"],
        "Truth-photon recovery",
        (0.0, 0.64),
        r"Truth photon $E_T$ [GeV]",
    ))
    family_leakage_csv = BASE / "target80_model_family_inclusive_jet_tight_fraction_1x3.csv"
    if family_leakage_csv.stat().st_size > 0:
        paths.append(redraw_points(
            "target80_model_family_inclusive_jet_tight_fraction_1x3.csv",
            "slide_jet_tight_fraction_target80_family_compact.png",
            ["box", "centInput", "minOpt", "cent3", "cent7", "ptCent3"],
            "Inclusive-jet tight fraction",
            (0.0, 0.58),
            r"Photon candidate $E_T$ [GeV]",
        ))
    paths.append(tradeoff("slide_signal_recovery_vs_jet_leakage_target80_compact.png"))
    paths.append(metric_table())
    print("[DONE] wrote clean target80 redraws:")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
