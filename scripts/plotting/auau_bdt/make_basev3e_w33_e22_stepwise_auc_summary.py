#!/usr/bin/env python3
"""Make a slide-ready stepwise AUC summary for base-v3E+w33+E22 ablations."""

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
from matplotlib.patches import FancyBboxPatch
import pandas as pd


ROOT = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
BASE_DIAG = ROOT / "validation/basev3e_controls_20260518_1110/validation_deep_diagnostics.json"
E22_METRICS = ROOT / "validation/basev3e_w33_e22ratio_20260518_192630/validation/validation_metrics.json"
OUT_DIR = ROOT / "slideReady/basev3e_w33_e22ratio_diagnostics"

CENTRALITY_KEY = "0_20"
CENTRALITY_LABEL = "0-20%"


def load_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def base_auc(base_diag: dict, product: str) -> tuple[float, int, int, int]:
    row = base_diag["products"][product]["auc_by_centrality"][CENTRALITY_KEY]
    return (
        float(row["auc"]),
        int(row["entries"]),
        int(row["signal_entries"]),
        int(row["background_entries"]),
    )


def e22_auc(e22_metrics: dict, product: str) -> float:
    return float(e22_metrics["products"][product]["auc_by_centrality"][CENTRALITY_KEY])


def add_sphenix_label(ax: plt.Axes) -> None:
    ax.text(
        0.885,
        0.975,
        "sPHENIX",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=18,
        fontstyle="italic",
        fontweight="bold",
        clip_on=False,
    )
    ax.text(
        0.888,
        0.975,
        " Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        clip_on=False,
    )


def build_rows() -> tuple[pd.DataFrame, tuple[int, int, int]]:
    base_diag = load_json(BASE_DIAG)
    e22_metrics = load_json(E22_METRICS)

    base_value, entries, sig_entries, bkg_entries = base_auc(base_diag, "baseBDT_v3E_withCentrality")
    width_value, _, _, _ = base_auc(base_diag, "baseBDT_v3E_withCentrality_w33")

    rows = [
        {
            "step": "reference",
            "model": "Base v3E + centrality",
            "short": "base",
            "inputs": 12,
            "auc": base_value,
            "increment_from": "none",
            "source": str(BASE_DIAG),
        },
        {
            "step": "width ablation",
            "model": "+ weta33/wphi33",
            "short": "w33",
            "inputs": 14,
            "auc": width_value,
            "increment_from": "base",
            "source": str(BASE_DIAG),
        },
        {
            "step": "E22 ratio ablation",
            "model": "+ E22/E37",
            "short": "E22/E37",
            "inputs": 15,
            "auc": e22_auc(e22_metrics, "baseBDT_v3E_withCentrality_w33_E22E37"),
            "increment_from": "w33",
            "source": str(E22_METRICS),
        },
        {
            "step": "E22 ratio ablation",
            "model": "+ E22/E53",
            "short": "E22/E53",
            "inputs": 15,
            "auc": e22_auc(e22_metrics, "baseBDT_v3E_withCentrality_w33_E22E53"),
            "increment_from": "w33",
            "source": str(E22_METRICS),
        },
        {
            "step": "E22 ratio ablation",
            "model": "+ E22/E37 + E22/E53",
            "short": "both ratios",
            "inputs": 16,
            "auc": e22_auc(e22_metrics, "baseBDT_v3E_withCentrality_w33_E22E37_E22E53"),
            "increment_from": "w33",
            "source": str(E22_METRICS),
        },
    ]
    df = pd.DataFrame(rows)
    df["delta_vs_base"] = df["auc"] - base_value
    df["delta_vs_base_percent"] = 100.0 * df["delta_vs_base"] / base_value
    df["delta_vs_w33"] = df["auc"] - width_value
    df["delta_vs_w33_percent"] = 100.0 * df["delta_vs_w33"] / width_value
    return df, (entries, sig_entries, bkg_entries)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df, counts = build_rows()
    entries, sig_entries, bkg_entries = counts
    base_value = float(df.loc[df["short"] == "base", "auc"].iloc[0])
    width_value = float(df.loc[df["short"] == "w33", "auc"].iloc[0])

    out_csv = OUT_DIR / "basev3e_w33_e22ratio_stepwise_auc_summary.csv"
    df.to_csv(out_csv, index=False)

    colors = {
        "base": "#6B7280",
        "w33": "#7C3AED",
        "E22/E37": "#2563EB",
        "E22/E53": "#F97316",
        "both ratios": "#16A34A",
    }
    edges = {
        "base": "#374151",
        "w33": "#5B21B6",
        "E22/E37": "#1E40AF",
        "E22/E53": "#C2410C",
        "both ratios": "#166534",
    }

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.labelsize": 16,
            "xtick.labelsize": 13.5,
            "ytick.labelsize": 14,
        }
    )

    fig = plt.figure(figsize=(16.0, 9.0), dpi=180)
    fig.patch.set_facecolor("white")

    # Slide-body layout: variable contract on the left, AUC response on the right.
    ax_vars = fig.add_axes([0.045, 0.185, 0.342, 0.655])
    ax = fig.add_axes([0.430, 0.185, 0.525, 0.655])
    for this_ax in (ax_vars, ax):
        this_ax.set_facecolor("white")

    plot_df = df.reset_index(drop=True)
    n_rows = len(plot_df)
    y_by_index = {i: n_rows - 1 - i for i in range(n_rows)}
    ymin, ymax = -0.52, n_rows - 0.10
    xmin = min(df["auc"]) - 0.0018
    xmax = max(df["auc"]) + 0.0064
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax_vars.set_xlim(0.0, 1.0)
    ax_vars.set_ylim(ymin, ymax)

    fig.text(
        0.045,
        0.925,
        f"{CENTRALITY_LABEL} centrality bin: targeted BDT input ablation",
        ha="left",
        va="bottom",
        fontsize=24,
        fontweight="bold",
        color="#0F172A",
    )
    fig.text(
        0.045,
        0.895,
        "Each row changes only the listed variables; AUC is evaluated on Photon12+20 signal vs Jet12+20+30 inclusive background.",
        ha="left",
        va="bottom",
        fontsize=13.8,
        color="#475569",
    )
    fig.text(
        0.045,
        0.865,
        "Read left to right: define the exact input change, then check whether signal/background ranking improves.",
        ha="left",
        va="bottom",
        fontsize=13.8,
        color="#475569",
    )

    ax_vars.text(
        0.00,
        ymax + 0.03,
        "Input contract",
        ha="left",
        va="bottom",
        fontsize=15,
        fontweight="bold",
        color="#0F172A",
        clip_on=False,
    )
    ax.text(
        xmin,
        ymax + 0.03,
        "Validation AUC response",
        ha="left",
        va="bottom",
        fontsize=15,
        fontweight="bold",
        color="#0F172A",
        clip_on=False,
    )

    variable_rows = {
        "base": {
            "heading": "Reference model",
            "variables": "base v3E + centrality",
            "note": "includes cluster_weta_cogx, cluster_wphi_cogx",
        },
        "w33": {
            "heading": "Add local 3x3 widths",
            "variables": "cluster_weta33_cogx, cluster_wphi33_cogx",
            "note": "same moment, restricted to the EMCal core",
        },
        "E22/E37": {
            "heading": "Add one E22 ratio",
            "variables": "e22_over_e37",
            "note": "energy sharing / shower compactness",
        },
        "E22/E53": {
            "heading": "Add one E22 ratio",
            "variables": "e22_over_e53",
            "note": "orthogonal E22 tail comparison",
        },
        "both ratios": {
            "heading": "Add both E22 ratios",
            "variables": "e22_over_e37 + e22_over_e53",
            "note": "best targeted ablation tested here",
        },
    }

    for i, (_, row) in enumerate(plot_df.iterrows()):
        y = y_by_index[i]
        short = str(row["short"])
        color = colors[short]
        edge = edges[short]
        band_color = "#F8FAFC" if i % 2 == 0 else "#FFFFFF"

        ax.axhspan(y - 0.39, y + 0.39, color=band_color, zorder=0)
        ax_vars.axhspan(y - 0.39, y + 0.39, color=band_color, zorder=0)

        box = FancyBboxPatch(
            (0.01, y - 0.305),
            0.955,
            0.61,
            boxstyle="round,pad=0.012,rounding_size=0.035",
            facecolor="white",
            edgecolor=color if short != "base" else "#94A3B8",
            linewidth=1.25 if short != "base" else 0.9,
            alpha=0.98,
            zorder=2,
        )
        ax_vars.add_patch(box)
        info = variable_rows[short]
        ax_vars.text(
            0.045,
            y + 0.155,
            info["heading"],
            ha="left",
            va="center",
            fontsize=12.9,
            fontweight="bold",
            color="#0F172A",
            zorder=3,
        )
        ax_vars.text(
            0.045,
            y - 0.035,
            info["variables"],
            ha="left",
            va="center",
            fontsize=10.8,
            color=edge if short != "base" else "#334155",
            fontfamily="DejaVu Sans Mono",
            zorder=3,
        )
        ax_vars.text(
            0.045,
            y - 0.205,
            info["note"],
            ha="left",
            va="center",
            fontsize=10.0,
            color="#64748B",
            zorder=3,
        )

    ax.axvline(
        base_value,
        color="#334155",
        linestyle=(0, (4, 4)),
        linewidth=2.0,
        alpha=0.9,
        zorder=1,
    )
    ax.axvline(
        width_value,
        color="#7C3AED",
        linestyle=(0, (2, 3)),
        linewidth=2.0,
        alpha=0.75,
        zorder=1,
    )
    ax.text(
        base_value,
        ymax - 0.10,
        f"base\n{base_value:.4f}",
        ha="center",
        va="top",
        fontsize=10.7,
        fontweight="bold",
        color="#334155",
        linespacing=0.95,
        bbox={
            "facecolor": "white",
            "edgecolor": "#CBD5E1",
            "boxstyle": "round,pad=0.22",
            "linewidth": 0.8,
            "alpha": 0.95,
        },
        zorder=5,
    )
    ax.text(
        width_value,
        ymax - 0.10,
        f"width baseline\n{width_value:.4f}",
        ha="center",
        va="top",
        fontsize=10.7,
        fontweight="bold",
        color="#6D28D9",
        linespacing=0.95,
        bbox={
            "facecolor": "white",
            "edgecolor": "#C4B5FD",
            "boxstyle": "round,pad=0.22",
            "linewidth": 0.8,
            "alpha": 0.95,
        },
        zorder=5,
    )

    for i, (_, row) in enumerate(plot_df.iterrows()):
        y = y_by_index[i]
        short = str(row["short"])
        auc = float(row["auc"])
        color = colors[short]
        edge = edges[short]

        if short == "base":
            start = auc
            delta_text = "reference"
        elif short == "w33":
            start = base_value
            delta_text = f"+{float(row['delta_vs_base_percent']):.2f}% from base"
        else:
            start = width_value
            delta_text = f"+{float(row['delta_vs_w33_percent']):.2f}% from width baseline"

        if short != "base":
            ax.hlines(
                y,
                min(start, auc),
                max(start, auc),
                color=color,
                linewidth=9.0 if short != "both ratios" else 10.5,
                alpha=0.28 if short != "both ratios" else 0.34,
                zorder=2,
            )
        ax.scatter(
            [auc],
            [y],
            s=325 if short != "both ratios" else 430,
            color=color,
            edgecolor=edge,
            linewidth=2.4,
            zorder=4,
        )

        x_text = auc + 0.00045
        label_ha = "left"
        if short == "both ratios":
            x_text = auc - 0.00055
            label_ha = "right"
        auc_font = 16.8 if short != "both ratios" else 18.4
        gain_font = 13.4 if short != "both ratios" else 14.2
        ax.text(
            x_text,
            y + 0.125,
            f"AUC {auc:.4f}",
            ha=label_ha,
            va="center",
            fontsize=auc_font,
            fontweight="bold",
            color="#111827",
            zorder=5,
        )
        ax.text(
            x_text,
            y - 0.175,
            delta_text,
            ha=label_ha,
            va="center",
            fontsize=gain_font,
            color=color if short != "base" else "#475569",
            zorder=5,
        )

    for this_ax in (ax_vars,):
        this_ax.set_xticks([])
        this_ax.set_yticks([])
        for spine in this_ax.spines.values():
            spine.set_visible(False)

    ax.set_yticks([])
    ax.set_xlabel(f"Validation AUC in {CENTRALITY_LABEL} centrality bin", labelpad=10)
    ax.grid(axis="x", color="#CBD5E1", linewidth=0.95, alpha=0.88)
    ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#0F172A")
    ax.spines["bottom"].set_linewidth(1.15)

    add_sphenix_label(ax)

    fig.text(
        0.045,
        0.105,
        (
            f"Centrality bin plotted: {CENTRALITY_LABEL}. Finite-score entries: {entries:,} "
            f"({sig_entries:,} signal, {bkg_entries:,} background)."
        ),
        ha="left",
        va="bottom",
        fontsize=10.8,
        color="#64748B",
    )
    fig.text(
        0.045,
        0.075,
        "Conclusion: local 3x3 widths make a cleaner baseline; both E22 ratios give the largest targeted AUC gain.",
        ha="left",
        va="bottom",
        fontsize=12.6,
        fontweight="bold",
        color="#0F172A",
    )

    out_png = OUT_DIR / "basev3e_w33_e22ratio_stepwise_auc_summary.png"
    fig.savefig(out_png, facecolor="white")
    plt.close(fig)

    print(out_png)
    print(out_csv)


if __name__ == "__main__":
    main()
