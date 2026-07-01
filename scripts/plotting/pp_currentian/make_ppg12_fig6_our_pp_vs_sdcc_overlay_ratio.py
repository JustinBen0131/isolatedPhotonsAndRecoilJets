#!/usr/bin/env python3
"""Render PPG12 Fig. 6 source data against our pp RecoilJets output.

This intentionally mirrors the cleaned DataThief-vs-SDCC validation plot, but
replaces the DataThief-export points with the THE-85 nominal-only pp analysis
output.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
EFF_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/efficiency_stage"
SOURCE_CSV = EFF_DIR / "ppg12_fig6_direct_vs_the85_nominalonly_overlay_10to35.csv"
SOURCE_MANIFEST = EFF_DIR / "ppg12_fig6_direct_vs_the85_nominalonly_overlay_10to35_manifest.json"
OUT_PNG = EFF_DIR / "ppg12_fig6_our_pp_analysis_vs_sdcc_overlay_ratio.png"
OUT_MANIFEST = EFF_DIR / "ppg12_fig6_our_pp_analysis_vs_sdcc_overlay_ratio_manifest.json"

COLORS = {
    "reco": "#111111",
    "reco_id": "#d62d91",
    "reco_id_iso": "#2f9638",
}
LABELS = {
    "reco": r"$\varepsilon_{\mathrm{reco}}$",
    "reco_id": r"$\varepsilon_{\mathrm{reco}}\times\varepsilon_{\mathrm{ID}}$",
    "reco_id_iso": r"$\varepsilon_{\mathrm{reco}}\times\varepsilon_{\mathrm{ID}}\times\varepsilon_{\mathrm{iso}}$",
}


def load_rows() -> dict[str, list[dict[str, float]]]:
    out: dict[str, list[dict[str, float]]] = defaultdict(list)
    with SOURCE_CSV.open() as f:
        for row in csv.DictReader(f):
            stage = row["stage"]
            out[stage].append(
                {
                    "pt_low": float(row["pt_low"]),
                    "pt_high": float(row["pt_high"]),
                    "pt_mid": float(row["pt_mid"]),
                    "ppg12_eff": float(row["ppg12_eff"]),
                    "our_eff": float(row["the85_nominalonly_eff"]),
                    "ratio": float(row["the85_over_ppg12"]),
                }
            )
    return {k: sorted(v, key=lambda r: r["pt_mid"]) for k, v in out.items()}


def render() -> dict[str, dict[str, float]]:
    rows = load_rows()
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=200,
        sharex=True,
        gridspec_kw={"height_ratios": [3.25, 1.0], "hspace": 0.05},
    )

    summary: dict[str, dict[str, float]] = {}
    for stage in ("reco", "reco_id", "reco_id_iso"):
        color = COLORS[stage]
        stage_rows = rows[stage]
        x = np.asarray([r["pt_mid"] for r in stage_rows], dtype=float)
        xerr = np.asarray(
            [
                [r["pt_mid"] - r["pt_low"] for r in stage_rows],
                [r["pt_high"] - r["pt_mid"] for r in stage_rows],
            ],
            dtype=float,
        )
        ppg12 = np.asarray([r["ppg12_eff"] for r in stage_rows], dtype=float)
        ours = np.asarray([r["our_eff"] for r in stage_rows], dtype=float)
        ratio = ours / ppg12

        ax.errorbar(
            x,
            ppg12,
            xerr=xerr,
            fmt="o",
            ms=6.0,
            mfc="white",
            mec=color,
            mew=1.55,
            ecolor=color,
            elinewidth=1.0,
            capsize=0,
            alpha=0.96,
            zorder=5,
        )
        ax.plot(
            x,
            ours,
            "o",
            color=color,
            markerfacecolor=color,
            markeredgewidth=0.9,
            ms=4.0,
            linestyle="None",
            zorder=6,
        )
        rax.plot(
            x,
            ratio,
            "o",
            color=color,
            markerfacecolor=color,
            markeredgewidth=0.8,
            ms=4.2,
            linestyle="None",
        )
        summary[stage] = {
            "min_ratio": float(np.nanmin(ratio)),
            "max_ratio": float(np.nanmax(ratio)),
            "mean_ratio": float(np.nanmean(ratio)),
        }

    ax.set_xlim(10.0, 35.0)
    ax.set_ylim(0.0, 1.15)
    ax.set_ylabel("Efficiency", fontsize=17)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
    ax.text(0.055, 0.93, "sPHENIX", transform=ax.transAxes, fontsize=15.0, fontstyle="italic", fontweight="bold")
    ax.text(0.255, 0.93, "Internal", transform=ax.transAxes, fontsize=15.0)
    ax.text(0.055, 0.865, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=12.5)
    ax.text(0.735, 0.93, "PYTHIA8", transform=ax.transAxes, fontsize=15.0)
    ax.text(0.735, 0.865, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=12.5)

    source_handles = [
        Line2D([0], [0], marker="o", color="0.15", markerfacecolor="0.15", markeredgecolor="0.15", lw=0, ms=6.0, label="Our pp analysis (filled)"),
        Line2D([0], [0], marker="o", color="0.15", markerfacecolor="white", markeredgecolor="0.15", markeredgewidth=1.55, lw=0, ms=6.0, label="PPG12 SDCC source (open)"),
    ]
    color_handles = [
        Line2D([0], [0], marker="o", color=COLORS["reco"], markerfacecolor=COLORS["reco"], lw=0, ms=5.8, label=LABELS["reco"]),
        Line2D([0], [0], marker="o", color=COLORS["reco_id"], markerfacecolor=COLORS["reco_id"], lw=0, ms=5.8, label=LABELS["reco_id"]),
        Line2D([0], [0], marker="o", color=COLORS["reco_id_iso"], markerfacecolor=COLORS["reco_id_iso"], lw=0, ms=5.8, label=LABELS["reco_id_iso"]),
    ]
    source_legend = ax.legend(
        handles=source_handles,
        loc="lower left",
        bbox_to_anchor=(0.02, 0.185),
        frameon=False,
        fontsize=12.1,
        handlelength=1.0,
        borderpad=0.2,
        labelspacing=0.55,
    )
    ax.add_artist(source_legend)
    ax.legend(
        handles=color_handles,
        loc="lower left",
        bbox_to_anchor=(0.02, 0.015),
        frameon=False,
        fontsize=12.1,
        handlelength=1.0,
        borderpad=0.2,
        labelspacing=0.55,
    )

    rax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylim(0.55, 1.45)
    rax.set_ylabel("Our / PPG12", fontsize=12.5)
    rax.set_xlabel(r"$E_{\mathrm{T}}^{\gamma,\mathrm{truth}}$ [GeV]", fontsize=15)
    rax.minorticks_on()
    rax.tick_params(which="both", direction="in", top=True, right=True, labelsize=11)
    fig.subplots_adjust(left=0.14, right=0.97, top=0.98, bottom=0.09)
    fig.savefig(OUT_PNG)
    plt.close(fig)
    return summary


def main() -> None:
    summary = render()
    upstream = json.loads(SOURCE_MANIFEST.read_text()) if SOURCE_MANIFEST.exists() else {}
    OUT_MANIFEST.write_text(
        json.dumps(
            {
                "artifact": str(OUT_PNG),
                "source_csv": str(SOURCE_CSV),
                "source_manifest": str(SOURCE_MANIFEST),
                "comparison": "filled markers are THE-85 nominal-only pp RecoilJets output; open markers are PPG12 SDCC source values",
                "ratio_panel": "our pp analysis / PPG12 SDCC source",
                "upstream_scope": upstream.get("the85", {}),
                "ratio_summary": summary,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(OUT_PNG)
    print(OUT_MANIFEST)


if __name__ == "__main__":
    main()
