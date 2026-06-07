#!/usr/bin/env python3
"""Render full-count total-calo event-veto histograms from score-cache summary."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np


OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
HIST_JSON = OUTDIR / "the32_low_calo_full_count_histograms_v1.json"
PNG = OUTDIR / "the32_low_calo_full_count_histograms_5pct_v16_20260606.png"
PHONE_PNG = OUTDIR / "the32_low_calo_full_count_histograms_5pct_v16_phone_refresh_20260606.png"
MANIFEST = OUTDIR / "the32_low_calo_full_count_histograms_5pct_v16_20260606.json"

BLUE = "#1f77b4"
RED = "#d62728"
INK = "#172033"
MUTED = "#5d6b7a"
BOX_FACE = "#F8FAFC"
BOX_EDGE = "#CBD5E1"
GOLD = "#B7791F"
GOLD_SOFT = "#FFF7ED"
GOLD_LINE = "#F1C27D"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def step_xy(edges: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.repeat(edges, 2)[1:-1], np.repeat(counts, 2)


def compact_count(value: int) -> str:
    return f"{value:,}"


def render() -> None:
    payload = json.loads(HIST_JSON.read_text())
    edges = np.asarray(payload["energy_edges"], dtype="float64")
    panels = payload["panels"]
    ymax = max(
        max(max(panel["retained_hist"]), max(panel["removed_hist"]))
        for panel in panels
    )

    fig, axes = plt.subplots(4, 4, figsize=(12.8, 7.2), dpi=200, sharex=True, sharey=True)
    fig.patch.set_facecolor("white")

    for ax, panel in zip(axes.ravel(), panels):
        retained = np.asarray(panel["retained_hist"], dtype="int64")
        removed = np.asarray(panel["removed_hist"], dtype="int64")
        bx, by = step_xy(edges, retained)
        rx, ry = step_xy(edges, removed)

        ax.set_facecolor("white")
        ax.set_yscale("log")
        ax.plot(bx, np.maximum(by, 0.9), color=BLUE, lw=1.9, label="retained")
        ax.plot(rx, np.maximum(ry, 0.9), color=RED, lw=1.9, label="removed")
        ax.fill_between(bx, np.maximum(by, 0.9), 0.9, color=BLUE, alpha=0.06)
        ax.fill_between(rx, np.maximum(ry, 0.9), 0.9, color=RED, alpha=0.08)
        ax.axvline(panel["threshold"], color="white", lw=4.2, ls=(0, (3.0, 3.0)), zorder=4)
        ax.axvline(panel["threshold"], color=GOLD, lw=2.15, ls=(0, (3.0, 3.0)), zorder=5)

        lo = int(panel["cent_lo"])
        hi = int(panel["cent_hi"])
        total = int(panel["event_total"])
        kept = int(panel["event_retained"])
        cut = int(panel["event_removed"])
        frac = 100.0 * float(panel["event_removed_fraction"])

        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0.8, max(12.0, ymax * 2.6))
        ax.set_title(f"{lo}-{hi}%", fontsize=10.8, fontweight="bold", color=INK, pad=14)
        ax.text(
            0.00,
            1.075,
            f"N={compact_count(total)}",
            transform=ax.transAxes,
            fontsize=7.55,
            color=INK,
            fontweight="bold",
            ha="left",
            va="bottom",
            clip_on=False,
            zorder=10,
        )
        ax.text(
            0.36,
            1.075,
            f"blue {compact_count(kept)}",
            transform=ax.transAxes,
            fontsize=7.35,
            color=BLUE,
            ha="left",
            va="bottom",
            clip_on=False,
            zorder=10,
        )
        ax.text(
            1.00,
            1.075,
            f"red {compact_count(cut)} ({frac:.2f}%)",
            transform=ax.transAxes,
            fontsize=7.35,
            color=RED if cut else MUTED,
            fontweight="bold" if cut else "normal",
            ha="right",
            va="bottom",
            clip_on=False,
            zorder=10,
        )
        ax.set_xticks([2.8, 3.0, 3.2, 3.4])
        ax.tick_params(axis="both", labelsize=8.0, colors=INK, length=3.2, width=0.85)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_linewidth(0.95)
            spine.set_color(INK)

    for row in range(4):
        axes[row, 0].set_ylabel("raw event counts\nper energy bin", fontsize=8.8, color=INK)

    total = sum(int(panel["event_total"]) for panel in panels)
    removed = sum(int(panel["event_removed"]) for panel in panels)

    fig.text(
        0.055,
        0.965,
        "Full-count total-calo energy counts by 5% centrality bin",
        fontsize=24.7,
        fontweight="bold",
        color=INK,
        va="top",
    )
    formula_box = FancyBboxPatch(
        (0.052, 0.697),
        0.905,
        0.183,
        boxstyle="round,pad=0.008,rounding_size=0.010",
        transform=fig.transFigure,
        facecolor=BOX_FACE,
        edgecolor=BOX_EDGE,
        linewidth=1.1,
        zorder=1,
    )
    fig.patches.append(formula_box)
    formula_chip = FancyBboxPatch(
        (0.064, 0.808),
        0.885,
        0.046,
        boxstyle="round,pad=0.004,rounding_size=0.006",
        transform=fig.transFigure,
        facecolor=GOLD_SOFT,
        edgecolor=GOLD_LINE,
        linewidth=0.9,
        zorder=2,
    )
    fig.patches.append(formula_chip)
    fig.text(
        0.066,
        0.840,
        "Applied event-energy veto",
        fontsize=14.1,
        fontweight="bold",
        color=GOLD,
        va="top",
        zorder=3,
    )
    fig.text(
        0.612,
        0.874,
        "dashed gold = threshold in that 5% bin; red = removed before training",
        fontsize=9.45,
        fontweight="bold",
        color=GOLD,
        ha="left",
        va="top",
        zorder=3,
    )
    for x0, x1 in [(0.550, 0.561), (0.570, 0.581), (0.590, 0.601)]:
        fig.add_artist(
            Line2D(
                [x0, x1],
                [0.868, 0.868],
                transform=fig.transFigure,
                color=GOLD,
                lw=2.2,
                solid_capstyle="butt",
                zorder=3,
            )
        )
    fig.text(
        0.506,
        0.832,
        r"$y=\log_{10}(E_{\rm CEMC}+E_{\rm IHCal}+E_{\rm OHCal}+1);\quad \mathrm{remove\ event\ if}\ y<T_{\rm bin}$",
        fontsize=12.4,
        fontweight="bold",
        color=INK,
        ha="center",
        va="center",
        zorder=3,
    )
    row_rule = FancyBboxPatch(
        (0.064, 0.772),
        0.625,
        0.025,
        boxstyle="round,pad=0.002,rounding_size=0.004",
        transform=fig.transFigure,
        facecolor="white",
        edgecolor="#E2E8F0",
        linewidth=0.65,
        zorder=2,
    )
    fig.patches.append(row_rule)
    row_definitions = FancyBboxPatch(
        (0.064, 0.740),
        0.625,
        0.025,
        boxstyle="round,pad=0.002,rounding_size=0.004",
        transform=fig.transFigure,
        facecolor="white",
        edgecolor="#E2E8F0",
        linewidth=0.65,
        zorder=2,
    )
    fig.patches.append(row_definitions)
    row_derivation = FancyBboxPatch(
        (0.064, 0.708),
        0.625,
        0.025,
        boxstyle="round,pad=0.002,rounding_size=0.004",
        transform=fig.transFigure,
        facecolor="white",
        edgecolor="#E2E8F0",
        linewidth=0.65,
        zorder=2,
    )
    fig.patches.append(row_derivation)
    stat_card = FancyBboxPatch(
        (0.712, 0.708),
        0.245,
        0.089,
        boxstyle="round,pad=0.004,rounding_size=0.006",
        transform=fig.transFigure,
        facecolor="#FFF5F5",
        edgecolor="#FECACA",
        linewidth=0.85,
        zorder=2,
    )
    fig.patches.append(stat_card)
    fig.text(
        0.074,
        0.790,
        "Threshold",
        fontsize=9.0,
        fontweight="bold",
        color=GOLD,
        va="top",
        zorder=4,
    )
    fig.text(
        0.182,
        0.790,
        r"Per 5% bin: $T_{\rm bin}=\max[\mathrm{median}(y)-5\cdot1.4826\,\mathrm{MAD}(y),\,q_{0.1\%}(y)]$",
        fontsize=9.25,
        fontweight="bold",
        color=INK,
        va="top",
        zorder=3,
    )
    fig.text(
        0.074,
        0.758,
        "Definitions",
        fontsize=9.0,
        fontweight="bold",
        color=GOLD,
        va="top",
        zorder=4,
    )
    fig.text(
        0.182,
        0.758,
        r"MAD = median absolute deviation; $q_{0.1\%}$ = lower 0.1% percentile",
        fontsize=9.15,
        color=INK,
        fontweight="bold",
        va="top",
        zorder=4,
    )
    fig.text(
        0.074,
        0.726,
        "Derivation",
        fontsize=9.0,
        fontweight="bold",
        color=GOLD,
        va="top",
        zorder=4,
    )
    fig.text(
        0.182,
        0.726,
        "Lower edge of the normal event-energy band; blind to merge weights and sample type.",
        fontsize=9.0,
        color=MUTED,
        va="top",
        zorder=4,
    )
    fig.text(
        0.835,
        0.784,
        "Full counts",
        fontsize=9.2,
        color=RED,
        fontweight="bold",
        ha="center",
        va="top",
        zorder=4,
    )
    fig.text(
        0.835,
        0.763,
        f"{total:,} events",
        fontsize=9.2,
        color=INK,
        fontweight="bold",
        ha="center",
        va="top",
        zorder=4,
    )
    fig.text(
        0.835,
        0.741,
        f"{removed:,} removed ({100.0 * removed / total:.2f}%)",
        fontsize=9.2,
        fontweight="bold",
        color=RED,
        ha="center",
        va="top",
        zorder=4,
    )
    fig.text(
        0.835,
        0.722,
        "raw event histograms",
        fontsize=8.0,
        color=MUTED,
        ha="center",
        va="top",
        zorder=4,
    )
    fig.text(
        0.535,
        0.055,
        r"Plotted event energy:  $\log_{10}(E_{\rm CEMC}+E_{\rm IHCal}+E_{\rm OHCal}+1)$",
        fontsize=13.0,
        color=INK,
        ha="center",
        va="bottom",
    )

    fig.subplots_adjust(left=0.075, right=0.982, top=0.620, bottom=0.140, wspace=0.125, hspace=0.880)
    fig.savefig(PNG, dpi=200)
    fig.savefig(PHONE_PNG, dpi=200)
    plt.close(fig)

    manifest = {
        "png": str(PNG),
        "phone_refresh_png": str(PHONE_PNG),
        "histogram_source": str(HIST_JSON),
        "source_cache_list": payload["source_cache_list"],
        "cache_count": payload["cache_count"],
        "centrality_domain": payload["centrality_domain"],
        "cut_variable": payload["cut_variable"],
        "event_total_0_80": total,
        "event_removed_0_80": removed,
        "event_removed_fraction_0_80": removed / total,
        "display_contract": "Full retained and full removed raw event-count histograms; no retained-event display subsampling.",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    render()
