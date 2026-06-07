#!/usr/bin/env python3
"""Make raw-count line histograms for the total-calo event-veto diagnostic."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np


OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
MARKER_JSON = OUTDIR / "the32_low_calo_marker_sample_v1.json"
CUT_JSON = OUTDIR / "the32_low_calo_cut_v1.json"
PNG = OUTDIR / "the32_low_calo_energy_raw_count_histograms_5pct_cutlines_v20_20260606.png"
PHONE_PNG = OUTDIR / "the32_low_calo_energy_raw_count_histograms_5pct_cutlines_v20_phone_refresh_20260606.png"
MANIFEST = OUTDIR / "the32_low_calo_energy_raw_count_histograms_5pct_cutlines_v20_20260606.json"

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


def load_payloads() -> tuple[dict, dict]:
    return json.loads(MARKER_JSON.read_text()), json.loads(CUT_JSON.read_text())


def arrays(payload: dict, key: str) -> tuple[np.ndarray, np.ndarray]:
    markers = payload["markers"][key]
    return (
        np.asarray(markers["centrality"], dtype=np.float32),
        np.asarray(markers["log10_calo"], dtype=np.float32),
    )


def step_xy(edges: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.repeat(edges, 2)[1:-1], np.repeat(counts, 2)


def render() -> None:
    payload, cut = load_payloads()
    kept_x, kept_y = arrays(payload, "retained")
    removed_x, removed_y = arrays(payload, "removed")
    envelope = cut["envelope"]
    centrality_bins = [(lo, lo + 5) for lo in range(0, 80, 5)]
    energy_edges = np.linspace(2.70, 3.42, 73)

    prepared: list[dict] = []
    ymax = 0
    for lo, hi in centrality_bins:
        kept = kept_y[(kept_x >= lo) & (kept_x < hi) & np.isfinite(kept_y)]
        removed = removed_y[(removed_x >= lo) & (removed_x < hi) & np.isfinite(removed_y)]
        kept_counts, _ = np.histogram(kept, bins=energy_edges)
        removed_counts, _ = np.histogram(removed, bins=energy_edges)
        thresholds = [
            float(row["threshold"])
            for row in envelope
            if float(row["cent_lo"]) >= lo and float(row["cent_hi"]) <= hi
        ]
        ymax = max(ymax, int(kept_counts.max(initial=0)), int(removed_counts.max(initial=0)))
        prepared.append(
            {
                "lo": lo,
                "hi": hi,
                "kept": kept,
                "removed": removed,
                "kept_counts": kept_counts,
                "removed_counts": removed_counts,
                "thresholds": thresholds,
            }
        )

    fig, axes = plt.subplots(4, 4, figsize=(12.8, 7.2), dpi=200, sharex=True, sharey=True)
    fig.patch.set_facecolor("white")
    for ax, panel in zip(axes.ravel(), prepared):
        ax.set_facecolor("white")
        kx, ky = step_xy(energy_edges, panel["kept_counts"])
        rx, ry = step_xy(energy_edges, panel["removed_counts"])
        ax.plot(kx, ky, color=BLUE, lw=2.0)
        ax.plot(rx, ry, color=RED, lw=2.0)
        ax.fill_between(kx, ky, color=BLUE, alpha=0.07)
        ax.fill_between(rx, ry, color=RED, alpha=0.07)
        for threshold in panel["thresholds"]:
            ax.axvline(threshold, color="white", lw=4.3, ls=(0, (3.2, 3.0)), zorder=4)
            ax.axvline(threshold, color=GOLD, lw=2.2, ls=(0, (3.2, 3.0)), zorder=5)
        ax.set_xlim(energy_edges[0], energy_edges[-1])
        ax.set_ylim(0, ymax * 1.12)
        ax.set_title(f"{panel['lo']}-{panel['hi']}%", fontsize=10.7, fontweight="bold", color=INK, pad=6)
        threshold = panel["thresholds"][0] if panel["thresholds"] else np.nan
        label_x = 0.955 if np.isfinite(threshold) and threshold < 3.0 else 0.045
        label_ha = "right" if label_x > 0.5 else "left"
        ax.text(
            label_x,
            0.895,
            f"blue sample: {panel['kept'].size:,}",
            transform=ax.transAxes,
            fontsize=7.8,
            color=BLUE,
            ha=label_ha,
            va="top",
            zorder=10,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.9},
        )
        ax.text(
            label_x,
            0.770,
            f"red removed: {panel['removed'].size:,}",
            transform=ax.transAxes,
            fontsize=7.8,
            color=RED if panel["removed"].size else MUTED,
            fontweight="bold" if panel["removed"].size else "normal",
            ha=label_ha,
            va="top",
            zorder=10,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.9},
        )
        ax.set_xticks([2.8, 3.0, 3.2, 3.4])
        ax.tick_params(axis="both", labelsize=8.0, colors=INK, length=3.2, width=0.85)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_linewidth(0.95)
            spine.set_color(INK)

    for row in range(4):
        axes[row, 0].set_ylabel("raw counts per\nenergy bin", fontsize=8.9, color=INK)
    fig.text(
        0.055,
        0.965,
        "Raw total-calo energy counts by 5% centrality bin",
        fontsize=25.8,
        fontweight="bold",
        color=INK,
        va="top",
    )
    formula_box = FancyBboxPatch(
        (0.052, 0.720),
        0.905,
        0.155,
        boxstyle="round,pad=0.008,rounding_size=0.010",
        transform=fig.transFigure,
        facecolor=BOX_FACE,
        edgecolor=BOX_EDGE,
        linewidth=1.1,
        zorder=1,
    )
    fig.patches.append(formula_box)
    formula_chip = FancyBboxPatch(
        (0.064, 0.795),
        0.885,
        0.045,
        boxstyle="round,pad=0.004,rounding_size=0.006",
        transform=fig.transFigure,
        facecolor=GOLD_SOFT,
        edgecolor=GOLD_LINE,
        linewidth=0.85,
        zorder=2,
    )
    fig.patches.append(formula_chip)
    fig.text(
        0.066,
        0.828,
        "Applied event-energy veto",
        fontsize=14.1,
        fontweight="bold",
        color=GOLD,
        va="top",
        zorder=3,
    )
    fig.text(
        0.944,
        0.858,
        "dashed gold = threshold in that 5% bin; red = removed before training",
        fontsize=11.6,
        fontweight="bold",
        color=GOLD,
        ha="right",
        va="top",
        zorder=3,
    )
    fig.text(
        0.506,
        0.818,
        r"$T_{\rm bin}=\max[\mathrm{median}(y)-5\cdot1.4826\,\mathrm{MAD}(y),\,q_{0.1\%}(y)]$",
        fontsize=13.2,
        fontweight="bold",
        color=INK,
        ha="center",
        va="center",
        zorder=3,
    )
    fig.text(
        0.066,
        0.775,
        r"Rule: red events are removed when the plotted event energy is below $T_{\rm bin}$",
        fontsize=11.4,
        fontweight="bold",
        color=RED,
        va="top",
        zorder=3,
    )
    definition_strip = FancyBboxPatch(
        (0.062, 0.727),
        0.885,
        0.035,
        boxstyle="round,pad=0.002,rounding_size=0.004",
        transform=fig.transFigure,
        facecolor=BOX_FACE,
        edgecolor="none",
        linewidth=0.0,
        zorder=2,
    )
    fig.patches.append(definition_strip)
    fig.text(
        0.066,
        0.752,
        r"$\mathrm{MAD}$ = median absolute deviation, scaled as event-band width",
        fontsize=10.5,
        color=MUTED,
        va="top",
        zorder=4,
    )
    fig.text(
        0.530,
        0.752,
        r"$q_{0.1\%}$ = lower 0.1% percentile floor, so rare tails do not set threshold",
        fontsize=10.5,
        color=MUTED,
        va="top",
        zorder=4,
    )
    fig.text(
        0.535,
        0.055,
        r"Plotted energy:  $\log_{10}(E_{\rm CEMC}+E_{\rm IHCal}+E_{\rm OHCal}+1)$",
        fontsize=13.6,
        color=INK,
        ha="center",
        va="bottom",
    )
    fig.subplots_adjust(left=0.075, right=0.982, top=0.670, bottom=0.140, wspace=0.125, hspace=0.560)
    fig.savefig(PNG, dpi=200)
    fig.savefig(PHONE_PNG, dpi=200)
    plt.close(fig)

    MANIFEST.write_text(
        json.dumps(
            {
                "png": str(PNG),
                "phone_refresh_png": str(PHONE_PNG),
                "marker_source": str(MARKER_JSON),
                "cut_source": str(CUT_JSON),
                "histogram_x": "log10(total calorimeter energy + 1)",
                "normalization": "none",
                "smoothing": "none",
                "y_axis": "raw counts per log-energy bin",
                "centrality_bins": [f"{lo}-{hi}%" for lo, hi in centrality_bins],
                "cut_formula": (
                    "For each 5% centrality bin: y=log10(E_CEMC+E_IHCal+E_OHCal+1); "
                    "T_bin=max[median(y)-5*1.4826*MAD(y), q_0.1%(y)]; reject event if y<T_bin."
                ),
                "cut_formula_terms": (
                    "MAD is the median absolute deviation, scaled by 1.4826 as a stable width estimate. "
                    "q_0.1%(y) is the lower 0.1% percentile floor, used so rare tails do not set the threshold."
                ),
                "cut_derivation_note": (
                    "The threshold is a per-bin low-energy floor derived from the bin's "
                    "total-calo event-energy distribution."
                ),
                "display_note": "Blue is the retained marker sample; red is the full removed-event marker set.",
                "panel_counts": [
                    {
                        "centrality_bin": f"{panel['lo']}-{panel['hi']}%",
                        "retained_marker_count": int(panel["kept"].size),
                        "removed_event_count": int(panel["removed"].size),
                        "retained_peak_bin_count": int(panel["kept_counts"].max(initial=0)),
                        "removed_peak_bin_count": int(panel["removed_counts"].max(initial=0)),
                        "thresholds_drawn": [float(x) for x in panel["thresholds"]],
                    }
                    for panel in prepared
                ],
            },
            indent=2,
        )
        + "\n"
    )


def main() -> int:
    render()
    print(PNG)
    print(PHONE_PNG)
    print(MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
