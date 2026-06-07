#!/usr/bin/env python3
"""Make line-histogram views of total-calo energy for the low-calo veto."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
MARKER_JSON = OUTDIR / "the32_low_calo_marker_sample_v1.json"
PNG = OUTDIR / "the32_low_calo_energy_line_histograms_10pct_audience_20260606.png"
PHONE_PNG = OUTDIR / "the32_low_calo_energy_line_histograms_10pct_audience_phone_refresh_20260606.png"
MANIFEST = OUTDIR / "the32_low_calo_energy_line_histograms_10pct_audience_20260606.json"

BLUE = "#1f77b4"
RED = "#d62728"
INK = "#172033"
MUTED = "#5d6b7a"
GRID = "#D0D5DD"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def load_payload() -> dict:
    return json.loads(MARKER_JSON.read_text())


def arrays(payload: dict, key: str) -> tuple[np.ndarray, np.ndarray]:
    markers = payload["markers"][key]
    return (
        np.asarray(markers["centrality"], dtype=np.float32),
        np.asarray(markers["log10_calo"], dtype=np.float32),
    )


def peak_shape_hist(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(values, bins=edges)
    if counts.sum() == 0:
        return np.zeros(edges.size - 1, dtype=np.float64)
    kernel = np.asarray([1, 4, 10, 16, 19, 16, 10, 4, 1], dtype=np.float64)
    kernel /= kernel.sum()
    smooth = np.convolve(counts.astype(np.float64), kernel, mode="same")
    peak = float(np.max(smooth))
    return smooth / peak if peak > 0 else smooth


def render() -> None:
    payload = load_payload()
    kept_x, kept_y = arrays(payload, "retained")
    removed_x, removed_y = arrays(payload, "removed")
    bins = [(lo, lo + 10) for lo in range(0, 80, 10)]
    hist_edges = np.linspace(2.70, 3.42, 73)
    centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])

    fig, axes = plt.subplots(2, 4, figsize=(12.8, 7.2), dpi=200, sharex=True, sharey=True)
    fig.patch.set_facecolor("white")
    axes_flat = axes.ravel()
    metadata: list[dict] = []

    prepared: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]] = []
    for lo, hi in bins:
        kept = kept_y[(kept_x >= lo) & (kept_x < hi) & np.isfinite(kept_y)]
        removed = removed_y[(removed_x >= lo) & (removed_x < hi) & np.isfinite(removed_y)]
        kept_hist = peak_shape_hist(kept, hist_edges)
        removed_hist = peak_shape_hist(removed, hist_edges)
        prepared.append((kept, removed, kept_hist, removed_hist, int(kept.size), int(removed.size)))

    for ax, (lo, hi), (kept, removed, kept_hist, removed_hist, kept_n, removed_n) in zip(axes_flat, bins, prepared):
        ax.set_facecolor("white")
        ax.plot(centers, kept_hist, color=BLUE, lw=2.7, alpha=0.95)
        if removed_n:
            ax.plot(centers, removed_hist, color=RED, lw=2.7, alpha=0.95)
        ax.fill_between(centers, kept_hist, color=BLUE, alpha=0.08)
        if removed_n:
            ax.fill_between(centers, removed_hist, color=RED, alpha=0.08)
        ax.set_xlim(hist_edges[0], hist_edges[-1])
        ax.set_ylim(0, 1.15)
        ax.set_title(f"{lo}-{hi}% centrality", fontsize=13.8, fontweight="bold", color=INK, pad=14)
        ax.text(
            0.50,
            1.015,
            f"{removed_n:,} removed",
            transform=ax.transAxes,
            fontsize=10.4,
            color=RED if removed_n else MUTED,
            fontweight="bold" if removed_n else "normal",
            ha="center",
            va="bottom",
            clip_on=False,
        )
        if hi == 50:
            pass
        ax.tick_params(axis="both", labelsize=10.2, colors=INK, length=4.0, width=0.9)
        ax.set_yticks([0.0, 0.5, 1.0])
        ax.set_xticks([2.8, 3.0, 3.2, 3.4])
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_linewidth(0.95)
            spine.set_color(INK)
        metadata.append(
            {
                "centrality_bin": f"{lo}-{hi}%",
                "retained_marker_count": kept_n,
                "removed_event_count": removed_n,
                "retained_y_median": float(np.median(kept)) if kept_n else None,
                "removed_y_median": float(np.median(removed)) if removed_n else None,
            }
        )

    axes[0, 0].set_ylabel("relative shape\n(same peak height)", fontsize=11.7, color=INK)
    axes[1, 0].set_ylabel("relative shape\n(same peak height)", fontsize=11.7, color=INK)
    for ax in axes[1, :]:
        ax.set_xlabel(r"$\log_{10}(\mathrm{total\ calorimeter\ energy}+1)$", fontsize=10.8, color=INK)

    fig.text(
        0.055,
        0.955,
        "Total calorimeter energy shapes by centrality",
        fontsize=25.0,
        fontweight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.055,
        0.908,
        "Curves are scaled to the same height to compare shape; red numbers show how many events were removed.",
        fontsize=12.4,
        color=MUTED,
        va="top",
    )
    fig.text(0.700, 0.955, "blue = not removed", fontsize=13.2, color=BLUE, va="top")
    fig.text(0.835, 0.955, "red = removed", fontsize=13.2, color=RED, va="top")
    fig.subplots_adjust(left=0.072, right=0.982, top=0.825, bottom=0.125, wspace=0.120, hspace=0.360)
    fig.savefig(PNG, dpi=200)
    fig.savefig(PHONE_PNG, dpi=200)
    plt.close(fig)

    MANIFEST.write_text(
        json.dumps(
            {
                "png": str(PNG),
                "phone_refresh_png": str(PHONE_PNG),
                "marker_source": str(MARKER_JSON),
                "histogram_x": "log10(E_calo_total + 1)",
                "normalization": "each color line is smoothed and scaled to peak height 1 within its centrality panel",
                "audience_note": "Red counts are printed above each panel to show prevalence because the line shapes are scaled to the same peak height.",
                "centrality_bins": [f"{lo}-{hi}%" for lo, hi in bins],
                "panel_metadata": metadata,
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
