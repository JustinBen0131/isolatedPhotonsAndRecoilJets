#!/usr/bin/env python3
"""Make a 0-80% centrality data-only view of the low-calo event veto."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np


OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
MARKER_JSON = OUTDIR / "the32_low_calo_marker_sample_v1.json"
REJECTION_AUDIT = OUTDIR / "audit" / "the32_low_calo_rejection_audit_v1.csv"
PNG = OUTDIR / "the32_low_calo_data_only_0to80_no_threshold_20260606.png"
PHONE_PNG = OUTDIR / "the32_low_calo_data_only_0to80_no_threshold_phone_refresh_20260606.png"
MANIFEST = OUTDIR / "the32_low_calo_data_only_0to80_no_threshold_20260606.json"

BLUE = "#1f77b4"
RED = "#d62728"
INK = "#172033"
MUTED = "#5d6b7a"
GRID = "#D0D5DD"
SOFT_BLUE = "#EEF6FF"
SOFT_RED = "#FFF4ED"
SOFT_GRAY = "#F8FAFC"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def load_payload() -> dict:
    return json.loads(MARKER_JSON.read_text())


def marker_arrays(payload: dict, key: str) -> tuple[np.ndarray, np.ndarray]:
    markers = payload["markers"][key]
    return (
        np.asarray(markers["centrality"], dtype=np.float32),
        np.asarray(markers["log10_calo"], dtype=np.float32),
    )


def sample_by_integer_cent(
    mask: np.ndarray,
    x: np.ndarray,
    *,
    per_bin: int,
    seed: int,
    xlim: tuple[float, float],
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for lo in np.arange(xlim[0], xlim[1], 1):
        idx = np.flatnonzero(mask & (x >= lo) & (x < lo + 1))
        if idx.size > per_bin:
            idx = rng.choice(idx, size=per_bin, replace=False)
            idx.sort()
        if idx.size:
            selected.append(idx)
    return np.concatenate(selected) if selected else np.array([], dtype=np.int64)


def jitter(x: np.ndarray, *, seed: int, width: float = 0.42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return x + rng.uniform(-width, width, size=x.size).astype(np.float32)


def region_rejection_summary() -> list[tuple[str, int, int, float]]:
    buckets = {
        "0-20%": [0, 0],
        "20-50%": [0, 0],
        "50-80%": [0, 0],
    }
    with REJECTION_AUDIT.open() as handle:
        for row in csv.DictReader(handle):
            label = row["centrality_bin"]
            if label not in buckets:
                continue
            buckets[label][0] += int(row["event_rejected"])
            buckets[label][1] += int(row["event_total"])
    return [
        (label, values[0], values[1], values[0] / values[1] if values[1] else float("nan"))
        for label, values in buckets.items()
    ]


def add_box(ax, x: float, y: float, w: float, h: float, facecolor: str) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            linewidth=1.15,
            edgecolor=GRID,
            facecolor=facecolor,
        )
    )


def render() -> None:
    payload = load_payload()
    kept_x, kept_y = marker_arrays(payload, "retained")
    cut_x, cut_y = marker_arrays(payload, "removed")
    counts = payload["event_counts"]
    total_events = int(counts["event_total_after_per_cache_dedup"])
    removed_events = int(counts["event_removed"])
    removed_pct = 100.0 * removed_events / total_events

    xlim = (0.0, 80.0)
    ylim = (2.70, 3.43)
    kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
    cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])
    kept_idx = sample_by_integer_cent(kept_mask, kept_x, per_bin=95, seed=6601, xlim=xlim)
    cut_idx = sample_by_integer_cent(cut_mask, cut_x, per_bin=150, seed=6602, xlim=xlim)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.patch.set_facecolor("white")
    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.axis("off")

    canvas.text(
        0.055,
        0.940,
        "Data-only view of the total-calo event veto",
        fontsize=27.0,
        fontweight="bold",
        color=INK,
        va="top",
    )
    canvas.text(
        0.055,
        0.878,
        "No threshold line is drawn here: color alone shows the applied event-level cut result.",
        fontsize=14.4,
        color=MUTED,
        va="top",
    )

    ax = fig.add_axes([0.060, 0.155, 0.610, 0.665])
    ax.set_facecolor("white")
    ax.scatter(
        jitter(kept_x[kept_idx], seed=6603),
        kept_y[kept_idx],
        s=10,
        c=BLUE,
        alpha=0.20,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        jitter(cut_x[cut_idx], seed=6604),
        cut_y[cut_idx],
        s=16,
        c=RED,
        alpha=0.62,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile", fontsize=16.2, color=INK)
    ax.set_ylabel(r"$\log_{10}(E_{\rm calo}^{\rm total}+1)$", fontsize=16.2, color=INK)
    ax.tick_params(axis="both", labelsize=12.5, colors=INK, length=4.8, width=1.05)
    for spine in ax.spines.values():
        spine.set_linewidth(1.08)
        spine.set_color(INK)
    ax.grid(False)
    label_box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.4}
    ax.text(0.020, 0.967, "blue = not removed", transform=ax.transAxes, color=BLUE, fontsize=14.4, va="top", bbox=label_box)
    ax.text(0.305, 0.967, "red = removed", transform=ax.transAxes, color=RED, fontsize=14.4, va="top", bbox=label_box)

    add_box(canvas, 0.705, 0.625, 0.245, 0.205, SOFT_BLUE)
    canvas.text(0.725, 0.795, "What the colors mean", fontsize=16.8, fontweight="bold", color=INK, va="top")
    canvas.text(0.725, 0.748, "Blue events passed the total-calo veto.", fontsize=12.2, color=BLUE, va="top")
    canvas.text(0.725, 0.710, "Red events failed it and were removed.", fontsize=12.2, color=RED, va="top")
    canvas.text(0.725, 0.670, "This plot intentionally hides the threshold line.", fontsize=11.6, color=MUTED, va="top")

    add_box(canvas, 0.705, 0.355, 0.245, 0.225, SOFT_RED)
    canvas.text(0.725, 0.545, "Why red fades near 45%", fontsize=16.0, fontweight="bold", color=INK, va="top")
    canvas.text(0.725, 0.503, "Blue band falls with centrality.", fontsize=11.7, color=INK, va="top")
    canvas.text(0.725, 0.466, "Low-calo band stays near y ~ 2.84.", fontsize=11.7, color=INK, va="top")
    canvas.text(0.725, 0.429, "Near 40-45%, the veto boundary reaches", fontsize=11.2, color=INK, va="top")
    canvas.text(0.725, 0.397, "that low band. Past that, only deeper", fontsize=11.2, color=INK, va="top")
    canvas.text(0.725, 0.365, "low-energy tails are still red.", fontsize=11.2, color=INK, va="top")

    add_box(canvas, 0.705, 0.075, 0.245, 0.235, SOFT_GRAY)
    canvas.text(0.725, 0.275, "Regional removal rate", fontsize=16.0, fontweight="bold", color=INK, va="top")
    y0 = 0.232
    for i, (label, removed, total, frac) in enumerate(region_rejection_summary()):
        canvas.text(
            0.725,
            y0 - i * 0.044,
            f"{label}: {removed:,}/{total:,} events ({100.0 * frac:.2f}%)",
            fontsize=11.4,
            color=RED if i < 2 else MUTED,
            fontweight="bold" if i < 2 else "normal",
            va="top",
        )
    canvas.text(
        0.725,
        0.100,
        f"Total: {removed_events:,}/{total_events:,} events ({removed_pct:.1f}%)",
        fontsize=11.5,
        color=INK,
        fontweight="bold",
        va="top",
    )

    fig.savefig(PNG, dpi=200)
    fig.savefig(PHONE_PNG, dpi=200)
    plt.close(fig)
    MANIFEST.write_text(
        json.dumps(
            {
                "png": str(PNG),
                "phone_refresh_png": str(PHONE_PNG),
                "marker_source": str(MARKER_JSON),
                "threshold_line_drawn": False,
                "color_semantics": {
                    "blue": "event was not removed by the event-level total-calo veto",
                    "red": "event was removed by the event-level total-calo veto",
                },
                "centrality_range": "0-80%",
                "y_range": list(ylim),
                "regional_rejection_summary": [
                    {
                        "centrality": label,
                        "removed_events": removed,
                        "total_events": total,
                        "removed_fraction": frac,
                    }
                    for label, removed, total, frac in region_rejection_summary()
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
