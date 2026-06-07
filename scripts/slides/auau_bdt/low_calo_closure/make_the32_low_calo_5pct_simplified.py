#!/usr/bin/env python3
"""Make a simplified display-binned view of the low-calo event-veto result."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
MARKER_JSON = OUTDIR / "the32_low_calo_marker_sample_v1.json"
CUT_JSON = OUTDIR / "the32_low_calo_cut_v1.json"
PNG = OUTDIR / "the32_low_calo_10pct_display_data_only_20260606.png"
PHONE_PNG = OUTDIR / "the32_low_calo_10pct_display_data_only_phone_refresh_20260606.png"
MANIFEST = OUTDIR / "the32_low_calo_10pct_display_data_only_20260606.json"
DISPLAY_BIN_WIDTH = 10.0

BLUE = "#1f77b4"
RED = "#d62728"
INK = "#172033"
MUTED = "#5d6b7a"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def arrays(payload: dict, key: str) -> tuple[np.ndarray, np.ndarray]:
    markers = payload["markers"][key]
    return (
        np.asarray(markers["centrality"], dtype=np.float32),
        np.asarray(markers["log10_calo"], dtype=np.float32),
    )


def bin_display_x(x: np.ndarray, *, seed: int) -> np.ndarray:
    """Spread integer-like centrality values inside wider display bins."""
    rng = np.random.default_rng(seed)
    bin_lo = np.floor(np.clip(x, 0, 79.999) / DISPLAY_BIN_WIDTH) * DISPLAY_BIN_WIDTH
    return bin_lo + rng.uniform(0.55, DISPLAY_BIN_WIDTH - 0.55, size=x.size).astype(np.float32)


def sample_retained_by_bin(x: np.ndarray, y: np.ndarray, *, per_bin: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for lo in np.arange(0, 80, DISPLAY_BIN_WIDTH):
        idx = np.flatnonzero((x >= lo) & (x < lo + DISPLAY_BIN_WIDTH) & (y >= 2.76) & (y <= 3.42))
        if idx.size > per_bin:
            idx = rng.choice(idx, size=per_bin, replace=False)
            idx.sort()
        if idx.size:
            selected.append(idx)
    return np.concatenate(selected) if selected else np.array([], dtype=np.int64)


def removed_counts_by_bin(x: np.ndarray, envelope: list[dict]) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    labels: list[str] = []
    counts: list[int] = []
    totals: list[int] = []
    fractions: list[float] = []
    for lo in np.arange(0, 80, DISPLAY_BIN_WIDTH):
        hi = lo + DISPLAY_BIN_WIDTH
        removed = int(((x >= lo) & (x < hi)).sum())
        total = int(
            sum(
                int(row["n_events"])
                for row in envelope
                if float(row["cent_lo"]) >= lo and float(row["cent_hi"]) <= hi
            )
        )
        labels.append(f"{int(lo):02d}-{int(hi):02d}")
        counts.append(removed)
        totals.append(total)
        fractions.append(removed / total if total else np.nan)
    return labels, np.asarray(counts), np.asarray(totals), np.asarray(fractions)


def render() -> None:
    payload = load_json(MARKER_JSON)
    cut = load_json(CUT_JSON)
    envelope = cut["envelope"]
    kept_x, kept_y = arrays(payload, "retained")
    removed_x, removed_y = arrays(payload, "removed")
    labels, removed_counts, totals, removed_frac = removed_counts_by_bin(removed_x, envelope)

    kept_idx = sample_retained_by_bin(kept_x, kept_y, per_bin=2300, seed=7011)
    removed_mask = (removed_x >= 0) & (removed_x < 80) & (removed_y >= 2.76) & (removed_y <= 3.42)
    removed_idx = np.flatnonzero(removed_mask)

    kept_x_plot = bin_display_x(kept_x[kept_idx], seed=7012)
    removed_x_plot = bin_display_x(removed_x[removed_idx], seed=7013)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=(4.9, 1.15),
        left=0.075,
        right=0.965,
        top=0.800,
        bottom=0.115,
        hspace=0.105,
    )
    ax = fig.add_subplot(gs[0])
    rate_ax = fig.add_subplot(gs[1], sharex=ax)

    ax.scatter(
        kept_x_plot,
        kept_y[kept_idx],
        s=8,
        c=BLUE,
        alpha=0.16,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        removed_x_plot,
        removed_y[removed_idx],
        s=12,
        c=RED,
        alpha=0.52,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xlim(0, 80)
    ax.set_ylim(2.76, 3.42)
    ax.set_ylabel(r"$\log_{10}(E_{\rm calo}^{\rm total}+1)$", fontsize=15.0, color=INK)
    ax.tick_params(axis="both", labelsize=11.8, colors=INK, length=4.5, width=1.0)
    ax.tick_params(axis="x", labelbottom=False)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(1.05)
        spine.set_color(INK)

    label_box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.70, "pad": 1.2}
    ax.text(0.018, 0.970, "blue = not removed", transform=ax.transAxes, fontsize=14.5, color=BLUE, va="top", bbox=label_box)
    ax.text(0.275, 0.970, "red = removed", transform=ax.transAxes, fontsize=14.5, color=RED, va="top", bbox=label_box)

    centers = np.arange(DISPLAY_BIN_WIDTH / 2.0, 80, DISPLAY_BIN_WIDTH)
    rate_ax.bar(centers, 100.0 * removed_frac, width=8.5, color=RED, alpha=0.35, edgecolor=RED, linewidth=1.0)
    rate_ax.set_xlim(0, 80)
    rate_ax.set_ylim(0, max(9.5, 100.0 * np.nanmax(removed_frac) * 1.15))
    rate_ax.set_ylabel("removed\n%", fontsize=12.5, color=INK)
    rate_ax.set_xlabel("Centrality percentile, rebinned to 10% for display only", fontsize=13.6, color=INK)
    rate_ax.set_xticks(np.arange(0, 81, 10))
    rate_ax.tick_params(axis="both", labelsize=11.2, colors=INK, length=4.3, width=1.0)
    rate_ax.grid(False)
    for spine in rate_ax.spines.values():
        spine.set_linewidth(1.05)
        spine.set_color(INK)
    for i in [3, 4, 5]:
        rate_ax.text(
            centers[i],
            100.0 * removed_frac[i] + 0.35,
            f"{removed_counts[i]:,}",
            fontsize=10.2,
            color=RED if i < 5 else MUTED,
            ha="center",
            va="bottom",
            fontweight="bold" if i < 5 else "normal",
        )

    fig.text(
        0.055,
        0.940,
        "Simplified display-binned view of the total-calo event veto",
        fontsize=25.0,
        fontweight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.055,
        0.885,
        "No threshold line: colors show the cut result; centrality is rebinned to 10% only for display clarity.",
        fontsize=13.6,
        color=MUTED,
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
                "cut_source": str(CUT_JSON),
                "threshold_line_drawn": False,
                "x_display": "centrality values are randomized within 10% display bins; cut labels are unchanged",
                "removed_counts_by_display_bin": [
                    {
                        "centrality_bin": label,
                        "removed_events": int(count),
                        "audit_total_events": int(total),
                        "removed_fraction": float(frac),
                    }
                    for label, count, total, frac in zip(labels, removed_counts, totals, removed_frac)
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
