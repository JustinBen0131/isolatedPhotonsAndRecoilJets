#!/usr/bin/env python3
"""Make marker-only THE-32 low-calo data views from raw score-cache rows."""

from __future__ import annotations

import json
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
MARKER_JSON = OUTDIR / "the32_low_calo_marker_sample_v1.json"
CUT_JSON = OUTDIR / "the32_low_calo_cut_v1.json"
REJECTION_AUDIT = OUTDIR / "audit" / "the32_low_calo_rejection_audit_v1.csv"

BLUE = "#1f77b4"
RED = "#d62728"
CUT = "#111827"
CUT_HIGHLIGHT = "#f59e0b"
INK = "#1b2638"
MUTED = "#5d6b7a"
CALO_LOG_LABEL = r"$\log_{10}(E_{\rm CEMC} + E_{\rm IHCal} + E_{\rm OHCal} + 1)$"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def load_markers() -> dict:
    return json.loads(MARKER_JSON.read_text())


def load_cut_envelope() -> list[dict]:
    return json.loads(CUT_JSON.read_text())["envelope"]


def load_region_rejection_fractions() -> dict[str, tuple[int, int, float]]:
    totals = {
        "0-20%": [0, 0],
        "20-50%": [0, 0],
        "50-80%": [0, 0],
    }
    if not REJECTION_AUDIT.exists():
        return {}
    with REJECTION_AUDIT.open() as handle:
        for row in csv.DictReader(handle):
            label = row["centrality_bin"]
            if label not in totals:
                continue
            totals[label][0] += int(row["event_rejected"])
            totals[label][1] += int(row["event_total"])
    return {
        label: (rejected, total, rejected / total if total else float("nan"))
        for label, (rejected, total) in totals.items()
    }


def arrays(payload: dict, kind: str) -> tuple[np.ndarray, np.ndarray]:
    markers = payload["markers"][kind]
    return (
        np.asarray(markers["centrality"], dtype=np.float32),
        np.asarray(markers["log10_calo"], dtype=np.float32),
    )


def deterministic_jitter(x: np.ndarray, *, width: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return x + rng.uniform(-width, width, size=x.size).astype(np.float32)


def stratified_sample(
    mask: np.ndarray,
    x: np.ndarray,
    *,
    target_per_bin: int,
    seed: int,
    x_hi: float = 35.0,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for lo in np.arange(0, x_hi, 1):
        in_bin = np.flatnonzero(mask & (x >= lo) & (x < lo + 1))
        if in_bin.size > target_per_bin:
            in_bin = rng.choice(in_bin, size=target_per_bin, replace=False)
            in_bin.sort()
        if in_bin.size:
            selected.append(in_bin)
    if not selected:
        return np.array([], dtype=np.int64)
    return np.concatenate(selected)


def cut_step_arrays(envelope: list[dict], *, xlim: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    for row in envelope:
        lo = max(float(row["cent_lo"]), xlim[0])
        hi = min(float(row["cent_hi"]), xlim[1])
        if hi <= xlim[0] or lo >= xlim[1]:
            continue
        xs.extend([lo, hi])
        ys.extend([float(row["threshold"]), float(row["threshold"])])
    return np.asarray(xs, dtype=np.float32), np.asarray(ys, dtype=np.float32)


def x_jitter_for_range(x: np.ndarray, *, seed: int, width: float = 0.36) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return x + rng.uniform(-width, width, size=x.size).astype(np.float32)


def region_sample(mask: np.ndarray, x: np.ndarray, *, xlim: tuple[float, float], target_per_bin: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for lo in np.arange(xlim[0], xlim[1], 1):
        in_bin = np.flatnonzero(mask & (x >= lo) & (x < lo + 1))
        if in_bin.size > target_per_bin:
            in_bin = rng.choice(in_bin, size=target_per_bin, replace=False)
            in_bin.sort()
        if in_bin.size:
            selected.append(in_bin)
    if not selected:
        return np.array([], dtype=np.int64)
    return np.concatenate(selected)


def qualitative_region_panel(payload: dict, out: Path) -> None:
    envelope = load_cut_envelope()
    fractions = load_region_rejection_fractions()
    kept_x, kept_y = arrays(payload, "retained")
    cut_x, cut_y = arrays(payload, "removed")
    ylim = (2.36, 3.36)
    regions = [
        ("0-20%", (0.0, 20.0), "central", 220, 320, 4101),
        ("20-50%", (20.0, 50.0), "mid-central", 150, 260, 4201),
        ("50-80%", (50.0, 80.0), "peripheral", 150, 260, 4301),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 7.8), dpi=200, sharey=True)
    fig.patch.set_facecolor("white")
    for ax, (label, xlim, short_label, kept_cap, cut_cap, seed) in zip(axes, regions):
        ax.set_facecolor("white")
        kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
        cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])
        kept_idx = region_sample(kept_mask, kept_x, xlim=xlim, target_per_bin=kept_cap, seed=seed)
        cut_idx = region_sample(cut_mask, cut_x, xlim=xlim, target_per_bin=cut_cap, seed=seed + 1)
        ax.scatter(
            x_jitter_for_range(kept_x[kept_idx], seed=seed + 2),
            kept_y[kept_idx],
            s=9,
            c=BLUE,
            alpha=0.18,
            linewidths=0,
            rasterized=True,
        )
        ax.scatter(
            x_jitter_for_range(cut_x[cut_idx], seed=seed + 3),
            cut_y[cut_idx],
            s=14,
            c=RED,
            alpha=0.58,
            linewidths=0,
            rasterized=True,
        )
        step_x, step_y = cut_step_arrays(envelope, xlim=xlim)
        ax.plot(step_x, step_y, color="white", lw=6.2, solid_capstyle="butt", zorder=8)
        ax.plot(step_x, step_y, color=CUT, lw=3.1, solid_capstyle="butt", zorder=9)
        ax.plot(step_x, step_y, color=CUT_HIGHLIGHT, lw=1.25, solid_capstyle="butt", zorder=10)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_title(f"{label} {short_label}", fontsize=18, weight="bold", color=INK, pad=12)
        ax.set_xlabel("Centrality percentile", fontsize=14.5, color=INK)
        if label in fractions:
            rejected, total, frac = fractions[label]
            frac_text = f"{100.0 * frac:.2f}% removed"
            count_text = f"{rejected:,} of {total:,} events"
        else:
            frac_text = "audit fraction unavailable"
            count_text = ""
        ax.text(
            0.04,
            0.055,
            frac_text + ("\n" + count_text if count_text else ""),
            transform=ax.transAxes,
            fontsize=12.2,
            color=RED if label != "50-80%" else MUTED,
            weight="bold" if label != "50-80%" else "normal",
            va="bottom",
            ha="left",
            bbox=dict(facecolor="white", edgecolor="#d0d5dd", boxstyle="round,pad=0.28", alpha=0.92),
        )
        ax.tick_params(axis="both", labelsize=12.2, colors=INK, length=4.5, width=1.0)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_linewidth(1.05)
            spine.set_color(INK)

    axes[0].set_ylabel(CALO_LOG_LABEL, fontsize=15.5, color=INK)
    fig.text(
        0.055,
        0.965,
        "Low-total-calo events are a central/mid-central excess",
        fontsize=25.5,
        weight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.055,
        0.918,
        "Same y-scale in all panels: the red excess is dense in central/mid-central events and becomes rare in 50-80%.",
        fontsize=13.8,
        color=MUTED,
        va="top",
    )
    fig.text(0.055, 0.878, "blue = passes cut", fontsize=12.6, color=BLUE, va="top")
    fig.text(0.205, 0.878, "red = fails cut", fontsize=12.6, color=RED, va="top")
    fig.text(0.345, 0.878, "black/gold = applied threshold", fontsize=12.6, color=INK, va="top")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.79, bottom=0.12, wspace=0.10)
    fig.savefig(out, dpi=200)
    plt.close(fig)


def audience_clean_panel(payload: dict, out: Path) -> None:
    kept_x, kept_y = arrays(payload, "retained")
    cut_x, cut_y = arrays(payload, "removed")
    xlim = (0.0, 35.0)
    ylim = (2.80, 3.36)
    kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
    cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])

    # The x-axis is integer-like in the source data. This jitter is display-only:
    # it prevents overplotting from reading like drawn vertical guide lines.
    kept_idx = stratified_sample(kept_mask, kept_x, target_per_bin=220, seed=3201)
    cut_idx = stratified_sample(cut_mask, cut_x, target_per_bin=260, seed=3202)
    kept_x_plot = deterministic_jitter(kept_x[kept_idx], width=0.50, seed=3203)
    cut_x_plot = deterministic_jitter(cut_x[cut_idx], width=0.50, seed=3204)

    fig, ax = plt.subplots(figsize=(16, 9), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.scatter(
        kept_x_plot,
        kept_y[kept_idx],
        s=10,
        c=BLUE,
        alpha=0.17,
        linewidths=0,
        label="retained event sample",
        rasterized=True,
    )
    ax.scatter(
        cut_x_plot,
        cut_y[cut_idx],
        s=17,
        c=RED,
        alpha=0.62,
        linewidths=0,
        label="removed low-calo event sample",
        rasterized=True,
    )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile", fontsize=18, color=INK)
    ax.set_ylabel(CALO_LOG_LABEL, fontsize=18, color=INK)
    ax.tick_params(axis="both", labelsize=14, colors=INK, length=5, width=1.1)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color(INK)
    ax.grid(False)
    ax.legend(
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.92,
        fontsize=16,
        markerscale=1.8,
        handletextpad=0.5,
        borderaxespad=0.9,
    )
    fig.text(0.08, 0.94, "Low-calo events separate below the retained data", ha="left", va="top", fontsize=25, weight="bold", color=INK)
    fig.subplots_adjust(left=0.095, right=0.975, top=0.875, bottom=0.12)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def audience_binned_panel(payload: dict, out: Path) -> None:
    kept_x, kept_y = arrays(payload, "retained")
    cut_x, cut_y = arrays(payload, "removed")
    xlim = (0.0, 35.0)
    ylim = (2.80, 3.36)
    kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
    cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])

    rng_kept = np.random.default_rng(3211)
    rng_cut = np.random.default_rng(3212)
    kept_parts: list[np.ndarray] = []
    cut_parts: list[np.ndarray] = []
    for lo in np.arange(0, 35, 5):
        kept_bin = np.flatnonzero(kept_mask & (kept_x >= lo) & (kept_x < lo + 5))
        cut_bin = np.flatnonzero(cut_mask & (cut_x >= lo) & (cut_x < lo + 5))
        if kept_bin.size > 700:
            kept_bin = rng_kept.choice(kept_bin, size=700, replace=False)
            kept_bin.sort()
        if cut_bin.size > 700:
            cut_bin = rng_cut.choice(cut_bin, size=700, replace=False)
            cut_bin.sort()
        kept_parts.append(kept_bin)
        cut_parts.append(cut_bin)
    kept_idx = np.concatenate([x for x in kept_parts if x.size])
    cut_idx = np.concatenate([x for x in cut_parts if x.size])

    # Display each integer centrality measurement inside its enclosing 5% bin so
    # the plot reads as event populations, not barcode-like integer columns.
    kept_bin_lo = np.floor(kept_x[kept_idx] / 5.0) * 5.0
    cut_bin_lo = np.floor(cut_x[cut_idx] / 5.0) * 5.0
    kept_x_plot = kept_bin_lo + rng_kept.uniform(0.35, 4.65, size=kept_idx.size)
    cut_x_plot = cut_bin_lo + rng_cut.uniform(0.35, 4.65, size=cut_idx.size)

    fig, ax = plt.subplots(figsize=(16, 9), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.scatter(
        kept_x_plot,
        kept_y[kept_idx],
        s=11,
        c=BLUE,
        alpha=0.18,
        linewidths=0,
        label="retained event sample",
        rasterized=True,
    )
    ax.scatter(
        cut_x_plot,
        cut_y[cut_idx],
        s=18,
        c=RED,
        alpha=0.60,
        linewidths=0,
        label="removed low-calo event sample",
        rasterized=True,
    )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile, shown in 5% display bins", fontsize=18, color=INK)
    ax.set_ylabel(CALO_LOG_LABEL, fontsize=18, color=INK)
    ax.set_xticks(np.arange(0, 36, 5))
    ax.tick_params(axis="both", labelsize=14, colors=INK, length=5, width=1.1)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color(INK)
    ax.grid(False)
    ax.legend(
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.92,
        fontsize=16,
        markerscale=1.8,
        handletextpad=0.5,
        borderaxespad=0.9,
    )
    fig.text(0.08, 0.94, "Low-calo events form a separated population", ha="left", va="top", fontsize=25, weight="bold", color=INK)
    fig.subplots_adjust(left=0.095, right=0.975, top=0.875, bottom=0.12)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def audience_matched_panel(payload: dict, out: Path) -> None:
    kept_x, kept_y = arrays(payload, "retained")
    cut_x, cut_y = arrays(payload, "removed")
    xlim = (0.0, 35.0)
    ylim = (2.80, 3.36)
    kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
    cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])
    kept_idx = stratified_sample(kept_mask, kept_x, target_per_bin=150, seed=3221)
    cut_idx = stratified_sample(cut_mask, cut_x, target_per_bin=150, seed=3222)
    kept_x_plot = deterministic_jitter(kept_x[kept_idx], width=0.50, seed=3223)
    cut_x_plot = deterministic_jitter(cut_x[cut_idx], width=0.50, seed=3224)

    fig, ax = plt.subplots(figsize=(16, 9), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.scatter(
        kept_x_plot,
        kept_y[kept_idx],
        s=12,
        c=BLUE,
        alpha=0.23,
        linewidths=0,
        label="retained event sample",
        rasterized=True,
    )
    ax.scatter(
        cut_x_plot,
        cut_y[cut_idx],
        s=16,
        c=RED,
        alpha=0.58,
        linewidths=0,
        label="removed low-calo event sample",
        rasterized=True,
    )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile", fontsize=18, color=INK)
    ax.set_ylabel(CALO_LOG_LABEL, fontsize=18, color=INK)
    ax.tick_params(axis="both", labelsize=14, colors=INK, length=5, width=1.1)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color(INK)
    ax.grid(False)
    ax.legend(
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.92,
        fontsize=16,
        markerscale=1.8,
        handletextpad=0.5,
        borderaxespad=0.9,
    )
    fig.text(0.08, 0.94, "Matched marker samples: removed events sit below retained events", ha="left", va="top", fontsize=25, weight="bold", color=INK)
    fig.subplots_adjust(left=0.095, right=0.975, top=0.875, bottom=0.12)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def blair_mattermost_panel(payload: dict, out: Path) -> None:
    kept_x, kept_y = arrays(payload, "retained")
    cut_x, cut_y = arrays(payload, "removed")
    xlim = (0.0, 35.0)
    ylim = (2.80, 3.36)
    kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
    cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])
    kept_idx = stratified_sample(kept_mask, kept_x, target_per_bin=145, seed=3231)
    cut_idx = stratified_sample(cut_mask, cut_x, target_per_bin=145, seed=3232)
    kept_x_plot = deterministic_jitter(kept_x[kept_idx], width=0.50, seed=3233)
    cut_x_plot = deterministic_jitter(cut_x[cut_idx], width=0.50, seed=3234)

    fig, ax = plt.subplots(figsize=(14.5, 8.4), dpi=190)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.scatter(
        kept_x_plot,
        kept_y[kept_idx],
        s=11,
        c=BLUE,
        alpha=0.24,
        linewidths=0,
        label="retained event sample",
        rasterized=True,
    )
    ax.scatter(
        cut_x_plot,
        cut_y[cut_idx],
        s=15,
        c=RED,
        alpha=0.60,
        linewidths=0,
        label="low-calo events removed by cut",
        rasterized=True,
    )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile", fontsize=17, color=INK)
    ax.set_ylabel(CALO_LOG_LABEL, fontsize=17, color=INK)
    ax.tick_params(axis="both", labelsize=13, colors=INK, length=5, width=1.1)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color(INK)
    ax.grid(False)
    ax.legend(
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.94,
        fontsize=14,
        markerscale=1.8,
        handletextpad=0.5,
        borderaxespad=0.9,
    )
    fig.text(0.08, 0.945, "Low-calo events form a separated population", ha="left", va="top", fontsize=23, weight="bold", color=INK)
    fig.text(
        0.08,
        0.902,
        "0-35% centrality; matched marker samples and display-only x-jitter are used only to avoid overplotting.",
        ha="left",
        va="top",
        fontsize=12.5,
        color=MUTED,
    )
    fig.subplots_adjust(left=0.105, right=0.975, top=0.84, bottom=0.13)
    fig.savefig(out, dpi=190)
    plt.close(fig)


def blair_canvas_panel(
    payload: dict,
    out: Path,
    *,
    x_hi: float = 35.0,
    ylim: tuple[float, float] = (2.80, 3.36),
) -> None:
    envelope = load_cut_envelope()
    kept_x, kept_y = arrays(payload, "retained")
    cut_x, cut_y = arrays(payload, "removed")
    xlim = (0.0, float(x_hi))
    kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
    cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])
    kept_idx = stratified_sample(kept_mask, kept_x, target_per_bin=135, seed=3241, x_hi=x_hi)
    cut_idx = stratified_sample(cut_mask, cut_x, target_per_bin=135, seed=3242, x_hi=x_hi)
    kept_x_plot = deterministic_jitter(kept_x[kept_idx], width=0.50, seed=3243)
    cut_x_plot = deterministic_jitter(cut_x[cut_idx], width=0.50, seed=3244)

    total = int(payload["event_counts"]["event_total_after_per_cache_dedup"])
    removed = int(payload["event_counts"]["event_removed"])
    removed_frac = 100.0 * removed / total
    kept_med = float(np.median(kept_y[kept_mask]))
    cut_med = float(np.median(cut_y[cut_mask]))

    fig = plt.figure(figsize=(16, 9), dpi=200)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes((0.070, 0.135, 0.600, 0.735))
    ax.set_facecolor("white")
    ax.scatter(
        kept_x_plot,
        kept_y[kept_idx],
        s=11,
        c=BLUE,
        alpha=0.24,
        linewidths=0,
        label="passes cut",
        rasterized=True,
    )
    ax.scatter(
        cut_x_plot,
        cut_y[cut_idx],
        s=15,
        c=RED,
        alpha=0.58,
        linewidths=0,
        label="fails cut",
        rasterized=True,
    )
    step_x, step_y = cut_step_arrays(envelope, xlim=xlim)
    ax.plot(step_x, step_y, color="white", lw=7.0, solid_capstyle="butt", zorder=8)
    ax.plot(step_x, step_y, color=CUT, lw=3.5, solid_capstyle="butt", zorder=9)
    ax.plot(step_x, step_y, color=CUT_HIGHLIGHT, lw=1.4, solid_capstyle="butt", zorder=10)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile", fontsize=18, color=INK)
    ax.set_ylabel(CALO_LOG_LABEL, fontsize=17, color=INK)
    ax.tick_params(axis="both", labelsize=14, colors=INK, length=5, width=1.1)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color(INK)
    ax.grid(False)
    fig.text(0.070, 0.955, "Per-event total calorimeter energy threshold on merged embedding", ha="left", va="top", fontsize=30, weight="bold", color=INK)

    info_ax = fig.add_axes((0.705, 0.135, 0.255, 0.735))
    info_ax.set_facecolor("white")
    for spine in info_ax.spines.values():
        spine.set_edgecolor("#cbd5e1")
        spine.set_linewidth(1.1)
    info_ax.set_xticks([])
    info_ax.set_yticks([])
    info_ax.set_xlim(0, 1)
    info_ax.set_ylim(0, 1)
    def separator(y: float) -> None:
        info_ax.plot([0.06, 0.94], [y, y], color="#d8e0ea", lw=1.0)

    info_ax.text(0.06, 0.940, "Data shown", fontsize=17.5, weight="bold", color=INK)
    info_ax.text(0.06, 0.890, "Unweighted event-level points", fontsize=12.8, color=INK)
    info_ax.text(0.06, 0.852, "from merged embedding score caches", fontsize=12.2, color=INK)
    info_ax.text(0.09, 0.806, "Photon12, Photon20", fontsize=12.2, color=INK)
    info_ax.text(0.09, 0.768, "Jet12, Jet20, Jet30, Jet40", fontsize=12.2, color=INK)
    separator(0.718)

    info_ax.text(0.06, 0.668, "Y-axis", fontsize=17.5, weight="bold", color=INK)
    info_ax.text(0.06, 0.620, "Per-event total calorimeter energy", fontsize=12.7, color=INK)
    info_ax.text(0.06, 0.580, r"$y=\log_{10}(E_{\rm calo}+1)$", fontsize=13.8, color=INK)
    info_ax.text(0.06, 0.538, r"$E_{\rm calo}=E_{\rm CEMC}+E_{\rm IHCal}+E_{\rm OHCal}$", fontsize=13.0, color=INK)
    info_ax.text(0.06, 0.498, "not photon-cluster energy", fontsize=12.1, color=MUTED)
    separator(0.448)

    info_ax.text(0.06, 0.398, "How to read it", fontsize=17.5, weight="bold", color=INK)
    info_ax.text(0.06, 0.350, "blue = retained event band", fontsize=12.7, color=BLUE)
    info_ax.text(0.06, 0.310, "red = low-energy excess removed", fontsize=12.7, color=RED)
    info_ax.plot([0.07, 0.25], [0.266, 0.266], color=CUT, lw=3.4, solid_capstyle="butt")
    info_ax.plot([0.07, 0.25], [0.266, 0.266], color=CUT_HIGHLIGHT, lw=1.4, solid_capstyle="butt")
    info_ax.text(0.29, 0.266, "applied threshold", va="center", fontsize=12.2, color=INK)
    separator(0.212)

    info_ax.text(0.06, 0.162, "Filter effect", fontsize=17.5, weight="bold", color=INK)
    info_ax.text(0.06, 0.108, f"{removed:,} events removed", fontsize=16.2, weight="bold", color=RED)
    info_ax.text(0.06, 0.064, f"{removed_frac:.1f}% of {total:,} diagnostic events", fontsize=12.4, color=INK)
    info_ax.text(0.06, 0.026, f"0-{int(x_hi)}% median y: red {cut_med:.3f}, blue {kept_med:.3f}", fontsize=10.4, color=MUTED)

    fig.savefig(out, dpi=200)
    plt.close(fig)


def smooth_threshold(cent: np.ndarray) -> np.ndarray:
    cent = np.asarray(cent, dtype=np.float64)
    threshold = np.full(cent.shape, np.nan, dtype=np.float64)
    in_range = (cent >= 0.0) & (cent < 80.0)
    central = in_range & (cent < 50.0)
    peripheral = in_range & (cent >= 50.0)
    threshold[central] = 3.1875 - 0.01226 * cent[central] + 0.0000983 * cent[central] ** 2
    threshold[peripheral] = 2.8233
    return threshold


def smooth_candidate_panel(payload: dict, out: Path) -> None:
    envelope = load_cut_envelope()
    kept_x, kept_y = arrays(payload, "retained")
    cut_x, cut_y = arrays(payload, "removed")
    all_x = np.concatenate([kept_x, cut_x]).astype(np.float32, copy=False)
    all_y = np.concatenate([kept_y, cut_y]).astype(np.float32, copy=False)
    validated_fail = np.concatenate([
        np.zeros(len(kept_x), dtype=bool),
        np.ones(len(cut_x), dtype=bool),
    ])
    xlim = (0.0, 80.0)
    ylim = (2.35, 3.36)
    threshold = smooth_threshold(all_x)
    smooth_fail = np.isfinite(threshold) & (all_y < threshold)
    in_view = (all_x >= xlim[0]) & (all_x <= xlim[1]) & (all_y >= ylim[0]) & (all_y <= ylim[1])
    pass_idx = stratified_sample(in_view & ~smooth_fail, all_x, target_per_bin=120, seed=3251, x_hi=80.0)
    fail_idx = stratified_sample(in_view & smooth_fail, all_x, target_per_bin=120, seed=3252, x_hi=80.0)

    mid = np.asarray([(row["cent_lo"] + row["cent_hi"]) / 2.0 for row in envelope], dtype=np.float64)
    actual = np.asarray([row["threshold"] for row in envelope], dtype=np.float64)
    approx = smooth_threshold(mid)
    max_diff = float(np.nanmax(np.abs(approx - actual)))
    displayed_agreement = float(np.mean(smooth_fail[in_view] == validated_fail[in_view]))
    fail_frac = 100.0 * float(np.mean(smooth_fail[in_view]))

    fig = plt.figure(figsize=(16, 9), dpi=200)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes((0.070, 0.135, 0.600, 0.735))
    ax.set_facecolor("white")
    ax.scatter(
        x_jitter_for_range(all_x[pass_idx], seed=3253, width=0.42),
        all_y[pass_idx],
        s=10,
        c=BLUE,
        alpha=0.22,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        x_jitter_for_range(all_x[fail_idx], seed=3254, width=0.42),
        all_y[fail_idx],
        s=15,
        c=RED,
        alpha=0.60,
        linewidths=0,
        rasterized=True,
    )
    curve_x = np.linspace(0.0, 80.0, 500)
    curve_y = smooth_threshold(curve_x)
    ax.plot(curve_x, curve_y, color="white", lw=7.0, zorder=8)
    ax.plot(curve_x, curve_y, color=CUT, lw=3.4, zorder=9)
    ax.plot(curve_x, curve_y, color=CUT_HIGHLIGHT, lw=1.35, zorder=10)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile", fontsize=18, color=INK)
    ax.set_ylabel(CALO_LOG_LABEL, fontsize=17, color=INK)
    ax.tick_params(axis="both", labelsize=14, colors=INK, length=5, width=1.1)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color(INK)
    ax.grid(False)
    fig.text(0.070, 0.955, "Smooth event-energy threshold candidate", ha="left", va="top", fontsize=30, weight="bold", color=INK)
    fig.text(
        0.070,
        0.910,
        "Same event variable; smooth curve replaces the 5%-bin step table for easier communication.",
        ha="left",
        va="top",
        fontsize=13.8,
        color=MUTED,
    )

    info_ax = fig.add_axes((0.705, 0.135, 0.255, 0.735))
    info_ax.set_facecolor("white")
    for spine in info_ax.spines.values():
        spine.set_edgecolor("#cbd5e1")
        spine.set_linewidth(1.1)
    info_ax.set_xticks([])
    info_ax.set_yticks([])
    info_ax.set_xlim(0, 1)
    info_ax.set_ylim(0, 1)

    def separator(y: float) -> None:
        info_ax.plot([0.06, 0.94], [y, y], color="#d8e0ea", lw=1.0)

    info_ax.text(0.06, 0.940, "Smooth candidate", fontsize=17.5, weight="bold", color=INK)
    info_ax.text(0.06, 0.890, r"$y=\log_{10}(E_{\rm calo}+1)$", fontsize=14.8, color=INK)
    info_ax.text(0.06, 0.842, r"$E_{\rm calo}=E_{\rm CEMC}+E_{\rm IHCal}+E_{\rm OHCal}$", fontsize=12.5, color=INK)
    info_ax.text(0.06, 0.798, "Remove event if y < T(C)", fontsize=13.3, color=RED, weight="bold")
    separator(0.750)

    info_ax.text(0.06, 0.700, "Function shown", fontsize=17.5, weight="bold", color=INK)
    info_ax.text(0.06, 0.650, r"$T(C)=3.1875-0.01226C$", fontsize=13.0, color=INK)
    info_ax.text(0.06, 0.610, r"$\quad +\,9.83\times10^{-5}C^2$", fontsize=13.0, color=INK)
    info_ax.text(0.06, 0.566, r"for $0\leq C<50$", fontsize=12.6, color=INK)
    info_ax.text(0.06, 0.518, r"$T(C)=2.8233$ for $50\leq C<80$", fontsize=12.6, color=INK)
    info_ax.text(0.06, 0.472, "C = centrality percentile", fontsize=11.9, color=MUTED)
    separator(0.420)

    info_ax.text(0.06, 0.370, "How it compares", fontsize=17.5, weight="bold", color=INK)
    info_ax.text(0.06, 0.320, f"max table difference: {max_diff:.3f}", fontsize=12.8, color=INK)
    info_ax.text(0.06, 0.278, f"displayed-label agreement: {100.0 * displayed_agreement:.1f}%", fontsize=12.8, color=INK)
    info_ax.text(0.06, 0.236, f"displayed smooth-fail fraction: {fail_frac:.1f}%", fontsize=12.8, color=INK)
    separator(0.188)

    info_ax.text(0.06, 0.138, "Plot colors", fontsize=17.5, weight="bold", color=INK)
    info_ax.text(0.06, 0.094, "blue = passes smooth candidate", fontsize=12.3, color=BLUE)
    info_ax.text(0.06, 0.054, "red = fails smooth candidate", fontsize=12.3, color=RED)
    info_ax.text(0.06, 0.016, "black/gold = smooth T(C)", fontsize=11.6, color=INK)

    fig.savefig(out, dpi=200)
    plt.close(fig)


def marker_panel(
    payload: dict,
    out: Path,
    *,
    title: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    kept_x, kept_y = arrays(payload, "retained")
    cut_x, cut_y = arrays(payload, "removed")
    kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
    cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])

    fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.scatter(
        kept_x[kept_mask],
        kept_y[kept_mask],
        s=8,
        c=BLUE,
        alpha=0.16,
        linewidths=0,
        label=f"retained event sample ({int(kept_mask.sum()):,})",
        rasterized=True,
    )
    ax.scatter(
        cut_x[cut_mask],
        cut_y[cut_mask],
        s=14,
        c=RED,
        alpha=0.76,
        linewidths=0,
        label=f"events removed by cut ({int(cut_mask.sum()):,})",
        rasterized=True,
    )

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile", fontsize=18, color=INK)
    ax.set_ylabel(CALO_LOG_LABEL, fontsize=18, color=INK)
    ax.tick_params(axis="both", labelsize=14, colors=INK, length=5, width=1.2)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color(INK)
    ax.grid(False)

    ax.legend(
        loc="lower left",
        frameon=False,
        fontsize=15,
        markerscale=1.8,
        handletextpad=0.5,
        borderaxespad=0.8,
    )
    fig.text(0.08, 0.94, title, ha="left", va="top", fontsize=25, weight="bold", color=INK)
    fig.subplots_adjust(left=0.095, right=0.975, top=0.875, bottom=0.12)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def main() -> None:
    payload = load_markers()
    full = OUTDIR / "the32_low_calo_marker_only_full_v1.png"
    central = OUTDIR / "the32_low_calo_marker_only_central_v1.png"
    clean = OUTDIR / "the32_low_calo_marker_audience_clean_v1.png"
    binned = OUTDIR / "the32_low_calo_marker_audience_5pct_bins_v1.png"
    matched = OUTDIR / "the32_low_calo_marker_matched_samples_v1.png"
    blair = OUTDIR / "the32_low_calo_blair_mattermost_v1.png"
    blair_canvas = OUTDIR / "the32_low_calo_blair_mattermost_canvas_v1.png"
    blair_canvas_full = OUTDIR / "the32_low_calo_blair_mattermost_canvas_0to80_v1.png"
    blair_canvas_full_ywide = OUTDIR / "the32_low_calo_blair_mattermost_canvas_0to80_ywide_v1.png"
    qualitative = OUTDIR / "the32_low_calo_centrality_regions_qualitative_v1.png"
    smooth_candidate = OUTDIR / "the32_low_calo_smooth_threshold_candidate_0to80_v1.png"
    marker_panel(
        payload,
        full,
        title="Low-calo tail in raw event markers",
        xlim=(0, 80),
        ylim=(2.70, 3.45),
    )
    marker_panel(
        payload,
        central,
        title="Central events: removed markers sit below the retained band",
        xlim=(0, 35),
        ylim=(2.70, 3.42),
    )
    audience_clean_panel(payload, clean)
    audience_binned_panel(payload, binned)
    audience_matched_panel(payload, matched)
    blair_mattermost_panel(payload, blair)
    blair_canvas_panel(payload, blair_canvas)
    blair_canvas_panel(payload, blair_canvas_full, x_hi=80.0)
    blair_canvas_panel(payload, blair_canvas_full_ywide, x_hi=80.0, ylim=(2.35, 3.36))
    qualitative_region_panel(payload, qualitative)
    smooth_candidate_panel(payload, smooth_candidate)
    manifest = {
        "marker_source": str(MARKER_JSON),
        "outputs": [
            str(full),
            str(central),
            str(clean),
            str(binned),
            str(matched),
            str(blair),
            str(blair_canvas),
            str(blair_canvas_full),
            str(blair_canvas_full_ywide),
            str(qualitative),
            str(smooth_candidate),
        ],
        "event_counts": payload["event_counts"],
        "plot_contract": "Markers only: no threshold line, no fit, no heatmap, no bar chart.",
        "audience_clean_note": "Audience-clean panel uses deterministic display-only centrality jitter and stratified marker samples to avoid overplotting.",
        "audience_binned_note": "5pct-bins panel displays integer centrality measurements within 5% bins to avoid barcode-like integer gaps.",
        "audience_matched_note": "Matched panel caps retained and removed markers to the same per-centrality count so the separation is visible without implying red is the majority.",
        "blair_mattermost_note": "Mattermost panel is a matched-marker event-level view for communication; text message should carry the quantitative construction details.",
        "blair_canvas_note": "Canvas panel places the construction and cut-effect statistics inside the plotting area while avoiding the main data bands.",
        "smooth_candidate_note": "Smooth candidate panel recolors the displayed marker sample by the smooth threshold approximation; it is not the exact validated training cut until separately rerun/revalidated.",
    }
    (OUTDIR / "the32_low_calo_marker_scatter_v1.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(full)
    print(central)


if __name__ == "__main__":
    main()
