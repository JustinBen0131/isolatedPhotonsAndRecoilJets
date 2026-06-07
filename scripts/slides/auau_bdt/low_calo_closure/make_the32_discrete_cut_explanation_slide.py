#!/usr/bin/env python3
"""Build a clear slide explaining the validated THE-32 discrete low-calo cut."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np


OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
CUT_JSON = OUTDIR / "the32_low_calo_cut_v1.json"
MARKER_JSON = OUTDIR / "the32_low_calo_marker_sample_v1.json"
SLIDE_DIR = OUTDIR / "slide_ready_discrete_cut_20260606"
PNG = SLIDE_DIR / "the32_discrete_low_calo_cut_explanation_slide_20260606.png"
SCRIPT = SLIDE_DIR / "the32_discrete_low_calo_cut_explanation_script_20260606.md"
MANIFEST = SLIDE_DIR / "the32_discrete_low_calo_cut_explanation_manifest_20260606.json"

BLUE = "#1f77b4"
RED = "#d62728"
INK = "#172033"
MUTED = "#566475"
GRID = "#D0D5DD"
CUT = "#111827"
GOLD = "#F59E0B"
SOFT_RED = "#FFF4ED"
SOFT_BLUE = "#EEF6FF"
SOFT_GREEN = "#F0FDF4"
SOFT_GRAY = "#F8FAFC"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def load_payloads() -> tuple[dict, dict]:
    return json.loads(CUT_JSON.read_text()), json.loads(MARKER_JSON.read_text())


def arrays(payload: dict, kind: str) -> tuple[np.ndarray, np.ndarray]:
    markers = payload["markers"][kind]
    return (
        np.asarray(markers["centrality"], dtype=np.float32),
        np.asarray(markers["log10_calo"], dtype=np.float32),
    )


def step_arrays(envelope: list[dict], xlim: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    for row in envelope:
        lo = max(float(row["cent_lo"]), xlim[0])
        hi = min(float(row["cent_hi"]), xlim[1])
        if hi <= xlim[0] or lo >= xlim[1]:
            continue
        xs.extend([lo, hi])
        ys.extend([float(row["threshold"]), float(row["threshold"])])
    return np.asarray(xs), np.asarray(ys)


def sample_by_bin(mask: np.ndarray, x: np.ndarray, *, per_bin: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for lo in np.arange(0, 80, 1):
        idx = np.flatnonzero(mask & (x >= lo) & (x < lo + 1))
        if idx.size > per_bin:
            idx = rng.choice(idx, size=per_bin, replace=False)
            idx.sort()
        if idx.size:
            selected.append(idx)
    return np.concatenate(selected) if selected else np.array([], dtype=np.int64)


def jitter(x: np.ndarray, *, seed: int, width: float = 0.40) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return x + rng.uniform(-width, width, size=x.size)


def add_round_box(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    facecolor: str,
    edgecolor: str = GRID,
    radius: float = 0.018,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        linewidth=1.2,
        edgecolor=edgecolor,
        facecolor=facecolor,
    )
    ax.add_patch(patch)
    return patch


def fmt_threshold(row: dict) -> str:
    return f"{int(row['cent_lo']):02d}-{int(row['cent_hi']):02d}: {float(row['threshold']):.3f}"


def render_slide() -> None:
    cut, markers = load_payloads()
    envelope = cut["envelope"]
    event_counts = markers["event_counts"]
    removed_events = int(event_counts["event_removed"])
    total_events = int(event_counts["event_total_after_per_cache_dedup"])
    removed_pct = 100.0 * removed_events / total_events
    kept_x, kept_y = arrays(markers, "retained")
    cut_x, cut_y = arrays(markers, "removed")
    xlim = (0.0, 35.0)
    ylim = (2.80, 3.36)
    visible_envelope = [row for row in envelope if float(row["cent_hi"]) <= xlim[1]]
    kept_mask = (kept_x >= xlim[0]) & (kept_x <= xlim[1]) & (kept_y >= ylim[0]) & (kept_y <= ylim[1])
    cut_mask = (cut_x >= xlim[0]) & (cut_x <= xlim[1]) & (cut_y >= ylim[0]) & (cut_y <= ylim[1])
    kept_idx = sample_by_bin(kept_mask, kept_x, per_bin=105, seed=6101)
    cut_idx = sample_by_bin(cut_mask, cut_x, per_bin=125, seed=6102)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.patch.set_facecolor("white")
    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.axis("off")

    canvas.text(
        0.055,
        0.942,
        "Per-event total-calo energy veto before BDT training",
        fontsize=27.0,
        fontweight="bold",
        color=INK,
        va="top",
    )
    canvas.text(
        0.055,
        0.878,
        "Definition: event-level veto on total CEMC+IHCal+OHCal energy, not a photon-cluster or BDT-score cut.",
        fontsize=14.2,
        color=MUTED,
        va="top",
    )

    ax = fig.add_axes([0.060, 0.250, 0.590, 0.585])
    ax.set_facecolor("white")
    ax.scatter(
        jitter(kept_x[kept_idx], seed=6103),
        kept_y[kept_idx],
        s=9,
        c=BLUE,
        alpha=0.20,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        jitter(cut_x[cut_idx], seed=6104),
        cut_y[cut_idx],
        s=14,
        c=RED,
        alpha=0.60,
        linewidths=0,
        rasterized=True,
    )
    sx, sy = step_arrays(visible_envelope, xlim)
    ax.plot(sx, sy, color="white", lw=7.0, solid_capstyle="butt", zorder=8)
    ax.plot(sx, sy, color=CUT, lw=3.4, solid_capstyle="butt", zorder=9)
    ax.plot(sx, sy, color=GOLD, lw=1.35, solid_capstyle="butt", zorder=10)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Centrality percentile", fontsize=15.3, color=INK)
    ax.set_ylabel(r"$\log_{10}(E_{\rm calo}^{\rm total}+1)$", fontsize=15.3, color=INK)
    ax.tick_params(axis="both", labelsize=12.2, colors=INK, length=4.8, width=1.05)
    for spine in ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(1.05)
    ax.grid(False)
    label_box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.4}
    ax.text(0.020, 0.967, "blue = kept", transform=ax.transAxes, color=BLUE, fontsize=14.4, va="top", bbox=label_box)
    ax.text(0.222, 0.967, "red = removed", transform=ax.transAxes, color=RED, fontsize=14.4, va="top", bbox=label_box)
    ax.text(0.472, 0.967, "black/gold = applied threshold", transform=ax.transAxes, color=INK, fontsize=14.4, va="top", bbox=label_box)

    add_round_box(canvas, 0.685, 0.655, 0.265, 0.180, facecolor=SOFT_BLUE)
    canvas.text(0.705, 0.805, "Data shown", fontsize=17.0, fontweight="bold", color=INK, va="top")
    canvas.text(0.705, 0.759, "Each marker is one merged-embedding event.", fontsize=12.8, color=INK, va="top")
    canvas.text(0.705, 0.717, "Samples: Photon12/20 and Jet12/20/30/40.", fontsize=12.0, color=INK, va="top")
    canvas.text(0.705, 0.676, "Unweighted diagnostic view, 0-35% centrality.", fontsize=11.9, color=MUTED, va="top")

    add_round_box(canvas, 0.685, 0.390, 0.265, 0.225, facecolor=SOFT_RED)
    canvas.text(0.705, 0.585, "How the thresholds were set", fontsize=15.8, fontweight="bold", color=INK, va="top")
    canvas.text(0.705, 0.545, r"Work separately in each 5% centrality bin.", fontsize=11.8, color=INK, va="top")
    canvas.text(0.705, 0.510, r"Center: median of  $y=\log_{10}(E_{\rm calo}^{\rm total}+1)$", fontsize=11.4, color=INK, va="top")
    canvas.text(0.705, 0.475, "Width: median event-to-event spread.", fontsize=11.3, color=INK, va="top")
    canvas.text(0.705, 0.440, r"$T_{\rm bin}=\max(\mathrm{lower\ 0.1\%\ floor},$", fontsize=11.5, color=RED, fontweight="bold", va="top")
    canvas.text(0.735, 0.412, r"$\mathrm{center}-5\times\mathrm{width})$", fontsize=11.5, color=RED, fontweight="bold", va="top")

    add_round_box(canvas, 0.685, 0.185, 0.265, 0.155, facecolor=SOFT_GREEN)
    canvas.text(0.705, 0.315, "Applied veto + check", fontsize=16.6, fontweight="bold", color=INK, va="top")
    canvas.text(0.705, 0.279, r"Reject whole event if  $y<T_{\rm bin}$", fontsize=12.8, color=RED, fontweight="bold", va="top")
    canvas.text(
        0.705,
        0.250,
        f"{removed_events:,} diagnostic events removed ({removed_pct:.1f}%)",
        fontsize=11.8,
        color=INK,
        fontweight="bold",
        va="top",
    )
    canvas.text(0.705, 0.221, "0 retained events/candidates below threshold", fontsize=11.2, color=INK, va="top")
    canvas.text(0.705, 0.196, "BDT validation AUC stable: 0.894693 -> 0.895321", fontsize=10.4, color=INK, va="top")

    table_rows = [fmt_threshold(row) for row in visible_envelope]
    add_round_box(canvas, 0.060, 0.040, 0.890, 0.112, facecolor=SOFT_GRAY)
    canvas.text(0.078, 0.128, "Frozen 5% centrality-bin thresholds shown for the plotted range", fontsize=14.0, fontweight="bold", color=INK, va="top")
    for i, text in enumerate(table_rows):
        canvas.text(0.078 + i * 0.121, 0.086, text, fontsize=11.1, color=INK, va="top")

    SLIDE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG, dpi=200)
    plt.close(fig)

    SCRIPT.write_text(
        "\n".join(
            [
                "# Discrete Total-Calo Event Veto Explanation",
                "",
                "This slide shows the event-level total calorimeter energy check that was applied before BDT training.",
                "Each point is one merged-embedding event from Photon12/20 and Jet12/20/30/40.",
                "The plotted view is restricted to 0-35% centrality because that is the cleanest visual range for the separated low-energy excess.",
                "The vertical axis is log10 of the total CEMC plus IHCal plus OHCal event energy plus one.",
                "The black and gold step line is the actual fixed threshold table used in the retraining.",
                "In each 5% centrality bin, the threshold is set from the event-energy distribution: median center, median event-to-event spread, and a threshold five widths below the center, bounded by the lower 0.1% tail.",
                "Events below that bin threshold are removed before the training matrix is built.",
                f"The event-level diagnostic sample removes {removed_events:,} of {total_events:,} events, or {removed_pct:.1f}%.",
                "This is not a cut on the photon candidate, the BDT score, the truth label, or the source sample.",
                "The closure check is that after upstream filtering and full-stat validation, zero retained events and zero retained candidates remain below the threshold, while the noIso AUC stays stable.",
                "",
            ]
        )
    )
    MANIFEST.write_text(
        json.dumps(
            {
                "slide_png": str(PNG),
                "speaker_script": str(SCRIPT),
                "cut_json": str(CUT_JSON),
                "marker_json": str(MARKER_JSON),
                "cut_variable": cut["cut_variable"],
                "centrality_bins": "5%-wide bins over 0-80%",
                "plotted_centrality_range": "0-35%",
                "plotted_threshold_bins": [fmt_threshold(row) for row in visible_envelope],
                "threshold_derivation": (
                    "For each 5% centrality bin, compute the median of y, "
                    "compute the median event-to-event spread around that center, "
                    "then set T_bin=max(lower 0.1% floor, center - 5*width)."
                ),
                "validated_training_cut": "event rejected if log10(CEMC+IHCal+OHCal+1) is below the fixed threshold for its 5% centrality bin",
                "diagnostic_events_removed": removed_events,
                "diagnostic_event_total": total_events,
                "diagnostic_events_removed_percent": removed_pct,
                "note": "This slide explains the discrete validated table, not the smooth approximation candidate.",
            },
            indent=2,
        )
        + "\n"
    )


def main() -> int:
    render_slide()
    print(PNG)
    print(SCRIPT)
    print(MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
