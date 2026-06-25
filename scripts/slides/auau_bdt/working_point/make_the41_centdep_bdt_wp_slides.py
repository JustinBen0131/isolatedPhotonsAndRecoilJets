#!/usr/bin/env python3
"""Extract and render THE-41 centrality-dependent BDT WP slide candidates."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


CENT_EDGES = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0]
ET_RANGE = (15.0, 35.0)
ET_EDGES = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0]
TARGETS = [0.90, 0.80, 0.70]
TARGET_COLORS = {0.90: "#D65F9E", 0.80: "#1B9E77", 0.70: "#2F6FB3"}
TARGET_NAMES = {0.90: "WP90", 0.80: "WP80", 0.70: "WP70"}
INK = "#111827"
MUTED = "#4B5563"
GRID = "#E5E7EB"
PANEL = "#F8FAFC"
PANEL_EDGE = "#CBD5E1"
YELLOW_TINT = "#FFF4C7"
YELLOW_EDGE = "#F4D35E"
GREEN_TINT = "#E8F6EF"
GREEN_EDGE = "#A7D8BF"
BLUE_TINT = "#EAF2FB"
BLUE_EDGE = "#9EC5E8"
PINK_TINT = "#FCE7F3"
PINK_EDGE = "#F3A6CE"


def score_cache_paths(path: Path) -> list[Path]:
    if path.suffix == ".npz":
        return [path]
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise FileNotFoundError(f"No score caches listed in {path}")
    return [Path(line) for line in lines]


def parse_edges(spec: str) -> list[float]:
    values = [float(tok) for tok in spec.split(",") if tok.strip()]
    if len(values) < 2 or any(values[i + 1] <= values[i] for i in range(len(values) - 1)):
        raise ValueError(f"Bad edge list: {spec}")
    return values


def threshold_for_target(sig, bkg, target: float) -> dict:
    import numpy as np

    if sig.size == 0 or bkg.size == 0:
        return {
            "threshold": math.nan,
            "threshold_stat_err": math.nan,
            "threshold_stat_err_low": math.nan,
            "threshold_stat_err_high": math.nan,
            "signal_efficiency": math.nan,
            "background_fake_rate": math.nan,
            "signal_entries": int(sig.size),
            "background_entries": int(bkg.size),
        }
    q = max(0.0, min(1.0, 1.0 - target))
    thr = float(np.quantile(sig, q))
    # One-sigma quantile uncertainty from binomial order-statistics:
    # sigma_q = sqrt(q(1-q)/N), mapped back to score through local quantiles.
    if sig.size > 1:
        sigma_q = math.sqrt(max(0.0, q * (1.0 - q)) / float(sig.size))
        q_low = max(0.0, q - sigma_q)
        q_high = min(1.0, q + sigma_q)
        thr_low = float(np.quantile(sig, q_low))
        thr_high = float(np.quantile(sig, q_high))
        err_low = max(0.0, thr - thr_low)
        err_high = max(0.0, thr_high - thr)
        err = 0.5 * (err_low + err_high)
    else:
        err_low = math.nan
        err_high = math.nan
        err = math.nan
    return {
        "threshold": thr,
        "threshold_stat_err": err,
        "threshold_stat_err_low": err_low,
        "threshold_stat_err_high": err_high,
        "signal_efficiency": float(np.mean(sig > thr)),
        "background_fake_rate": float(np.mean(bkg > thr)),
        "signal_entries": int(sig.size),
        "background_entries": int(bkg.size),
    }


def extract(args: argparse.Namespace) -> dict:
    import numpy as np

    cent_edges = parse_edges(args.cent_edges)
    et_edges = parse_edges(args.et_edges)
    et_min, et_max = [float(x) for x in args.et_range.split(",")]
    scan_thresholds = np.linspace(0.0, 0.95, 96)
    by_cent = {(cent_edges[i], cent_edges[i + 1], cls): [] for i in range(len(cent_edges) - 1) for cls in (0, 1)}
    by_cell = {
        (clo, chi, elo, ehi, cls): []
        for clo, chi in zip(cent_edges[:-1], cent_edges[1:])
        for elo, ehi in zip(et_edges[:-1], et_edges[1:])
        for cls in (0, 1)
    }
    files_loaded = 0
    rows_loaded = 0
    missing_files: list[str] = []

    for cache in score_cache_paths(args.manifest):
        if not cache.exists():
            missing_files.append(str(cache))
            continue
        with np.load(cache, allow_pickle=True) as z:
            missing = [key for key in ("is_signal", "cluster_Et", "centrality", args.score_key) if key not in z.files]
            if missing:
                raise KeyError(f"{cache} missing keys {missing}")
            y = np.asarray(z["is_signal"], dtype=np.int8)
            et = np.asarray(z["cluster_Et"], dtype=np.float32)
            cent = np.asarray(z["centrality"], dtype=np.float32)
            score = np.asarray(z[args.score_key], dtype=np.float32)
            base = (
                np.isfinite(et)
                & np.isfinite(cent)
                & np.isfinite(score)
                & np.isin(y, [0, 1])
                & (et >= et_min)
                & (et < et_max)
                & (cent >= cent_edges[0])
                & (cent < cent_edges[-1])
            )
            files_loaded += 1
            rows_loaded += int(np.sum(base))
            if not np.any(base):
                continue
            for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
                cmask = base & (cent >= clo) & (cent < chi)
                if not np.any(cmask):
                    continue
                for cls in (0, 1):
                    arr = score[cmask & (y == cls)]
                    if arr.size:
                        by_cent[(clo, chi, cls)].append(arr.astype(np.float32, copy=False))
                for elo, ehi in zip(et_edges[:-1], et_edges[1:]):
                    emask = cmask & (et >= elo) & (et < ehi)
                    if not np.any(emask):
                        continue
                    for cls in (0, 1):
                        arr = score[emask & (y == cls)]
                        if arr.size:
                            by_cell[(clo, chi, elo, ehi, cls)].append(arr.astype(np.float32, copy=False))

    if missing_files:
        raise FileNotFoundError(f"Missing score-cache files: {missing_files[:3]} (n={len(missing_files)})")

    rows = []
    scans = []
    et_cells = []
    for target in TARGETS:
        for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
            sig = np.concatenate(by_cent[(clo, chi, 1)]) if by_cent[(clo, chi, 1)] else np.array([], dtype=np.float32)
            bkg = np.concatenate(by_cent[(clo, chi, 0)]) if by_cent[(clo, chi, 0)] else np.array([], dtype=np.float32)
            item = threshold_for_target(sig, bkg, target)
            rows.append(
                {
                    "target_signal_efficiency": target,
                    "wp_label": TARGET_NAMES[target],
                    "centrality_min": clo,
                    "centrality_max": chi,
                    "centrality_center": 0.5 * (clo + chi),
                    "centrality_label": f"{int(clo)}-{int(chi)}%",
                    **item,
                }
            )
            for elo, ehi in zip(et_edges[:-1], et_edges[1:]):
                sig_cell = np.concatenate(by_cell[(clo, chi, elo, ehi, 1)]) if by_cell[(clo, chi, elo, ehi, 1)] else np.array([], dtype=np.float32)
                bkg_cell = np.concatenate(by_cell[(clo, chi, elo, ehi, 0)]) if by_cell[(clo, chi, elo, ehi, 0)] else np.array([], dtype=np.float32)
                cell_item = threshold_for_target(sig_cell, bkg_cell, target)
                et_cells.append(
                    {
                        "target_signal_efficiency": target,
                        "wp_label": TARGET_NAMES[target],
                        "centrality_min": clo,
                        "centrality_max": chi,
                        "centrality_center": 0.5 * (clo + chi),
                        "centrality_label": f"{int(clo)}-{int(chi)}%",
                        "et_min": elo,
                        "et_max": ehi,
                        "et_center": 0.5 * (elo + ehi),
                        **cell_item,
                    }
                )
    for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
        sig = np.concatenate(by_cent[(clo, chi, 1)]) if by_cent[(clo, chi, 1)] else np.array([], dtype=np.float32)
        bkg = np.concatenate(by_cent[(clo, chi, 0)]) if by_cent[(clo, chi, 0)] else np.array([], dtype=np.float32)
        for thr in scan_thresholds:
            scans.append(
                {
                    "centrality_min": clo,
                    "centrality_max": chi,
                    "centrality_center": 0.5 * (clo + chi),
                    "centrality_label": f"{int(clo)}-{int(chi)}%",
                    "threshold": float(thr),
                    "signal_efficiency": float(np.mean(sig > thr)) if sig.size else math.nan,
                    "background_fake_rate": float(np.mean(bkg > thr)) if bkg.size else math.nan,
                }
            )

    return {
        "schema": "THE41_CENTRALITY_DEPENDENT_BDT_WP_V1",
        "manifest": str(args.manifest),
        "score_key": args.score_key,
        "et_range": [et_min, et_max],
        "et_edges": et_edges,
        "cent_edges": cent_edges,
        "targets": TARGETS,
        "files_loaded": files_loaded,
        "rows_loaded": rows_loaded,
        "source_label": args.source_label,
        "model_label": args.model_label,
        "training_sample": args.training_sample,
        "training_inputs": args.training_inputs,
        "rows": rows,
        "et_cells": et_cells,
        "efficiency_scan": scans,
    }


def flat_fit_rows(payload: dict) -> list[dict]:
    import numpy as np

    if payload.get("flat_rows"):
        return list(payload["flat_rows"])

    cells = payload.get("et_cells") or []
    if not cells:
        return payload["rows"]
    cent_labels = [f"{int(lo)}-{int(hi)}%" for lo, hi in zip(payload["cent_edges"][:-1], payload["cent_edges"][1:])]
    rows = []
    for target in TARGETS:
        for label in cent_labels:
            sub = [
                c
                for c in cells
                if c["centrality_label"] == label
                and abs(float(c["target_signal_efficiency"]) - target) < 1e-9
                and math.isfinite(float(c["threshold"]))
            ]
            if not sub:
                continue
            y = np.array([float(c["threshold"]) for c in sub], dtype=float)
            const = float(np.mean(y))
            residuals = y - const
            rows.append(
                {
                    "target_signal_efficiency": target,
                    "wp_label": TARGET_NAMES[target],
                    "centrality_min": float(sub[0]["centrality_min"]),
                    "centrality_max": float(sub[0]["centrality_max"]),
                    "centrality_center": float(sub[0]["centrality_center"]),
                    "centrality_label": label,
                    "threshold": const,
                    "signal_efficiency": target,
                    "background_fake_rate": math.nan,
                    "signal_entries": int(sum(int(c["signal_entries"]) for c in sub)),
                    "background_entries": int(sum(int(c["background_entries"]) for c in sub)),
                    "n_et_points": int(len(sub)),
                    "threshold_stat_err": math.nan,
                    "threshold_stat_err_low": math.nan,
                    "threshold_stat_err_high": math.nan,
                    "flat_max_abs_residual": float(np.max(np.abs(residuals))),
                    "flat_rms_residual": float(np.sqrt(np.mean(residuals**2))),
                }
            )
    return rows


def fit_lines(rows: list[dict]) -> dict[float, dict]:
    import numpy as np

    fits = {}
    for target in TARGETS:
        sub = [r for r in rows if abs(float(r["target_signal_efficiency"]) - target) < 1e-9]
        x = np.array([float(r["centrality_center"]) for r in sub], dtype=float)
        y = np.array([float(r["threshold"]) for r in sub], dtype=float)
        slope, intercept = np.polyfit(x, y, 1)
        pred = slope * x + intercept
        fits[target] = {
            "slope": float(slope),
            "intercept": float(intercept),
            "max_abs_residual": float(np.max(np.abs(y - pred))),
            "rms_residual": float(np.sqrt(np.mean((y - pred) ** 2))),
        }
    return fits


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.0,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )
    return plt


def add_box(fig, rect, facecolor, edgecolor, title, body, title_color=INK, body_size=11.0):
    import matplotlib.pyplot as plt

    ax = fig.add_axes(rect)
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=facecolor, edgecolor=edgecolor, linewidth=1.0))
    ax.text(0.045, 0.78, title, fontsize=14.5, fontweight="bold", color=title_color, ha="left", va="center")
    ax.text(0.055, 0.39, body, fontsize=body_size, color=INK, ha="left", va="center", linespacing=1.18)
    return ax


def rows_by_target(rows: list[dict], target: float) -> list[dict]:
    return [r for r in rows if abs(float(r["target_signal_efficiency"]) - target) < 1e-9]


def audience_sample_label(payload: dict) -> str:
    sample = str(payload.get("training_sample", "")).strip()
    if "Branch A" in sample or not sample:
        return "Photon+Jet Embedded 12+20, Inclusive Jet Embedded 12+20+30+40"
    return sample


def draw_slide_points(payload: dict, flat_rows: list[dict], out: Path) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    setup_matplotlib()
    cells = payload.get("et_cells", [])
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(0.045, 0.965, "Flat fits define one BDT cut per centrality bin", fontsize=27.0, fontweight="bold", color=INK, ha="left", va="top")

    subtitle = fig.add_axes([0.045, 0.825, 0.560, 0.085])
    subtitle.axis("off")
    subtitle.add_patch(
        plt.matplotlib.patches.FancyBboxPatch(
            (0.008, 0.060),
            0.984,
            0.880,
            transform=subtitle.transAxes,
            boxstyle="round,pad=0.018,rounding_size=0.045",
            facecolor=BLUE_TINT,
            edgecolor=BLUE_EDGE,
            linewidth=1.4,
        )
    )
    subtitle.text(
        0.030,
        0.54,
        "Flat fits versus $E_T$ give one threshold per centrality bin for the centrality-dependent working point.",
        transform=subtitle.transAxes,
        fontsize=15.4,
        color=INK,
        ha="left",
        va="center",
    )

    legend = fig.add_axes([0.635, 0.825, 0.325, 0.085])
    legend.axis("off")
    legend.add_patch(plt.Rectangle((0, 0), 1, 1, transform=legend.transAxes, facecolor=PANEL, edgecolor=PANEL_EDGE, linewidth=1.0))
    legend.text(0.040, 0.74, "Legend", transform=legend.transAxes, fontsize=14.0, fontweight="bold", color=INK, ha="left", va="center", zorder=3)
    legend.plot([0.045, 0.145], [0.34, 0.34], color=MUTED, lw=1.5, transform=legend.transAxes, clip_on=False, zorder=3)
    legend.plot([0.095], [0.34], "o", ms=5.3, color=MUTED, mec="white", mew=0.6, transform=legend.transAxes, clip_on=False, zorder=4)
    legend.text(0.175, 0.34, "point + error + fit", transform=legend.transAxes, fontsize=12.2, color=MUTED, ha="left", va="center", zorder=3)
    legend.plot([0.500, 0.500], [0.18, 0.86], color=PANEL_EDGE, lw=1.0, transform=legend.transAxes, clip_on=False, zorder=3)
    for y, target in zip([0.72, 0.48, 0.24], TARGETS):
        legend.plot([0.555], [y], "o", ms=5.4, color=TARGET_COLORS[target], mec="white", mew=0.6, transform=legend.transAxes, clip_on=False, zorder=4)
        legend.text(0.590, y, f"{TARGET_NAMES[target]} = {int(target*100)}%", transform=legend.transAxes, fontsize=13.0, fontweight="bold", color=TARGET_COLORS[target], ha="left", va="center", zorder=3)

    cent_labels = [r["centrality_label"] for r in rows_by_target(flat_rows, 0.80)]
    panel_rects = [
        [0.065, 0.535, 0.205, 0.235],
        [0.295, 0.535, 0.205, 0.235],
        [0.525, 0.535, 0.205, 0.235],
        [0.755, 0.535, 0.205, 0.235],
        [0.105, 0.265, 0.225, 0.235],
        [0.390, 0.265, 0.225, 0.235],
        [0.675, 0.265, 0.225, 0.235],
    ]
    finite_thresholds = [float(c["threshold"]) for c in cells if math.isfinite(float(c["threshold"]))]
    finite_thresholds += [float(r["threshold"]) for r in flat_rows if math.isfinite(float(r["threshold"]))]
    ymin = max(0.20, min(finite_thresholds) - 0.055)
    ymax = min(0.90, max(finite_thresholds) + 0.055)
    for idx, label in enumerate(cent_labels):
        ax = fig.add_axes(panel_rects[idx])
        ax.set_xlim(payload["et_edges"][0] - 0.7, payload["et_edges"][-1] + 0.7)
        ax.set_ylim(ymin, ymax)
        ax.grid(True, color=GRID, lw=0.7, alpha=0.85)
        ax.set_title(label, fontsize=15.0, fontweight="bold", pad=7)
        ax.tick_params(labelsize=13.0, length=3)
        if idx in (0, 4):
            ax.set_ylabel("BDT score cut", fontsize=13.4)
        else:
            ax.set_yticklabels([])
        if idx >= 4:
            ax.set_xlabel("cluster $E_T$ bin center [GeV]", fontsize=13.4)
        else:
            ax.set_xticklabels([])
        for target in TARGETS:
            sub = [
                c
                for c in cells
                if c["centrality_label"] == label and abs(float(c["target_signal_efficiency"]) - target) < 1e-9
            ]
            sub = sorted(sub, key=lambda c: float(c["et_center"]))
            x = np.array([float(c["et_center"]) for c in sub], dtype=float)
            y = np.array([float(c["threshold"]) for c in sub], dtype=float)
            err_low = np.array([float(c.get("threshold_stat_err_low", math.nan)) for c in sub], dtype=float)
            err_high = np.array([float(c.get("threshold_stat_err_high", math.nan)) for c in sub], dtype=float)
            if np.all(np.isfinite(err_low)) and np.all(np.isfinite(err_high)):
                ax.errorbar(
                    x,
                    y,
                    yerr=np.vstack([err_low, err_high]),
                    fmt="o",
                    ms=4.8,
                    color=TARGET_COLORS[target],
                    mec="white",
                    mew=0.6,
                    elinewidth=1.0,
                    capsize=2.0,
                    capthick=0.9,
                    alpha=0.96,
                    zorder=3,
                )
            else:
                ax.plot(x, y, "o", ms=4.8, color=TARGET_COLORS[target], mec="white", mew=0.6, alpha=0.96, zorder=3)
            row = next(r for r in rows_by_target(flat_rows, target) if r["centrality_label"] == label)
            const = float(row["threshold"])
            ax.axhline(const, color=TARGET_COLORS[target], lw=1.45, alpha=0.95)
        if idx == 6:
            ax.text(0.96, 0.07, f"{int(payload['et_range'][0])} <= $E_T$ < {int(payload['et_range'][1])} GeV", transform=ax.transAxes, fontsize=13.0, color=MUTED, ha="right", va="bottom")

    tab = fig.add_axes([0.080, 0.025, 0.840, 0.185])
    tab.axis("off")
    tab.add_patch(plt.Rectangle((0, 0), 1, 1, transform=tab.transAxes, facecolor=PANEL, edgecolor=PANEL_EDGE, linewidth=1.0))
    tab.text(0.030, 0.78, "Flat-fit constants used as centrality-fit points", fontsize=15.0, fontweight="bold", color=INK, ha="left", va="center")
    labels = [r["centrality_label"] for r in rows_by_target(flat_rows, 0.80)]
    xs = np.linspace(0.385, 0.950, len(labels))
    for x, label in zip(xs, labels):
        tab.text(x, 0.78, label, fontsize=13.0, color=MUTED, ha="center", va="center")
    for j, target in enumerate(TARGETS):
        y = 0.55 - j * 0.19
        tab.text(0.055, y, f"{TARGET_NAMES[target]} ({int(target*100)}%)", fontsize=13.8, fontweight="bold", color=TARGET_COLORS[target], ha="left", va="center")
        for x, r in zip(xs, rows_by_target(flat_rows, target)):
            tab.text(x, y, f"{float(r['threshold']):.3f}", fontsize=13.6, color=INK, ha="center", va="center")
    tab.text(
        0.500,
        0.065,
        "Statistical note: error bars are included on the points; they are smaller than the markers on this scale.",
        fontsize=12.4,
        color=MUTED,
        ha="center",
        va="center",
    )

    fig.savefig(out, dpi=160)
    plt.close(fig)


def draw_slide_fits(payload: dict, fits: dict[float, dict], out: Path) -> None:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import FancyBboxPatch

    setup_matplotlib()
    rows = flat_fit_rows(payload)
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    nominal = 0.80
    nominal_fit = fits[nominal]
    fig.text(0.045, 0.965, "The 80% efficiency fit defines the tight BDT selection", fontsize=29.0, fontweight="bold", color=INK, ha="left", va="top")

    subtitle = fig.add_axes([0.045, 0.835, 0.690, 0.055])
    subtitle.axis("off")
    subtitle.add_patch(
        FancyBboxPatch(
            (0, 0),
            1,
            1,
            transform=subtitle.transAxes,
            boxstyle="round,pad=0.010,rounding_size=0.025",
            facecolor=BLUE_TINT,
            edgecolor=BLUE_EDGE,
            linewidth=1.3,
        )
    )
    subtitle.text(
        0.030,
        0.52,
        "Seven centrality-bin thresholds are fit versus centrality; WP80 sets the nominal split.",
        transform=subtitle.transAxes,
        fontsize=15.8,
        color=INK,
        ha="left",
        va="center",
    )

    ax = fig.add_axes([0.070, 0.245, 0.610, 0.545])
    grid = np.linspace(0, 80, 300)
    nominal_line = nominal_fit["slope"] * grid + nominal_fit["intercept"]
    nominal_band = nominal_fit["rms_residual"]
    ax.fill_between(
        grid,
        nominal_line - nominal_band,
        nominal_line + nominal_band,
        facecolor="#D6F3E5",
        edgecolor="none",
        alpha=0.90,
        label=f"WP80 RMS residual band ($\\pm${nominal_band:.4f})",
        zorder=1,
    )
    for target in TARGETS:
        sub = rows_by_target(rows, target)
        x = np.array([float(r["centrality_center"]) for r in sub])
        y = np.array([float(r["threshold"]) for r in sub])
        fit = fits[target]
        is_nominal = abs(target - nominal) < 1e-9
        alpha = 1.0 if is_nominal else 0.25
        lw = 3.2 if is_nominal else 1.55
        ms = 10.2 if is_nominal else 6.8
        zorder = 5 if is_nominal else 3
        label = "WP80 nominal fit" if is_nominal else f"{TARGET_NAMES[target]} context"
        ax.plot(
            grid,
            fit["slope"] * grid + fit["intercept"],
            lw=lw,
            color=TARGET_COLORS[target],
            alpha=alpha,
            label=label,
            zorder=zorder,
            solid_capstyle="round",
        )
        ax.plot(x, y, "o", ms=ms, color=TARGET_COLORS[target], alpha=alpha, mec="white", mew=1.1, zorder=zorder + 1)
    ax.set_xlim(0, 80)
    ymin = min(float(r["threshold"]) for r in rows) - 0.050
    ymax = max(float(r["threshold"]) for r in rows) + 0.125
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("centrality percentile c", fontsize=14.8)
    ax.set_ylabel("BDT score cut", fontsize=14.8)
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80])
    ax.tick_params(labelsize=13.6)
    ax.grid(True, color=GRID, lw=0.8, alpha=0.90)
    legend = ax.legend(loc="upper left", frameon=True, fontsize=13.4, facecolor="white", edgecolor=PANEL_EDGE, framealpha=0.92)
    legend.get_frame().set_linewidth(0.8)
    ax.set_title("Centrality-dependent BDT threshold from flat-fit constants", fontsize=17.0, fontweight="bold", pad=10)
    ax.text(
        0.965,
        0.115,
        f"Band: WP80 fit RMS residual = {nominal_band:.4f}",
        transform=ax.transAxes,
        fontsize=12.8,
        fontweight="bold",
        color=TARGET_COLORS[nominal],
        ha="right",
        va="bottom",
    )
    ax.text(
        0.965,
        0.060,
        f"Nominal: T$_{{80}}$(c) = {nominal_fit['intercept']:.4f} {nominal_fit['slope']:+.5f} c",
        transform=ax.transAxes,
        fontsize=14.2,
        fontweight="bold",
        color=TARGET_COLORS[nominal],
        ha="right",
        va="bottom",
    )

    decision = fig.add_axes([0.715, 0.485, 0.265, 0.330])
    decision.axis("off")
    decision.add_patch(
        FancyBboxPatch(
            (0, 0),
            1,
            1,
            transform=decision.transAxes,
            boxstyle="round,pad=0.010,rounding_size=0.025",
            facecolor=GREEN_TINT,
            edgecolor=GREEN_EDGE,
            linewidth=1.4,
        )
    )
    decision.add_patch(plt.Rectangle((0.000, 0.000), 0.018, 1.000, transform=decision.transAxes, facecolor=TARGET_COLORS[nominal], edgecolor="none"))
    decision.text(0.070, 0.865, "Final working definition", transform=decision.transAxes, fontsize=17.0, fontweight="bold", color=INK, ha="left", va="center")
    decision.text(0.070, 0.700, "Nominal tight-BDT cut:", transform=decision.transAxes, fontsize=14.4, color=INK, ha="left", va="center")
    decision.text(0.070, 0.570, f"T$_{{80}}$(c) = {nominal_fit['intercept']:.4f} {nominal_fit['slope']:+.5f} c", transform=decision.transAxes, fontsize=17.4, fontweight="bold", color=TARGET_COLORS[nominal], ha="left", va="center")
    decision.text(0.070, 0.455, f"Fit residuals: RMS={nominal_fit['rms_residual']:.4f}, max={nominal_fit['max_abs_residual']:.4f}", transform=decision.transAxes, fontsize=13.4, color=MUTED, ha="left", va="center")
    decision.text(0.070, 0.305, "Tight BDT", transform=decision.transAxes, fontsize=14.8, fontweight="bold", color=INK, ha="left", va="center")
    decision.text(0.395, 0.305, "score > T$_{80}$(c)", transform=decision.transAxes, fontsize=14.8, color=INK, ha="left", va="center")
    decision.text(0.070, 0.160, "Non-tight", transform=decision.transAxes, fontsize=14.8, fontweight="bold", color=INK, ha="left", va="center")
    decision.text(0.395, 0.160, "score <= T$_{80}$(c)", transform=decision.transAxes, fontsize=14.8, color=INK, ha="left", va="center")

    context = fig.add_axes([0.715, 0.255, 0.265, 0.190])
    context.axis("off")
    context.add_patch(
        FancyBboxPatch(
            (0, 0),
            1,
            1,
            transform=context.transAxes,
            boxstyle="round,pad=0.010,rounding_size=0.025",
            facecolor=PANEL,
            edgecolor=PANEL_EDGE,
            linewidth=1.0,
        )
    )
    context.text(0.060, 0.765, "Context curves", transform=context.transAxes, fontsize=14.8, fontweight="bold", color=INK, ha="left", va="center")
    context.text(0.060, 0.505, f"WP90: T(c) = {fits[0.90]['intercept']:.4f} {fits[0.90]['slope']:+.5f} c", transform=context.transAxes, fontsize=13.4, color=TARGET_COLORS[0.90], fontweight="bold", ha="left", va="center")
    context.text(0.060, 0.275, f"WP70: T(c) = {fits[0.70]['intercept']:.4f} {fits[0.70]['slope']:+.5f} c", transform=context.transAxes, fontsize=13.4, color=TARGET_COLORS[0.70], fontweight="bold", ha="left", va="center")

    foot = fig.add_axes([0.070, 0.030, 0.910, 0.150])
    foot.axis("off")
    foot.add_patch(
        FancyBboxPatch(
            (0, 0),
            1,
            1,
            transform=foot.transAxes,
            boxstyle="round,pad=0.008,rounding_size=0.018",
            facecolor=PANEL,
            edgecolor=PANEL_EDGE,
            linewidth=1.0,
        )
    )
    total_sig = sum(int(r["signal_entries"]) for r in rows_by_target(rows, 0.80))
    total_bkg = sum(int(r["background_entries"]) for r in rows_by_target(rows, 0.80))
    foot.text(0.030, 0.670, "Calibration sample", fontsize=15.6, fontweight="bold", color=INK, ha="left", va="center")
    foot.text(0.220, 0.670, audience_sample_label(payload), fontsize=14.8, color=INK, ha="left", va="center")
    foot.text(0.030, 0.315, "Inputs + rows", fontsize=15.6, fontweight="bold", color=INK, ha="left", va="center")
    input_line = f"{payload['training_inputs']}"
    count_line = (
        f"S={total_sig:,}, Incl.={total_bkg:,}; "
        f"{int(payload['et_range'][0])} <= cluster $E_T$ < {int(payload['et_range'][1])} GeV"
    )
    if payload.get("full_matrix_rows") is not None:
        count_line += f"; full matrix rows={int(payload['full_matrix_rows']):,}"
    foot.text(
        0.220,
        0.405,
        input_line,
        fontsize=13.9,
        color=MUTED,
        ha="left",
        va="center",
    )
    foot.text(0.220, 0.205, count_line, fontsize=13.9, color=MUTED, ha="left", va="center")

    fig.savefig(out, dpi=160)
    plt.close(fig)


def write_csv(payload: dict, path: Path, fits: dict[float, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "row_type",
        "target_signal_efficiency",
        "wp_label",
        "centrality_min",
        "centrality_max",
        "centrality_center",
        "centrality_label",
        "et_min",
        "et_max",
        "et_center",
        "threshold",
        "threshold_stat_err",
        "threshold_stat_err_low",
        "threshold_stat_err_high",
        "signal_efficiency",
        "background_fake_rate",
        "signal_entries",
        "background_entries",
        "signal_weight_sum",
        "background_weight_sum",
        "signal_effective_entries",
        "background_effective_entries",
        "n_et_points",
        "flat_max_abs_residual",
        "flat_rms_residual",
        "fit_intercept",
        "fit_slope",
        "fit_max_abs_residual",
        "fit_rms_residual",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in flat_fit_rows(payload):
            fit = fits[float(row["target_signal_efficiency"])]
            out = dict(row)
            out.update(
                {
                    "row_type": "flat_fit_constant",
                    "fit_intercept": fit["intercept"],
                    "fit_slope": fit["slope"],
                    "fit_max_abs_residual": fit["max_abs_residual"],
                    "fit_rms_residual": fit["rms_residual"],
                }
            )
            writer.writerow(out)
        for row in payload.get("et_cells", []):
            out = dict(row)
            out.update(
                {
                    "row_type": "et_bin_point",
                    "n_et_points": "",
                    "flat_max_abs_residual": "",
                    "flat_rms_residual": "",
                    "fit_intercept": "",
                    "fit_slope": "",
                    "fit_max_abs_residual": "",
                    "fit_rms_residual": "",
                }
            )
            writer.writerow(out)


def write_scripts(outdir: Path, payload: dict, fits: dict[float, dict]) -> dict[str, str]:
    script1 = outdir / "the41_centdep_bdt_wp_flat_points_script.md"
    script2 = outdir / "the41_centdep_bdt_wp_linear_fits_script.md"
    script1.write_text(
        "# THE-41 Slide Script - Centrality-Panel Working-Point Derivation\n\n"
        "First, I want to show the intermediate step before the centrality-dependent BDT cut is defined. "
        "Each panel is one fine centrality bin. Within that panel, the points are the BDT score cuts extracted independently in cluster E_T bins for the same target signal efficiencies.\n\n"
        "The color convention follows the sliding-isolation slides: pink is 90 percent, green is 80 percent, and blue is 70 percent signal efficiency. "
        "For each color, the horizontal line is the flat fit across E_T inside that centrality bin. "
        "Those flat-fit constants are the numbers carried forward to the centrality-fit slide.\n\n"
        "So the logic is deliberately two-step. First we avoid choosing an E_T-dependent shape by fitting a flat score cut in each centrality environment. "
        "Then, on the next slide, we fit those seven flat constants as a function of centrality to define the proposed runtime BDT threshold.\n"
    )
    wp80_fit = fits[0.80]
    script2.write_text(
        "# THE-41 Slide Script - Linear Centrality Fits\n\n"
        "Now I take the seven flat-fit constants from the previous slide and fit each efficiency target with a linear function of centrality. "
        "This is the direct analogue of the sliding-isolation cutoff procedure, but applied to the BDT score threshold.\n\n"
        "The main decision on this slide is that I am taking the 80 percent signal-efficiency line as the nominal working point. "
        f"That gives T80 of centrality equals {wp80_fit['intercept']:.4f} plus {wp80_fit['slope']:.5f} times the centrality percentile.\n\n"
        f"The shaded green band around that line is the RMS residual of the WP80 centrality fit, which is {wp80_fit['rms_residual']:.4f} in BDT score. "
        f"The largest WP80 point-to-line residual is {wp80_fit['max_abs_residual']:.4f}, so this is a compact visual check of how well the linear centrality model describes the seven flat-fit constants.\n\n"
        "The selection definition is then explicit. "
        "A candidate is tight BDT if its BDT score is above T80 of centrality. "
        "A candidate is non-tight if its BDT score is less than or equal to that threshold. "
        "The 90 and 70 percent curves are kept on the plot as context, but the green WP80 line is the cut definition I would carry into the shower-shape overlay step.\n"
    )
    return {"points_script": str(script1), "fits_script": str(script2)}


def render(args: argparse.Namespace) -> dict:
    payload = json.loads(args.input.read_text())
    flats = flat_fit_rows(payload)
    fits = fit_lines(flats)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    points_png = outdir / "the41_centdep_bdt_wp_flat_points_slide.png"
    fits_png = outdir / "the41_centdep_bdt_wp_linear_fits_slide.png"
    draw_slide_points(payload, flats, points_png)
    draw_slide_fits(payload, fits, fits_png)
    summary_csv = outdir / "the41_centdep_bdt_wp_summary.csv"
    write_csv(payload, summary_csv, fits)
    scripts = write_scripts(outdir, payload, fits)
    manifest = {
        "schema": "THE41_CENTDEP_BDT_WP_SLIDE_MANIFEST_V1",
        "input": str(args.input),
        "points_png": str(points_png),
        "fits_png": str(fits_png),
        "summary_csv": str(summary_csv),
        "speaker_scripts": scripts,
        "fits": {TARGET_NAMES[k]: v for k, v in fits.items()},
        "n_et_cells": len(payload.get("et_cells", [])),
        "n_flat_fit_constants": len(flats),
        "stat_uncertainty_method": payload.get(
            "stat_uncertainty_method",
            "One-sigma binomial order-statistic quantile uncertainty is drawn for each E_T-bin threshold point. Flat constants and centrality lines are simple unweighted fits.",
        ),
        "source_label": payload.get("source_label"),
        "model_label": payload.get("model_label"),
        "files_loaded": payload.get("files_loaded"),
        "rows_loaded": payload.get("rows_loaded"),
        "caveat": "Local full-slide candidates; Google Slides was not mutated.",
    }
    manifest_path = outdir / "the41_centdep_bdt_wp_slide_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["extract", "render"], required=True)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--input", type=Path)
    ap.add_argument("--json-out", type=Path)
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--score-key", default="score_globalEtCent1535_bdt_noIso")
    ap.add_argument("--cent-edges", default=",".join(str(x) for x in CENT_EDGES))
    ap.add_argument("--et-range", default="15,35")
    ap.add_argument("--et-edges", default=",".join(str(x) for x in ET_EDGES))
    ap.add_argument("--source-label", default="THE8_branchA_jet12_20_30_40_scorecache_fullstat_20260527")
    ap.add_argument("--model-label", default="global EtCent1535 no-iso BDT")
    ap.add_argument("--training-sample", default="Photon+Jet Embedded 12+20, Inclusive Jet Embedded 12+20+30+40")
    ap.add_argument("--training-inputs", default="baseV3E + centrality + weta33/wphi33")
    args = ap.parse_args()
    if args.mode == "extract" and (args.manifest is None or args.json_out is None):
        ap.error("--manifest and --json-out are required for --mode extract")
    if args.mode == "render" and (args.input is None or args.outdir is None):
        ap.error("--input and --outdir are required for --mode render")
    return args


def main() -> None:
    args = parse_args()
    if args.mode == "extract":
        payload = extract(args)
        text = "__THE41_JSON_BEGIN__\n" + json.dumps(payload, separators=(",", ":")) + "\n__THE41_JSON_END__\n"
        if str(args.json_out) == "-":
            print(text, end="")
        else:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            print(args.json_out)
    else:
        manifest = render(args)
        print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
