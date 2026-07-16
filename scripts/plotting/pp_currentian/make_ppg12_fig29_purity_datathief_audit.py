#!/usr/bin/env python3
"""DataThief audit for PPG12 IAN Fig. 29 / paper Fig. 5 purity points.

The SDCC source contract is the PPG12 nominal final-yield ROOT file:

    /sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root

with `gpurity` and `gpurity_leak`. This script compares those graph points to
the visible marker layer from the published Fig. 5 PNG, using the local
DataThief jar for the image-coordinate transform.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

from make_ppg12_datathief_validation_overlays import (
    DATATHIEF_JAR,
    jar_md5,
    load_datathief_csv,
    run_datathief_export,
)


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024"
    / "fig29_purity_datathief_audit"
)

PAPER_FIG5 = REPO / "usefulDocs/ppg12/sPH-JET-2026-02_public_figures_20260620/sPH-JET-2026-02-Fig5.png"
IAN_FIG29_CROP = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_current_ian_fig29_left_purity.png"
)
SDCC_EXTRACT_CSV = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_ratio_diagnostic/ppg12_photon_final_bdt_nom_extract.csv"
)

PPG12_SOURCE_ROOT = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root"
PPG12_PLOTTING_MACRO = "ppg12codeGit/plotting/plot_purity_selection.C"
DATATHIEF_INTERNAL_Y_SCALE = 20.0


def load_sdcc_points() -> dict[str, dict[str, np.ndarray]]:
    out: dict[str, list[dict[str, float]]] = {"raw": [], "leakage_corrected": []}
    name_to_series = {"gpurity": "raw", "gpurity_leak": "leakage_corrected"}
    with SDCC_EXTRACT_CSV.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("kind") != "GRAPH" or row.get("source") != "data":
                continue
            series = name_to_series.get(row.get("name", ""))
            if series is None:
                continue
            out[series].append(
                {
                    "x": float(row["x"]),
                    "y": float(row["y"]),
                    "ex_low": float(row["ex_low"]),
                    "ex_high": float(row["ex_high"]),
                    "ey_low": float(row["ey_low"]),
                    "ey_high": float(row["ey_high"]),
                }
            )
    converted: dict[str, dict[str, np.ndarray]] = {}
    for series, rows in out.items():
        rows = sorted(rows, key=lambda r: r["x"])
        converted[series] = {k: np.asarray([row[k] for row in rows], dtype=float) for k in rows[0]}
    return converted


def grouped_centers(indices: np.ndarray, values: np.ndarray) -> list[tuple[float, float, float]]:
    if len(indices) == 0:
        return []
    groups: list[tuple[int, int]] = []
    start = prev = int(indices[0])
    for raw_idx in indices[1:]:
        idx = int(raw_idx)
        if idx == prev + 1:
            prev = idx
        else:
            groups.append((start, prev))
            start = prev = idx
    groups.append((start, prev))
    return [(0.5 * (lo + hi), lo, hi) for lo, hi in groups]


def choose_evenly_spaced_ticks(
    candidates: list[float],
    *,
    n_ticks: int,
    force_first: float | None = None,
    force_last: float | None = None,
) -> list[float]:
    """Choose an evenly spaced major-tick sequence from noisy tick candidates."""
    unique = sorted({round(float(c), 3) for c in candidates})
    if force_first is not None and force_first not in unique:
        unique.append(force_first)
    if force_last is not None and force_last not in unique:
        unique.append(force_last)
    unique = sorted(unique)

    required = {force_first, force_last} - {None}
    optional = [c for c in unique if c not in required]
    n_optional = n_ticks - len(required)
    best: tuple[float, list[float]] | None = None
    for combo in itertools.combinations(optional, n_optional):
        ticks = sorted(list(required) + list(combo))
        if len(ticks) != n_ticks:
            continue
        idx = np.arange(n_ticks, dtype=float)
        vals = np.asarray(ticks, dtype=float)
        coeff = np.polyfit(idx, vals, 1)
        pred = coeff[0] * idx + coeff[1]
        rms = float(np.sqrt(np.mean((vals - pred) ** 2)))
        # Strongly prefer sequences with the intended endpoints.
        if force_first is not None and not math.isclose(ticks[0], force_first, abs_tol=2.0):
            rms += 1000.0
        if force_last is not None and not math.isclose(ticks[-1], force_last, abs_tol=2.0):
            rms += 1000.0
        if best is None or rms < best[0]:
            best = (rms, ticks)
    if best is None:
        raise RuntimeError(f"Could not choose {n_ticks} evenly spaced ticks from {candidates}")
    return best[1]


def axis_from_visible_ticks(image: Path) -> dict[str, object]:
    """Infer the paper Fig. 5 axis references from visible major ticks."""
    arr = np.asarray(Image.open(image).convert("RGB"))
    gray = np.dot(arr[..., :3], [0.299, 0.587, 0.114])
    dark = gray < 70

    # The full-height frame columns and full-width frame rows are unambiguous
    # in the public figure PNG. They also bound the tick search windows.
    col_counts = dark.sum(axis=0)
    row_counts = dark.sum(axis=1)
    x_frame_candidates = np.where(col_counts > 250)[0]
    y_frame_candidates = np.where(row_counts > 250)[0]
    if len(x_frame_candidates) < 2 or len(y_frame_candidates) < 2:
        raise RuntimeError(f"Could not identify frame lines in {image}")
    x_left = float(x_frame_candidates[0])
    x_right = float(x_frame_candidates[-1])
    y_top_frame = float(y_frame_candidates[0])
    y_bottom = float(y_frame_candidates[-1])

    bottom_strip = dark[int(y_bottom) - 35 : int(y_bottom) + 1, :].sum(axis=0)
    x_groups = grouped_centers(
        np.where(
            (bottom_strip > 25)
            & (np.arange(len(bottom_strip)) >= x_left)
            & (np.arange(len(bottom_strip)) <= x_right)
        )[0],
        bottom_strip,
    )
    x_major_candidates = [center for center, lo, hi in x_groups if hi - lo <= 4]
    try:
        x_major = choose_evenly_spaced_ticks(
            x_major_candidates,
            n_ticks=6,
            force_first=float(x_major_candidates[0]) if x_major_candidates else 164.5,
            force_last=float(x_major_candidates[-1]) if x_major_candidates else 1098.0,
        )
    except RuntimeError:
        # Stable fallback from the public Fig. 5 PNG; retained in the manifest.
        x_major = [164.5, 351.5, 538.0, 724.5, 911.0, 1098.0]

    left_strip = dark[:, int(x_left) : int(x_left) + 32].sum(axis=1)
    y_groups = grouped_centers(
        np.where(
            (left_strip > 25)
            & (np.arange(len(left_strip)) >= y_top_frame)
            & (np.arange(len(left_strip)) <= y_bottom)
        )[0],
        left_strip,
    )
    # Major y ticks are 0.0, 0.2, ..., 1.0. The top frame is above 1.0 and is
    # not used as the y-axis value reference.
    y_major_candidates = [center for center, lo, hi in y_groups if hi - lo <= 4 and center > 50]
    try:
        y_major_ascending = choose_evenly_spaced_ticks(
            y_major_candidates,
            n_ticks=6,
            force_first=float(min(y_major_candidates)) if y_major_candidates else 138.5,
            force_last=float(max(y_major_candidates)) if y_major_candidates else 904.0,
        )
        y_major = sorted(y_major_ascending, reverse=True)
    except RuntimeError:
        y_major = [904.0, 751.0, 598.0, 445.0, 292.0, 138.5]

    return {
        "x_10": float(x_major[0]),
        "x_15": float(x_major[1]),
        "x_20": float(x_major[2]),
        "x_25": float(x_major[3]),
        "x_30": float(x_major[4]),
        "x_35": float(x_major[5]),
        "y_0": float(y_major[0]),
        "y_0p2": float(y_major[1]),
        "y_0p4": float(y_major[2]),
        "y_0p6": float(y_major[3]),
        "y_0p8": float(y_major[4]),
        "y_1p0": float(y_major[5]),
        "y_top_frame": y_top_frame,
        "x_left_frame": x_left,
        "x_right_frame": x_right,
    }


def x_to_pixel(x: float, axis: dict[str, object]) -> float:
    return float(axis["x_10"]) + (x - 10.0) / 25.0 * (float(axis["x_35"]) - float(axis["x_10"]))


def y_to_pixel(y: float, axis: dict[str, object]) -> float:
    return float(axis["y_0"]) - y * (float(axis["y_0"]) - float(axis["y_1p0"]))


def pixel_to_y(py: float, axis: dict[str, object]) -> float:
    return (float(axis["y_0"]) - py) / (float(axis["y_0"]) - float(axis["y_1p0"]))


def detect_leakage_points(
    arr: np.ndarray,
    axis: dict[str, object],
    sdcc: dict[str, np.ndarray],
) -> list[tuple[str, float, float]]:
    """Find filled blue leakage-corrected markers in the paper Fig. 5 image."""
    blue = (arr[:, :, 2] > 135) & (arr[:, :, 0] < 90) & (arr[:, :, 1] < 160)
    points: list[tuple[str, float, float]] = []
    for x, y in zip(sdcc["x"], sdcc["y"]):
        xp = x_to_pixel(float(x), axis)
        yp = y_to_pixel(float(y), axis)
        best: tuple[float, float, float] | None = None
        for cy in range(int(round(yp)) - 34, int(round(yp)) + 35):
            for cx in range(int(round(xp)) - 24, int(round(xp)) + 25):
                y0 = max(0, cy - 17)
                y1 = min(blue.shape[0], cy + 18)
                x0 = max(0, cx - 17)
                x1 = min(blue.shape[1], cx + 18)
                yy, xx = np.ogrid[y0:y1, x0:x1]
                rr = (xx - cx) ** 2 + (yy - cy) ** 2
                patch = blue[y0:y1, x0:x1]
                disk = rr <= 10**2
                annulus = (rr > 11**2) & (rr <= 17**2)
                score = (
                    float(patch[disk].sum())
                    - 0.12 * float(patch[annulus].sum())
                    - 0.05 * abs(cx - xp)
                    - 0.03 * abs(cy - yp)
                )
                if best is None or score > best[0]:
                    best = (score, float(cx), float(cy))
        if best is None:
            raise RuntimeError(f"Could not find leakage-corrected marker near E_T={x}")
        points.append(("leakage_corrected", best[1], best[2]))
    return points


def detect_raw_points(
    arr: np.ndarray,
    axis: dict[str, object],
    sdcc: dict[str, np.ndarray],
) -> list[tuple[str, float, float]]:
    """Find raw open-circle marker rows in the paper Fig. 5 image.

    The raw markers are open black circles drawn on top of black error bars.
    The most stable center estimator is the horizontal error-bar row at the
    known PPG12 bin center.
    """
    black = (arr[:, :, 0] < 80) & (arr[:, :, 1] < 80) & (arr[:, :, 2] < 80)
    points: list[tuple[str, float, float]] = []
    for x, y in zip(sdcc["x"], sdcc["y"]):
        xp = x_to_pixel(float(x), axis)
        yp = y_to_pixel(float(y), axis)
        ylo = max(0, int(round(yp)) - 45)
        yhi = min(black.shape[0] - 1, int(round(yp)) + 45)
        xlo = max(0, int(round(xp)) - 35)
        xhi = min(black.shape[1] - 1, int(round(xp)) + 35)
        counts = black[ylo : yhi + 1, xlo : xhi + 1].sum(axis=1).astype(float)
        smoothed = np.convolve(counts, np.ones(3) / 3.0, mode="same")
        candidate_rows = np.arange(ylo, yhi + 1, dtype=float)
        score = smoothed - 0.12 * np.abs(candidate_rows - yp)
        cy = float(candidate_rows[int(np.argmax(score))])
        # Use the known binned x-position for assignment; the y-coordinate is
        # what this audit is checking against the SDCC ROOT graph.
        points.append(("raw", xp, cy))
    return points


def run_datathief(points: list[tuple[str, float, float]], axis: dict[str, object]) -> Path:
    raw_csv = OUTDIR / "ppg12_fig5_purity_marker_pixels_datathief_raw_export.csv"
    # DataThief's legacy parser is more reliable with a normalized x-axis and a
    # larger y-axis span. The exported values are converted back after export.
    refs = [
        (0, float(axis["x_10"]), float(axis["y_0"]), 0.0, 0.0),
        (1, float(axis["x_35"]), float(axis["y_0"]), 1.0, 0.0),
        (2, float(axis["x_10"]), float(axis["y_1p0"]), 0.0, DATATHIEF_INTERNAL_Y_SCALE),
    ]
    run_datathief_export(PAPER_FIG5, refs, points, raw_csv)
    return raw_csv


def write_local_tick_transform(points: list[tuple[str, float, float]], axis: dict[str, object]) -> Path:
    """Fallback for the same 3-reference affine transform DataThief should do."""
    out = OUTDIR / "ppg12_fig5_purity_marker_pixels_local_tick_transform.csv"
    x10 = float(axis["x_10"])
    x35 = float(axis["x_35"])
    y0 = float(axis["y_0"])
    y1 = float(axis["y_1p0"])
    with out.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "point_index", "datathief_x", "datathief_y"])
        counts: dict[str, int] = {}
        for series, px, py in points:
            idx = counts.get(series, 0)
            counts[series] = idx + 1
            xval = 10.0 + 25.0 * (px - x10) / (x35 - x10)
            yval = (y0 - py) / (y0 - y1)
            writer.writerow([series, idx, f"{xval:.12g}", f"{yval:.12g}"])
    return out


def load_valid_digitized_points(
    pixel_points: list[tuple[str, float, float]],
    axis: dict[str, object],
) -> tuple[Path, dict[str, list[tuple[float, float]]], dict[str, object]]:
    status: dict[str, object] = {
        "datathief_jar_attempted": True,
        "datathief_internal_y_scale": DATATHIEF_INTERNAL_Y_SCALE,
    }
    try:
        raw_dt_csv = run_datathief(pixel_points, axis)
        dt_norm = load_datathief_csv(raw_dt_csv)
        dt = {
            series: [(10.0 + 25.0 * xv, yv / DATATHIEF_INTERNAL_Y_SCALE) for xv, yv in rows]
            for series, rows in dt_norm.items()
        }
        flat = [coord for rows in dt.values() for pair in rows for coord in pair]
        if not flat or not np.all(np.isfinite(np.asarray(flat, dtype=float))):
            raise RuntimeError(f"DataThief jar returned non-finite coordinates in {raw_dt_csv}")
        expected_csv = write_local_tick_transform(pixel_points, axis)
        expected = load_datathief_csv(expected_csv)
        for series in ("raw", "leakage_corrected"):
            if len(dt.get(series, [])) != len(expected.get(series, [])):
                raise RuntimeError(
                    f"DataThief jar returned {len(dt.get(series, []))} {series} points; "
                    f"expected {len(expected.get(series, []))}"
                )
            for idx, ((xv, yv), (x_exp, y_exp)) in enumerate(zip(dt[series], expected[series])):
                if abs(xv - x_exp) > 0.1 or abs(yv - y_exp) > 0.02:
                    raise RuntimeError(
                        f"DataThief jar coordinate sanity check failed for {series} point {idx}: "
                        f"got ({xv:.6g}, {yv:.6g}), expected approximately "
                        f"({x_exp:.6g}, {y_exp:.6g})"
                    )
        status.update(
            {
                "coordinate_transform_method": "official_datathief_jar",
                "datathief_jar_status": "success",
            }
        )
        return raw_dt_csv, dt, status
    except Exception as exc:
        fallback_csv = write_local_tick_transform(pixel_points, axis)
        dt = load_datathief_csv(fallback_csv)
        status.update(
            {
                "coordinate_transform_method": "local_tick_affine_transform_fallback",
                "datathief_jar_status": "failed",
                "datathief_jar_error": repr(exc),
                "fallback_note": (
                    "The local transform uses the same three axis reference points "
                    "that were supplied to DataThief: x=10 at left, x=35 at right, "
                    "and y=1.0 at the visible major tick. The official jar attempt "
                    "is preserved but not used if it fails, returns NaN, or returns "
                    "coordinates that fail the reference-point sanity check."
                ),
            }
        )
        return fallback_csv, dt, status


def write_tables(
    sdcc: dict[str, dict[str, np.ndarray]],
    dt: dict[str, list[tuple[float, float]]],
    pixel_points: list[tuple[str, float, float]],
) -> tuple[Path, Path, Path]:
    sdcc_csv = OUTDIR / "ppg12_fig29_purity_sdcc_points.csv"
    dt_csv = OUTDIR / "ppg12_fig5_purity_datathief_points.csv"
    ratio_csv = OUTDIR / "ppg12_fig5_datathief_over_sdcc_purity_ratios.csv"

    with sdcc_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "point_index", "x", "purity", "ex_low", "ex_high", "ey_low", "ey_high"])
        for series in ("raw", "leakage_corrected"):
            n = len(sdcc[series]["x"])
            for i in range(n):
                writer.writerow(
                    [
                        series,
                        i,
                        f"{sdcc[series]['x'][i]:.12g}",
                        f"{sdcc[series]['y'][i]:.12g}",
                        f"{sdcc[series]['ex_low'][i]:.12g}",
                        f"{sdcc[series]['ex_high'][i]:.12g}",
                        f"{sdcc[series]['ey_low'][i]:.12g}",
                        f"{sdcc[series]['ey_high'][i]:.12g}",
                    ]
                )

    pixel_by_series: dict[str, list[tuple[float, float]]] = {"raw": [], "leakage_corrected": []}
    for series, px, py in pixel_points:
        pixel_by_series[series].append((px, py))

    with dt_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "point_index", "datathief_x", "datathief_purity", "pixel_x", "pixel_y"])
        for series in ("raw", "leakage_corrected"):
            for i, ((xv, yv), (px, py)) in enumerate(zip(dt[series], pixel_by_series[series])):
                writer.writerow([series, i, f"{xv:.12g}", f"{yv:.12g}", f"{px:.6f}", f"{py:.6f}"])

    with ratio_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "series",
                "point_index",
                "x_sdcc",
                "x_datathief",
                "purity_sdcc",
                "purity_datathief",
                "datathief_over_sdcc",
                "delta_abs",
            ]
        )
        for series in ("raw", "leakage_corrected"):
            for i, (xv, yv) in enumerate(dt[series]):
                y_sdcc = float(sdcc[series]["y"][i])
                ratio = yv / y_sdcc if y_sdcc else math.nan
                writer.writerow(
                    [
                        series,
                        i,
                        f"{sdcc[series]['x'][i]:.12g}",
                        f"{xv:.12g}",
                        f"{y_sdcc:.12g}",
                        f"{yv:.12g}",
                        f"{ratio:.12g}",
                        f"{(yv - y_sdcc):.12g}",
                    ]
                )
    return sdcc_csv, dt_csv, ratio_csv


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
        }
    )


def plot_overlay(sdcc: dict[str, dict[str, np.ndarray]], dt: dict[str, list[tuple[float, float]]]) -> Path:
    setup_style()
    out = OUTDIR / "ppg12_fig5_fig29_purity_datathief_vs_sdcc_overlay_ratio.png"
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.2, 7.8),
        dpi=200,
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.04},
    )

    colors = {"raw": "black", "leakage_corrected": "#1267c9"}
    labels = {"raw": "w/o sig. leak. corr.", "leakage_corrected": "w/ sig. leak. corr."}
    for series in ("leakage_corrected", "raw"):
        ax.errorbar(
            sdcc[series]["x"],
            sdcc[series]["y"],
            xerr=[sdcc[series]["ex_low"], sdcc[series]["ex_high"]],
            yerr=[sdcc[series]["ey_low"], sdcc[series]["ey_high"]],
            fmt="o",
            ms=5.2,
            lw=1.1,
            color=colors[series],
            mfc=colors[series],
            mec=colors[series],
            label=f"SDCC {labels[series]}",
            zorder=3,
        )
        x_dt = np.asarray([p[0] for p in dt[series]], dtype=float)
        y_dt = np.asarray([p[1] for p in dt[series]], dtype=float)
        ax.plot(
            x_dt,
            y_dt,
            "o",
            ms=8.0,
            mfc="none",
            mec=colors[series],
            mew=1.6,
            label=f"DataThief {labels[series]}",
            zorder=4,
        )
        ratio = y_dt / sdcc[series]["y"]
        rax.plot(
            sdcc[series]["x"],
            ratio,
            "o",
            ms=5.2,
            mfc="none",
            mec=colors[series],
            mew=1.2,
            label=labels[series],
            zorder=3,
        )

    ax.text(0.07, 0.92, "sPHENIX", transform=ax.transAxes, fontsize=17, fontstyle="italic", fontweight="bold")
    ax.text(0.285, 0.92, "Internal", transform=ax.transAxes, fontsize=17)
    ax.text(0.07, 0.845, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=14.5)
    ax.text(0.07, 0.775, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=14.5)
    ax.text(0.07, 0.705, "PPG12 paper Fig. 5 / IAN Fig. 29", transform=ax.transAxes, fontsize=12.5)
    ax.set_ylabel("Purity", fontsize=18)
    ax.set_xlim(10, 35)
    ax.set_ylim(0.0, 1.15)
    ax.minorticks_on()
    ax.tick_params(which="both", top=True, right=True, labelsize=13)

    source_handles = [
        Line2D([0], [0], marker="o", color="black", mfc="black", mec="black", linestyle="None", ms=6, label="PPG12 SDCC ROOT"),
        Line2D([0], [0], marker="o", color="black", mfc="none", mec="black", linestyle="None", ms=6, label="Fig. 5 digitized markers"),
    ]
    series_handles = [
        Line2D([0], [0], marker="o", color=colors["leakage_corrected"], linestyle="-", ms=6, label="w/ sig. leak. corr."),
        Line2D([0], [0], marker="o", color=colors["raw"], linestyle="-", ms=6, label="w/o sig. leak. corr."),
    ]
    leg1 = ax.legend(handles=series_handles, loc="lower left", bbox_to_anchor=(0.10, 0.15), frameon=False, fontsize=12.5)
    ax.add_artist(leg1)
    ax.legend(handles=source_handles, loc="lower left", bbox_to_anchor=(0.10, 0.03), frameon=False, fontsize=11.5)

    rax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylabel("DataThief / SDCC", fontsize=14)
    rax.set_xlabel(r"$E_T^\gamma$ [GeV]", fontsize=18)
    rax.set_ylim(0.985, 1.015)
    rax.minorticks_on()
    rax.tick_params(which="both", top=True, right=True, labelsize=13)
    rax.legend(loc="lower left", frameon=False, fontsize=11.5, ncol=2, handlelength=1.2)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.97, bottom=0.10)
    fig.savefig(out)
    plt.close(fig)
    return out


def write_manifest(
    *,
    axis: dict[str, object],
    raw_datathief_csv: Path,
    sdcc_csv: Path,
    dt_csv: Path,
    ratio_csv: Path,
    png: Path,
    transform_status: dict[str, object],
) -> Path:
    ratios: dict[str, list[float]] = {"raw": [], "leakage_corrected": []}
    with ratio_csv.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            ratios[row["series"]].append(float(row["datathief_over_sdcc"]))
    summary = {
        series: {
            "n": len(vals),
            "mean_ratio": float(np.nanmean(np.asarray(vals, dtype=float))),
            "max_abs_ratio_minus_one": float(np.nanmax(np.abs(np.asarray(vals, dtype=float) - 1.0))),
        }
        for series, vals in ratios.items()
    }
    manifest = {
        "artifact": str(png),
        "comparison": "PPG12 visible Fig. 5/Fig. 29 purity markers exported through DataThief vs PPG12 SDCC final ROOT graph values",
        "ppg12_source_root": PPG12_SOURCE_ROOT,
        "ppg12_objects": ["gpurity", "gpurity_leak"],
        "ppg12_plotting_macro": PPG12_PLOTTING_MACRO,
        "source_paper_fig5_png": str(PAPER_FIG5),
        "source_ian_fig29_crop": str(IAN_FIG29_CROP),
        "local_sdcc_graph_extract_csv": str(SDCC_EXTRACT_CSV),
        "output_sdcc_points_csv": str(sdcc_csv),
        "output_datathief_points_csv": str(dt_csv),
        "output_ratio_csv": str(ratio_csv),
        "raw_datathief_export_csv": str(raw_datathief_csv),
        "datathief_jar": str(DATATHIEF_JAR),
        "datathief_jar_md5": jar_md5(),
        "transform_status": transform_status,
        "axis_reference_pixels": axis,
        "axis_reference_values": [
            {"pixel": [axis["x_10"], axis["y_0"]], "value": [10.0, 0.0]},
            {"pixel": [axis["x_35"], axis["y_0"]], "value": [35.0, 0.0]},
            {"pixel": [axis["x_10"], axis["y_1p0"]], "value": [10.0, 1.0]},
        ],
        "point_assignment_note": (
            "Known PPG12 bin centers are used to assign each visible marker. "
            "Leakage-corrected y points use filled-blue marker template centers; "
            "raw y points use the open-circle horizontal error-bar row, because "
            "the open marker shares black pixels with its error bars. The final "
            "data values are the DataThief axis transform of those pixels."
        ),
        "ratio_summary": summary,
    }
    out = OUTDIR / "ppg12_fig5_fig29_purity_datathief_vs_sdcc_manifest.json"
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if not PAPER_FIG5.exists():
        raise FileNotFoundError(PAPER_FIG5)
    if not SDCC_EXTRACT_CSV.exists():
        raise FileNotFoundError(SDCC_EXTRACT_CSV)

    shutil.copy2(PAPER_FIG5, OUTDIR / "source_paper_fig5_used_for_datathief.png")
    if IAN_FIG29_CROP.exists():
        shutil.copy2(IAN_FIG29_CROP, OUTDIR / "source_ian_fig29_crop_reference.png")

    sdcc = load_sdcc_points()
    axis = axis_from_visible_ticks(PAPER_FIG5)
    arr = np.asarray(Image.open(PAPER_FIG5).convert("RGB"))
    pixel_points = detect_raw_points(arr, axis, sdcc["raw"])
    pixel_points += detect_leakage_points(arr, axis, sdcc["leakage_corrected"])
    raw_dt_csv, dt, transform_status = load_valid_digitized_points(pixel_points, axis)
    sdcc_csv, dt_csv, ratio_csv = write_tables(sdcc, dt, pixel_points)
    png = plot_overlay(sdcc, dt)
    manifest = write_manifest(
        axis=axis,
        raw_datathief_csv=raw_dt_csv,
        sdcc_csv=sdcc_csv,
        dt_csv=dt_csv,
        ratio_csv=ratio_csv,
        png=png,
        transform_status=transform_status,
    )
    print(f"Wrote {png}")
    print(f"Wrote {ratio_csv}")
    print(f"Wrote {manifest}")


if __name__ == "__main__":
    main()
