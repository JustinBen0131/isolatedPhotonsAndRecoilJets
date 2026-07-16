#!/usr/bin/env python3
"""Strict local THE-76 PPG12 SIM parity artifacts.

Inputs are deliberately separated:
  * exact IAN raster crops, digitized through the local DataThief jar;
  * live SDCC PPG12 source CSVs extracted from the source ROOT objects;
  * July 1 RecoilJets SIM aggregate CSVs extracted from raw campaign ROOTs.

This script does not submit jobs, merge ROOT files, edit SDCC state, or mutate
Google Slides. It only writes local PNG/CSV/JSON artifacts.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

from make_ppg12_datathief_validation_overlays import (  # noqa: E402
    DATATHIEF_JAR,
    jar_md5,
    load_datathief_csv,
    run_datathief_export,
)


TAG = "the76_ppg12_parity_full_20260701_003024"
OUT = REPO / "dataOutput/ppg12Parity" / TAG
SOURCE_CSV = OUT / "reference_audit/sdcc_pull/live_ppg12_ian_source_points.csv"
CURRENT_CSV = OUT / "july1_sim_aggregate/july1_strict_stitch_points.csv"
JET_COMPONENT_CSV = OUT / "july1_sim_aggregate/july1_jet_component_stitch_points.csv"
TEFF_CSV = OUT / "strict_iso_efficiency/ppg12_teff_reference_points.csv"

FIG5_CROP = OUT / "reference_audit/exact_ian_crops/fig5_photon_exact_ian_pad.png"
FIG6_CROP = OUT / "reference_audit/exact_ian_crops/fig6_inclusivejet_exact_ian_pad.png"


COLORS = {
    "photon5": "#e7298a",
    "photon10": "#1b9e37",
    "photon20": "#1296f3",
    "jet8": "#e7298a",
    "jet12": "#1b9e37",
    "jet20": "#1296f3",
    "jet30": "#ff6f00",
    "jet40": "#e7298a",
}


@dataclass(frozen=True)
class AxisSpec:
    image: Path
    group: str
    samples: tuple[str, ...]
    xlim: tuple[float, float]
    top_xlim: tuple[float, float]
    top_ylim: tuple[float, float]
    top_x_refs: tuple[tuple[float, float], tuple[float, float]]
    top_log_y_refs: tuple[tuple[float, float], tuple[float, float]]
    ratio_x_refs: tuple[tuple[float, float], tuple[float, float]]
    ratio_y_refs: tuple[tuple[float, float], tuple[float, float]]
    xlabel: str
    ylabel: str
    title: str
    source_note: str


SPECS = {
    "photon": AxisSpec(
        image=FIG5_CROP,
        group="photon",
        samples=("photon5", "photon10", "photon20"),
        xlim=(10.0, 40.0),
        top_xlim=(10.0, 40.0),
        top_ylim=(3.0e-2, 2.0e4),
        top_x_refs=((98.0, 10.0), (476.0, 40.0)),
        top_log_y_refs=((296.0, -1.0), (50.0, 4.0)),
        ratio_x_refs=((98.0, 10.0), (476.0, 40.0)),
        ratio_y_refs=((452.0, 0.85), (301.0, 1.15)),
        xlabel=r"Leading $E_T^\gamma$ [GeV]",
        ylabel=r"$d\sigma/dE_T^\gamma$ [pb / GeV]",
        title="PPG12 IAN Fig. 5 photon+jet stitch",
        source_note="Fig.5 source: photon_max_pT_uncut.root, PPG12 no-bin-width photon stitch convention",
    ),
    "jet": AxisSpec(
        image=FIG6_CROP,
        group="jet",
        samples=("jet8", "jet12", "jet20", "jet30", "jet40"),
        xlim=(9.0, 50.0),
        top_xlim=(9.0, 50.0),
        top_ylim=(1.0e4, 2.0e12),
        top_x_refs=((108.0, 9.0), (486.0, 50.0)),
        top_log_y_refs=((309.0, 4.0), (9.0, 12.0)),
        ratio_x_refs=((108.0, 9.0), (486.0, 50.0)),
        ratio_y_refs=((459.0, 0.85), (309.0, 1.15)),
        xlabel=r"Leading $p_T^\mathrm{jet}$ [GeV]",
        ylabel="counts",
        title="PPG12 IAN Fig. 6 inclusive-jet stitch",
        source_note="Fig.6 source: MC_efficiency_jet*_bdt_nom.root h_max_truth_jet_pT, 1 GeV rebinned counts",
    ),
}


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def read_rows(path: Path) -> list[dict[str, str]]:
    require(path)
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def hex_rgb(color: str) -> np.ndarray:
    color = color.lstrip("#")
    return np.asarray([int(color[i : i + 2], 16) for i in (0, 2, 4)], dtype=float)


def data_to_pixel_linear(
    x: float,
    y: float,
    x_refs: tuple[tuple[float, float], tuple[float, float]],
    y_refs: tuple[tuple[float, float], tuple[float, float]],
) -> tuple[float, float]:
    (px0, x0), (px1, x1) = x_refs
    (py0, y0), (py1, y1) = y_refs
    px = px0 + (x - x0) / (x1 - x0) * (px1 - px0)
    py = py0 + (y - y0) / (y1 - y0) * (py1 - py0)
    return px, py


def data_to_pixel_logy(spec: AxisSpec, x: float, y: float) -> tuple[float, float]:
    return data_to_pixel_linear(x, math.log10(y), spec.top_x_refs, spec.top_log_y_refs)


def find_colored_marker_pixels(spec: AxisSpec, rows: list[dict[str, str]]) -> list[dict[str, object]]:
    arr = np.asarray(Image.open(spec.image).convert("RGB"), dtype=float)
    out: list[dict[str, object]] = []
    for row in rows:
        sample = row["sample"]
        x = float(row["bin_center"])
        y = float(row["value"])
        if sample not in spec.samples or y <= 0 or not (spec.xlim[0] <= x <= spec.xlim[1]):
            continue
        px0, py0 = data_to_pixel_logy(spec, x, y)
        target = hex_rgb(COLORS[sample])
        dist = np.sqrt(np.sum((arr - target) ** 2, axis=2))
        color_mask = dist < 115.0
        xlo = max(0, int(round(px0)) - 10)
        xhi = min(arr.shape[1] - 1, int(round(px0)) + 10)
        ylo = max(0, int(round(py0)) - 10)
        yhi = min(arr.shape[0] - 1, int(round(py0)) + 10)
        yy, xx = np.where(color_mask[ylo : yhi + 1, xlo : xhi + 1])
        mode = "sample_color"
        if len(xx) < 2:
            # The Fig. 6 fit line overlaps some markers. Fall back to any
            # saturated colored pixels near the predicted location.
            xlo = max(0, int(round(px0)) - 16)
            xhi = min(arr.shape[1] - 1, int(round(px0)) + 16)
            ylo = max(0, int(round(py0)) - 16)
            yhi = min(arr.shape[0] - 1, int(round(py0)) + 16)
            crop = arr[ylo : yhi + 1, xlo : xhi + 1]
            sat = (crop.max(axis=2) - crop.min(axis=2) > 40.0) & (crop.max(axis=2) > 110.0)
            yy, xx = np.where(sat)
            mode = "saturated_local"
        if len(xx) < 2:
            continue
        abs_x = xlo + xx.astype(float)
        abs_y = ylo + yy.astype(float)
        r2 = (abs_x - px0) ** 2 + (abs_y - py0) ** 2
        weights = np.exp(-0.5 * r2 / (3.0**2))
        out.append(
            {
                **row,
                "pixel_x": float(np.average(abs_x, weights=weights)),
                "pixel_y": float(np.average(abs_y, weights=weights)),
                "pred_pixel_x": float(px0),
                "pred_pixel_y": float(py0),
                "pixels_used": int(len(xx)),
                "match_mode": mode,
            }
        )
    return out


def connected_components(mask: np.ndarray) -> list[dict[str, float]]:
    h, w = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    out: list[dict[str, float]] = []
    for y in range(h):
        xs = np.where(mask[y] & (~seen[y]))[0]
        for x0 in xs:
            x0 = int(x0)
            if seen[y, x0] or not mask[y, x0]:
                continue
            stack = [(y, x0)]
            seen[y, x0] = True
            pts: list[tuple[int, int]] = []
            while stack:
                cy, cx = stack.pop()
                pts.append((cy, cx))
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ny = cy + dy
                    nx = cx + dx
                    if 0 <= ny < h and 0 <= nx < w and mask[ny, nx] and not seen[ny, nx]:
                        seen[ny, nx] = True
                        stack.append((ny, nx))
            if len(pts) < 3:
                continue
            yy = np.asarray([p[0] for p in pts], dtype=float)
            xx = np.asarray([p[1] for p in pts], dtype=float)
            out.append(
                {
                    "area": float(len(pts)),
                    "x": float(xx.mean()),
                    "y": float(yy.mean()),
                    "xmin": float(xx.min()),
                    "xmax": float(xx.max()),
                    "ymin": float(yy.min()),
                    "ymax": float(yy.max()),
                    "width": float(xx.max() - xx.min() + 1),
                    "height": float(yy.max() - yy.min() + 1),
                }
            )
    return out


def find_black_ratio_marker_pixels(spec: AxisSpec, rows: list[dict[str, str]]) -> list[dict[str, object]]:
    arr = np.asarray(Image.open(spec.image).convert("RGB"), dtype=float)
    black = (arr[:, :, 0] < 85) & (arr[:, :, 1] < 85) & (arr[:, :, 2] < 85)
    panel_xlo = int(round(min(spec.ratio_x_refs[0][0], spec.ratio_x_refs[1][0]))) + 3
    panel_xhi = int(round(max(spec.ratio_x_refs[0][0], spec.ratio_x_refs[1][0]))) - 3
    panel_ylo = int(round(min(spec.ratio_y_refs[0][0], spec.ratio_y_refs[1][0]))) + 4
    panel_yhi = int(round(max(spec.ratio_y_refs[0][0], spec.ratio_y_refs[1][0]))) - 4
    black_filtered = np.zeros_like(black)
    black_filtered[panel_ylo : panel_yhi + 1, panel_xlo : panel_xhi + 1] = black[
        panel_ylo : panel_yhi + 1, panel_xlo : panel_xhi + 1
    ]
    out: list[dict[str, object]] = []
    for row in rows:
        sample = row["sample"]
        x = float(row["bin_center"])
        ratio = float(row["fit_ratio"])
        if sample not in spec.samples or ratio <= 0 or not (spec.xlim[0] <= x <= spec.xlim[1]):
            continue
        if not (0.82 <= ratio <= 1.18):
            continue
        px0, py0 = data_to_pixel_linear(x, ratio, spec.ratio_x_refs, spec.ratio_y_refs)
        xlo = max(panel_xlo, int(round(px0)) - 8)
        xhi = min(panel_xhi, int(round(px0)) + 8)
        ylo = panel_ylo
        yhi = panel_yhi
        local = black_filtered[ylo : yhi + 1, xlo : xhi + 1]
        components = []
        for comp in connected_components(local):
            if not (4 <= comp["area"] <= 220):
                continue
            if not (2 <= comp["width"] <= 15 and 2 <= comp["height"] <= 24):
                continue
            abs_x = xlo + comp["x"]
            abs_y = ylo + comp["y"]
            score = abs(abs_x - px0) + 0.08 * abs(abs_y - py0)
            components.append((score, abs_x, abs_y, comp))
        if not components:
            continue
        components.sort(key=lambda item: item[0])
        _, abs_x, abs_y, comp = components[0]
        out.append(
            {
                **row,
                "pixel_x": float(abs_x),
                "pixel_y": float(abs_y),
                "pred_pixel_x": float(px0),
                "pred_pixel_y": float(py0),
                "pixels_used": int(comp["area"]),
                "match_mode": "connected_black_marker_component",
            }
        )
    return out


def finite_datathief(raw: dict[str, list[tuple[float, float]]], expected: int) -> bool:
    pts = [pt for series in raw.values() for pt in series]
    if len(pts) != expected:
        return False
    return all(math.isfinite(x) and math.isfinite(y) for x, y in pts)


def export_log_datathief(spec: AxisSpec, pixel_rows: list[dict[str, object]], out_csv: Path) -> dict[str, list[tuple[float, float]]]:
    out_csv.unlink(missing_ok=True)
    if out_csv.exists():
        raw = load_datathief_csv(out_csv)
        if sum(len(v) for v in raw.values()) == len(pixel_rows):
            converted: dict[str, list[tuple[float, float]]] = {}
            log0 = spec.top_log_y_refs[0][1]
            log1 = spec.top_log_y_refs[1][1]
            for sample, pts in raw.items():
                converted[sample] = [(x, 10.0 ** (log0 + y / 100.0 * (log1 - log0))) for x, y in pts]
            return converted
    (py0, log0), (py1, log1) = spec.top_log_y_refs
    (px0, x0), (px1, x1) = spec.top_x_refs
    refs = [(0, px0, py0, x0, 0.0), (1, px1, py0, x1, 0.0), (2, px0, py1, x0, 100.0)]
    points = [(str(r["sample"]), float(r["pixel_x"]), float(r["pixel_y"])) for r in pixel_rows]
    try:
        run_datathief_export(spec.image, refs, points, out_csv)
        raw = load_datathief_csv(out_csv)
    except Exception:
        raw = {}
    if finite_datathief(raw, len(pixel_rows)):
        converted = {}
        for sample, pts in raw.items():
            converted[sample] = [(x, 10.0 ** (log0 + y / 100.0 * (log1 - log0))) for x, y in pts]
        return converted
    converted: dict[str, list[tuple[float, float]]] = {}
    for r in pixel_rows:
        sample = str(r["sample"])
        px = float(r["pixel_x"])
        py = float(r["pixel_y"])
        x = x0 + (px - px0) / (px1 - px0) * (x1 - x0)
        logy = log0 + (py - py0) / (py1 - py0) * (log1 - log0)
        converted.setdefault(sample, []).append((x, 10.0**logy))
    return converted


def export_linear_datathief(spec: AxisSpec, pixel_rows: list[dict[str, object]], out_csv: Path) -> dict[str, list[tuple[float, float]]]:
    out_csv.unlink(missing_ok=True)
    if out_csv.exists():
        raw = load_datathief_csv(out_csv)
        if sum(len(v) for v in raw.values()) == len(pixel_rows):
            return raw
    (px0, x0), (px1, x1) = spec.ratio_x_refs
    (py0, y0), (py1, y1) = spec.ratio_y_refs
    refs = [(0, px0, py0, 0.0, 0.0), (1, px1, py0, 100.0, 0.0), (2, px0, py1, 0.0, 100.0)]
    points = [(str(r["sample"]), float(r["pixel_x"]), float(r["pixel_y"])) for r in pixel_rows]
    try:
        run_datathief_export(spec.image, refs, points, out_csv)
        raw = load_datathief_csv(out_csv)
    except Exception:
        raw = {}
    if not finite_datathief(raw, len(pixel_rows)):
        raw = {}
        for r in pixel_rows:
            sample = str(r["sample"])
            px = float(r["pixel_x"])
            py = float(r["pixel_y"])
            x_scaled = (px - px0) / (px1 - px0) * 100.0
            y_scaled = (py - py0) / (py1 - py0) * 100.0
            raw.setdefault(sample, []).append((x_scaled, y_scaled))
    converted: dict[str, list[tuple[float, float]]] = {}
    for sample, pts in raw.items():
        converted[sample] = [
            (
                x0 + x_scaled / 100.0 * (x1 - x0),
                y0 + y_scaled / 100.0 * (y1 - y0),
            )
            for x_scaled, y_scaled in pts
        ]
    return converted


def make_datathief_reference_audit(group: str, source_rows: list[dict[str, str]]) -> tuple[Path, Path]:
    spec = SPECS[group]
    out_dir = OUT / "reference_audit/ian_datathief"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [r for r in source_rows if r["group"] == group and r["sample"] in spec.samples]

    ratio_pixels = find_black_ratio_marker_pixels(spec, rows)
    ratio_export = export_linear_datathief(spec, ratio_pixels, out_dir / f"{group}_ian_ratio_datathief_points.csv")

    sample_indices_ratio = {s: 0 for s in spec.samples}
    audit_rows: list[dict[str, object]] = []
    for r in ratio_pixels:
        sample = str(r["sample"])
        idx = sample_indices_ratio[sample]
        sample_indices_ratio[sample] += 1
        dt_x, dt_y = ratio_export[sample][idx]
        src = float(r["fit_ratio"])
        audit_rows.append(
            {
                "panel": "lower_ratio",
                "group": group,
                "sample": sample,
                "source_x": float(r["bin_center"]),
                "source_y": src,
                "datathief_x": dt_x,
                "datathief_y": dt_y,
                "ratio_datathief_over_source": dt_y / src if src > 0 else math.nan,
                "pixel_x": float(r["pixel_x"]),
                "pixel_y": float(r["pixel_y"]),
                "match_mode": r["match_mode"],
            }
        )

    csv_path = out_dir / f"{group}_ian_datathief_vs_sdcc_source_audit.csv"
    fields = [
        "panel",
        "group",
        "sample",
        "source_x",
        "source_y",
        "datathief_x",
        "datathief_y",
        "ratio_datathief_over_source",
        "pixel_x",
        "pixel_y",
        "match_mode",
    ]
    write_rows(csv_path, audit_rows, fields)

    plt.rcParams.update(style_rc())
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=180,
        sharex=False,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.09},
    )
    ax.set_title(spec.title + ": IAN lower-ratio DataThief audit", fontsize=15.5)
    ax.axhline(1.0, color="0.55", lw=0.9, ls=(0, (4, 4)))
    rax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
    ax.set_ylabel("IAN lower-panel ratio", fontsize=12.5)
    rax.set_ylabel("DataThief / SDCC fit-ratio", fontsize=11.8)
    ax.set_xlim(*spec.xlim)
    rax.set_xlim(*spec.xlim)
    ax.set_ylim(0.84, 1.16)
    rax.set_ylim(0.96, 1.05 if group == "jet" else 1.04)
    for sample in spec.samples:
        pts = [r for r in audit_rows if r["sample"] == sample]
        if not pts:
            continue
        color = COLORS[sample]
        ax.plot(
            [float(r["source_x"]) for r in pts],
            [float(r["source_y"]) for r in pts],
            "o",
            ms=4.0,
            mfc="white",
            mec=color,
            mew=1.0,
            linestyle="none",
            label=f"{sample} SDCC source",
        )
        ax.plot(
            [float(r["datathief_x"]) for r in pts],
            [float(r["datathief_y"]) for r in pts],
            "s",
            ms=3.4,
            color=color,
            linestyle="none",
            label=f"{sample} IAN DataThief" if sample == spec.samples[0] else None,
        )
        rax.plot(
            [float(r["source_x"]) for r in pts],
            [float(r["ratio_datathief_over_source"]) for r in pts],
            "o",
            ms=3.6,
            color=color,
            linestyle="none",
            label=sample,
        )
    rax.set_xlabel(spec.xlabel, fontsize=14)
    ax.legend(frameon=False, ncol=2 if group == "jet" else 1, fontsize=9.2, loc="best")
    png_path = out_dir / f"{group}_ian_datathief_vs_sdcc_source_audit.png"
    fig.savefig(png_path)
    plt.close(fig)

    ratios = [float(r["ratio_datathief_over_source"]) for r in audit_rows if math.isfinite(float(r["ratio_datathief_over_source"]))]
    write_json(
        out_dir / f"{group}_ian_datathief_vs_sdcc_source_manifest.json",
        {
            "status": "ok" if ratios else "no_datathief_points",
            "figure": spec.title,
            "source_image": str(spec.image),
            "source_csv": str(SOURCE_CSV),
            "validated_panel": "lower_fit_ratio",
            "audit_csv": str(csv_path),
            "audit_png": str(png_path),
            "datathief_jar": str(DATATHIEF_JAR),
            "datathief_jar_md5": jar_md5(),
            "n_audit_rows": len(audit_rows),
            "mean_abs_ratio_minus_one": float(np.mean(np.abs(np.asarray(ratios) - 1.0))) if ratios else None,
            "max_abs_ratio_minus_one": float(np.max(np.abs(np.asarray(ratios) - 1.0))) if ratios else None,
            "axis_spec": {
                "top_x_refs": spec.top_x_refs,
                "top_log_y_refs": spec.top_log_y_refs,
                "ratio_x_refs": spec.ratio_x_refs,
                "ratio_y_refs": spec.ratio_y_refs,
            },
        },
    )
    return csv_path, png_path


def make_sdcc_source_reproduction(group: str, source_rows: list[dict[str, str]]) -> tuple[Path, Path]:
    spec = SPECS[group]
    out_dir = OUT / "reference_audit/sdcc_reproduction"
    rows = [
        r
        for r in source_rows
        if r["group"] == group
        and r["sample"] in spec.samples
        and spec.xlim[0] <= float(r["bin_center"]) <= spec.xlim[1]
    ]
    rows.sort(key=lambda r: (spec.samples.index(r["sample"]), float(r["bin_center"])))
    csv_path = out_dir / f"{group}_sdcc_source_reproduction_points.csv"
    fields = [
        "group",
        "sample",
        "bin_low",
        "bin_high",
        "bin_center",
        "value",
        "error",
        "fit_value",
        "fit_ratio",
        "used_in_stitch",
        "source_file",
        "source_object",
    ]
    write_rows(csv_path, rows, fields)

    plt.rcParams.update(style_rc())
    fig = plt.figure(figsize=(7.72, 9.98), dpi=180)
    gs = fig.add_gridspec(2, 1, height_ratios=(3.0, 1.0), hspace=0.04, left=0.15, right=0.98, top=0.95, bottom=0.085)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(*spec.xlim)
    ax.set_ylim(*spec.top_ylim)
    ax.set_ylabel(spec.ylabel, fontsize=15)
    rax.set_xlabel(spec.xlabel, fontsize=15)
    rax.set_ylabel("Source / Fit", fontsize=12.5)
    rax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylim(0.84, 1.16)
    ax.tick_params(which="both", labelbottom=False, labelsize=11.5, length=5)
    rax.tick_params(which="both", labelsize=11.5, length=5)
    for sample in spec.samples:
        pts = [r for r in rows if r["sample"] == sample]
        if not pts:
            continue
        color = COLORS[sample]
        x = [float(r["bin_center"]) for r in pts]
        y = [float(r["value"]) for r in pts]
        yerr = [float(r["error"]) for r in pts]
        ratio = [float(r["fit_ratio"]) for r in pts]
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            ms=3.8,
            color=color,
            ecolor=color,
            elinewidth=0.45,
            linestyle="none",
            label=sample,
        )
        rax.plot(x, ratio, "o", ms=3.4, color=color, linestyle="none")
    fit_pts = sorted(rows, key=lambda r: float(r["bin_center"]))
    if fit_pts:
        ax.plot(
            [float(r["bin_center"]) for r in fit_pts],
            [float(r["fit_value"]) for r in fit_pts],
            color="0.15",
            lw=1.0,
            alpha=0.7,
            label="PPG12 fit",
        )
    ax.legend(frameon=False, ncol=2 if group == "jet" else 1, fontsize=10.5, loc="upper right")
    ax.text(0.03, 0.94, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=14)
    ax.text(0.03, 0.89, "SDCC source reproduction, no DataThief gate", transform=ax.transAxes, ha="left", va="top", fontsize=10.5)
    png_path = out_dir / f"{group}_sdcc_source_reproduction.png"
    fig.savefig(png_path)
    plt.close(fig)
    write_json(
        out_dir / f"{group}_sdcc_source_reproduction_manifest.json",
        {
            "status": "ok",
            "plot_png": str(png_path),
            "points_csv": str(csv_path),
            "group": group,
            "ppg12_reference_csv": str(SOURCE_CSV),
            "source_note": spec.source_note,
            "no_datathief_required": True,
        },
    )
    return png_path, csv_path


def best_constant_scale(rows: list[dict[str, object]], ratio_key: str) -> float:
    ratios = [float(r[ratio_key]) for r in rows if math.isfinite(float(r[ratio_key])) and float(r[ratio_key]) > 0]
    if not ratios:
        return 1.0
    return 1.0 / float(np.median(ratios))


def make_current_vs_sdcc_shape_overlay(
    group: str,
    source_rows: list[dict[str, str]],
    current_rows: list[dict[str, str]],
) -> tuple[Path, Path]:
    spec = SPECS[group]
    out_dir = OUT / ("strict_stitched_photonjet" if group == "photon" else "strict_stitched_inclusivejet")
    src = {
        key(r["sample"], float(r["bin_center"])): r
        for r in source_rows
        if r["group"] == group and r["sample"] in spec.samples and spec.xlim[0] <= float(r["bin_center"]) <= spec.xlim[1]
    }
    cur = {
        key(r["sample"], float(r["bin_center"])): r
        for r in current_rows
        if r["group"] == group and r["sample"] in spec.samples and spec.xlim[0] <= float(r["bin_center"]) <= spec.xlim[1]
    }
    rows: list[dict[str, object]] = []
    for k in sorted(src, key=lambda t: (spec.samples.index(t[0]), t[1])):
        srow = src[k]
        crow = cur.get(k)
        if crow is None:
            continue
        sy = float(srow["value"])
        cy = float(crow["value"])
        if sy <= 0 or cy <= 0:
            continue
        rows.append(
            {
                "group": group,
                "sample": k[0],
                "bin_center": k[1],
                "bin_low": srow["bin_low"],
                "bin_high": srow["bin_high"],
                "ppg12_source_value": sy,
                "ppg12_source_error": float(srow["error"]),
                "ppg12_fit_value": float(srow["fit_value"]),
                "ppg12_source_over_fit": float(srow["fit_ratio"]),
                "current_value_unscaled": cy,
                "current_error_unscaled": float(crow["error"]),
                "current_over_sdcc_unscaled": cy / sy,
            }
        )
    scale = best_constant_scale(rows, "current_over_sdcc_unscaled")
    for row in rows:
        current_scaled = float(row["current_value_unscaled"]) * scale
        current_error_scaled = float(row["current_error_unscaled"]) * scale
        fit = float(row["ppg12_fit_value"])
        source = float(row["ppg12_source_value"])
        row["current_scale_constant"] = scale
        row["current_value_scaled"] = current_scaled
        row["current_error_scaled"] = current_error_scaled
        row["scaled_current_over_sdcc"] = current_scaled / source if source > 0 else math.nan
        row["scaled_current_over_fit"] = current_scaled / fit if fit > 0 else math.nan

    csv_path = out_dir / f"{group}_current_vs_sdcc_shape_overlay_points.csv"
    fields = [
        "group",
        "sample",
        "bin_center",
        "bin_low",
        "bin_high",
        "ppg12_source_value",
        "ppg12_source_error",
        "ppg12_fit_value",
        "ppg12_source_over_fit",
        "current_value_unscaled",
        "current_error_unscaled",
        "current_over_sdcc_unscaled",
        "current_scale_constant",
        "current_value_scaled",
        "current_error_scaled",
        "scaled_current_over_sdcc",
        "scaled_current_over_fit",
    ]
    write_rows(csv_path, rows, fields)

    plt.rcParams.update(style_rc())
    fig = plt.figure(figsize=(7.2, 9.0), dpi=180)
    gs = fig.add_gridspec(2, 1, height_ratios=(3.35, 1.1), hspace=0.04, left=0.14, right=0.985, top=0.965, bottom=0.09)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(*spec.xlim)
    ax.set_ylim(*spec.top_ylim)
    ax.set_ylabel(spec.ylabel, fontsize=16)
    rax.set_xlabel(spec.xlabel, fontsize=15)
    rax.set_ylabel("Current / SDCC", fontsize=12.5)
    rax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylim(0.84, 1.12 if group == "photon" else 1.18)
    ax.tick_params(which="both", labelbottom=False, labelsize=11.5, length=5)
    rax.tick_params(which="both", labelsize=11.5, length=5)

    for sample in spec.samples:
        color = COLORS[sample]
        pts = [r for r in rows if r["sample"] == sample]
        if not pts:
            continue
        x = [float(r["bin_center"]) for r in pts]
        ax.errorbar(
            x,
            [float(r["ppg12_source_value"]) for r in pts],
            yerr=[float(r["ppg12_source_error"]) for r in pts],
            fmt="o",
            ms=3.6,
            mfc="white",
            mec=color,
            mew=1.0,
            ecolor=color,
            elinewidth=0.45,
            linestyle="none",
            zorder=3,
        )
        ax.errorbar(
            x,
            [float(r["current_value_scaled"]) for r in pts],
            yerr=[float(r["current_error_scaled"]) for r in pts],
            fmt="s",
            ms=3.2,
            mfc=color,
            mec=color,
            mew=0.6,
            ecolor=color,
            elinewidth=0.4,
            linestyle="none",
            zorder=4,
        )
        rax.errorbar(
            x,
            [float(r["scaled_current_over_sdcc"]) for r in pts],
            yerr=[
                float(r["current_error_scaled"]) / float(r["ppg12_source_value"])
                if float(r["ppg12_source_value"]) > 0
                else math.nan
                for r in pts
            ],
            fmt="s",
            ms=3.0,
            mfc=color,
            mec=color,
            mew=0.5,
            ecolor=color,
            elinewidth=0.35,
            linestyle="none",
            zorder=4,
        )

    fit_pts = sorted(rows, key=lambda r: float(r["bin_center"]))
    if fit_pts:
        ax.plot(
            [float(r["bin_center"]) for r in fit_pts],
            [float(r["ppg12_fit_value"]) for r in fit_pts],
            color="red",
            lw=1.1,
            label="PPG12 fit",
        )

    ax.text(0.98, 0.94, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=14.5)
    ax.text(0.98, 0.86, f"{group.replace('jet','inclusive+jet').replace('photon','photon+jet')} stitch comparison", transform=ax.transAxes, ha="right", va="top", fontsize=12.5)
    ax.text(0.98, 0.80, f"current output scaled by one constant = {scale:.6g}", transform=ax.transAxes, ha="right", va="top", fontsize=10.8)
    handles = [
        Line2D([0], [0], marker="o", lw=0, color="black", mfc="white", mec="black", label="PPG12 SDCC source"),
        Line2D([0], [0], marker="s", lw=0, color="black", mfc="black", mec="black", label="Current pp output"),
        Line2D([0], [0], lw=1.1, color="red", label="PPG12 fit"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower left", fontsize=10.5)
    png_path = out_dir / f"{group}_current_vs_sdcc_shape_overlay.png"
    fig.savefig(png_path)
    plt.close(fig)

    ratios = [float(r["scaled_current_over_sdcc"]) for r in rows if math.isfinite(float(r["scaled_current_over_sdcc"]))]
    unscaled_ratios = [float(r["current_over_sdcc_unscaled"]) for r in rows if math.isfinite(float(r["current_over_sdcc_unscaled"]))]
    write_json(
        out_dir / f"{group}_current_vs_sdcc_shape_overlay_manifest.json",
        {
            "status": "ok_shape_overlay",
            "plot_png": str(png_path),
            "points_csv": str(csv_path),
            "current_scale_constant": scale,
            "scale_definition": "1 / median(current_unscaled / PPG12_SDCC_source) over matched plotted bins",
            "ratio_panel": "scaled_current / PPG12_SDCC_source",
            "unscaled_current_over_sdcc_summary": {
                "median": float(np.median(unscaled_ratios)) if unscaled_ratios else None,
                "min": float(np.min(unscaled_ratios)) if unscaled_ratios else None,
                "max": float(np.max(unscaled_ratios)) if unscaled_ratios else None,
            },
            "scaled_current_over_sdcc_summary": {
                "median": float(np.median(ratios)) if ratios else None,
                "mean": float(np.mean(ratios)) if ratios else None,
                "min": float(np.min(ratios)) if ratios else None,
                "max": float(np.max(ratios)) if ratios else None,
                "rms_to_one": float(math.sqrt(np.mean((np.asarray(ratios) - 1.0) ** 2))) if ratios else None,
            },
            "note": "This is the eyeball/shape plot Justin requested; the strict no-scale plot remains separately written.",
        },
    )
    return png_path, csv_path


def make_ian_crop_vs_sdcc_side_by_side(group: str, sdcc_png: Path) -> Path:
    spec = SPECS[group]
    out_dir = OUT / "reference_audit/sdcc_reproduction"
    crop = np.asarray(Image.open(spec.image).convert("RGB"))
    repro = np.asarray(Image.open(sdcc_png).convert("RGB"))
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 8.0), dpi=180, gridspec_kw={"width_ratios": [1.0, 1.15]})
    for ax, img, title in [
        (axes[0], crop, "Exact PPG12 IAN crop"),
        (axes[1], repro, "Live SDCC source reproduction"),
    ]:
        ax.imshow(img)
        ax.set_title(title, fontsize=15)
        ax.axis("off")
    fig.suptitle(spec.title + " visual source check", fontsize=17)
    fig.tight_layout()
    png_path = out_dir / f"{group}_ian_crop_vs_sdcc_reproduction_side_by_side.png"
    fig.savefig(png_path)
    plt.close(fig)
    write_json(
        out_dir / f"{group}_ian_crop_vs_sdcc_reproduction_side_by_side_manifest.json",
        {
            "status": "ok",
            "plot_png": str(png_path),
            "ian_crop": str(spec.image),
            "sdcc_reproduction_png": str(sdcc_png),
            "purpose": "visual eyeball check replacing DataThief as the primary source-validation gate",
        },
    )
    return png_path


def style_rc() -> dict[str, object]:
    return {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }


def key(sample: str, x: float) -> tuple[str, float]:
    return sample, round(float(x), 6)


def make_strict_stitch_plot(group: str, source_rows: list[dict[str, str]], current_rows: list[dict[str, str]]) -> tuple[Path, Path]:
    spec = SPECS[group]
    out_dir = OUT / ("strict_stitched_photonjet" if group == "photon" else "strict_stitched_inclusivejet")
    src = {
        key(r["sample"], float(r["bin_center"])): r
        for r in source_rows
        if r["group"] == group and r["sample"] in spec.samples and spec.xlim[0] <= float(r["bin_center"]) <= spec.xlim[1]
    }
    cur = {
        key(r["sample"], float(r["bin_center"])): r
        for r in current_rows
        if r["group"] == group and r["sample"] in spec.samples and spec.xlim[0] <= float(r["bin_center"]) <= spec.xlim[1]
    }
    rows: list[dict[str, object]] = []
    missing = []
    for k in sorted(src, key=lambda t: (spec.samples.index(t[0]), t[1])):
        srow = src[k]
        crow = cur.get(k)
        if crow is None:
            missing.append({"sample": k[0], "x": k[1]})
            continue
        sy = float(srow["value"])
        cy = float(crow["value"])
        sey = float(srow["error"])
        cey = float(crow["error"])
        rows.append(
            {
                "group": group,
                "sample": k[0],
                "bin_center": k[1],
                "bin_low": srow["bin_low"],
                "bin_high": srow["bin_high"],
                "ppg12_source_value": sy,
                "ppg12_source_error": sey,
                "july1_value": cy,
                "july1_error": cey,
                "ratio_july1_over_ppg12_source": cy / sy if sy > 0 else math.nan,
                "ratio_error": cey / sy if sy > 0 else math.nan,
                "july1_value_mode": crow["value_mode"],
                "july1_source_files": crow["source_files"],
                "july1_events_processed_metadata": crow["events_processed_metadata"],
            }
        )

    csv_path = out_dir / f"{group}_strict_july1_over_ppg12_source_points.csv"
    fields = [
        "group",
        "sample",
        "bin_center",
        "bin_low",
        "bin_high",
        "ppg12_source_value",
        "ppg12_source_error",
        "july1_value",
        "july1_error",
        "ratio_july1_over_ppg12_source",
        "ratio_error",
        "july1_value_mode",
        "july1_source_files",
        "july1_events_processed_metadata",
    ]
    write_rows(csv_path, rows, fields)

    plt.rcParams.update(style_rc())
    fig = plt.figure(figsize=(7.72, 9.98), dpi=180)
    gs = fig.add_gridspec(2, 1, height_ratios=(3.25, 1.0), hspace=0.04, left=0.16, right=0.98, top=0.97, bottom=0.085)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(*spec.xlim)
    ax.set_ylim(*spec.top_ylim)
    ax.set_ylabel(spec.ylabel, fontsize=16)
    rax.set_xlabel(spec.xlabel, fontsize=16)
    rax.set_ylabel("July1 / PPG12", fontsize=13)
    rax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylim(0.75, 1.25)
    if group == "jet":
        rax.set_ylim(0.0, 2.6)
    ax.tick_params(which="both", labelsize=12.5, labelbottom=False, length=6)
    rax.tick_params(which="both", labelsize=12.5, length=6)

    for sample in spec.samples:
        color = COLORS[sample]
        spts = [r for r in src.values() if r["sample"] == sample]
        cpts = [r for r in rows if r["sample"] == sample]
        if spts:
            ax.errorbar(
                [float(r["bin_center"]) for r in spts],
                [float(r["value"]) for r in spts],
                yerr=[float(r["error"]) for r in spts],
                fmt="o",
                ms=4.4,
                mfc="white",
                mec=color,
                mew=1.0,
                ecolor=color,
                elinewidth=0.55,
                linestyle="none",
                zorder=3,
            )
        if cpts:
            ax.errorbar(
                [float(r["bin_center"]) for r in cpts],
                [float(r["july1_value"]) for r in cpts],
                yerr=[float(r["july1_error"]) for r in cpts],
                fmt="s",
                ms=3.9,
                mfc=color,
                mec=color,
                mew=0.8,
                ecolor=color,
                elinewidth=0.55,
                linestyle="none",
                zorder=4,
            )
            rax.errorbar(
                [float(r["bin_center"]) for r in cpts],
                [float(r["ratio_july1_over_ppg12_source"]) for r in cpts],
                yerr=[float(r["ratio_error"]) for r in cpts],
                fmt="s",
                ms=3.4,
                mfc=color,
                mec=color,
                mew=0.6,
                ecolor=color,
                elinewidth=0.5,
                linestyle="none",
                zorder=4,
            )

    handles = [
        Line2D([0], [0], marker="o", color="black", lw=0, mfc="white", mec="black", label="PPG12 SDCC source"),
        Line2D([0], [0], marker="s", color="black", lw=0, mfc="black", mec="black", label="July1 RecoilJets"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower left", fontsize=11.5)
    ax.text(0.98, 0.96, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=15)
    ax.text(0.98, 0.90, spec.source_note, transform=ax.transAxes, ha="right", va="top", fontsize=9.8)
    png_path = out_dir / f"{group}_strict_july1_over_ppg12_source.png"
    fig.savefig(png_path)
    plt.close(fig)

    ratios = [float(r["ratio_july1_over_ppg12_source"]) for r in rows if math.isfinite(float(r["ratio_july1_over_ppg12_source"]))]
    status = "ok" if not missing else "missing_bins"
    strict_status_note = "direct bin-by-bin overlay"
    if group == "jet" and ratios:
        outside = [r for r in ratios if r < 0.8 or r > 1.2]
        if outside:
            status = "weighted_count_exposure_or_component_mismatch"
            strict_status_note = (
                "PPG12 Fig.6 is a weighted-count spectrum, not a normalized cross-section; "
                "the July1 all-component current output is shown without any fitted scale and "
                "therefore exposes the remaining exposure/component mismatch."
            )
    write_json(
        out_dir / f"{group}_strict_july1_over_ppg12_source_manifest.json",
        {
            "status": status,
            "strict_status_note": strict_status_note,
            "plot_png": str(png_path),
            "points_csv": str(csv_path),
            "group": group,
            "campaign_tag": TAG,
            "ppg12_reference_csv": str(SOURCE_CSV),
            "july1_current_csv": str(CURRENT_CSV),
            "missing_bins": missing,
            "n_reference_bins": len(src),
            "n_matched_bins": len(rows),
            "no_global_fitted_scale_applied": True,
            "normalization_note": "Photon uses raw*xsec/events no-bin-width PPG12 Fig.5 contract; jet uses PPG12 legacy xsec/jet50 weighted-count convention with no fitted scale.",
            "ratio_summary": {
                "mean": float(np.mean(ratios)) if ratios else None,
                "median": float(np.median(ratios)) if ratios else None,
                "rms_to_one": float(math.sqrt(np.mean((np.asarray(ratios) - 1.0) ** 2))) if ratios else None,
                "min": float(np.min(ratios)) if ratios else None,
                "max": float(np.max(ratios)) if ratios else None,
            },
        },
    )
    return png_path, csv_path


def sum_component_rows(rows: list[dict[str, str]], modes: dict[str, set[tuple[str, str]]]) -> dict[str, dict[tuple[str, float], dict[str, float]]]:
    out: dict[str, dict[tuple[str, float], dict[str, float]]] = {mode: {} for mode in modes}
    for row in rows:
        if row["group"] != "jet":
            continue
        sample = row["sample"]
        x = round(float(row["bin_center"]), 6)
        component_key = (row["period"], row["component"])
        for mode, accepted in modes.items():
            if component_key not in accepted:
                continue
            dest = out[mode].setdefault(
                (sample, x),
                {
                    "value": 0.0,
                    "err2": 0.0,
                    "source_files": 0.0,
                    "events_processed_metadata": 0.0,
                },
            )
            dest["value"] += float(row["value"])
            dest["err2"] += float(row["error"]) ** 2
            dest["source_files"] += float(row["source_files"])
            dest["events_processed_metadata"] += float(row["events_processed_metadata"])
    return out


def make_jet_component_audit(source_rows: list[dict[str, str]], component_rows: list[dict[str, str]]) -> tuple[Path, Path]:
    spec = SPECS["jet"]
    out_dir = OUT / "strict_stitched_inclusivejet"
    modes = {
        "all_components": {("0mrad", "single"), ("0mrad", "double"), ("1p5mrad", "single"), ("1p5mrad", "double")},
        "single_all_periods": {("0mrad", "single"), ("1p5mrad", "single")},
        "double_all_periods": {("0mrad", "double"), ("1p5mrad", "double")},
        "0mrad_single": {("0mrad", "single")},
        "1p5mrad_single": {("1p5mrad", "single")},
        "0mrad_double": {("0mrad", "double")},
        "1p5mrad_double": {("1p5mrad", "double")},
    }
    src = {
        key(r["sample"], float(r["bin_center"])): r
        for r in source_rows
        if r["group"] == "jet" and r["sample"] in spec.samples and spec.xlim[0] <= float(r["bin_center"]) <= spec.xlim[1]
    }
    summed = sum_component_rows(component_rows, modes)
    rows: list[dict[str, object]] = []
    for mode in modes:
        for k in sorted(src, key=lambda t: (spec.samples.index(t[0]), t[1])):
            srow = src[k]
            crow = summed[mode].get(k)
            if crow is None:
                continue
            sy = float(srow["value"])
            rows.append(
                {
                    "group": "jet",
                    "mode": mode,
                    "sample": k[0],
                    "bin_center": k[1],
                    "bin_low": srow["bin_low"],
                    "bin_high": srow["bin_high"],
                    "ppg12_source_value": sy,
                    "ppg12_source_error": float(srow["error"]),
                    "july1_value": crow["value"],
                    "july1_error": math.sqrt(crow["err2"]),
                    "ratio_july1_over_ppg12_source": crow["value"] / sy if sy > 0 else math.nan,
                    "ratio_error": math.sqrt(crow["err2"]) / sy if sy > 0 else math.nan,
                    "july1_source_files": int(crow["source_files"]),
                    "july1_events_processed_metadata": crow["events_processed_metadata"],
                }
            )

    csv_path = out_dir / "jet_component_exposure_audit_points.csv"
    fields = [
        "group",
        "mode",
        "sample",
        "bin_center",
        "bin_low",
        "bin_high",
        "ppg12_source_value",
        "ppg12_source_error",
        "july1_value",
        "july1_error",
        "ratio_july1_over_ppg12_source",
        "ratio_error",
        "july1_source_files",
        "july1_events_processed_metadata",
    ]
    write_rows(csv_path, rows, fields)

    plt.rcParams.update(style_rc())
    fig, axes = plt.subplots(3, 1, figsize=(7.72, 9.98), dpi=180, sharex=True, gridspec_kw={"hspace": 0.08})
    panels = [
        ("Raw all-component current / PPG12 Fig.6", ("all_components",), (0.0, 2.6)),
        ("Single-interaction components", ("single_all_periods", "0mrad_single", "1p5mrad_single"), (0.0, 2.1)),
        ("Double-interaction components", ("double_all_periods", "0mrad_double", "1p5mrad_double"), (0.0, 0.65)),
    ]
    mode_styles = {
        "all_components": ("black", "s"),
        "single_all_periods": ("#1f77b4", "o"),
        "double_all_periods": ("#9467bd", "o"),
        "0mrad_single": ("#2ca02c", "^"),
        "1p5mrad_single": ("#ff7f0e", "v"),
        "0mrad_double": ("#17becf", "^"),
        "1p5mrad_double": ("#d62728", "v"),
    }
    for ax, (title, panel_modes, ylim) in zip(axes, panels):
        ax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
        ax.set_title(title, fontsize=12.5, loc="left")
        ax.set_xlim(*spec.xlim)
        ax.set_ylim(*ylim)
        ax.set_ylabel("July1 / PPG12", fontsize=11.5)
        for mode in panel_modes:
            pts = [r for r in rows if r["mode"] == mode and math.isfinite(float(r["ratio_july1_over_ppg12_source"]))]
            if not pts:
                continue
            color, marker = mode_styles[mode]
            ax.errorbar(
                [float(r["bin_center"]) for r in pts],
                [float(r["ratio_july1_over_ppg12_source"]) for r in pts],
                yerr=[float(r["ratio_error"]) for r in pts],
                fmt=marker,
                ms=2.7,
                color=color,
                ecolor=color,
                elinewidth=0.35,
                linestyle="none",
                alpha=0.72,
                label=mode,
            )
        ax.legend(frameon=False, fontsize=9.3, ncol=3, loc="upper right")
        ax.tick_params(which="both", labelsize=10.5, length=5)
    axes[-1].set_xlabel(spec.xlabel, fontsize=13.5)
    fig.text(
        0.02,
        0.015,
        "No fitted scale applied. This audit explains why Fig.6 raw weighted-count parity is not yet bin-by-bin exact for the July1 component mix/exposure.",
        fontsize=9.5,
    )
    png_path = out_dir / "jet_component_exposure_audit.png"
    fig.savefig(png_path)
    plt.close(fig)

    summary = {}
    for mode in modes:
        mode_ratios = [
            float(r["ratio_july1_over_ppg12_source"])
            for r in rows
            if r["mode"] == mode and math.isfinite(float(r["ratio_july1_over_ppg12_source"]))
        ]
        summary[mode] = {
            "n": len(mode_ratios),
            "median": float(np.median(mode_ratios)) if mode_ratios else None,
            "min": float(np.min(mode_ratios)) if mode_ratios else None,
            "max": float(np.max(mode_ratios)) if mode_ratios else None,
        }
    write_json(
        out_dir / "jet_component_exposure_audit_manifest.json",
        {
            "status": "diagnostic_exposure_component_mismatch_not_a_strict_pass",
            "plot_png": str(png_path),
            "points_csv": str(csv_path),
            "ppg12_reference_csv": str(SOURCE_CSV),
            "july1_component_csv": str(JET_COMPONENT_CSV),
            "no_global_fitted_scale_applied": True,
            "ratio_summary_by_mode": summary,
            "interpretation": (
                "PPG12 Fig.6 source is a weighted-count diagnostic from MC_efficiency_jet*_bdt_nom.root. "
                "The July1 current campaign has period/component-split outputs; the raw all-component sum "
                "does not bin-by-bin equal the IAN weighted-count exposure. This artifact exposes the component "
                "structure before any decision to reproduce the PPG12 no-suffix exposure contract exactly."
            ),
        },
    )
    return png_path, csv_path


def make_iso_artifacts() -> tuple[Path, Path]:
    rows = read_rows(TEFF_CSV)
    out_dir = OUT / "strict_iso_efficiency"
    curves = {r["object"]: [] for r in rows}
    for r in rows:
        curves[r["object"]].append(r)
    by_obj = {obj: sorted(vals, key=lambda r: int(r["bin"])) for obj, vals in curves.items()}

    # Exact absolute curves available from PPG12 TEfficiency.
    reco = by_obj.get("eff_reco_eta_0", [])
    all_eff = by_obj.get("eff_all_eta_0", [])
    iso_cond = by_obj.get("eff_iso_eta_0", [])
    id_cond = by_obj.get("eff_id_eta_0", [])

    derived_rows: list[dict[str, object]] = []
    for r in reco:
        derived_rows.append({**r, "derived_curve": "reco_absolute", "derived_eff": r["eff"], "strict_status": "stored_ppg12_tefficiency"})
    for r in all_eff:
        derived_rows.append({**r, "derived_curve": "reco_x_iso_x_id_absolute", "derived_eff": r["eff"], "strict_status": "stored_ppg12_tefficiency"})
    for r_iso, r_id in zip(iso_cond, id_cond):
        # This product is a useful stored-conditional diagnostic, not the same
        # as an independently stored absolute reco x eID curve.
        derived = float(r_iso["eff"]) * float(r_id["eff"])
        derived_rows.append(
            {
                **r_iso,
                "object": "eff_iso_eta_0_times_eff_id_eta_0",
                "curve": "iso_cond_x_id_cond",
                "derived_curve": "conditional_iso_times_conditional_id",
                "derived_eff": derived,
                "strict_status": "diagnostic_not_absolute_reco_x_eid",
            }
        )

    csv_path = out_dir / "strict_isolation_efficiency_reference_points.csv"
    fields = list(derived_rows[0].keys()) if derived_rows else []
    write_rows(csv_path, derived_rows, fields)

    plt.rcParams.update(style_rc())
    fig, ax = plt.subplots(figsize=(7.72, 5.3), dpi=180)
    for label, obj, color, marker in [
        ("PPG12 reco", "eff_reco_eta_0", "#1b4f9c", "o"),
        ("PPG12 all = reco+iso+id", "eff_all_eta_0", "#c43c39", "s"),
        ("conditional iso x conditional id", "conditional", "#4a8f3a", "^"),
    ]:
        if obj == "conditional":
            pts = [r for r in derived_rows if r["derived_curve"] == "conditional_iso_times_conditional_id"]
            x = [float(r["bin_center"]) for r in pts]
            y = [float(r["derived_eff"]) for r in pts]
            yerr = None
        else:
            pts = by_obj.get(obj, [])
            x = [float(r["bin_center"]) for r in pts]
            y = [float(r["eff"]) for r in pts]
            yerr = [[float(r["err_low"]) for r in pts], [float(r["err_high"]) for r in pts]]
        if x:
            ax.errorbar(x, y, yerr=yerr, fmt=marker, ms=4.0, linestyle="none", color=color, label=label)
    ax.set_xlim(8, 45)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(r"Truth $E_T^\gamma$ [GeV]", fontsize=14)
    ax.set_ylabel("Efficiency", fontsize=14)
    ax.text(0.03, 0.07, "July1 current strict TEfficiency-equivalent objects: missing in raw outputs", transform=ax.transAxes, fontsize=10.5)
    ax.legend(frameon=False, fontsize=10.8, loc="lower right")
    png_path = out_dir / "strict_isolation_efficiency_reference_missing_current.png"
    fig.savefig(png_path)
    plt.close(fig)

    write_json(
        out_dir / "strict_isolation_efficiency_manifest.json",
        {
            "status": "missing_current_tefficiency_equivalent",
            "reference_csv": str(TEFF_CSV),
            "points_csv": str(csv_path),
            "plot_png": str(png_path),
            "ppg12_reference_objects": ["eff_reco_eta_0", "eff_id_eta_0", "eff_iso_eta_0", "eff_all_eta_0"],
            "current_required_objects": [
                "eff_reco_eta_0 equivalent numerator/denominator",
                "eff_id_eta_0 equivalent numerator/denominator",
                "eff_iso_eta_0 equivalent numerator/denominator",
                "eff_all_eta_0 equivalent numerator/denominator",
            ],
            "caveat": "PPG12 stores reco and all as absolute TEfficiencies, while iso/id are conditional; an exact absolute reco x eID curve is not independently recoverable from these four objects alone.",
            "required_code_fix_next_sim_pass": "Add PPG12-style truth-photon efficiency numerator/denominator histograms in RecoilJets for reco, iso conditional, id conditional, and all absolute, with the same eta=0 pT binning and weighted-event semantics.",
        },
    )
    return png_path, csv_path


def write_provenance(outputs: dict[str, str]) -> Path:
    path = OUT / "strict_ppg12_sim_parity_provenance.md"
    lines = [
        "# THE-76 Strict July 1 PPG12 SIM Parity Provenance",
        "",
        f"Campaign tag: `{TAG}`",
        "",
        "PPG12 reference:",
        f"- IAN PDF crops: `{FIG5_CROP}` and `{FIG6_CROP}`",
        f"- Live SDCC source CSV: `{SOURCE_CSV}`",
        "- Fig. 5 source: `/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root`",
        "- Fig. 6 source: `/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet*_bdt_nom.root`",
        "",
        "July 1 current source:",
        f"- Raw SIM aggregate CSV: `{CURRENT_CSV}`",
        f"- Inclusive-jet period/component aggregate CSV: `{JET_COMPONENT_CSV}`",
        "",
        "Generated artifacts:",
    ]
    for key_name, value in sorted(outputs.items()):
        lines.append(f"- {key_name}: `{value}`")
    lines.extend(
        [
            "",
            "Caveats:",
            "- DataThief is not used as the primary gate here; the primary source check is exact IAN crop versus live SDCC reproduction side by side.",
            "- Stitched spectra apply no arbitrary fitted global normalization.",
            "- Photon normalization uses the PPG12 Fig. 5 no-bin-width contract: raw*xsec/events.",
            "- Inclusive-jet normalization uses the PPG12 Fig. 6 weighted-count contract: raw*legacy_xsec/jet50.",
            "- The inclusive-jet Fig. 6 raw all-component overlay is intentionally not called a strict pass when the weighted-count exposure/component mix differs from the PPG12 no-suffix source.",
            "- Strict isolation efficiency cannot be claimed from the July 1 raw outputs unless TEfficiency-equivalent numerator/denominator objects are present.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> None:
    source_rows = read_rows(SOURCE_CSV)
    current_rows = read_rows(CURRENT_CSV)
    jet_component_rows = read_rows(JET_COMPONENT_CSV)
    outputs: dict[str, str] = {}
    for group in ("photon", "jet"):
        sdcc_repro_png, sdcc_repro_csv = make_sdcc_source_reproduction(group, source_rows)
        side_by_side_png = make_ian_crop_vs_sdcc_side_by_side(group, sdcc_repro_png)
        strict_png, strict_csv = make_strict_stitch_plot(group, source_rows, current_rows)
        shape_png, shape_csv = make_current_vs_sdcc_shape_overlay(group, source_rows, current_rows)
        outputs[f"{group}_sdcc_reproduction_png"] = str(sdcc_repro_png)
        outputs[f"{group}_sdcc_reproduction_csv"] = str(sdcc_repro_csv)
        outputs[f"{group}_ian_vs_sdcc_side_by_side_png"] = str(side_by_side_png)
        outputs[f"{group}_strict_png"] = str(strict_png)
        outputs[f"{group}_strict_csv"] = str(strict_csv)
        outputs[f"{group}_shape_overlay_png"] = str(shape_png)
        outputs[f"{group}_shape_overlay_csv"] = str(shape_csv)
    jet_component_png, jet_component_csv = make_jet_component_audit(source_rows, jet_component_rows)
    outputs["jet_component_audit_png"] = str(jet_component_png)
    outputs["jet_component_audit_csv"] = str(jet_component_csv)
    iso_png, iso_csv = make_iso_artifacts()
    outputs["iso_png"] = str(iso_png)
    outputs["iso_csv"] = str(iso_csv)
    provenance = write_provenance(outputs)
    outputs["provenance"] = str(provenance)
    print(json.dumps(outputs, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
