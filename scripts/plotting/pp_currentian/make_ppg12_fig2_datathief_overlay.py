#!/usr/bin/env python3
"""Digitize PPG12 IAN Fig. 2 and regenerate the two isolation-fraction pads.

The source image is the rendered page containing Fig. 2 in the current PPG12
IAN.  Colored curve pixels are detected near the integer isolation-threshold
bin centers, then exported through the local DataThief jar coordinate transform
using explicit axis reference points.  The direct-photon pad has essentially
degenerate colored curves in the rendered IAN image, so the visible common
curve is used for all three direct-photon pT intervals and this caveat is
recorded in the CSV/manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_SOURCE = REPO / "tmp/ppg12_ian_fig2/page-006.png"
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_fig2_datathief_overlay_20260630"
)
HELPER_DIR = REPO / "scripts/plotting/pp_currentian"
DATATHIEF_JAR = REPO / "local_tools/datathief/Datathief.jar"

# Reuse the existing local DataThief bridge rather than duplicating the Java
# reflection shim in this focused Fig. 2 script.
sys.path.insert(0, str(HELPER_DIR))
from make_ppg12_datathief_validation_overlays import run_datathief_export  # noqa: E402


@dataclass(frozen=True)
class AxisSpec:
    name: str
    x0: float
    x1: float
    y_top: float
    y_bottom: float
    y_min: float
    y_max: float
    source_crop: tuple[int, int, int, int]


@dataclass(frozen=True)
class SeriesSpec:
    key: str
    label: str
    color: str
    mask_name: str


AXES = {
    "direct": AxisSpec(
        name="direct",
        x0=240.0,
        x1=637.0,
        y_top=728.0,
        y_bottom=999.0,
        y_min=0.90,
        y_max=1.05,
        source_crop=(185, 640, 690, 1130),
    ),
    "fragmentation": AxisSpec(
        name="fragmentation",
        x0=748.0,
        x1=1143.0,
        y_top=728.0,
        y_bottom=999.0,
        y_min=0.80,
        y_max=1.10,
        source_crop=(690, 640, 1200, 1130),
    ),
}

SERIES = [
    SeriesSpec("pt10_15", r"$10 < p_T^\gamma < 15$ GeV", "#cc33cc", "magenta"),
    SeriesSpec("pt15_20", r"$15 < p_T^\gamma < 20$ GeV", "#2ca02c", "green"),
    SeriesSpec("pt25_30", r"$25 < p_T^\gamma < 30$ GeV", "#1f77b4", "blue"),
]

X_VALUES = np.asarray([i + 0.5 for i in range(20)], dtype=float)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return ap.parse_args()


def jar_md5() -> str | None:
    if not DATATHIEF_JAR.exists():
        return None
    return hashlib.md5(DATATHIEF_JAR.read_bytes()).hexdigest()


def mask_for_color(rgb: np.ndarray, name: str) -> np.ndarray:
    r = rgb[..., 0].astype(np.int16)
    g = rgb[..., 1].astype(np.int16)
    b = rgb[..., 2].astype(np.int16)
    if name == "magenta":
        saturated = (r > 145) & (g < 95) & (b > 105)
        antialias = (r > 120) & (b > 100) & (g < 210) & ((r - g) > 20) & ((b - g) > 20)
        return saturated | antialias
    if name == "green":
        return (g > 70) & (r < 135) & (b < 130)
    if name == "blue":
        return (b > 120) & (r < 135) & (g < 175)
    raise ValueError(f"Unknown color mask: {name}")


def data_to_pixel_x(x_value: float, axis: AxisSpec) -> float:
    return axis.x0 + (x_value / 20.0) * (axis.x1 - axis.x0)


def approx_pixel_to_data_y(pixel_y: float, axis: AxisSpec) -> float:
    return axis.y_min + (axis.y_bottom - pixel_y) / (axis.y_bottom - axis.y_top) * (
        axis.y_max - axis.y_min
    )


def y_clusters(ys: np.ndarray) -> list[np.ndarray]:
    if len(ys) == 0:
        return []
    ys = np.sort(ys.astype(int))
    groups: list[list[int]] = [[int(ys[0])]]
    for val in ys[1:]:
        val = int(val)
        if val <= groups[-1][-1] + 3:
            groups[-1].append(val)
        else:
            groups.append([val])
    return [np.asarray(g, dtype=int) for g in groups if len(g) >= 3]


def extract_pixels_for_series(
    rgb: np.ndarray,
    axis: AxisSpec,
    series: SeriesSpec,
    *,
    max_data_y: float,
    min_pixels: int = 8,
) -> list[tuple[str, float, float]]:
    mask = mask_for_color(rgb, series.mask_name)
    x_min = int(round(axis.x0))
    x_max = int(round(axis.x1))
    y_min = int(round(axis.y_top))
    y_max = int(round(axis.y_bottom))
    mask[:y_min, :] = False
    mask[y_max + 1 :, :] = False
    mask[:, :x_min] = False
    mask[:, x_max + 1 :] = False

    pixels: list[tuple[str, float, float]] = []
    for x_value in X_VALUES:
        px = data_to_pixel_x(float(x_value), axis)
        found = None
        for half_width in (4, 6, 9, 12):
            xlo = max(x_min, int(round(px)) - half_width)
            xhi = min(x_max, int(round(px)) + half_width)
            yy, xx = np.where(mask[y_min : y_max + 1, xlo : xhi + 1])
            if len(yy) < min_pixels:
                continue
            yy = yy + y_min
            xx = xx + xlo
            candidates = []
            for group in y_clusters(yy):
                select = np.isin(yy, group)
                if int(select.sum()) < min_pixels:
                    continue
                py = float(np.mean(yy[select]))
                py_data = approx_pixel_to_data_y(py, axis)
                if axis.y_min - 0.04 <= py_data <= max_data_y:
                    candidates.append(
                        {
                            "n": int(select.sum()),
                            "px": float(np.mean(xx[select])),
                            "py": py,
                            "data_y": py_data,
                        }
                    )
            if candidates:
                # Legend strokes sit above the data curve in both Fig. 2 pads.
                # Choose the lowest valid data value when multiple colored
                # clusters are present in the same x window.
                found = sorted(candidates, key=lambda item: item["data_y"])[0]
                break
        if found is not None:
            pixels.append((series.key, found["px"], found["py"]))
    return pixels


def load_datathief_csv(path: Path) -> dict[str, list[tuple[float, float]]]:
    out: dict[str, list[tuple[float, float]]] = {}
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            out.setdefault(row["series"], []).append(
                (float(row["datathief_x"]), float(row["datathief_y"]))
            )
    return out


def export_with_retry(
    source: Path,
    refs: list[tuple[int, float, float, float, float]],
    pixel_points: list[tuple[str, float, float]],
    raw_csv: Path,
    *,
    expected_rows: int,
    max_attempts: int = 6,
) -> dict[str, list[tuple[float, float]]]:
    last_error = "not run"
    for attempt in range(1, max_attempts + 1):
        if raw_csv.exists():
            raw_csv.unlink()
        try:
            run_datathief_export(source, refs, pixel_points, raw_csv)
            raw = load_datathief_csv(raw_csv)
            flat = [point for rows in raw.values() for point in rows]
            valid = (
                len(flat) == expected_rows
                and all(np.isfinite(x) and np.isfinite(y) for x, y in flat)
                and (len(flat) == 0 or np.nanmax([abs(y) for _, y in flat]) > 1.0)
            )
            if valid:
                return raw
            last_error = f"invalid DataThief export on attempt {attempt}: {raw_csv}"
        except Exception as exc:  # DataThief occasionally starts before image load is stable.
            last_error = f"{type(exc).__name__}: {exc}"
        time.sleep(0.35)
    raise RuntimeError(last_error)


def export_pad_points(
    source: Path,
    rgb: np.ndarray,
    outdir: Path,
    axis: AxisSpec,
) -> tuple[dict[str, list[dict[str, float | str]]], dict[str, object]]:
    outdir.mkdir(parents=True, exist_ok=True)
    source_crop = outdir / f"ppg12_fig2_source_{axis.name}_crop.png"
    Image.open(source).convert("RGB").crop(axis.source_crop).save(source_crop)

    if axis.name == "direct":
        blue = SERIES[-1]
        common_pixels = extract_pixels_for_series(
            rgb,
            axis,
            blue,
            max_data_y=1.006,
            min_pixels=6,
        )
        pixel_points: list[tuple[str, float, float]] = []
        for series in SERIES:
            pixel_points.extend((series.key, px, py) for _, px, py in common_pixels)
        caveat = (
            "direct pad: colored pT curves are visually degenerate/overpainted "
            "in the IAN PNG; visible common curve digitized once and assigned "
            "to all three pT intervals"
        )
    else:
        pixel_points = []
        for series in SERIES:
            min_pixels = 2 if series.mask_name == "magenta" else 6
            pixel_points.extend(
                extract_pixels_for_series(
                    rgb,
                    axis,
                    series,
                    max_data_y=1.025,
                    min_pixels=min_pixels,
                )
            )
        caveat = "fragmentation pad: three colored pT curves digitized separately"

    raw_csv = outdir / f"ppg12_fig2_{axis.name}_datathief_raw_shifted.csv"
    # Keep the shifted DataThief y range small.  The Java exporter handles the
    # 15/30 unit shifted ranges reliably; a 300 unit range collapsed the y
    # transform to zero for the fragmentation pad in this figure.
    y_shift_scale = 100.0
    refs = [
        (0, axis.x0, axis.y_bottom, 0.0, 0.0),
        (1, axis.x1, axis.y_bottom, 20.0, 0.0),
        (2, axis.x0, axis.y_top, 0.0, (axis.y_max - axis.y_min) * y_shift_scale),
    ]
    raw = export_with_retry(
        source,
        refs,
        pixel_points,
        raw_csv,
        expected_rows=len(pixel_points),
    )

    rows_by_series: dict[str, list[dict[str, float | str]]] = {}
    for series in SERIES:
        src_points = [(px, py) for key, px, py in pixel_points if key == series.key]
        rows = []
        for idx, ((x_raw, y_raw), (px, py)) in enumerate(zip(raw.get(series.key, []), src_points)):
            rows.append(
                {
                    "series": series.key,
                    "label": series.label.replace("$", ""),
                    "point_index": idx,
                    "cutoff_gev": float(x_raw),
                    "fraction": float(axis.y_min + y_raw / y_shift_scale),
                    "pixel_x": float(px),
                    "pixel_y": float(py),
                    "source_pad": axis.name,
                    "note": caveat,
                }
            )
        rows_by_series[series.key] = rows

    metadata = {
        "source_crop": str(source_crop),
        "raw_shifted_datathief_csv": str(raw_csv),
        "refs": refs,
        "y_shift_scale": y_shift_scale,
        "axis_pixels": {
            "x0": axis.x0,
            "x1": axis.x1,
            "y_top": axis.y_top,
            "y_bottom": axis.y_bottom,
            "y_min": axis.y_min,
            "y_max": axis.y_max,
        },
        "caveat": caveat,
        "n_points": {key: len(rows) for key, rows in rows_by_series.items()},
    }
    return rows_by_series, metadata


def write_combined_csv(path: Path, pads: dict[str, dict[str, list[dict[str, float | str]]]]) -> None:
    with path.open("w", newline="") as f:
        fieldnames = [
            "source_pad",
            "series",
            "label",
            "point_index",
            "cutoff_gev",
            "fraction",
            "pixel_x",
            "pixel_y",
            "note",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for pad_name in ("direct", "fragmentation"):
            for series in SERIES:
                for row in pads[pad_name][series.key]:
                    writer.writerow(row)


def plot_pad(
    rows_by_series: dict[str, list[dict[str, float | str]]],
    axis: AxisSpec,
    output: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.05,
        }
    )
    fig, ax = plt.subplots(figsize=(7.96, 7.22), dpi=100)
    marker_styles = {
        "pt10_15": {"marker": "o", "s": 54, "facecolors": "none", "linewidths": 1.35},
        "pt15_20": {"marker": "s", "s": 42, "facecolors": "none", "linewidths": 1.25},
        "pt25_30": {"marker": "^", "s": 28, "facecolors": None, "linewidths": 1.0},
    }
    for series in SERIES:
        rows = rows_by_series[series.key]
        xs = np.asarray([float(row["cutoff_gev"]) for row in rows], dtype=float)
        ys = np.asarray([float(row["fraction"]) for row in rows], dtype=float)
        order = np.argsort(xs)
        style = marker_styles[series.key]
        facecolors = style["facecolors"] if style["facecolors"] is not None else series.color
        ax.scatter(
            xs[order],
            ys[order],
            marker=style["marker"],
            s=style["s"],
            facecolors=facecolors,
            edgecolors=series.color,
            linewidths=style["linewidths"],
            label=series.label,
            zorder=3,
        )

    ax.set_xlim(0.0, 20.0)
    ax.set_ylim(axis.y_min, axis.y_max)
    ax.set_xlabel(r"$E_T^{iso}$ Cutoff [GeV]", fontsize=16)
    ax.set_ylabel("Fraction of Events", fontsize=16)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
    ax.grid(False)
    title = "Direct photons" if axis.name == "direct" else "Fragmentation photons"
    ax.text(0.06, 0.93, title, transform=ax.transAxes, fontsize=15)
    ax.text(0.06, 0.875, r"$R=0.3,\ |\eta^\gamma|<0.7,\ |z_{vtx}|<30$ cm", transform=ax.transAxes, fontsize=11)
    ax.text(0.06, 0.825, "DataThief export from PPG12 IAN Fig. 2 PNG", transform=ax.transAxes, fontsize=10.5)
    if axis.name == "direct":
        ax.text(
            0.06,
            0.775,
            "curves overlap in source; common visible curve used",
            transform=ax.transAxes,
            fontsize=10.5,
        )
    ax.legend(loc="lower right", frameon=False, fontsize=11.2, handlelength=2.0)
    fig.subplots_adjust(left=0.15, right=0.965, top=0.965, bottom=0.13)
    fig.savefig(output)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not args.source.exists():
        raise FileNotFoundError(args.source)
    args.outdir.mkdir(parents=True, exist_ok=True)

    rgb = np.asarray(Image.open(args.source).convert("RGB"))
    pads: dict[str, dict[str, list[dict[str, float | str]]]] = {}
    metadata: dict[str, object] = {}
    for name, axis in AXES.items():
        rows, meta = export_pad_points(args.source, rgb.copy(), args.outdir, axis)
        pads[name] = rows
        metadata[name] = meta

    combined_csv = args.outdir / "ppg12_fig2_datathief_points.csv"
    write_combined_csv(combined_csv, pads)

    direct_png = args.outdir / "ppg12_fig2_direct_datathief_overlay_796x722.png"
    frag_png = args.outdir / "ppg12_fig2_fragmentation_datathief_overlay_796x722.png"
    plot_pad(pads["direct"], AXES["direct"], direct_png)
    plot_pad(pads["fragmentation"], AXES["fragmentation"], frag_png)

    manifest = {
        "source_pdf": str(REPO / "usefulDocs/ppg12/current/PPG12_CURRENT_IAN_2026-05-21_v4.pdf"),
        "source_image": str(args.source),
        "datathief_jar": str(DATATHIEF_JAR),
        "datathief_jar_md5": jar_md5(),
        "combined_points_csv": str(combined_csv),
        "outputs": {
            "direct_png": str(direct_png),
            "fragmentation_png": str(frag_png),
        },
        "pads": metadata,
        "note": (
            "This is a plot-only digitization artifact for side-by-side IAN "
            "checking. Direct-photon pT curves are not independently separable "
            "from the rendered source image."
        ),
    }
    manifest_path = args.outdir / "ppg12_fig2_datathief_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
