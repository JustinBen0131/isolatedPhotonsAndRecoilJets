#!/usr/bin/env python3
"""Digitize the PPG12 IAN Fig. 3 screenshot and redraw its Fig. 3 contract.

This is a screenshot/DataThief reference only. The authoritative ROOT source is
still PPG12's MC_efficiency_noiso.root. The purpose here is to verify that the
visible top-panel points reproduce the Direct/Total lower-panel shape.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, LogLocator, MultipleLocator, NullFormatter
import numpy as np
from PIL import Image


REPO = Path(__file__).resolve().parents[3]
DEFAULT_SCREENSHOT = Path(
    "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/"
    "NSIRD_screencaptureui_FaGLMM/Screenshot 2026-07-02 at 12.51.34\u202fPM.png"
)
DEFAULT_OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217"
    / "fig3_datathief_reference"
)


@dataclass(frozen=True)
class Series:
    key: str
    label: str
    color: str
    marker: str
    min_area: int
    max_area: int


SERIES = [
    Series("total", "Total", "#e15b9a", "o", 10, 80),
    Series("direct", "Direct", "#2c7a19", "s", 16, 80),
    Series("frag", "Fragmentation", "#2b83ff", "^", 8, 80),
]

# PPG12 Fig. 3 crop calibration from the supplied top-panel screenshot.
X_PIXEL_MIN = 58.0
X_PIXEL_MAX = 383.0
Y_PIXEL_TOP = 3.0
Y_PIXEL_BOTTOM = 228.0
X_PHYS_MIN = 10.0
X_PHYS_MAX = 35.0
LOG10_Y_TOP = 8.0
LOG10_Y_BOTTOM = 4.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screenshot", type=Path, default=DEFAULT_SCREENSHOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def pixel_to_x(xpix: float) -> float:
    return X_PHYS_MIN + (xpix - X_PIXEL_MIN) / (X_PIXEL_MAX - X_PIXEL_MIN) * (X_PHYS_MAX - X_PHYS_MIN)


def pixel_to_y(ypix: float) -> float:
    log10_y = LOG10_Y_TOP - (ypix - Y_PIXEL_TOP) / (Y_PIXEL_BOTTOM - Y_PIXEL_TOP) * (LOG10_Y_TOP - LOG10_Y_BOTTOM)
    return 10.0**log10_y


def connected_components(mask: np.ndarray) -> list[dict[str, float]]:
    height, width = mask.shape
    seen = np.zeros(mask.shape, dtype=bool)
    components: list[dict[str, float]] = []
    ys, xs = np.where(mask)
    for y0, x0 in zip(ys, xs):
        if seen[y0, x0]:
            continue
        stack = [(int(y0), int(x0))]
        seen[y0, x0] = True
        pts: list[tuple[int, int]] = []
        while stack:
            y, x = stack.pop()
            pts.append((y, x))
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dy == 0 and dx == 0:
                        continue
                    yy = y + dy
                    xx = x + dx
                    if 0 <= yy < height and 0 <= xx < width and mask[yy, xx] and not seen[yy, xx]:
                        seen[yy, xx] = True
                        stack.append((yy, xx))
        arr = np.asarray(pts)
        components.append(
            {
                "area": float(len(pts)),
                "xpix": float(arr[:, 1].mean()),
                "ypix": float(arr[:, 0].mean()),
                "xmin": float(arr[:, 1].min()),
                "xmax": float(arr[:, 1].max()),
                "ymin": float(arr[:, 0].min()),
                "ymax": float(arr[:, 0].max()),
            }
        )
    return components


def color_masks(rgb: np.ndarray) -> dict[str, np.ndarray]:
    height, width = rgb.shape[:2]
    plot = np.zeros((height, width), dtype=bool)
    plot[int(Y_PIXEL_TOP) : int(Y_PIXEL_BOTTOM) + 1, int(X_PIXEL_MIN) : int(X_PIXEL_MAX) + 1] = True
    legend = np.zeros((height, width), dtype=bool)
    legend[:115, 200:330] = True
    r = rgb[:, :, 0]
    g = rgb[:, :, 1]
    b = rgb[:, :, 2]
    usable = plot & ~legend
    return {
        "total": (r > 170) & (g > 45) & (g < 170) & (b > 90) & usable,
        "direct": (g > 70) & (r < 140) & (b < 140) & usable,
        "frag": (b > 140) & (r < 140) & (g > 60) & usable,
    }


def extract_points(screenshot: Path) -> tuple[list[dict[str, float | str]], dict[str, str]]:
    image = Image.open(screenshot).convert("RGB")
    rgb = np.asarray(image)
    masks = color_masks(rgb)
    rows: list[dict[str, float | str]] = []
    expected_centers = np.arange(10.5, 35.0, 1.0)
    qa: dict[str, str] = {}
    for series in SERIES:
        raw = [
            comp
            for comp in connected_components(masks[series.key])
            if series.min_area <= comp["area"] <= series.max_area
        ]
        used: set[int] = set()
        found = 0
        for expected in expected_centers:
            candidates = []
            for idx, comp in enumerate(raw):
                if idx in used:
                    continue
                x = pixel_to_x(comp["xpix"])
                if abs(x - expected) > 0.35:
                    continue
                # Reject axis/tick fragments at the very bottom/right.
                if comp["ypix"] > Y_PIXEL_BOTTOM - 0.5 or comp["xpix"] >= X_PIXEL_MAX - 0.1:
                    continue
                candidates.append((abs(x - expected), comp))
            if not candidates:
                continue
            _, comp = min(candidates, key=lambda item: (item[0], -item[1]["area"]))
            used.add(raw.index(comp))
            x = pixel_to_x(comp["xpix"])
            y = pixel_to_y(comp["ypix"])
            rows.append(
                {
                    "source": "ppg12_ian_fig3_datathief",
                    "series": series.key,
                    "label": series.label,
                    "expected_xcenter_gev": float(expected),
                    "xcenter_gev": x,
                    "counts": y,
                    "log10_counts": math.log10(y),
                    "xpix": comp["xpix"],
                    "ypix": comp["ypix"],
                    "component_area_px": comp["area"],
                }
            )
            found += 1
        qa[f"{series.key}_points"] = str(found)
    return rows, qa


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def grouped(rows: list[dict[str, float | str]]) -> dict[str, list[dict[str, float | str]]]:
    out = {series.key: [] for series in SERIES}
    for row in rows:
        out[str(row["series"])].append(row)
    for key in out:
        out[key].sort(key=lambda row: float(row["expected_xcenter_gev"]))
    return out


def direct_total_ratio(rows: list[dict[str, float | str]]) -> list[dict[str, float]]:
    by = grouped(rows)
    total_by_x = {round(float(row["expected_xcenter_gev"]), 3): float(row["counts"]) for row in by["total"]}
    direct_by_x = {round(float(row["expected_xcenter_gev"]), 3): float(row["counts"]) for row in by["direct"]}
    ratio = []
    for x in sorted(set(total_by_x) & set(direct_by_x)):
        ratio.append({"xcenter_gev": x, "direct_over_total": direct_by_x[x] / total_by_x[x]})
    return ratio


def draw(path: Path, rows: list[dict[str, float | str]]) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 15,
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
        }
    )
    fig = plt.figure(figsize=(8.0, 8.89), dpi=100)
    ax = fig.add_axes([0.13, 0.421, 0.79, 0.507])
    rax = fig.add_axes([0.13, 0.100, 0.79, 0.292], sharex=ax)
    for axis in (ax, rax):
        axis.tick_params(which="both", top=True, right=True, labelsize=13)
        axis.minorticks_on()

    group = grouped(rows)
    for series in SERIES:
        series_rows = group[series.key]
        ax.errorbar(
            [float(row["expected_xcenter_gev"]) for row in series_rows],
            [float(row["counts"]) for row in series_rows],
            fmt=series.marker,
            markersize=5.4,
            markerfacecolor=series.color,
            markeredgecolor=series.color,
            ecolor=series.color,
            elinewidth=0.75,
            capsize=0,
            linestyle="none",
            label=series.label,
        )

    ratio = direct_total_ratio(rows)
    rax.plot(
        [row["xcenter_gev"] for row in ratio],
        [row["direct_over_total"] for row in ratio],
        color="#2c7a19",
        linewidth=1.4,
        drawstyle="steps-mid",
    )

    ax.set_yscale("log")
    ax.set_xlim(10.0, 35.0)
    ax.set_ylim(8.0e3, 1.3e8)
    rax.set_ylim(0.50, 1.00)
    ax.tick_params(labelbottom=False)
    ax.xaxis.set_major_locator(FixedLocator([10, 15, 20, 25, 30, 35]))
    ax.xaxis.set_minor_locator(MultipleLocator(1.0))
    rax.xaxis.set_major_locator(FixedLocator([10, 15, 20, 25, 30, 35]))
    rax.xaxis.set_minor_locator(MultipleLocator(1.0))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=60))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ratio_ticks = np.arange(0.50, 1.001, 0.05)
    rax.yaxis.set_major_locator(FixedLocator(ratio_ticks))
    rax.yaxis.set_major_formatter(
        FixedFormatter([f"{tick:.2f}".rstrip("0").rstrip(".") for tick in ratio_ticks])
    )
    rax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.set_ylabel("Counts", fontsize=15, labelpad=12)
    rax.set_ylabel("Direct/Total", fontsize=15, labelpad=11)
    rax.set_xlabel(r"$p_T$ [GeV]", fontsize=15, labelpad=6)
    rax.axhline(1.0, color="0.4", linestyle=(0, (4, 4)), linewidth=0.9)

    ax.text(0.065, 0.97, r"$\bf{\it{sPHENIX}}$ Simulation", transform=ax.transAxes, fontsize=11, va="top")
    ax.text(0.065, 0.92, "Photon Jet Samples", transform=ax.transAxes, fontsize=11, va="top")
    ax.text(0.25, 0.85, r"Pythia, $\sqrt{s}$=200 GeV", transform=ax.transAxes, fontsize=10, va="top")
    ax.text(0.25, 0.80, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=10, va="top")
    ax.text(0.25, 0.75, r"$R = 0.3,\ E_T^{iso} < 4$ GeV", transform=ax.transAxes, fontsize=10, va="top")
    ax.text(
        0.37,
        0.095,
        "DataThief Equivalence Check under same plotting constraints\n"
        "PPG12 SDCC data to reproduce not found",
        transform=ax.transAxes,
        fontsize=8.0,
        ha="center",
        va="bottom",
        color="0.18",
    )
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.60, 0.82),
        frameon=False,
        fontsize=10.5,
        borderaxespad=0.0,
        labelspacing=0.35,
        handletextpad=0.5,
    )
    fig.savefig(path)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows, qa = extract_points(args.screenshot)
    points_csv = args.out_dir / "ppg12_fig3_ian_datathief_points.csv"
    ratio_csv = args.out_dir / "ppg12_fig3_ian_datathief_direct_over_total.csv"
    png = args.out_dir / "ppg12_fig3_ian_datathief_reproduction.png"
    manifest = args.out_dir / "ppg12_fig3_ian_datathief_manifest.json"
    write_csv(points_csv, rows)
    ratio = direct_total_ratio(rows)
    with ratio_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["xcenter_gev", "direct_over_total"])
        writer.writeheader()
        writer.writerows(ratio)
    draw(png, rows)
    manifest.write_text(
        json.dumps(
            {
                "artifact": "PPG12 IAN Fig.3 DataThief reproduction from screenshot crop",
                "screenshot": str(args.screenshot),
                "axis_calibration": {
                    "x_pixel_min": X_PIXEL_MIN,
                    "x_pixel_max": X_PIXEL_MAX,
                    "y_pixel_top": Y_PIXEL_TOP,
                    "y_pixel_bottom": Y_PIXEL_BOTTOM,
                    "x_phys_min": X_PHYS_MIN,
                    "x_phys_max": X_PHYS_MAX,
                    "log10_y_top": LOG10_Y_TOP,
                    "log10_y_bottom": LOG10_Y_BOTTOM,
                },
                "qa": qa,
                "points_csv": str(points_csv),
                "ratio_csv": str(ratio_csv),
                "png": str(png),
                "caveat": "Screenshot/DataThief reference only; exact PPG12 SDCC ROOT remains authoritative.",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(json.dumps({"png": str(png), "points_csv": str(points_csv), "ratio_csv": str(ratio_csv), "manifest": str(manifest)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
