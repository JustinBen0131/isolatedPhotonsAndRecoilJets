#!/usr/bin/env python3
"""Digitize PPG12 Fig. 13 black data markers from a PNG and compare to SDCC ROOT data.

This is a "data thief" cross-check: the PNG is treated as the reference image,
the top-panel axes are calibrated from the rendered tick marks, and only the
black filled data-marker centers are digitized.  The SDCC ROOT-extracted JSON is
then overlaid to verify that the local ROOT source corresponds to the IAN plot.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


DEFAULT_DIR = Path(
    "dataOutput/ppg12PhotonYield/"
    "ppg12_photon_yield_v1_data_20260620/"
    "shower_shape_reference_validation/fig13_e11_e33"
)
DEFAULT_JSON = DEFAULT_DIR / "ppg12_sdcc_fig13_e11_to_e33_histograms.json"
DEFAULT_PNG = Path(
    "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/"
    "NSIRD_screencaptureui_UNWpDZ/Screenshot 2026-06-29 at 10.11.06\u202fPM.png"
)
DEFAULT_OUT = DEFAULT_DIR / "ppg12_ian_png_digitized_black_data_vs_sdcc_overlay.png"
DEFAULT_SLIDE_OUT = DEFAULT_DIR / "ppg12_ian_png_digitized_black_data_vs_sdcc_overlay_slidefit_772x998.png"
DEFAULT_CSV = DEFAULT_DIR / "ppg12_ian_png_digitized_black_data_vs_sdcc_points.csv"
DEFAULT_SUMMARY = DEFAULT_DIR / "ppg12_ian_png_digitized_black_data_vs_sdcc_summary.json"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--image", type=Path, default=DEFAULT_PNG)
    ap.add_argument("--sdcc-json", type=Path, default=DEFAULT_JSON)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--slide-output", type=Path, default=DEFAULT_SLIDE_OUT)
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    ap.add_argument(
        "--no-slide-output",
        action="store_true",
        help="Only write --output, not the 830x815 slide-fit PNG.",
    )
    return ap.parse_args()


def clusters(indices: np.ndarray) -> list[np.ndarray]:
    if len(indices) == 0:
        return []
    out = []
    cur = [int(indices[0])]
    for idx in indices[1:]:
        idx = int(idx)
        if idx <= cur[-1] + 1:
            cur.append(idx)
        else:
            out.append(np.asarray(cur, dtype=float))
            cur = [idx]
    out.append(np.asarray(cur, dtype=float))
    return out


def weighted_cluster_center(values: np.ndarray, weights: np.ndarray) -> float:
    total = float(weights.sum())
    if total <= 0:
        return float(values.mean())
    return float((values * weights).sum() / total)


def calibrate_axes(rgb: np.ndarray) -> dict[str, float]:
    """Return pixel-to-data calibration for the top panel.

    ROOT-style tick marks make the axis calibration more reliable than using
    the frame top: the top frame is above the 0.18 y tick, so use the 0.18 and
    0.02 major tick rows plus the bottom axis row inferred from the tick train.
    """

    gray = np.dot(rgb[..., :3], [0.299, 0.587, 0.114])
    dark = gray < 40
    col_counts = dark.sum(axis=0)
    row_counts = dark.sum(axis=1)

    vertical_cols = np.where(col_counts > 900)[0]
    x_groups = clusters(vertical_cols)
    x_left = weighted_cluster_center(x_groups[0], col_counts[x_groups[0].astype(int)])
    x_right = weighted_cluster_center(x_groups[-1], col_counts[x_groups[-1].astype(int)])

    horizontal_rows = np.where(row_counts > 500)[0]
    y_groups = clusters(horizontal_rows)
    y_top_frame = float(y_groups[0].mean())
    y_bottom_frame = float(y_groups[1].mean())

    # Major y ticks on the left axis. Major ticks are wider than minor ticks.
    tick_region = dark[:, int(round(x_left)) + 1 : int(round(x_left)) + 24]
    tick_counts = tick_region.sum(axis=1)
    tick_rows = np.where((tick_counts >= 20) & (np.arange(len(tick_counts)) > y_top_frame))[0]
    tick_rows = tick_rows[tick_rows < y_bottom_frame + 1]
    tick_groups = clusters(tick_rows)
    tick_centers = [float(g.mean()) for g in tick_groups]

    # The usable major tick sequence is 0.18, 0.16, ..., 0.02.  The 0.00 tick
    # is the frame bottom.  Minor ticks and a tiny near-axis data marker can add
    # extra rows near the bottom, so fit the dominant major-tick spacing.
    majors = tick_centers[:9]
    y_at_018 = majors[0]
    y_at_002 = majors[-1]
    slope_pix_per_y = (y_at_002 - y_at_018) / (0.02 - 0.18)
    y_at_zero = y_at_002 - slope_pix_per_y * 0.02

    return {
        "x_left": float(x_left),
        "x_right": float(x_right),
        "y_top_frame": y_top_frame,
        "y_bottom_frame": y_bottom_frame,
        "y_at_018": float(y_at_018),
        "y_at_zero": float(y_at_zero),
    }


def y_to_data(pixel_y: float, cal: dict[str, float]) -> float:
    return (cal["y_at_zero"] - pixel_y) / (cal["y_at_zero"] - cal["y_at_018"]) * 0.18


def x_to_pixel(x_value: float, cal: dict[str, float]) -> float:
    return cal["x_left"] + x_value * (cal["x_right"] - cal["x_left"])


def filled_circle_score(dark: np.ndarray, cx: int, cy: int, radius: int = 9) -> tuple[float, float]:
    yy, xx = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    disk = xx * xx + yy * yy <= radius * radius
    patch = dark[cy - radius : cy + radius + 1, cx - radius : cx + radius + 1]
    if patch.shape != disk.shape:
        return -1.0, -1.0
    central = dark[cy - 4 : cy + 5, cx - 4 : cx + 5]
    return float(patch[disk].mean()), float(central.mean())


def digitize_black_markers(
    rgb: np.ndarray, centers: list[float], cal: dict[str, float]
) -> list[dict[str, float | str]]:
    """Find black filled data-marker centers near the known histogram bin centers."""

    gray = np.dot(rgb[..., :3], [0.299, 0.587, 0.114])
    dark = gray < 60

    out = []
    for center in centers:
        expected_x = x_to_pixel(center, cal)
        xlo = int(max(cal["x_left"] + 4, round(expected_x) - 15))
        xhi = int(min(cal["x_right"] - 4, round(expected_x) + 15))
        ylo = int(cal["y_at_018"] + 180)
        yhi = int(cal["y_at_zero"] - 2)

        best = None
        for cy in range(ylo, yhi + 1):
            for cx in range(xlo, xhi + 1):
                density, central = filled_circle_score(dark, cx, cy)
                if density < 0:
                    continue
                # Filled black markers score near 2 here.  Text glyph holes and
                # axis/tick fragments fail the central-density term.
                score = density + central - 0.006 * abs(cx - expected_x)
                if best is None or score > best["score"]:
                    best = {
                        "score": score,
                        "pixel_x": float(cx),
                        "pixel_y": float(cy),
                        "radius_px": 9.0,
                        "density": density,
                        "central_density": central,
                        "method": "disk_filled_marker",
                    }

        if best and best["density"] > 0.75 and best["central_density"] > 0.80:
            chosen = best
        else:
            # Tail bins sit on the x-axis in this figure, so the bottom half of
            # the marker can be hidden by the frame.  Use a constrained local
            # black-pixel component as an explicit axis-edge estimate.
            x_c = int(round(expected_x))
            xlo = max(0, x_c - 16)
            xhi = min(dark.shape[1] - 1, x_c + 16)
            ylo = int(max(cal["y_at_zero"] - 55, cal["y_at_018"]))
            yhi = int(min(cal["y_at_zero"] - 1, dark.shape[0] - 1))
            sub = dark[ylo : yhi + 1, xlo : xhi + 1]
            yy, xx = np.where(sub)
            if len(xx) == 0:
                raise RuntimeError(f"No marker candidate found for x={center:.2f}")
            # Pick the densest local disk if available; otherwise use the lower
            # edge of the visible dark component.  This path is only expected
            # for the last two tail markers close to y=0.
            local_best = None
            for cy in range(max(ylo + 8, int(cal["y_at_zero"] - 20)), yhi + 1):
                for cx in range(xlo + 8, xhi - 7):
                    density, central = filled_circle_score(dark, cx, cy)
                    score = density + 0.4 * central - 0.006 * abs(cx - expected_x)
                    if local_best is None or score > local_best["score"]:
                        local_best = {
                            "score": score,
                            "pixel_x": float(cx),
                            "pixel_y": float(cy),
                            "density": density,
                            "central_density": central,
                        }
            if local_best and local_best["density"] > 0.45:
                clipped_center_y = float(local_best["pixel_y"])
                clipped_center_x = float(local_best["pixel_x"])
                density = float(local_best["density"])
                central = float(local_best["central_density"])
            else:
                clipped_center_y = min(cal["y_at_zero"] - 1.0, ylo + float(yy.max()) - 0.5)
                clipped_center_x = xlo + float(xx.mean())
                density = float(len(xx))
                central = -1.0
            chosen = {
                "pixel_x": clipped_center_x,
                "pixel_y": float(clipped_center_y),
                "radius_px": 9.0,
                "density": density,
                "central_density": central,
                "method": "axis_edge_component_estimate",
            }
        out.append(
            {
                "bin_center": float(center),
                "pixel_x": float(chosen["pixel_x"]),
                "pixel_y": float(chosen["pixel_y"]),
                "digitized_value": float(y_to_data(float(chosen["pixel_y"]), cal)),
                "method": str(chosen["method"]),
            }
        )
    return out


def write_tables(
    rows: list[dict[str, float | str]],
    sdcc_data: dict,
    cal: dict[str, float],
    args: argparse.Namespace,
) -> None:
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "bin_center",
                "sdcc_value",
                "sdcc_error",
                "digitized_value",
                "digitized_over_sdcc",
                "pixel_x",
                "pixel_y",
                "method",
            ],
        )
        writer.writeheader()
        for row, sdcc_y, sdcc_e in zip(rows, sdcc_data["values"], sdcc_data["errors"]):
            sdcc_y = float(sdcc_y)
            writer.writerow(
                {
                    "bin_center": f"{float(row['bin_center']):.6g}",
                    "sdcc_value": f"{sdcc_y:.10g}",
                    "sdcc_error": f"{float(sdcc_e):.10g}",
                    "digitized_value": f"{float(row['digitized_value']):.10g}",
                    "digitized_over_sdcc": f"{float(row['digitized_value']) / sdcc_y:.10g}"
                    if sdcc_y > 0
                    else "",
                    "pixel_x": f"{float(row['pixel_x']):.3f}",
                    "pixel_y": f"{float(row['pixel_y']):.3f}",
                    "method": row["method"],
                }
            )

    summary = {
        "source_image": str(args.image),
        "sdcc_json": str(args.sdcc_json),
        "output_png": str(args.output),
        "output_csv": str(args.csv),
        "axis_calibration": cal,
        "rows": rows,
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def make_plot(
    rows: list[dict[str, float | str]],
    sdcc_data: dict,
    output: Path,
    *,
    slide_fit: bool = False,
) -> None:
    x = np.asarray(sdcc_data["centers"], dtype=float)
    sdcc_y = np.asarray(sdcc_data["values"], dtype=float)
    sdcc_err = np.asarray(sdcc_data["errors"], dtype=float)
    digit_y = np.asarray([float(r["digitized_value"]) for r in rows], dtype=float)
    ratio = digit_y / sdcc_y

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 15 if slide_fit else 16,
            "axes.linewidth": 1.2 if slide_fit else 1.4,
            "xtick.major.size": 6 if slide_fit else 7,
            "ytick.major.size": 6 if slide_fit else 7,
            "xtick.minor.size": 3 if slide_fit else 4,
            "ytick.minor.size": 3 if slide_fit else 4,
        }
    )
    fig_size = (7.72, 9.98) if slide_fit else (10.5, 8.2)
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=fig_size,
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.04},
    )

    ax.errorbar(
        x,
        sdcc_y,
        yerr=sdcc_err,
        fmt="o",
        color="black",
        ms=4.2 if slide_fit else 5.8,
        lw=1.0 if slide_fit else 1.2,
        capsize=0,
        label="SDCC ROOT data",
        zorder=3,
    )
    ax.plot(
        x,
        digit_y,
        "s",
        color="#d62728",
        markerfacecolor="none",
        markeredgewidth=1.35 if slide_fit else 1.7,
        ms=5.0 if slide_fit else 6.4,
        label="Data-thief from IAN PNG",
        zorder=4,
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.18)
    ax.set_ylabel("normalized counts", fontsize=17 if slide_fit else None)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    header_fs = 15 if slide_fit else 18
    label_fs = 11 if slide_fit else 14
    ax.text(
        0.050,
        0.92,
        "sPHENIX",
        transform=ax.transAxes,
        fontsize=header_fs,
        fontstyle="italic",
        fontweight="bold",
    )
    ax.text(0.235, 0.92, "Internal", transform=ax.transAxes, fontsize=header_fs)
    ax.text(0.050, 0.845, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=label_fs)
    ax.text(0.050, 0.780, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=label_fs)
    ax.text(0.050, 0.715, r"$22<p_T<28$ GeV", transform=ax.transAxes, fontsize=label_fs)
    ax.text(0.050, 0.650, "w/o nbkg cut", transform=ax.transAxes, fontsize=label_fs)
    ax.legend(
        loc="upper right",
        frameon=False,
        fontsize=14.0 if slide_fit else 15.5,
        handlelength=1.4 if slide_fit else 1.6,
    )

    rax.axhline(1.0, color="black", lw=1.0, ls=(0, (4, 4)))
    rax.plot(x, ratio, "o", color="#d62728", ms=4.0 if slide_fit else 5.0)
    rax.set_xlim(0.0, 1.0)
    rax.set_ylim(0.85, 1.15)
    rax.set_ylabel("PNG / SDCC", fontsize=15 if slide_fit else 14)
    rax.set_xlabel("e11_to_e33", fontsize=17 if slide_fit else None)
    rax.minorticks_on()
    rax.tick_params(which="both", direction="in", top=True, right=True)
    if slide_fit:
        ax.tick_params(axis="both", which="major", labelsize=15)
        rax.tick_params(axis="both", which="major", labelsize=15)

    if slide_fit:
        fig.subplots_adjust(left=0.145, right=0.985, top=0.985, bottom=0.115)
        dpi = 100
        bbox = None
    else:
        dpi = 220
        bbox = "tight"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches=bbox)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    rgb = np.asarray(Image.open(args.image).convert("RGB"))
    with args.sdcc_json.open() as f:
        payload = json.load(f)
    sdcc_data = payload["data"]
    centers = [float(x) for x in sdcc_data["centers"]]
    cal = calibrate_axes(rgb)
    rows = digitize_black_markers(rgb, centers, cal)
    write_tables(rows, sdcc_data, cal, args)
    make_plot(rows, sdcc_data, args.output)
    if not args.no_slide_output:
        make_plot(rows, sdcc_data, args.slide_output, slide_fit=True)
    print(args.output)
    if not args.no_slide_output:
        print(args.slide_output)
    print(args.csv)
    print(args.summary)


if __name__ == "__main__":
    main()
