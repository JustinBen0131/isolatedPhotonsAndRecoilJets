#!/usr/bin/env python3
"""Reproduce the PPG12 IAN Fig. 37 iteration-stability plot.

The preferred input is a local CSV extracted read-only from Shuhang's SDCC ROOT
files using the same calculation as ppg12codeGit/plotting/plot_unfold_iter.C.
If that CSV is absent, the script falls back to the older screenshot-digitized
points and labels the output as such.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[3]
OUTDIR = REPO / "dataOutput/ppg12Parity/the76_ppg12_fig37_iteration_stability_reference"
SDCC_ROOT_CSV = OUTDIR / "ppg12_fig37_iteration_stability_sdcc_nomold_root_points.csv"


def load_root_points(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows: list[dict[str, str]] = []
    with path.open() as f:
        reader = csv.DictReader(f)
        rows.extend(reader)
    iteration = np.array([float(row["iteration"]) for row in rows], dtype=float)
    stat_uncertainty = np.array([float(row["total_relative_stat_uncertainty"]) for row in rows], dtype=float)
    relative_deviation = np.array([float(row["total_relative_deviation"]) for row in rows], dtype=float)
    stat_plus_dev = np.array([float(row["stat_plus_dev_quadrature"]) for row in rows], dtype=float)
    return iteration, stat_uncertainty, relative_deviation, stat_plus_dev


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    if SDCC_ROOT_CSV.exists():
        iteration, stat_uncertainty, relative_deviation, stat_plus_dev = load_root_points(SDCC_ROOT_CSV)
        csv_path = SDCC_ROOT_CSV
        png_path = OUTDIR / "ppg12_fig37_iteration_stability_sdcc_nomold_root_reproduction.png"
        numeric_source = "SDCC ROOT extraction from Photon_final_bdt_nomold.root + MC_response_bdt_nomold.root"
        source_status = "ROOT-derived; matches supplied IAN Fig. 37 iteration-1 scale"
    else:
        # Digitized from the supplied PPG12 IAN Fig. 37 screenshot.
        # The PPG12 macro draws:
        #   black: quadrature-summed relative statistical uncertainty
        #   blue: quadrature-summed relative deviation
        #   red:  quadrature of the two above
        iteration = np.arange(1, 11, dtype=float)
        relative_deviation = np.array([0.495, 0.088, 0.044, 0.030, 0.024, 0.020, 0.017, 0.014, 0.012, 0.011])
        stat_plus_dev = np.array([0.515, 0.218, 0.241, 0.272, 0.297, 0.316, 0.337, 0.355, 0.370, 0.386])
        stat_uncertainty = np.sqrt(np.maximum(stat_plus_dev**2 - relative_deviation**2, 0.0))

        csv_path = OUTDIR / "ppg12_fig37_iteration_stability_digitized_points.csv"
        with csv_path.open("w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "iteration",
                    "total_relative_stat_uncertainty",
                    "total_relative_deviation",
                    "stat_plus_dev_quadrature",
                    "source",
                ],
            )
            writer.writeheader()
            for i, stat, dev, quad in zip(iteration, stat_uncertainty, relative_deviation, stat_plus_dev):
                writer.writerow(
                    {
                        "iteration": int(i),
                        "total_relative_stat_uncertainty": f"{stat:.8g}",
                        "total_relative_deviation": f"{dev:.8g}",
                        "stat_plus_dev_quadrature": f"{quad:.8g}",
                        "source": "digitized_from_supplied_ppg12_ian_fig37_screenshot",
                    }
                )
        png_path = OUTDIR / "ppg12_fig37_iteration_stability_reproduction.png"
        numeric_source = "digitized from supplied PPG12 IAN Fig. 37 screenshot"
        source_status = "fallback only; SDCC ROOT extraction CSV was not present"

    # Minimal ROOT-like renderer using Pillow so this remains dependency-light.
    width = height = 720
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
        candidates = [
            "/System/Library/Fonts/Helvetica.ttc",
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
            "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
        ]
        for path in candidates:
            try:
                return ImageFont.truetype(path, size=size)
            except Exception:
                pass
        return ImageFont.load_default()

    f_axis = font(25)
    f_tick = font(23)
    f_text = font(22)
    f_small = font(18)
    f_bold = font(23, bold=True)
    f_sphenix = font(22, bold=True)

    left, right, top, bottom = 92, 675, 62, 612
    xmin, xmax = 0.0, 10.5
    ymin, ymax = 0.0, 0.7

    def xp(x: float) -> int:
        return int(round(left + (x - xmin) / (xmax - xmin) * (right - left)))

    def yp(y: float) -> int:
        return int(round(bottom - (y - ymin) / (ymax - ymin) * (bottom - top)))

    # Frame.
    draw.rectangle([left, top, right, bottom], outline="black", width=2)

    # Major/minor ticks. ROOT-like inward ticks on all sides.
    for x in np.arange(0, 10.1, 1.0):
        px = xp(float(x))
        draw.line([px, bottom, px, bottom - 18], fill="black", width=2)
        draw.line([px, top, px, top + 18], fill="black", width=2)
        label = str(int(x))
        tw = draw.textlength(label, font=f_tick)
        draw.text((px - tw / 2, bottom + 8), label, fill="black", font=f_tick)
    for x in np.arange(0.5, 10.1, 0.5):
        px = xp(float(x))
        draw.line([px, bottom, px, bottom - 9], fill="black", width=1)
        draw.line([px, top, px, top + 9], fill="black", width=1)

    for y in np.arange(0, 0.71, 0.1):
        py = yp(float(y))
        draw.line([left, py, left + 18, py], fill="black", width=2)
        draw.line([right, py, right - 18, py], fill="black", width=2)
        label = "0" if abs(y) < 1e-12 else f"{y:.1f}"
        tw = draw.textlength(label, font=f_tick)
        draw.text((left - tw - 12, py - 13), label, fill="black", font=f_tick)
    for y in np.arange(0.02, 0.70, 0.02):
        py = yp(float(y))
        draw.line([left, py, left + 9, py], fill="black", width=1)
        draw.line([right, py, right - 9, py], fill="black", width=1)

    # Axis labels.
    xlabel = "Iteration"
    draw.text((right - draw.textlength(xlabel, font=f_axis), bottom + 50), xlabel, fill="black", font=f_axis)
    ylabel = "√δ"
    y_label_img = Image.new("RGBA", (190, 44), (255, 255, 255, 0))
    y_label_draw = ImageDraw.Draw(y_label_img)
    y_label_draw.text((0, 0), ylabel, fill="black", font=f_axis)
    y_label_img = y_label_img.rotate(90, expand=True)
    image.paste(y_label_img, (22, top - 3), y_label_img)

    def marker(x: float, y: float, color: str, r: int = 4) -> None:
        px, py = xp(x), yp(y)
        draw.ellipse([px - r, py - r, px + r, py + r], fill=color, outline=color)

    for x, y in zip(iteration, stat_uncertainty):
        marker(float(x), float(y), "black", 4)
    for x, y in zip(iteration, relative_deviation):
        marker(float(x), float(y), "blue", 4)
    for x, y in zip(iteration, stat_plus_dev):
        marker(float(x), float(y), "red", 4)

    # Text block and legend.
    tx = xp(5.0)
    draw.text((tx, yp(0.675)), "sPHENIX", fill="black", font=f_sphenix)
    draw.text((tx + 118, yp(0.675)), "Internal", fill="black", font=f_text)
    draw.text((tx, yp(0.625)), "p+p √s=200 GeV", fill="black", font=f_text)
    draw.text((tx, yp(0.575)), "|η^γ| < 0.7", fill="black", font=f_text)

    legend_x, legend_y = xp(5.0), yp(0.44)
    legend_rows = [
        ("black", "total relative stat. uncertainty"),
        ("blue", "total relative deviation"),
        ("red", "stat. + dev. (quadrature)"),
    ]
    for i, (color, label) in enumerate(legend_rows):
        yy = legend_y + i * 29
        draw.line([legend_x - 8, yy, legend_x + 8, yy], fill=color, width=1)
        draw.ellipse([legend_x - 4, yy - 4, legend_x + 4, yy + 4], fill=color, outline=color)
        draw.text((legend_x + 20, yy - 13), label, fill="black", font=f_small)

    image.save(png_path)

    manifest = {
        "plot": str(png_path),
        "csv": str(csv_path),
        "provenance": {
            "ppg12_macro": "ppg12codeGit/plotting/plot_unfold_iter.C",
            "macro_source_root": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nomold.root",
            "macro_response_root": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_response_bdt_nomold.root",
            "macro_histograms": [
                "h_unfold_sub_leak_1",
                "h_unfold_sub_leak_2",
                "h_unfold_sub_leak_3",
                "h_unfold_sub_leak_4",
                "h_unfold_sub_leak_5",
                "h_unfold_sub_leak_6",
                "h_unfold_sub_leak_7",
                "h_unfold_sub_leak_8",
                "h_unfold_sub_leak_9",
                "h_unfold_sub_leak_10",
                "h_pT_truth_response_0",
            ],
            "source_status": source_status,
            "numeric_source": numeric_source,
            "calculation": "plot_unfold_iter.C equivalent: calcDelta=(h_this-h_ref)/h_ref; quadrature over bins 2..N-2",
        },
    }
    manifest_path = png_path.with_name(png_path.stem + "_manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(png_path)
    print(csv_path)
    print(manifest_path)


if __name__ == "__main__":
    main()
