#!/usr/bin/env python3
"""Build a full-slide checkpoint for the PPG12 Fig.29 purity ratio diagnostic."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DIAG_DIR = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_ratio_diagnostic"
)
RATIO_CSV = DIAG_DIR / "current_vs_ppg12_fig29_purity_ratio_points.csv"
OUT_DIR = DIAG_DIR / "slide_checkpoint"
OUT_PNG = OUT_DIR / "ppg12_fig29_purity_ratio_checkpoint_slide.png"
OUT_SCRIPT = OUT_DIR / "ppg12_fig29_purity_ratio_checkpoint_speaker_script.md"
OUT_MANIFEST = OUT_DIR / "ppg12_fig29_purity_ratio_checkpoint_manifest.json"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def read_rows() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with RATIO_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append({key: float(value) for key, value in row.items()})
    return rows


def render_plot(rows: list[dict[str, float]], out: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
            "axes.titlesize": 18,
            "axes.labelsize": 20,
            "xtick.labelsize": 17,
            "ytick.labelsize": 17,
            "legend.fontsize": 14,
        }
    )
    x = np.array([r["pt_center"] for r in rows])
    xerr = np.array([(r["pt_hi"] - r["pt_lo"]) / 2.0 for r in rows])

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11.0, 8.35),
        dpi=220,
        sharex=True,
        gridspec_kw={"height_ratios": [1.9, 1.05], "hspace": 0.055},
    )
    ax = axes[0]
    ax.errorbar(
        x - 0.06,
        [r["ppg12_raw"] for r in rows],
        xerr=xerr,
        yerr=[r["ppg12_raw_err"] for r in rows],
        fmt="o",
        color="black",
        ms=7.4,
        lw=1.8,
        capsize=0,
        label="PPG12 raw",
    )
    ax.errorbar(
        x + 0.06,
        [r["current_raw"] for r in rows],
        xerr=xerr,
        yerr=[r["current_raw_err"] for r in rows],
        fmt="D",
        mfc="white",
        mec="#0072B2",
        mew=1.8,
        ecolor="#0072B2",
        color="#0072B2",
        ms=7.2,
        lw=1.8,
        capsize=0,
        label="Current pp raw",
    )
    ax.errorbar(
        x - 0.06,
        [r["ppg12_corrected"] for r in rows],
        xerr=xerr,
        yerr=[r["ppg12_corrected_err"] for r in rows],
        fmt="^",
        color="#D55E00",
        ms=7.8,
        lw=1.8,
        capsize=0,
        label="PPG12 leakage corrected",
    )
    ax.errorbar(
        x + 0.06,
        [r["current_corrected"] for r in rows],
        xerr=xerr,
        yerr=[r["current_corrected_err"] for r in rows],
        fmt="s",
        mfc="white",
        mec="#009E73",
        mew=1.9,
        ecolor="#009E73",
        color="#009E73",
        ms=7.4,
        lw=1.8,
        capsize=0,
        label="Current pp leakage corrected",
    )
    ax.set_ylabel("Purity")
    ax.set_ylim(0.0, 1.16)
    ax.grid(True, color="0.80", alpha=0.38, lw=0.7)
    ax.legend(loc="upper left", frameon=False, ncol=2, columnspacing=1.0, handletextpad=0.5)
    ax.text(
        0.982,
        0.105,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p{+}p$, $\sqrt{s}=200$ GeV",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=15,
    )

    ratio_ax = axes[1]
    ratio_ax.axhline(1.0, color="0.15", lw=1.3)
    ratio_ax.plot(
        x,
        [r["current_raw_over_ppg12"] for r in rows],
        "D-",
        color="black",
        lw=2.25,
        ms=7.0,
        label="Raw current / PPG12",
    )
    ratio_ax.plot(
        x,
        [r["current_corrected_over_ppg12"] for r in rows],
        "^-",
        color="#D55E00",
        lw=2.25,
        ms=7.4,
        label="Leakage corrected current / PPG12",
    )
    ratio_ax.set_ylabel("Ratio")
    ratio_ax.set_xlabel(r"Cluster $E_T$ [GeV]")
    ratio_ax.set_xlim(9.4, 36.6)
    ratio_ax.set_ylim(0.2, 1.18)
    ratio_ax.grid(True, color="0.80", alpha=0.38, lw=0.7)
    ratio_ax.legend(loc="lower left", frameon=False)
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def draw_wrapped_text(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_width: int,
    text_font: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int],
    line_gap: int = 10,
) -> int:
    x, y = xy
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if draw.textbbox((0, 0), trial, font=text_font)[2] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    cursor = y
    for line in lines:
        draw.text((x, cursor), line, font=text_font, fill=fill)
        cursor += text_font.size + line_gap
    return cursor


def draw_arrow_bullet(
    draw: ImageDraw.ImageDraw,
    text: str,
    y: int,
    max_width: int,
    body_font: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int],
) -> int:
    x = 1804
    arrow_x = x
    arrow_y = y + 14
    draw.polygon(
        [(arrow_x, arrow_y), (arrow_x + 38, arrow_y + 21), (arrow_x, arrow_y + 42)],
        fill=(28, 91, 150),
    )
    return draw_wrapped_text(draw, text, (x + 62, y), max_width - 62, body_font, fill, line_gap=14) + 34


def crop_white_margin(image: Image.Image, margin: int = 12) -> Image.Image:
    arr = np.asarray(image)
    non_white = np.any(arr < 248, axis=2)
    rows = np.where(non_white.any(axis=1))[0]
    cols = np.where(non_white.any(axis=0))[0]
    if not len(rows) or not len(cols):
        return image
    left = max(int(cols[0]) - margin, 0)
    top = max(int(rows[0]) - margin, 0)
    right = min(int(cols[-1]) + margin, image.width - 1)
    bottom = min(int(rows[-1]) + margin, image.height - 1)
    return image.crop((left, top, right + 1, bottom + 1))


def compose_slide(rows: list[dict[str, float]], plot_path: Path) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    canvas = Image.new("RGB", (2560, 1440), "white")
    draw = ImageDraw.Draw(canvas)

    title_font = font(82, bold=True)
    body_font = font(43)

    title = "PPG12 Purity Checkpoint"
    draw.text((110, 54), title, font=title_font, fill=(20, 20, 20))

    plot = crop_white_margin(Image.open(plot_path).convert("RGB"), margin=8)
    plot.thumbnail((1715, 1215), Image.Resampling.LANCZOS)
    plot_x, plot_y = 60, 175
    canvas.paste(plot, (plot_x, plot_y))

    bullets = [
        "Raw purity is close: current/PPG12 = 0.80-0.97 through 10-26 GeV.",
        "Leakage is the main open issue: the current leakage shift is only 0.32-0.41 of PPG12 in stable bins.",
        "B and D signal-sideband leakage are near zero relative to PPG12; C is much closer.",
        "Next: reproduce PPG12 sideband semantics and toy-throw errors before claiming closure.",
    ]
    for y, bullet in zip((185, 435, 715, 1000), bullets):
        draw_arrow_bullet(draw, bullet, y, 665, body_font, (22, 22, 22))

    canvas.save(OUT_PNG)


def write_script(rows: list[dict[str, float]]) -> None:
    stable = [r for r in rows if r["pt_hi"] <= 26]
    raw_min = min(r["current_raw_over_ppg12"] for r in stable)
    raw_max = max(r["current_raw_over_ppg12"] for r in stable)
    shift_min = min(r["current_shift_over_ppg12_shift"] for r in stable)
    shift_max = max(r["current_shift_over_ppg12_shift"] for r in stable)
    OUT_SCRIPT.write_text(
        f"""# Speaker Script: PPG12 Fig.29 Purity Checkpoint

This is the current checkpoint against the exact PPG12 Fig.29 source graphs, not a final closure claim.

The main thing to say is that the raw purity is not wildly off. In the stable 10 to 26 GeV bins, our raw points are about {raw_min:.2f} to {raw_max:.2f} of the PPG12 values, which means the data selection and object chain are in the right broad neighborhood.

The important remaining difference is the leakage correction. PPG12 moves the points up much more strongly than our current correction does; our leakage shift is only about {shift_min:.2f} to {shift_max:.2f} of PPG12 in the stable bins. The fraction diagnostic points specifically to the B and D signal-sideband leakage terms being near zero compared with PPG12, while the C leakage term is much closer.

So the next debugging target is not another blind production pass. It is to reproduce the PPG12 leakage-sideband filling and the PPG12 uncertainty procedure: the isolation-sideband/sliding-window implementation details, plus their twenty-thousand-toy error propagation and Gaussian fit procedure.

The useful takeaway is that the infrastructure is now producing the correct object family and binning, and the discrepancy has been localized enough to drive the next code comparison.
""",
        encoding="utf-8",
    )


def write_manifest(rows: list[dict[str, float]], plot_path: Path) -> None:
    stable = [r for r in rows if r["pt_hi"] <= 26]
    OUT_MANIFEST.write_text(
        json.dumps(
            {
                "slide_png": str(OUT_PNG),
                "speaker_script": str(OUT_SCRIPT),
                "input_ratio_csv": str(RATIO_CSV),
                "plot_panel": str(plot_path),
                "ppg12_source": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root",
                "current_source": str(
                    REPO
                    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/purity_fig29_comparison/current_pp_photon_yield_purity_points.csv"
                ),
                "stable_bin_hi_max": 26,
                "raw_ratio_range_stable": [
                    min(r["current_raw_over_ppg12"] for r in stable),
                    max(r["current_raw_over_ppg12"] for r in stable),
                ],
                "leakage_shift_ratio_range_stable": [
                    min(r["current_shift_over_ppg12_shift"] for r in stable),
                    max(r["current_shift_over_ppg12_shift"] for r in stable),
                ],
                "interpretation": "checkpoint diagnostic; current pp output is not numerically closed to PPG12 Fig.29 yet",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_rows()
    plot_path = OUT_DIR / "ppg12_fig29_purity_ratio_two_panel.png"
    render_plot(rows, plot_path)
    compose_slide(rows, plot_path)
    write_script(rows)
    write_manifest(rows, plot_path)
    print(OUT_PNG)
    print(OUT_SCRIPT)


if __name__ == "__main__":
    main()
