#!/usr/bin/env python3
"""Build a follow-up slide for the reco-cluster ET leakage diagnostic."""

from __future__ import annotations

import csv
import json
import math
from io import BytesIO
from pathlib import Path
from textwrap import wrap

import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[4]
OUT = REPO / "dataOutput/slides/wp_gammajets_6_1_26/reco_cluster_et_leakage_20260604"
CSV_PATH = (
    REPO
    / "dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_jet1234_20260526"
    / "inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50.csv"
)
SUMMARY_PATH = (
    REPO
    / "dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_jet1234_20260526"
    / "inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50_summary.json"
)
SOURCE_SLIDE_PNG = (
    REPO
    / "dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_jet1234_20260526"
    / "inclusive_jet1234_reco_cluster_et_leakage_slide7_followup_12to50.png"
)

W, H = 2560, 1440
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
TIMES_BOLD_ITALIC = FONT_DIR / "Times New Roman Bold Italic.ttf"

INK = (18, 24, 38)
MUTED = (78, 86, 101)
LINE = (208, 216, 228)
BLUE = (33, 76, 195)
ORANGE = (238, 114, 18)
MAGENTA = (204, 43, 168)
GREEN = (35, 150, 55)
SIGNAL_ACCENT = (26, 111, 116)
SIGNAL_SOFT = (235, 248, 248)
SIGNAL_LINE = (174, 214, 216)
BACKGROUND_ACCENT = (92, 78, 152)
BACKGROUND_SOFT = (244, 241, 250)
BACKGROUND_LINE = (204, 198, 230)
RESULT_ACCENT = (42, 91, 117)
RESULT_SOFT = (241, 246, 249)

SAMPLES = ("Jet12", "Jet20", "Jet30", "Jet40")
COLORS = {
    "Jet12": BLUE,
    "Jet20": ORANGE,
    "Jet30": MAGENTA,
    "Jet40": GREEN,
    "Sum": INK,
}
MARKERS = {"Jet12": "o", "Jet20": "s", "Jet30": "^", "Jet40": "v", "Sum": "D"}


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    path = TIMES
    if bold and italic:
        path = TIMES_BOLD_ITALIC
    elif bold:
        path = TIMES_BOLD
    elif italic:
        path = TIMES_ITALIC
    try:
        return ImageFont.truetype(str(path), size)
    except OSError:
        return ImageFont.load_default()


F = {
    "title": font(60, bold=True),
    "subtitle": font(31),
    "panel_title": font(52, bold=True),
    "panel_subtitle": font(32),
    "card_title": font(45, bold=True),
    "card_body": font(36),
    "card_body_bold": font(36, bold=True),
    "body": font(31),
    "body_bold": font(31, bold=True),
    "small": font(27),
    "small_bold": font(27, bold=True),
}


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill, outline=LINE, width=2, radius=16) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], msg: str, fnt, fill=INK, anchor=None) -> None:
    draw.text(xy, msg, font=fnt, fill=fill, anchor=anchor)


def rich_line(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    parts: list[tuple[str, ImageFont.FreeTypeFont, tuple[int, int, int]]],
) -> int:
    cursor = x
    max_h = 0
    for msg, fnt, fill in parts:
        bbox = draw.textbbox((cursor, y), msg, font=fnt)
        draw.text((cursor, y), msg, font=fnt, fill=fill)
        cursor += bbox[2] - bbox[0]
        max_h = max(max_h, bbox[3] - bbox[1])
    return y + max_h


def wrapped(draw: ImageDraw.ImageDraw, xy: tuple[int, int], msg: str, fnt, fill=INK, width_chars=44, spacing=8) -> int:
    x, y = xy
    for para in msg.split("\n"):
        for line in wrap(para, width=width_chars) or [""]:
            draw.text((x, y), line, font=fnt, fill=fill)
            y += fnt.size + spacing
    return y


def paste_contained(canvas: Image.Image, img: Image.Image, box: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    bw, bh = x1 - x0, y1 - y0
    scale = min(bw / img.width, bh / img.height)
    nw, nh = int(img.width * scale), int(img.height * scale)
    resized = img.resize((nw, nh), Image.Resampling.LANCZOS)
    bg = Image.new("RGB", (bw, bh), "white")
    ox, oy = (bw - nw) // 2, (bh - nh) // 2
    bg.paste(resized, (ox, oy))
    canvas.paste(bg, (x0, y0))
    return (x0 + ox, y0 + oy, x0 + ox + nw, y0 + oy + nh)


def load_rows() -> list[dict]:
    rows: list[dict] = []
    with CSV_PATH.open(newline="") as handle:
        for row in csv.DictReader(handle):
            row = dict(row)
            for key in ("x_low", "x_high", "x_center", "x_err_low", "x_err_high", "weighted_entries", "weighted_error"):
                row[key] = float(row[key])
            rows.append(row)
    return rows


def pivot_rows(rows: list[dict]) -> list[dict]:
    by_bin: dict[tuple[float, float, float], dict] = {}
    for row in rows:
        key = (row["x_low"], row["x_high"], row["x_center"])
        target = by_bin.setdefault(
            key,
            {
                "x_low": row["x_low"],
                "x_high": row["x_high"],
                "x_center": row["x_center"],
                "x_err_low": row["x_err_low"],
                "x_err_high": row["x_err_high"],
            },
        )
        sample = row["sample"]
        target[f"weighted_entries_{sample}"] = row["weighted_entries"]
        target[f"weighted_error_{sample}"] = row["weighted_error"]

    out: list[dict] = []
    for _, row in sorted(by_bin.items()):
        sum_y = 0.0
        err2 = 0.0
        for sample in SAMPLES:
            y = float(row.get(f"weighted_entries_{sample}", 0.0))
            ey = float(row.get(f"weighted_error_{sample}", 0.0))
            row[f"weighted_entries_{sample}"] = y
            row[f"weighted_error_{sample}"] = ey
            sum_y += y
            err2 += ey * ey
        row["weighted_entries_Sum"] = sum_y
        row["weighted_error_Sum"] = math.sqrt(err2)
        for sample in SAMPLES:
            row[f"fraction_{sample}"] = row[f"weighted_entries_{sample}"] / sum_y if sum_y > 0.0 else 0.0
        out.append(row)
    return out


def pct(value: float) -> str:
    return f"{100.0 * value:.0f}%" if value == 0.0 else f"{100.0 * value:.1f}%"


def fraction(summary: dict, label: str, sample: str) -> float:
    for entry in summary["bins_of_interest"]:
        if entry["x_range"] == label:
            return float(entry["fractions"].get(sample, 0.0))
    return float("nan")


def render_plot(pivot: list[dict]) -> Image.Image:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(12.0, 8.25), dpi=120, facecolor="white")
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[3.3, 1.15],
        left=0.105,
        right=0.965,
        bottom=0.115,
        top=0.885,
        hspace=0.055,
    )
    ax = fig.add_subplot(gs[0])
    fax = fig.add_subplot(gs[1], sharex=ax)

    xs = np.array([r["x_center"] for r in pivot], dtype=float)
    xerr = np.array([[r["x_err_low"] for r in pivot], [r["x_err_high"] for r in pivot]], dtype=float)

    for sample in SAMPLES + ("Sum",):
        y = np.array([r[f"weighted_entries_{sample}"] for r in pivot], dtype=float)
        ey = np.array([r[f"weighted_error_{sample}"] for r in pivot], dtype=float)
        mask = y > 0.0
        color = np.array(COLORS[sample]) / 255.0
        label = sample if sample != "Sum" else "Weighted sum"
        ax.errorbar(
            xs[mask],
            y[mask],
            yerr=ey[mask],
            xerr=xerr[:, mask],
            fmt=MARKERS[sample],
            markersize=5.3 if sample != "Sum" else 5.9,
            markerfacecolor=color if sample != "Sum" else "white",
            markeredgecolor=color,
            markeredgewidth=1.0,
            linestyle="none",
            color=color,
            ecolor=color,
            elinewidth=0.9,
            capsize=1.6,
            label=label,
            zorder=5 if sample == "Sum" else 4,
        )

    for sample in SAMPLES:
        color = np.array(COLORS[sample]) / 255.0
        vals = np.array([r[f"fraction_{sample}"] for r in pivot], dtype=float)
        fax.plot(
            xs,
            vals,
            marker=MARKERS[sample],
            markersize=5.3,
            linestyle="none",
            color=color,
            label=sample,
            zorder=3,
        )

    ax.set_yscale("log")
    ax.set_xlim(12.0, 50.0)
    positive = []
    for row in pivot:
        for sample in SAMPLES + ("Sum",):
            y = row[f"weighted_entries_{sample}"]
            if y > 0.0:
                positive.append(y)
    ax.set_ylim(max(1.0, min(positive) * 0.45), max(positive) * 2.4)
    ax.set_ylabel("Weighted entries / 1 GeV bin", fontsize=15.0)
    ax.tick_params(labelbottom=False, labelsize=12.5, length=6)
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(which="major", color="0.84", linestyle=":", linewidth=0.75)
    ax.grid(which="minor", color="0.91", linestyle=":", linewidth=0.45)

    fax.set_ylim(-0.035, 1.05)
    fax.set_yticks([0.0, 0.25, 0.50, 0.75, 1.0])
    fax.set_ylabel("Fraction of\nweighted sum", fontsize=13.5)
    fax.set_xlabel(r"Reco photon-cluster $E_T$ ($p_T^\gamma$) [GeV]", fontsize=15.0, loc="right")
    fax.tick_params(labelsize=12.5, length=6)
    fax.grid(which="major", color="0.84", linestyle=":", linewidth=0.75)

    ax.text(
        0.965,
        0.980,
        r"$\it{\bf{sPHENIX}}$ Internal" + "\nPYTHIA8 embedded inclusive jet, 0-80%",
        transform=ax.transAxes,
        fontsize=14.0,
        ha="right",
        va="top",
        linespacing=1.05,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 2.0},
        zorder=20,
    )
    ax.text(
        0.040,
        0.110,
        "Jet12: 12-21 GeV     Jet20: 21-31 GeV\nJet30: 31-41 GeV     Jet40: >=41 GeV",
        transform=ax.transAxes,
        fontsize=17.0,
        ha="left",
        va="bottom",
        color=np.array(INK) / 255.0,
        linespacing=1.18,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "0.82", "alpha": 0.90, "linewidth": 0.8},
        zorder=20,
    )

    legend = ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.962, 0.785),
        frameon=False,
        fontsize=14.0,
        handlelength=1.8,
        labelspacing=0.30,
        borderaxespad=0.2,
        markerscale=1.05,
    )
    for item in legend.get_texts():
        item.set_fontsize(14.0)

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120, facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def add_plot_panel(canvas: Image.Image, draw: ImageDraw.ImageDraw, plot: Image.Image) -> None:
    box = (72, 184, 1654, 1348)
    x0, y0, x1, y1 = box
    rounded(draw, box, (255, 255, 255), outline=BACKGROUND_LINE, width=3, radius=18)
    draw.rounded_rectangle((x0 + 24, y0 + 22, x0 + 37, y0 + 100), radius=6, fill=BACKGROUND_ACCENT)
    text(draw, (x0 + 68, y0 + 14), "Embedded inclusive-jet source composition", F["panel_title"], INK)
    text(
        draw,
        (x0 + 70, y0 + 72),
        "Weighted reco-cluster energy in 1 GeV bins, with sample fractions below",
        F["panel_subtitle"],
        MUTED,
    )
    paste_contained(canvas, plot, (x0 + 50, y0 + 120, x1 - 50, y1 - 34))


def sample_chip(draw: ImageDraw.ImageDraw, x: int, y: int, name: str, window: str, color: tuple[int, int, int]) -> None:
    rounded(draw, (x, y, x + 304, y + 76), (255, 255, 255), outline=(224, 230, 239), radius=10)
    draw.ellipse((x + 20, y + 29, x + 38, y + 47), fill=color)
    text(draw, (x + 52, y + 12), name, F["small_bold"], color)
    text(draw, (x + 156, y + 15), window, F["small"], INK)


def add_info_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title_msg: str,
    accent: tuple[int, int, int],
    fill: tuple[int, int, int],
    outline: tuple[int, int, int],
    lines: list[tuple[str, str]],
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, fill, outline=outline, width=3, radius=18)
    draw.rounded_rectangle((x0 + 20, y0 + 24, x0 + 31, y1 - 24), radius=5, fill=accent)
    text(draw, (x0 + 54, y0 + 22), title_msg, F["card_title"], INK)
    y = y0 + 82
    for label, body in lines:
        y = rich_line(draw, x0 + 54, y, [(label, F["card_body_bold"], accent), (body, F["card_body"], MUTED)])
        y += 12


def add_fraction_card(draw: ImageDraw.ImageDraw, summary: dict) -> None:
    box = (1710, 800, 2488, 1098)
    x0, y0, x1, y1 = box
    rounded(draw, box, RESULT_SOFT, outline=(190, 208, 218), width=3, radius=18)
    draw.rounded_rectangle((x0 + 20, y0 + 24, x0 + 31, y1 - 24), radius=5, fill=RESULT_ACCENT)
    text(draw, (x0 + 54, y0 + 22), "Composition readout", F["card_title"], INK)

    entries = [
        ("22-24 GeV", [("Jet20", fraction(summary, "22-24", "Jet20"), ORANGE), ("Jet30", fraction(summary, "22-24", "Jet30"), MAGENTA)]),
        ("26-35 GeV", [("Jet30", fraction(summary, "26-35", "Jet30"), MAGENTA), ("Jet40", fraction(summary, "26-35", "Jet40"), GREEN)]),
        ("40-50 GeV", [("Jet30", fraction(summary, "40-50", "Jet30"), MAGENTA), ("Jet40", fraction(summary, "40-50", "Jet40"), GREEN)]),
    ]
    y = y0 + 88
    for label, pieces in entries:
        text(draw, (x0 + 54, y), label, F["body_bold"], INK)
        cursor = x0 + 230
        for name, value, color in pieces:
            rounded(draw, (cursor, y - 4, cursor + 214, y + 42), (255, 255, 255), outline=(224, 230, 239), radius=9)
            draw.ellipse((cursor + 14, y + 11, cursor + 28, y + 25), fill=color)
            text(draw, (cursor + 38, y + 3), f"{name} {pct(value)}", F["small_bold"], color)
            cursor += 232
        y += 58


def build() -> tuple[Path, Path, Path]:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = json.loads(SUMMARY_PATH.read_text())
    pivot = pivot_rows(load_rows())
    plot = render_plot(pivot)

    plot_png = OUT / "standardized_reco_cluster_et_leakage_plot.png"
    plot.save(plot_png)

    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)

    title = "Reco-cluster energy leakage check follows the stitched inclusive-jet ownership"
    text(draw, (72, 46), title, F["title"], INK)
    subtitle = "The diagnostic asks which embedded inclusive-jet sample populates the reconstructed photon-cluster-energy tail after weighting."
    text(draw, (76, 126), subtitle, F["subtitle"], MUTED)

    add_plot_panel(canvas, draw, plot)

    add_info_card(
        draw,
        (1710, 184, 2488, 414),
        "What is plotted",
        BACKGROUND_ACCENT,
        BACKGROUND_SOFT,
        BACKGROUND_LINE,
        [
            ("Object: ", "Reco photon-cluster energy",),
            ("Binning: ", "1 GeV bins, 12-50 GeV",),
            ("Weights: ", "Counted cross sections",),
        ],
    )

    add_info_card(
        draw,
        (1710, 444, 2488, 770),
        "Main readout",
        SIGNAL_ACCENT,
        SIGNAL_SOFT,
        SIGNAL_LINE,
        [
            ("Jet12: ", "leaves the high-energy tail.",),
            ("Jet20: ", "carries the 22-26 GeV transition.",),
            ("Jet30/40: ", "share the high-tail handoff.",),
        ],
    )

    add_fraction_card(draw, summary)

    rounded(draw, (1710, 1128, 2488, 1398), (255, 255, 255), outline=BACKGROUND_LINE, width=3, radius=18)
    text(draw, (1744, 1150), "Ownership windows", F["card_title"], INK)
    sample_chip(draw, 1744, 1210, "Jet12", "12-21", BLUE)
    sample_chip(draw, 2078, 1210, "Jet20", "21-31", ORANGE)
    sample_chip(draw, 1744, 1294, "Jet30", "31-41", MAGENTA)
    sample_chip(draw, 2078, 1294, "Jet40", ">=41", GREEN)

    png = OUT / "reco_cluster_et_leakage_followup_slide.png"
    script = OUT / "reco_cluster_et_leakage_followup_slide_script.md"
    manifest = OUT / "reco_cluster_et_leakage_followup_slide_manifest.json"
    canvas.save(png)

    script.write_text(
        """# WP GammaJets Follow-up Script - Reco-cluster Energy Leakage

This is the follow-up to the embedded stitching closure slide. The previous slide showed that the truth-filter ownership windows and cross-section weights give continuous stitched spectra. Here I am checking a related but different object: the reconstructed photon-cluster energy distribution in the embedded inclusive-jet background.

The plot shows the ABCD-summed reconstructed photon-cluster energy distribution from the inclusive Jet12, Jet20, Jet30, and Jet40 samples after applying the counted effective-cross-section weights. The upper panel is the weighted yield in true 1 GeV reconstructed-cluster-energy bins from 12 to 50 GeV. The lower panel shows the fraction of the weighted sum coming from each generator-filter sample.

The main thing to look for is whether a lower-threshold sample is leaking into a region where it should not dominate. Jet12 falls out above the low-energy region, which is what we want for the high-energy tail. Jet20 carries the transition around 22 to 26 GeV. Then Jet30 takes over through the 31 to 41 GeV ownership window, and Jet40 begins to contribute in the high-energy tail.

So this slide is not claiming a final physics result by itself. It is a sanity check on the embedded inclusive-jet background bookkeeping. Together with the stitched truth-filter spectrum, it supports the statement that the four inclusive embedded samples have non-overlapping ownership and a sensible reconstructed-energy composition after weighting.
""",
        encoding="utf-8",
    )

    manifest.write_text(
        json.dumps(
            {
                "png": str(png),
                "script": str(script),
                "standardized_plot_png": str(plot_png),
                "google_slides_mutation": False,
                "source_google_slide": {
                    "presentation_id": "1GAoEcN9UGOkUVxTgT0m1vRhghl3O1jd28LwSHg9klBs",
                    "slide_object_id": "g3a89b5423f1_0_192",
                    "image_title_readback": "inclusive_jet1234_reco_cluster_et_leakage_slide7_followup_12to50.png",
                },
                "source_png_from_slide": str(SOURCE_SLIDE_PNG),
                "source_csv": str(CSV_PATH),
                "source_summary": str(SUMMARY_PATH),
                "style_relation": "Follow-up to embedded stitching closure slide; uses matching Times typography, panel accents, and Jet12/20/30/40 color semantics.",
                "summary_bins_of_interest": summary["bins_of_interest"],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return png, script, manifest


if __name__ == "__main__":
    for path in build():
        print(path)
