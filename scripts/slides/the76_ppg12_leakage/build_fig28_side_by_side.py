#!/usr/bin/env python3
"""Build a 1x2 slide candidate comparing PPG12 Fig. 28 with current output."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
    "ppg12_photon_yield_v1_data_20260620/purity_fig29_comparison/"
    "ppg12_ratio_diagnostic/leakage_components"
)
CSV_PATH = BASE / "current_vs_ppg12_fig29_leakage_components_bcd.csv"
REF_DIR = BASE / "ppg12_reference_fig28_sdcc"
PPG12_FIG28_PNG = REF_DIR / "analysis_note_leakage_fraction_et_sim_p1.png"
OUT_DIR = BASE / "slide_candidates"
RHS_PNG = OUT_DIR / "the76_recoiljets_current_leakage_ppg12_style_20260629.png"
SLIDE_PNG = OUT_DIR / "the76_fig28_ppg12_vs_recoiljets_current_20260629.png"
MANIFEST = OUT_DIR / "the76_fig28_ppg12_vs_recoiljets_current_20260629_manifest.json"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
TIMES_BOLD_ITALIC = FONT_DIR / "Times New Roman Bold Italic.ttf"


def font(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(path), size=size)


F_TITLE = font(TIMES_BOLD, 76)
F_PANEL = font(TIMES_BOLD, 34)
F_NOTE = font(TIMES, 38)
F_SMALL = font(TIMES, 25)
F_PLOT_LARGE = font(TIMES, 70)
F_PLOT_MED = font(TIMES, 48)
F_PLOT_SMALL = font(TIMES, 44)
F_PLOT_SMALL_BOLD = font(TIMES_BOLD, 44)
F_PLOT_SMALL_BOLD_ITALIC = font(TIMES_BOLD_ITALIC, 44)


def read_rows() -> list[dict[str, float]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def sha256(path: Path) -> str:
    return subprocess.check_output(["shasum", "-a", "256", str(path)], text=True).split()[0]


def map_xy(
    x: float,
    y: float,
    rect: tuple[int, int, int, int],
    xmin: float = 10.0,
    xmax: float = 35.0,
    ymin: float = 0.0,
    ymax: float = 1.3,
) -> tuple[int, int]:
    x0, y0, x1, y1 = rect
    px = x0 + (x - xmin) / (xmax - xmin) * (x1 - x0)
    py = y1 - (y - ymin) / (ymax - ymin) * (y1 - y0)
    return round(px), round(py)


def draw_rotated_text(
    img: Image.Image,
    text: str,
    font_obj: ImageFont.FreeTypeFont,
    center: tuple[int, int],
    *,
    angle: int = 90,
) -> None:
    tmp = Image.new("RGBA", (700, 110), (255, 255, 255, 0))
    d = ImageDraw.Draw(tmp)
    bbox = d.textbbox((0, 0), text, font=font_obj)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    d.text(((tmp.width - tw) // 2, (tmp.height - th) // 2 - bbox[1]), text, font=font_obj, fill=(0, 0, 0, 255))
    rot = tmp.rotate(angle, expand=True)
    img.alpha_composite(rot, (center[0] - rot.width // 2, center[1] - rot.height // 2))


def line_with_width(draw: ImageDraw.ImageDraw, points: list[tuple[int, int]], color: tuple[int, int, int], width: int) -> None:
    if len(points) >= 2:
        draw.line(points, fill=color, width=width, joint="curve")


def step_points(edges: list[float], values: list[float], rect: tuple[int, int, int, int]) -> list[tuple[int, int]]:
    pts: list[tuple[int, int]] = []
    for idx, val in enumerate(values):
        lo = max(edges[idx], 10.0)
        hi = min(edges[idx + 1], 35.0)
        if hi <= 10.0 or lo >= 35.0:
            continue
        p0 = map_xy(lo, val, rect)
        p1 = map_xy(hi, val, rect)
        if not pts:
            pts.append(p0)
        else:
            pts.append(p0)
        pts.append(p1)
    return pts


def draw_current_ppg12_style(rows: list[dict[str, float]]) -> None:
    """Render current RecoilJets leakage with a PPG12-like step-line plot."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (1701, 1632), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)

    plot = (276, 85, 1628, 1325)
    edges = [rows[0]["pt_lo"]] + [r["pt_hi"] for r in rows]
    b = [r["current_f_b"] for r in rows]
    c = [r["current_f_c"] for r in rows]
    d = [r["current_f_d"] for r in rows]

    # Frame and ROOT-like ticks.
    draw.rectangle(plot, outline=(0, 0, 0), width=3)
    for xt in range(10, 36):
        x, _ = map_xy(float(xt), 0.0, plot)
        major = xt % 5 == 0
        tick = 52 if major else 28
        width = 3 if major else 2
        draw.line([(x, plot[3]), (x, plot[3] - tick)], fill=(0, 0, 0), width=width)
        draw.line([(x, plot[1]), (x, plot[1] + tick)], fill=(0, 0, 0), width=width)
        if major:
            draw.text((x, plot[3] + 36), f"{xt}", font=F_PLOT_LARGE, fill=(0, 0, 0), anchor="ma")

    for i in range(0, 27):
        yt = i * 0.05
        if yt > 1.3:
            continue
        _, y = map_xy(10.0, yt, plot)
        major = abs((yt * 10) % 2) < 1e-6
        tick = 48 if major else 24
        width = 3 if major else 2
        draw.line([(plot[0], y), (plot[0] + tick, y)], fill=(0, 0, 0), width=width)
        draw.line([(plot[2], y), (plot[2] - tick, y)], fill=(0, 0, 0), width=width)
        if major and yt > 0:
            label = "1" if abs(yt - 1.0) < 1e-6 else f"{yt:.1f}"
            draw.text((plot[0] - 18, y), label, font=F_PLOT_LARGE, fill=(0, 0, 0), anchor="rm")

    # Physics text and legend.
    draw.text((310, 125), "sPHENIX", font=F_PLOT_SMALL_BOLD_ITALIC, fill=(0, 0, 0), anchor="la")
    draw.text((535, 125), "Internal", font=F_PLOT_SMALL, fill=(0, 0, 0), anchor="la")
    draw.text((310, 210), "p+p sqrt(s) = 200 GeV", font=F_PLOT_SMALL, fill=(0, 0, 0), anchor="la")
    draw.text((310, 290), "|eta^gamma| < 0.7", font=F_PLOT_SMALL, fill=(0, 0, 0), anchor="la")
    draw.text((310, 370), "RecoilJets current", font=F_PLOT_SMALL, fill=(0, 0, 0), anchor="la")

    leg_x, leg_y = 760, 145
    labels = [
        ((0, 0, 0), "B/A tight noniso"),
        ((255, 55, 25), "C/A nontight iso"),
        ((40, 80, 255), "D/A nontight noniso"),
    ]
    for idx, (color, label) in enumerate(labels):
        y = leg_y + idx * 118
        draw.line([(leg_x, y), (leg_x + 72, y)], fill=color, width=4)
        draw.text((leg_x + 95, y - 45), label, font=F_PLOT_MED, fill=(0, 0, 0), anchor="la")

    line_with_width(draw, step_points(edges, b, plot), (0, 0, 0), 4)
    line_with_width(draw, step_points(edges, c, plot), (255, 55, 25), 4)
    line_with_width(draw, step_points(edges, d, plot), (40, 80, 255), 4)

    draw_rotated_text(img, "Signal leakage", F_PLOT_LARGE, (90, 710), angle=90)
    draw.text((985, 1515), "E_T^gamma,rec [GeV]", font=F_PLOT_LARGE, fill=(0, 0, 0), anchor="la")
    img.convert("RGB").save(RHS_PNG)


def fit_inside(
    img: Image.Image,
    box: tuple[int, int, int, int],
    *,
    bg: tuple[int, int, int] = (255, 255, 255),
) -> tuple[Image.Image, tuple[int, int]]:
    x0, y0, x1, y1 = box
    max_w, max_h = x1 - x0, y1 - y0
    scale = min(max_w / img.width, max_h / img.height)
    new_size = (round(img.width * scale), round(img.height * scale))
    resized = img.convert("RGB").resize(new_size, Image.Resampling.LANCZOS)
    x = x0 + (max_w - new_size[0]) // 2
    y = y0 + (max_h - new_size[1]) // 2
    canvas = Image.new("RGB", (max_w, max_h), bg)
    canvas.paste(resized, ((max_w - new_size[0]) // 2, (max_h - new_size[1]) // 2))
    return canvas, (x, y)


def build_slide() -> None:
    rows = read_rows()
    draw_current_ppg12_style(rows)
    if not PPG12_FIG28_PNG.exists():
        raise FileNotFoundError(f"Missing rendered PPG12 Fig. 28 source: {PPG12_FIG28_PNG}")

    lhs = Image.open(PPG12_FIG28_PNG).convert("RGB")
    rhs = Image.open(RHS_PNG).convert("RGB")
    slide = Image.new("RGB", (2560, 1440), (255, 255, 255))
    draw = ImageDraw.Draw(slide)

    title_color = (28, 38, 54)
    bullet_blue = (24, 103, 179)
    border_blue = (169, 200, 238)
    left_tab = (82, 91, 109)
    right_tab = (32, 120, 188)

    draw.text((82, 48), "Leakage comparison to PPG12 and this analysis", font=F_TITLE, fill=title_color, anchor="la")
    bullet_line = (
        "Check whether pp leakage from this analysis is consistent with the PPG12 working point; "
        "small discrepancies remain, but the pattern is overall coherent."
    )
    bullet_x = 82
    bullet_center_y = 184
    tri = [
        (bullet_x, bullet_center_y - 15),
        (bullet_x, bullet_center_y + 15),
        (bullet_x + 28, bullet_center_y),
    ]
    draw.polygon(tri, fill=bullet_blue)
    draw.text((bullet_x + 52, bullet_center_y), bullet_line, font=F_NOTE, fill=title_color, anchor="lm")

    left_card = (82, 262, 1248, 1312)
    right_card = (1360, 262, 2526, 1312)
    for card in [left_card, right_card]:
        draw.rounded_rectangle(card, radius=12, fill=(255, 255, 255), outline=border_blue, width=3)

    def draw_tab(card: tuple[int, int, int, int], text: str, color: tuple[int, int, int], width: int) -> None:
        x0, y0, _, _ = card
        tab = (x0 + 18, y0 + 14, x0 + 18 + width, y0 + 66)
        draw.rounded_rectangle(tab, radius=8, fill=color, outline=color, width=1)
        # Square off the lower corners so the tab reads like the deck reference.
        draw.rectangle((tab[0], tab[3] - 15, tab[2], tab[3]), fill=color, outline=color)
        draw.text((tab[0] + 16, tab[1] + 10), text, font=F_PANEL, fill=(255, 255, 255), anchor="la")

    draw_tab(left_card, "PPG12 IAN Fig. 28 source", left_tab, 500)
    draw_tab(right_card, "This analysis base pp output", right_tab, 560)

    left_plot = (135, 350, 1200, 1258)
    right_plot = (1413, 350, 2478, 1258)
    lhs_canvas, _ = fit_inside(lhs, left_plot)
    rhs_canvas, _ = fit_inside(rhs, right_plot)
    slide.paste(lhs_canvas, (left_plot[0], left_plot[1]))
    slide.paste(rhs_canvas, (right_plot[0], right_plot[1]))
    draw.text((2516, 1378), "3", font=font(TIMES_BOLD, 28), fill=(0, 0, 0), anchor="ra")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    slide.save(SLIDE_PNG)

    manifest = {
        "schema": "the76_fig28_side_by_side_v1",
        "slide_png": str(SLIDE_PNG),
        "rhs_png": str(RHS_PNG),
        "ppg12_fig28_png": str(PPG12_FIG28_PNG),
        "csv": str(CSV_PATH),
        "ppg12_fig28_sdcc_path": "/sphenix/user/shuhangli/ppg12/PPG12-analysis-note/Figures/analysis/leakage_fraction_et_sim.pdf",
        "later_nonmatching_sdcc_path": "/sphenix/user/shuhangli/ppg12/plotting/figures/leakage_fraction_et_sim_bdt_nom.pdf",
        "sha256": {
            "slide_png": sha256(SLIDE_PNG),
            "rhs_png": sha256(RHS_PNG),
            "ppg12_fig28_png": sha256(PPG12_FIG28_PNG),
            "csv": sha256(CSV_PATH),
        },
        "current_values": rows,
        "note": "The actual analysis-note Fig. 28 source differs from the later regenerated bdt_nom leakage artifact at high ET.",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"slide_png={SLIDE_PNG}")
    print(f"rhs_png={RHS_PNG}")
    print(f"manifest={MANIFEST}")


if __name__ == "__main__":
    build_slide()
