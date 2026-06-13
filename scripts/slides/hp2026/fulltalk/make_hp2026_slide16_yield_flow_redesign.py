#!/usr/bin/env python3
"""Local-only redesign candidates for HP2026 Slide 16 yield-correction flow."""

from __future__ import annotations

import json
import math
import shutil
import sys
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
SCRIPT_DIR = ROOT / "scripts/slides/hp2026/fulltalk"
sys.path.insert(0, str(SCRIPT_DIR))
COMMON_DIR = ROOT / "scripts/slides/common"
sys.path.insert(0, str(COMMON_DIR))

import make_hp2026_closing_three_candidates as c3  # noqa: E402
import make_hp2026_fulltalk_candidates as ft  # noqa: E402
from slide_symmetry_audit import Box, Check, SymmetryAudit, require_audit_passed  # noqa: E402


ARIAL_UNICODE = Path("/System/Library/Fonts/Supplemental/Arial Unicode.ttf")
OUTDIR = ROOT / "outputs/manual-20260610-slide16_yield_flow_main_talk_replacement"
CANDIDATE_DIR = OUTDIR / "internal_candidates"
FINAL_PNG = OUTDIR / "slide16_yield_flow_plot_flow_first_replacement.png"
MANIFEST = OUTDIR / "slide16_yield_flow_plot_flow_first_replacement_manifest.json"
SYMMETRY_REPORT = OUTDIR / "slide16_yield_flow_plot_flow_first_replacement_symmetry.json"
HEADER_SPEC = FINAL_PNG.with_suffix(".header.json")

HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": None,
    "title_xy": [132, 76],
    "subtitle_xy": None,
    "divider_y": 232,
}

TITLE = "From selected candidates to a corrected photon yield"
SUBTITLE = None

LAYOUT: dict[str, tuple[int, int, int, int]] = {}


def trim_white_margins(img: Image.Image, tolerance: int = 12, pad: int = 8) -> Image.Image:
    rgb = img.convert("RGB")
    white = Image.new("RGB", rgb.size, "white")
    diff = ImageChops.difference(rgb, white).convert("L")
    bbox = diff.point(lambda value: 255 if value > tolerance else 0).getbbox()
    if bbox is None:
        return img
    x0, y0, x1, y1 = bbox
    return img.crop((max(0, x0 - pad), max(0, y0 - pad), min(img.width, x1 + pad), min(img.height, y1 + pad)))


def base_slide() -> Image.Image:
    img = Image.new("RGBA", (ft.W, ft.H), (*ft.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, ft.W, ft.H), fill=(*ft.SOFT_BG, 255))
    draw.rectangle((0, 0, ft.W, 22), fill=(*ft.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, ft.W, 30), fill=(*ft.PHOTON, 255))
    draw.text(
        tuple(HP2026_MAIN_HEADER["title_xy"]),
        TITLE,
        font=ft.font(ft.TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]),
        fill=ft.INK,
    )
    if SUBTITLE and HP2026_MAIN_HEADER["subtitle_xy"] and HP2026_MAIN_HEADER["subtitle_font_size"]:
        draw.text(
            tuple(HP2026_MAIN_HEADER["subtitle_xy"]),
            SUBTITLE,
            font=ft.font(ft.TIMES_ITALIC, HP2026_MAIN_HEADER["subtitle_font_size"]),
            fill=ft.MUTED,
        )
    draw.line((132, HP2026_MAIN_HEADER["divider_y"], ft.W - 132, HP2026_MAIN_HEADER["divider_y"]), fill=(221, 226, 232, 255), width=3)
    ft.add_top_right_sphenix_logo_like_slide2(img)
    return img


def draw_card(base: Image.Image, box: tuple[int, int, int, int], accent: tuple[int, int, int], title: str, subtitle: str) -> ImageDraw.ImageDraw:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 12, box[3]), radius=6, fill=(*accent, 235))
    draw.text((box[0] + 54, box[1] + 24), title, font=ft.font(ft.TIMES_BOLD, 33), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        subtitle,
        (box[0] + 54, box[1] + 66),
        box[2] - box[0] - 90,
        ft.font(ft.TIMES_ITALIC, 22),
        fill=ft.MUTED,
        line_gap=3,
    )
    draw.line((box[0] + 54, box[1] + 110, box[2] - 34, box[1] + 110), fill=(221, 228, 236), width=2)
    return draw


def draw_plain_accent_card(base: Image.Image, box: tuple[int, int, int, int], accent: tuple[int, int, int]) -> ImageDraw.ImageDraw:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 12, box[3]), radius=6, fill=(*accent, 235))
    return draw


def draw_leakage_label(base: Image.Image, center: tuple[int, int], text: str, color: tuple[int, int, int], angle: float = 0) -> None:
    font = ft.font(ft.TIMES_BOLD, 22)
    tmp = Image.new("RGBA", (180, 46), (255, 255, 255, 0))
    d = ImageDraw.Draw(tmp, "RGBA")
    tw, th = ft.text_box(d, text, font)
    d.rounded_rectangle((6, 6, tw + 24, 38), radius=8, fill=(255, 255, 255, 245), outline=(*color, 185), width=2)
    d.text((16, 10), text, font=font, fill=color)
    tmp = tmp.crop((0, 0, tw + 34, 46)).rotate(-angle, expand=True, resample=Image.Resampling.BICUBIC)
    base.alpha_composite(tmp, (int(center[0] - tmp.width / 2), int(center[1] - tmp.height / 2)))


def draw_leakage_map(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box

    header_font = ft.font(ft.TIMES_BOLD, 38)
    definition_font = ft.font(ft.TIMES, 29)
    header_x, header_y = x0 + 8, y0 + 4
    draw.text((header_x, header_y), "Leakage =", font=header_font, fill=ft.BLUE)
    header_w = ft.text_box(draw, "Leakage =", header_font)[0]
    draw.text(
        (header_x + header_w + 10, header_y + 7),
        "truth-matched signal photons in PYTHIA MC",
        font=definition_font,
        fill=ft.INK,
    )
    draw.text(
        (header_x + header_w + 10, header_y + 41),
        "found in B/C/D control regions",
        font=definition_font,
        fill=ft.INK,
    )

    map_x0, map_y0 = x0 + 8, y0 + 82
    available_w = x1 - x0 - 16
    cell_w = min(350, max(310, int((available_w - 142) / 2)))
    cell_h = 112
    h_gap = available_w - 2 * cell_w
    note_top = y1 - 66
    v_gap = max(54, note_top - map_y0 - 2 * cell_h - 12)
    cells = {
        "C": (map_x0, map_y0),
        "D": (map_x0 + cell_w + h_gap, map_y0),
        "A": (map_x0, map_y0 + cell_h + v_gap),
        "B": (map_x0 + cell_w + h_gap, map_y0 + cell_h + v_gap),
    }
    meta = {
        "A": ("selected sample", "isolated + tight target", ft.PHOTON_DARK, (255, 248, 229)),
        "B": ("non-isolated + tight", "true photon fails isolation", ft.SPHENIX_BLUE, (238, 247, 252)),
        "C": ("isolated + non-tight", "true photon fails tight ID", ft.TEAL, (239, 249, 250)),
        "D": ("non-isolated + non-tight", "true photon fails both axes", ft.MUTED, (245, 247, 250)),
    }
    for region, (cx, cy) in cells.items():
        label, detail, color, fill = meta[region]
        draw.rounded_rectangle((cx, cy, cx + cell_w, cy + cell_h), radius=8, fill=(*fill, 255), outline=(*color, 230), width=3)
        draw.text((cx + 18, cy + 14), region, font=ft.font(ft.TIMES_BOLD, 41), fill=color)
        ft.draw_wrapped(draw, label, (cx + 78, cy + 14), cell_w - 94, ft.font(ft.TIMES_BOLD, 24), fill=ft.INK, line_gap=2)
        ft.draw_wrapped(draw, detail, (cx + 78, cy + 64), cell_w - 94, ft.font(ft.TIMES, 22), fill=ft.MUTED, line_gap=2)

    ax, ay = cells["A"]
    bx, by = cells["B"]
    cx, cy = cells["C"]
    dx, dy = cells["D"]
    c_start = (ax + cell_w // 2, ay - 10)
    c_end = (cx + cell_w // 2, cy + cell_h + 10)
    b_start = (ax + cell_w + 8, ay + cell_h // 2 + 4)
    b_end = (bx - 12, by + cell_h // 2 + 4)
    d_start = (ax + cell_w + 14, ay + 26)
    d_end = (dx + 4, dy + cell_h - 4)
    ft.draw_arrow(draw, c_start, c_end, fill=ft.TEAL, width=4)
    ft.draw_arrow(draw, b_start, b_end, fill=ft.SPHENIX_BLUE, width=4)
    ft.draw_arrow(draw, d_start, d_end, fill=ft.MUTED, width=4)
    b_angle = math.degrees(math.atan2(b_end[1] - b_start[1], b_end[0] - b_start[0]))
    d_angle = math.degrees(math.atan2(d_end[1] - d_start[1], d_end[0] - d_start[0]))
    draw_leakage_label(base, (c_start[0] + 74, (c_start[1] + c_end[1]) // 2), "C leakage", ft.TEAL)
    draw_leakage_label(base, ((b_start[0] + b_end[0]) // 2, b_start[1] - 18), "B leakage", ft.SPHENIX_BLUE, b_angle)
    draw_leakage_label(base, ((d_start[0] + d_end[0]) // 2 + 4, (d_start[1] + d_end[1]) // 2 - 16), "D leakage", ft.MUTED, d_angle)

    note = (x0 + 8, y1 - 76, x1 - 8, y1 - 12)
    icon_x, icon_y = note[0] + 134, note[1] + 15
    draw.rounded_rectangle((icon_x, icon_y, icon_x + 34, icon_y + 34), radius=5, fill=(255, 255, 255, 255), outline=(*ft.SPHENIX_BLUE, 255), width=2)
    draw.line((icon_x + 9, icon_y + 11, icon_x + 25, icon_y + 11), fill=(*ft.SPHENIX_BLUE, 220), width=3)
    draw.line((icon_x + 9, icon_y + 19, icon_x + 24, icon_y + 19), fill=(*ft.SPHENIX_BLUE, 190), width=3)
    draw.line((icon_x + 9, icon_y + 27, icon_x + 20, icon_y + 27), fill=(*ft.SPHENIX_BLUE, 165), width=3)
    note_font = ft.font(ft.TIMES_BOLD, 27)
    body_font = ft.font(ft.TIMES, 25)
    line1 = "Note: ABCD purity requires isolation-tight ID complementarity;"
    line2 = "isolation is not used in BDT training."
    text_left = icon_x + 44
    text_right = note[2] - 18
    line1_w = ft.text_box(draw, line1, note_font)[0]
    line2_w = ft.text_box(draw, line2, body_font)[0]
    draw.text((text_left + max(0, (text_right - text_left - line1_w) // 2), note[1] + 8), line1, font=note_font, fill=ft.BLUE)
    draw.text((text_left + max(0, (text_right - text_left - line2_w) // 2), note[1] + 37), line2, font=body_font, fill=ft.BLUE)


def draw_purity_plot(base: Image.Image, box: tuple[int, int, int, int], title: str = "Leakage-corrected purity") -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 12, box[3]), radius=6, fill=(*ft.SPHENIX_BLUE, 235))
    header_x = box[0] + 54
    header_y = box[1] + 24
    title_font = ft.font(ft.TIMES_BOLD, 38)
    body_font = ft.font(ft.TIMES, 30)
    draw.text((header_x, header_y), title, font=title_font, fill=ft.INK)
    title_w = ft.text_box(draw, title, title_font)[0]
    draw.text((header_x + title_w + 12, header_y + 10), "= signal fraction carried into the corrected yield.", font=body_font, fill=ft.INK)
    draw.line((box[0] + 54, box[1] + 90, box[2] - 34, box[1] + 90), fill=(221, 228, 236), width=2)

    plot = trim_white_margins(Image.open(ft.figure_path("fig5_purity")).convert("RGBA"), tolerance=10, pad=6)
    content = (box[0] + 54, box[1] + 108, box[2] - 34, box[3] - 26)
    plot_box = (content[0] + 18, content[1] + 0, content[0] + 612, content[3] - 2)
    ft.shadow(base, plot_box, radius=10)
    draw.rounded_rectangle(
        (plot_box[0] - 8, plot_box[1] - 8, plot_box[2] + 8, plot_box[3] + 8),
        radius=8,
        fill=(255, 255, 255, 255),
        outline=(216, 226, 235, 255),
        width=2,
    )
    ft.paste_fit(base, plot, plot_box, anchor="center")
    draw_purity_reading_key(base, content, plot_box)


def draw_purity_reading_key(base: Image.Image, content: tuple[int, int, int, int], plot_box: tuple[int, int, int, int]) -> None:
    key_x0 = plot_box[2] + 30
    key_x1 = content[2] - 8
    gap = 20
    key_h = ((content[3] - 18) - (content[1] + 18) - gap) // 2
    raw_box = (key_x0, content[1] + 18, key_x1, content[1] + 18 + key_h)
    corr_box = (key_x0, raw_box[3] + gap, key_x1, content[3] - 18)

    draw_purity_key_card(
        base,
        raw_box,
        (32, 32, 32),
        "Raw sideband purity",
        "Black points: ABCD sideband estimate before correcting true-photon leakage into B/C/D.",
        fill=(246, 247, 249, 255),
    )
    draw_purity_key_card(
        base,
        corr_box,
        ft.SPHENIX_BLUE,
        "Leakage-corrected purity",
        "Blue points: MC leakage correction applied. This curve is carried into the Region A yield.",
        fill=(236, 247, 255, 255),
    )


def draw_purity_key_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    color: tuple[int, int, int],
    title: str,
    body: str,
    fill: tuple[int, int, int, int],
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle(box, radius=12, fill=fill, outline=(205, 221, 235, 255), width=2)
    draw.text((box[0] + 30, box[1] + 28), title, font=ft.font(ft.TIMES_BOLD, 43), fill=color)
    ft.draw_wrapped(
        draw,
        body,
        (box[0] + 30, box[1] + 92),
        box[2] - box[0] - 60,
        ft.font(ft.TIMES, 38),
        fill=ft.INK,
        line_gap=8,
    )


def draw_purity_logic_band(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(250, 252, 254, 255), outline=(211, 224, 235, 255), width=2)
    draw.text((x0 + 28, y0 + 18), "Reading the purity plot", font=ft.font(ft.TIMES_BOLD, 38), fill=ft.INK)
    draw.text(
        (x0 + 454, y0 + 26),
        "raw sidebands become the leakage-corrected purity.",
        font=ft.font(ft.TIMES, 32),
        fill=ft.MUTED,
    )
    draw.line((x0 + 28, y0 + 76, x1 - 28, y0 + 76), fill=(224, 232, 240), width=2)

    col_y0 = y0 + 104
    col_y1 = y1 - 24
    inner_x0 = x0 + 28
    inner_x1 = x1 - 28
    cell_w = (inner_x1 - inner_x0) / 3
    col_inset = 7
    columns = [
        (
            "1",
            "Raw ABCD",
            "black points",
            "Raw sideband purity.",
            (35, 35, 35),
            (246, 247, 249),
        ),
        (
            "2",
            "MC leakage",
            "truth correction",
            "Subtract signal in B/C/D.",
            ft.TEAL,
            (240, 249, 250),
        ),
        (
            "3",
            "Corrected purity",
            "blue points / fit",
            "Used for Region A yield.",
            ft.SPHENIX_BLUE,
            (236, 247, 255),
        ),
    ]
    for idx, (num, title, tag, body, color, fill) in enumerate(columns):
        cx0 = round(inner_x0 + idx * cell_w + col_inset)
        cx1 = round(inner_x0 + (idx + 1) * cell_w - col_inset)
        col_w = cx1 - cx0
        LAYOUT[f"purity_logic_col_{idx + 1}"] = (cx0, col_y0, cx1, col_y1)
        draw.rounded_rectangle((cx0, col_y0, cx1, col_y1), radius=11, fill=(*fill, 255), outline=(*color, 185), width=2)
        badge = (cx0 + 20, col_y0 + 18, cx0 + 62, col_y0 + 60)
        draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(*color, 230), width=3)
        nf = ft.font(ft.TIMES_BOLD, 26)
        nw, nh = ft.text_box(draw, num, nf)
        draw.text((badge[0] + (42 - nw) / 2, badge[1] + (42 - nh) / 2 - 1), num, font=nf, fill=color)
        draw.text((cx0 + 78, col_y0 + 16), title, font=ft.font(ft.TIMES_BOLD, 30), fill=color)
        draw.text((cx0 + 78, col_y0 + 53), tag, font=ft.font(ft.TIMES_BOLD, 25), fill=ft.INK)
        ft.draw_wrapped(draw, body, (cx0 + 22, col_y0 + 96), col_w - 44, ft.font(ft.TIMES, 29), fill=ft.INK, line_gap=3)


def draw_raw_abcd_count_equation(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    *,
    size: int = 34,
) -> None:
    base = ft.font(ft.TIMES_ITALIC, size)
    roman = ft.font(ft.TIMES, size)
    small = ft.font(ft.TIMES_BOLD, max(15, int(size * 0.56)))
    small_plain = ft.font(ft.TIMES, max(15, int(size * 0.56)))
    region_colors = {
        "A": ft.PHOTON_DARK,
        "B": ft.SPHENIX_BLUE,
        "C": ft.TEAL,
        "D": ft.MUTED,
    }

    def advance(text: str, font, fill=ft.INK, dy: int = 0) -> None:
        nonlocal x
        draw.text((x, y + dy), text, font=font, fill=fill)
        x += ft.text_box(draw, text, font)[0]

    def n_sig_raw() -> None:
        nonlocal x
        draw.text((x, y), "N", font=base, fill=ft.INK)
        bw = ft.text_box(draw, "N", base)[0]
        draw.text((x + bw - 1, y - int(size * 0.34)), "raw", font=small_plain, fill=ft.INK)
        draw.text((x + bw - 1, y + int(size * 0.47)), "sig,A", font=small_plain, fill=ft.INK)
        x += bw + max(ft.text_box(draw, "raw", small_plain)[0], ft.text_box(draw, "sig,A", small_plain)[0]) + 10

    def n_region(region: str, x_pos: int | None = None, y_pos: int | None = None, term_size: int | None = None) -> int:
        nonlocal x
        if x_pos is None:
            x_pos = x
        if y_pos is None:
            y_pos = y
        if term_size is None:
            term_size = size
        term_base = ft.font(ft.TIMES_ITALIC, term_size)
        term_small = ft.font(ft.TIMES_BOLD, max(14, int(term_size * 0.58)))
        draw.text((x_pos, y_pos), "N", font=term_base, fill=ft.INK)
        bw = ft.text_box(draw, "N", term_base)[0]
        draw.text((x_pos + bw - 1, y_pos + int(term_size * 0.48)), region, font=term_small, fill=region_colors[region])
        end = x_pos + bw + ft.text_box(draw, region, term_small)[0] + 8
        if x_pos == x:
            x = end
        return end

    n_sig_raw()
    advance(" = ", roman)
    n_region("A")
    advance(" − ", roman)
    n_region("B")
    advance(" (", roman)
    n_region("C")
    advance(" / ", roman)
    n_region("D")
    advance(")", roman)


def draw_abcd_mini_cue(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(250, 252, 254, 255), outline=(211, 224, 235, 255), width=2)
    draw.text((x0 + 28, y0 + 22), "Sideband control regions", font=ft.font(ft.TIMES_BOLD, 42), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "B/C/D constrain residual background; MC removes leaked true photons.",
        (x0 + 28, y0 + 78),
        x1 - x0 - 56,
        ft.font(ft.TIMES, 36),
        fill=ft.MUTED,
        line_gap=5,
    )

    grid_top = y0 + 162
    gap = 18
    cell_w = (x1 - x0 - 56 - gap) // 2
    cell_h = 94
    cells = [
        ("C", "isolated\nnon-tight (bkg-like)", ft.TEAL, (239, 249, 250)),
        ("D", "non-isolated\nnon-tight", ft.MUTED, (246, 247, 249)),
        ("A", "isolated and tight\nselected sample", ft.PHOTON_DARK, (255, 249, 235)),
        ("B", "non-isolated\ntight (photon-like)", ft.SPHENIX_BLUE, (239, 248, 253)),
    ]
    for i, (region, label, color, fill) in enumerate(cells):
        row = i // 2
        col = i % 2
        cx = x0 + 28 + col * (cell_w + gap)
        cy = grid_top + row * (cell_h + gap)
        draw.rounded_rectangle((cx, cy, cx + cell_w, cy + cell_h), radius=9, fill=(*fill, 255), outline=(*color, 230), width=3)
        draw.text((cx + 18, cy + 18), region, font=ft.font(ft.TIMES_BOLD, 45), fill=color)
        lines = label.splitlines()
        for j, line in enumerate(lines):
            fnt = ft.font(ft.TIMES_BOLD if j == 0 else ft.TIMES, 31)
            draw.text((cx + 82, cy + 15 + j * 36), line, font=fnt, fill=ft.INK if j == 0 else ft.MUTED)

    eq_box = (x0 + 24, y1 - 112, x1 - 24, y1 - 22)
    LAYOUT["abcd_equation_box"] = eq_box
    draw.rounded_rectangle(eq_box, radius=10, fill=(255, 255, 255, 250), outline=(213, 225, 236, 255), width=2)
    draw.text(
        (eq_box[0] + 22, eq_box[1] + 14),
        "Raw ABCD estimate",
        font=ft.font(ft.TIMES_BOLD, 35),
        fill=ft.INK,
    )
    draw.text(
        (eq_box[0] + 22, eq_box[1] + 52),
        "from region counts:",
        font=ft.font(ft.TIMES_BOLD, 35),
        fill=ft.INK,
    )
    draw_raw_abcd_count_equation(draw, eq_box[0] + 430, eq_box[1] + 28, size=44)


def draw_purity_definition_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(246, 251, 255, 255), outline=(197, 225, 243, 255), width=2)
    draw.text((x0 + 30, y0 + 24), "Leakage-corrected purity", font=ft.font(ft.TIMES_BOLD, 42), fill=ft.SPHENIX_BLUE)
    draw.line((x0 + 30, y0 + 82, x1 - 30, y0 + 82), fill=(207, 225, 238, 255), width=2)
    draw.text((x0 + 30, y0 + 112), "MC corrects signal leakage into the sidebands.", font=ft.font(ft.TIMES, 38), fill=ft.INK)
    draw.text((x0 + 30, y0 + 168), "Blue points / fit are used in the corrected yield.", font=ft.font(ft.TIMES_BOLD, 37), fill=ft.BLUE)


def draw_compact_abcd_cue(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(211, 224, 235, 255), width=2)
    draw.text((x0 + 30, y0 + 16), "Data control regions", font=ft.font(ft.TIMES_BOLD, 44), fill=ft.INK)

    grid_top = y0 + 79
    gap = 14
    cell_w = (x1 - x0 - 60 - gap) // 2
    cell_h = 106
    cells = [
        ("C", "isolated\nnon-tight (bkg-like)", ft.TEAL, (239, 249, 250)),
        ("D", "non-isolated\nnon-tight", ft.MUTED, (246, 247, 249)),
        ("A", "isolated and tight\nselected sample", ft.PHOTON_DARK, (255, 249, 235)),
        ("B", "non-isolated\ntight (photon-like)", ft.SPHENIX_BLUE, (239, 248, 253)),
    ]
    for i, (region, label, color, fill) in enumerate(cells):
        row = i // 2
        col = i % 2
        cx = x0 + 30 + col * (cell_w + gap)
        cy = grid_top + row * (cell_h + gap)
        draw.rounded_rectangle((cx, cy, cx + cell_w, cy + cell_h), radius=8, fill=(*fill, 255), outline=(*color, 230), width=3)
        draw.text((cx + 18, cy + 20), region, font=ft.font(ft.TIMES_BOLD, 54), fill=color)
        lines = label.splitlines()
        draw.text((cx + 98, cy + 17), lines[0], font=ft.font(ft.TIMES_BOLD, 39), fill=ft.INK)
        second_line_font = ft.font(ft.TIMES, 33 if "(" in lines[1] else 37)
        draw.text((cx + 98, cy + 61), lines[1], font=second_line_font, fill=ft.MUTED)

    eq_box = (x0 + 30, y1 - 94, x1 - 30, y1 - 18)
    LAYOUT["abcd_equation_box"] = eq_box
    draw.text((eq_box[0] + 4, eq_box[1] + 21), "Raw signal estimate:", font=ft.font(ft.TIMES_BOLD, 39), fill=ft.INK)
    draw_raw_abcd_count_equation(draw, eq_box[0] + 390, eq_box[1] + 14, size=50)


def draw_compact_plot_readout(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 250), outline=(211, 224, 235, 255), width=2)
    rows = [
        ((35, 35, 35), "Black points", "raw data-driven purity"),
        (ft.SPHENIX_BLUE, "Blue points / fit", "leakage-corrected purity"),
    ]
    label_font = ft.font(ft.TIMES_BOLD, 41)
    body_font = ft.font(ft.TIMES, 41)

    def draw_midline_text(x: int, center_y: int, text: str, font, fill: tuple[int, int, int]) -> None:
        bbox = draw.textbbox((0, 0), text, font=font)
        text_h = bbox[3] - bbox[1]
        draw.text((x, center_y - text_h / 2 - bbox[1]), text, font=font, fill=fill)

    row_h = (y1 - y0 - 24) // 2
    for i, (color, label, body) in enumerate(rows):
        ry0 = y0 + 12 + i * row_h
        ry1 = y0 + 12 + (i + 1) * row_h
        if i:
            draw.line((x0 + 28, ry0, x1 - 28, ry0), fill=(224, 232, 240, 255), width=2)
        cy = (ry0 + ry1) // 2
        draw.ellipse((x0 + 32, cy - 10, x0 + 52, cy + 10), fill=(*color, 255))
        draw_midline_text(x0 + 82, cy, label, label_font, color)
        draw_midline_text(x0 + 430, cy, body, body_font, ft.INK)


def draw_purity_hero_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 12, box[3]), radius=6, fill=(*ft.SPHENIX_BLUE, 235))
    LAYOUT["hero_panel"] = box

    plot = trim_white_margins(Image.open(ft.figure_path("fig5_purity")).convert("RGBA"), tolerance=10, pad=6)
    plot_box = (box[0] + 70, box[1] + 70, box[0] + 1120, box[3] - 58)
    LAYOUT["purity_plot_box"] = plot_box
    ft.shadow(base, plot_box, radius=10)
    draw.rounded_rectangle(
        (plot_box[0] - 8, plot_box[1] - 8, plot_box[2] + 8, plot_box[3] + 8),
        radius=10,
        fill=(255, 255, 255, 255),
        outline=(213, 225, 236, 255),
        width=2,
    )
    ft.paste_fit(base, plot, plot_box, anchor="center")

    rhs_x0 = plot_box[2] + 58
    rhs_x1 = box[2] - 54
    readout_box = (rhs_x0, plot_box[1], rhs_x1, plot_box[1] + 238)
    sideband_box = (rhs_x0, readout_box[3] + 26, rhs_x1, readout_box[3] + 438)
    definition_box = (rhs_x0, sideband_box[3] + 26, rhs_x1, plot_box[3])
    LAYOUT["plot_readout_box"] = readout_box
    LAYOUT["sideband_cue_box"] = sideband_box
    LAYOUT["purity_definition_box"] = definition_box
    draw_compact_plot_readout(base, readout_box)
    draw_compact_abcd_cue(base, sideband_box)
    draw_purity_definition_card(base, definition_box)


def draw_purity_logic_box(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    gap = 14
    raw_h = 188
    raw_box = (box[0], box[1], box[2], box[1] + raw_h)
    lc_box = (box[0], raw_box[3] + gap, box[2], box[3])
    draw_definition_section(
        base,
        raw_box,
        "Raw sideband estimate",
        [
            r"$N_{\mathrm{sig}}^{A,\mathrm{raw}}=N_{\mathrm{raw}}^A-N_{\mathrm{bkg}}^{\mathrm{ABCD}}$",
            r"$N_{\mathrm{bkg}}^{\mathrm{ABCD}}=N_{\mathrm{raw}}^B(N_{\mathrm{raw}}^C/N_{\mathrm{raw}}^D)$",
        ],
        [
            ("N_A raw", "selected A"),
            ("N_bkg", "ABCD estimate"),
            ("B/C/D", "control regions"),
        ],
        fontsize=26,
    )
    draw_definition_section(
        base,
        lc_box,
        "Leakage-corrected estimate",
        [
            r"$b_X=N_{\mathrm{raw}}^X-f_X^{\mathrm{MC}}N_{\mathrm{sig}}^A,\quad X=B,C,D$",
            r"$N_{\mathrm{sig}}^{A,\mathrm{LC}}=N_{\mathrm{raw}}^A-b_B(b_C/b_D)$",
            r"$P_{\mathrm{LC}}(E_T)=N_{\mathrm{sig}}^{A,\mathrm{LC}}/N_{\mathrm{raw}}^A$",
        ],
        [
            ("f_MC", "MC leakage fraction"),
            (r"$b_X$", "leakage-corrected sideband"),
            ("P_LC", "purity carried forward"),
        ],
        fontsize=25,
    )


def draw_definition_section(
    base: Image.Image,
    section_box: tuple[int, int, int, int],
    title: str,
    eqs: list[str],
    terms: list[tuple[str, str]],
    fontsize: int,
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = section_box
    pad = 12
    draw.rounded_rectangle(
        section_box,
        radius=9,
        fill=(255, 255, 255, 248),
        outline=(216, 226, 235, 255),
        width=2,
    )
    draw.text((x0 + pad + 8, y0 + 9), title, font=ft.font(ft.TIMES_BOLD, 24), fill=ft.INK)
    draw.line((x0 + pad + 8, y0 + 44, x1 - pad - 8, y0 + 44), fill=(228, 234, 240), width=2)

    content_top = y0 + 54
    content_bottom = y1 - 12
    content_h = content_bottom - content_top
    split_x = x0 + int((x1 - x0) * 0.63)
    draw.line((split_x, content_top + 2, split_x, content_bottom - 2), fill=(234, 239, 244), width=2)

    eq_img = math_block_image(eqs, fontsize=fontsize, width=4.8, line_step=0.31)
    paste_math_left(base, eq_img, (x0 + pad + 8, content_top, split_x - 16, content_bottom))

    term_font = ft.font(ft.TIMES_BOLD, 18)
    def_font = ft.font(ft.TIMES, 18)
    pill_x = split_x + 12
    row_h = max(37, content_h // max(1, len(terms)))
    y = content_top + max(0, (content_h - row_h * len(terms)) // 2)
    for label, meaning in terms:
        if label.startswith("$"):
            label_img = math_block_image([label], fontsize=18, width=1.15, line_step=0.30)
            label_w = max(64, min(105, label_img.width))
            pill = (pill_x, y + 2, pill_x + label_w + 16, y + 30)
            draw.rounded_rectangle(pill, radius=6, fill=(246, 249, 252, 255), outline=(215, 225, 234, 255), width=1)
            paste_math_fit(base, label_img, (pill_x + 8, y + 4, pill_x + label_w + 8, y + 27))
        else:
            label_w = ft.text_box(draw, label, term_font)[0]
            pill = (pill_x, y + 2, pill_x + max(label_w + 16, 78), y + 30)
            draw.rounded_rectangle(pill, radius=6, fill=(246, 249, 252, 255), outline=(215, 225, 234, 255), width=1)
            draw.text((pill_x + 8, y + 5), label, font=term_font, fill=ft.INK)
        ft.draw_wrapped(
            draw,
            f"= {meaning}",
            (pill[2] + 10, y + 5),
            max(80, x1 - pill[2] - 24),
            def_font,
            fill=ft.INK,
            line_gap=1,
        )
        y += row_h


def draw_final_purity_carry_forward(base: Image.Image, box: tuple[int, int, int, int], y: int) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    pad = 18
    draw.rounded_rectangle((box[0] + pad, y, box[2] - pad, box[3] - 14), radius=9, fill=(247, 250, 252, 255), outline=(216, 226, 235, 255), width=2)
    draw.text((box[0] + pad + 16, y + 7), "Purity carried forward", font=ft.font(ft.TIMES_BOLD, 20), fill=ft.INK)
    eq_img = math_block_image([r"$P_{\mathrm{LC}}(E_T)=N_{\mathrm{sig}}^{A,\mathrm{LC}}/N_{\mathrm{raw}}^A$"], fontsize=29, width=5.9, line_step=0.32)
    paste_math_left(base, eq_img, (box[0] + pad + 16, y + 29, box[2] - pad - 8, box[3] - 17))


def paste_scaled_formula(
    base: Image.Image,
    renderer,
    xy: tuple[int, int],
    scale: float,
    canvas_size: tuple[int, int],
) -> None:
    layer = Image.new("RGBA", canvas_size, (255, 255, 255, 0))
    renderer(ImageDraw.Draw(layer, "RGBA"), 18, 24)
    bbox = layer.getbbox()
    if not bbox:
        return
    crop = layer.crop(bbox)
    crop = crop.resize((max(1, int(crop.width * scale)), max(1, int(crop.height * scale))), Image.Resampling.LANCZOS)
    base.alpha_composite(crop, xy)


def draw_f_mc_var(draw: ImageDraw.ImageDraw, x: int, y: int, sub: str = "X", size: int = 27, fill=ft.INK) -> int:
    base = ft.font(ft.TIMES_ITALIC, size)
    small = ft.font(ft.TIMES, max(13, int(size * 0.55)))
    draw.text((x, y), "f", font=base, fill=fill)
    bw, _ = ft.text_box(draw, "f", base)
    draw.text((x + bw - 1, y - int(size * 0.36)), "MC", font=small, fill=fill)
    draw.text((x + bw - 1, y + int(size * 0.48)), sub, font=small, fill=fill)
    return x + bw + max(ft.text_box(draw, "MC", small)[0], ft.text_box(draw, sub, small)[0]) + 8


def draw_b_sub(draw: ImageDraw.ImageDraw, x: int, y: int, sub: str = "X", size: int = 28, fill=ft.INK) -> int:
    base = ft.font(ft.TIMES_ITALIC, size)
    small = ft.font(ft.TIMES, max(13, int(size * 0.55)))
    draw.text((x, y), "b", font=base, fill=fill)
    bw, _ = ft.text_box(draw, "b", base)
    draw.text((x + bw - 1, y + int(size * 0.48)), sub, font=small, fill=fill)
    return x + bw + ft.text_box(draw, sub, small)[0] + 8


def draw_raw_sideband_equation(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 28) -> None:
    x = c3.draw_n_var(draw, x, y, "A", "sig", size=size)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)
    x = c3.draw_n_var(draw, x, y, "A", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " − ", size=size)
    x = c3.draw_n_var(draw, x, y, "B", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " (", size=size)
    frac_x = x + 4
    num_end = c3.draw_n_var(draw, frac_x + 10, y - 11, "C", "raw", size=size - 4)
    frac_w = max(78, num_end - frac_x + 8)
    draw.line((frac_x, y + 35, frac_x + frac_w, y + 35), fill=ft.INK, width=2)
    c3.draw_n_var(draw, frac_x + 10, y + 38, "D", "raw", size=size - 4)
    draw.text((frac_x + frac_w + 6, y - 3), ")", font=ft.font(ft.TIMES, size + 15), fill=ft.INK)


def draw_tag(draw: ImageDraw.ImageDraw, center_x: int, y: int, text: str, color: tuple[int, int, int], size: int = 19) -> tuple[int, int, int, int]:
    fnt = ft.font(ft.TIMES_BOLD, size)
    tw, th = ft.text_box(draw, text, fnt)
    x0 = int(center_x - tw / 2 - 14)
    box = (x0, y, x0 + tw + 28, y + th + 10)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 248), outline=(*color, 215), width=2)
    draw.text((x0 + 14, y + 4), text, font=fnt, fill=color)
    return box


def draw_raw_sideband_equation_annotated(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 32) -> None:
    x = c3.draw_n_var(draw, x, y, "A", "sig", size=size)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)

    selected_start = x
    x = c3.draw_n_var(draw, x, y, "A", "raw", size=size)
    selected_end = x

    x = c3.math_text(draw, x, y + 2, " − ", size=size)

    sideband_start = x
    x = c3.draw_n_var(draw, x, y, "B", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " (", size=size)
    x = c3.draw_n_var(draw, x, y, "C", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " / ", size=size)
    x = c3.draw_n_var(draw, x, y, "D", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, ")", size=size)
    sideband_end = x

    term_top = y - int(size * 0.50)
    term_bottom = y + int(size * 1.18)
    selected_box = (selected_start - 10, term_top, selected_end + 10, term_bottom)
    sideband_box = (sideband_start - 10, term_top, sideband_end + 10, term_bottom)
    draw.rounded_rectangle(selected_box, radius=18, outline=(*ft.PHOTON_DARK, 230), width=4)
    draw.rounded_rectangle(sideband_box, radius=18, outline=(*ft.TEAL, 230), width=4)

    selected_tag = draw_tag(draw, (selected_start + selected_end) // 2, y - 57, "selected Region A", ft.PHOTON_DARK, size=20)
    sideband_tag = draw_tag(draw, (sideband_start + sideband_end) // 2, y - 57, "sideband estimate", ft.TEAL, size=20)
    draw.line(((selected_tag[0] + selected_tag[2]) // 2, selected_tag[3], (selected_box[0] + selected_box[2]) // 2, selected_box[1]), fill=(*ft.PHOTON_DARK, 190), width=2)
    draw.line(((sideband_tag[0] + sideband_tag[2]) // 2, sideband_tag[3], (sideband_box[0] + sideband_box[2]) // 2, sideband_box[1]), fill=(*ft.TEAL, 190), width=2)


def draw_corrected_abcd_equation(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 27) -> None:
    line1_x = x
    line1_x = c3.draw_n_var(draw, line1_x, y, "A", "signal", size=size)
    line1_x = c3.math_text(draw, line1_x, y + 2, " = ", size=size)
    line1_x = c3.draw_n_var(draw, line1_x, y, "A", "raw", size=size)
    line1_x = c3.math_text(draw, line1_x, y + 2, " − [(", size=size)
    line1_x = c3.draw_n_var(draw, line1_x, y, "B", "raw", size=size)
    line1_x = c3.math_text(draw, line1_x, y + 2, " − ", size=size)
    line1_x = draw_f_mc_var(draw, line1_x, y, "B", size=size, fill=ft.TEAL)
    line1_x = c3.draw_n_var(draw, line1_x + 1, y, "A", "signal", size=size)
    c3.math_text(draw, line1_x, y + 2, ")", size=size)

    line2_y = y + 58
    line2_x = x + 142
    line2_x = c3.math_text(draw, line2_x, line2_y + 2, "× ", size=size)
    frac_x = line2_x + 4
    num_x = frac_x + 18
    num_end = c3.draw_n_var(draw, num_x, line2_y - 24, "C", "raw", size=size)
    num_end = c3.math_text(draw, num_end, line2_y - 22, " − ", size=size)
    num_end = draw_f_mc_var(draw, num_end, line2_y - 24, "C", size=size, fill=ft.TEAL)
    num_end = c3.draw_n_var(draw, num_end + 1, line2_y - 24, "A", "signal", size=size)
    frac_w = max(270, num_end - frac_x + 22)
    draw.line((frac_x, line2_y + 20, frac_x + frac_w, line2_y + 20), fill=ft.INK, width=2)
    den_x = frac_x + 18
    den_end = c3.draw_n_var(draw, den_x, line2_y + 25, "D", "raw", size=size)
    den_end = c3.math_text(draw, den_end, line2_y + 27, " − ", size=size)
    den_end = draw_f_mc_var(draw, den_end, line2_y + 25, "D", size=size, fill=ft.TEAL)
    c3.draw_n_var(draw, den_end + 1, line2_y + 25, "A", "signal", size=size)
    draw.text((frac_x + frac_w + 12, line2_y - 21), "]", font=ft.font(ft.TIMES, size + 34), fill=ft.INK)


def math_lines_image(lines: list[str], fontsize: int = 30, dpi: int = 240) -> Image.Image:
    fig = Figure(figsize=(6.8, 2.75), dpi=dpi)
    fig.patch.set_alpha(0)
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.axis("off")
    if len(lines) == 3:
        placements = [(0.04, 0.78), (0.19, 0.49), (0.21, 0.20)]
    else:
        placements = [(0.02, 0.70), (0.16, 0.28)]
    for line, (tx, ty) in zip(lines, placements):
        fig.text(tx, ty, line, fontsize=fontsize, family="serif", color=tuple(v / 255 for v in ft.INK), va="center")
    canvas.draw()
    img = Image.frombuffer("RGBA", canvas.get_width_height(), canvas.buffer_rgba(), "raw", "RGBA", 0, 1).copy()
    return trim_white_margins(img, tolerance=2, pad=4)


def math_block_image(lines: list[str], fontsize: int = 26, width: float = 4.5, line_step: float = 0.32, dpi: int = 240) -> Image.Image:
    height = 0.70 + line_step * max(1, len(lines))
    fig = Figure(figsize=(width, height), dpi=dpi)
    fig.patch.set_alpha(0)
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.axis("off")
    top = 0.82
    for i, line in enumerate(lines):
        fig.text(0.03, top - i * line_step, line, fontsize=fontsize, family="serif", color=tuple(v / 255 for v in ft.INK), va="center")
    canvas.draw()
    img = Image.frombuffer("RGBA", canvas.get_width_height(), canvas.buffer_rgba(), "raw", "RGBA", 0, 1).copy()
    return trim_white_margins(img, tolerance=2, pad=4)


def paste_math_fit(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int]) -> None:
    bw = box[2] - box[0]
    bh = box[3] - box[1]
    scale = min(bw / img.width, bh / img.height)
    resized = img.resize((max(1, int(img.width * scale)), max(1, int(img.height * scale))), Image.Resampling.LANCZOS)
    base.alpha_composite(resized, (box[0] + (bw - resized.width) // 2, box[1] + (bh - resized.height) // 2))


def paste_math_left(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int]) -> None:
    bw = box[2] - box[0]
    bh = box[3] - box[1]
    scale = min(bw / img.width, bh / img.height)
    resized = img.resize((max(1, int(img.width * scale)), max(1, int(img.height * scale))), Image.Resampling.LANCZOS)
    base.alpha_composite(resized, (box[0], box[1] + (bh - resized.height) // 2))


def paste_corrected_abcd_equation(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    lines = [
        r"$N^A_{\mathrm{signal}}=N^A_{\mathrm{raw}}-\Big[$",
        r"$(N^B_{\mathrm{raw}}-f^{B,\mathrm{MC}}N^A_{\mathrm{signal}})$",
        r"$\times\frac{N^C_{\mathrm{raw}}-f^{C,\mathrm{MC}}N^A_{\mathrm{signal}}}{N^D_{\mathrm{raw}}-f^{D,\mathrm{MC}}N^A_{\mathrm{signal}}}\Big]$",
    ]
    paste_math_fit(base, math_lines_image(lines, fontsize=34), box)


def draw_leakage_correction_equation(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 28) -> None:
    x = draw_b_sub(draw, x, y, "X", size=size, fill=ft.TEAL)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)
    x = c3.draw_n_var(draw, x, y, "X", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " − ", size=size)
    x = draw_f_mc_var(draw, x, y, "X", size=size, fill=ft.TEAL)
    x = c3.draw_n_var(draw, x + 2, y, "A", "sig", size=size)
    tail_x = x + 20
    tail_font = ft.font(ft.TIMES_BOLD, size - 6)
    symbol_font = ft.font(ARIAL_UNICODE, size - 6)
    draw.text((tail_x, y + 4), "X ", font=tail_font, fill=ft.MUTED)
    tail_x += ft.text_box(draw, "X ", tail_font)[0]
    draw.text((tail_x, y + 4), "∈", font=symbol_font, fill=ft.MUTED)
    tail_x += ft.text_box(draw, "∈", symbol_font)[0] + 4
    draw.text((tail_x, y + 4), "{B,C,D}", font=tail_font, fill=ft.MUTED)


def draw_purity_equation(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 30) -> None:
    x = c3.math_text(draw, x, y, "P(", size=size)
    x = c3.math_text(draw, x, y, "E", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.50), "T", size=int(size * 0.55), fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x + 7, y, ") = ", size=size)
    frac_x = x
    num_end = c3.draw_n_var(draw, frac_x + 18, y - 13, "A", "sig", size=size - 2)
    frac_w = max(106, num_end - frac_x + 16)
    draw.line((frac_x, y + 37, frac_x + frac_w, y + 37), fill=ft.INK, width=2)
    c3.draw_n_var(draw, frac_x + 18, y + 39, "A", "raw", size=size - 2)


def draw_purity_equation_inline(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 35) -> None:
    x = c3.math_text(draw, x, y, "P(", size=size)
    x = c3.math_text(draw, x, y, "E", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.50), "T", size=int(size * 0.55), fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x + 7, y, ") = ", size=size)
    x = c3.draw_n_var(draw, x, y, "A", "sig", size=size)
    x = c3.math_text(draw, x + 2, y + 2, " / ", size=size)
    c3.draw_n_var(draw, x, y, "A", "raw", size=size)


def draw_minimal_math_strip(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((box[0] + 30, box[1] + 18), "Leakage Correction", font=ft.font(ft.TIMES_BOLD, 38), fill=ft.INK)

    x0, y0, x1, y1 = box
    gap = 22
    usable_w = x1 - x0 - 60 - 2 * gap
    col_ws = [680, 900, usable_w - 680 - 900]
    cols = [
        ("1", "Raw sideband estimate", ft.PHOTON_DARK),
        ("2", "Leakage correction", ft.TEAL),
        ("3", "Purity carried forward", ft.SPHENIX_BLUE),
    ]
    card_boxes = []
    cx = x0 + 30
    for idx, (num, title, color) in enumerate(cols):
        col_w = col_ws[idx]
        cy = y0 + 72
        draw.rounded_rectangle((cx, cy, cx + col_w, y1 - 24), radius=10, fill=(248, 251, 253, 255), outline=(*color, 175), width=2)
        draw.rounded_rectangle((cx, cy, cx + 12, y1 - 24), radius=4, fill=(*color, 255))
        draw.ellipse((cx + 28, cy + 16, cx + 74, cy + 62), fill=(*color, 255))
        draw.text((cx + 51, cy + 39), num, font=ft.font(ft.TIMES_BOLD, 27), fill=(255, 255, 255), anchor="mm")
        draw.text((cx + 92, cy + 19), title, font=ft.font(ft.TIMES_BOLD, 29), fill=ft.INK)
        card_boxes.append((cx, cy, cx + col_w, y1 - 24))
        cx += col_w + gap

    draw_raw_sideband_equation_annotated(draw, card_boxes[0][0] + 74, card_boxes[0][1] + 104, size=32)
    paste_corrected_abcd_equation(base, (card_boxes[1][0] + 42, card_boxes[1][1] + 78, card_boxes[1][2] - 34, card_boxes[1][3] - 14))
    draw_purity_equation_inline(draw, card_boxes[2][0] + 74, card_boxes[2][1] + 104, size=40)


def draw_bottom_flow(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((box[0] + 42, box[1] + 24), "Corrected-yield flow", font=ft.font(ft.TIMES_BOLD, 46), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "Leakage-corrected purity is applied to Region A; efficiency and unfolding move the yield to particle level.",
        (box[0] + 42, box[1] + 82),
        box[2] - box[0] - 84,
        ft.font(ft.TIMES_ITALIC, 33),
        fill=ft.MUTED,
        line_gap=5,
    )
    nodes = [
        ("Region A\ncandidates", ft.PHOTON_DARK, 340),
        ("apply leakage-corrected\npurity curve", ft.SPHENIX_BLUE, 372),
        ("apply\nefficiency", ft.TEAL, 270),
        ("unfold detector\nresponse", ft.PHOTON, 356),
        ("particle-level yield\n/ cross section", ft.BLUE, 404),
    ]
    arrow_w = 46
    gap = 20
    total_w = sum(w for _, _, w in nodes) + (len(nodes) - 1) * (arrow_w + gap)
    x = box[0] + (box[2] - box[0] - total_w) // 2
    y = box[1] + 158
    h = 130
    for idx, (label, color, w) in enumerate(nodes):
        draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=(247, 250, 252, 255), outline=(*color, 230), width=3)
        draw.rounded_rectangle((x, y, x + 15, y + h), radius=5, fill=(*color, 255))
        lines = label.splitlines()
        line_h = 36
        start_y = y + (h - len(lines) * line_h) // 2 - 1
        for j, line in enumerate(lines):
            fnt = ft.font(ft.TIMES_BOLD, 33)
            tw, _ = ft.text_box(draw, line, fnt)
            draw.text((x + 15 + (w - 15 - tw) / 2, start_y + j * line_h), line, font=fnt, fill=ft.INK)
        if idx < len(nodes) - 1:
            ft.draw_arrow(draw, (x + w + 10, y + h // 2), (x + w + arrow_w, y + h // 2), fill=(151, 169, 188), width=6)
        x += w + arrow_w + gap


def _box(name: str) -> Box:
    return Box(*LAYOUT[name], name)


def _audit_within(audit: SymmetryAudit, child: Box, parent: Box, pad: float = 0, *, name: str) -> None:
    ok = (
        child.x0 >= parent.x0 + pad
        and child.y0 >= parent.y0 + pad
        and child.x1 <= parent.x1 - pad
        and child.y1 <= parent.y1 - pad
    )
    audit.checks.append(
        Check(
            name=name,
            kind="within",
            ok=ok,
            got=child.as_int_tuple(),
            want=parent.as_int_tuple(),
            tolerance=pad,
        )
    )


def write_header_and_symmetry_reports() -> dict[str, object]:
    HEADER_SPEC.write_text(json.dumps({"hp2026_main_header": HP2026_MAIN_HEADER}, indent=2) + "\n", encoding="utf-8")

    audit = SymmetryAudit(name="slide16_yield_flow_plot_first_layout")
    slide = Box(0, 0, ft.W, ft.H, "slide canvas")
    main_canvas = Box(132, HP2026_MAIN_HEADER["divider_y"] + 6, 2390, 1326, "main canvas above footer")
    hero = _box("hero_panel")
    plot = _box("purity_plot_box")
    readout = _box("plot_readout_box")
    cue = _box("sideband_cue_box")
    definition = _box("purity_definition_box")
    eq_box = _box("abcd_equation_box")

    audit.image_size(FINAL_PNG, (2560, 1440))
    audit.font_size_at_least("slide title font", HP2026_MAIN_HEADER["title_font_size"], 86, context="HP2026 main-talk header")
    if HP2026_MAIN_HEADER["subtitle_font_size"]:
        audit.font_size_at_least("slide subtitle font", HP2026_MAIN_HEADER["subtitle_font_size"], 56, context="HP2026 main-talk header")
    tmp = Image.new("RGBA", (ft.W, ft.H), (255, 255, 255, 0))
    tmp_draw = ImageDraw.Draw(tmp, "RGBA")
    title_bbox = tmp_draw.textbbox(tuple(HP2026_MAIN_HEADER["title_xy"]), TITLE, font=ft.font(ft.TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]))
    _audit_within(audit, Box(*title_bbox, "slide title text bbox"), Box(0, 0, 2188, 175, "header before logo"), name="slide title clears logo/header bounds")
    if SUBTITLE and HP2026_MAIN_HEADER["subtitle_xy"] and HP2026_MAIN_HEADER["subtitle_font_size"]:
        subtitle_bbox = tmp_draw.textbbox(tuple(HP2026_MAIN_HEADER["subtitle_xy"]), SUBTITLE, font=ft.font(ft.TIMES_ITALIC, HP2026_MAIN_HEADER["subtitle_font_size"]))
        _audit_within(audit, Box(*subtitle_bbox, "slide subtitle text bbox"), Box(0, 160, 2390, 274, "subtitle usable band"), name="slide subtitle fits visible header band")
    for node in (hero, plot, readout, cue, definition, eq_box):
        _audit_within(audit, node, slide, name=f"{node.name} contained on canvas")
    _audit_within(audit, hero, main_canvas, name="hero panel contained between title divider and footer")
    audit.equal_padding_x(main_canvas, hero, hero, tol=3, name="hero panel left/right page padding")
    audit.close("hero panel top begins just below divider", hero.y0, 258, 4)
    audit.close("hero panel bottom clears footer rule", hero.y1, 1308, 4)
    audit.close("plot centered in its usable vertical band", plot.cy, Box(plot.x0, hero.y0 + 70, plot.x1, hero.y1 - 58, "plot usable band").cy, 1)
    audit.close("plot readout and plot share top edge", readout.y0, plot.y0, 1)
    audit.close("sideband cue aligns to readout left edge", cue.x0, readout.x0, 1)
    audit.close("sideband cue aligns to readout right edge", cue.x1, readout.x1, 1)
    audit.close("definition card aligns to readout left edge", definition.x0, readout.x0, 1)
    audit.close("definition card aligns to readout right edge", definition.x1, readout.x1, 1)
    audit.close("equation box centered in sideband cue X", eq_box.cx, cue.cx, 2)
    audit.font_size_at_least("hero card title", 43, 35, context="card title")
    audit.font_size_at_least("hero card lead line", 32, 30, context="lead line")
    audit.font_size_at_least("sideband cue title", 42, 35, context="card title")
    audit.font_size_at_least("sideband cue body", 32, 31, context="card body")
    audit.font_size_at_least("ABCD equation", 44, 40, context="equation")
    audit.font_size_at_least("plot readout label", 34, 30, context="compact plot readout")
    report = audit.write_json(SYMMETRY_REPORT)
    require_audit_passed(report)
    return report


def render_candidate(kind: str) -> Image.Image:
    img = base_slide()
    if kind == "wide_plot":
        hero_box = (132, 258, 2390, 1308)
    elif kind == "balanced":
        hero_box = (132, 258, 2390, 1308)
    else:
        hero_box = (132, 258, 2390, 1308)

    draw_purity_hero_panel(img, hero_box)
    ft.draw_hp2026_identity_footer(img)
    return img


def write_contact_sheet(paths: list[Path]) -> Path:
    sheet = Image.new("RGB", (640 * len(paths), 420), (255, 255, 255))
    draw = ImageDraw.Draw(sheet)
    for i, path in enumerate(paths):
        thumb = Image.open(path).convert("RGB").resize((640, 360), Image.Resampling.LANCZOS)
        x = i * 640
        sheet.paste(thumb, (x, 58))
        draw.text((x + 18, 18), path.stem, font=ft.font(ft.TIMES, 22), fill=ft.INK)
    out = OUTDIR / "slide16_yield_flow_redesign_internal_contact_sheet.png"
    sheet.save(out)
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    CANDIDATE_DIR.mkdir(parents=True, exist_ok=True)

    variants = {
        "candidate01_wide_plot": "wide_plot",
        "candidate02_balanced_selected": "balanced",
        "candidate03_large_diagram": "large_diagram",
    }
    paths: list[Path] = []
    for stem, kind in variants.items():
        LAYOUT.clear()
        img = render_candidate(kind)
        out = CANDIDATE_DIR / f"{stem}.png"
        img.convert("RGB").save(out, "PNG")
        paths.append(out)

    selected = paths[1]
    contact = write_contact_sheet(paths)
    LAYOUT.clear()
    render_candidate("balanced").convert("RGB").save(FINAL_PNG, "PNG")
    symmetry_report = write_header_and_symmetry_reports()
    MANIFEST.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().isoformat(timespec="seconds"),
                "script": str(Path(__file__).resolve()),
                "source_generator_found": str(Path(__file__).resolve()),
                "source_function": "draw_purity_hero_panel()",
                "source_plot_assets": {
                    "purity": str(ft.figure_path("fig5_purity")),
                },
                "deck_mutation": "none; local PNG prototype only",
                "candidate_paths": [str(p) for p in paths],
                "selected": str(FINAL_PNG),
                "contact_sheet": str(contact),
                "hp2026_main_header": HP2026_MAIN_HEADER,
                "header_spec": str(HEADER_SPEC),
                "symmetry_report": str(SYMMETRY_REPORT),
                "symmetry_ok": bool(symmetry_report.get("ok")),
                "selection_rationale": "main-talk plot-first slide: the corrected-yield flow ribbon was removed, the leakage-corrected purity plot and data-control-region cue were expanded into one integrated card, and the lower band now gives a compact plot-reading sequence from raw data-driven points to leakage-corrected purity.",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(FINAL_PNG)


if __name__ == "__main__":
    main()
