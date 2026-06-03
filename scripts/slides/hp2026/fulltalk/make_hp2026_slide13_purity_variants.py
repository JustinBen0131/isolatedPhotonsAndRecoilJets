#!/usr/bin/env python3
"""Render rough layout variants for the HP2026 purity/corrections slide.

These are local PNG candidates only. They reuse current PPG12 paper figure
crops from the full-talk generator and hand-render compact equation/sideband
explanations for visual comparison.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw

import make_hp2026_closing_three_candidates as c3
import make_hp2026_fulltalk_candidates as ft


VARIANT_DIR = ft.OUTPUT / "slide13_purity_layout_variants"
SCRIPT_DIR = ft.SCRIPT_DIR / "slide13_purity_layout_variants"
CONTACT_SHEET = VARIANT_DIR / "hp2026_slide13_purity_variants_contact_sheet.png"
MANIFEST = VARIANT_DIR / "hp2026_slide13_purity_variants_manifest.json"
EQUATION_VARIANT_DIR = VARIANT_DIR / "equation_box_eye_exam"
EQUATION_CONTACT_SHEET = EQUATION_VARIANT_DIR / "hp2026_slide13_equation_box_eye_exam_contact_sheet.png"
EQUATION_MANIFEST = EQUATION_VARIANT_DIR / "hp2026_slide13_equation_box_eye_exam_manifest.json"


def save_script(stem: str, title: str, body: str) -> Path:
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    path = SCRIPT_DIR / f"{stem}_script.md"
    path.write_text(f"# HP2026 Slide 13 Variant Script - {title}\n\n{body.strip()}\n", encoding="utf-8")
    return path


def base() -> Image.Image:
    img = ft.base_slide(
        "From selected candidates to a corrected photon yield",
        "The sideband method measures purity; efficiency and unfolding convert that yield to particle level.",
    )
    ft.add_top_right_sphenix_logo_like_slide2(img)
    return img


def panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], title: str | None = None) -> None:
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    if title:
        draw.text((box[0] + 28, box[1] + 22), title, font=ft.font(ft.TIMES_BOLD, 30), fill=ft.INK)


def shadow_panel(img: Image.Image, box: tuple[int, int, int, int], title: str | None = None) -> ImageDraw.ImageDraw:
    draw = ImageDraw.Draw(img, "RGBA")
    ft.shadow(img, box)
    panel(draw, box, title)
    return draw


def trim_white_margins(img: Image.Image, tolerance: int = 12, pad: int = 8) -> Image.Image:
    rgb = img.convert("RGB")
    white = Image.new("RGB", rgb.size, "white")
    diff = ImageChops.difference(rgb, white).convert("L")
    bbox = diff.point(lambda value: 255 if value > tolerance else 0).getbbox()
    if bbox is None:
        return img
    x0, y0, x1, y1 = bbox
    x0 = max(0, x0 - pad)
    y0 = max(0, y0 - pad)
    x1 = min(img.width, x1 + pad)
    y1 = min(img.height, y1 + pad)
    return img.crop((x0, y0, x1, y1))


def plot_in_panel(
    img: Image.Image,
    key: str,
    box: tuple[int, int, int, int],
    title: str,
    subtitle: str,
    inset_top: int = 96,
    inset: int = 24,
) -> None:
    draw = shadow_panel(img, box, title)
    ft.draw_wrapped(
        draw,
        subtitle,
        (box[0] + 28, box[1] + 60),
        box[2] - box[0] - 56,
        ft.font(ft.TIMES_ITALIC, 21),
        fill=ft.MUTED,
        line_gap=3,
    )
    plot = Image.open(ft.figure_path(key)).convert("RGBA")
    ft.paste_fit(img, plot, (box[0] + inset, box[1] + inset_top, box[2] - inset, box[3] - inset), anchor="center")


def draw_bcd_mini(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], compact: bool = False) -> None:
    x0, y0, x1, y1 = box
    labels = [
        ("A", "iso + tight", ft.PHOTON_DARK, (0, 0)),
        ("B", "noniso + tight", ft.SPHENIX_BLUE, (1, 0)),
        ("C", "iso + nontight", ft.TEAL, (0, 1)),
        ("D", "noniso + nontight", ft.MUTED, (1, 1)),
    ]
    cell_w = (x1 - x0) // 2
    cell_h = (y1 - y0) // 2
    for letter, label, color, (cx, cy) in labels:
        bx0, by0 = x0 + cx * cell_w, y0 + cy * cell_h
        fill = (255, 248, 229) if letter == "A" else (247, 250, 252)
        draw.rounded_rectangle((bx0, by0, bx0 + cell_w - 4, by0 + cell_h - 4), radius=7, fill=(*fill, 255), outline=(*color, 230), width=2)
        draw.text((bx0 + 10, by0 + 8), letter, font=ft.font(ft.TIMES_BOLD, 28 if not compact else 23), fill=color)
        ft.draw_wrapped(
            draw,
            label,
            (bx0 + (42 if not compact else 34), by0 + 12),
            cell_w - (52 if not compact else 42),
            ft.font(ft.TIMES_BOLD, 16 if not compact else 13),
            fill=ft.INK,
            line_gap=2,
        )


def leakage_card(img: Image.Image, box: tuple[int, int, int, int], title: str = "What leakage means") -> None:
    draw = shadow_panel(img, box, title)
    x0, y0, x1, y1 = box
    ft.draw_wrapped(
        draw,
        "Signal in B/C/D means true photons escaped Region A before the sidebands were used as background controls.",
        (x0 + 28, y0 + 62),
        x1 - x0 - 56,
        ft.font(ft.TIMES_ITALIC, 21),
        fill=ft.BLUE,
        line_gap=3,
    )
    draw_bcd_mini(draw, (x0 + 30, y0 + 126, x0 + 310, y0 + 278))
    bullets = [
        ("B", "tight but non-isolated", ft.SPHENIX_BLUE),
        ("C", "isolated but non-tight", ft.TEAL),
        ("D", "fails both axes", ft.MUTED),
    ]
    yy = y0 + 134
    for letter, text, color in bullets:
        draw.text((x0 + 350, yy), letter, font=ft.font(ft.TIMES_BOLD, 27), fill=color)
        ft.draw_wrapped(draw, text, (x0 + 390, yy + 2), x1 - x0 - 414, ft.font(ft.TIMES, 24), fill=ft.INK, line_gap=2)
        yy += 58
    ft.draw_wrapped(
        draw,
        "Eq 3 subtracts leaked signal so B/C/D remain valid background handles.",
        (x0 + 30, y1 - 74),
        x1 - x0 - 60,
        ft.font(ft.TIMES_BOLD, 22),
        fill=ft.BLUE,
        line_gap=4,
    )


def leakage_compact_card(img: Image.Image, box: tuple[int, int, int, int], title: str = "Leakage definitions") -> None:
    draw = shadow_panel(img, box, title)
    x0, y0, x1, y1 = box
    ft.draw_wrapped(
        draw,
        "Signal leakage means true photons escaped Region A before B/C/D were used as background controls.",
        (x0 + 28, y0 + 58),
        x1 - x0 - 56,
        ft.font(ft.TIMES_ITALIC, 20),
        fill=ft.BLUE,
        line_gap=3,
    )
    entries = [
        ("B", "tight, non-isolated", ft.SPHENIX_BLUE),
        ("C", "isolated, non-tight", ft.TEAL),
        ("D", "non-isolated, non-tight", ft.MUTED),
    ]
    yy = y0 + 128
    for letter, text, color in entries:
        draw.rounded_rectangle((x0 + 30, yy - 2, x0 + 68, yy + 36), radius=5, fill=(247, 250, 252, 255), outline=(*color, 220), width=2)
        draw.text((x0 + 42, yy + 2), letter, font=ft.font(ft.TIMES_BOLD, 24), fill=color)
        draw.text((x0 + 84, yy + 2), text, font=ft.font(ft.TIMES, 24), fill=ft.INK)
        yy += 50
    ft.draw_wrapped(
        draw,
        "Eq 3 removes leaked signal before purity is evaluated.",
        (x0 + 30, y1 - 58),
        x1 - x0 - 60,
        ft.font(ft.TIMES_BOLD, 22),
        fill=ft.BLUE,
        line_gap=3,
    )


def draw_arrow_label(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    color: tuple[int, int, int],
    angle_degrees: float = 0.0,
) -> None:
    x, y = xy
    label_font = ft.font(ft.TIMES_BOLD, 21)
    tw, th = ft.text_box(draw, text, label_font)
    pad_x, pad_y = 10, 5
    label = Image.new("RGBA", (int(tw + 2 * pad_x + 6), int(th + 2 * pad_y + 6)), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label, "RGBA")
    label_draw.rounded_rectangle(
        (2, 2, label.width - 2, label.height - 2),
        radius=7,
        fill=(255, 255, 255, 246),
        outline=(*color, 190),
        width=2,
    )
    label_draw.text((pad_x + 3, pad_y + 3), text, font=label_font, fill=color)
    rotated = label.rotate(-angle_degrees, expand=True, resample=Image.Resampling.BICUBIC)
    draw._image.alpha_composite(rotated, (int(x - rotated.width / 2), int(y - rotated.height / 2)))


def draw_leakage_map(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(248, 251, 253, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    header_font = ft.font(ft.TIMES_BOLD, 34)
    definition_font = ft.font(ft.TIMES, 23)
    header_x, header_y = x0 + 26, y0 + 18
    draw.text((header_x, header_y), "Leakage", font=header_font, fill=ft.BLUE)
    header_w = ft.text_box(draw, "Leakage", header_font)[0]
    draw.text(
        (header_x + header_w + 9, header_y + 8),
        "= truth-matched signal photons in PYTHIA MC found in regions B/C/D",
        font=definition_font,
        fill=ft.INK,
    )

    map_x0, map_y0 = x0 + 34, y0 + 96
    available_w = x1 - x0 - 68
    cell_w = min(322, max(280, int((available_w - 160) / 2)))
    cell_h = 104
    h_gap = available_w - 2 * cell_w
    v_gap = 72
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
    for letter, (cx, cy) in cells.items():
        label, detail, color, fill = meta[letter]
        draw.rounded_rectangle((cx, cy, cx + cell_w, cy + cell_h), radius=8, fill=(*fill, 255), outline=(*color, 230), width=3)
        draw.text((cx + 18, cy + 14), letter, font=ft.font(ft.TIMES_BOLD, 35), fill=color)
        ft.draw_wrapped(
            draw,
            label,
            (cx + 70, cy + 14),
            cell_w - 88,
            ft.font(ft.TIMES_BOLD, 22),
            fill=ft.INK,
            line_gap=2,
        )
        ft.draw_wrapped(
            draw,
            detail,
            (cx + 70, cy + 58),
            cell_w - 88,
            ft.font(ft.TIMES, 18),
            fill=ft.MUTED,
            line_gap=2,
        )

    ax, ay = cells["A"]
    bx, by = cells["B"]
    cx, cy = cells["C"]
    dx, dy = cells["D"]
    c_start = (ax + cell_w // 2, ay - 10)
    c_end = (cx + cell_w // 2, cy + cell_h + 10)
    ft.draw_arrow(draw, c_start, c_end, fill=ft.TEAL, width=4)

    b_start = (ax + cell_w + 8, ay + cell_h // 2 + 4)
    b_end = (bx - 12, by + cell_h // 2 + 4)
    ft.draw_arrow(draw, b_start, b_end, fill=ft.SPHENIX_BLUE, width=4)

    d_start = (ax + cell_w + 14, ay + 26)
    d_end = (dx + 4, dy + cell_h - 4)
    ft.draw_arrow(draw, d_start, d_end, fill=ft.MUTED, width=4)

    c_mid = ((c_start[0] + c_end[0]) // 2 + 74, (c_start[1] + c_end[1]) // 2)
    b_mid = ((b_start[0] + b_end[0]) // 2, b_start[1] - 18)
    d_mid = ((d_start[0] + d_end[0]) // 2 + 4, (d_start[1] + d_end[1]) // 2 - 16)
    b_angle = math.degrees(math.atan2(b_end[1] - b_start[1], b_end[0] - b_start[0]))
    d_angle = math.degrees(math.atan2(d_end[1] - d_start[1], d_end[0] - d_start[0]))
    draw_arrow_label(draw, c_mid, "C leakage", ft.TEAL)
    draw_arrow_label(draw, b_mid, "B leakage", ft.SPHENIX_BLUE, angle_degrees=b_angle)
    draw_arrow_label(draw, d_mid, "D leakage", ft.MUTED, angle_degrees=d_angle)

    note_box = (x0 + 34, y1 - 72, x1 - 34, y1 - 22)
    draw.rounded_rectangle(note_box, radius=9, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    icon_x, icon_y = note_box[0] + 20, note_box[1] + 10
    draw.rounded_rectangle(
        (icon_x, icon_y, icon_x + 30, icon_y + 30),
        radius=4,
        fill=(255, 255, 255, 0),
        outline=(*ft.SPHENIX_BLUE, 255),
        width=2,
    )
    draw.line((icon_x + 8, icon_y + 10, icon_x + 22, icon_y + 10), fill=(*ft.SPHENIX_BLUE, 210), width=2)
    draw.line((icon_x + 8, icon_y + 17, icon_x + 21, icon_y + 17), fill=(*ft.SPHENIX_BLUE, 180), width=2)
    draw.line((icon_x + 8, icon_y + 24, icon_x + 17, icon_y + 24), fill=(*ft.SPHENIX_BLUE, 150), width=2)
    note_font = ft.font(ft.TIMES_BOLD, 20)
    note_x = icon_x + 42
    note_y = note_box[1] + 9
    draw.text((note_x, note_y), "Note:", font=note_font, fill=ft.BLUE)
    note_w, _ = ft.text_box(draw, "Note:", note_font)
    body_font = ft.font(ft.TIMES, 20)
    body_x = note_x + note_w + 10
    draw.text((body_x, note_y), "ABCD purity requires isolation-tight ID complementarity", font=body_font, fill=ft.BLUE)
    draw.text((body_x, note_y + 22), "isolation is not used in BDT training.", font=body_font, fill=ft.BLUE)


def draw_compact_eq3_leakage(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    size = 22
    x = c3.math_text(draw, x, y, "For X = B, C, D:  ", size=20, fill=ft.MUTED, bold=True)
    x = c3.math_text(draw, x, y, "X", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x, y + 13, "corr", size=14, fill=ft.BLUE)
    x = c3.math_text(draw, x + 8, y, " = ", size=size)
    x = c3.draw_n_var(draw, x, y - 2, "X", "raw", size=size)
    x = c3.math_text(draw, x, y, " − ", size=size)
    x = c3.draw_f_var(draw, x, y - 2, "X,MC", size=size)
    c3.draw_n_var(draw, x, y - 2, "A", "signal", size=size)


def draw_eq4_compact(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    size = 23
    x = c3.math_text(draw, x, y, "Purity(", size=size)
    x = c3.math_text(draw, x, y, "P", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x, y, ") = ", size=size)
    frac_x = x
    num_end = c3.draw_n_var(draw, frac_x + 16, y - 12, "A", "signal", size=21)
    frac_w = max(86, num_end - frac_x + 14)
    draw.line((frac_x, y + 27, frac_x + frac_w, y + 27), fill=ft.INK, width=2)
    c3.draw_n_var(draw, frac_x + 16, y + 27, "A", "raw", size=21)


def draw_corr_var(draw: ImageDraw.ImageDraw, x: int, y: int, letter: str, size: int = 29, fill=ft.INK) -> int:
    base = ft.font(ft.TIMES_BOLD, size)
    sub = ft.font(ft.TIMES, max(14, int(size * 0.55)))
    draw.text((x, y), letter, font=base, fill=fill)
    bw = ft.text_box(draw, letter, base)[0]
    draw.text((x + bw + 2, y + int(size * 0.50)), "corr", font=sub, fill=fill)
    return x + bw + ft.text_box(draw, "corr", sub)[0] + 10


def draw_eq3_rail(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    size = 30
    x = c3.draw_n_var(draw, x, y - 2, "A", "signal", size=size)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)
    x = c3.draw_n_var(draw, x, y - 2, "A", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " − [", size=size)
    x = draw_corr_var(draw, x, y + 2, "B", size=28)
    x = c3.math_text(draw, x, y + 2, " × ", size=size)
    draw.text((x, y - 9), "(", font=ft.font(ft.TIMES, 58), fill=ft.INK)
    x += 26
    frac_x = x
    num_end = draw_corr_var(draw, frac_x + 14, y - 10, "C", size=26)
    frac_w = max(88, num_end - frac_x + 12)
    draw.line((frac_x, y + 38, frac_x + frac_w, y + 38), fill=ft.INK, width=2)
    draw_corr_var(draw, frac_x + 14, y + 42, "D", size=26)
    draw.text((frac_x + frac_w + 8, y - 9), ")]", font=ft.font(ft.TIMES, 58), fill=ft.INK)


def draw_sideband_corr_definition(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    size = 23
    x = c3.math_text(draw, x, y, "For X = B, C, D:   ", size=20, fill=ft.MUTED, bold=True)
    x = draw_corr_var(draw, x, y, "X", size=23, fill=ft.BLUE)
    x = c3.math_text(draw, x, y, " = ", size=size)
    x = c3.draw_n_var(draw, x, y - 2, "X", "raw", size=size)
    x = c3.math_text(draw, x, y, " − ", size=size)
    x = c3.draw_f_var(draw, x, y - 2, "X,MC", size=size)
    c3.draw_n_var(draw, x, y - 2, "A", "signal", size=size)


def draw_var(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    symbol: str,
    sub: str | None = None,
    sup: str | None = None,
    size: int = 30,
    fill: tuple[int, int, int] = ft.INK,
    bold: bool = False,
) -> int:
    base = ft.font(ft.TIMES_BOLD if bold else ft.TIMES_ITALIC, size)
    small = ft.font(ft.TIMES, max(13, int(size * 0.54)))
    draw.text((x, y), symbol, font=base, fill=fill)
    bw = ft.text_box(draw, symbol, base)[0]
    extra_w = 0
    if sup:
        draw.text((x + bw + 1, y - int(size * 0.28)), sup, font=small, fill=fill)
        extra_w = max(extra_w, ft.text_box(draw, sup, small)[0])
    if sub:
        draw.text((x + bw + 1, y + int(size * 0.48)), sub, font=small, fill=fill)
        extra_w = max(extra_w, ft.text_box(draw, sub, small)[0])
    return x + bw + extra_w + 8


def draw_region_raw(draw: ImageDraw.ImageDraw, x: int, y: int, region: str, size: int = 28) -> int:
    return draw_var(draw, x, y, region, sup="raw", size=size, bold=True)


def draw_s_a(draw: ImageDraw.ImageDraw, x: int, y: int, sup: str | None = None, size: int = 30) -> int:
    return draw_var(draw, x, y, "S", sub="A", sup=sup, size=size)


def draw_raw_abcd_formula(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 28) -> None:
    x = draw_s_a(draw, x, y, sup="raw", size=size)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)
    x = draw_region_raw(draw, x, y, "A", size=size)
    x = c3.math_text(draw, x, y + 2, " − ", size=size)
    x = draw_region_raw(draw, x, y, "B", size=size)
    x = c3.math_text(draw, x, y + 2, " × ", size=size)
    draw.text((x, y - 9), "(", font=ft.font(ft.TIMES, int(size * 2.0)), fill=ft.INK)
    x += int(size * 0.9)
    frac_x = x
    num_end = draw_region_raw(draw, frac_x + 10, y - 12, "C", size=max(22, size - 3))
    frac_w = max(72, num_end - frac_x + 10)
    draw.line((frac_x, y + int(size * 1.22), frac_x + frac_w, y + int(size * 1.22)), fill=ft.INK, width=2)
    draw_region_raw(draw, frac_x + 10, y + int(size * 1.23), "D", size=max(22, size - 3))
    draw.text((frac_x + frac_w + 7, y - 9), ")", font=ft.font(ft.TIMES, int(size * 2.0)), fill=ft.INK)


def draw_leakage_factor_formula(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 28) -> None:
    items = [("B", ft.SPHENIX_BLUE), ("C", ft.TEAL), ("D", ft.MUTED)]
    for idx, (letter, color) in enumerate(items):
        x = c3.math_text(draw, x, y, "f", size=size, fill=color, bold=True)
        x = c3.math_text(draw, x - 2, y + int(size * 0.45), letter, size=int(size * 0.56), fill=color, bold=True)
        x = c3.math_text(draw, x + 6, y, " = ", size=size)
        frac_x = x
        x = c3.math_text(draw, frac_x + 6, y - 12, letter, size=size - 1, fill=color, bold=True)
        x = c3.math_text(draw, x, y - 2, "sig", size=int(size * 0.52))
        frac_w = max(68, x - frac_x + 8)
        draw.line((frac_x, y + int(size * 1.10), frac_x + frac_w, y + int(size * 1.10)), fill=ft.INK, width=2)
        x = c3.math_text(draw, frac_x + 8, y + int(size * 1.10), "A", size=size - 1, fill=ft.PHOTON_DARK, bold=True)
        x = c3.math_text(draw, x, y + int(size * 1.20), "sig", size=int(size * 0.52))
        x = frac_x + frac_w + (24 if idx < len(items) - 1 else 0)


def draw_background_sideband_formula(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 23) -> None:
    start_x = x
    rows = [
        ("B", ft.SPHENIX_BLUE, "B"),
        ("C", ft.TEAL, "C"),
        ("D", ft.MUTED, "D"),
    ]
    for idx, (letter, color, raw) in enumerate(rows):
        yy = y + idx * 38
        xx = c3.math_text(draw, start_x, yy, "b", size=size, fill=color, bold=True)
        xx = c3.math_text(draw, xx - 2, yy + int(size * 0.48), letter, size=int(size * 0.56), fill=color, bold=True)
        xx = c3.math_text(draw, xx + 6, yy, " = ", size=size)
        xx = draw_region_raw(draw, xx, yy, raw, size=size)
        xx = c3.math_text(draw, xx, yy, " − f", size=size)
        xx = c3.math_text(draw, xx - 2, yy + int(size * 0.48), letter, size=int(size * 0.56), fill=color, bold=True)
        draw_s_a(draw, xx + 6, yy, size=size)

    xx = c3.math_text(draw, start_x + 248, y + 34, "b", size=size + 1, fill=ft.PHOTON_DARK, bold=True)
    xx = c3.math_text(draw, xx - 2, y + 34 + int(size * 0.48), "A", size=int(size * 0.56), fill=ft.PHOTON_DARK, bold=True)
    xx = c3.math_text(draw, xx + 6, y + 34, " = ", size=size + 1)
    xx = c3.math_text(draw, xx, y + 34, "b", size=size + 1, fill=ft.SPHENIX_BLUE, bold=True)
    xx = c3.math_text(draw, xx - 2, y + 34 + int(size * 0.48), "B", size=int(size * 0.56), fill=ft.SPHENIX_BLUE, bold=True)
    xx = c3.math_text(draw, xx + 6, y + 34, " × b", size=size + 1)
    xx = c3.math_text(draw, xx - 2, y + 34 + int(size * 0.48), "C", size=int(size * 0.56), fill=ft.TEAL, bold=True)
    xx = c3.math_text(draw, xx + 8, y + 34, "/ b", size=size + 1)
    c3.math_text(draw, xx - 2, y + 34 + int(size * 0.48), "D", size=int(size * 0.56), fill=ft.MUTED, bold=True)


def draw_corrected_purity_formula(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 27) -> None:
    x = draw_s_a(draw, x, y, size=size)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)
    x = draw_region_raw(draw, x, y, "A", size=size)
    x = c3.math_text(draw, x, y + 2, " − b", size=size)
    x = c3.math_text(draw, x - 2, y + int(size * 0.48), "A", size=int(size * 0.56), fill=ft.PHOTON_DARK, bold=True)
    x = c3.math_text(draw, x + 22, y + 2, "     P", size=size, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.48), "corr", size=int(size * 0.56), fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x + 8, y + 2, " = ", size=size)
    frac_x = x
    end = draw_s_a(draw, frac_x + 18, y - 12, size=size - 1)
    frac_w = max(88, end - frac_x + 16)
    draw.line((frac_x, y + int(size * 1.18), frac_x + frac_w, y + int(size * 1.18)), fill=ft.INK, width=2)
    draw_region_raw(draw, frac_x + 18, y + int(size * 1.16), "A", size=size - 1)


def draw_b_var(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    sub: str,
    sup: str | None = None,
    size: int = 28,
    fill: tuple[int, int, int] = ft.INK,
    bold: bool = False,
) -> int:
    return draw_var(draw, x, y, "b", sub=sub, sup=sup, size=size, fill=fill, bold=bold)


def draw_leakage_factor_compact(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 25) -> None:
    x = c3.math_text(draw, x, y, "f", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.48), "X", size=int(size * 0.56), fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x + 6, y, " = ", size=size)
    frac_x = x
    x = c3.math_text(draw, frac_x + 12, y - 12, "X", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x, y - 1, "sig", size=int(size * 0.52))
    frac_w = max(86, x - frac_x + 12)
    draw.line((frac_x, y + int(size * 1.10), frac_x + frac_w, y + int(size * 1.10)), fill=ft.INK, width=2)
    x = c3.math_text(draw, frac_x + 12, y + int(size * 1.10), "A", size=size, fill=ft.PHOTON_DARK, bold=True)
    c3.math_text(draw, x, y + int(size * 1.20), "sig", size=int(size * 0.52))
    c3.math_text(draw, frac_x + frac_w + 24, y, "for X = B, C, D", size=size - 2, fill=ft.MUTED, bold=True)


def draw_background_compact(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 24) -> None:
    x0 = x
    x = c3.math_text(draw, x, y, "b", size=size, fill=ft.TEAL, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.48), "X", size=int(size * 0.56), fill=ft.TEAL, bold=True)
    x = c3.math_text(draw, x + 6, y, " = ", size=size)
    x = c3.math_text(draw, x, y, "X", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x, y - 6, "raw", size=int(size * 0.52))
    x = c3.math_text(draw, x + 8, y, " − f", size=size)
    x = c3.math_text(draw, x - 2, y + int(size * 0.48), "X", size=int(size * 0.56), fill=ft.BLUE, bold=True)
    x = draw_s_a(draw, x + 8, y, size=size)
    x = c3.math_text(draw, x + 24, y, "for X = B, C, D", size=size - 2, fill=ft.MUTED, bold=True)
    x = c3.math_text(draw, x0 + 442, y, "b", size=size, fill=ft.PHOTON_DARK, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.48), "A", size=int(size * 0.56), fill=ft.PHOTON_DARK, bold=True)
    x = c3.math_text(draw, x + 6, y, " = ", size=size)
    x = c3.math_text(draw, x, y, "b", size=size, fill=ft.SPHENIX_BLUE, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.48), "B", size=int(size * 0.56), fill=ft.SPHENIX_BLUE, bold=True)
    x = c3.math_text(draw, x + 6, y, "b", size=size, fill=ft.TEAL, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.48), "C", size=int(size * 0.56), fill=ft.TEAL, bold=True)
    x = c3.math_text(draw, x + 6, y, "/b", size=size)
    c3.math_text(draw, x - 2, y + int(size * 0.48), "D", size=int(size * 0.56), fill=ft.MUTED, bold=True)


def draw_n_fraction(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    num_sup: str,
    num_sub: str,
    den_sup: str,
    den_sub: str,
    size: int = 22,
    fill: tuple[int, int, int] = ft.INK,
) -> int:
    num_end = c3.draw_n_var(draw, x + 12, y - 10, num_sup, num_sub, size=size, fill=fill)
    frac_w = max(84, num_end - x + 14)
    draw.line((x, y + 30, x + frac_w, y + 30), fill=ft.INK, width=2)
    c3.draw_n_var(draw, x + 12, y + 33, den_sup, den_sub, size=size, fill=fill)
    return x + frac_w + 8


def draw_b_fraction(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    num_sub: str,
    den_sub: str,
    size: int = 24,
) -> int:
    num_end = draw_b_var(draw, x + 14, y - 9, num_sub, size=size, fill=ft.TEAL, bold=True)
    frac_w = max(62, num_end - x + 10)
    draw.line((x, y + 29, x + frac_w, y + 29), fill=ft.INK, width=2)
    draw_b_var(draw, x + 14, y + 31, den_sub, size=size, fill=ft.MUTED, bold=True)
    return x + frac_w + 8


def draw_eq2_physical(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 23) -> None:
    x = c3.draw_n_var(draw, x, y, "A,RAW", "signal", size=size)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)
    x = c3.draw_n_var(draw, x, y, "A", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " − ", size=size)
    x = c3.draw_n_var(draw, x, y, "B", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " × ", size=size)
    x = c3.draw_n_var(draw, x, y, "C", "raw", size=size)
    x = c3.math_text(draw, x + 2, y + 2, "/", size=size)
    c3.draw_n_var(draw, x + 2, y, "D", "raw", size=size)


def draw_eq3_step_form(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 21) -> None:
    x1 = c3.draw_f_var(draw, x, y, "X,MC", size=size + 1, fill=ft.BLUE)
    x1 = c3.math_text(draw, x1, y + 1, " = ", size=size + 1)
    x1 = c3.draw_n_var(draw, x1, y, "X", "signal", size=size + 1, fill=ft.BLUE)
    x1 = c3.math_text(draw, x1 + 2, y + 1, "/", size=size + 1)
    x1 = c3.draw_n_var(draw, x1 + 2, y, "A", "signal", size=size + 1, fill=ft.BLUE)
    c3.math_text(draw, x1 + 18, y + 1, "X = B, C, D", size=size - 2, fill=ft.MUTED, bold=True)

    y2 = y + 42
    x2 = draw_b_var(draw, x, y2, "X", size=size + 1, fill=ft.TEAL, bold=True)
    x2 = c3.math_text(draw, x2, y2 + 1, " = ", size=size + 1)
    x2 = c3.draw_n_var(draw, x2, y2, "X", "raw", size=size + 1)
    x2 = c3.math_text(draw, x2, y2 + 1, " − ", size=size + 1)
    x2 = c3.draw_f_var(draw, x2, y2, "X,MC", size=size + 1, fill=ft.BLUE)
    c3.draw_n_var(draw, x2 + 8, y2, "A", "signal", size=size + 1)


def draw_eq3_eq4_final(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 23) -> None:
    x1 = c3.draw_n_var(draw, x, y, "A", "signal", size=size)
    x1 = c3.math_text(draw, x1, y + 2, " = ", size=size)
    x1 = c3.draw_n_var(draw, x1, y, "A", "raw", size=size)
    x1 = c3.math_text(draw, x1, y + 2, " − ", size=size)
    x1 = draw_b_var(draw, x1, y + 2, "B", size=size, fill=ft.SPHENIX_BLUE, bold=True)
    x1 = c3.math_text(draw, x1, y + 2, " × ", size=size)
    x1 = draw_b_var(draw, x1, y + 2, "C", size=size, fill=ft.TEAL, bold=True)
    x1 = c3.math_text(draw, x1 + 2, y + 2, "/", size=size)
    x1 = draw_b_var(draw, x1 + 2, y + 2, "D", size=size, fill=ft.MUTED, bold=True)

    x2 = c3.math_text(draw, x1 + 26, y + 2, "Purity(", size=size)
    x2 = c3.math_text(draw, x2, y + 2, "P", size=size, fill=ft.BLUE, bold=True)
    x2 = c3.math_text(draw, x2, y + 2, ") = ", size=size)
    x2 = c3.draw_n_var(draw, x2, y, "A", "signal", size=size, fill=ft.BLUE)
    x2 = c3.math_text(draw, x2 + 2, y + 2, "/", size=size)
    c3.draw_n_var(draw, x2 + 2, y, "A", "raw", size=size, fill=ft.BLUE)


def draw_symbol_key_math(draw: ImageDraw.ImageDraw, key: str, x: int, y: int) -> None:
    if key == "raw":
        c3.draw_n_var(draw, x, y, "X", "raw", size=22)
    elif key == "signal":
        c3.draw_n_var(draw, x, y, "A", "signal", size=22, fill=ft.PHOTON_DARK)
    elif key == "leakage":
        c3.draw_f_var(draw, x, y, "X,MC", size=24, fill=ft.BLUE)
    elif key == "background":
        draw_b_var(draw, x, y + 2, "X", size=24, fill=ft.TEAL, bold=True)


def draw_symbol_key(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(210, 222, 233, 255), width=2)
    draw.text((x0 + 22, y0 + 14), "Symbol key", font=ft.font(ft.TIMES_BOLD, 24), fill=ft.INK)
    draw.text((x0 + 22, y0 + 42), "X = A, B, C, or D", font=ft.font(ft.TIMES_ITALIC, 18), fill=ft.MUTED)

    entries = [
        ("raw", "raw candidates", "counted in region X", ft.INK),
        ("signal", "signal in A", "after background subtraction", ft.PHOTON_DARK),
        ("leakage", "leakage fraction", "truth-matched signal in SIM", ft.BLUE),
        ("background", "background-only", "sideband count after leakage", ft.TEAL),
    ]
    yy = y0 + 72
    for key, lead, detail, color in entries:
        draw.rounded_rectangle((x0 + 18, yy, x1 - 18, yy + 44), radius=8, fill=(248, 251, 253, 255), outline=(*color, 135), width=2)
        draw_symbol_key_math(draw, key, x0 + 34, yy + 8)
        draw.text((x0 + 122, yy + 5), lead, font=ft.font(ft.TIMES_BOLD, 18), fill=ft.INK)
        draw.text((x0 + 122, yy + 26), detail, font=ft.font(ft.TIMES, 15), fill=ft.MUTED)
        yy += 54


def draw_physical_calculation(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(210, 222, 233, 255), width=2)
    draw.text((x0 + 22, y0 + 14), "Physical read of the algebra", font=ft.font(ft.TIMES_BOLD, 24), fill=ft.INK)
    draw.text((x0 + 22, y0 + 42), "count A, clean B/C/D, then form the signal fraction", font=ft.font(ft.TIMES_ITALIC, 18), fill=ft.MUTED)

    rows = [
        {
            "badge": "Eq 2",
            "title": "raw ABCD estimate",
            "note": "subtract a sideband-predicted background from A",
            "color": ft.PHOTON_DARK,
            "height": 62,
            "renderer": draw_eq2_physical,
            "size": 21,
        },
        {
            "badge": "Eq 3",
            "title": "leakage-corrected sidebands",
            "note": "remove true photons from B/C/D before using them as controls",
            "color": ft.TEAL,
            "height": 88,
            "renderer": draw_eq3_step_form,
            "size": 19,
        },
        {
            "badge": "Eq 4",
            "title": "corrected signal and purity",
            "note": "the purity curve is this fraction versus photon energy",
            "color": ft.BLUE,
            "height": 66,
            "renderer": draw_eq3_eq4_final,
            "size": 20,
        },
    ]
    ry = y0 + 68
    for idx, row in enumerate(rows):
        color = row["color"]
        height = row["height"]
        draw.rounded_rectangle((x0 + 20, ry, x1 - 20, ry + height), radius=8, fill=(248, 251, 253, 255), outline=(*color, 150), width=2)
        draw.rounded_rectangle((x0 + 20, ry, x0 + 34, ry + height), radius=4, fill=(*color, 255))
        draw.text((x0 + 50, ry + 9), row["badge"], font=ft.font(ft.TIMES_BOLD, 21), fill=color)
        draw.text((x0 + 128, ry + 8), row["title"], font=ft.font(ft.TIMES_BOLD, 19), fill=ft.INK)
        draw.text((x0 + 128, ry + 31), row["note"], font=ft.font(ft.TIMES_ITALIC, 15), fill=ft.MUTED)
        row["renderer"](draw, x0 + 514, ry + (15 if idx != 1 else 8), size=row["size"])
        if idx < len(rows) - 1:
            ft.draw_arrow(draw, (x0 + 448, ry + height + 1), (x0 + 448, ry + height + 11), fill=(151, 169, 188), width=3)
        ry += height + 8


def paste_scaled_equation(
    img: Image.Image,
    renderer,
    xy: tuple[int, int],
    scale: float,
    canvas_size: tuple[int, int] = (760, 190),
) -> tuple[int, int, int, int]:
    canvas = Image.new("RGBA", canvas_size, (255, 255, 255, 0))
    renderer(ImageDraw.Draw(canvas, "RGBA"), 18, 28)
    bbox = canvas.getbbox()
    if not bbox:
        return (xy[0], xy[1], xy[0], xy[1])
    crop = canvas.crop(bbox)
    fitted = crop.resize((int(crop.width * scale), int(crop.height * scale)), Image.Resampling.LANCZOS)
    img.alpha_composite(fitted, xy)
    return (xy[0], xy[1], xy[0] + fitted.width, xy[1] + fitted.height)


def draw_equation_method_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=10, fill=(248, 251, 253, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 22, y0 + 18), "Equations used to correct Region A", font=ft.font(ft.TIMES_BOLD, 25), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "Eq 2 is the ABCD estimate; Eq 3 subtracts truth-tagged signal leakage from B/C/D; Eq 4 is the purity applied to the yield.",
        (x0 + 22, y0 + 52),
        x1 - x0 - 44,
        ft.font(ft.TIMES_ITALIC, 18),
        fill=ft.BLUE,
        line_gap=2,
    )
    sections = [
        ("Eq 2: raw sidebands estimate the background", ft.PHOTON_DARK, y0 + 112, 96, c3.draw_eq2, 42),
        ("Eq 3: leakage-corrected sidebands", ft.SPHENIX_BLUE, y0 + 220, 98, draw_compact_eq3_leakage, 42),
        ("Eq 4: purity carried into the result", ft.TEAL, y0 + 324, 96, draw_eq4_compact, 36),
    ]
    for heading, color, sy, height, renderer, render_offset in sections:
        draw.rounded_rectangle((x0 + 22, sy, x1 - 22, sy + height), radius=8, fill=(255, 255, 255, 255), outline=(*color, 175), width=2)
        draw.rounded_rectangle((x0 + 22, sy, x0 + 34, sy + height), radius=4, fill=(*color, 255))
        draw.text((x0 + 48, sy + 10), heading, font=ft.font(ft.TIMES_BOLD, 19), fill=color)
        renderer(draw, x0 + 48, sy + render_offset)


def draw_wide_equation_panel(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(248, 251, 253, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 28, y0 + 20), "Equations used to correct Region A", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "Read the paper equations as a physical bookkeeping chain: raw candidates, simulated signal leakage, background-only sidebands, then purity.",
        (x0 + 28, y0 + 60),
        x1 - x0 - 56,
        ft.font(ft.TIMES_ITALIC, 22),
        fill=ft.BLUE,
        line_gap=3,
    )

    content = (x0 + 28, y0 + 116, x1 - 28, y1 - 26)
    legend = (content[0], content[1], content[0] + 430, content[3])
    calculation = (legend[2] + 24, content[1], content[2], content[3])
    draw_symbol_key(draw, legend)
    draw_physical_calculation(draw, calculation)


def draw_eq2_abcd_inline(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 27) -> None:
    draw_eq2_physical(draw, x + 118, y - 1, size=32)


def draw_eq3_sideband_inline(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 25) -> None:
    start_x = x + 58
    size = 26
    x = start_x
    x = draw_b_var(draw, x, y, "X", size=size, fill=ft.INK, bold=True)
    x = c3.math_text(draw, x, y + 1, " = ", size=size)
    x = c3.draw_n_var(draw, x, y, "X", "raw", size=size)
    x = c3.math_text(draw, x, y + 1, " − ", size=size)
    x = c3.draw_f_var(draw, x, y, "X,MC", size=size, fill=ft.BLUE)
    x = c3.draw_n_var(draw, x + 8, y, "A,RAW", "signal", size=size)
    c3.math_text(draw, x + 10, y + 1, ", X = B,C,D", size=17, fill=ft.MUTED, bold=True)

    x = start_x + 482
    x = c3.draw_n_var(draw, x, y, "A", "signal", size=size)
    x = c3.math_text(draw, x, y + 1, " = ", size=size)
    x = c3.draw_n_var(draw, x, y, "A,RAW", "signal", size=size)
    x = c3.math_text(draw, x, y + 1, " − ", size=size)
    x = draw_b_var(draw, x, y + 1, "B", size=size, fill=ft.INK, bold=True)
    x = c3.math_text(draw, x, y + 1, " × ", size=size)
    x = draw_b_var(draw, x, y + 1, "C", size=size, fill=ft.INK, bold=True)
    x = c3.math_text(draw, x + 2, y + 1, "/", size=size)
    draw_b_var(draw, x + 2, y + 1, "D", size=size, fill=ft.INK, bold=True)


def draw_eq3_signal_inline(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 26) -> None:
    x = c3.draw_n_var(draw, x, y, "A", "signal", size=size)
    x = c3.math_text(draw, x, y + 1, " = ", size=size)
    x = c3.draw_n_var(draw, x, y, "A", "raw", size=size)
    x = c3.math_text(draw, x, y + 1, " − ", size=size)
    x = draw_b_var(draw, x, y + 1, "B", size=size, fill=ft.SPHENIX_BLUE, bold=True)
    x = c3.math_text(draw, x, y + 1, " × ", size=size)
    x = draw_b_var(draw, x, y + 1, "C", size=size, fill=ft.TEAL, bold=True)
    x = c3.math_text(draw, x + 2, y + 1, "/", size=size)
    draw_b_var(draw, x + 2, y + 1, "D", size=size, fill=ft.MUTED, bold=True)


def draw_eq4_purity_inline(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 27) -> None:
    x = c3.math_text(draw, x, y + 1, "Purity(", size=size)
    x = c3.math_text(draw, x, y + 1, "P", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x, y + 1, ") = ", size=size)
    x = c3.draw_n_var(draw, x, y, "A", "signal", size=size, fill=ft.BLUE)
    x = c3.math_text(draw, x + 2, y + 1, "/", size=size)
    c3.draw_n_var(draw, x + 2, y, "A", "raw", size=size, fill=ft.BLUE)


def draw_purity_comparison_inline(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 24) -> None:
    corrected_blue = (88, 118, 255)
    size = 31
    raw_y = y - 2
    corr_y = raw_y
    x = x + 116
    corr_x = x + 430

    xx = c3.math_text(draw, x, raw_y + 2, "P", size=size, fill=ft.INK, bold=True)
    xx = c3.math_text(draw, xx - 2, raw_y + int(size * 0.54), "raw", size=max(13, int(size * 0.50)), fill=ft.INK)
    xx = c3.math_text(draw, xx + 10, raw_y + 2, " = ", size=size)
    xx = c3.draw_n_var(draw, xx, raw_y, "A,RAW", "signal", size=size, fill=ft.INK)
    xx = c3.math_text(draw, xx + 2, raw_y + 2, "/", size=size)
    c3.draw_n_var(draw, xx + 2, raw_y, "A", "raw", size=size, fill=ft.INK)

    xx = c3.math_text(draw, corr_x, corr_y + 2, "P", size=size, fill=corrected_blue, bold=True)
    xx = c3.math_text(draw, xx - 2, corr_y + int(size * 0.54), "corr", size=max(13, int(size * 0.50)), fill=corrected_blue)
    xx = c3.math_text(draw, xx + 10, corr_y + 2, " = ", size=size, fill=corrected_blue)
    xx = c3.draw_n_var(draw, xx, corr_y, "A", "signal", size=size, fill=corrected_blue)
    xx = c3.math_text(draw, xx + 2, corr_y + 2, "/", size=size, fill=corrected_blue)
    c3.draw_n_var(draw, xx + 2, corr_y, "A", "raw", size=size, fill=corrected_blue)


def draw_definition_ribbon_large_rows(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=9, fill=(255, 255, 255, 255), outline=(210, 222, 233, 255), width=2)
    sections = [
        (x0 + 22, "raw", ft.INK),
        (x0 + 346, "signal", ft.INK),
        (x0 + 690, "f", ft.BLUE),
        (x0 + 1218, "b", ft.INK),
    ]

    c3.draw_n_var(draw, sections[0][0], y0 + 24, "X", "raw", size=26, fill=ft.INK)
    draw.text((sections[0][0] + 110, y0 + 31), "counted candidates", font=ft.font(ft.TIMES, 23), fill=ft.MUTED)

    c3.draw_n_var(draw, sections[1][0], y0 + 24, "A", "signal", size=26, fill=ft.INK)
    draw.text((sections[1][0] + 118, y0 + 31), "corrected signal in A", font=ft.font(ft.TIMES, 23), fill=ft.MUTED)

    xx = c3.draw_f_var(draw, sections[2][0], y0 + 24, "X,MC", size=27, fill=ft.BLUE)
    xx = c3.math_text(draw, xx, y0 + 29, " = ", size=24)
    xx = c3.draw_n_var(draw, xx, y0 + 24, "X", "leak", size=24, fill=ft.BLUE)
    xx = c3.math_text(draw, xx + 1, y0 + 29, "/", size=24)
    xx = c3.draw_n_var(draw, xx + 1, y0 + 24, "A", "signal", size=24, fill=ft.BLUE)
    draw.text((xx + 12, y0 + 31), "signal leakage fraction", font=ft.font(ft.TIMES, 23), fill=ft.MUTED)

    xx = draw_b_var(draw, sections[3][0], y0 + 27, "X", size=29, fill=ft.INK, bold=True)
    draw.text((xx + 12, y0 + 31), "corrected sideband count", font=ft.font(ft.TIMES, 23), fill=ft.MUTED)


def draw_equation_box_frame(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], subtitle: str) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(248, 251, 253, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 28, y0 + 20), "Equations used to correct Region A", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    return (x0 + 28, y0 + 76, x1 - 28, y1 - 26)


def draw_definition_chip(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    lead: str,
    detail: str,
    color: tuple[int, int, int],
) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=8, fill=(255, 255, 255, 255), outline=(*color, 160), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 10, y1), radius=4, fill=(*color, 255))
    draw.text((x0 + 22, y0 + 10), lead, font=ft.font(ft.TIMES_BOLD, 19), fill=ft.INK)
    ft.draw_wrapped(draw, detail, (x0 + 22, y0 + 34), x1 - x0 - 44, ft.font(ft.TIMES, 16), fill=ft.MUTED, line_gap=1)


def draw_equation_row(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    badge: str,
    title: str,
    note: str,
    color: tuple[int, int, int],
    renderer,
    formula_size: int,
) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=9, fill=(255, 255, 255, 255), outline=(*color, 170), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 14, y1), radius=4, fill=(*color, 255))
    draw.text((x0 + 26, y0 + 9), badge, font=ft.font(ft.TIMES_BOLD, 27), fill=color)
    draw.text((x0 + 118, y0 + 9), title, font=ft.font(ft.TIMES_BOLD, 24), fill=ft.INK)
    draw.text((x0 + 118, y0 + 36), note, font=ft.font(ft.TIMES_ITALIC, 22), fill=ft.MUTED)
    renderer(draw, x0 + 618, y0 + 18, size=formula_size)


def draw_equation_box_simple_stack(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = draw_equation_box_frame(draw, box, "Simplest read: estimate the signal in A, correct sideband leakage, then convert to purity.")
    ribbon_h = 50
    draw.rounded_rectangle((x0, y0, x1, y0 + ribbon_h), radius=9, fill=(255, 255, 255, 255), outline=(210, 222, 233, 255), width=2)
    definitions = [
        ("A", "selected region", ft.PHOTON_DARK),
        ("B/C/D", "control regions", ft.SPHENIX_BLUE),
        ("f", "signal leakage fraction", ft.BLUE),
        ("b", "background-only sideband", ft.TEAL),
    ]
    xx = x0 + 26
    for lead, detail, color in definitions:
        draw.text((xx, y0 + 12), lead, font=ft.font(ft.TIMES_BOLD, 21), fill=color)
        draw.text((xx + 70, y0 + 14), detail, font=ft.font(ft.TIMES, 17), fill=ft.MUTED)
        xx += 320

    rows = [
        ("Eq 2", "raw ABCD estimate", "background in A estimated from B, C, D", ft.PHOTON_DARK, draw_eq2_abcd_inline, 24),
        ("Eq 3", "remove leaked signal from sidebands", "make B/C/D background-like before using them", ft.TEAL, draw_eq3_sideband_inline, 24),
        ("Eq 4", "purity", "corrected signal fraction in A", ft.BLUE, draw_eq4_purity_inline, 25),
    ]
    row_h, gap = 66, 13
    yy = y0 + ribbon_h + 18
    for args in rows:
        draw_equation_row(draw, (x0, yy, x1, yy + row_h), *args)
        yy += row_h + gap


def draw_equation_box_large_rows(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = draw_equation_box_frame(draw, box, "")
    def_y = y0
    ribbon_h = 78
    draw_definition_ribbon_large_rows(draw, (x0, def_y, x1, def_y + ribbon_h))
    rows = [
        ("1", "Raw estimate", "candidate count minus ABCD background estimate", ft.PHOTON_DARK, draw_eq2_abcd_inline, 28),
        ("2", "Leakage correction", "sideband count after removing true-photon leakage", ft.TEAL, draw_eq3_sideband_inline, 28),
        ("3", "Purity", "raw sidebands; blue uses leakage-corrected sidebands", ft.BLUE, draw_purity_comparison_inline, 24),
    ]
    yy = def_y + ribbon_h + 10
    for args in rows:
        draw_equation_row(draw, (x0, yy, x1, yy + 76), *args)
        yy += 84


def draw_equation_box_three_columns(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = draw_equation_box_frame(draw, box, "Three columns: each equation gets one job, one formula, and one variable cue.")
    col_gap = 20
    col_w = (x1 - x0 - 2 * col_gap) // 3
    columns = [
        ("Eq 2", "Estimate A background", "Use B, C, D sidebands before leakage correction.", ft.PHOTON_DARK, draw_eq2_abcd_inline, 19, "A is selected; B/C/D are controls."),
        ("Eq 3", "Clean sidebands", "Subtract true-photon leakage from each sideband.", ft.TEAL, draw_eq3_sideband_inline, 19, "X means B, C, or D."),
        ("Eq 4", "Apply purity", "Use corrected signal fraction in the yield.", ft.BLUE, draw_eq4_purity_inline, 23, "P is applied to the Region A yield."),
    ]
    for idx, (badge, title, note, color, renderer, size, cue) in enumerate(columns):
        cx0 = x0 + idx * (col_w + col_gap)
        cx1 = cx0 + col_w
        draw.rounded_rectangle((cx0, y0, cx1, y1), radius=10, fill=(255, 255, 255, 255), outline=(*color, 170), width=2)
        draw.rounded_rectangle((cx0, y0, cx0 + 14, y1), radius=4, fill=(*color, 255))
        draw.text((cx0 + 30, y0 + 24), badge, font=ft.font(ft.TIMES_BOLD, 31), fill=color)
        ft.draw_wrapped(draw, title, (cx0 + 30, y0 + 66), col_w - 60, ft.font(ft.TIMES_BOLD, 24), fill=ft.INK, line_gap=2)
        ft.draw_wrapped(draw, note, (cx0 + 30, y0 + 126), col_w - 60, ft.font(ft.TIMES_ITALIC, 18), fill=ft.MUTED, line_gap=3)
        renderer(draw, cx0 + 30, y0 + 204, size=size)
        draw.rounded_rectangle((cx0 + 26, y1 - 72, cx1 - 26, y1 - 22), radius=8, fill=(248, 251, 253, 255), outline=(*color, 120), width=2)
        ft.draw_wrapped(draw, cue, (cx0 + 44, y1 - 60), col_w - 88, ft.font(ft.TIMES_BOLD, 17), fill=ft.BLUE, line_gap=1)


def draw_equation_box_symbol_left(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = draw_equation_box_frame(draw, box, "Definitions stay fixed on the left; the right side is only the calculation chain.")
    left = (x0, y0, x0 + 438, y1)
    right = (left[2] + 24, y0, x1, y1)
    draw_symbol_key(draw, left)
    rows = [
        ("Eq 2", "raw estimate", "ABCD background estimate", ft.PHOTON_DARK, draw_eq2_abcd_inline, 24),
        ("Eq 3", "leakage removal", "clean B/C/D first", ft.TEAL, draw_eq3_sideband_inline, 23),
        ("Eq 4", "purity", "signal fraction in A", ft.BLUE, draw_eq4_purity_inline, 25),
    ]
    row_h, gap = 78, 17
    yy = right[1]
    for args in rows:
        draw_equation_row(draw, (right[0], yy, right[2], yy + row_h), *args)
        yy += row_h + gap


def draw_equation_box_bookkeeping_flow(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = draw_equation_box_frame(draw, box, "Bookkeeping flow: raw A, corrected sidebands, corrected signal, purity.")
    steps = [
        ("1", "Raw ABCD", ft.PHOTON_DARK, draw_eq2_abcd_inline, 19, "Start from the sideband estimate."),
        ("2", "Clean B/C/D", ft.TEAL, draw_eq3_sideband_inline, 18, "Remove signal leakage from controls."),
        ("3", "Correct A", ft.PHOTON, draw_eq3_signal_inline, 18, "Subtract the corrected background."),
        ("4", "Purity", ft.BLUE, draw_eq4_purity_inline, 20, "Signal fraction applied to A."),
    ]
    step_gap = 18
    step_w = (x1 - x0 - 3 * step_gap) // 4
    yy = y0 + 8
    for idx, (badge, title, color, renderer, size, note) in enumerate(steps):
        sx = x0 + idx * (step_w + step_gap)
        draw.rounded_rectangle((sx, yy, sx + step_w, y1), radius=10, fill=(255, 255, 255, 255), outline=(*color, 170), width=2)
        draw.ellipse((sx + 24, yy + 20, sx + 64, yy + 60), fill=(*color, 255))
        tw = ft.text_box(draw, badge, ft.font(ft.TIMES_BOLD, 22))[0]
        draw.text((sx + 44 - tw / 2, yy + 25), badge, font=ft.font(ft.TIMES_BOLD, 22), fill=(255, 255, 255))
        draw.text((sx + 82, yy + 23), title, font=ft.font(ft.TIMES_BOLD, 23), fill=ft.INK)
        ft.draw_wrapped(draw, note, (sx + 24, yy + 78), step_w - 48, ft.font(ft.TIMES_ITALIC, 17), fill=ft.MUTED, line_gap=2)
        renderer(draw, sx + 24, yy + 150, size=size)
        if idx < len(steps) - 1:
            ft.draw_arrow(draw, (sx + step_w + 2, yy + 146), (sx + step_w + step_gap - 4, yy + 146), fill=(151, 169, 188), width=3)


def draw_refined_method_column(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = shadow_panel(img, box, "1. Define and correct the sidebands")
    x0, y0, x1, y1 = box
    ft.draw_wrapped(
        draw,
        "The key risk is that true photons can leave A and populate the sideband controls.",
        (x0 + 28, y0 + 62),
        x1 - x0 - 56,
        ft.font(ft.TIMES_ITALIC, 20),
        fill=ft.MUTED,
        line_gap=3,
    )
    draw_leakage_map(draw, (x0 + 26, y0 + 116, x1 - 26, y0 + 506))
    draw_equation_method_panel(draw, (x0 + 26, y0 + 534, x1 - 26, y1 - 28))


def draw_refined_purity_panel(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    plot_in_panel(
        img,
        "fig5_purity",
        box,
        "2. Measure purity after leakage correction",
        "blue points and fit: signal-leakage-corrected purity used for the yield",
        inset_top=92,
        inset=18,
    )


def draw_purity_panel_inside_method(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(248, 251, 253, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    label = "Leakage Corrected Purity"
    label_font = ft.font(ft.TIMES_BOLD, 26)
    label_w = ft.text_box(draw, label, label_font)[0]
    draw.text((x0 + (x1 - x0 - label_w) / 2, y0 + 18), label, font=label_font, fill=ft.BLUE)
    plot = trim_white_margins(Image.open(ft.figure_path("fig5_purity")).convert("RGBA"))
    ft.paste_fit(img, plot, (x0 + 28, y0 + 58, x1 - 28, y1 - 24), anchor="center")


def draw_refined_result_panel(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = shadow_panel(img, box, "Build the final result")
    x0, y0, x1, y1 = box
    ft.draw_wrapped(
        draw,
        "Apply purity, efficiency, and unfolding corrections to turn selected candidates into a particle-level yield.",
        (x0 + 30, y0 + 62),
        x1 - x0 - 60,
        ft.font(ft.TIMES_ITALIC, 20),
        fill=ft.MUTED,
        line_gap=3,
    )
    node_w = x1 - x0 - 84
    node_h = 46
    nodes = [
        ("Region A candidates", ft.PHOTON_DARK, x0 + 42, y0 + 126),
        ("apply purity correction", ft.SPHENIX_BLUE, x0 + 42, y0 + 184),
        ("apply efficiency correction", ft.TEAL, x0 + 42, y0 + 242),
        ("unfold to particle level", ft.PHOTON, x0 + 42, y0 + 300),
    ]
    for label, color, nx, ny in nodes:
        bx = (nx, ny, nx + node_w, ny + node_h)
        draw.rounded_rectangle(bx, radius=8, fill=(247, 250, 252, 255), outline=(*color, 220), width=2)
        draw.rounded_rectangle((bx[0], bx[1], bx[0] + 10, bx[3]), radius=4, fill=(*color, 255))
        font = ft.font(ft.TIMES_BOLD, 20)
        tw = ft.text_box(draw, label, font)[0]
        draw.text((bx[0] + 10 + (node_w - 10 - tw) / 2, bx[1] + 12), label, font=font, fill=ft.INK)
    draw.rounded_rectangle((x0 + 34, y1 - 86, x1 - 34, y1 - 24), radius=9, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    target_font = ft.font(ft.TIMES_BOLD, 19)
    body_font = ft.font(ft.TIMES, 19)
    target_x, target_y = x0 + 56, y1 - 70
    draw.text((target_x, target_y), "Target result:", font=target_font, fill=ft.BLUE)
    target_w, _ = ft.text_box(draw, "Target result:", target_font)
    body_x = target_x + target_w + 10
    draw.text((body_x, target_y), "particle-level inclusive isolated prompt photons,", font=body_font, fill=ft.BLUE)
    draw.text((body_x, target_y + 22), "checked against PHENIX and pQCD theory.", font=body_font, fill=ft.BLUE)


def equation_card(img: Image.Image, box: tuple[int, int, int, int], title: str = "Purity equations") -> None:
    draw = shadow_panel(img, box, title)
    x0, y0, x1, y1 = box
    ft.draw_wrapped(
        draw,
        "Eq 2 estimates the background in A; Eq 3 corrects sideband leakage; Eq 4 converts corrected signal count to purity.",
        (x0 + 28, y0 + 60),
        x1 - x0 - 56,
        ft.font(ft.TIMES_ITALIC, 20),
        fill=ft.BLUE,
        line_gap=3,
    )
    sections = [
        ("Eq 2: no leakage", y0 + 118, ft.PHOTON_DARK, c3.draw_eq2, 112),
        ("Eq 3: leakage-corrected", y0 + 252, ft.SPHENIX_BLUE, c3.draw_eq3_compact, 192),
        ("Eq 4: purity", y0 + 474, ft.TEAL, c3.draw_eq4, 116),
    ]
    for heading, sy, color, renderer, height in sections:
        draw.rounded_rectangle((x0 + 24, sy, x1 - 24, sy + height), radius=9, fill=(247, 250, 252, 255), outline=(*color, 180), width=2)
        draw.rounded_rectangle((x0 + 24, sy, x0 + 36, sy + height), radius=4, fill=(*color, 255))
        draw.text((x0 + 50, sy + 14), heading, font=ft.font(ft.TIMES_BOLD, 20), fill=color)
        renderer(draw, x0 + 50, sy + 48)


def equation_ribbon(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = shadow_panel(img, box, "ABCD purity logic")
    x0, y0, x1, y1 = box
    ft.draw_wrapped(
        draw,
        "Region A is selected; B/C/D are sideband controls. Eq 3 removes true-photon leakage before the purity is evaluated.",
        (x0 + 30, y0 + 58),
        x1 - x0 - 60,
        ft.font(ft.TIMES_ITALIC, 21),
        fill=ft.BLUE,
        line_gap=3,
    )
    draw_bcd_mini(draw, (x0 + 30, y0 + 124, x0 + 380, y0 + 256), compact=True)
    c3.draw_eq2(draw, x0 + 438, y0 + 122)
    c3.draw_eq4(draw, x0 + 1040, y0 + 126)
    ft.draw_wrapped(
        draw,
        "B: tight/nonisolated   C: isolated/non-tight   D: nonisolated/non-tight",
        (x0 + 438, y0 + 224),
        x1 - x0 - 468,
        ft.font(ft.TIMES_BOLD, 21),
        fill=ft.MUTED,
        line_gap=2,
    )


def efficiency_mini(img: Image.Image, box: tuple[int, int, int, int], title: str = "Efficiency correction") -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(img, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    title_font = ft.font(ft.TIMES_BOLD, 30)
    title_w = ft.text_box(draw, title, title_font)[0]
    draw.text((x0 + (x1 - x0 - title_w) / 2, y0 + 24), title, font=title_font, fill=ft.INK)
    plot = trim_white_margins(Image.open(ft.figure_path("fig6_efficiencies")).convert("RGBA"))
    ft.paste_fit(img, plot, (x0 + 30, y0 + 78, x1 - 30, y1 - 26), anchor="center")


def final_bridge(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    labels = [
        ("Region A", ft.PHOTON_DARK),
        ("× P", ft.SPHENIX_BLUE),
        ("÷ ε_tot", ft.TEAL),
        ("unfold", ft.PHOTON),
        ("cross section", ft.BLUE),
    ]
    w = (x1 - x0 - 4 * 42) // 5
    x = x0
    for idx, (label, color) in enumerate(labels):
        draw.rounded_rectangle((x, y0, x + w, y1), radius=8, fill=(247, 250, 252, 255), outline=(*color, 230), width=2)
        tw = ft.text_box(draw, label, ft.font(ft.TIMES_BOLD, 22))[0]
        draw.text((x + (w - tw) / 2, y0 + 17), label, font=ft.font(ft.TIMES_BOLD, 22), fill=ft.INK)
        if idx < len(labels) - 1:
            ft.draw_arrow(draw, (x + w + 8, (y0 + y1) // 2), (x + w + 34, (y0 + y1) // 2), fill=(151, 169, 188), width=3)
        x += w + 42


def variant_a() -> tuple[Path, Path]:
    img = base()
    plot_in_panel(img, "fig5_purity", (132, 302, 1390, 1184), "Purity curve is the main correction", "blue points: signal-leakage-corrected purity entering the yield", inset_top=90, inset=22)
    leakage_card(img, (1444, 302, 2390, 666), "Sidebands and leakage")
    equation_card(img, (1444, 704, 2390, 1294), "Equations 2-4")
    draw = ImageDraw.Draw(img, "RGBA")
    final_bridge(draw, (154, 1208, 1368, 1282))
    ft.draw_hp2026_identity_footer(img)
    png = VARIANT_DIR / "hp2026_slide13_purity_variant_A_purity_hero.png"
    img.convert("RGB").save(png)
    script = save_script("hp2026_slide13_purity_variant_A_purity_hero", "A - Purity Hero", SCRIPT_BODY)
    return png, script


def variant_b() -> tuple[Path, Path]:
    img = base()
    equation_ribbon(img, (132, 288, 2390, 560))
    plot_in_panel(img, "fig5_purity", (132, 604, 1366, 1294), "Purity after leakage correction", "the fitted blue curve corrects the selected yield", inset_top=88, inset=22)
    efficiency_mini(img, (1418, 604, 2390, 902), "Efficiency after selection")
    leakage_compact_card(img, (1418, 936, 2390, 1294), "Leakage definitions")
    ft.draw_hp2026_identity_footer(img)
    png = VARIANT_DIR / "hp2026_slide13_purity_variant_B_top_logic_ribbon.png"
    img.convert("RGB").save(png)
    script = save_script("hp2026_slide13_purity_variant_B_top_logic_ribbon", "B - Top Logic Ribbon", SCRIPT_BODY)
    return png, script


def variant_c() -> tuple[Path, Path]:
    img = base()
    leakage_card(img, (132, 304, 914, 660), "ABCD regions and leakage")
    equation_card(img, (132, 696, 914, 1294), "Equations that correct A")
    plot_in_panel(img, "fig5_purity", (954, 304, 1714, 1294), "Measured purity", "blue curve is used for the yield correction", inset_top=86, inset=18)
    efficiency_mini(img, (1752, 304, 2390, 962), "Efficiency")
    draw = shadow_panel(img, (1752, 1000, 2390, 1294), "From selected to result")
    ft.draw_wrapped(draw, "Region A yield × purity ÷ total efficiency, then unfold detector response.", (1784, 1062), 574, ft.font(ft.TIMES_BOLD, 27), fill=ft.BLUE, line_gap=6)
    ft.draw_hp2026_identity_footer(img)
    png = VARIANT_DIR / "hp2026_slide13_purity_variant_C_three_column.png"
    img.convert("RGB").save(png)
    script = save_script("hp2026_slide13_purity_variant_C_three_column", "C - Three Column", SCRIPT_BODY)
    return png, script


def variant_c_refined() -> tuple[Path, Path]:
    img = base()
    draw = ImageDraw.Draw(img, "RGBA")
    method_box = (132, 304, 1758, 1294)
    ft.shadow(img, method_box)
    draw.rounded_rectangle(method_box, radius=14, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw_leakage_map(draw, (160, 332, 1110, 794))
    draw_purity_panel_inside_method(img, (1138, 332, 1730, 794))
    efficiency_mini(img, (1792, 304, 2390, 812), "Apply efficiency")
    draw_wide_equation_panel(img, (160, 826, 1730, 1266))
    draw_refined_result_panel(img, (1792, 844, 2390, 1294))
    ft.draw_hp2026_identity_footer(img)
    png = VARIANT_DIR / "hp2026_slide13_purity_variant_C2_leakage_map_refined.png"
    img.convert("RGB").save(png)
    script = save_script("hp2026_slide13_purity_variant_C2_leakage_map_refined", "C2 - Leakage Map Refined", C2_SCRIPT_BODY)
    return png, script


def equation_box_variant(stem: str, title: str, renderer) -> tuple[Path, Path]:
    img = base()
    draw = ImageDraw.Draw(img, "RGBA")
    method_box = (132, 304, 1758, 1294)
    ft.shadow(img, method_box)
    draw.rounded_rectangle(method_box, radius=14, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw_leakage_map(draw, (160, 332, 1110, 794))
    draw_purity_panel_inside_method(img, (1138, 332, 1730, 794))
    efficiency_mini(img, (1792, 304, 2390, 812), "Apply efficiency")
    renderer(img, (160, 826, 1730, 1266))
    draw_refined_result_panel(img, (1792, 844, 2390, 1294))
    ft.draw_hp2026_identity_footer(img)
    EQUATION_VARIANT_DIR.mkdir(parents=True, exist_ok=True)
    png = EQUATION_VARIANT_DIR / f"{stem}.png"
    img.convert("RGB").save(png)
    script = save_script(stem, title, C2_SCRIPT_BODY)
    return png, script


def equation_box_eye_exam_variants() -> list[tuple[Path, Path]]:
    variants = [
        ("hp2026_slide13_eqbox_A_simple_stack", "Eq Box A - Simple Stack", draw_equation_box_simple_stack),
        ("hp2026_slide13_eqbox_B_large_rows", "Eq Box B - Large Rows", draw_equation_box_large_rows),
        ("hp2026_slide13_eqbox_C_three_columns", "Eq Box C - Three Columns", draw_equation_box_three_columns),
        ("hp2026_slide13_eqbox_D_symbol_left", "Eq Box D - Symbol Left", draw_equation_box_symbol_left),
        ("hp2026_slide13_eqbox_E_bookkeeping_flow", "Eq Box E - Bookkeeping Flow", draw_equation_box_bookkeeping_flow),
    ]
    return [equation_box_variant(stem, title, renderer) for stem, title, renderer in variants]


def variant_d() -> tuple[Path, Path]:
    img = base()
    draw = ImageDraw.Draw(img, "RGBA")
    draw_bcd_mini(draw, (1768, 92, 2058, 176), compact=True)
    plot_in_panel(img, "fig5_purity", (132, 304, 2390, 906), "Purity is the measured signal fraction in Region A", "black: no signal-leakage correction; blue: corrected purity used for the yield", inset_top=84, inset=20)
    equation_ribbon(img, (132, 940, 1542, 1294))
    efficiency_mini(img, (1582, 940, 2390, 1294), "Efficiency is the next correction")
    ft.draw_hp2026_identity_footer(img)
    png = VARIANT_DIR / "hp2026_slide13_purity_variant_D_plot_pair_ribbon.png"
    img.convert("RGB").save(png)
    script = save_script("hp2026_slide13_purity_variant_D_plot_pair_ribbon", "D - Plot With Logic Ribbon", SCRIPT_BODY)
    return png, script


def variant_e() -> tuple[Path, Path]:
    img = base()
    plot_in_panel(img, "fig5_purity", (132, 304, 1220, 1294), "1. Measure purity", "sidebands plus leakage correction give the signal fraction", inset_top=88, inset=18)
    equation_card(img, (1260, 304, 1904, 968), "2. Correct the sidebands")
    leakage_compact_card(img, (1938, 304, 2390, 968), "3. Interpret leakage")
    efficiency_mini(img, (1260, 1002, 1904, 1294), "4. Apply efficiency")
    draw = shadow_panel(img, (1938, 1002, 2390, 1294), "5. Build the yield")
    ft.draw_wrapped(draw, "Region A yield × Purity ÷ ε_tot, unfolded to particle level.", (1968, 1064), 392, ft.font(ft.TIMES_BOLD, 25), fill=ft.BLUE, line_gap=6)
    ft.draw_hp2026_identity_footer(img)
    png = VARIANT_DIR / "hp2026_slide13_purity_variant_E_numbered_flow.png"
    img.convert("RGB").save(png)
    script = save_script("hp2026_slide13_purity_variant_E_numbered_flow", "E - Numbered Flow", SCRIPT_BODY)
    return png, script


SCRIPT_BODY = """With the ID and isolation regions defined, this slide explains how the selected Region A candidates become a corrected photon yield.

The first point is purity. Region A is the tight-and-isolated selected sample, while B, C, and D are the sideband regions. B is tight but non-isolated, C is isolated but non-tight, and D fails both axes. If true prompt photons leak into those sidebands, the sidebands are not purely background-like, so Eq 3 subtracts that leakage before the sidebands are used to estimate the background in A.

Then Eq 4 is the definition that matters for the yield: purity is the corrected signal count in A divided by the raw candidate count in A. That is why the blue purity curve is the curve that gets carried into the result.

The efficiency plot is the second accounting step. Once the selected yield is corrected for purity, it still has to be corrected for reconstruction, photon ID, and isolation efficiency, and then unfolded to particle level. So the take-home message is that the next cross-section plot is not a selected-candidate spectrum. It is a purity-corrected, efficiency-corrected, unfolded photon yield."""


C2_SCRIPT_BODY = """At this point in the talk, we have already defined tight ID and isolation. Now I want to be very explicit about what the purity correction is protecting us from.

The selected sample is Region A: isolated and tight. The sidebands are B, C, and D. In PYTHIA MC, with truth matching, leakage means truth-matched signal photons that are found in B, C, or D instead of staying in A.

If a true photon leaks into B, it still looks tight, but it fails isolation. If it leaks into C, it is isolated, but it fails the tight ID requirement. If it leaks into D, it fails both axes. So the leakage map is telling us whether true photons are being pushed out by the isolation boundary, the shower-shape boundary, or both.

The important assumption is that the two axes are complementary. Isolation is not used as a BDT training feature, so tight ID and isolation are separate handles in this ABCD construction rather than two copies of the same information.

The bottom box is the bookkeeping version of that statement. The symbol key is there so the notation is not a black box: raw N is the raw count in a region, signal N in A is the corrected signal count in the selected region, f from MC is the simulated truth-matched leakage fraction, and b is the background-only sideband count after that leakage is removed.

Eq 2 is the raw ABCD estimate. It subtracts the background predicted from B, C, and D from the candidates in A. The problem is that this only works cleanly if B, C, and D are background-like. Eq 3 is the correction for that: first measure the leakage fractions in simulation, then subtract the leaked true-photon contribution from B, C, and D before using them as sidebands. Eq 4 is then the quantity that enters the yield: purity is the corrected signal count in A divided by the raw candidate count in A.

So when I show the purity plot, the blue points are the object to keep in mind. They are Eq 4 after the Eq 3 leakage correction. The black points show the same purity construction before that leakage correction. That is why the blue curve is the one that gets carried into the corrected yield. The efficiency plot is the next accounting step: reconstruction, ID, and isolation survival before unfolding to particle level."""


def write_contact_sheet(outputs: list[tuple[Path, Path]]) -> Path:
    thumb_w, thumb_h = 512, 288
    label_h = 54
    sheet = Image.new("RGB", (thumb_w * len(outputs), thumb_h + label_h), "white")
    draw = ImageDraw.Draw(sheet)
    for idx, (png, _) in enumerate(outputs):
        x = idx * thumb_w
        draw.text((x + 12, 14), png.stem.replace("hp2026_slide13_purity_variant_", ""), font=ft.font(ft.TIMES, 20), fill=ft.INK)
        thumb = Image.open(png).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(thumb, (x, label_h))
    sheet.save(CONTACT_SHEET)
    return CONTACT_SHEET


def write_equation_contact_sheet(outputs: list[tuple[Path, Path]]) -> Path:
    thumb_w, thumb_h = 512, 288
    label_h = 58
    sheet = Image.new("RGB", (thumb_w * len(outputs), thumb_h + label_h), "white")
    draw = ImageDraw.Draw(sheet)
    for idx, (png, _) in enumerate(outputs):
        x = idx * thumb_w
        label = png.stem.replace("hp2026_slide13_eqbox_", "")
        draw.text((x + 12, 14), label, font=ft.font(ft.TIMES_BOLD, 20), fill=ft.INK)
        thumb = Image.open(png).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(thumb, (x, label_h))
    sheet.save(EQUATION_CONTACT_SHEET)
    return EQUATION_CONTACT_SHEET


def write_manifest(outputs: list[tuple[Path, Path]], contact_sheet: Path) -> Path:
    data = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "target_live_slide": "HPslides_v1 slide object g3e74ae09204_0_85, displayed as slide 13",
        "source_pdf": str(ft.PAPER.relative_to(ft.ROOT)),
        "paper_figures": {
            "fig5_purity": {
                "page": ft.FIGURES["fig5_purity"].page,
                "crop_box": list(ft.FIGURES["fig5_purity"].box),
                "asset": str(ft.figure_path("fig5_purity").relative_to(ft.ROOT)),
            },
            "fig6_efficiencies": {
                "page": ft.FIGURES["fig6_efficiencies"].page,
                "crop_box": list(ft.FIGURES["fig6_efficiencies"].box),
                "asset": str(ft.figure_path("fig6_efficiencies").relative_to(ft.ROOT)),
            },
        },
        "outputs": [
            {
                "png": str(png.relative_to(ft.ROOT)),
                "speaker_script": str(script.relative_to(ft.ROOT)),
                "size": [ft.W, ft.H],
                "mode": "RGB",
            }
            for png, script in outputs
        ],
        "contact_sheet": str(contact_sheet.relative_to(ft.ROOT)),
        "notes": [
            "Equations are hand-rendered as slide graphics, not screenshot crops.",
            "B/C/D leakage definitions are integrated with leakage interpretation.",
            "No slide numbers or provenance footers are baked into the variants.",
            "Current paper plot crops still show sPHENIX Internal and need public/approved labels before final HP insertion.",
        ],
    }
    MANIFEST.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return MANIFEST


def write_equation_manifest(outputs: list[tuple[Path, Path]], contact_sheet: Path) -> Path:
    data = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "purpose": "Five eye-exam variants that keep Slide 13/C2 fixed except for the Equations used to correct Region A box.",
        "equation_box_variants": [
            {
                "png": str(png.relative_to(ft.ROOT)),
                "speaker_script": str(script.relative_to(ft.ROOT)),
                "size": [ft.W, ft.H],
                "mode": "RGB",
            }
            for png, script in outputs
        ],
        "contact_sheet": str(contact_sheet.relative_to(ft.ROOT)),
        "notes": [
            "Only the bottom-left equation box organization varies across these five outputs.",
            "Equations are hand-rendered with visual superscripts/subscripts rather than caret notation.",
            "The surrounding leakage map, purity plot, efficiency plot, footer, and logo are held fixed for comparison.",
            "No Google Slides mutation was performed.",
        ],
    }
    EQUATION_MANIFEST.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return EQUATION_MANIFEST


def main() -> None:
    VARIANT_DIR.mkdir(parents=True, exist_ok=True)
    ft.prepare_figures()
    outputs = [renderer() for renderer in (variant_a, variant_b, variant_c, variant_c_refined, variant_d, variant_e)]
    equation_outputs = equation_box_eye_exam_variants()
    contact = write_contact_sheet(outputs)
    equation_contact = write_equation_contact_sheet(equation_outputs)
    manifest = write_manifest(outputs, contact)
    equation_manifest = write_equation_manifest(equation_outputs, equation_contact)
    print(manifest)
    print(equation_manifest)
    for png, script in outputs:
        print(png)
        print(script)
    for png, script in equation_outputs:
        print(png)
        print(script)
    print(contact)
    print(equation_contact)


if __name__ == "__main__":
    main()
