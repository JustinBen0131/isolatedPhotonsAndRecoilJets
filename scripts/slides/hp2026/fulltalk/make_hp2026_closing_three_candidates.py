#!/usr/bin/env python3
"""Render the three-slide HP2026 closing sequence after tight photon ID.

The sequence compresses the remaining PPG12 story into:
  9. isolation + sideband geometry,
 10. purity + corrections,
 11. final cross section + baseline payoff.

All main plots are cropped from the current PPG12 paper draft by the shared
full-talk generator. Generated drawings are explanatory slide scaffolding only.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_fulltalk_candidates as ft


OUTPUT = ft.OUTPUT
SCRIPT_DIR = ft.SCRIPT_DIR
CONTACT_SHEET = OUTPUT / "hp2026_slides09_11_closing_three_contact_sheet.png"
MANIFEST = OUTPUT / "hp2026_slides09_11_closing_three_manifest.json"
DERIVED_ASSETS = ft.ASSETS / "closing_three"


def text_w(draw: ImageDraw.ImageDraw, text: str, size: int, font_path: Path = ft.TIMES_BOLD) -> int:
    return ft.text_box(draw, text, ft.font(font_path, size))[0]


def save_script(slide_no: int, title: str, body: str, stem: str) -> Path:
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    path = SCRIPT_DIR / f"{stem}.md"
    path.write_text(f"# HP2026 Slide {slide_no} Script - {title}\n\n{body.strip()}\n", encoding="utf-8")
    return path


def draw_formula_parts(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    parts: list[tuple[str, int, int, Path, tuple[int, int, int]]],
) -> int:
    x, y = xy
    for text, size, dy, font_path, fill in parts:
        fnt = ft.font(font_path, size)
        draw.text((x, y + dy), text, font=fnt, fill=fill)
        x += ft.text_box(draw, text, fnt)[0]
    return x


def draw_plot_card(
    base: Image.Image,
    key: str,
    box: tuple[int, int, int, int],
    title: str,
    subtitle: str | None = None,
    inset_top: int = 104,
    inset: int = 28,
    anchor: str = "center",
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((box[0] + 30, box[1] + 22), title, font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    if subtitle:
        ft.draw_wrapped(
            draw,
            subtitle,
            (box[0] + 30, box[1] + 62),
            box[2] - box[0] - 60,
            ft.font(ft.TIMES_ITALIC, 22),
            fill=ft.MUTED,
            line_gap=3,
        )
    img = Image.open(ft.figure_path(key)).convert("RGBA")
    return ft.paste_fit(base, img, (box[0] + inset, box[1] + inset_top, box[2] - inset, box[3] - inset), anchor=anchor)


def derived_fig8_top_ratio() -> Path:
    DERIVED_ASSETS.mkdir(parents=True, exist_ok=True)
    out = DERIVED_ASSETS / "fig8_cross_section_top_plus_theory_ratio.png"
    if out.exists():
        return out
    img = Image.open(ft.figure_path("fig8_cross_section")).convert("RGB")
    crop = img.crop((0, 0, img.width, 1320))
    crop.save(out, "PNG")
    return out


def draw_derived_plot_card(
    base: Image.Image,
    image_path: Path,
    box: tuple[int, int, int, int],
    title: str,
    subtitle: str,
    inset_top: int = 92,
    inset: int = 26,
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((box[0] + 30, box[1] + 22), title, font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        subtitle,
        (box[0] + 30, box[1] + 62),
        box[2] - box[0] - 60,
        ft.font(ft.TIMES_ITALIC, 22),
        fill=ft.MUTED,
        line_gap=3,
    )
    img = Image.open(image_path).convert("RGBA")
    return ft.paste_fit(base, img, (box[0] + inset, box[1] + inset_top, box[2] - inset, box[3] - inset), anchor="center")


def draw_iso_equation_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 24), "Isolation definition", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    draw_wrapped = ft.draw_wrapped
    draw_wrapped(
        draw,
        "sum nearby topo-cluster transverse energy, then subtract the candidate",
        (x0 + 30, y0 + 66),
        x1 - x0 - 60,
        ft.font(ft.TIMES_ITALIC, 22),
        fill=ft.MUTED,
        line_gap=3,
    )
    eq_box = (x0 + 30, y0 + 128, x1 - 30, y0 + 214)
    draw.rounded_rectangle(eq_box, radius=10, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    fx = eq_box[0] + 26
    fy = eq_box[1] + 20
    draw.text((fx, fy), "E", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    draw.text((fx + 18, fy + 14), "T", font=ft.font(ft.TIMES_BOLD, 18), fill=ft.INK)
    draw.text((fx + 26, fy - 12), "iso,reco", font=ft.font(ft.TIMES_BOLD, 17), fill=ft.INK)
    x = fx + 116
    draw.text((x, fy), "< 0.49 + 0.037 E", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    x += ft.text_box(draw, "< 0.49 + 0.037 E", ft.font(ft.TIMES_BOLD, 31))[0]
    draw.text((x, fy + 14), "T", font=ft.font(ft.TIMES_BOLD, 18), fill=ft.INK)
    x += 30
    draw.text((x, fy), "[GeV]", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    draw.rounded_rectangle((x0 + 30, y0 + 238, x1 - 30, y1 - 28), radius=10, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_wrapped(
        draw,
        "small isolation energy means the photon candidate is experimentally quiet in its neighborhood",
        (x0 + 56, y0 + 262),
        x1 - x0 - 112,
        ft.font(ft.TIMES_BOLD, 24),
        fill=ft.BLUE,
        line_gap=5,
    )


def draw_abcd_logic_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 24), "The two axes form the sidebands", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "tight/non-tight ID and isolated/non-isolated regions estimate the residual background",
        (x0 + 30, y0 + 66),
        x1 - x0 - 60,
        ft.font(ft.TIMES_ITALIC, 22),
        fill=ft.MUTED,
        line_gap=3,
    )
    gx0, gy0, gx1, gy1 = x0 + 48, y0 + 128, x1 - 46, y1 - 42
    row_w = 168
    col_y = gy0 + 64
    cell_x0 = gx0 + row_w
    cell_w = (gx1 - cell_x0) // 2
    cell_h = (gy1 - col_y) // 2
    grid_edge = (145, 112, 112, 230)
    draw.rounded_rectangle((gx0, gy0, gx1, gy1), radius=8, fill=(253, 253, 252, 255), outline=grid_edge, width=3)
    draw.rectangle((gx0, gy0, gx1, col_y), fill=(248, 250, 252, 255), outline=grid_edge, width=2)
    draw.rectangle((gx0, col_y, cell_x0, gy1), fill=(248, 250, 252, 255), outline=grid_edge, width=2)
    draw.line((cell_x0, gy0, cell_x0, gy1), fill=grid_edge, width=2)
    draw.line((cell_x0 + cell_w, gy0, cell_x0 + cell_w, gy1), fill=grid_edge, width=2)
    draw.line((gx0, col_y, gx1, col_y), fill=grid_edge, width=2)
    draw.line((gx0, col_y + cell_h, gx1, col_y + cell_h), fill=grid_edge, width=2)
    draw.text((gx0 + 20, gy0 + 22), "γ ID", font=ft.font(ft.TIMES_BOLD, 25), fill=ft.MUTED)
    headers = [
        (cell_x0, cell_x0 + cell_w, "isolated"),
        (cell_x0 + cell_w, gx1, "non-isolated"),
    ]
    for hx0, hx1, label in headers:
        tw = text_w(draw, label, 25, ft.TIMES_BOLD)
        draw.text((hx0 + (hx1 - hx0 - tw) // 2, gy0 + 20), label, font=ft.font(ft.TIMES_BOLD, 25), fill=ft.INK)
    row_labels = [
        (col_y, col_y + cell_h, "non-tight\nID"),
        (col_y + cell_h, gy1, "tight\nID"),
    ]
    for ry0, ry1, label in row_labels:
        lines = label.splitlines()
        line_h = 28
        start_y = ry0 + (ry1 - ry0 - line_h * len(lines)) // 2
        for i, line in enumerate(lines):
            tw = text_w(draw, line, 25, ft.TIMES_BOLD)
            draw.text((gx0 + (row_w - tw) // 2, start_y + i * line_h), line, font=ft.font(ft.TIMES_BOLD, 25), fill=ft.INK)
    regions = [
        ((cell_x0, col_y, cell_x0 + cell_w, col_y + cell_h), "C", "non-tight,\nisolated", ft.TEAL),
        ((cell_x0 + cell_w, col_y, gx1, col_y + cell_h), "D", "non-tight,\nnon-isolated", ft.MUTED),
        ((cell_x0, col_y + cell_h, cell_x0 + cell_w, gy1), "A", "tight,\nisolated", ft.PHOTON_DARK),
        ((cell_x0 + cell_w, col_y + cell_h, gx1, gy1), "B", "tight,\nnon-isolated", ft.SPHENIX_BLUE),
    ]
    for rbox, letter, label, color in regions:
        rx0, ry0, rx1, ry1 = rbox
        if letter == "A":
            draw.rounded_rectangle((rx0 + 10, ry0 + 10, rx1 - 10, ry1 - 10), radius=8, fill=(255, 248, 229, 255), outline=(*color, 230), width=3)
        tx = rx0 + 26
        ty = ry0 + 28
        draw.text((tx, ty), f"{letter}:", font=ft.font(ft.TIMES_BOLD, 34), fill=color)
        draw.multiline_text((tx + 56, ty + 4), label, font=ft.font(ft.TIMES_BOLD, 28), fill=ft.INK, spacing=5)


def draw_correction_flow(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 32, y0 + 18), "Correction accounting bridge", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    steps = [
        ("Region A\nselected yield", ft.PHOTON_DARK, 322),
        ("× purity\nsignal fraction", ft.SPHENIX_BLUE, 284),
        ("/ efficiency\nselection survival", ft.TEAL, 322),
        ("unfold\nresponse", ft.PHOTON, 250),
        ("particle-level\ncross section", ft.BLUE, 332),
    ]
    arrow_w = 52
    gap = 18
    total_w = sum(w for _, _, w in steps) + (len(steps) - 1) * (arrow_w + gap)
    sx = x0 + (x1 - x0 - total_w) // 2
    sy = y0 + 74
    node_h = 74
    for i, (label, color, w) in enumerate(steps):
        draw.rounded_rectangle((sx, sy, sx + w, sy + node_h), radius=10, fill=(247, 250, 252, 255), outline=(*color, 230), width=3)
        draw.rounded_rectangle((sx, sy, sx + 13, sy + node_h), radius=5, fill=(*color, 255))
        lines = label.splitlines()
        line_font = ft.font(ft.TIMES_BOLD, 24)
        line_h = 27
        start_y = sy + (node_h - line_h * len(lines)) // 2 - 1
        for j, line in enumerate(lines):
            tw = ft.text_box(draw, line, line_font)[0]
            draw.text((sx + 13 + (w - 13 - tw) // 2, start_y + j * line_h), line, font=line_font, fill=ft.INK)
        if i < len(steps) - 1:
            ft.draw_arrow(draw, (sx + w + 12, sy + node_h // 2), (sx + w + arrow_w + 4, sy + node_h // 2), fill=(151, 169, 188), width=4)
        sx += w + arrow_w + gap


def draw_n_var(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    sup: str,
    sub: str,
    size: int = 28,
    fill: tuple[int, int, int] = ft.INK,
) -> int:
    base = ft.font(ft.TIMES_ITALIC, size)
    small = ft.font(ft.TIMES, max(13, int(size * 0.55)))
    draw.text((x, y), "N", font=base, fill=fill)
    bw, _ = ft.text_box(draw, "N", base)
    draw.text((x + bw - 1, y - int(size * 0.34)), sup, font=small, fill=fill)
    draw.text((x + bw - 1, y + int(size * 0.48)), sub, font=small, fill=fill)
    sub_w, _ = ft.text_box(draw, sub, small)
    sup_w, _ = ft.text_box(draw, sup, small)
    return x + bw + max(sub_w, sup_w) + 9


def draw_f_var(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    sup: str,
    size: int = 27,
    fill: tuple[int, int, int] = ft.INK,
) -> int:
    base = ft.font(ft.TIMES_ITALIC, size)
    small = ft.font(ft.TIMES, max(13, int(size * 0.55)))
    draw.text((x, y), "f", font=base, fill=fill)
    bw, _ = ft.text_box(draw, "f", base)
    draw.text((x + bw - 1, y - int(size * 0.34)), sup, font=small, fill=fill)
    sup_w, _ = ft.text_box(draw, sup, small)
    return x + bw + sup_w + 8


def math_text(draw: ImageDraw.ImageDraw, x: int, y: int, text: str, size: int = 27, fill=ft.INK, bold: bool = False) -> int:
    fnt = ft.font(ft.TIMES_BOLD if bold else ft.TIMES, size)
    draw.text((x, y), text, font=fnt, fill=fill)
    return x + ft.text_box(draw, text, fnt)[0]


def draw_eq2(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    size = 25
    x = draw_n_var(draw, x, y, "A", "signal", size=size)
    x = math_text(draw, x, y + 2, " = ", size=size)
    x = draw_n_var(draw, x, y, "A", "raw", size=size)
    x = math_text(draw, x, y + 2, " − ", size=size)
    x = draw_n_var(draw, x, y, "B", "raw", size=size)
    x = math_text(draw, x, y + 2, " × ", size=size)
    draw.text((x, y - 9), "(", font=ft.font(ft.TIMES, 52), fill=ft.INK)
    x += 24
    frac_x = x
    num_end = draw_n_var(draw, frac_x + 10, y - 10, "C", "raw", size=22)
    frac_w = max(70, num_end - frac_x + 8)
    draw.line((frac_x, y + 33, frac_x + frac_w, y + 33), fill=ft.INK, width=2)
    draw_n_var(draw, frac_x + 10, y + 38, "D", "raw", size=22)
    draw.text((frac_x + frac_w + 6, y - 9), ")", font=ft.font(ft.TIMES, 52), fill=ft.INK)


def draw_eq3_compact(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    size = 22
    start_x = x
    x = math_text(draw, x, y, "For X = B, C, D:   ", size=20, fill=ft.MUTED, bold=True)
    x = math_text(draw, x, y, "X", size=size, fill=ft.BLUE, bold=True)
    x = math_text(draw, x, y + 2, "corr", size=15, fill=ft.BLUE)
    x = math_text(draw, x + 8, y, " = ", size=size)
    x = draw_n_var(draw, x, y - 2, "X", "raw", size=size)
    x = math_text(draw, x, y, " − ", size=size)
    x = draw_f_var(draw, x, y - 2, "X,MC", size=size)
    draw_n_var(draw, x, y - 2, "A", "signal", size=size)

    x = math_text(draw, start_x, y + 54, "Then:   ", size=20, fill=ft.MUTED, bold=True)
    x = draw_n_var(draw, x, y + 52, "A", "signal", size=23)
    x = math_text(draw, x, y + 54, " = ", size=23)
    x = draw_n_var(draw, x, y + 52, "A", "raw", size=23)
    x = math_text(draw, x, y + 54, " − B", size=23)
    x = math_text(draw, x, y + 68, "corr", size=15)
    x = math_text(draw, x + 8, y + 54, " × ", size=23)
    draw.text((x, y + 43), "(", font=ft.font(ft.TIMES, 48), fill=ft.INK)
    x += 23
    frac_x = x
    num_font = ft.font(ft.TIMES, 22)
    draw.text((frac_x + 12, y + 44), "C", font=ft.font(ft.TIMES_BOLD, 21), fill=ft.INK)
    draw.text((frac_x + 30, y + 57), "corr", font=ft.font(ft.TIMES, 14), fill=ft.INK)
    draw.line((frac_x, y + 82, frac_x + 72, y + 82), fill=ft.INK, width=2)
    draw.text((frac_x + 12, y + 84), "D", font=ft.font(ft.TIMES_BOLD, 21), fill=ft.INK)
    draw.text((frac_x + 30, y + 97), "corr", font=ft.font(ft.TIMES, 14), fill=ft.INK)
    draw.text((frac_x + 78, y + 43), ")", font=ft.font(ft.TIMES, 48), fill=ft.INK)


def draw_eq4(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    size = 27
    x = math_text(draw, x, y, "Purity(", size=size)
    x = math_text(draw, x, y, "P", size=size, fill=ft.BLUE, bold=True)
    x = math_text(draw, x, y, ") = ", size=size)
    frac_x = x
    num_end = draw_n_var(draw, frac_x + 20, y - 14, "A", "signal", size=24)
    frac_w = max(96, num_end - frac_x + 18)
    draw.line((frac_x, y + 33, frac_x + frac_w, y + 33), fill=ft.INK, width=2)
    draw_n_var(draw, frac_x + 20, y + 36, "A", "raw", size=24)


def draw_purity_equation_stack(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=11, fill=(247, 250, 252, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 24, y0 + 20), "ABCD purity equations", font=ft.font(ft.TIMES_BOLD, 27), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "Eq 2 estimates background in A from B, C, D; Eq 3 subtracts true-photon leakage from the sidebands; Eq 4 turns the corrected signal count into purity.",
        (x0 + 24, y0 + 58),
        x1 - x0 - 48,
        ft.font(ft.TIMES_ITALIC, 19),
        fill=ft.BLUE,
        line_gap=3,
    )

    sections = [
        ("Eq 2: no signal leakage", y0 + 136, ft.PHOTON_DARK, draw_eq2),
        ("Eq 3: leakage-corrected sidebands", y0 + 282, ft.SPHENIX_BLUE, draw_eq3_compact),
        ("Eq 4: purity used for the yield", y0 + 512, ft.TEAL, draw_eq4),
    ]
    for heading, sy, color, renderer in sections:
        draw.rounded_rectangle((x0 + 22, sy - 12, x1 - 22, sy + (126 if "Eq 3" not in heading else 206)), radius=9, fill=(255, 255, 255, 255), outline=(*color, 170), width=2)
        draw.rounded_rectangle((x0 + 22, sy - 12, x0 + 34, sy + (126 if "Eq 3" not in heading else 206)), radius=4, fill=(*color, 255))
        draw.text((x0 + 50, sy + 2), heading, font=ft.font(ft.TIMES_BOLD, 21), fill=color)
        renderer(draw, x0 + 52, sy + 40)


def draw_purity_equation_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 22), "1. Purity from ABCD sidebands", font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "blue points are the signal-leakage-corrected purity values used to correct the selected yield",
        (x0 + 30, y0 + 62),
        x1 - x0 - 60,
        ft.font(ft.TIMES_ITALIC, 21),
        fill=ft.MUTED,
        line_gap=3,
    )
    img = Image.open(ft.figure_path("fig5_purity")).convert("RGBA")
    ft.paste_fit(base, img, (x0 + 30, y0 + 120, x0 + 820, y1 - 34), anchor="center")
    draw_purity_equation_stack(draw, (x0 + 850, y0 + 120, x1 - 30, y1 - 34))


def draw_leakage_meaning_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 28, y0 + 22), "What leakage means", font=ft.font(ft.TIMES_BOLD, 30), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "Signal in B/C/D means a true photon escaped the selected region A.",
        (x0 + 28, y0 + 62),
        x1 - x0 - 56,
        ft.font(ft.TIMES_ITALIC, 20),
        fill=ft.BLUE,
        line_gap=3,
    )
    gx0, gy0 = x0 + 30, y0 + 112
    cell_w, cell_h = 82, 62
    labels = [
        ("A", "iso +\ntight", ft.PHOTON_DARK, (0, 0)),
        ("B", "noniso +\ntight", ft.SPHENIX_BLUE, (1, 0)),
        ("C", "iso +\nnontight", ft.TEAL, (0, 1)),
        ("D", "noniso +\nnontight", ft.MUTED, (1, 1)),
    ]
    for letter, label, color, (cx, cy) in labels:
        bx0 = gx0 + cx * cell_w
        by0 = gy0 + cy * cell_h
        fill = (255, 248, 229) if letter == "A" else (247, 250, 252)
        draw.rounded_rectangle((bx0, by0, bx0 + cell_w, by0 + cell_h), radius=7, fill=(*fill, 255), outline=(*color, 220), width=2)
        draw.text((bx0 + 8, by0 + 8), letter, font=ft.font(ft.TIMES_BOLD, 26), fill=color)
        draw.multiline_text((bx0 + 32, by0 + 12), label, font=ft.font(ft.TIMES_BOLD, 13), fill=ft.INK, spacing=2)

    bullets = [
        ("B", "passes tight ID, fails isolation", ft.SPHENIX_BLUE),
        ("C", "isolated, fails tight ID", ft.TEAL),
        ("D", "fails both axes", ft.MUTED),
    ]
    yy = y0 + 116
    for letter, text, color in bullets:
        lx = x0 + 226
        draw.text((lx, yy), letter, font=ft.font(ft.TIMES_BOLD, 24), fill=color)
        draw.text((lx + 34, yy + 1), text, font=ft.font(ft.TIMES, 21), fill=ft.INK)
        yy += 58
    ft.draw_wrapped(
        draw,
        "Eq 3 removes that leaked signal before B/C/D are used as background controls.",
        (x0 + 28, y1 - 72),
        x1 - x0 - 56,
        ft.font(ft.TIMES_BOLD, 21),
        fill=ft.BLUE,
        line_gap=4,
    )


def draw_compact_correction_bridge(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    pieces = [
        ("Region A yield", ft.PHOTON_DARK, 220),
        ("× Purity", ft.SPHENIX_BLUE, 170),
        ("÷ total efficiency", ft.TEAL, 252),
        ("unfold response", ft.PHOTON, 238),
        ("particle-level cross section", ft.BLUE, 360),
    ]
    arrow_w, gap = 42, 18
    total_w = sum(w for _, _, w in pieces) + (len(pieces) - 1) * (arrow_w + gap)
    x = x0 + (x1 - x0 - total_w) // 2
    y = y0 + 28
    h = 64
    for idx, (label, color, w) in enumerate(pieces):
        draw.rounded_rectangle((x, y, x + w, y + h), radius=9, fill=(247, 250, 252, 255), outline=(*color, 225), width=3)
        draw.rounded_rectangle((x, y, x + 12, y + h), radius=4, fill=(*color, 255))
        tw = ft.text_box(draw, label, ft.font(ft.TIMES_BOLD, 23))[0]
        draw.text((x + 14 + (w - 14 - tw) / 2, y + 20), label, font=ft.font(ft.TIMES_BOLD, 23), fill=ft.INK)
        if idx < len(pieces) - 1:
            ft.draw_arrow(draw, (x + w + 10, y + h // 2), (x + w + arrow_w, y + h // 2), fill=(151, 169, 188), width=4)
        x += w + arrow_w + gap


def draw_result_takeaway_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 24), "Final takeaways", font=ft.font(ft.TIMES_BOLD, 34), fill=ft.INK)
    bullets = [
        ("first sPHENIX p+p isolated prompt-photon baseline", ft.PHOTON_DARK),
        ("BDT ID + isolation + data-driven purity", ft.SPHENIX_BLUE),
        ("agreement with NLO pQCD within uncertainties", ft.TEAL),
        ("foundation for future photon + jet and heavy-ion work", ft.BLUE),
    ]
    y = y0 + 92
    for text, color in bullets:
        draw.ellipse((x0 + 36, y + 8, x0 + 56, y + 28), fill=(*color, 255))
        ft.draw_wrapped(draw, text, (x0 + 76, y), x1 - x0 - 112, ft.font(ft.TIMES_BOLD, 26), fill=ft.INK, line_gap=4)
        y += 68


def slide09() -> tuple[Path, Path]:
    img = ft.base_slide(
        "Isolation gives the second axis for purity",
        "Tight ID selects photon-like showers; isolation tests whether the neighborhood around the candidate is quiet.",
    )
    ft.add_top_right_sphenix_logo(img)
    draw_plot_card(
        img,
        "fig3_isolation",
        (132, 318, 1294, 1192),
        "Isolation-energy shape",
        "tight-ID data, non-tight background-enriched data, and signal MC",
        inset_top=106,
        inset=34,
    )
    draw_iso_equation_card(img, (1350, 318, 2390, 622))
    draw_abcd_logic_card(img, (1350, 662, 2390, 1192))
    ft.draw_hp2026_identity_footer(img)
    png = OUTPUT / "hp2026_slide09_isolation_sideband_axis_closing.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        9,
        "Isolation gives the second axis for purity",
        "Now that the tight photon-ID score has selected photon-like showers, the second axis is isolation. The question becomes: is the region around the photon candidate quiet, or does it look like the candidate is sitting inside nearby jet activity?\n\nOn the left, the isolation-energy distribution shows the key separation. The tight-ID data sit near the signal-MC peak at low isolation energy, while the non-tight data provide a background-enriched shape with a longer high-isolation tail. That is exactly the behavior we need: isolation is not just another cut, it is a second experimental handle that distinguishes prompt-like isolated photons from background-rich nearby activity.\n\nThe working definition is shown on the upper right. The analysis uses the reconstructed isolation energy and applies the current PPG12 threshold, 0.49 plus 0.037 times the candidate transverse energy. In words, small isolation energy means the photon candidate is experimentally quiet in its neighborhood.\n\nThe bottom-right panel is the purity logic. Tight versus non-tight ID gives one axis, and isolated versus non-isolated gives the other. Region A is the selected sample, but the background remaining in A is constrained by the other three regions. So the slide sets up the central measurement idea: after BDT photon ID and isolation, purity can be measured from data rather than assumed from simulation.",
        "hp2026_slide09_isolation_sideband_axis_closing_script",
    )
    return png, script


def slide10() -> tuple[Path, Path]:
    img = ft.base_slide(
        "From selected candidates to a corrected photon yield",
        "The ABCD sidebands measure purity; efficiency and unfolding turn the corrected yield into the physics result.",
    )
    ft.add_top_right_sphenix_logo(img)
    draw_purity_equation_panel(img, (132, 314, 1588, 1294))
    draw_plot_card(
        img,
        "fig6_efficiencies",
        (1628, 314, 2390, 844),
        "2. How much prompt signal survives?",
        "reconstruction, ID, isolation, and combined efficiency are explicit corrections",
        inset_top=94,
        inset=18,
    )
    draw_leakage_meaning_card(img, (1628, 876, 2390, 1294))
    ft.draw_hp2026_identity_footer(img)
    png = OUTPUT / "hp2026_slide10_purity_corrections_closing.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        13,
        "From selected candidates to a corrected photon yield",
        "With the ID and isolation regions defined, the analysis has a selected sample, but it does not yet have a physics result. This slide is the accounting bridge between the candidate sample and the corrected photon yield.\n\nThe left side is the purity story. Region A is the tight and isolated sample. In the ideal ABCD picture, Eq 2 says that the background in A can be estimated from the three sideband regions B, C, and D. But the important refinement is Eq 3. If true prompt photons leak into B, C, or D, those sidebands are not pure background controls anymore, so the analysis subtracts the expected signal leakage before using the sidebands. That is why the blue points in the purity plot are the right curve to carry forward.\n\nThe leakage map on the right is how I would explain this physically. Leakage into B means a true photon passes tight ID but fails isolation, so the isolation boundary pushed real signal into the non-isolated region. Leakage into C means a true isolated photon fails tight ID, so the non-tight sideband contains real signal. Leakage into D means the photon fails both axes. The point is not that these failures dominate; the point is that the correction explicitly accounts for them instead of pretending the sidebands are perfectly background-like.\n\nEq 4 then turns the corrected signal count in A into the purity. After that, the efficiency plot answers the second question: how much true prompt signal survives reconstruction, photon ID, and isolation. The final bridge is simple: start from the selected Region A yield, multiply by purity, divide by the total efficiency, unfold the detector response, and then the next slide can honestly be interpreted as the particle-level cross section.",
        "hp2026_slide10_purity_corrections_closing_script",
    )
    return png, script


def slide11() -> tuple[Path, Path]:
    img = ft.base_slide(
        "sPHENIX p+p isolated photons: baseline established",
        "The corrected spectrum agrees with NLO pQCD within uncertainties and anchors the future hard-probes program.",
    )
    ft.add_top_right_sphenix_logo(img)
    draw_derived_plot_card(
        img,
        derived_fig8_top_ratio(),
        (132, 286, 1596, 1296),
        "Main result: isolated prompt-photon cross section",
        "current PPG12 paper figure: cross section and primary theory/data comparison",
        inset_top=88,
        inset=26,
    )
    draw_plot_card(
        img,
        "fig9_phenix",
        (1650, 286, 2390, 840),
        "RHIC context",
        "comparison with PHENIX at the same collision energy",
        inset_top=92,
        inset=26,
    )
    draw_result_takeaway_card(img, (1650, 884, 2390, 1296))
    ft.draw_hp2026_identity_footer(img)
    png = OUTPUT / "hp2026_slide11_result_baseline_closing.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        11,
        "sPHENIX p+p isolated photons: baseline established",
        "This is the closing result slide. The large plot on the left is the corrected isolated prompt-photon cross section in p+p collisions at 200 GeV. By this point in the talk, the audience has seen how the points were built: compact EMCal showers, BDT photon ID, isolation, data-driven purity, efficiency corrections, and unfolding.\n\nThe main physics statement is that the measured spectrum is consistent with NLO pQCD calculations within the quoted uncertainties. I would keep the narration centered there, because this is the money plot of the talk. It shows that sPHENIX can reconstruct and measure isolated prompt photons in p+p at RHIC with a result that is quantitatively comparable to theory.\n\nThe smaller plot on the right puts the measurement in the RHIC context by comparing to PHENIX at the same collision energy. The point is not to overclaim a dramatic discrepancy. The point is that sPHENIX now has its own p+p isolated-photon baseline, with the analysis machinery needed for future photon and photon-plus-jet work.\n\nSo I would end by saying: this measurement establishes the first sPHENIX p+p isolated prompt-photon baseline, using BDT photon ID, isolation, and data-driven purity. It agrees with perturbative QCD within uncertainties, and it provides the baseline that the future heavy-ion hard-probes photon program will build on.",
        "hp2026_slide11_result_baseline_closing_script",
    )
    return png, script


def write_contact_sheet(outputs: list[tuple[int, Path, Path]]) -> Path:
    thumb_w, thumb_h = 640, 360
    label_h = 56
    sheet = Image.new("RGB", (len(outputs) * thumb_w, thumb_h + label_h), "white")
    draw = ImageDraw.Draw(sheet)
    for i, (slide_no, png, _) in enumerate(outputs):
        x = i * thumb_w
        draw.text((x + 16, 16), f"Slide {slide_no}: {png.stem}", font=ft.font(ft.TIMES, 22), fill=ft.INK)
        thumb = Image.open(png).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(thumb, (x, label_h))
    sheet.save(CONTACT_SHEET, "PNG")
    return CONTACT_SHEET


def write_manifest(outputs: list[tuple[int, Path, Path]], contact_sheet: Path) -> Path:
    figure_keys = ["fig3_isolation", "fig5_purity", "fig6_efficiencies", "fig8_cross_section", "fig9_phenix"]
    embedded = []
    for key in figure_keys:
        spec = ft.FIGURES[key]
        embedded.append(
            {
                "key": key,
                "label": spec.label,
                "source_type": spec.source,
                "source_pdf": str(ft.PAPER.relative_to(ft.ROOT)),
                "pdf_page": spec.page,
                "crop_box_px_at_3000_long_edge": list(spec.box),
                "asset": str(ft.figure_path(key).relative_to(ft.ROOT)),
            }
        )
    embedded.append(
        {
            "source_type": "generated explanatory graphic",
            "description": "Clean isolation definition card, sideband map, correction flow, and result takeaway card.",
        }
    )
    embedded.append(
        {
            "source_type": "paper plot derived crop",
            "description": "Presentation crop from PPG12 paper Fig. 8 showing the cross section and primary theory/data ratio.",
            "asset": str(derived_fig8_top_ratio().relative_to(ft.ROOT)),
            "source_asset": str(ft.figure_path("fig8_cross_section").relative_to(ft.ROOT)),
        }
    )
    data = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "sequence": "HP2026 closing three after tight photon ID",
        "slides": [
            {
                "slide": slide_no,
                "png": str(png.relative_to(ft.ROOT)),
                "speaker_script": str(script.relative_to(ft.ROOT)),
                "size": [ft.W, ft.H],
                "mode": "RGB",
            }
            for slide_no, png, script in outputs
        ],
        "qa_contact_sheet": str(contact_sheet.relative_to(ft.ROOT)),
        "embedded_assets": embedded,
        "notes": [
            "Main physics plots are cropped from usefulDocs/sPHENIX_PPG12_Paper_2026-05-21_current_draft.pdf.",
            "No slide numbers or provenance footers are baked into the PNGs.",
            "No Google Slides mutation was performed.",
        ],
    }
    MANIFEST.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return MANIFEST


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    ft.prepare_figures()
    derived_fig8_top_ratio()
    renderers = [(9, slide09), (10, slide10), (11, slide11)]
    outputs = []
    for slide_no, renderer in renderers:
        png, script = renderer()
        outputs.append((slide_no, png, script))
    contact_sheet = write_contact_sheet(outputs)
    manifest = write_manifest(outputs, contact_sheet)
    print(manifest)
    for _, png, _ in outputs:
        print(png)
    print(contact_sheet)


if __name__ == "__main__":
    main()
