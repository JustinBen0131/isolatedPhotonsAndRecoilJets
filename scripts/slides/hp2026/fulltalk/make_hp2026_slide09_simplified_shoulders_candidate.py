#!/usr/bin/env python3
"""Local PNG-only candidate for simplified HP2026 Slide 9.

This renders a review candidate only. It does not mutate Google Slides.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_fulltalk_candidates as full
import make_hp2026_slide08_simplified_core_candidate as slide8


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260610-slide9-shoulders-simplified"
PNG_PATH = OUT_DIR / "slide09_shoulders_simplified_candidate.png"
SCRIPT_PATH = OUT_DIR / "slide09_shoulders_simplified_candidate_script.md"
MANIFEST_PATH = OUT_DIR / "manifest.json"

W, H = full.W, full.H
PHOTON_RED = slide8.PHOTON_LABEL_RED
GEOM_CORE_PURPLE = slide8.GEOM_CORE_PURPLE
GEOM_TOWER_GREEN = slide8.GEOM_TOWER_GREEN
WIDTH_BLUE = full.SPHENIX_BLUE
WIDTH_RED = (231, 78, 52)
WIDTH_INK = full.INK

HEADER_SPEC = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": None,
    "title_xy": [132, 76],
    "subtitle_xy": None,
    "divider_y": 232,
}


def draw_centered(draw: ImageDraw.ImageDraw, text: str, center: tuple[float, float], font, fill) -> None:
    bbox = draw.textbbox((0, 0), text, font=font)
    draw.text(
        (center[0] - (bbox[0] + bbox[2]) / 2, center[1] - (bbox[1] + bbox[3]) / 2),
        text,
        font=font,
        fill=fill,
    )


def draw_arrow(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    fill: tuple[int, int, int],
    width: int = 5,
    head: int = 16,
) -> None:
    draw.line((start[0], start[1], end[0], end[1]), fill=(*fill, 245), width=width)
    angle = math.atan2(end[1] - start[1], end[0] - start[0])
    pts = [
        end,
        (end[0] - head * math.cos(angle - 0.46), end[1] - head * math.sin(angle - 0.46)),
        (end[0] - head * math.cos(angle + 0.46), end[1] - head * math.sin(angle + 0.46)),
    ]
    draw.polygon(pts, fill=(*fill, 245))


def draw_double_arrow(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    fill: tuple[int, int, int],
    width: int = 5,
    head: int = 16,
) -> None:
    draw_arrow(draw, start, end, fill=fill, width=width, head=head)
    draw_arrow(draw, end, start, fill=fill, width=width, head=head)


def draw_base() -> Image.Image:
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text(tuple(HEADER_SPEC["title_xy"]), "Reading photon-like shower shapes", font=full.font(full.TIMES_BOLD, 86), fill=full.INK)
    y = HEADER_SPEC["divider_y"]
    draw.line((132, y, W - 132, y), fill=(221, 226, 232, 255), width=3)
    full.add_top_right_sphenix_logo_like_slide2(img)
    return img


def panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], accent: tuple[int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 230))


def badge(draw: ImageDraw.ImageDraw, xy: tuple[int, int], label: str, color: tuple[int, int, int]) -> None:
    x, y = xy
    r = 34
    draw.ellipse((x, y, x + 2 * r, y + 2 * r), fill=(255, 255, 255, 255), outline=(*color, 235), width=4)
    draw_centered(draw, label, (x + r, y + r), full.font(full.TIMES_BOLD, 37), color)


def draw_e11_e33_line(draw: ImageDraw.ImageDraw, xy: tuple[int, int], size: int, fill=full.BLUE) -> int:
    x, y = xy
    x = slide8.draw_e11_e33(draw, (x, y), size, fill=fill)
    eq_f = full.font(full.TIMES_BOLD, size)
    draw.text((x, y), ": ", font=eq_f, fill=full.INK)
    x += full.text_box(draw, ": ", eq_f)[0]
    body = "center vs local 3×3 energy"
    draw.text((x, y), body, font=full.font(full.TIMES, size), fill=full.INK)
    return x + full.text_box(draw, body, full.font(full.TIMES, size))[0]


def draw_w_symbol(draw: ImageDraw.ImageDraw, xy: tuple[int, int], subscript: str, size: int, fill=full.INK) -> int:
    sub = max(17, round(size * 0.62))
    return full.draw_formula_run(
        draw,
        xy,
        [
            ("w", size, 0, full.TIMES_BOLD),
            (subscript, sub, round(size * 0.36), full.TIMES_BOLD),
        ],
        fill=fill,
    )


def w_symbol_width(draw: ImageDraw.ImageDraw, subscript: str, size: int) -> int:
    sub = max(17, round(size * 0.62))
    return (
        full.text_box(draw, "w", full.font(full.TIMES_BOLD, size))[0]
        + full.text_box(draw, subscript, full.font(full.TIMES_BOLD, sub))[0]
    )


def draw_mini_core_summary(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    panel(draw, box, full.PHOTON_DARK)
    badge(draw, (x0 + 44, y0 + 38), "1", full.PHOTON_DARK)

    title = "Is the core compact?"
    tf = full.font(full.TIMES_BOLD, 45)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 22, y0 + 54), title, font=tf, fill=full.INK)

    values = [
        [0.04, 0.06, 0.08, 0.06, 0.04],
        [0.06, 0.18, 0.30, 0.18, 0.06],
        [0.08, 0.36, 0.95, 0.42, 0.08],
        [0.06, 0.20, 0.38, 0.19, 0.06],
        [0.04, 0.06, 0.08, 0.06, 0.04],
    ]
    diffuse = [
        [0.13, 0.20, 0.28, 0.20, 0.13],
        [0.22, 0.37, 0.52, 0.40, 0.22],
        [0.30, 0.56, 0.78, 0.59, 0.32],
        [0.22, 0.40, 0.54, 0.43, 0.23],
        [0.13, 0.20, 0.30, 0.20, 0.13],
    ]
    cell = 62
    gap = 54
    pair_w = 2 * 5 * cell + gap
    gx = x0 + (x1 - x0 - pair_w) // 2 + 8
    gy = y0 + 166
    left = full.draw_small_tower_grid(draw, (gx, gy), cell, values, highlight=(2, 2, 3, 3), highlight_color=GEOM_TOWER_GREEN, secondary=(1, 1, 4, 4), secondary_color=GEOM_CORE_PURPLE)
    right = full.draw_small_tower_grid(draw, (gx + 5 * cell + gap, gy), cell, diffuse, highlight=(2, 2, 3, 3), highlight_color=GEOM_TOWER_GREEN, secondary=(1, 1, 4, 4), secondary_color=GEOM_CORE_PURPLE)
    for grid in (left, right):
        gx0, gy0, _, _ = grid
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*GEOM_CORE_PURPLE, 255), width=6)
        draw.rectangle((gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell), outline=(*GEOM_TOWER_GREEN, 255), width=6)
        core_f = full.font(full.TIMES_BOLD, 20)
        core_label = "3×3 local core"
        ctw, cth = full.text_box(draw, core_label, core_f)
        draw.text((gx0 + 2.5 * cell - ctw / 2, gy0 + cell - cth - 7), core_label, font=core_f, fill=GEOM_CORE_PURPLE)
        tower_f = full.font(full.TIMES_BOLD, 20)
        tower_label = "1×1 center tower"
        ttw, _ = full.text_box(draw, tower_label, tower_f)
        draw.text((gx0 + 2.5 * cell - ttw / 2, gy0 + 3 * cell + 9), tower_label, font=tower_f, fill=GEOM_TOWER_GREEN)
    cf = full.font(full.TIMES_BOLD, 32)
    for text, grid, color in (("photon-like", left, PHOTON_RED), ("background-like", right, full.SPHENIX_BLUE)):
        gx0, _, gx1, gy1 = grid
        tw, _ = full.text_box(draw, text, cf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, gy1 + 22), text, font=cf, fill=color)

    eq_box = (x0 + 42, y0 + 610, x1 - 24, y0 + 718)
    eq_size = 40
    eq_w = slide8.e11_e33_width(draw, eq_size) + full.text_box(draw, ": center vs local 3×3 energy", full.font(full.TIMES, eq_size))[0]
    draw_e11_e33_line(draw, (round(eq_box[0] + (eq_box[2] - eq_box[0] - eq_w) / 2), eq_box[1] + 28), eq_size)

    take = (x0 + 64, y1 - 190, x1 - 44, y1 - 64)
    draw.rounded_rectangle(take, radius=11, fill=(253, 251, 246, 255), outline=(218, 207, 185, 255), width=2)
    prefix = "Photon-like:"
    pf = full.font(full.TIMES_BOLD, 42)
    bf = full.font(full.TIMES, 42)
    body = " compact center tower."
    total_w = full.text_box(draw, prefix, pf)[0] + full.text_box(draw, body, bf)[0]
    tx = take[0] + (take[2] - take[0] - total_w) / 2
    ty = take[1] + (take[3] - take[1] - full.text_box(draw, prefix, pf)[1]) / 2 - 1
    draw.text((tx, ty), prefix, font=pf, fill=full.INK)
    tx += full.text_box(draw, prefix, pf)[0]
    draw.text((tx, ty), body, font=bf, fill=full.INK)


def shoulder_values(narrow: bool) -> list[list[float]]:
    if narrow:
        return [
            [0.03, 0.04, 0.05, 0.04, 0.03],
            [0.04, 0.08, 0.16, 0.08, 0.04],
            [0.05, 0.16, 0.00, 0.17, 0.05],
            [0.04, 0.08, 0.15, 0.08, 0.04],
            [0.03, 0.04, 0.05, 0.04, 0.03],
        ]
    return [
        [0.12, 0.18, 0.24, 0.18, 0.12],
        [0.18, 0.34, 0.46, 0.36, 0.18],
        [0.25, 0.48, 0.00, 0.50, 0.28],
        [0.18, 0.36, 0.46, 0.38, 0.20],
        [0.12, 0.18, 0.26, 0.18, 0.12],
    ]


def draw_omitted_seed(draw: ImageDraw.ImageDraw, grid: tuple[int, int, int, int], cell: int) -> None:
    gx0, gy0, _, _ = grid
    c0 = gx0 + 2 * cell
    r0 = gy0 + 2 * cell
    draw.rectangle((c0, r0, c0 + cell, r0 + cell), fill=(255, 255, 255, 255), outline=(188, 198, 208, 255), width=3)
    draw.line((c0 + 12, r0 + 12, c0 + cell - 12, r0 + cell - 12), fill=(108, 116, 126, 230), width=4)
    draw.line((c0 + cell - 12, r0 + 12, c0 + 12, r0 + cell - 12), fill=(108, 116, 126, 230), width=4)
    draw.ellipse((c0 + cell / 2 - 9, r0 + cell / 2 - 9, c0 + cell / 2 + 9, r0 + cell / 2 + 9), fill=(0, 0, 0, 255))


def draw_width_labels(draw: ImageDraw.ImageDraw, grid: tuple[int, int, int, int], cell: int, *, long: bool) -> None:
    gx0, gy0, gx1, gy1 = grid
    center_x = (gx0 + gx1) / 2
    if long:
        x_a, x_b = gx0 + 46, gx1 - 46
        y_a, y_b = gy0 + 34, gy1 - 54
    else:
        x_a, x_b = center_x - 82, center_x + 82
        y_a, y_b = gy0 + 78, gy1 - 94
    draw_double_arrow(draw, (x_a, gy0 + 20), (x_b, gy0 + 20), fill=WIDTH_INK, width=5, head=15)
    label_size = 33
    tw = w_symbol_width(draw, "φ", label_size)
    draw_w_symbol(draw, (round((x_a + x_b - tw) / 2), gy0 - 30), "φ", label_size, fill=WIDTH_INK)
    draw_double_arrow(draw, (gx0 - 32, y_a), (gx0 - 32, y_b), fill=WIDTH_INK, width=5, head=15)
    tw = w_symbol_width(draw, "η", label_size)
    th = full.text_box(draw, "w", full.font(full.TIMES_BOLD, label_size))[1]
    draw_w_symbol(draw, (round(gx0 - 78), round((y_a + y_b - th) / 2)), "η", label_size, fill=WIDTH_INK)


def draw_shoulders_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    panel(draw, box, full.SPHENIX_BLUE)
    badge(draw, (x0 + 44, y0 + 38), "2", full.SPHENIX_BLUE)

    title = "Are the shoulders narrow?"
    tf = full.font(full.TIMES_BOLD, 58)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2, y0 + 48), title, font=tf, fill=full.INK)
    subtitle = "Energy map with the center tower omitted"
    sf = full.font(full.TIMES_BOLD, 42)
    sw, _ = full.text_box(draw, subtitle, sf)
    draw.text((x0 + (x1 - x0 - sw) / 2, y0 + 122), subtitle, font=sf, fill=full.BLUE)

    cell = 88
    gap = 190
    pair_w = 2 * 5 * cell + gap
    gx = x0 + (x1 - x0 - pair_w) // 2 + 44
    gy = y0 + 224
    left = full.draw_small_tower_grid(draw, (gx, gy), cell, shoulder_values(True), highlight=None, secondary=None)
    right = full.draw_small_tower_grid(draw, (gx + 5 * cell + gap, gy), cell, shoulder_values(False), highlight=None, secondary=None)
    for grid, is_long in ((left, False), (right, True)):
        draw_omitted_seed(draw, grid, cell)
        draw_width_labels(draw, grid, cell, long=is_long)

    cf = full.font(full.TIMES_BOLD, 43)
    for text, grid, color in (("photon-like", left, PHOTON_RED), ("background-like", right, full.SPHENIX_BLUE)):
        gx0, _, gx1, gy1 = grid
        tw, _ = full.text_box(draw, text, cf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, gy1 + 34), text, font=cf, fill=color)

    definition_y = y0 + 805
    definition_size = 45
    f2 = full.font(full.TIMES, definition_size)
    total = (
        w_symbol_width(draw, "η", definition_size)
        + full.text_box(draw, ", ", f2)[0]
        + w_symbol_width(draw, "φ", definition_size)
        + full.text_box(draw, " = seed-excluded spread of surrounding energy", f2)[0]
    )
    tx = x0 + (x1 - x0 - total) / 2 + 18
    tx = draw_w_symbol(draw, (round(tx), definition_y), "η", definition_size, fill=WIDTH_INK)
    draw.text((tx, definition_y), ", ", font=f2, fill=full.INK)
    tx += full.text_box(draw, ", ", f2)[0]
    tx = draw_w_symbol(draw, (round(tx), definition_y), "φ", definition_size, fill=WIDTH_INK)
    draw.text((tx, definition_y), " = seed-excluded spread of surrounding energy", font=f2, fill=full.INK)

    take = (x0 + 110, y0 + 888, x1 - 78, y1 - 38)
    draw.rounded_rectangle(take, radius=11, fill=(239, 247, 252, 255), outline=(182, 216, 238, 255), width=3)
    prefix = "Photon-like:"
    pf = full.font(full.TIMES_BOLD, 49)
    bf = full.font(full.TIMES, 49)
    body = " narrow shoulders around the compact core."
    total_w = full.text_box(draw, prefix, pf)[0] + full.text_box(draw, body, bf)[0]
    tx = take[0] + (take[2] - take[0] - total_w) / 2
    ty = take[1] + (take[3] - take[1] - full.text_box(draw, prefix, pf)[1]) / 2 - 2
    draw.text((tx, ty), prefix, font=pf, fill=full.INK)
    tx += full.text_box(draw, prefix, pf)[0]
    draw.text((tx, ty), body, font=bf, fill=full.INK)


def build_slide() -> Image.Image:
    img = draw_base()
    draw = ImageDraw.Draw(img, "RGBA")
    left = (132, 258, 870, 1282)
    right = (912, 258, W - 132, 1282)
    draw_mini_core_summary(draw, left)
    draw_shoulders_panel(draw, right)
    full.draw_hp2026_identity_footer(img)
    return img.convert("RGB")


def write_script() -> None:
    SCRIPT_PATH.write_text(
        """# Slide 9 script - shoulder width

After the compact-core check, the next question is whether the surrounding energy stays narrow once the seed tower is omitted.

The photon-like example keeps the seed-excluded energy close to the core in both eta and phi. The background-like example spreads that surrounding energy more broadly, which is what the shoulder-width variables capture.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img = build_slide()
    img.save(PNG_PATH, "PNG")
    write_script()
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "google_slides_mutation": False,
                "hp2026_main_header": HEADER_SPEC,
                "source_baseline": "scripts/slides/hp2026/fulltalk/make_hp2026_shower_no_blur_variants.py",
                "output_png": str(PNG_PATH.relative_to(ROOT)),
                "companion_script": str(SCRIPT_PATH.relative_to(ROOT)),
                "slide_intent": "Simplified Slide 9 local PNG candidate: keep the compact-core context, then explain seed-excluded shoulder width with pictures and one bottom takeaway.",
                "advisor_comment_addressed": "Slides 8-10 need less detail; keep pictures and bottom boxes, remove dense definitions.",
                "removed_from_old_slide9": [
                    "slide subtitle under the title",
                    "core fraction E2x2/Ecluster detail in the context card",
                    "dense shoulder-width prose rows",
                    "old small-card typography and cramped definitions",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG_PATH)
    print(SCRIPT_PATH)


if __name__ == "__main__":
    main()
