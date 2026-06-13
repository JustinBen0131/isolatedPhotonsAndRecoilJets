#!/usr/bin/env python3
"""Local PNG-only candidate for simplified HP2026 Slide 10.

This renders a review candidate only. It does not mutate Google Slides.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_fulltalk_candidates as full
import make_hp2026_slide08_simplified_core_candidate as slide8
import make_hp2026_slide09_simplified_shoulders_candidate as slide9


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260610-slide10-shape-summary-simplified"
PNG_PATH = OUT_DIR / "slide10_shape_summary_simplified_candidate.png"
SCRIPT_PATH = OUT_DIR / "slide10_shape_summary_simplified_candidate_script.md"
MANIFEST_PATH = OUT_DIR / "manifest.json"

W, H = full.W, full.H
PHOTON_RED = slide8.PHOTON_LABEL_RED
BACKGROUND_BLUE = full.SPHENIX_BLUE
CORE_PURPLE = slide8.GEOM_CORE_PURPLE
TOWER_GREEN = slide8.GEOM_TOWER_GREEN
STRIP_GREEN = TOWER_GREEN
REGION_PURPLE = CORE_PURPLE

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


def draw_formula_run_centered(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    parts: list[tuple[str, int, int, Path]],
    *,
    fill: tuple[int, int, int] = full.INK,
) -> None:
    width = 0
    max_h = 0
    for text, size, dy, font_path in parts:
        fnt = full.font(font_path, size)
        tw, th = full.text_box(draw, text, fnt)
        width += tw
        max_h = max(max_h, th + abs(dy))
    x = box[0] + (box[2] - box[0] - width) / 2
    y = box[1] + (box[3] - box[1] - max_h) / 2 - 2
    full.draw_formula_run(draw, (round(x), round(y)), parts, fill=fill)


def draw_base() -> Image.Image:
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text(tuple(HEADER_SPEC["title_xy"]), "Reading photon-like shower shapes", font=full.font(full.TIMES_BOLD, 86), fill=full.INK)
    draw.line((132, HEADER_SPEC["divider_y"], W - 132, HEADER_SPEC["divider_y"]), fill=(221, 226, 232, 255), width=3)
    full.add_top_right_sphenix_logo_like_slide2(img)
    return img


def panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], accent: tuple[int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 232))


def badge(draw: ImageDraw.ImageDraw, xy: tuple[int, int], label: str, color: tuple[int, int, int]) -> None:
    x, y = xy
    r = 29
    draw.ellipse((x, y, x + 2 * r, y + 2 * r), fill=(255, 255, 255, 255), outline=(*color, 235), width=4)
    draw_centered(draw, label, (x + r, y + r), full.font(full.TIMES_BOLD, 32), color)


def caption_pair(draw: ImageDraw.ImageDraw, left: tuple[int, int, int, int], right: tuple[int, int, int, int], y: int, size: int = 31) -> None:
    cf = full.font(full.TIMES_BOLD, size)
    for text, grid, color in (("photon-like", left, PHOTON_RED), ("background-like", right, BACKGROUND_BLUE)):
        gx0, _, gx1, _ = grid
        tw, _ = full.text_box(draw, text, cf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y), text, font=cf, fill=color)


def draw_slide10_width_labels(
    draw: ImageDraw.ImageDraw,
    grid: tuple[int, int, int, int],
    cell: int,
    *,
    long: bool,
    eta_side: str,
) -> None:
    """Compact w_eta / w_phi labels for the narrower Slide 10 shoulder card."""
    gx0, gy0, gx1, gy1 = grid
    center_x = (gx0 + gx1) / 2
    if long:
        x_a, x_b = gx0 + 42, gx1 - 42
        y_a, y_b = gy0 + 34, gy1 - 54
    else:
        x_a, x_b = center_x - 72, center_x + 72
        y_a, y_b = gy0 + 72, gy1 - 82

    ink = full.INK
    label_size = 27
    slide9.draw_double_arrow(draw, (x_a, gy0 + 18), (x_b, gy0 + 18), fill=ink, width=4, head=13)
    tw = slide9.w_symbol_width(draw, "φ", label_size)
    slide9.draw_w_symbol(draw, (round((x_a + x_b - tw) / 2), gy0 - 28), "φ", label_size, fill=ink)

    if eta_side == "right":
        arrow_x = gx1 + 30
        label_x = gx1 + 42
    else:
        arrow_x = gx0 - 30
        label_x = gx0 - 76
    slide9.draw_double_arrow(draw, (arrow_x, y_a), (arrow_x, y_b), fill=ink, width=4, head=13)
    th = full.text_box(draw, "w", full.font(full.TIMES_BOLD, label_size))[1]
    slide9.draw_w_symbol(draw, (round(label_x), round((y_a + y_b - th) / 2)), "η", label_size, fill=ink)


def draw_core_card(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    panel(draw, box, full.PHOTON_DARK)
    badge(draw, (x0 + 38, y0 + 32), "1", full.PHOTON_DARK)
    title = "Core compact?"
    tf = full.font(full.TIMES_BOLD, 40)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 14, y0 + 42), title, font=tf, fill=full.INK)

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
    cell, gap = 48, 50
    pair_w = 2 * 5 * cell + gap
    gx = x0 + (x1 - x0 - pair_w) // 2 + 10
    gy = y0 + 152
    left = full.draw_small_tower_grid(draw, (gx, gy), cell, values, highlight=(2, 2, 3, 3), highlight_color=TOWER_GREEN, secondary=(1, 1, 4, 4), secondary_color=CORE_PURPLE)
    right = full.draw_small_tower_grid(draw, (gx + 5 * cell + gap, gy), cell, diffuse, highlight=(2, 2, 3, 3), highlight_color=TOWER_GREEN, secondary=(1, 1, 4, 4), secondary_color=CORE_PURPLE)
    for grid in (left, right):
        gx0, gy0, _, _ = grid
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*CORE_PURPLE, 255), width=5)
        draw.rectangle((gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell), outline=(*TOWER_GREEN, 255), width=5)
    caption_pair(draw, left, right, y0 + 414, size=29)

    definition = (x0 + 48, y1 - 244, x1 - 36, y1 - 168)
    draw_formula_run_centered(
        draw,
        definition,
        [
            ("E", 30, 0, full.TIMES_BOLD),
            ("11", 19, 11, full.TIMES_BOLD),
            (" / E", 30, 0, full.TIMES_BOLD),
            ("33", 19, 11, full.TIMES_BOLD),
            (" = center-tower dominance", 30, 0, full.TIMES),
        ],
        fill=full.INK,
    )

    take = (x0 + 48, y1 - 132, x1 - 36, y1 - 34)
    draw.rounded_rectangle(take, radius=9, fill=(253, 251, 246, 255), outline=(218, 207, 185, 255), width=2)
    draw_centered(draw, "Photon-like: compact center tower.", ((take[0] + take[2]) / 2, (take[1] + take[3]) / 2), full.font(full.TIMES_BOLD, 36), full.INK)


def draw_shoulders_card(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    panel(draw, box, full.SPHENIX_BLUE)
    badge(draw, (x0 + 38, y0 + 32), "2", full.SPHENIX_BLUE)
    title = "Shoulders narrow?"
    tf = full.font(full.TIMES_BOLD, 38)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 16, y0 + 42), title, font=tf, fill=full.INK)

    cell, gap = 48, 50
    pair_w = 2 * 5 * cell + gap
    gx = x0 + (x1 - x0 - pair_w) // 2 + 10
    gy = y0 + 152
    left = full.draw_small_tower_grid(draw, (gx, gy), cell, slide9.shoulder_values(True), highlight=None, secondary=None)
    right = full.draw_small_tower_grid(draw, (gx + 5 * cell + gap, gy), cell, slide9.shoulder_values(False), highlight=None, secondary=None)
    for grid in (left, right):
        slide9.draw_omitted_seed(draw, grid, cell)
    draw_slide10_width_labels(draw, left, cell, long=False, eta_side="left")
    draw_slide10_width_labels(draw, right, cell, long=True, eta_side="right")
    caption_pair(draw, left, right, y0 + 414, size=29)

    definition = (x0 + 48, y1 - 244, x1 - 36, y1 - 168)
    draw_formula_run_centered(
        draw,
        definition,
        [
            ("w", 31, 0, full.TIMES_BOLD),
            ("η", 20, 12, full.TIMES_BOLD),
            (", w", 31, 0, full.TIMES_BOLD),
            ("φ", 20, 12, full.TIMES_BOLD),
            (" = seed-excluded spread", 31, 0, full.TIMES),
        ],
        fill=full.INK,
    )

    take = (x0 + 48, y1 - 132, x1 - 36, y1 - 34)
    draw.rounded_rectangle(take, radius=9, fill=(239, 247, 252, 255), outline=(182, 216, 238, 255), width=2)
    draw_centered(draw, "Photon-like: narrow surrounding energy.", ((take[0] + take[2]) / 2, (take[1] + take[3]) / 2), full.font(full.TIMES_BOLD, 36), full.INK)


def stretch_values(stretched: bool) -> list[list[float]]:
    if not stretched:
        return [
            [0.03, 0.04, 0.05, 0.04, 0.03],
            [0.05, 0.10, 0.18, 0.11, 0.05],
            [0.07, 0.24, 0.72, 0.25, 0.08],
            [0.05, 0.11, 0.20, 0.12, 0.05],
            [0.03, 0.04, 0.05, 0.04, 0.03],
        ]
    return [
        [0.05, 0.08, 0.11, 0.10, 0.07],
        [0.08, 0.18, 0.28, 0.31, 0.20],
        [0.12, 0.32, 0.62, 0.66, 0.44],
        [0.08, 0.18, 0.30, 0.34, 0.24],
        [0.05, 0.08, 0.12, 0.11, 0.07],
    ]


def draw_stretch_card(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    panel(draw, box, full.TEAL)
    badge(draw, (x0 + 38, y0 + 32), "3", full.TEAL)
    title = "Stretched or split?"
    tf = full.font(full.TIMES_BOLD, 42)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 16, y0 + 42), title, font=tf, fill=full.INK)
    bridge = "Narrow strip compared with wider local region"
    bf = full.font(full.TIMES_BOLD, 29)
    btw, _ = full.text_box(draw, bridge, bf)
    draw.text((x0 + (x1 - x0 - btw) / 2 + 10, y0 + 96), bridge, font=bf, fill=full.BLUE)

    cell, gap = 56, 52
    pair_w = 2 * 5 * cell + gap
    gx = x0 + (x1 - x0 - pair_w) // 2 + 10
    gy = y0 + 162
    left = full.draw_small_tower_grid(draw, (gx, gy), cell, stretch_values(False), highlight=None, secondary=None)
    right = full.draw_small_tower_grid(draw, (gx + 5 * cell + gap, gy), cell, stretch_values(True), highlight=None, secondary=None)
    for grid in (left, right):
        gx0, gy0, _, _ = grid
        draw.rectangle((gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell), outline=(*REGION_PURPLE, 255), width=5)
        strip_box = (gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 4 * cell)
        draw.rectangle(strip_box, outline=(*STRIP_GREEN, 255), width=5)
        region_f = full.font(full.TIMES_BOLD, 24)
        region_label = "3×5 local region"
        rtw, rth = full.text_box(draw, region_label, region_f)
        draw.text((gx0 + 2.5 * cell - rtw / 2, gy0 - rth - 4), region_label, font=region_f, fill=REGION_PURPLE)
        strip_f = full.font(full.TIMES_BOLD, 25)
        strip_label = "3×2 strip"
        stw, sth = full.text_box(draw, strip_label, strip_f)
        draw.text(
            (gx0 + 2.5 * cell - stw / 2, (strip_box[1] + strip_box[3] - sth) / 2 + 5),
            strip_label,
            font=strip_f,
            fill=STRIP_GREEN,
        )
    caption_pair(draw, left, right, gy + 5 * cell + 26, size=32)

    definition = (x0 + 48, y1 - 200, x1 - 36, y1 - 134)
    draw_formula_run_centered(
        draw,
        definition,
        [
            ("E", 34, 0, full.TIMES_BOLD),
            ("3×2", 21, 13, full.TIMES_BOLD),
            (" / E", 34, 0, full.TIMES_BOLD),
            ("3×5", 21, 13, full.TIMES_BOLD),
            (" = narrow-strip share", 34, 0, full.TIMES),
        ],
        fill=full.INK,
    )

    take = (x0 + 48, y1 - 132, x1 - 36, y1 - 34)
    draw.rounded_rectangle(take, radius=9, fill=(240, 250, 250, 255), outline=(179, 222, 224, 255), width=2)
    draw_centered(draw, "Photon-like: little elongation or splitting.", ((take[0] + take[2]) / 2, (take[1] + take[3]) / 2), full.font(full.TIMES_BOLD, 36), full.INK)


def draw_formula_chip(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], parts: list[tuple[str, int, int, Path]]) -> None:
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(215, 224, 234, 255), width=2)
    width = 0
    for text, size, _dy, font_path in parts:
        width += full.text_box(draw, text, full.font(font_path, size))[0]
    x = box[0] + (box[2] - box[0] - width) / 2
    y = box[1] + 27
    full.draw_formula_run(draw, (round(x), y), parts, fill=full.INK)


def draw_classical_cut_band(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=13, fill=(249, 251, 253, 255), outline=(217, 226, 235, 255), width=2)
    title = "Classical box cuts"
    tf = full.font(full.TIMES_BOLD, 52)
    draw.text((x0 + 34, y0 + 28), title, font=tf, fill=full.INK)
    subtitle = "If cuts for"
    sf = full.font(full.TIMES, 54)
    sx = x0 + 34
    sy = y0 + 112
    draw.text((sx, sy), subtitle, font=sf, fill=full.INK)
    sx += full.text_box(draw, subtitle, sf)[0] + 18
    colors = [full.PHOTON_DARK, full.SPHENIX_BLUE, full.TEAL]
    for idx, color in enumerate(colors, start=1):
        r = 34
        draw.ellipse((sx, sy - 4, sx + 2 * r, sy + 2 * r - 4), fill=(255, 255, 255, 255), outline=(*color, 235), width=5)
        draw_centered(draw, str(idx), (sx + r, sy + r - 4), full.font(full.TIMES_BOLD, 37), color)
        sx += 2 * r + 16
        if idx < 3:
            plus_f = full.font(full.TIMES_BOLD, 45)
            draw.text((sx, sy + 5), "+", font=plus_f, fill=full.LIGHT_MUTED)
            sx += full.text_box(draw, "+", plus_f)[0] + 16
    phrase = "pass  →  accept candidate"
    draw.text((sx + 10, sy), phrase, font=full.font(full.TIMES_BOLD, 54), fill=full.BLUE)

    divider_x = x0 + 1300
    draw.line((divider_x, y0 + 48, divider_x, y1 - 48), fill=(214, 224, 234, 255), width=3)
    bullets = [
        "Transparent but rigid",
        "Hard to exploit interdependencies",
    ]
    bx = divider_x + 68
    by = y0 + 68
    hf = full.font(full.TIMES_BOLD, 50)
    for head in bullets:
        draw.ellipse((bx, by + 17, bx + 24, by + 41), fill=(*full.INK, 255))
        draw.text((bx + 48, by), head, font=hf, fill=full.INK)
        by += 116


def build_slide() -> Image.Image:
    img = draw_base()
    draw = ImageDraw.Draw(img, "RGBA")
    top = 258
    cards_bottom = 966
    gap = 40
    card_w = (W - 2 * 132 - 2 * gap) // 3
    cards = [
        (132, top, 132 + card_w, cards_bottom),
        (132 + card_w + gap, top, 132 + 2 * card_w + gap, cards_bottom),
        (132 + 2 * (card_w + gap), top, W - 132, cards_bottom),
    ]
    draw_core_card(draw, cards[0])
    draw_shoulders_card(draw, cards[1])
    draw_stretch_card(draw, cards[2])
    draw_classical_cut_band(draw, (132, 998, W - 132, 1282))
    full.draw_hp2026_identity_footer(img)
    return img.convert("RGB")


def write_script() -> None:
    SCRIPT_PATH.write_text(
        """# Slide 10 script - shower-shape checks to box cuts

Putting the three shower-shape checks together, the variables are transparent: compact core, narrow surrounding energy, and little elongation or splitting.

Classically, those handles become rectangular cuts. That is easy to explain, but it treats the thresholds mostly independently. This motivates the next step: a learned classifier can use the same shower-shape information while learning dependencies between variables.
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
                "slide_intent": "Simplified Slide 10: three visual shower-shape checks plus a bottom rectangular-cut bridge toward BDT motivation.",
                "advisor_comment_addressed": "Slides 8-10 need less detail; keep pictures and bottom boxes, remove dense definitions.",
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
