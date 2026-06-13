#!/usr/bin/env python3
"""Local BDT yes/no build sequence for HP2026.

This generator prototypes Slides 11-13 as local PNGs only.  It imports the
current HP2026 full-talk drawing primitives so typography, footer, logo, and
deck chrome match the existing Hard Probes slide system.  It does not mutate
Google Slides.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFilter

import make_hp2026_fulltalk_candidates as full


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260609-slides11-13-bdt-yesno-sequence"
FUNNEL_ASSET = OUT_DIR / "assets/32_funelo_ink_transparent.png"
W, H = full.W, full.H

TITLE = "From fixed cuts to a BDT score"
SIMPLIFIED_TITLE = "From fixed cuts to a photon-ID score"
SIMPLIFIED_SUBTITLE = "The BDT uses the same shower-shape handles, but learns how to combine them."
COMPRESSED_SLIDE_A_TITLE = "Shower shape gives interpretable photon-ID handles"
COMPRESSED_SLIDE_A_SUBTITLE = "Representative geometry handles connect the EMCal pictures to photon ID."
COMPRESSED_SLIDE_A_V2_TITLE = "Interpretable shower-shape handles for photon ID"
COMPRESSED_SLIDE_B_TITLE = "From cut logic to a BDT score"
COMPRESSED_SLIDE_B_SUBTITLE = "The BDT asks similar shower-shape questions, but learns how to combine them."
SUBTITLES = [
    "Fixed yes/no gates are transparent, but each variable is thresholded separately.",
    "The same shower-shape handles feed one learned conditional tree.",
    "A boosted ensemble combines many learned trees into one photon-ID ranking score.",
]
HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": None,
    "title_xy": [132, 76],
    "subtitle_xy": None,
    "divider_y": 232,
}

SCRIPT_MD = """# Slides 11-13 BDT yes/no sequence

## Slide 11 - Fixed rectangular cuts

After defining the shower-shape variables, the most transparent classical approach is to ask a set of fixed yes/no questions.

For example, we can require the core to be compact, the surrounding widths to be narrow, and the deposit not to look stretched or split. Each row here is a rectangular cut: the threshold is fixed, and the candidate either passes or fails.

That is clean and interpretable, but it is rigid. It does not learn the interdependencies between variables; each variable is thresholded separately.

## Slide 12 - One learned decision tree

The BDT starts from the same familiar shower-shape handles, but instead of applying them as independent fixed cuts, it learns a conditional path.

In one tree, the first question can be something like whether the shoulders are narrow. Depending on that answer, the next question can change: one branch may ask whether the core is compact, while another asks whether the deposit is stretched or split.

So the tree keeps the yes/no language, but learns the order of the questions from the data.

## Slide 13 - Boosted decision tree ensemble

The boosted decision tree repeats that idea many times. Each shallow tree gives a weighted vote, and the ensemble combines those votes into one score.

The score is a ranking axis: candidates closer to one are more photon-like, and candidates closer to zero are more background-like.

Importantly, the BDT score is not itself the final purity. The purity is calibrated later with the sideband method.
"""

SIMPLIFIED_SCRIPT_MD = """# Slide 11 - From fixed cuts to a photon-ID score

I want to make clear what the BDT is doing without turning this into a machine-learning lecture.

The inputs are the same physical shower-shape handles we just defined: whether the core is compact, whether the shoulders are narrow, and whether the cluster looks split or elongated. A rectangular-cut analysis applies those questions one by one with fixed thresholds. That is transparent, but rigid.

A decision tree keeps the same yes/no language, but learns a conditional order. Depending on the first answer, the next shower-shape question can change.

The boosted decision tree repeats that idea many times and combines the shallow trees into one photon-ID score. Low score is more background-like; high score is more photon-like. The score is the identification axis, not the final purity; the purity is calibrated later with sidebands.
"""

COMPRESSED_SCRIPT_MD = """# Two-slide compression replacing detailed Slides 8-13

## Slide A - Shower shape gives interpretable photon-ID handles

Before introducing the BDT, I want to show the kind of information it uses. Shower shape gives us physically interpretable handles. A photon-like EMCal deposit tends to have a compact core, narrow shoulders, and no split or elongated local structure.

These are representative examples, not the full list of inputs. The important point is that photon ID is built from geometry questions about the shower, not arbitrary numbers. A classical rectangular-cut analysis answers questions like these with fixed thresholds.

## Slide B - From cut logic to a BDT score

Classical cuts apply fixed thresholds to these shower-shape handles independently. A decision tree asks the same kind of questions, but the next question can depend on the previous answer.

The BDT repeats that logic with many shallow trees and combines their contributions into one photon-ID score. Low score is more background-like; high score is more photon-like. The score is the identification axis that feeds the later purity correction.
"""

ACCENTS = {
    "cuts": full.PHOTON_DARK,
    "tree": full.SPHENIX_BLUE,
    "boost": full.TEAL,
    "fail": (198, 82, 64),
}


def centered_text(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    text: str,
    fnt: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int],
    *,
    stroke_width: int = 0,
    stroke_fill: tuple[int, int, int] | None = None,
) -> None:
    bbox = draw.textbbox((0, 0), text, font=fnt, stroke_width=stroke_width)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    x = (box[0] + box[2] - tw) / 2 - bbox[0]
    y = (box[1] + box[3] - th) / 2 - bbox[1]
    kwargs = {"font": fnt, "fill": fill}
    if stroke_width:
        kwargs.update({"stroke_width": stroke_width, "stroke_fill": stroke_fill or (255, 255, 255)})
    draw.text((x, y), text, **kwargs)


def draw_lines_centered(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    lines: list[tuple[str, ImageFont.FreeTypeFont, tuple[int, int, int]]],
    *,
    line_gap: int = 8,
) -> None:
    heights = []
    widths = []
    for text, fnt, _ in lines:
        bbox = draw.textbbox((0, 0), text, font=fnt)
        widths.append(bbox[2] - bbox[0])
        heights.append(bbox[3] - bbox[1])
    total = sum(heights) + line_gap * (len(lines) - 1)
    y = (box[1] + box[3] - total) / 2
    for (text, fnt, fill), width, height in zip(lines, widths, heights):
        draw.text(((box[0] + box[2] - width) / 2, y), text, font=fnt, fill=fill)
        y += height + line_gap


def draw_centered_segments(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    segments: list[tuple[str, ImageFont.FreeTypeFont, tuple[int, int, int]]],
) -> None:
    widths = []
    heights = []
    bboxes = []
    for text, fnt, _ in segments:
        bbox = draw.textbbox((0, 0), text, font=fnt)
        bboxes.append(bbox)
        widths.append(bbox[2] - bbox[0])
        heights.append(bbox[3] - bbox[1])
    total_w = sum(widths)
    max_h = max(heights)
    x = (box[0] + box[2] - total_w) / 2
    y_mid = (box[1] + box[3]) / 2
    for (text, fnt, fill), width, _, bbox in zip(segments, widths, heights, bboxes):
        y = y_mid - max_h / 2 - bbox[1]
        draw.text((x, y), text, font=fnt, fill=fill)
        x += width


def pill(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    text: str,
    accent: tuple[int, int, int],
    *,
    fill: tuple[int, int, int] = (255, 255, 255),
    size: int = 30,
    bold: bool = True,
    text_fill: tuple[int, int, int] | None = None,
) -> None:
    draw.rounded_rectangle(box, radius=(box[3] - box[1]) // 2, fill=(*fill, 255), outline=(*accent, 235), width=3)
    centered_text(draw, box, text, full.font(full.TIMES_BOLD if bold else full.TIMES, size), text_fill or accent)


def draw_slide_shell(subtitle: str | None, *, title: str = TITLE) -> Image.Image:
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text(tuple(HP2026_MAIN_HEADER["title_xy"]), title, font=full.font(full.TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]), fill=full.INK)
    if subtitle and HP2026_MAIN_HEADER["subtitle_xy"] and HP2026_MAIN_HEADER["subtitle_font_size"]:
        draw.text(tuple(HP2026_MAIN_HEADER["subtitle_xy"]), subtitle, font=full.font(full.TIMES_ITALIC, HP2026_MAIN_HEADER["subtitle_font_size"]), fill=full.MUTED)
    draw.line((132, HP2026_MAIN_HEADER["divider_y"], W - 132, HP2026_MAIN_HEADER["divider_y"]), fill=(221, 226, 232, 255), width=3)
    full.add_top_right_sphenix_logo_like_slide2(img)
    return img


def draw_number_badge(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], number: str, accent: tuple[int, int, int]) -> None:
    x0, y0, _, y1 = box
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 220))
    badge = (x0 + 38, y0 + 34, x0 + 106, y0 + 102)
    draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(*accent, 235), width=4)
    centered_text(draw, badge, number, full.font(full.TIMES_BOLD, 38), accent)


def panel(
    base: Image.Image,
    box: tuple[int, int, int, int],
    accent: tuple[int, int, int],
    number: str | None = None,
    *,
    show_sidebar: bool = False,
) -> ImageDraw.ImageDraw:
    full.shadow(base, box)
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    if show_sidebar:
        x0, y0, _, y1 = box
        draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 230))
    if number is not None:
        draw_number_badge(draw, box, number, accent)
    return draw


QUESTION_CHIPS = [
    ("core compact?", full.PHOTON_DARK),
    ("shoulders narrow?", full.SPHENIX_BLUE),
    ("not stretched/split?", full.TEAL),
]

CUT_EXAMPLES = [
    ("core compact?", "(E11/E33 > a)", full.PHOTON_DARK, True, 0.62),
    ("shoulders narrow?", "(wη, wφ < b)", full.SPHENIX_BLUE, False, 0.44),
    ("not stretched/split?", "(E32/E35 > c)", full.TEAL, True, 0.72),
]

HANDLE_ROWS = [
    {
        "title": "Compact core",
        "tag": "E11/E33",
        "meaning": "energy concentrated near the seed/core",
        "question": "Is the core compact?",
        "accent": full.PHOTON_DARK,
        "kind": "core",
    },
    {
        "title": "Narrow shoulders",
        "tag": "wη, wφ",
        "meaning": "less lateral spread around the shower center",
        "question": "Are the shoulders narrow?",
        "accent": full.SPHENIX_BLUE,
        "kind": "shoulders",
    },
    {
        "title": "Unsplit local deposit",
        "tag": "E3x2/E3x5",
        "meaning": "not elongated, stretched, or two-lobed",
        "question": "Is the shower unsplit?",
        "accent": full.TEAL,
        "kind": "split",
    },
]


def draw_question_math_pill(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    question: str,
    math_label: str,
    accent: tuple[int, int, int],
    *,
    question_size: int = 31,
    math_size: int = 20,
    fill: tuple[int, int, int] = (255, 255, 255),
    math_bold: bool = True,
) -> None:
    draw.rounded_rectangle(box, radius=(box[3] - box[1]) // 2, fill=(*fill, 255), outline=(*accent, 235), width=4)
    q_font = full.font(full.TIMES_BOLD, question_size)
    m_font = full.font(full.TIMES_BOLD if math_bold else full.TIMES, math_size)
    qb = draw.textbbox((0, 0), question, font=q_font)
    mb = draw.textbbox((0, 0), math_label, font=m_font)
    gap = max(18, int((box[2] - box[0]) * 0.04))
    total_w = (qb[2] - qb[0]) + gap + (mb[2] - mb[0])
    cy = (box[1] + box[3]) / 2
    x = (box[0] + box[2] - total_w) / 2
    qy = cy - (qb[3] - qb[1]) / 2 - qb[1]
    my = cy - (mb[3] - mb[1]) / 2 - mb[1]
    draw.text((x, qy), question, font=q_font, fill=accent)
    draw.text((x + (qb[2] - qb[0]) + gap, my), math_label, font=m_font, fill=full.INK)


def draw_question_chips(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], *, large: bool = True) -> None:
    x0, y0, x1, y1 = box
    gap = 24 if large else 16
    chip_h = 62 if large else 48
    chip_w = (x1 - x0 - 2 * gap) // 3
    for idx, (label, accent) in enumerate(QUESTION_CHIPS):
        cx0 = x0 + idx * (chip_w + gap)
        pill(
            draw,
            (cx0, y0, cx0 + chip_w, y0 + chip_h),
            label,
            accent,
            fill=(248, 251, 253),
            size=30 if large else 22,
        )


def draw_gate_check(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, accent: tuple[int, int, int], detail: str) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=10, fill=(248, 251, 253, 255), outline=(218, 226, 235, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=5, fill=(*accent, 240))
    check = (x0 + 38, y0 + 19, x0 + 78, y0 + 59)
    draw.ellipse(check, outline=(*accent, 235), width=4, fill=(255, 255, 255, 255))
    draw.line((check[0] + 10, check[1] + 21, check[0] + 19, check[1] + 31), fill=(*full.TEAL, 255), width=5)
    draw.line((check[0] + 19, check[1] + 31, check[0] + 32, check[1] + 12), fill=(*full.TEAL, 255), width=5)
    draw.text((x0 + 100, y0 + 15), label, font=full.font(full.TIMES_BOLD, 31), fill=full.INK)
    draw.text((x0 + 100, y0 + 53), detail, font=full.font(full.TIMES_ITALIC, 23), fill=full.MUTED)


def draw_threshold_gate(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    question: str,
    math_label: str,
    threshold_label: str,
    rule: str,
    accent: tuple[int, int, int],
    *,
    pass_high: bool,
    threshold_frac: float = 0.55,
) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=15, fill=(248, 251, 253, 255), outline=(208, 221, 232, 255), width=3)

    qbox = (x0 + 38, y0 + 18, x0 + 840, y0 + 102)
    draw_question_math_pill(draw, qbox, question, math_label, accent, question_size=41, math_size=38)
    draw.text((x0 + 54, y0 + 114), rule, font=full.font(full.TIMES_ITALIC, 32), fill=full.BLUE)

    bar = (x0 + 900, y0 + 58, x1 - 58, y0 + 112)
    draw.rounded_rectangle(bar, radius=25, fill=(229, 235, 241, 255))
    threshold = int(bar[0] + threshold_frac * (bar[2] - bar[0]))
    if pass_high:
        pass_region = (threshold, bar[1], bar[2], bar[3])
        fail_label_xy = (bar[0], bar[1] - 50)
        pass_label_xy = (threshold + 18, bar[1] - 50)
    else:
        pass_region = (bar[0], bar[1], threshold, bar[3])
        pass_label_xy = (bar[0], bar[1] - 50)
        fail_label_xy = (threshold + 18, bar[1] - 50)
    draw.rounded_rectangle(pass_region, radius=25, fill=(*accent, 92))
    draw.line((threshold, bar[1] - 24, threshold, bar[3] + 14), fill=(*accent, 245), width=7)
    cut_font = full.font(full.TIMES_BOLD, 48)
    cut_w, cut_h = full.text_box(draw, threshold_label, cut_font)
    if threshold_label == "b":
        cut_x = threshold + 10
    else:
        cut_x = threshold - cut_w - 10
    draw.text((cut_x, bar[3] - 4), threshold_label, font=cut_font, fill=accent)
    draw.text(pass_label_xy, "pass", font=full.font(full.TIMES_BOLD, 40), fill=accent)
    draw.text(fail_label_xy, "fail", font=full.font(full.TIMES_BOLD, 40), fill=full.MUTED)


def draw_rigid_combination(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=15, fill=(255, 255, 255, 255), outline=(208, 221, 232, 255), width=3)
    centered_text(draw, (x0 + 28, y0 + 28, x1 - 28, y0 + 92), "Rigid combination rule", full.font(full.TIMES_BOLD, 43), full.INK)
    centered_text(draw, (x0 + 42, y0 + 98, x1 - 42, y0 + 142), "all fixed gates must pass", full.font(full.TIMES_ITALIC, 34), full.BLUE)

    row_y = y0 + 182
    for idx, (label, accent) in enumerate(QUESTION_CHIPS):
        cy = row_y + idx * 82
        draw.ellipse((x0 + 62, cy, x0 + 110, cy + 48), fill=(255, 255, 255, 255), outline=(*accent, 235), width=5)
        draw.line((x0 + 74, cy + 27, x0 + 86, cy + 38), fill=(*full.TEAL, 255), width=6)
        draw.line((x0 + 86, cy + 38, x0 + 104, cy + 13), fill=(*full.TEAL, 255), width=6)
        draw.text((x0 + 132, cy + 4), label, font=full.font(full.TIMES_BOLD, 34), fill=full.INK)
        if idx < 2:
            centered_text(draw, (x0 + 282, cy + 50, x0 + 370, cy + 88), "AND", full.font(full.TIMES_BOLD, 25), full.LIGHT_MUTED)

    out = (x0 + 48, y1 - 118, x1 - 48, y1 - 24)
    draw.rounded_rectangle(out, radius=14, fill=(248, 251, 253, 255), outline=(208, 221, 232, 255), width=2)
    draw_lines_centered(
        draw,
        (out[0] + 22, out[1] + 6, out[2] - 22, out[3] - 6),
        [
            ("survives only if", full.font(full.TIMES_BOLD, 34), full.BLUE),
            ("all fixed cuts pass", full.font(full.TIMES_BOLD, 34), full.BLUE),
        ],
        line_gap=2,
    )


def draw_rectangular_cut_plot(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(247, 250, 252, 255), outline=(218, 226, 235, 255), width=2)
    plot = (x0 + 70, y0 + 60, x1 - 54, y1 - 74)
    draw.rectangle(plot, fill=(255, 255, 255, 255), outline=(163, 176, 190, 255), width=3)
    # Light point cloud.
    photon_pts = [(0.24, 0.72), (0.30, 0.76), (0.35, 0.68), (0.41, 0.79), (0.48, 0.71), (0.39, 0.61)]
    bkg_pts = [(0.58, 0.36), (0.69, 0.47), (0.75, 0.31), (0.52, 0.22), (0.82, 0.58), (0.64, 0.20)]
    px0, py0, px1, py1 = plot
    for fx, fy in photon_pts:
        x = px0 + fx * (px1 - px0)
        y = py1 - fy * (py1 - py0)
        draw.ellipse((x - 9, y - 9, x + 9, y + 9), fill=(*full.SPHENIX_BLUE, 210), outline=(255, 255, 255, 230), width=2)
    for fx, fy in bkg_pts:
        x = px0 + fx * (px1 - px0)
        y = py1 - fy * (py1 - py0)
        draw.ellipse((x - 9, y - 9, x + 9, y + 9), fill=(205, 75, 55, 205), outline=(255, 255, 255, 230), width=2)
    # Accepted rectangular region.
    rect = (
        px0 + 0.16 * (px1 - px0),
        py1 - 0.86 * (py1 - py0),
        px0 + 0.53 * (px1 - px0),
        py1 - 0.52 * (py1 - py0),
    )
    draw.rounded_rectangle(rect, radius=8, fill=(245, 181, 34, 42), outline=(*full.PHOTON_DARK, 245), width=5)
    draw.line((rect[2], py0, rect[2], py1), fill=(*full.PHOTON_DARK, 150), width=3)
    draw.line((px0, rect[3], px1, rect[3]), fill=(*full.PHOTON_DARK, 150), width=3)
    draw.text((rect[0] + 20, rect[1] + 16), "pass region", font=full.font(full.TIMES_BOLD, 31), fill=full.PHOTON_DARK)
    draw.text((px0 + 10, py1 + 24), "shoulder width", font=full.font(full.TIMES_BOLD, 28), fill=full.MUTED)
    draw.text((px0 - 54, py0 + 18), "core", font=full.font(full.TIMES_BOLD, 27), fill=full.MUTED)
    draw.text((px0 - 54, py0 + 52), "fraction", font=full.font(full.TIMES_BOLD, 27), fill=full.MUTED)


def draw_funnel_transition(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=20, fill=(255, 252, 243, 255), outline=(238, 220, 172, 255), width=2)
    draw.text((x0 + 56, y0 + 34), "same yes/no handles", font=full.font(full.TIMES_BOLD, 41), fill=full.BLUE)
    draw.text((x0 + 56, y0 + 86), "but now learned together", font=full.font(full.TIMES_ITALIC, 31), fill=full.MUTED)
    left_x = x0 + 58
    chip_y = y0 + 150
    chip_w = 480
    chip_h = 74
    funnel_labels = [
        ("core", full.PHOTON_DARK),
        ("width", full.SPHENIX_BLUE),
        ("split", full.TEAL),
    ]
    for idx, (label, accent) in enumerate(funnel_labels):
        cy = chip_y + idx * 88
        draw.rounded_rectangle(
            (left_x, cy, left_x + chip_w, cy + chip_h),
            radius=17,
            fill=(255, 255, 255, 255),
            outline=(*accent, 235),
            width=3,
        )
        centered_text(
            draw,
            (left_x + 18, cy + 4, left_x + chip_w - 18, cy + chip_h - 4),
            label,
            full.font(full.TIMES_BOLD, 43),
            accent,
        )

    # Draw a clean funnel/combiner: literal enough to read as a funnel, restrained
    # enough that it does not become a cartoon centerpiece.
    funnel_cx = x0 + 878
    top_y = y0 + 98
    mouth_w = 362
    mouth_h = 74
    throat_y = y0 + 296
    stem_bottom = y0 + 382
    accent = full.PHOTON_DARK
    fill = (245, 181, 34, 42)
    # Body fill and side walls.
    body = [
        (funnel_cx - mouth_w // 2 + 16, top_y + mouth_h // 2),
        (funnel_cx + mouth_w // 2 - 16, top_y + mouth_h // 2),
        (funnel_cx + 42, throat_y),
        (funnel_cx + 24, stem_bottom),
        (funnel_cx - 24, stem_bottom),
        (funnel_cx - 42, throat_y),
    ]
    draw.polygon(body, fill=fill, outline=(*accent, 225))
    draw.line(body + [body[0]], fill=(*accent, 230), width=5, joint="curve")
    # Mouth rim.
    outer = (funnel_cx - mouth_w // 2, top_y, funnel_cx + mouth_w // 2, top_y + mouth_h)
    inner = (funnel_cx - mouth_w // 2 + 24, top_y + 10, funnel_cx + mouth_w // 2 - 24, top_y + mouth_h - 10)
    draw.ellipse(outer, fill=(255, 248, 229, 255), outline=(*accent, 240), width=5)
    draw.ellipse(inner, outline=(*accent, 190), width=4)
    # Collar and stem hint.
    draw.arc((funnel_cx - 48, throat_y - 10, funnel_cx + 48, throat_y + 24), start=0, end=180, fill=(*accent, 220), width=4)
    draw.line((funnel_cx - 24, throat_y + 8, funnel_cx - 20, stem_bottom), fill=(*accent, 230), width=5)
    draw.line((funnel_cx + 24, throat_y + 8, funnel_cx + 20, stem_bottom), fill=(*accent, 230), width=5)
    draw.arc((funnel_cx - 23, stem_bottom - 13, funnel_cx + 24, stem_bottom + 15), start=0, end=180, fill=(*accent, 220), width=4)
    centered_text(
        draw,
        (funnel_cx - 132, top_y + 116, funnel_cx + 132, top_y + 170),
        "combine",
        full.font(full.TIMES_BOLD, 37),
        full.PHOTON_DARK,
    )

    for idx, (_, accent) in enumerate(funnel_labels):
        sy = chip_y + idx * 88 + chip_h // 2
        full.draw_arrow(draw, (left_x + chip_w + 30, sy), (funnel_cx - mouth_w // 2 - 18, top_y + mouth_h // 2 + (idx - 1) * 14), fill=accent, width=5)

    bdt = (x1 - 620, y0 + 144, x1 - 58, y1 - 92)
    draw.rounded_rectangle(bdt, radius=16, fill=(239, 246, 250, 255), outline=(*full.SPHENIX_BLUE, 230), width=3)
    draw_lines_centered(
        draw,
        bdt,
        [
            ("next: one learned tree", full.font(full.TIMES_BOLD, 38), full.BLUE),
            ("same variables,", full.font(full.TIMES_BOLD, 31), full.INK),
            ("learned order", full.font(full.TIMES_BOLD, 31), full.INK),
        ],
        line_gap=7,
    )
    full.draw_arrow(draw, (funnel_cx + mouth_w // 2 + 34, y0 + 248), (bdt[0] - 24, y0 + 248), fill=(151, 169, 188), width=6)


def make_silver_funnel_asset(width: int = 420, height: int = 500) -> Image.Image:
    """Build a black-and-white reference-style funnel sketch asset."""
    scale = 4
    w, h = width * scale, height * scale
    img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img, "RGBA")

    def pt(x: float, y: float) -> tuple[int, int]:
        return int(x * scale), int(y * scale)

    def box(x0: float, y0: float, x1: float, y1: float) -> tuple[int, int, int, int]:
        return int(x0 * scale), int(y0 * scale), int(x1 * scale), int(y1 * scale)

    def cubic_points(
        p0: tuple[float, float],
        p1: tuple[float, float],
        p2: tuple[float, float],
        p3: tuple[float, float],
        steps: int = 34,
    ) -> list[tuple[float, float]]:
        pts = []
        for i in range(steps + 1):
            t = i / steps
            mt = 1 - t
            x = mt**3 * p0[0] + 3 * mt**2 * t * p1[0] + 3 * mt * t**2 * p2[0] + t**3 * p3[0]
            y = mt**3 * p0[1] + 3 * mt**2 * t * p1[1] + 3 * mt * t**2 * p2[1] + t**3 * p3[1]
            pts.append((x, y))
        return pts

    ink = (14, 17, 20)
    gray = (128, 138, 148)
    light = (255, 255, 255)

    cx = width / 2
    mouth_y = 62
    mouth_w = 344
    mouth_h = 82
    body_top_y = mouth_y + mouth_h * 0.62
    neck_y = 300
    stem_bottom = 438
    left_lip = (cx - mouth_w / 2 + 10, body_top_y)
    right_lip = (cx + mouth_w / 2 - 10, body_top_y)
    left_neck = (cx - 19, neck_y)
    right_neck = (cx + 19, neck_y)
    left_stem = (cx - 14, stem_bottom)
    right_stem = (cx + 14, stem_bottom)
    left_cone = cubic_points(left_lip, (cx - 148, body_top_y + 82), (cx - 50, neck_y - 18), left_neck, steps=44)
    right_cone = cubic_points(right_lip, (cx + 148, body_top_y + 82), (cx + 50, neck_y - 18), right_neck, steps=44)
    left_stem_curve = cubic_points(left_neck, (cx - 17, neck_y + 52), (cx - 15, stem_bottom - 58), left_stem, steps=24)[1:]
    right_stem_curve = cubic_points(right_neck, (cx + 17, neck_y + 52), (cx + 15, stem_bottom - 58), right_stem, steps=24)[1:]
    left_curve = left_cone + left_stem_curve
    right_curve = right_cone + right_stem_curve
    body = left_curve + list(reversed(right_curve))

    # Soft shadow behind the object.
    shadow = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow, "RGBA")
    sd.polygon([pt(x + 5, y + 7) for x, y in body], fill=(35, 48, 62, 20))
    sd.ellipse(box(cx - mouth_w / 2 + 5, mouth_y + 6, cx + mouth_w / 2 + 5, mouth_y + mouth_h + 6), fill=(35, 48, 62, 18))
    shadow = shadow.filter(ImageFilter.GaussianBlur(5 * scale))
    img.alpha_composite(shadow)
    draw = ImageDraw.Draw(img, "RGBA")

    # Continuous funnel body: reference-like cone tapering directly into a
    # narrow long stem, with no collar/nozzle block.
    draw.polygon([pt(*p) for p in body], fill=(*light, 252))

    body_mask = Image.new("L", (w, h), 0)
    ImageDraw.Draw(body_mask).polygon([pt(*p) for p in body], fill=255)
    hatch_layer = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    hd = ImageDraw.Draw(hatch_layer, "RGBA")

    # Long internal contour lines in the same visual style as the reference.
    hd.line([pt(cx - 122, body_top_y + 20), pt(cx - 70, neck_y - 36)], fill=(*ink, 100), width=2 * scale)
    hd.line([pt(cx + 100, body_top_y + 18), pt(cx + 45, neck_y - 36)], fill=(*ink, 128), width=2 * scale)
    hd.line([pt(cx + 128, body_top_y + 15), pt(cx + 54, neck_y - 20)], fill=(*ink, 178), width=3 * scale)
    hd.line([pt(cx + 132, body_top_y + 28), pt(cx + 63, neck_y - 6)], fill=(*ink, 112), width=2 * scale)

    # Reference-like engraving: dense diagonal hatch on right, lighter short
    # horizontal strokes on the cone and stem. This layer is clipped to the
    # body mask so sketch strokes cannot cross the funnel outline.
    for j in range(34):
        y0h = body_top_y + 14 + j * 7
        x0h = cx + 126 - j * 3.6
        y1h = y0h + 36
        x1h = max(cx + 18, x0h - 54)
        if y1h < neck_y + 18:
            hd.line((pt(x0h, y0h), pt(x1h, y1h)), fill=(*ink, 86), width=1 * scale)
    for j in range(22):
        y0h = body_top_y + 30 + j * 7.5
        x0h = cx + 46 + j * 1.6
        y1h = y0h + 32
        x1h = cx + 118 - j * 2.7
        if y1h < neck_y - 2 and x1h > x0h:
            hd.line((pt(x0h, y0h), pt(x1h, y1h)), fill=(*ink, 42), width=1 * scale)
    for j in range(44):
        y = body_top_y + 20 + j * 7
        span = max(8, 96 - j * 2.35)
        x_start = cx - span
        x_end = cx - span + 34 + (j % 3) * 9
        if y < neck_y + 10:
            hd.line((pt(x_start, y), pt(x_end, y + 1)), fill=(*ink, 58), width=1 * scale)
    for j in range(16):
        y = body_top_y + 38 + j * 8
        x_start = cx + 18 + j * 1.1
        x_end = x_start + 34 - min(18, j)
        if y < neck_y - 10:
            hd.line((pt(x_start, y), pt(x_end, y + 1)), fill=(*ink, 58), width=1 * scale)
    for j in range(15):
        y = neck_y + 22 + j * 10
        hd.line((pt(cx + 3, y), pt(cx + 15, y + 2)), fill=(*ink, 70), width=1 * scale)

    hatch_alpha = ImageChops.multiply(hatch_layer.getchannel("A"), body_mask)
    hatch_layer.putalpha(hatch_alpha)
    img.alpha_composite(hatch_layer)
    draw = ImageDraw.Draw(img, "RGBA")

    # Crisp sidewalls on top of the clipped engraving.
    draw.line([pt(*p) for p in left_curve], fill=(*ink, 255), width=5 * scale, joint="curve")
    draw.line([pt(*p) for p in right_curve], fill=(*ink, 255), width=5 * scale, joint="curve")

    for j in range(10):
        x0h = cx - 116 + j * 16
        draw.line((pt(x0h, mouth_y + 20), pt(x0h + 5, mouth_y + 48 - j * 1.4)), fill=(*ink, 116), width=2 * scale)

    # Rolled elliptical rim with inner opening.
    outer = box(cx - mouth_w / 2, mouth_y, cx + mouth_w / 2, mouth_y + mouth_h)
    mid = box(cx - mouth_w / 2 + 15, mouth_y + 8, cx + mouth_w / 2 - 15, mouth_y + mouth_h - 8)
    inner = box(cx - mouth_w / 2 + 46, mouth_y + 21, cx + mouth_w / 2 - 46, mouth_y + mouth_h - 21)
    draw.ellipse(outer, fill=(*light, 255), outline=(*ink, 255), width=6 * scale)
    draw.ellipse(mid, fill=(*light, 255), outline=(*ink, 245), width=4 * scale)
    draw.ellipse(inner, fill=(*light, 255), outline=(*gray, 170), width=2 * scale)
    rim_hatch = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    rh = ImageDraw.Draw(rim_hatch, "RGBA")
    rim_mask = Image.new("L", (w, h), 0)
    ImageDraw.Draw(rim_mask).ellipse(inner, fill=255)
    # Engraved top-opening strokes, clipped to the inner ellipse.
    for j in range(19):
        x = cx - 110 + j * 10
        y0r = mouth_y + 20 + (j % 4) * 1.3
        y1r = mouth_y + 47 - j * 0.75
        if x < cx + 52:
            rh.line((pt(x, y0r), pt(x + 3, y1r)), fill=(*ink, 116), width=2 * scale)
    for j in range(17):
        y = mouth_y + 22 + j * 2.4
        x0r = cx - 100 + j * 5.2
        x1r = cx + 88 - j * 2.0
        if x1r > x0r + 24:
            rh.line((pt(x0r, y), pt(x1r, y + 0.5)), fill=(*ink, 48), width=1 * scale)
    rim_alpha = ImageChops.multiply(rim_hatch.getchannel("A"), rim_mask)
    rim_hatch.putalpha(rim_alpha)
    img.alpha_composite(rim_hatch)
    draw = ImageDraw.Draw(img, "RGBA")
    # Dark crescent shadows on the left/right inner rim as in the reference.
    draw.arc(box(cx - mouth_w / 2 + 18, mouth_y + 9, cx + mouth_w / 2 - 18, mouth_y + mouth_h - 6), start=166, end=203, fill=(*ink, 255), width=11 * scale)
    draw.arc(box(cx - mouth_w / 2 + 18, mouth_y + 9, cx + mouth_w / 2 - 18, mouth_y + mouth_h - 6), start=337, end=14, fill=(*ink, 255), width=11 * scale)
    draw.arc(box(cx - mouth_w / 2 + 22, mouth_y + 5, cx + mouth_w / 2 - 22, mouth_y + mouth_h - 15), start=188, end=352, fill=(*ink, 74), width=2 * scale)

    # Reference-style rounded tube end, integrated into the same outline.
    draw.arc(box(cx - 18, stem_bottom - 11, cx + 18, stem_bottom + 10), start=0, end=180, fill=(*ink, 255), width=4 * scale)

    return img.resize((width, height), Image.Resampling.LANCZOS)


def draw_tree_input_funnel(
    base: Image.Image,
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    target: tuple[int, int],
) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=14, fill=(248, 251, 253, 255), outline=(208, 221, 232, 255), width=2)

    funnel_cx = (x0 + x1) // 2
    asset_w = 372
    asset_h = 490
    asset_x = funnel_cx - asset_w // 2
    asset_y = y0 + 98
    mouth_y = asset_y + 128
    mouth_w = 310

    falling_inputs = [
        ("E11/E33 > a", full.PHOTON_DARK, (x0 + 122, y0 + 62), (funnel_cx - 54, mouth_y - 10)),
        ("wη, wφ < b", full.SPHENIX_BLUE, (funnel_cx, y0 + 32), (funnel_cx, mouth_y - 4)),
        ("E32/E35 > c", full.TEAL, (x1 - 122, y0 + 62), (funnel_cx + 54, mouth_y - 10)),
    ]
    funnel_asset = make_silver_funnel_asset(asset_w, asset_h)
    base.alpha_composite(funnel_asset, (asset_x, asset_y))

    bdt_box = (funnel_cx - 172, asset_y + int(462 * asset_h / 500), funnel_cx + 172, asset_y + int(520 * asset_h / 500))
    draw.rounded_rectangle(bdt_box, radius=10, fill=(*full.BLUE, 245), outline=(255, 255, 255, 230), width=2)
    centered_text(draw, bdt_box, "Single Decision Tree", full.font(full.TIMES_BOLD, 32), (255, 255, 255))

    for label, accent, label_center, end in falling_inputs:
        fnt = full.font(full.TIMES_BOLD, 33)
        tw, th = full.text_box(draw, label, fnt)
        pos = (int(label_center[0] - tw / 2), int(label_center[1] - th / 2))
        start = (int(label_center[0]), int(pos[1] + th + 12))
        full.draw_arrow(draw, start, end, fill=accent, width=6)
        draw.ellipse((end[0] - 8, end[1] - 8, end[0] + 8, end[1] + 8), fill=(*accent, 225), outline=(255, 255, 255, 230), width=2)
        draw.text(pos, label, font=fnt, fill=accent)

    # The tree begins adjacent to the BDT output; avoid a long diagonal arrow
    # that competes with the decision-tree branches.


def draw_manual_focus(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = panel(base, box, ACCENTS["cuts"], None, show_sidebar=True)
    x0, y0, x1, y1 = box
    title_box = (x0 + 138, y0 + 34, x1 - 70, y0 + 102)
    centered_text(draw, title_box, "Fixed rectangular cuts", full.font(full.TIMES_BOLD, 62), full.INK)
    draw.text((x0 + 126, y0 + 120), "Three fixed gates; all must pass.", font=full.font(full.TIMES_ITALIC, 39), fill=full.BLUE)

    gates = [
        ("core compact?", "(E11/E33 > a)", "a", "require concentrated local core", full.PHOTON_DARK, True, 0.62),
        ("shoulders narrow?", "(wη, wφ < b)", "b", "require small surrounding width", full.SPHENIX_BLUE, False, 0.44),
        ("not stretched/split?", "(E32/E35 > c)", "c", "reject elongated or split deposits", full.TEAL, True, 0.72),
    ]
    left = (x0 + 84, y0 + 182, x0 + 1590, y0 + 836)
    for i, (question, math_label, threshold_label, rule, accent, pass_high, threshold_frac) in enumerate(gates):
        y = left[1] + i * 202
        draw_threshold_gate(
            draw,
            (left[0], y, left[2], y + 168),
            question,
            math_label,
            threshold_label,
            rule,
            accent,
            pass_high=pass_high,
            threshold_frac=threshold_frac,
        )

    right = (x0 + 1648, left[1], x1 - 82, left[1] + 2 * 202 + 168)
    draw_rigid_combination(draw, right)

    takeaway = (x0 + 176, y1 - 136, x1 - 176, y1 - 42)
    draw.rounded_rectangle(takeaway, radius=15, fill=(248, 251, 253, 255), outline=(208, 221, 232, 255), width=2)
    draw_lines_centered(
        draw,
        takeaway,
        [
            ("Transparent, but rigid: interdependencies between variables are not learned.", full.font(full.TIMES_BOLD, 42), full.BLUE),
        ],
        line_gap=7,
    )


def draw_manual_funnel_focus(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = panel(base, box, ACCENTS["cuts"], None)
    x0, y0, x1, y1 = box
    centered_text(draw, (x0 + 128, y0 + 42, x1 - 70, y0 + 106), "From fixed cuts to a learned combination", full.font(full.TIMES_BOLD, 58), full.INK)
    draw.text((x0 + 142, y0 + 118), "Keep the same shower-shape handles, but let the classifier learn how they work together.", font=full.font(full.TIMES_ITALIC, 34), fill=full.BLUE)

    top_note = (x0 + 190, y0 + 184, x1 - 190, y0 + 282)
    draw.rounded_rectangle(top_note, radius=15, fill=(248, 251, 253, 255), outline=(218, 226, 235, 255), width=2)
    draw_lines_centered(
        draw,
        top_note,
        [
            ("The physics questions stay familiar; the combination becomes conditional and learned.", full.font(full.TIMES_BOLD, 38), full.BLUE),
        ],
        line_gap=6,
    )

    draw_funnel_transition(draw, (x0 + 116, y0 + 340, x1 - 116, y1 - 62))


def draw_manual_summary(base: Image.Image, box: tuple[int, int, int, int], *, muted: bool = False) -> None:
    draw = panel(base, box, ACCENTS["cuts"], "1")
    if muted:
        pass
    x0, y0, x1, y1 = box
    centered_text(draw, (x0 + 112, y0 + 36, x1 - 40, y0 + 94), "Manual cuts", full.font(full.TIMES_BOLD, 42), full.INK)
    draw_question_chips(draw, (x0 + 72, y0 + 132, x1 - 62, y0 + 180), large=False)
    gate_rows = [
        ("fixed core gate", full.PHOTON_DARK),
        ("fixed width gate", full.SPHENIX_BLUE),
        ("fixed strip gate", full.TEAL),
    ]
    for i, (gate_label, accent) in enumerate(gate_rows):
        y = y0 + 238 + i * 98
        draw_gate_check(draw, (x0 + 74, y, x1 - 64, y + 70), gate_label, accent, "yes/no threshold")
    note = (x0 + 74, y1 - 138, x1 - 64, y1 - 48)
    draw.rounded_rectangle(note, radius=11, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    draw_lines_centered(draw, note, [("clear, but rigid", full.font(full.TIMES_BOLD, 31), full.BLUE)])


def draw_tree_node(
    draw: ImageDraw.ImageDraw,
    center: tuple[int, int],
    label: str,
    accent: tuple[int, int, int],
    *,
    width: int = 300,
    height: int = 82,
    sublabel: str | None = None,
    label_size: int = 28,
    sublabel_size: int = 19,
) -> tuple[int, int, int, int]:
    cx, cy = center
    box = (cx - width // 2, cy - height // 2, cx + width // 2, cy + height // 2)
    draw.rounded_rectangle(box, radius=13, fill=(255, 255, 255, 255), outline=(*accent, 235), width=4)
    if sublabel:
        draw_lines_centered(
            draw,
            box,
            [
                (label, full.font(full.TIMES_BOLD, label_size), full.INK),
                (sublabel, full.font(full.TIMES_ITALIC, sublabel_size), full.MUTED),
            ],
            line_gap=2,
        )
    else:
        centered_text(draw, box, label, full.font(full.TIMES_BOLD, label_size), full.INK)
    return box


def draw_leaf(
    draw: ImageDraw.ImageDraw,
    center: tuple[int, int],
    label: str,
    accent: tuple[int, int, int],
    *,
    width: int = 230,
    label_size: int = 24,
) -> None:
    cx, cy = center
    box = (cx - width // 2, cy - 40, cx + width // 2, cy + 44)
    fill = (239, 246, 250) if accent == full.SPHENIX_BLUE else (255, 248, 239) if accent == full.PHOTON_DARK else (250, 242, 240)
    draw.rounded_rectangle(box, radius=12, fill=(*fill, 255), outline=(*accent, 235), width=3)
    draw_lines_centered(draw, box, [(label, full.font(full.TIMES_BOLD, label_size), full.INK)])


def draw_tree_diagram(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], *, compact: bool = False) -> None:
    x0, y0, x1, y1 = box
    cx = (x0 + x1) // 2
    if compact:
        root = (cx, y0 + 92)
        left = (cx - 120, y0 + 232)
        right = (cx + 120, y0 + 232)
        ll = (x0 + 78, y0 + 374)
        lr = (x0 + 214, y0 + 374)
        rl = (x1 - 214, y0 + 374)
        rr = (x1 - 78, y0 + 374)
        node_h = 60
        node_w = 196
        leaf_w = 108
        label_size = 22
        sublabel_size = 18
        leaf_size = 18
    else:
        root = (cx, y0 + 122)
        left = (cx - 390, y0 + 360)
        right = (cx + 390, y0 + 360)
        ll = (cx - 630, y0 + 620)
        lr = (cx - 225, y0 + 620)
        rl = (cx + 225, y0 + 620)
        rr = (cx + 630, y0 + 620)
        node_h = 118
        node_w = 520
        leaf_w = 350
        label_size = 44
        sublabel_size = 28
        leaf_size = 36

    links = [
        (root, left, "yes"),
        (root, right, "no"),
        (left, ll, "yes"),
        (left, lr, "no"),
        (right, rl, "yes"),
        (right, rr, "no"),
    ]
    branch_color = (130, 149, 168, 255)
    for start, end, label in links:
        draw.line((start[0], start[1] + node_h // 2, end[0], end[1] - node_h // 2), fill=branch_color, width=6 if not compact else 4)
        mx = (start[0] + end[0]) // 2
        my = (start[1] + end[1]) // 2
        label_x = mx + (-44 if label == "yes" else 16)
        label_y = my + (-42 if not compact else -24)
        full.callout_label(draw, (label_x, label_y), label, full.TEAL if label == "yes" else full.MUTED, bg=(255, 255, 255), size=33 if not compact else 17)

    draw_tree_node(draw, root, "shoulders narrow?", full.SPHENIX_BLUE, width=node_w, height=node_h, sublabel=None, label_size=label_size + (4 if not compact else 0), sublabel_size=sublabel_size)
    draw_tree_node(draw, left, "core compact?", full.PHOTON_DARK, width=node_w, height=node_h, sublabel=None, label_size=label_size + (4 if not compact else 0), sublabel_size=sublabel_size)
    right_label = "split?" if compact else "stretched or split?"
    draw_tree_node(draw, right, right_label, full.TEAL, width=node_w, height=node_h, sublabel=None, label_size=label_size + (4 if not compact else 0), sublabel_size=sublabel_size)
    if compact:
        draw_leaf(draw, ll, "photon", full.SPHENIX_BLUE, width=leaf_w + 8, label_size=leaf_size)
        draw_leaf(draw, lr, "mixed", full.PHOTON_DARK, width=leaf_w, label_size=leaf_size)
        draw_leaf(draw, rl, "bkg-like", ACCENTS["fail"], width=leaf_w + 24, label_size=leaf_size)
        draw_leaf(draw, rr, "mixed", full.PHOTON_DARK, width=leaf_w, label_size=leaf_size)
    else:
        draw_leaf(draw, ll, "photon-like leaf", full.SPHENIX_BLUE, width=leaf_w, label_size=leaf_size)
        draw_leaf(draw, lr, "mixed leaf", full.PHOTON_DARK, width=leaf_w - 36, label_size=leaf_size)
        draw_leaf(draw, rl, "background-like leaf", ACCENTS["fail"], width=leaf_w + 34, label_size=leaf_size)
        draw_leaf(draw, rr, "mixed leaf", full.PHOTON_DARK, width=leaf_w - 36, label_size=leaf_size)


def draw_tree_focus(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = panel(base, box, ACCENTS["tree"], "1")
    x0, y0, x1, y1 = box
    centered_text(draw, (x0 + 128, y0 + 34, x1 - 60, y0 + 100), "One learned decision tree", full.font(full.TIMES_BOLD, 58), full.INK)
    draw.text((x0 + 138, y0 + 114), "Same yes/no inputs; conditional order.", font=full.font(full.TIMES, 37), fill=full.BLUE)
    draw_tree_diagram(draw, (x0 + 88, y0 + 160, x1 - 88, y0 + 850), compact=False)
    takeaway = (x0 + 176, y1 - 118, x1 - 176, y1 - 40)
    draw.rounded_rectangle(takeaway, radius=13, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_lines_centered(
        draw,
        takeaway,
        [
            ("A tree learns the order of the same shower-shape questions.", full.font(full.TIMES_BOLD, 36), full.BLUE),
        ],
        line_gap=8,
    )


def draw_funnel_tree_focus(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = panel(base, box, ACCENTS["tree"], None)
    x0, y0, x1, y1 = box
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*ACCENTS["tree"], 230))
    centered_text(draw, (x0 + 128, y0 + 30, x1 - 60, y0 + 96), "One learned decision tree", full.font(full.TIMES_BOLD, 60), full.INK)
    draw_lines_centered(
        draw,
        (x0 + 176, y0 + 100, x1 - 176, y0 + 190),
        [
            ("Use the same familiar shower-shape cuts and yes/no decisions;", full.font(full.TIMES, 40), full.BLUE),
            ("choose an ML strategy that exploits interdependencies between variables.", full.font(full.TIMES, 40), full.BLUE),
        ],
        line_gap=8,
    )

    left = (x0 + 88, y0 + 208, x0 + 630, y1 - 124)
    right = (x0 + 650, y0 + 190, x1 - 60, y1 - 158)
    tree_cx = (right[0] + right[2]) // 2
    tree_root_target = (tree_cx - 262, right[1] + 164)
    draw_tree_input_funnel(base, draw, left, target=tree_root_target)
    draw_tree_diagram(draw, right, compact=False)

    takeaway = (x0 + 176, y1 - 96, x1 - 176, y1 - 28)
    draw.rounded_rectangle(takeaway, radius=13, fill=(232, 243, 250, 255), outline=(161, 198, 221, 255), width=2)
    draw_centered_segments(
        draw,
        takeaway,
        [
            ("Same variables, learned order: ", full.font(full.TIMES_BOLD, 41), full.BLUE),
            ("fixed cuts become conditional yes/no logic.", full.font(full.TIMES, 41), full.BLUE),
        ],
    )


def draw_tree_summary(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = panel(base, box, ACCENTS["tree"], "1")
    x0, y0, x1, y1 = box
    centered_text(draw, (x0 + 112, y0 + 36, x1 - 40, y0 + 94), "One tree", full.font(full.TIMES_BOLD, 43), full.INK)
    draw_tree_diagram(draw, (x0 + 40, y0 + 120, x1 - 40, y0 + 560), compact=True)
    note = (x0 + 74, y1 - 138, x1 - 64, y1 - 48)
    draw.rounded_rectangle(note, radius=11, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_lines_centered(draw, note, [("conditional order", full.font(full.TIMES_BOLD, 31), full.BLUE)])


def draw_mini_tree_custom(draw: ImageDraw.ImageDraw, origin: tuple[int, int], *, scale: float = 1.0, alpha: int = 210) -> None:
    full.draw_mini_tree(draw, origin, scale=scale, alpha=alpha)


def draw_boost_focus(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = panel(base, box, ACCENTS["boost"], None)
    x0, y0, x1, y1 = box
    centered_text(draw, (x0 + 112, y0 + 36, x1 - 40, y0 + 102), "Boosted ensemble", full.font(full.TIMES_BOLD, 56), full.INK)
    draw.text((x0 + 96, y0 + 116), "many shallow trees combine", font=full.font(full.TIMES_ITALIC, 36), fill=full.BLUE)

    tree_y = y0 + 280
    tree_gap = (x1 - x0 - 390) // 4
    origins = [(x0 + 190 + i * tree_gap, tree_y - (20 if i % 2 else 0)) for i in range(5)]
    for i, origin in enumerate(origins):
        draw_mini_tree_custom(draw, origin, scale=1.0, alpha=228 - 10 * i)
        centered_text(draw, (origin[0] - 36, origin[1] + 154, origin[0] + 36, origin[1] + 194), f"T{i+1}", full.font(full.TIMES_BOLD, 27), full.LIGHT_MUTED)
        if i < len(origins) - 1:
            draw.text((origin[0] + tree_gap // 2 - 10, origin[1] + 144), "+", font=full.font(full.TIMES_BOLD, 38), fill=full.LIGHT_MUTED)

    sum_center = ((x0 + x1) // 2, y0 + 585)
    for origin in origins:
        draw.line((origin[0], origin[1] + 170, sum_center[0], sum_center[1] - 72), fill=(168, 184, 199, 135), width=4)
    draw.ellipse((sum_center[0] - 64, sum_center[1] - 64, sum_center[0] + 64, sum_center[1] + 64), fill=(239, 246, 250, 255), outline=(*full.SPHENIX_BLUE, 230), width=4)
    centered_text(draw, (sum_center[0] - 64, sum_center[1] - 64, sum_center[0] + 64, sum_center[1] + 64), "Σ", full.font(full.TIMES_BOLD, 72), full.BLUE)
    centered_text(draw, (sum_center[0] - 150, sum_center[1] + 72, sum_center[0] + 150, sum_center[1] + 112), "weighted sum", full.font(full.TIMES_BOLD, 28), full.LIGHT_MUTED)

    score_label = (x0 + 86, y0 + 704, x1 - 86, y0 + 748)
    centered_text(draw, score_label, "BDT score", full.font(full.TIMES_BOLD, 35), full.BLUE)
    bar = (x0 + 86, y0 + 764, x1 - 86, y0 + 834)
    draw.rounded_rectangle(bar, radius=24, fill=(232, 236, 240, 255))
    for i in range(bar[2] - bar[0]):
        t = i / max(1, bar[2] - bar[0] - 1)
        r = int(210 * (1 - t) + full.SPHENIX_BLUE[0] * t)
        g = int(79 * (1 - t) + full.SPHENIX_BLUE[1] * t)
        b = int(64 * (1 - t) + full.SPHENIX_BLUE[2] * t)
        draw.line((bar[0] + i, bar[1], bar[0] + i, bar[3]), fill=(r, g, b, 220), width=1)
    draw.text((bar[0], bar[1] - 50), "0", font=full.font(full.TIMES_BOLD, 36), fill=full.MUTED)
    tw1, _ = full.text_box(draw, "1", full.font(full.TIMES_BOLD, 36))
    draw.text((bar[2] - tw1, bar[1] - 50), "1", font=full.font(full.TIMES_BOLD, 36), fill=full.MUTED)
    draw.text((bar[0] + 20, bar[1] + 16), "background-like", font=full.font(full.TIMES_BOLD, 31), fill=(255, 255, 255))
    tw, _ = full.text_box(draw, "photon-like", full.font(full.TIMES_BOLD, 31))
    draw.text((bar[2] - tw - 20, bar[1] + 16), "photon-like", font=full.font(full.TIMES_BOLD, 31), fill=(255, 255, 255))

    note = (x0 + 70, y1 - 126, x1 - 70, y1 - 42)
    draw.rounded_rectangle(note, radius=11, fill=(239, 249, 250, 255), outline=(200, 226, 232, 255), width=2)
    draw_lines_centered(draw, note, [("BDT score ranks candidates, not final purity.", full.font(full.TIMES_BOLD, 31), full.BLUE)])


def draw_boost_full_focus(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = panel(base, box, ACCENTS["boost"], None)
    x0, y0, x1, y1 = box
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*ACCENTS["boost"], 220))
    centered_text(draw, (x0 + 112, y0 + 32, x1 - 80, y0 + 96), "Boosted decision tree ensemble", full.font(full.TIMES_BOLD, 62), full.INK)
    draw_lines_centered(
        draw,
        (x0 + 180, y0 + 102, x1 - 180, y0 + 166),
        [("Repeat the learned yes/no tree; combine many weighted votes into one score.", full.font(full.TIMES, 40), full.BLUE)],
    )

    ensemble = (x0 + 94, y0 + 184, x1 - 94, y0 + 548)
    draw.rounded_rectangle(ensemble, radius=15, fill=(248, 251, 253, 255), outline=(205, 220, 232, 255), width=2)
    draw.text((ensemble[0] + 28, ensemble[1] + 22), "many shallow learned trees", font=full.font(full.TIMES_BOLD, 37), fill=full.INK)

    centers = [ensemble[0] + 260 + i * 405 for i in range(5)]
    card_w, card_h = 270, 236
    tree_accents = [full.TEAL, full.PHOTON_DARK, (116, 132, 148), full.TEAL, full.PHOTON_DARK]
    tree_fills = [
        (248, 251, 253, 255),
        (252, 250, 246, 255),
        (249, 250, 251, 255),
        (248, 251, 253, 255),
        (252, 250, 246, 255),
    ]
    origins: list[tuple[int, int]] = []
    for i, cx_tree in enumerate(centers):
        card = (cx_tree - card_w // 2, ensemble[1] + 96, cx_tree + card_w // 2, ensemble[1] + 96 + card_h)
        fill = tree_fills[i]
        outline = tree_accents[i]
        draw.rounded_rectangle(card, radius=14, fill=fill, outline=(205, 220, 232, 255), width=2)
        draw.text((card[0] + 20, card[1] + 24), f"tree {i + 1}", font=full.font(full.TIMES_BOLD, 34), fill=full.BLUE)
        draw.text((card[2] - 70, card[1] + 26), f"w{i + 1}", font=full.font(full.TIMES_ITALIC, 30), fill=full.MUTED)
        origin = (cx_tree, card[1] + 76)
        origins.append(origin)
        draw_mini_tree_custom(draw, origin, scale=1.28, alpha=242)

    bus_y = ensemble[3] + 20
    bus_x0, bus_x1 = centers[0], centers[-1]
    draw.line((bus_x0, bus_y, bus_x1, bus_y), fill=(118, 139, 158, 205), width=5)
    for cx_tree in centers:
        full.draw_arrow(draw, (cx_tree, ensemble[3] - 4), (cx_tree, bus_y - 5), fill=(118, 139, 158), width=4)
    sum_center = (centers[2], y0 + 646)
    full.draw_arrow(draw, (sum_center[0], bus_y + 8), (sum_center[0], sum_center[1] - 74), fill=(118, 139, 158), width=6)
    draw.ellipse((sum_center[0] - 68, sum_center[1] - 68, sum_center[0] + 68, sum_center[1] + 68), fill=(239, 246, 250, 255), outline=(*ACCENTS["boost"], 240), width=6)
    centered_text(draw, (sum_center[0] - 68, sum_center[1] - 72, sum_center[0] + 68, sum_center[1] + 66), "Σ", full.font(full.TIMES_BOLD, 88), full.BLUE)

    score_label = (x0 + 170, y0 + 714, x1 - 170, y0 + 780)
    main_fnt = full.font(full.TIMES_BOLD, 39)
    sub_fnt = full.font(full.TIMES_BOLD, 25)
    score_parts = [
        ("BDT score = weighted sum of tree votes = w", main_fnt, 0, full.BLUE),
        ("1", sub_fnt, 13, full.BLUE),
        ("T", main_fnt, 0, full.BLUE),
        ("1", sub_fnt, 13, full.BLUE),
        (" + w", main_fnt, 0, full.BLUE),
        ("2", sub_fnt, 13, full.BLUE),
        ("T", main_fnt, 0, full.BLUE),
        ("2", sub_fnt, 13, full.BLUE),
        (" + w", main_fnt, 0, full.BLUE),
        ("3", sub_fnt, 13, full.BLUE),
        ("T", main_fnt, 0, full.BLUE),
        ("3", sub_fnt, 13, full.BLUE),
        (" + ...", main_fnt, 0, full.BLUE),
    ]
    total_w = sum(full.text_box(draw, text, fnt)[0] for text, fnt, _, _ in score_parts)
    max_h = max(full.text_box(draw, text, fnt)[1] + max(0, dy) for text, fnt, dy, _ in score_parts)
    sx = (score_label[0] + score_label[2] - total_w) / 2
    sy = (score_label[1] + score_label[3] - max_h) / 2
    for text, fnt, dy, fill in score_parts:
        draw.text((sx, sy + dy), text, font=fnt, fill=fill)
        sx += full.text_box(draw, text, fnt)[0]
    bar = (x0 + 170, y0 + 802, x1 - 170, y0 + 890)
    draw.rounded_rectangle(bar, radius=34, fill=(232, 236, 240, 255))
    for i in range(bar[2] - bar[0]):
        t = i / max(1, bar[2] - bar[0] - 1)
        r = int(210 * (1 - t) + full.SPHENIX_BLUE[0] * t)
        g = int(79 * (1 - t) + full.SPHENIX_BLUE[1] * t)
        b = int(64 * (1 - t) + full.SPHENIX_BLUE[2] * t)
        draw.line((bar[0] + i, bar[1], bar[0] + i, bar[3]), fill=(r, g, b, 224), width=1)
    draw.text((bar[0], bar[1] - 58), "0", font=full.font(full.TIMES_BOLD, 43), fill=full.MUTED)
    one_w, _ = full.text_box(draw, "1", full.font(full.TIMES_BOLD, 43))
    draw.text((bar[2] - one_w, bar[1] - 58), "1", font=full.font(full.TIMES_BOLD, 43), fill=full.MUTED)
    draw.text((bar[0] + 34, bar[1] + 19), "background-like", font=full.font(full.TIMES_BOLD, 39), fill=(255, 255, 255))
    p_w, _ = full.text_box(draw, "photon-like", full.font(full.TIMES_BOLD, 39))
    draw.text((bar[2] - p_w - 34, bar[1] + 19), "photon-like", font=full.font(full.TIMES_BOLD, 39), fill=(255, 255, 255))

    note = (x0 + 210, y1 - 76, x1 - 210, y1 - 18)
    draw.rounded_rectangle(note, radius=14, fill=(238, 247, 247, 255), outline=(176, 211, 214, 255), width=2)
    draw_lines_centered(
        draw,
        note,
        [("Boosting reduces dependence on one initial split by averaging many learned tree votes.", full.font(full.TIMES_BOLD, 35), full.BLUE)],
    )


def slide11_manual() -> Image.Image:
    img = draw_slide_shell(None, title="From shower-shape cuts to a BDT score")
    draw = ImageDraw.Draw(img, "RGBA")
    top_card = (132, 258, W - 132, 996)
    boost = (132, 1030, W - 132, 1282)
    full.shadow(img, top_card)
    draw.rounded_rectangle(top_card, radius=13, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((top_card[0], top_card[1], top_card[0] + 12, top_card[3]), radius=6, fill=(*ACCENTS["tree"], 230))

    centered_text(
        draw,
        (top_card[0] + 88, top_card[1] + 22, top_card[2] - 88, top_card[1] + 66),
        "One decision tree",
        full.font(full.TIMES_BOLD, 56),
        full.INK,
    )
    emphasis = "Conditional splits exploit interdependencies that fixed cuts miss."
    emphasis_font = full.font(full.TIMES_ITALIC, 46)
    emphasis_w, emphasis_h = full.text_box(draw, emphasis, emphasis_font)
    emphasis_x = (top_card[0] + top_card[2] - emphasis_w) / 2
    emphasis_y = top_card[1] + 88
    draw.text((emphasis_x, emphasis_y), emphasis, font=emphasis_font, fill=full.BLUE)

    def numbered_node(
        center: tuple[int, int],
        number: str,
        text: str,
        accent: tuple[int, int, int],
        *,
        width: int = 520,
        height: int = 106,
        text_size: int = 42,
    ) -> tuple[int, int, int, int]:
        cx, cy = center
        box = (cx - width // 2, cy - height // 2, cx + width // 2, cy + height // 2)
        draw.rounded_rectangle(box, radius=16, fill=(255, 255, 255, 255), outline=(*accent, 238), width=4)
        chip = (box[0] + 22, cy - 24, box[0] + 70, cy + 24)
        draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 238), width=3)
        centered_text(draw, chip, number, full.font(full.TIMES_BOLD, 25), accent)
        centered_text(draw, (box[0] + 78, box[1] + 8, box[2] - 18, box[3] - 8), text, full.font(full.TIMES_BOLD, text_size), full.INK)
        return box

    def outcome_leaf(
        center: tuple[int, int],
        text: str,
        accent: tuple[int, int, int],
        *,
        width: int = 340,
    ) -> tuple[int, int, int, int]:
        cx, cy = center
        box = (cx - width // 2, cy - 42, cx + width // 2, cy + 46)
        if accent == full.SPHENIX_BLUE:
            fill = (239, 247, 252, 255)
        elif accent == ACCENTS["fail"]:
            fill = (253, 244, 242, 255)
        else:
            fill = (255, 250, 241, 255)
        draw.rounded_rectangle(box, radius=13, fill=fill, outline=(*accent, 235), width=3)
        centered_text(draw, box, text, full.font(full.TIMES_BOLD, 32), full.INK)
        return box

    tree = (top_card[0] + 250, top_card[1] + 126, top_card[2] - 250, top_card[1] + 682)
    cx = (tree[0] + tree[2]) // 2
    root = (cx, tree[1] + 105)
    left = (cx - 410, tree[1] + 315)
    right = (cx + 410, tree[1] + 315)
    leaves = [
        (cx - 660, tree[1] + 545, "more photon-like", ACCENTS["fail"], 378),
        (cx - 245, tree[1] + 545, "mixed", full.PHOTON_DARK, 270),
        (cx + 245, tree[1] + 545, "more background-like", full.SPHENIX_BLUE, 430),
        (cx + 660, tree[1] + 545, "mixed", full.PHOTON_DARK, 270),
    ]

    root_box = numbered_node(root, "1", "Core compact?", full.PHOTON_DARK, text_size=44)
    left_box = numbered_node(left, "2", "Shoulders narrow?", full.SPHENIX_BLUE, text_size=44)
    right_box = numbered_node(right, "3", "Stretched or split?", full.TEAL, text_size=42)

    branch_color = (132, 151, 170, 255)
    links = [
        (root, left, "yes"),
        (root, right, "no"),
        (left, leaves[0][0:2], "yes"),
        (left, leaves[1][0:2], "no"),
        (right, leaves[2][0:2], "yes"),
        (right, leaves[3][0:2], "no"),
    ]
    for start, end, label in links:
        start_y = start[1] + 54
        end_y = end[1] - 54 if end[1] <= tree[1] + 330 else end[1] - 48
        draw.line((start[0], start_y, end[0], end_y), fill=branch_color, width=6)
        mx = (start[0] + end[0]) // 2
        my = (start_y + end_y) // 2
        full.callout_label(
            draw,
            (mx + (-44 if label == "yes" else 14), my - 42),
            label,
            full.TEAL if label == "yes" else full.MUTED,
            bg=(255, 255, 255),
            size=31,
        )
    for leaf in leaves:
        outcome_leaf((leaf[0], leaf[1]), leaf[2], leaf[3], width=leaf[4])

    full.shadow(img, boost)
    draw.rounded_rectangle(boost, radius=15, fill=(248, 251, 253, 255), outline=(210, 222, 233, 255), width=2)
    draw.text((boost[0] + 40, boost[1] + 18), "Boosted Decision Tree (BDT):", font=full.font(full.TIMES_BOLD, 48), fill=full.INK)

    card_y = boost[1] + 102
    card_w, card_h = 222, 76
    card_step = 290
    card_x0 = boost[0] + 40
    tree_cards: list[tuple[int, int, int, int]] = []
    for i in range(3):
        x = card_x0 + i * card_step
        card = (x, card_y, x + card_w, card_y + card_h)
        tree_cards.append(card)
        fill = (239, 247, 252, 255) if i % 2 == 0 else (255, 250, 241, 255)
        draw.rounded_rectangle(card, radius=12, fill=fill, outline=(198, 213, 225, 255), width=2)
        draw_lines_centered(
            draw,
            card,
            [
                (f"tree {i + 1}", full.font(full.TIMES_BOLD, 35), full.BLUE),
            ],
        )

    sigma = (tree_cards[-1][2] + 152, card_y + card_h // 2)
    plus_font = full.font(full.TIMES_BOLD, 37)
    for left_card, right_card in zip(tree_cards, tree_cards[1:]):
        plus_w, plus_h = full.text_box(draw, "+", plus_font)
        px = (left_card[2] + right_card[0] - plus_w) / 2
        py = card_y + (card_h - plus_h) / 2 - 2
        draw.text((px, py), "+", font=plus_font, fill=full.LIGHT_MUTED)
    ellipsis_text = "+ ..."
    ellipsis_font = full.font(full.TIMES_BOLD, 41)
    ellipsis_w, ellipsis_h = full.text_box(draw, ellipsis_text, ellipsis_font)
    gap_left = tree_cards[-1][2]
    gap_right = sigma[0] - 48
    draw.text(((gap_left + gap_right - ellipsis_w) / 2, card_y + (card_h - ellipsis_h) / 2 - 1), ellipsis_text, font=ellipsis_font, fill=(118, 137, 156))
    full.draw_arrow(draw, (gap_right + 12, sigma[1]), (sigma[0] - 46, sigma[1]), fill=(126, 145, 164), width=5)
    draw.ellipse((sigma[0] - 42, sigma[1] - 42, sigma[0] + 42, sigma[1] + 42), fill=(239, 246, 250, 255), outline=(*ACCENTS["boost"], 230), width=4)
    centered_text(draw, (sigma[0] - 42, sigma[1] - 47, sigma[0] + 42, sigma[1] + 37), "Σ", full.font(full.TIMES_BOLD, 60), full.BLUE)

    divider = boost[0] + 1660
    score_title = "photon-ID score"
    score_font = full.font(full.TIMES_BOLD, 33)
    score_w, _ = full.text_box(draw, score_title, score_font)
    score_center_x = (sigma[0] + 96 + divider - 58) // 2
    draw.text((score_center_x - score_w // 2, boost[1] + 70), score_title, font=score_font, fill=full.INK)

    axis_left = sigma[0] + 118
    axis_right = divider - 72
    axis_y = sigma[1]
    full.draw_arrow(draw, (sigma[0] + 44, axis_y), (axis_left - 28, axis_y), fill=(126, 145, 164), width=5)

    # A bare score axis reads cleaner than a nested capsule here.
    draw.line((axis_left, axis_y, axis_right, axis_y), fill=(204, 214, 224, 255), width=13)
    axis_w = max(1, axis_right - axis_left)
    for i in range(axis_w + 1):
        t = i / axis_w
        r = int(full.SPHENIX_BLUE[0] * (1 - t) + ACCENTS["fail"][0] * t)
        g = int(full.SPHENIX_BLUE[1] * (1 - t) + ACCENTS["fail"][1] * t)
        b = int(full.SPHENIX_BLUE[2] * (1 - t) + ACCENTS["fail"][2] * t)
        draw.line((axis_left + i, axis_y - 3, axis_left + i, axis_y + 3), fill=(r, g, b, 230), width=1)
    for x, color in [(axis_left, full.SPHENIX_BLUE), (axis_right, ACCENTS["fail"])]:
        draw.ellipse((x - 17, axis_y - 17, x + 17, axis_y + 17), fill=(255, 255, 255, 255), outline=(*color, 245), width=4)
        draw.ellipse((x - 10, axis_y - 10, x + 10, axis_y + 10), fill=(*color, 255))

    endpoint_font = full.font(full.TIMES_BOLD, 45)
    zero_w, zero_h = full.text_box(draw, "0", endpoint_font)
    one_w, _ = full.text_box(draw, "1", endpoint_font)
    draw.text((axis_left - zero_w / 2, axis_y - 62), "0", font=endpoint_font, fill=full.SPHENIX_BLUE)
    draw.text((axis_right - one_w / 2, axis_y - 62), "1", font=endpoint_font, fill=ACCENTS["fail"])
    label_font = full.font(full.TIMES_BOLD, 26)
    draw.text((axis_left + 16, axis_y + 32), "background-like", font=label_font, fill=full.SPHENIX_BLUE)
    p_w, _ = full.text_box(draw, "photon-like", label_font)
    draw.text((axis_right - p_w - 16, axis_y + 32), "photon-like", font=label_font, fill=ACCENTS["fail"])

    draw.line((divider, boost[1] + 26, divider, boost[3] - 26), fill=(216, 226, 236, 255), width=3)
    bullet_x = divider + 42
    bullet_font = full.font(full.TIMES_BOLD, 40)
    bullet_body = full.font(full.TIMES, 39)
    bullet_rows = [
        ("Uses", "correlated information"),
        ("Returns", "one photon-ID score"),
    ]
    for i, (lead, body_text) in enumerate(bullet_rows):
        y = boost[1] + 52 + i * 78
        draw.ellipse((bullet_x, y + 16, bullet_x + 17, y + 33), fill=full.INK)
        draw.text((bullet_x + 34, y), lead, font=bullet_font, fill=full.INK)
        lead_w, _ = full.text_box(draw, lead, bullet_font)
        draw.text((bullet_x + 48 + lead_w, y), body_text, font=bullet_body, fill=full.BLUE)

    full.draw_hp2026_identity_footer(img)
    return img


def slide12_tree() -> Image.Image:
    img = draw_slide_shell(None, title="From fixed cuts to a single decision tree")
    draw_funnel_tree_focus(img, (132, 304, W - 132, 1282))
    full.draw_hp2026_identity_footer(img)
    return img


def slide13_boost() -> Image.Image:
    img = draw_slide_shell(None, title="From one tree to a photon-ID score")
    draw_boost_full_focus(img, (132, 304, W - 132, 1282))
    full.draw_hp2026_identity_footer(img)
    return img


def draw_simple_gate(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    accent: tuple[int, int, int],
) -> None:
    draw.rounded_rectangle(box, radius=18, fill=(248, 251, 253, 255), outline=(*accent, 220), width=4)
    check = (box[0] + 28, box[1] + 22, box[0] + 82, box[1] + 76)
    draw.ellipse(check, fill=(255, 255, 255, 255), outline=(*accent, 235), width=5)
    draw.line((check[0] + 13, check[1] + 30, check[0] + 25, check[1] + 42), fill=(*full.TEAL, 255), width=7)
    draw.line((check[0] + 25, check[1] + 42, check[0] + 43, check[1] + 17), fill=(*full.TEAL, 255), width=7)
    centered_text(draw, (box[0] + 108, box[1] + 8, box[2] - 22, box[3] - 8), label, full.font(full.TIMES_BOLD, 36), full.INK)


def draw_bridge_tree(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
) -> None:
    x0, y0, x1, y1 = box
    cx = (x0 + x1) // 2
    root = (cx, y0 + 74)
    left = (cx - 360, y0 + 230)
    right = (cx + 360, y0 + 230)
    leaves = [
        (cx - 545, y0 + 386, "photon-like", full.SPHENIX_BLUE, 286),
        (cx - 205, y0 + 386, "mixed", full.PHOTON_DARK, 210),
        (cx + 205, y0 + 386, "mixed", full.PHOTON_DARK, 210),
        (cx + 545, y0 + 386, "background-like", ACCENTS["fail"], 320),
    ]
    node_w, node_h = 420, 90
    leaf_h = 70
    branch = (130, 149, 168, 255)
    links = [
        (root, left, "yes"),
        (root, right, "no"),
        (left, leaves[0][0:2], "yes"),
        (left, leaves[1][0:2], "no"),
        (right, leaves[2][0:2], "no"),
        (right, leaves[3][0:2], "yes"),
    ]
    for start, end, yn in links:
        draw.line((start[0], start[1] + node_h // 2, end[0], end[1] - (node_h // 2 if end in [left, right] else leaf_h // 2)), fill=branch, width=5)
        mx = (start[0] + end[0]) // 2
        my = (start[1] + end[1]) // 2
        full.callout_label(draw, (mx - 38, my - 40), yn, full.TEAL if yn == "yes" else full.BLUE, bg=(255, 255, 255), size=32)

    draw_tree_node(draw, root, "shoulders narrow?", full.SPHENIX_BLUE, width=node_w, height=node_h, label_size=37)
    draw_tree_node(draw, left, "core compact?", full.PHOTON_DARK, width=node_w, height=node_h, label_size=37)
    draw_tree_node(draw, right, "stretched or split?", full.TEAL, width=node_w, height=node_h, label_size=37)
    for lx, ly, label, accent, leaf_w in leaves:
        fill = (239, 246, 250) if accent == full.SPHENIX_BLUE else (255, 248, 239) if accent == full.PHOTON_DARK else (250, 242, 240)
        leaf = (lx - leaf_w // 2, ly - leaf_h // 2, lx + leaf_w // 2, ly + leaf_h // 2)
        draw.rounded_rectangle(leaf, radius=13, fill=(*fill, 255), outline=(*accent, 235), width=3)
        centered_text(draw, leaf, label, full.font(full.TIMES_BOLD, 31), full.INK)


def draw_score_bar(draw: ImageDraw.ImageDraw, bar: tuple[int, int, int, int]) -> None:
    draw.rounded_rectangle(bar, radius=(bar[3] - bar[1]) // 2, fill=(232, 236, 240, 255))
    for i in range(bar[2] - bar[0]):
        t = i / max(1, bar[2] - bar[0] - 1)
        r = int(205 * (1 - t) + full.SPHENIX_BLUE[0] * t)
        g = int(74 * (1 - t) + full.SPHENIX_BLUE[1] * t)
        b = int(62 * (1 - t) + full.SPHENIX_BLUE[2] * t)
        draw.line((bar[0] + i, bar[1], bar[0] + i, bar[3]), fill=(r, g, b, 226), width=1)
    draw.text((bar[0], bar[1] - 54), "0", font=full.font(full.TIMES_BOLD, 42), fill=full.MUTED)
    one_w, _ = full.text_box(draw, "1", full.font(full.TIMES_BOLD, 42))
    draw.text((bar[2] - one_w, bar[1] - 54), "1", font=full.font(full.TIMES_BOLD, 42), fill=full.MUTED)
    draw.text((bar[0] + 28, bar[1] + 18), "background-like", font=full.font(full.TIMES_BOLD, 36), fill=(255, 255, 255))
    p_w, _ = full.text_box(draw, "photon-like", full.font(full.TIMES_BOLD, 36))
    draw.text((bar[2] - p_w - 28, bar[1] + 18), "photon-like", font=full.font(full.TIMES_BOLD, 36), fill=(255, 255, 255))


def wrap_lines_by_width(draw: ImageDraw.ImageDraw, text: str, max_width: int, fnt) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        candidate = " ".join(current + [word])
        if current and full.text_box(draw, candidate, fnt)[0] > max_width:
            lines.append(" ".join(current))
            current = [word]
        else:
            current.append(word)
    if current:
        lines.append(" ".join(current))
    return lines


def draw_wrapped_left(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    fnt,
    fill: tuple[int, int, int],
    max_width: int,
    *,
    line_gap: int = 6,
) -> int:
    x, y = xy
    for line in wrap_lines_by_width(draw, text, max_width, fnt):
        draw.text((x, y), line, font=fnt, fill=fill)
        y += full.text_box(draw, line, fnt)[1] + line_gap
    return y


def draw_tag_pill(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    accent: tuple[int, int, int],
    *,
    size: int = 29,
) -> tuple[int, int, int, int]:
    fnt = full.font(full.TIMES_BOLD, size)
    tw, th = full.text_box(draw, text, fnt)
    box = (xy[0], xy[1], xy[0] + tw + 36, xy[1] + th + 20)
    draw.rounded_rectangle(box, radius=(box[3] - box[1]) // 2, fill=(255, 255, 255, 255), outline=(*accent, 220), width=3)
    centered_text(draw, box, text, fnt, accent)
    return box


def shower_values(kind: str, photon_like: bool) -> list[list[float]]:
    if kind == "core" and photon_like:
        return [
            [0.03, 0.05, 0.07, 0.05, 0.03],
            [0.05, 0.16, 0.34, 0.17, 0.05],
            [0.07, 0.38, 0.98, 0.42, 0.07],
            [0.05, 0.18, 0.36, 0.18, 0.05],
            [0.03, 0.05, 0.07, 0.05, 0.03],
        ]
    if kind == "core":
        return [
            [0.12, 0.20, 0.28, 0.20, 0.12],
            [0.22, 0.40, 0.55, 0.40, 0.22],
            [0.30, 0.58, 0.80, 0.60, 0.32],
            [0.22, 0.42, 0.56, 0.42, 0.22],
            [0.12, 0.20, 0.30, 0.20, 0.12],
        ]
    if kind == "shoulders" and photon_like:
        return [
            [0.01, 0.02, 0.04, 0.02, 0.01],
            [0.02, 0.08, 0.20, 0.08, 0.02],
            [0.04, 0.18, 0.95, 0.18, 0.04],
            [0.02, 0.08, 0.20, 0.08, 0.02],
            [0.01, 0.02, 0.04, 0.02, 0.01],
        ]
    if kind == "shoulders":
        return [
            [0.10, 0.18, 0.25, 0.18, 0.10],
            [0.18, 0.34, 0.46, 0.34, 0.18],
            [0.25, 0.50, 0.82, 0.50, 0.25],
            [0.18, 0.34, 0.46, 0.34, 0.18],
            [0.10, 0.18, 0.25, 0.18, 0.10],
        ]
    if photon_like:
        return [
            [0.02, 0.04, 0.06, 0.04, 0.02],
            [0.03, 0.10, 0.20, 0.10, 0.03],
            [0.05, 0.32, 0.96, 0.34, 0.05],
            [0.03, 0.10, 0.20, 0.10, 0.03],
            [0.02, 0.04, 0.06, 0.04, 0.02],
        ]
    return [
        [0.02, 0.06, 0.10, 0.08, 0.04],
        [0.06, 0.18, 0.30, 0.25, 0.11],
        [0.18, 0.62, 0.98, 0.82, 0.46],
        [0.16, 0.54, 0.74, 0.60, 0.36],
        [0.04, 0.12, 0.18, 0.14, 0.06],
    ]


def draw_handle_grid(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    *,
    cell: int,
    kind: str,
    photon_like: bool,
    accent: tuple[int, int, int],
) -> tuple[int, int, int, int]:
    grid = full.draw_small_tower_grid(draw, xy, cell, shower_values(kind, photon_like))
    gx0, gy0, gx1, gy1 = grid
    if kind == "core":
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*full.SPHENIX_BLUE, 230), width=4)
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 3 * cell, gy0 + 3 * cell), outline=(*full.PHOTON_DARK, 245), width=4)
    elif kind == "shoulders":
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 215), outline=(155, 164, 174, 215), width=2)
        draw.line((center_cell[0] + 8, center_cell[1] + 8, center_cell[2] - 8, center_cell[3] - 8), fill=(155, 164, 174, 190), width=3)
        draw.line((center_cell[0] + 8, center_cell[3] - 8, center_cell[2] - 8, center_cell[1] + 8), fill=(155, 164, 174, 190), width=3)
        span = (gx1 - gx0) * (0.50 if photon_like else 0.86)
        y = gy0 + 14
        x_start = (gx0 + gx1 - span) / 2
        x_end = (gx0 + gx1 + span) / 2
        draw.line((x_start, y, x_end, y), fill=(*accent, 245), width=5)
        draw.polygon([(x_start, y), (x_start + 14, y - 10), (x_start + 14, y + 10)], fill=(*accent, 245))
        draw.polygon([(x_end, y), (x_end - 14, y - 10), (x_end - 14, y + 10)], fill=(*accent, 245))
    else:
        draw.rectangle((gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell), outline=(*full.SPHENIX_BLUE, 230), width=4)
        draw.rectangle((gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*full.PHOTON_DARK, 245), width=4)
    draw.ellipse((gx0 + 2.5 * cell - 8, gy0 + 2.5 * cell - 8, gx0 + 2.5 * cell + 8, gy0 + 2.5 * cell + 8), fill=(0, 0, 0, 255))
    return grid


def draw_embedded_label(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    accent: tuple[int, int, int],
    *,
    size: int = 24,
    pad_x: int = 12,
    pad_y: int = 5,
) -> tuple[int, int, int, int]:
    fnt = full.font(full.TIMES_BOLD, size)
    tw, th = full.text_box(draw, text, fnt)
    box = (xy[0], xy[1], xy[0] + tw + 2 * pad_x, xy[1] + th + 2 * pad_y)
    draw.rounded_rectangle(box, radius=(box[3] - box[1]) // 2, fill=(255, 255, 255, 235), outline=(*accent, 230), width=2)
    centered_text(draw, box, text, fnt, accent)
    return box


def draw_double_arrow(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    fill: tuple[int, int, int],
    width: int = 5,
    head: int = 13,
) -> None:
    draw.line((start[0], start[1], end[0], end[1]), fill=(*fill, 245), width=width)
    if abs(end[0] - start[0]) >= abs(end[1] - start[1]):
        y = start[1]
        draw.polygon([(start[0], y), (start[0] + head, y - head * 0.7), (start[0] + head, y + head * 0.7)], fill=(*fill, 245))
        draw.polygon([(end[0], y), (end[0] - head, y - head * 0.7), (end[0] - head, y + head * 0.7)], fill=(*fill, 245))
    else:
        x = start[0]
        draw.polygon([(x, start[1]), (x - head * 0.7, start[1] + head), (x + head * 0.7, start[1] + head)], fill=(*fill, 245))
        draw.polygon([(x, end[1]), (x - head * 0.7, end[1] - head), (x + head * 0.7, end[1] - head)], fill=(*fill, 245))


def draw_handle_grid_v4(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    *,
    cell: int,
    kind: str,
    photon_like: bool,
    accent: tuple[int, int, int],
) -> tuple[int, int, int, int]:
    grid = full.draw_small_tower_grid(draw, xy, cell, shower_values(kind, photon_like))
    gx0, gy0, gx1, gy1 = grid
    grid_w = gx1 - gx0
    if kind == "core":
        rect_3x3 = (gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell)
        rect_1x1 = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(rect_3x3, outline=(*full.SPHENIX_BLUE, 245), width=5)
        draw.rectangle(rect_1x1, outline=(*full.PHOTON_DARK, 255), width=6)
        draw_embedded_label(draw, (gx0 + cell + 7, gy0 + cell - 30), "3x3", full.SPHENIX_BLUE, size=22, pad_x=9, pad_y=3)
        draw_embedded_label(draw, (gx0 + 3 * cell - 2, gy0 + 2 * cell + 7), "1x1", full.PHOTON_DARK, size=22, pad_x=8, pad_y=3)
        draw_embedded_label(draw, (gx0 + 10, gy1 - 38), "E11/E33", full.PHOTON_DARK, size=25, pad_x=12, pad_y=4)
    elif kind == "shoulders":
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 218), outline=(155, 164, 174, 220), width=2)
        draw.line((center_cell[0] + 9, center_cell[1] + 9, center_cell[2] - 9, center_cell[3] - 9), fill=(155, 164, 174, 200), width=3)
        draw.line((center_cell[0] + 9, center_cell[3] - 9, center_cell[2] - 9, center_cell[1] + 9), fill=(155, 164, 174, 200), width=3)
        h_span = grid_w * (0.50 if photon_like else 0.86)
        hx0 = (gx0 + gx1 - h_span) / 2
        hx1 = (gx0 + gx1 + h_span) / 2
        hy = gy0 + 24
        draw_double_arrow(draw, (hx0, hy), (hx1, hy), fill=full.SPHENIX_BLUE, width=5, head=13)
        phi_label = "wφ"
        phi_font = full.font(full.TIMES_BOLD, 24)
        phi_w, phi_h = full.text_box(draw, phi_label, phi_font)
        draw.text(((hx0 + hx1 - phi_w) / 2, hy + 8), phi_label, font=phi_font, fill=full.SPHENIX_BLUE)
        v_span = grid_w * (0.48 if photon_like else 0.82)
        vx = gx1 - 22
        vy0 = (gy0 + gy1 - v_span) / 2
        vy1 = (gy0 + gy1 + v_span) / 2
        draw_double_arrow(draw, (vx, vy0), (vx, vy1), fill=full.TEAL, width=5, head=13)
        eta_label = "wη"
        eta_font = full.font(full.TIMES_BOLD, 24)
        draw.text((vx - 52, (vy0 + vy1) / 2 - 18), eta_label, font=eta_font, fill=full.TEAL)
        draw_embedded_label(draw, (gx0 + 10, gy1 - 38), "wη, wφ", full.SPHENIX_BLUE, size=25, pad_x=12, pad_y=4)
    else:
        rect_3x5 = (gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell)
        rect_3x2 = (gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 4 * cell)
        draw.rectangle(rect_3x5, outline=(*full.SPHENIX_BLUE, 245), width=5)
        draw.rectangle(rect_3x2, outline=(*full.PHOTON_DARK, 255), width=6)
        draw_embedded_label(draw, (gx0 + cell + 8, gy0 + 8), "3x5", full.SPHENIX_BLUE, size=22, pad_x=9, pad_y=3)
        draw_embedded_label(draw, (gx0 + cell + 8, gy0 + 2 * cell + 8), "3x2", full.PHOTON_DARK, size=22, pad_x=9, pad_y=3)
        draw_embedded_label(draw, (gx0 + 10, gy1 - 38), "E3x2/E3x5", full.TEAL, size=24, pad_x=10, pad_y=4)
    draw.ellipse((gx0 + 2.5 * cell - 9, gy0 + 2.5 * cell - 9, gx0 + 2.5 * cell + 9, gy0 + 2.5 * cell + 9), fill=(0, 0, 0, 255))
    return grid


def draw_handle_grid_v5(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    *,
    cell: int,
    kind: str,
    photon_like: bool,
    accent: tuple[int, int, int],
) -> tuple[int, int, int, int]:
    """Large tower-grid cue with geometry only; variable definitions live in the question boxes."""
    grid = full.draw_small_tower_grid(draw, xy, cell, shower_values(kind, photon_like))
    gx0, gy0, gx1, gy1 = grid
    grid_w = gx1 - gx0
    if kind == "core":
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*full.SPHENIX_BLUE, 245), width=5)
        draw.rectangle((gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell), outline=(*full.PHOTON_DARK, 255), width=7)
    elif kind == "shoulders":
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 218), outline=(155, 164, 174, 220), width=2)
        draw.line((center_cell[0] + 10, center_cell[1] + 10, center_cell[2] - 10, center_cell[3] - 10), fill=(155, 164, 174, 200), width=3)
        draw.line((center_cell[0] + 10, center_cell[3] - 10, center_cell[2] - 10, center_cell[1] + 10), fill=(155, 164, 174, 200), width=3)
        h_span = grid_w * (0.50 if photon_like else 0.86)
        hx0 = (gx0 + gx1 - h_span) / 2
        hx1 = (gx0 + gx1 + h_span) / 2
        draw_double_arrow(draw, (hx0, gy0 + 28), (hx1, gy0 + 28), fill=full.SPHENIX_BLUE, width=6, head=15)
        v_span = grid_w * (0.48 if photon_like else 0.82)
        vx = gx1 - 26
        vy0 = (gy0 + gy1 - v_span) / 2
        vy1 = (gy0 + gy1 + v_span) / 2
        draw_double_arrow(draw, (vx, vy0), (vx, vy1), fill=full.TEAL, width=6, head=15)
    else:
        draw.rectangle((gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell), outline=(*full.SPHENIX_BLUE, 245), width=5)
        draw.rectangle((gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*full.PHOTON_DARK, 255), width=7)
    draw.ellipse((gx0 + 2.5 * cell - 10, gy0 + 2.5 * cell - 10, gx0 + 2.5 * cell + 10, gy0 + 2.5 * cell + 10), fill=(0, 0, 0, 255))
    return grid


def draw_handle_grid_v10(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    *,
    cell: int,
    kind: str,
    photon_like: bool,
) -> tuple[int, int, int, int]:
    """Tower-grid cue using non-category colors for internal geometry annotations."""
    grid = full.draw_small_tower_grid(draw, xy, cell, shower_values(kind, photon_like))
    gx0, gy0, gx1, gy1 = grid
    grid_w = gx1 - gx0
    local_edge = (93, 105, 119)
    focus_edge = (188, 119, 18)
    eta_edge = (129, 91, 164)
    phi_edge = (24, 129, 137)
    if kind == "core":
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*local_edge, 245), width=5)
        draw.rectangle((gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell), outline=(*focus_edge, 255), width=7)
    elif kind == "shoulders":
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 218), outline=(155, 164, 174, 220), width=2)
        draw.line((center_cell[0] + 10, center_cell[1] + 10, center_cell[2] - 10, center_cell[3] - 10), fill=(155, 164, 174, 200), width=3)
        draw.line((center_cell[0] + 10, center_cell[3] - 10, center_cell[2] - 10, center_cell[1] + 10), fill=(155, 164, 174, 200), width=3)
        h_span = grid_w * (0.50 if photon_like else 0.86)
        hx0 = (gx0 + gx1 - h_span) / 2
        hx1 = (gx0 + gx1 + h_span) / 2
        draw_double_arrow(draw, (hx0, gy0 + 28), (hx1, gy0 + 28), fill=eta_edge, width=6, head=15)
        v_span = grid_w * (0.48 if photon_like else 0.82)
        vx = gx1 - 26
        vy0 = (gy0 + gy1 - v_span) / 2
        vy1 = (gy0 + gy1 + v_span) / 2
        draw_double_arrow(draw, (vx, vy0), (vx, vy1), fill=phi_edge, width=6, head=15)
    else:
        draw.rectangle((gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell), outline=(*local_edge, 245), width=5)
        draw.rectangle((gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*focus_edge, 255), width=7)
    draw.ellipse((gx0 + 2.5 * cell - 10, gy0 + 2.5 * cell - 10, gx0 + 2.5 * cell + 10, gy0 + 2.5 * cell + 10), fill=(0, 0, 0, 255))
    return grid


def draw_question_box_v5(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    number: int,
    question: str,
    definition: str,
) -> None:
    draw.rounded_rectangle(box, radius=18, fill=(250, 252, 253, 255), outline=(201, 212, 224, 255), width=3)
    badge = (box[0] + 28, box[1] + 47, box[0] + 86, box[1] + 105)
    draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(124, 137, 150, 235), width=3)
    centered_text(draw, badge, str(number), full.font(full.TIMES_BOLD, 31), full.INK)
    q_font = full.font(full.TIMES_BOLD, 40)
    d_font = full.font(full.TIMES, 31)
    draw.text((box[0] + 112, box[1] + 30), question, font=q_font, fill=full.INK)
    draw.text((box[0] + 112, box[1] + 88), definition, font=d_font, fill=full.INK)


def draw_handle_row(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    spec: dict[str, object],
    *,
    show_col_labels: bool = False,
) -> None:
    x0, y0, x1, y1 = box
    accent = spec["accent"]  # type: ignore[index]
    kind = spec["kind"]  # type: ignore[index]
    draw.rounded_rectangle(box, radius=12, fill=(250, 252, 253, 255), outline=(220, 228, 235, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 8, y1), radius=4, fill=(*accent, 220))
    draw.text((x0 + 28, y0 + 22), str(spec["title"]), font=full.font(full.TIMES_BOLD, 37), fill=full.INK)
    draw_tag_pill(draw, (x0 + 30, y0 + 78), str(spec["tag"]), accent, size=26)
    draw_wrapped_left(draw, (x0 + 30, y0 + 132), str(spec["meaning"]), full.font(full.TIMES, 28), full.BLUE, 430, line_gap=4)

    cell = 32
    gy = y0 + 28
    photon_x = x0 + 640
    bkg_x = x0 + 1025
    if show_col_labels:
        for label, lx, color in [
            ("Photon-like", photon_x, full.SPHENIX_BLUE),
            ("Background-like", bkg_x, ACCENTS["fail"]),
        ]:
            fnt = full.font(full.TIMES_BOLD, 35)
            tw, _ = full.text_box(draw, label, fnt)
            draw.text((lx + (5 * cell - tw) / 2, y0 - 44), label, font=fnt, fill=color)
    draw_handle_grid(draw, (photon_x, gy), cell=cell, kind=str(kind), photon_like=True, accent=accent)  # type: ignore[arg-type]
    draw_handle_grid(draw, (bkg_x, gy), cell=cell, kind=str(kind), photon_like=False, accent=accent)  # type: ignore[arg-type]


def slide08_shower_shape_handles() -> Image.Image:
    img = draw_slide_shell(COMPRESSED_SLIDE_A_SUBTITLE, title=COMPRESSED_SLIDE_A_TITLE)
    draw = ImageDraw.Draw(img, "RGBA")

    left = (132, 304, 1708, 1282)
    right = (1742, 304, W - 132, 1282)
    panel(img, left, full.SPHENIX_BLUE, show_sidebar=False)
    panel(img, right, full.PHOTON_DARK, show_sidebar=False)

    draw.text((left[0] + 36, left[1] + 28), "Representative EMCal shower-geometry handles", font=full.font(full.TIMES_BOLD, 45), fill=full.INK)
    draw.text((left[0] + 38, left[1] + 86), "Examples, not an exhaustive list of photon-ID variables", font=full.font(full.TIMES_ITALIC, 31), fill=full.MUTED)

    row_x0 = left[0] + 34
    row_x1 = left[2] - 34
    row_y = [left[1] + 188, left[1] + 418, left[1] + 648]
    row_h = 205
    for idx, spec in enumerate(HANDLE_ROWS):
        draw_handle_row(draw, (row_x0, row_y[idx], row_x1, row_y[idx] + row_h), spec, show_col_labels=idx == 0)

    bridge = (left[0] + 36, left[3] - 74, left[2] - 36, left[3] - 24)
    draw.rounded_rectangle(bridge, radius=12, fill=(255, 248, 229, 255), outline=(236, 218, 171, 255), width=2)
    centered_text(draw, bridge, "Box cuts answer questions like these with fixed thresholds.", full.font(full.TIMES_BOLD, 34), full.BLUE)

    draw.text((right[0] + 44, right[1] + 42), "Geometry questions", font=full.font(full.TIMES_BOLD, 54), fill=full.INK)
    draw_wrapped_left(
        draw,
        (right[0] + 46, right[1] + 112),
        "The BDT score is grounded in physically interpretable shower pictures.",
        full.font(full.TIMES_ITALIC, 31),
        full.MUTED,
        right[2] - right[0] - 92,
        line_gap=5,
    )
    q_y = right[1] + 220
    for idx, spec in enumerate(HANDLE_ROWS, 1):
        qbox = (right[0] + 44, q_y + (idx - 1) * 162, right[2] - 44, q_y + (idx - 1) * 162 + 118)
        accent = spec["accent"]  # type: ignore[index]
        draw.rounded_rectangle(qbox, radius=18, fill=(248, 251, 253, 255), outline=(*accent, 210), width=3)
        badge = (qbox[0] + 24, qbox[1] + 29, qbox[0] + 84, qbox[1] + 89)
        draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(*accent, 235), width=4)
        centered_text(draw, badge, str(idx), full.font(full.TIMES_BOLD, 35), accent)
        draw.text((qbox[0] + 112, qbox[1] + 34), str(spec["question"]), font=full.font(full.TIMES_BOLD, 38), fill=full.INK)
    cue = (right[0] + 44, right[3] - 206, right[2] - 44, right[3] - 52)
    draw.rounded_rectangle(cue, radius=18, fill=(239, 246, 250, 255), outline=(207, 224, 236, 255), width=2)
    draw_lines_centered(
        draw,
        (cue[0] + 28, cue[1] + 14, cue[2] - 28, cue[3] - 14),
        [
            ("Representative handles:", full.font(full.TIMES_BOLD, 34), full.BLUE),
            ("compactness, widths, and local topology", full.font(full.TIMES, 33), full.INK),
        ],
        line_gap=8,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v2() -> Image.Image:
    """Integrated Slide 8 candidate: large shower comparisons with attached questions."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    label_font = full.font(full.TIMES_BOLD, 36)
    ph_label = "Photon-like"
    bg_label = "Background-like"
    cell = 45
    grid_w = 5 * cell
    label_y = outer[1] + 36
    photon_x = 640
    bkg_x = 1000
    ph_w, _ = full.text_box(draw, ph_label, label_font)
    bg_w, _ = full.text_box(draw, bg_label, label_font)
    draw.text((photon_x + (grid_w - ph_w) / 2, label_y), ph_label, font=label_font, fill=full.SPHENIX_BLUE)
    draw.text((bkg_x + (grid_w - bg_w) / 2, label_y), bg_label, font=label_font, fill=ACCENTS["fail"])

    row_tops = [outer[1] + 112, outer[1] + 362, outer[1] + 612]
    row_h = 226
    row_label_x = outer[0] + 58
    ribbon_x0 = 1390
    ribbon_x1 = outer[2] - 56
    for idx, (spec, y0) in enumerate(zip(HANDLE_ROWS, row_tops)):
        y1 = y0 + row_h
        accent = spec["accent"]  # type: ignore[index]
        kind = str(spec["kind"])
        if idx:
            draw.line((outer[0] + 36, y0 - 12, outer[2] - 36, y0 - 12), fill=(224, 231, 237, 255), width=2)
        draw.rounded_rectangle((outer[0] + 34, y0 + 10, outer[0] + 44, y1 - 10), radius=5, fill=(*accent, 225))
        draw.text((row_label_x, y0 + 34), str(spec["title"]), font=full.font(full.TIMES_BOLD, 40), fill=full.INK)
        draw_tag_pill(draw, (row_label_x, y0 + 94), str(spec["tag"]), accent, size=28)

        grid_y = y0 + (row_h - grid_w) // 2 + 2
        draw_handle_grid(draw, (photon_x, grid_y), cell=cell, kind=kind, photon_like=True, accent=accent)  # type: ignore[arg-type]
        draw_handle_grid(draw, (bkg_x, grid_y), cell=cell, kind=kind, photon_like=False, accent=accent)  # type: ignore[arg-type]

        arrow_y = y0 + row_h // 2
        full.draw_arrow(draw, (bkg_x + grid_w + 48, arrow_y), (ribbon_x0 - 36, arrow_y), fill=(139, 156, 173), width=5)
        qbox = (ribbon_x0, y0 + 53, ribbon_x1, y0 + 173)
        draw.rounded_rectangle(qbox, radius=18, fill=(248, 251, 253, 255), outline=(*accent, 230), width=4)
        centered_text(draw, qbox, str(spec["question"]), full.font(full.TIMES_BOLD, 43), full.INK)

    bridge = (outer[0] + 36, outer[3] - 104, outer[2] - 36, outer[3] - 28)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    mid_x = bridge[0] + (bridge[2] - bridge[0]) // 2
    draw.line((mid_x, bridge[1] + 16, mid_x, bridge[3] - 16), fill=(229, 205, 144, 255), width=2)
    draw_lines_centered(
        draw,
        (bridge[0] + 26, bridge[1] + 8, mid_x - 30, bridge[3] - 8),
        [("Example box cuts threshold handles like these independently.", full.font(full.TIMES_BOLD, 33), full.BLUE)],
    )
    draw_lines_centered(
        draw,
        (mid_x + 30, bridge[1] + 8, bridge[2] - 26, bridge[3] - 8),
        [("Correlated handles are harder to tune with independent cuts.", full.font(full.TIMES_BOLD, 33), full.BLUE)],
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v3() -> Image.Image:
    """Composition-polished Slide 8 candidate with one continuous visual argument."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((164, 340, 178, 1112), radius=7, fill=(*full.SPHENIX_BLUE, 230))

    cell = 50
    grid_w = 5 * cell
    photon_x = 555
    bkg_x = 895
    label_y = outer[1] + 34
    label_font = full.font(full.TIMES_BOLD, 38)
    ph_label = "Photon-like"
    bg_label = "Background-like"
    ph_w, _ = full.text_box(draw, ph_label, label_font)
    bg_w, _ = full.text_box(draw, bg_label, label_font)
    draw.text((photon_x + (grid_w - ph_w) / 2, label_y), ph_label, font=label_font, fill=full.SPHENIX_BLUE)
    draw.text((bkg_x + (grid_w - bg_w) / 2, label_y), bg_label, font=label_font, fill=ACCENTS["fail"])

    row_tops = [350, 602, 854]
    label_x = 214
    qbox_x0 = 1324
    qbox_x1 = outer[2] - 56
    for idx, (spec, y0) in enumerate(zip(HANDLE_ROWS, row_tops)):
        accent = spec["accent"]  # type: ignore[index]
        kind = str(spec["kind"])
        if idx:
            draw.line((outer[0] + 52, y0 - 22, outer[2] - 52, y0 - 22), fill=(226, 233, 239, 255), width=2)

        draw.text((label_x, y0 + 64), str(spec["title"]), font=full.font(full.TIMES_BOLD, 39), fill=full.INK)
        draw_tag_pill(draw, (label_x, y0 + 122), str(spec["tag"]), accent, size=28)

        grid_y = y0
        draw_handle_grid(draw, (photon_x, grid_y), cell=cell, kind=kind, photon_like=True, accent=accent)  # type: ignore[arg-type]
        draw_handle_grid(draw, (bkg_x, grid_y), cell=cell, kind=kind, photon_like=False, accent=accent)  # type: ignore[arg-type]

        arrow_y = y0 + grid_w // 2
        full.draw_arrow(draw, (bkg_x + grid_w + 32, arrow_y), (qbox_x0 - 26, arrow_y), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, y0 + 66, qbox_x1, y0 + 164)
        draw.rounded_rectangle(qbox, radius=17, fill=(248, 251, 253, 255), outline=(*accent, 230), width=4)
        centered_text(draw, qbox, str(spec["question"]), full.font(full.TIMES_BOLD, 42), full.INK)

    bridge = (outer[0] + 36, outer[3] - 100, outer[2] - 36, outer[3] - 30)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 38),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v4() -> Image.Image:
    """Slide 8 candidate with row labels removed and variable cues embedded in images."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((164, 338, 178, 1124), radius=7, fill=(*full.SPHENIX_BLUE, 230))

    cell = 50
    grid_w = 5 * cell
    photon_x = 318
    bkg_x = 700
    label_y = outer[1] + 32
    label_font = full.font(full.TIMES_BOLD, 39)
    for label, lx, color in [
        ("Photon-like", photon_x, full.SPHENIX_BLUE),
        ("Background-like", bkg_x, ACCENTS["fail"]),
    ]:
        tw, _ = full.text_box(draw, label, label_font)
        draw.text((lx + (grid_w - tw) / 2, label_y), label, font=label_font, fill=color)

    row_tops = [350, 616, 882]
    qbox_x0 = 1104
    qbox_x1 = outer[2] - 56
    for idx, (spec, y0) in enumerate(zip(HANDLE_ROWS, row_tops)):
        accent = spec["accent"]  # type: ignore[index]
        kind = str(spec["kind"])
        if idx:
            draw.line((outer[0] + 52, y0 - 22, outer[2] - 52, y0 - 22), fill=(226, 233, 239, 255), width=2)

        draw_handle_grid_v4(draw, (photon_x, y0), cell=cell, kind=kind, photon_like=True, accent=accent)  # type: ignore[arg-type]
        draw_handle_grid_v4(draw, (bkg_x, y0), cell=cell, kind=kind, photon_like=False, accent=accent)  # type: ignore[arg-type]

        arrow_y = y0 + grid_w // 2
        full.draw_arrow(draw, (bkg_x + grid_w + 34, arrow_y), (qbox_x0 - 28, arrow_y), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, y0 + 82, qbox_x1, y0 + 180)
        draw.rounded_rectangle(qbox, radius=17, fill=(248, 251, 253, 255), outline=(*accent, 225), width=4)
        centered_text(draw, qbox, str(spec["question"]), full.font(full.TIMES_BOLD, 42), full.INK)

    bridge = (outer[0] + 36, outer[3] - 98, outer[2] - 36, outer[3] - 30)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 38),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v5() -> Image.Image:
    """Slide 8 candidate with larger visual columns and variable definitions in question boxes."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    cell = 54
    grid_w = 5 * cell
    col_top = outer[1] + 26
    col_bottom = outer[3] - 120
    photon_col = (220, col_top, 540, col_bottom)
    bkg_col = (592, col_top, 912, col_bottom)
    draw.rounded_rectangle(photon_col, radius=16, fill=(254, 244, 242, 255), outline=(236, 207, 202, 255), width=2)
    draw.rounded_rectangle(bkg_col, radius=16, fill=(239, 246, 252, 255), outline=(205, 220, 234, 255), width=2)

    label_font = full.font(full.TIMES_BOLD, 39)
    for label, col in [("Photon-like", photon_col), ("Background-like", bkg_col)]:
        tw, _ = full.text_box(draw, label, label_font)
        draw.text(((col[0] + col[2] - tw) / 2, col[1] + 22), label, font=label_font, fill=full.INK)

    row_tops = [386, 642, 898]
    question_specs = [
        ("Is the core compact?", "center tower / local 3x3 energy: E11/E33"),
        ("Are the shoulders narrow?", "lateral widths around the seed: wη, wφ"),
        ("Is the shower unsplit?", "3x2 strip relative to 3x5 window: E3x2/E3x5"),
    ]
    qbox_x0 = 1044
    qbox_x1 = outer[2] - 56
    for idx, (spec, y0, qspec) in enumerate(zip(HANDLE_ROWS, row_tops, question_specs)):
        kind = str(spec["kind"])
        accent = spec["accent"]  # type: ignore[index]
        if idx:
            draw.line((outer[0] + 52, y0 - 22, outer[2] - 52, y0 - 22), fill=(226, 233, 239, 255), width=2)

        photon_x = int((photon_col[0] + photon_col[2] - grid_w) / 2)
        bkg_x = int((bkg_col[0] + bkg_col[2] - grid_w) / 2)
        draw_handle_grid_v5(draw, (photon_x, y0), cell=cell, kind=kind, photon_like=True, accent=accent)  # type: ignore[arg-type]
        draw_handle_grid_v5(draw, (bkg_x, y0), cell=cell, kind=kind, photon_like=False, accent=accent)  # type: ignore[arg-type]

        arrow_y = y0 + grid_w // 2
        full.draw_arrow(draw, (bkg_col[2] + 32, arrow_y), (qbox_x0 - 28, arrow_y), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, y0 + 58, qbox_x1, y0 + 206)
        draw_question_box_v5(draw, qbox, number=idx + 1, question=qspec[0], definition=qspec[1])

    bridge = (outer[0] + 36, outer[3] - 98, outer[2] - 36, outer[3] - 30)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 38),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v6() -> Image.Image:
    """Slide 8 candidate with separated image tiles and compact tint legend."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    cell = 54
    grid_w = 5 * cell
    tile_w, tile_h = 350, 278
    photon_tiles_x = 176
    bkg_tiles_x = 590
    row_tops = [326, 612, 898]
    qbox_x0 = 1038
    qbox_x1 = outer[2] - 56
    question_specs = [
        ("Is the core compact?", "center tower / local 3x3 energy: E11/E33"),
        ("Are the shoulders narrow?", "lateral widths around the seed: wη, wφ"),
        ("Is the shower unsplit?", "3x2 strip relative to 3x5 window: E3x2/E3x5"),
    ]

    legend = (qbox_x1 - 500, outer[1] + 28, qbox_x1, outer[1] + 86)
    draw.rounded_rectangle(legend, radius=14, fill=(255, 255, 255, 240), outline=(218, 226, 234, 255), width=2)
    sw = 34
    lx = legend[0] + 22
    ly = legend[1] + 16
    legend_font = full.font(full.TIMES_BOLD, 27)
    draw.rounded_rectangle((lx, ly, lx + sw, ly + 26), radius=7, fill=(254, 231, 227, 255), outline=(235, 195, 188, 255), width=2)
    draw.text((lx + sw + 10, legend[1] + 14), "photon-like", font=legend_font, fill=full.INK)
    lx2 = legend[0] + 248
    draw.rounded_rectangle((lx2, ly, lx2 + sw, ly + 26), radius=7, fill=(226, 239, 251, 255), outline=(190, 211, 231, 255), width=2)
    draw.text((lx2 + sw + 10, legend[1] + 14), "background-like", font=legend_font, fill=full.INK)

    for idx, (spec, y0, qspec) in enumerate(zip(HANDLE_ROWS, row_tops, question_specs)):
        kind = str(spec["kind"])
        if idx:
            draw.line((outer[0] + 52, y0 - 16, outer[2] - 52, y0 - 16), fill=(226, 233, 239, 255), width=2)

        photon_tile = (photon_tiles_x, y0, photon_tiles_x + tile_w, y0 + tile_h)
        bkg_tile = (bkg_tiles_x, y0, bkg_tiles_x + tile_w, y0 + tile_h)
        draw.rounded_rectangle(photon_tile, radius=16, fill=(254, 244, 242, 255), outline=(236, 207, 202, 255), width=2)
        draw.rounded_rectangle(bkg_tile, radius=16, fill=(239, 246, 252, 255), outline=(205, 220, 234, 255), width=2)

        grid_y = y0 + (tile_h - grid_w) // 2
        photon_x = photon_tile[0] + (tile_w - grid_w) // 2
        bkg_x = bkg_tile[0] + (tile_w - grid_w) // 2
        draw_handle_grid_v5(draw, (photon_x, grid_y), cell=cell, kind=kind, photon_like=True, accent=full.PHOTON_DARK)  # type: ignore[arg-type]
        draw_handle_grid_v5(draw, (bkg_x, grid_y), cell=cell, kind=kind, photon_like=False, accent=full.PHOTON_DARK)  # type: ignore[arg-type]

        arrow_y = y0 + tile_h // 2
        full.draw_arrow(draw, (bkg_tile[2] + 30, arrow_y), (qbox_x0 - 30, arrow_y), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, y0 + 68, qbox_x1, y0 + 216)
        draw_question_box_v5(draw, qbox, number=idx + 1, question=qspec[0], definition=qspec[1])

    bridge = (outer[0] + 36, outer[3] - 98, outer[2] - 36, outer[3] - 30)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 38),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v7() -> Image.Image:
    """Slide 8 candidate with larger symmetric comparison columns and a compact column key."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    cell = 56
    grid_w = 5 * cell
    tile_w, tile_h = 370, 290
    photon_tiles_x = 160
    bkg_tiles_x = 570
    row_tops = [306, 604, 902]
    qbox_x0 = 1018
    qbox_x1 = outer[2] - 56
    q_h = 154
    question_specs = [
        ("Is the core compact?", "center tower / local 3x3 energy: E11/E33"),
        ("Are the shoulders narrow?", "lateral widths around the seed: wη, wφ"),
        ("Is the shower unsplit?", "3x2 strip relative to 3x5 window: E3x2/E3x5"),
    ]

    key = (qbox_x1 - 575, outer[1] + 24, qbox_x1, outer[1] + 82)
    draw.rounded_rectangle(key, radius=14, fill=(255, 255, 255, 242), outline=(218, 226, 234, 255), width=2)
    key_font = full.font(full.TIMES_BOLD, 27)
    sw = 34
    ky = key[1] + 16
    kx = key[0] + 22
    draw.rounded_rectangle((kx, ky, kx + sw, ky + 26), radius=7, fill=(254, 231, 227, 255), outline=(235, 195, 188, 255), width=2)
    draw.text((kx + sw + 10, key[1] + 14), "left: photon-like", font=key_font, fill=full.INK)
    kx2 = key[0] + 286
    draw.rounded_rectangle((kx2, ky, kx2 + sw, ky + 26), radius=7, fill=(226, 239, 251, 255), outline=(190, 211, 231, 255), width=2)
    draw.text((kx2 + sw + 10, key[1] + 14), "right: background-like", font=key_font, fill=full.INK)

    for idx, (spec, y0, qspec) in enumerate(zip(HANDLE_ROWS, row_tops, question_specs)):
        kind = str(spec["kind"])
        photon_tile = (photon_tiles_x, y0, photon_tiles_x + tile_w, y0 + tile_h)
        bkg_tile = (bkg_tiles_x, y0, bkg_tiles_x + tile_w, y0 + tile_h)
        draw.rounded_rectangle(photon_tile, radius=17, fill=(254, 244, 242, 255), outline=(235, 199, 193, 255), width=2)
        draw.rounded_rectangle(bkg_tile, radius=17, fill=(239, 246, 252, 255), outline=(197, 216, 234, 255), width=2)

        grid_y = y0 + (tile_h - grid_w) // 2
        photon_x = photon_tile[0] + (tile_w - grid_w) // 2
        bkg_x = bkg_tile[0] + (tile_w - grid_w) // 2
        draw_handle_grid_v5(draw, (photon_x, grid_y), cell=cell, kind=kind, photon_like=True, accent=full.PHOTON_DARK)  # type: ignore[arg-type]
        draw_handle_grid_v5(draw, (bkg_x, grid_y), cell=cell, kind=kind, photon_like=False, accent=full.PHOTON_DARK)  # type: ignore[arg-type]

        row_center = y0 + tile_h // 2
        full.draw_arrow(draw, (bkg_tile[2] + 26, row_center), (qbox_x0 - 28, row_center), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, row_center - q_h // 2, qbox_x1, row_center + q_h // 2)
        draw_question_box_v5(draw, qbox, number=idx + 1, question=qspec[0], definition=qspec[1])

    bridge = (outer[0] + 36, outer[3] - 82, outer[2] - 36, outer[3] - 24)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 37),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v8() -> Image.Image:
    """Slide 8 candidate with white comparison tiles, thick category borders, and local column key."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    cell = 54
    grid_w = 5 * cell
    tile_w, tile_h = 368, 282
    photon_tiles_x = 154
    bkg_tiles_x = 574
    row_tops = [346, 630, 914]
    qbox_x0 = 1024
    qbox_x1 = outer[2] - 56
    q_h = 154
    photon_edge = (221, 89, 75)
    bkg_edge = (51, 141, 207)
    question_specs = [
        ("Is the core compact?", "center tower / local 3x3 energy: E11/E33"),
        ("Are the shoulders narrow?", "lateral widths around the seed: wη, wφ"),
        ("Is the shower unsplit?", "3x2 strip relative to 3x5 window: E3x2/E3x5"),
    ]

    key_font = full.font(full.TIMES_BOLD, 31)
    key_y0, key_y1 = outer[1] + 24, outer[1] + 76
    for label, x0, edge in [
        ("photon-like", photon_tiles_x, photon_edge),
        ("background-like", bkg_tiles_x, bkg_edge),
    ]:
        key_box = (x0 + 12, key_y0, x0 + tile_w - 12, key_y1)
        draw.rounded_rectangle(key_box, radius=16, fill=(255, 255, 255, 245), outline=(*edge, 255), width=4)
        centered_text(draw, key_box, label, key_font, full.INK)

    for idx, (spec, y0, qspec) in enumerate(zip(HANDLE_ROWS, row_tops, question_specs)):
        kind = str(spec["kind"])
        photon_tile = (photon_tiles_x, y0, photon_tiles_x + tile_w, y0 + tile_h)
        bkg_tile = (bkg_tiles_x, y0, bkg_tiles_x + tile_w, y0 + tile_h)
        draw.rounded_rectangle(photon_tile, radius=18, fill=(255, 255, 255, 252), outline=(*photon_edge, 250), width=6)
        draw.rounded_rectangle(bkg_tile, radius=18, fill=(255, 255, 255, 252), outline=(*bkg_edge, 250), width=6)

        grid_y = y0 + (tile_h - grid_w) // 2
        photon_x = photon_tile[0] + (tile_w - grid_w) // 2
        bkg_x = bkg_tile[0] + (tile_w - grid_w) // 2
        draw_handle_grid_v5(draw, (photon_x, grid_y), cell=cell, kind=kind, photon_like=True, accent=full.PHOTON_DARK)  # type: ignore[arg-type]
        draw_handle_grid_v5(draw, (bkg_x, grid_y), cell=cell, kind=kind, photon_like=False, accent=full.PHOTON_DARK)  # type: ignore[arg-type]

        row_center = y0 + tile_h // 2
        full.draw_arrow(draw, (bkg_tile[2] + 26, row_center), (qbox_x0 - 28, row_center), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, row_center - q_h // 2, qbox_x1, row_center + q_h // 2)
        draw_question_box_v5(draw, qbox, number=idx + 1, question=qspec[0], definition=qspec[1])

    bridge = (outer[0] + 36, outer[3] - 70, outer[2] - 36, outer[3] - 18)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 36),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v9() -> Image.Image:
    """Slide 8 candidate with category color on the outermost grid boundary itself."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    cell = 58
    grid_w = 5 * cell
    photon_x = 190
    bkg_x = 585
    row_tops = [320, 610, 900]
    qbox_x0 = 1018
    qbox_x1 = outer[2] - 56
    q_h = 154
    photon_edge = (221, 89, 75)
    bkg_edge = (51, 141, 207)
    question_specs = [
        ("Is the core compact?", "center tower / local 3x3 energy: E11/E33"),
        ("Are the shoulders narrow?", "lateral widths around the seed: wη, wφ"),
        ("Is the shower unsplit?", "3x2 strip relative to 3x5 window: E3x2/E3x5"),
    ]

    for idx, (spec, y0, qspec) in enumerate(zip(HANDLE_ROWS, row_tops, question_specs)):
        kind = str(spec["kind"])
        p_grid = draw_handle_grid_v5(draw, (photon_x, y0), cell=cell, kind=kind, photon_like=True, accent=full.PHOTON_DARK)  # type: ignore[arg-type]
        b_grid = draw_handle_grid_v5(draw, (bkg_x, y0), cell=cell, kind=kind, photon_like=False, accent=full.PHOTON_DARK)  # type: ignore[arg-type]
        draw.rectangle(p_grid, outline=(*photon_edge, 255), width=7)
        draw.rectangle(b_grid, outline=(*bkg_edge, 255), width=7)

        row_center = y0 + grid_w // 2
        full.draw_arrow(draw, (b_grid[2] + 28, row_center), (qbox_x0 - 28, row_center), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, row_center - q_h // 2, qbox_x1, row_center + q_h // 2)
        draw_question_box_v5(draw, qbox, number=idx + 1, question=qspec[0], definition=qspec[1])

    bridge = (outer[0] + 36, outer[3] - 70, outer[2] - 36, outer[3] - 18)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 36),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v10() -> Image.Image:
    """Slide 8 candidate with separated tinted columns and category color only on grid borders."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    cell = 56
    grid_w = 5 * cell
    photon_x = 184
    bkg_x = 568
    row_tops = [322, 616, 910]
    qbox_x0 = 1010
    qbox_x1 = outer[2] - 56
    q_h = 154
    photon_edge = (221, 89, 75)
    bkg_edge = (51, 141, 207)
    question_specs = [
        ("Is the core compact?", "center tower / local 3x3 energy: E11/E33"),
        ("Are the shoulders narrow?", "lateral widths around the seed: wη, wφ"),
        ("Is the shower unsplit?", "3x2 strip relative to 3x5 window: E3x2/E3x5"),
    ]

    band_top, band_bottom = outer[1] + 28, outer[3] - 80
    photon_band = (outer[0] + 38, band_top, photon_x + grid_w + 38, band_bottom)
    bkg_band = (bkg_x - 38, band_top, bkg_x + grid_w + 38, band_bottom)
    question_band = (qbox_x0 - 36, band_top, qbox_x1 + 18, band_bottom)
    draw.rounded_rectangle(photon_band, radius=18, fill=(255, 249, 248, 255), outline=(247, 224, 221, 255), width=2)
    draw.rounded_rectangle(bkg_band, radius=18, fill=(247, 251, 255, 255), outline=(220, 233, 245, 255), width=2)
    draw.rounded_rectangle(question_band, radius=18, fill=(249, 250, 248, 255), outline=(224, 229, 232, 255), width=2)

    for idx, (spec, y0, qspec) in enumerate(zip(HANDLE_ROWS, row_tops, question_specs)):
        kind = str(spec["kind"])
        p_grid = draw_handle_grid_v10(draw, (photon_x, y0), cell=cell, kind=kind, photon_like=True)  # type: ignore[arg-type]
        b_grid = draw_handle_grid_v10(draw, (bkg_x, y0), cell=cell, kind=kind, photon_like=False)  # type: ignore[arg-type]
        draw.rectangle(p_grid, outline=(*photon_edge, 255), width=8)
        draw.rectangle(b_grid, outline=(*bkg_edge, 255), width=8)

        row_center = y0 + grid_w // 2
        full.draw_arrow(draw, (b_grid[2] + 30, row_center), (qbox_x0 - 30, row_center), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, row_center - q_h // 2, qbox_x1, row_center + q_h // 2)
        draw_question_box_v5(draw, qbox, number=idx + 1, question=qspec[0], definition=qspec[1])

    bridge = (outer[0] + 36, outer[3] - 70, outer[2] - 36, outer[3] - 18)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 36),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide08_shower_shape_handles_v11() -> Image.Image:
    """Slide 8 candidate with variable lanes and category color only on grid boundaries."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    cell = 56
    grid_w = 5 * cell
    photon_x = 184
    bkg_x = 568
    row_tops = [330, 620, 910]
    qbox_x0 = 1010
    qbox_x1 = outer[2] - 56
    q_h = 154
    photon_edge = (221, 89, 75)
    bkg_edge = (51, 141, 207)
    lane_specs = [
        ((255, 251, 241, 255), (239, 225, 189, 255), "Is the core compact?", "center tower / local 3x3 energy: E11/E33"),
        ((244, 251, 249, 255), (205, 228, 224, 255), "Are the shoulders narrow?", "lateral widths around the seed: wη, wφ"),
        ((249, 248, 253, 255), (222, 216, 235, 255), "Is the shower unsplit?", "3x2 strip relative to 3x5 window: E3x2/E3x5"),
    ]

    for idx, (spec, y0, lane) in enumerate(zip(HANDLE_ROWS, row_tops, lane_specs)):
        lane_fill, lane_edge, question, definition = lane
        kind = str(spec["kind"])
        lane_box = (outer[0] + 36, y0 - 6, outer[2] - 36, y0 + grid_w + 6)
        draw.rounded_rectangle(lane_box, radius=18, fill=lane_fill, outline=lane_edge, width=2)

        p_grid = draw_handle_grid_v10(draw, (photon_x, y0), cell=cell, kind=kind, photon_like=True)  # type: ignore[arg-type]
        b_grid = draw_handle_grid_v10(draw, (bkg_x, y0), cell=cell, kind=kind, photon_like=False)  # type: ignore[arg-type]
        draw.rectangle(p_grid, outline=(*photon_edge, 255), width=8)
        draw.rectangle(b_grid, outline=(*bkg_edge, 255), width=8)

        row_center = y0 + grid_w // 2
        full.draw_arrow(draw, (b_grid[2] + 30, row_center), (qbox_x0 - 30, row_center), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, row_center - q_h // 2, qbox_x1, row_center + q_h // 2)
        draw_question_box_v5(draw, qbox, number=idx + 1, question=question, definition=definition)

    bridge = (outer[0] + 36, outer[3] - 70, outer[2] - 36, outer[3] - 18)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 36),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def redraw_unsplit_regions(
    draw: ImageDraw.ImageDraw,
    grid: tuple[int, int, int, int],
    *,
    cell: int,
) -> None:
    """Redraw split-window annotations after the category border so they stay on top."""
    gx0, gy0, _, _ = grid
    local_edge = (93, 105, 119)
    focus_edge = (188, 119, 18)
    draw.rectangle((gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell), outline=(*local_edge, 255), width=6)
    draw.rectangle((gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*focus_edge, 255), width=8)


def slide08_shower_shape_handles_v12() -> Image.Image:
    """Slide 8 candidate with full-width variable lanes and strict row geometry."""
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 78), COMPRESSED_SLIDE_A_V2_TITLE, font=full.font(full.TIMES_BOLD, 82), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232, 255), width=3)

    outer = (132, 282, W - 132, 1280)
    draw.rounded_rectangle(outer, radius=14, fill=(255, 255, 255, 250), outline=(*full.PANEL_EDGE, 255), width=2)

    cell = 56
    grid_w = 5 * cell
    lane_h = 288
    lane_tops = [304, 602, 900]
    grid_offset_y = 4
    photon_x = 172
    bkg_x = 556
    qbox_x0 = 998
    qbox_x1 = outer[2] - 56
    q_h = 154
    photon_edge = (221, 89, 75)
    bkg_edge = (51, 141, 207)
    lane_specs = [
        ((255, 251, 241, 255), (239, 225, 189, 255), "Is the core compact?", "center tower / local 3x3 energy: E11/E33"),
        ((244, 251, 249, 255), (205, 228, 224, 255), "Are the shoulders narrow?", "lateral widths around the seed: wη, wφ"),
        ((249, 248, 253, 255), (222, 216, 235, 255), "Is the shower unsplit?", "3x2 strip relative to 3x5 window: E3x2/E3x5"),
    ]

    for idx, (spec, lane_top, lane) in enumerate(zip(HANDLE_ROWS, lane_tops, lane_specs)):
        lane_fill, lane_edge, question, definition = lane
        kind = str(spec["kind"])
        lane_box = (outer[0] + 36, lane_top, outer[2] - 36, lane_top + lane_h)
        draw.rounded_rectangle(lane_box, radius=18, fill=lane_fill, outline=lane_edge, width=2)

        grid_y = lane_top + grid_offset_y
        p_grid = draw_handle_grid_v10(draw, (photon_x, grid_y), cell=cell, kind=kind, photon_like=True)  # type: ignore[arg-type]
        b_grid = draw_handle_grid_v10(draw, (bkg_x, grid_y), cell=cell, kind=kind, photon_like=False)  # type: ignore[arg-type]
        draw.rectangle(p_grid, outline=(*photon_edge, 255), width=8)
        draw.rectangle(b_grid, outline=(*bkg_edge, 255), width=8)
        if kind == "split":
            redraw_unsplit_regions(draw, p_grid, cell=cell)
            redraw_unsplit_regions(draw, b_grid, cell=cell)

        row_center = grid_y + grid_w // 2
        full.draw_arrow(draw, (b_grid[2] + 30, row_center), (qbox_x0 - 30, row_center), fill=(132, 150, 168), width=5)
        qbox = (qbox_x0, row_center - q_h // 2, qbox_x1, row_center + q_h // 2)
        draw_question_box_v5(draw, qbox, number=idx + 1, question=question, definition=definition)

    bridge = (outer[0] + 36, outer[3] - 70, outer[2] - 36, outer[3] - 18)
    draw.rounded_rectangle(bridge, radius=16, fill=(255, 248, 229, 255), outline=(235, 216, 166, 255), width=2)
    centered_text(
        draw,
        bridge,
        "Box cuts threshold handles independently; correlated handles motivate the BDT.",
        full.font(full.TIMES_BOLD, 36),
        full.INK,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def slide09_cut_logic_to_bdt_score() -> Image.Image:
    img = slide11_simplified_bridge()
    draw = ImageDraw.Draw(img, "RGBA")
    # Repaint the header text only; the body layout is the cleaned vertical progression.
    draw.rectangle((0, 46, W, 245), fill=(*full.SOFT_BG, 255))
    draw.text(tuple(HP2026_MAIN_HEADER["title_xy"]), COMPRESSED_SLIDE_B_TITLE, font=full.font(full.TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)
    return img


def slide11_simplified_bridge() -> Image.Image:
    img = draw_slide_shell(SIMPLIFIED_SUBTITLE, title=SIMPLIFIED_TITLE)
    draw = ImageDraw.Draw(img, "RGBA")

    left_margin, right_margin = 132, W - 132
    # Top band: fixed cuts as independent gates.
    top = (left_margin, 304, right_margin, 500)
    full.shadow(img, top)
    draw.rounded_rectangle(top, radius=13, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((top[0], top[1], top[0] + 12, top[3]), radius=6, fill=(*full.PHOTON_DARK, 230))
    draw.text((top[0] + 42, top[1] + 26), "Fixed cuts test each handle separately", font=full.font(full.TIMES_BOLD, 43), fill=full.INK)
    gates_y = top[1] + 86
    gate_w = 415
    gate_gap = 28
    gate_x = top[0] + 46
    for i, (label, accent) in enumerate([
        ("compact core?", full.PHOTON_DARK),
        ("narrow shoulders?", full.SPHENIX_BLUE),
        ("not split?", full.TEAL),
    ]):
        draw_simple_gate(draw, (gate_x + i * (gate_w + gate_gap), gates_y, gate_x + i * (gate_w + gate_gap) + gate_w, gates_y + 94), label, accent)
    takeaway = (top[0] + 1460, top[1] + 72, top[2] - 44, top[3] - 38)
    draw.rounded_rectangle(takeaway, radius=16, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    draw_lines_centered(
        draw,
        (takeaway[0] + 28, takeaway[1] + 10, takeaway[2] - 28, takeaway[3] - 10),
        [
            ("Transparent, but rigid:", full.font(full.TIMES_BOLD, 31), full.BLUE),
            ("independent thresholds", full.font(full.TIMES_BOLD, 31), full.BLUE),
        ],
        line_gap=5,
    )

    # Middle card: the decision tree remains the visual bridge.
    mid = (left_margin, 524, right_margin, 1024)
    full.shadow(img, mid)
    draw.rounded_rectangle(mid, radius=13, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((mid[0], mid[1], mid[0] + 12, mid[3]), radius=6, fill=(*full.SPHENIX_BLUE, 230))
    centered_text(draw, (mid[0] + 42, mid[1] + 24, mid[2] - 42, mid[1] + 82), "One learned decision tree", full.font(full.TIMES_BOLD, 48), full.INK)
    draw_bridge_tree(draw, (mid[0] + 170, mid[1] + 94, mid[2] - 170, mid[3] - 58))
    tree_note = (mid[2] - 680, mid[1] + 28, mid[2] - 48, mid[1] + 116)
    draw.rounded_rectangle(tree_note, radius=14, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_lines_centered(
        draw,
        (tree_note[0] + 22, tree_note[1] + 8, tree_note[2] - 22, tree_note[3] - 8),
        [("Same inputs, learned conditional order.", full.font(full.TIMES_BOLD, 32), full.BLUE)],
    )

    # Bottom band: tree ensemble to one score.
    bot = (left_margin, 1048, right_margin, 1282)
    full.shadow(img, bot)
    draw.rounded_rectangle(bot, radius=13, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((bot[0], bot[1], bot[0] + 12, bot[3]), radius=6, fill=(*full.TEAL, 225))
    draw.text((bot[0] + 42, bot[1] + 24), "Boosting turns many trees into one score", font=full.font(full.TIMES_BOLD, 43), fill=full.INK)
    icon_y = bot[1] + 94
    icon_xs = [bot[0] + 450, bot[0] + 620, bot[0] + 790, bot[0] + 960]
    for idx, ix in enumerate(icon_xs):
        draw_mini_tree_custom(draw, (ix, icon_y), scale=0.82, alpha=230)
        if idx < len(icon_xs) - 1:
            draw.text((ix + 84, icon_y + 58), "+", font=full.font(full.TIMES_BOLD, 44), fill=full.LIGHT_MUTED)
    full.draw_arrow(draw, (bot[0] + 1140, bot[1] + 142), (bot[0] + 1310, bot[1] + 142), fill=(126, 145, 164), width=7)
    score = (bot[0] + 1345, bot[1] + 94, bot[2] - 350, bot[1] + 178)
    draw_score_bar(draw, score)
    score_note = (bot[2] - 318, bot[1] + 58, bot[2] - 44, bot[3] - 44)
    draw.rounded_rectangle(score_note, radius=14, fill=(238, 247, 247, 255), outline=(176, 211, 214, 255), width=2)
    draw_lines_centered(
        draw,
        (score_note[0] + 20, score_note[1] + 8, score_note[2] - 20, score_note[3] - 8),
        [
            ("A stable ensemble", full.font(full.TIMES_BOLD, 29), full.BLUE),
            ("ranks each candidate", full.font(full.TIMES_BOLD, 29), full.BLUE),
            ("photon-ID score.", full.font(full.TIMES_BOLD, 29), full.BLUE),
        ],
        line_gap=3,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def make_contact_sheet(
    paths: list[Path],
    *,
    name: str = "slide11_13_bdt_yesno_sequence_contact_sheet.png",
    labels: tuple[str, ...] | None = None,
) -> Path:
    thumb_w, thumb_h = 768, 432
    label_h = 86
    cols = max(1, len(paths))
    sheet = Image.new("RGB", (cols * thumb_w, label_h + thumb_h), (246, 249, 252))
    draw = ImageDraw.Draw(sheet, "RGBA")
    if labels is None:
        labels = tuple(f"slide {idx + 1}" for idx in range(len(paths)))
    for col, (path, label) in enumerate(zip(paths, labels)):
        x = col * thumb_w
        lf = full.font(full.TIMES_BOLD, 30)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((x + (thumb_w - tw) / 2, 28), label, font=lf, fill=full.BLUE)
        thumb = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(thumb, (x, label_h))
        draw.rectangle((x, label_h, x + thumb_w - 1, label_h + thumb_h - 1), outline=(204, 216, 228, 255), width=2)
    out = OUT_DIR / name
    sheet.save(out, "PNG")
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    slides = [
        ("slide08_shower_shape_interpretable_handles.png", slide08_shower_shape_handles),
        ("slide08_shower_shape_interpretable_handles_v2.png", slide08_shower_shape_handles_v2),
        ("slide08_shower_shape_interpretable_handles_v3.png", slide08_shower_shape_handles_v3),
        ("slide08_shower_shape_interpretable_handles_v4.png", slide08_shower_shape_handles_v4),
        ("slide08_shower_shape_interpretable_handles_v5.png", slide08_shower_shape_handles_v5),
        ("slide08_shower_shape_interpretable_handles_v6.png", slide08_shower_shape_handles_v6),
        ("slide08_shower_shape_interpretable_handles_v7.png", slide08_shower_shape_handles_v7),
        ("slide08_shower_shape_interpretable_handles_v8.png", slide08_shower_shape_handles_v8),
        ("slide08_shower_shape_interpretable_handles_v9.png", slide08_shower_shape_handles_v9),
        ("slide08_shower_shape_interpretable_handles_v10.png", slide08_shower_shape_handles_v10),
        ("slide08_shower_shape_interpretable_handles_v11.png", slide08_shower_shape_handles_v11),
        ("slide08_shower_shape_interpretable_handles_v12.png", slide08_shower_shape_handles_v12),
        ("slide09_cut_logic_to_bdt_score.png", slide09_cut_logic_to_bdt_score),
        ("slide11_bdt_fixed_cuts_to_photon_id_score_simplified.png", slide11_simplified_bridge),
        ("slide11_bdt_manual_rectangular_cuts.png", slide11_manual),
        ("slide12_bdt_funnel_to_single_tree.png", slide12_tree),
        ("slide13_bdt_boosted_ensemble_score.png", slide13_boost),
    ]
    paths: list[Path] = []
    for name, fn in slides:
        path = OUT_DIR / name
        fn().convert("RGB").save(path, "PNG")
        paths.append(path)
    contact = make_contact_sheet(
        [paths[11], paths[12]],
        name="slide08_09_shower_bdt_compression_contact_sheet.png",
        labels=("8 / shower handles", "9 / cut logic to BDT"),
    )
    backup_contact = make_contact_sheet(
        paths[14:],
        name="slide11_13_bdt_yesno_sequence_contact_sheet.png",
        labels=("11 / fixed gates", "12 / funnel to tree", "13 / boosted score"),
    )
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "source_generator": str(Path(__file__).resolve()),
        "style_source": str((ROOT / "scripts/slides/hp2026/fulltalk/make_hp2026_fulltalk_candidates.py").resolve()),
        "design_contract": "Compress old main Slides 8-13 into two main-talk slides: representative shower-shape handles -> fixed cuts/tree/BDT score, with detailed originals available for backup.",
        "hp2026_main_header": HP2026_MAIN_HEADER,
        "outputs": [str(p) for p in paths],
        "main_deck_replacements": [str(paths[11]), str(paths[12])],
        "slide08_review_candidate_v12": str(paths[11]),
        "slide08_review_candidate_v11": str(paths[10]),
        "slide08_review_candidate_v10": str(paths[9]),
        "slide08_review_candidate_v9": str(paths[8]),
        "slide08_review_candidate_v8": str(paths[7]),
        "slide08_review_candidate_v7": str(paths[6]),
        "slide08_review_candidate_v6": str(paths[5]),
        "slide08_review_candidate_v5": str(paths[4]),
        "slide08_review_candidate_v4": str(paths[3]),
        "slide08_review_candidate_v3": str(paths[2]),
        "slide08_review_candidate_v2": str(paths[1]),
        "superseded_slide08_candidate": str(paths[0]),
        "superseded_one_slide_replacement": str(paths[13]),
        "backup_detail_outputs": [str(p) for p in paths[14:]],
        "contact_sheet": str(contact),
        "backup_contact_sheet": str(backup_contact),
        "size": [W, H],
    }
    script_path = OUT_DIR / "slide11_13_bdt_yesno_sequence_script.md"
    script_path.write_text(SCRIPT_MD, encoding="utf-8")
    simplified_script_path = OUT_DIR / "slide11_bdt_fixed_cuts_to_photon_id_score_simplified_script.md"
    simplified_script_path.write_text(SIMPLIFIED_SCRIPT_MD, encoding="utf-8")
    compressed_script_path = OUT_DIR / "slide08_09_shower_bdt_compression_script.md"
    compressed_script_path.write_text(COMPRESSED_SCRIPT_MD, encoding="utf-8")
    manifest["script"] = str(script_path)
    manifest["main_deck_replacement_script"] = str(simplified_script_path)
    manifest["two_slide_compression_script"] = str(compressed_script_path)
    manifest_path = OUT_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(manifest_path)
    for path in paths:
        print(path)
    print(contact)
    print(backup_contact)


if __name__ == "__main__":
    main()
