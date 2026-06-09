#!/usr/bin/env python3
"""Local no-blur strategy prototypes for the HP2026 shower-shape build.

This is a local design comparison only. It imports the current HP2026
full-talk drawing primitives and writes PNG previews; it does not mutate the
Google Slides deck or the primary full-talk generator.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_fulltalk_candidates as full


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260608-slides8-10-no-blur-variants"
W, H = full.W, full.H

PANELS = [
    {
        "short": "Core compactness",
        "title": "Is the core compact?",
        "next": "next check",
        "takeaway": "Core energy stays concentrated near the seed.",
        "accent": full.PHOTON_DARK,
        "draw": full.draw_compact_core_panel,
    },
    {
        "short": "Shoulder width",
        "title": "Are the shoulders narrow?",
        "next": "next check",
        "takeaway": "Seed-excluded shoulders stay narrow.",
        "accent": full.SPHENIX_BLUE,
        "draw": full.draw_narrow_shoulders_panel,
    },
    {
        "short": "Stretch / split",
        "title": "Is the shower stretched or split?",
        "next": "final check",
        "takeaway": "Energy remains confined rather than elongated or split.",
        "accent": full.TEAL,
        "draw": full.draw_elongation_panel,
    },
]

SUBTITLES = [
    "First, check whether the EMCal energy is concentrated in the local core.",
    "Next, check whether the surrounding energy remains narrow around the core.",
    "A photon-like EMCal cluster has a concentrated core, narrow shoulders, and little evidence of elongation or splitting.",
]

BOTTOM_STEP_TEXT = [
    "Core compactness asks whether the candidate starts as one concentrated photon-like deposit.",
    "Adding shoulder width checks whether nearby energy stays narrow after the seed tower is removed.",
    "Together, these variables ask whether the EMCal energy looks like one compact photon, rather than a broad, elongated, or split decay-like shower.",
]


def draw_concept_badge(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    index: int,
    *,
    size: int = 66,
    offset: tuple[int, int] = (36, 34),
    bar: bool = True,
) -> None:
    """Draw a consistent numbered concept badge for the shower-shape sequence."""
    x0, y0, _, y1 = box
    accent = PANELS[index]["accent"]
    if bar:
        draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 220))
    cx0 = x0 + offset[0]
    cy0 = y0 + offset[1]
    chip = (cx0, cy0, cx0 + size, cy0 + size)
    draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 235), width=4)
    num = str(index + 1)
    nf = full.font(full.TIMES_BOLD, int(size * 0.55))
    bbox = draw.textbbox((0, 0), num, font=nf)
    cx = (chip[0] + chip[2]) / 2
    cy = (chip[1] + chip[3]) / 2
    draw.text((cx - (bbox[0] + bbox[2]) / 2, cy - (bbox[1] + bbox[3]) / 2), num, font=nf, fill=accent)


def draw_core_fraction_eq(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    size: int,
    *,
    prefix: str = "Core fraction = ",
    fill: tuple[int, int, int] = full.INK,
    prefix_bold: bool = True,
) -> int:
    small = max(15, round(size * 0.58))
    sub = max(14, round(size * 0.64))
    parts: list[tuple[str, int, int, Path]] = []
    if prefix:
        parts.append((prefix, size, 0, full.TIMES_BOLD if prefix_bold else full.TIMES))
    parts.extend(
        [
            ("E", size, 0, full.TIMES_BOLD),
            ("2x2", small, round(size * 0.34), full.TIMES_BOLD),
            (" / ", size, 0, full.TIMES),
            ("E", size, 0, full.TIMES_BOLD),
            ("cluster", small, round(size * 0.34), full.TIMES_BOLD),
        ]
    )
    return full.draw_formula_run(draw, xy, parts, fill=fill)


def draw_e1x1_e3x3_eq(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    size: int,
    *,
    prefix: str = "",
    fill: tuple[int, int, int] = full.INK,
) -> int:
    sub = max(14, round(size * 0.62))
    parts: list[tuple[str, int, int, Path]] = []
    if prefix:
        parts.append((prefix, size, 0, full.TIMES_BOLD))
    parts.extend(
        [
            ("E", size, 0, full.TIMES_BOLD),
            ("1x1", sub, round(size * 0.34), full.TIMES_BOLD),
            (" / E", size, 0, full.TIMES_BOLD),
            ("3x3", sub, round(size * 0.34), full.TIMES_BOLD),
        ]
    )
    return full.draw_formula_run(draw, xy, parts, fill=fill)


def draw_e3x2_e3x5_eq(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    size: int,
    *,
    fill: tuple[int, int, int] = full.INK,
) -> int:
    sub = max(14, round(size * 0.62))
    return full.draw_formula_run(
        draw,
        xy,
        [
            ("E", size, 0, full.TIMES_BOLD),
            ("3x2", sub, round(size * 0.34), full.TIMES_BOLD),
            (" / E", size, 0, full.TIMES_BOLD),
            ("3x5", sub, round(size * 0.34), full.TIMES_BOLD),
            (" = strip fraction", size, 0, full.TIMES),
        ],
        fill=fill,
    )


def panel_boxes() -> list[tuple[int, int, int, int]]:
    panel_top = 324
    panel_bottom = 1138
    gap = 38
    panel_w = (W - 2 * 132 - 2 * gap) // 3
    return [
        (132, panel_top, 132 + panel_w, panel_bottom),
        (132 + panel_w + gap, panel_top, 132 + 2 * panel_w + gap, panel_bottom),
        (132 + 2 * (panel_w + gap), panel_top, W - 132, panel_bottom),
    ]


def new_base(frame: int, subtitle: str | None = None) -> Image.Image:
    img = full.base_slide(
        "What the shower-shape variables measure",
        subtitle or SUBTITLES[frame],
    )
    full.add_top_right_sphenix_logo_like_slide2(img)
    full.draw_header_shower_icon(img)
    return img


def bottom_takeaway(img: Image.Image, text: str) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rounded_rectangle(
        (132, 1192, W - 132, 1308),
        radius=8,
        fill=(239, 246, 250, 255),
        outline=(213, 226, 235, 255),
        width=2,
    )
    full.draw_wrapped(
        draw,
        text,
        (176, 1218),
        W - 352,
        full.font(full.TIMES_ITALIC, 34),
        fill=full.BLUE,
        line_gap=6,
    )


def draw_placeholder(
    img: Image.Image,
    box: tuple[int, int, int, int],
    index: int,
    label: str = "next check",
    dim_completed: bool = False,
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    accent = PANELS[index]["accent"]
    fill = (248, 251, 253, 255) if not dim_completed else (252, 253, 254, 255)
    outline = (213, 226, 235, 255)
    draw.rounded_rectangle(box, radius=12, fill=fill, outline=outline, width=2)
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 220))

    chip = (x0 + 42, y0 + 38, x0 + 106, y0 + 102)
    draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 235), width=3)
    num = str(index + 1)
    nf = full.font(full.TIMES_BOLD, 38)
    tw, th = full.text_box(draw, num, nf)
    draw.text(
        (chip[0] + (chip[2] - chip[0] - tw) / 2, chip[1] + (chip[3] - chip[1] - th) / 2 - 1),
        num,
        font=nf,
        fill=accent,
    )

    title_font = full.font(full.TIMES_BOLD, 35)
    muted = (81, 91, 105)
    draw.text((x0 + 132, y0 + 45), f"Step {index + 1}", font=title_font, fill=muted)
    body = "Revealed on next click" if label == "next check" else label
    full.draw_wrapped(
        draw,
        body,
        (x0 + 132, y0 + 96),
        x1 - x0 - 190,
        full.font(full.TIMES_ITALIC, 27),
        fill=full.LIGHT_MUTED,
        line_gap=4,
    )

    # A quiet glyph gives the future card structure without leaking content.
    cx, cy = (x0 + x1) // 2, y0 + 392
    draw.rounded_rectangle(
        (cx - 154, cy - 114, cx + 154, cy + 114),
        radius=18,
        fill=(255, 255, 255, 210),
        outline=(222, 230, 238, 255),
        width=2,
    )
    for dx in (-72, 0, 72):
        alpha = 70 if not dim_completed else 52
        draw.ellipse((cx + dx - 28, cy - 28, cx + dx + 28, cy + 28), outline=(*accent, alpha), width=6)
    draw.line((cx - 100, cy + 88, cx + 100, cy + 88), fill=(215, 225, 235, 230), width=3)


def draw_completed_summary(img: Image.Image, box: tuple[int, int, int, int], index: int) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    accent = PANELS[index]["accent"]
    draw.rounded_rectangle(box, radius=12, fill=(251, 253, 254, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 210))

    chip = (x0 + 40, y0 + 36, x0 + 102, y0 + 98)
    draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 225), width=3)
    draw.line((chip[0] + 17, chip[1] + 33, chip[0] + 29, chip[1] + 45), fill=(*accent, 255), width=5)
    draw.line((chip[0] + 29, chip[1] + 45, chip[0] + 47, chip[1] + 20), fill=(*accent, 255), width=5)
    draw.text((x0 + 128, y0 + 43), PANELS[index]["short"], font=full.font(full.TIMES_BOLD, 35), fill=full.INK)
    full.draw_wrapped(
        draw,
        PANELS[index]["takeaway"],
        (x0 + 128, y0 + 98),
        x1 - x0 - 172,
        full.font(full.TIMES_ITALIC, 28),
        fill=full.BLUE,
        line_gap=5,
    )
    # Keep the rest intentionally calm; this is a completed memory trace.
    mid_y = y0 + 415
    draw.rounded_rectangle((x0 + 70, mid_y - 54, x1 - 70, mid_y + 54), radius=12, fill=(240, 247, 251, 255), outline=(218, 226, 235, 255), width=2)
    draw.text((x0 + 102, mid_y - 18), "covered", font=full.font(full.TIMES_ITALIC, 29), fill=full.LIGHT_MUTED)


def draw_progress_rail(draw: ImageDraw.ImageDraw, frame: int, y: int = 302) -> None:
    labels = ["Core", "Shoulders", "Stretch / split"]
    xs = [760, 1280, 1800]
    for i in range(2):
        draw.line((xs[i] + 48, y, xs[i + 1] - 48, y), fill=(178, 194, 210, 255), width=5)
        draw.polygon([(xs[i + 1] - 48, y), (xs[i + 1] - 68, y - 12), (xs[i + 1] - 68, y + 12)], fill=(178, 194, 210, 255))
    for i, (x, label) in enumerate(zip(xs, labels)):
        active = i == frame
        done = i < frame
        accent = PANELS[i]["accent"]
        fill = accent if active else ((255, 255, 255) if done else (243, 248, 251))
        outline = accent if (active or done) else full.PANEL_EDGE
        draw.ellipse((x - 42, y - 42, x + 42, y + 42), fill=(*fill, 255), outline=(*outline, 255), width=4)
        if done:
            draw.line((x - 16, y + 2, x - 3, y + 16), fill=(*accent, 255), width=5)
            draw.line((x - 3, y + 16, x + 20, y - 18), fill=(*accent, 255), width=5)
        else:
            text = str(i + 1)
            tf = full.font(full.TIMES_BOLD, 36)
            tw, th = full.text_box(draw, text, tf)
            draw.text((x - tw / 2, y - th / 2 - 2), text, font=tf, fill=(255, 255, 255) if active else accent)
        rail_label = label if i <= frame else f"Step {i + 1}"
        lf = full.font(full.TIMES_BOLD if active else full.TIMES, 25)
        tw, _ = full.text_box(draw, rail_label, lf)
        draw.text((x - tw / 2, y + 54), rail_label, font=lf, fill=full.INK if active else full.MUTED)


def variant_a(frame: int) -> Image.Image:
    img = new_base(frame)
    boxes = panel_boxes()
    revealed = set(range(frame + 1))
    for idx, box in enumerate(boxes):
        if idx in revealed:
            PANELS[idx]["draw"](img, box)
        else:
            draw_placeholder(img, box, idx)
    bottom_takeaway(img, BOTTOM_STEP_TEXT[frame])
    full.draw_hp2026_identity_footer(img)
    return img


def variant_b(frame: int) -> Image.Image:
    img = new_base(frame)
    boxes = panel_boxes()
    for idx, box in enumerate(boxes):
        if idx < frame:
            draw_completed_summary(img, box, idx)
        elif idx == frame:
            PANELS[idx]["draw"](img, box)
        else:
            draw_placeholder(img, box, idx)
    bottom_takeaway(
        img,
        [
            "Current check: core compactness. Upcoming checks stay informationally absent until the click.",
            "Current check: seed-excluded shoulder width. The core check remains as a concise memory trace.",
            "Current check: elongation or splitting. The three checks now support the full photon-like-shower takeaway.",
        ][frame],
    )
    full.draw_hp2026_identity_footer(img)
    return img


def variant_c(frame: int) -> Image.Image:
    img = new_base(frame, subtitle=SUBTITLES[frame])
    draw = ImageDraw.Draw(img, "RGBA")
    draw_progress_rail(draw, frame, y=294)
    # One active card uses the middle/lower canvas so the audience never reads ahead.
    active_box = (348, 396, W - 348, 1138)
    PANELS[frame]["draw"](img, active_box)
    bottom_takeaway(
        img,
        [
            "Step 1 isolates the local-core question before any later variables appear.",
            "Step 2 adds the seed-excluded shoulder-width question without showing the final check yet.",
            "Step 3 completes the sequence: compact core, narrow shoulders, and no elongation or splitting.",
        ][frame],
    )
    full.draw_hp2026_identity_footer(img)
    return img


def draw_simple_opaque_cover(img: Image.Image, box: tuple[int, int, int, int], index: int) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    accent = PANELS[index]["accent"]
    draw.rounded_rectangle(
        box,
        radius=12,
        fill=(248, 251, 253, 255),
        outline=(213, 226, 235, 255),
        width=2,
    )
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 215))
    chip = (x0 + 40, y0 + 34, x0 + 96, y0 + 90)
    draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 220), width=3)
    step = str(index + 1)
    sf = full.font(full.TIMES_BOLD, 32)
    tw, th = full.text_box(draw, step, sf)
    draw.text(
        (chip[0] + (chip[2] - chip[0] - tw) / 2, chip[1] + (chip[3] - chip[1] - th) / 2 - 1),
        step,
        font=sf,
        fill=accent,
    )


def variant_d_simple_opaque(frame: int) -> Image.Image:
    # This intentionally preserves the current full slide language and layout;
    # only the inactive-panel treatment changes from blur/veil to opaque cover.
    img = full.build_slide05_shower_shape_variables(set())
    boxes = panel_boxes()
    hidden = {1, 2} if frame == 0 else ({2} if frame == 1 else set())
    for idx in sorted(hidden):
        draw_simple_opaque_cover(img, boxes[idx], idx)
    return img.convert("RGB")


def draw_core_summary_card(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 0, size=60, offset=(34, 34), bar=True)
    title = "Is the core compact?"
    tf = full.font(full.TIMES_BOLD, 36)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 18, y0 + 42), title, font=tf, fill=full.INK)

    values = [
        [0.03, 0.05, 0.07, 0.05, 0.03],
        [0.05, 0.18, 0.32, 0.18, 0.05],
        [0.07, 0.36, 0.95, 0.42, 0.07],
        [0.05, 0.20, 0.38, 0.19, 0.05],
        [0.03, 0.05, 0.07, 0.05, 0.03],
    ]
    cell = 60
    gx = x0 + (x1 - x0 - 5 * cell) // 2
    gy = y0 + 134
    full.draw_small_tower_grid(
        draw,
        (gx, gy),
        cell,
        values,
        highlight=(1, 1, 3, 3),
        highlight_color=full.PHOTON_DARK,
        secondary=(1, 1, 4, 4),
        secondary_color=full.SPHENIX_BLUE,
    )
    full.callout_label(draw, (gx + cell + 4, gy + cell + 14), "2x2 core", full.PHOTON_DARK, size=18)
    full.callout_label(draw, (gx + 3 * cell - 1, gy + 3 * cell + 8), "3x3 local core", full.SPHENIX_BLUE, size=18)
    draw.text((gx + 2 * cell + 16, gy + 2 * cell + 39), "seed", font=full.font(full.TIMES_ITALIC, 19), fill=full.INK)
    draw.line((gx + 2 * cell - 12, gy + cell + 48, gx + cell + 16, gy + cell + 18), fill=(*full.PHOTON_DARK, 180), width=2)
    draw.line((gx + 3 * cell + 19, gy + 3 * cell + 18, gx + 4 * cell - 8, gy + 4 * cell - 8), fill=(*full.SPHENIX_BLUE, 170), width=2)

    draw_core_fraction_eq(draw, (x0 + 72, y0 + 510), 30, prefix="Core fraction = ")
    draw_e1x1_e3x3_eq(draw, (x0 + 72, y0 + 584), 30, prefix="")
    draw.text((x0 + 226, y0 + 584), " = center-tower dominance", font=full.font(full.TIMES, 30), fill=full.INK)
    full.draw_wrapped(
        draw,
        "Photon-like: energy stays concentrated in the core.",
        (x0 + 72, y0 + 696),
        x1 - x0 - 122,
        full.font(full.TIMES_ITALIC, 30),
        fill=full.BLUE,
        line_gap=8,
    )


def draw_large_core_focus_panel(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 0, size=72, offset=(76, 48), bar=True)
    title = "Is the core compact?"
    tf = full.font(full.TIMES_BOLD, 43)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2, y0 + 44), title, font=tf, fill=full.INK)
    draw.text((x0 + 178, y0 + 132), "Local EMCal energy map", font=full.font(full.TIMES_BOLD, 31), fill=full.BLUE)

    values = [
        [0.03, 0.05, 0.07, 0.05, 0.03],
        [0.05, 0.18, 0.32, 0.18, 0.05],
        [0.07, 0.36, 0.95, 0.42, 0.07],
        [0.05, 0.20, 0.38, 0.19, 0.05],
        [0.03, 0.05, 0.07, 0.05, 0.03],
    ]
    cell = 80
    gx = x0 + 300
    gy = y0 + 178
    grid = full.draw_small_tower_grid(
        draw,
        (gx, gy),
        cell,
        values,
        highlight=(1, 1, 3, 3),
        highlight_color=full.PHOTON_DARK,
        secondary=(1, 1, 4, 4),
        secondary_color=full.SPHENIX_BLUE,
    )
    full.callout_label(draw, (gx + 1 * cell + 14, gy + cell + 22), "2x2 core", full.PHOTON_DARK, size=24)
    full.callout_label(draw, (gx + 3 * cell + 2, gy + 3 * cell + 20), "3x3 local core", full.SPHENIX_BLUE, size=24)
    draw.text((gx + 2 * cell + 23, gy + 2 * cell + 54), "seed", font=full.font(full.TIMES_ITALIC, 25), fill=full.INK)
    draw.line((gx + 2 * cell - 14, gy + cell + 68, gx + cell + 22, gy + cell + 24), fill=(*full.PHOTON_DARK, 190), width=3)
    draw.line((gx + 3 * cell + 28, gy + 3 * cell + 30, gx + 4 * cell - 12, gy + 4 * cell - 12), fill=(*full.SPHENIX_BLUE, 180), width=3)

    draw_core_fraction_eq(draw, (x0 + 80, y0 + 620), 38, prefix="Core fraction = ")
    draw_e1x1_e3x3_eq(draw, (x0 + 80, y0 + 704), 38)
    draw.text((x0 + 276, y0 + 704), " = center-tower dominance in the local core", font=full.font(full.TIMES, 38), fill=full.INK)
    full.draw_wrapped(
        draw,
        "Photon-like: energy stays concentrated in the core.",
        (x0 + 80, y0 + 790),
        960,
        full.font(full.TIMES_ITALIC, 37),
        fill=full.BLUE,
        line_gap=8,
    )

    note = (x0 + 1184, y0 + 176, x1 - 76, y0 + 658)
    draw.rounded_rectangle(note, radius=12, fill=(247, 251, 253, 255), outline=(218, 226, 235, 255), width=2)
    draw.rounded_rectangle((note[0], note[1], note[0] + 12, note[3]), radius=6, fill=(*full.PHOTON_DARK, 235))
    draw.text((note[0] + 42, note[1] + 32), "What this first check establishes", font=full.font(full.TIMES_BOLD, 36), fill=full.INK)
    rows = [
        ("local core", "energy concentrated near the cluster seed"),
        ("core fraction", "2x2 core energy divided by total cluster energy"),
        ("E1x1/E3x3", "center-tower dominance inside the local 3x3"),
    ]
    y = note[1] + 104
    for label, body in rows:
        if label == "E1x1/E3x3":
            draw_e1x1_e3x3_eq(draw, (note[0] + 44, y), 30, fill=full.PHOTON_DARK)
        else:
            draw.text((note[0] + 44, y), label, font=full.font(full.TIMES_BOLD, 30), fill=full.PHOTON_DARK)
        full.draw_wrapped(draw, body, (note[0] + 246, y), note[2] - note[0] - 286, full.font(full.TIMES, 30), fill=full.MUTED, line_gap=4)
        y += 92
    takeaway = (note[0], note[3] + 26, note[2], note[3] + 146)
    draw.rounded_rectangle(takeaway, radius=10, fill=(255, 250, 239, 255), outline=(234, 206, 143, 255), width=2)
    full.draw_wrapped(
        draw,
        "A prompt-photon-like cluster should begin as one concentrated EMCal deposit.",
        (takeaway[0] + 34, takeaway[1] + 27),
        takeaway[2] - takeaway[0] - 68,
        full.font(full.TIMES_BOLD, 30),
        fill=full.BLUE,
        line_gap=5,
    )


def draw_large_shoulders_panel(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 1, size=66, offset=(40, 36), bar=True)
    title = "Are the shoulders narrow?"
    tf = full.font(full.TIMES_BOLD, 43)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2, y0 + 42), title, font=tf, fill=full.INK)

    compact = [
        [0.01, 0.02, 0.04, 0.02, 0.01],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.04, 0.18, 0.92, 0.18, 0.04],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.01, 0.02, 0.04, 0.02, 0.01],
    ]
    broad = [
        [0.10, 0.18, 0.25, 0.18, 0.10],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.25, 0.48, 0.80, 0.48, 0.25],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.10, 0.18, 0.25, 0.18, 0.10],
    ]
    cell = 76
    total_grid_w = 2 * 5 * cell + 160
    left_x = x0 + (x1 - x0 - total_grid_w) // 2
    gy = y0 + 152
    grids = [
        full.draw_small_tower_grid(draw, (left_x, gy), cell, compact),
        full.draw_small_tower_grid(draw, (left_x + 5 * cell + 160, gy), cell, broad),
    ]
    for idx, grid in enumerate(grids):
        gx0, gy0, gx1, _ = grid
        if idx == 0:
            start, end = gx0 + 138, gx1 - 138
        else:
            start, end = gx0 + 34, gx1 - 34
        draw.line((start, gy0 - 28, end, gy0 - 28), fill=(230, 70, 45, 230), width=5)
        draw.line((start, gy0 - 40, start, gy0 - 16), fill=(230, 70, 45, 230), width=5)
        draw.line((end, gy0 - 40, end, gy0 - 16), fill=(230, 70, 45, 230), width=5)
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 215), outline=(160, 168, 176, 215), width=2)
        draw.line((center_cell[0] + 10, center_cell[1] + 10, center_cell[2] - 10, center_cell[3] - 10), fill=(160, 168, 176, 190), width=3)
        draw.line((center_cell[0] + 10, center_cell[3] - 10, center_cell[2] - 10, center_cell[1] + 10), fill=(160, 168, 176, 190), width=3)
        draw.ellipse((gx0 + 2.5 * cell - 10, gy0 + 2.5 * cell - 10, gx0 + 2.5 * cell + 10, gy0 + 2.5 * cell + 10), fill=(0, 0, 0, 255))
    for label, grid in (("photon-like", grids[0]), ("broad / multi-tower", grids[1])):
        gx0, _, gx1, _ = grid
        lf = full.font(full.TIMES_ITALIC, 31)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y0 + 562), label, font=lf, fill=full.MUTED)
    full.callout_label(draw, (x0 + (x1 - x0) // 2 - 78, y0 + 495), "seed removed", full.LIGHT_MUTED, size=22)
    full.draw_formula_run(
        draw,
        (x0 + 78, y0 + 636),
        [("Shoulder width", 36, 0, full.TIMES_BOLD), (" = seed-excluded spread", 36, 0, full.TIMES)],
        fill=full.INK,
    )
    full.draw_formula_run(
        draw,
        (x0 + 78, y0 + 684),
        [("w", 28, 0, full.TIMES_BOLD), ("η", 20, 11, full.TIMES_BOLD), ("cogX", 17, -11, full.TIMES_BOLD), (" and w", 28, 0, full.TIMES_BOLD), ("φ", 20, 11, full.TIMES_BOLD), ("cogX", 17, -11, full.TIMES_BOLD)],
        fill=full.MUTED,
    )
    full.draw_wrapped(
        draw,
        "Tests the surrounding shower after the hottest tower is removed.",
        (x0 + 78, y0 + 732),
        x1 - x0 - 156,
        full.font(full.TIMES, 32),
        fill=full.MUTED,
        line_gap=8,
    )
    full.draw_wrapped(
        draw,
        "Photon-like: narrow shoulders around the core.",
        (x0 + 78, y0 + 790),
        x1 - x0 - 156,
        full.font(full.TIMES_ITALIC, 36),
        fill=full.BLUE,
        line_gap=8,
    )


def draw_shoulders_summary_card(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 1, size=58, offset=(30, 30), bar=True)
    title = "Are the shoulders narrow?"
    tf = full.font(full.TIMES_BOLD, 36)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 18, y0 + 40), title, font=tf, fill=full.INK)

    compact = [
        [0.01, 0.02, 0.04, 0.02, 0.01],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.04, 0.18, 0.92, 0.18, 0.04],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.01, 0.02, 0.04, 0.02, 0.01],
    ]
    broad = [
        [0.10, 0.18, 0.25, 0.18, 0.10],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.25, 0.48, 0.80, 0.48, 0.25],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.10, 0.18, 0.25, 0.18, 0.10],
    ]
    cell = 58
    gap = 46
    total = 2 * 5 * cell + gap
    left_x = x0 + (x1 - x0 - total) // 2
    gy = y0 + 132
    grids = [
        full.draw_small_tower_grid(draw, (left_x, gy), cell, compact),
        full.draw_small_tower_grid(draw, (left_x + 5 * cell + gap, gy), cell, broad),
    ]
    for idx, grid in enumerate(grids):
        gx0, gy0, gx1, _ = grid
        start, end = (gx0 + 106, gx1 - 106) if idx == 0 else (gx0 + 26, gx1 - 26)
        draw.line((start, gy0 - 20, end, gy0 - 20), fill=(230, 70, 45, 220), width=4)
        draw.line((start, gy0 - 30, start, gy0 - 10), fill=(230, 70, 45, 220), width=4)
        draw.line((end, gy0 - 30, end, gy0 - 10), fill=(230, 70, 45, 220), width=4)
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 215), outline=(160, 168, 176, 215), width=2)
        draw.line((center_cell[0] + 8, center_cell[1] + 8, center_cell[2] - 8, center_cell[3] - 8), fill=(160, 168, 176, 190), width=3)
        draw.line((center_cell[0] + 8, center_cell[3] - 8, center_cell[2] - 8, center_cell[1] + 8), fill=(160, 168, 176, 190), width=3)
        draw.ellipse((gx0 + 2.5 * cell - 8, gy0 + 2.5 * cell - 8, gx0 + 2.5 * cell + 8, gy0 + 2.5 * cell + 8), fill=(0, 0, 0, 255))
    for label, grid in (("photon-like", grids[0]), ("broad / multi-tower", grids[1])):
        gx0, _, gx1, _ = grid
        lf = full.font(full.TIMES_ITALIC, 27)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y0 + 438), label, font=lf, fill=full.MUTED)
    full.draw_formula_run(
        draw,
        (x0 + 54, y0 + 530),
        [("Shoulder width", 30, 0, full.TIMES_BOLD), (" = seed-excluded spread", 30, 0, full.TIMES)],
        fill=full.INK,
    )
    full.draw_wrapped(
        draw,
        "Photon-like: narrow shoulders around the core.",
        (x0 + 54, y0 + 682),
        x1 - x0 - 96,
        full.font(full.TIMES_ITALIC, 30),
        fill=full.BLUE,
        line_gap=7,
    )


def draw_stretch_summary_card(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 2, size=58, offset=(30, 30), bar=True)
    title = "Is the shower stretched or split?"
    tf = full.font(full.TIMES_BOLD, 34)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 22, y0 + 40), title, font=tf, fill=full.INK)

    compact = [
        [0.02, 0.04, 0.06, 0.04, 0.02],
        [0.03, 0.10, 0.20, 0.10, 0.03],
        [0.05, 0.32, 0.95, 0.34, 0.05],
        [0.03, 0.10, 0.20, 0.10, 0.03],
        [0.02, 0.04, 0.06, 0.04, 0.02],
    ]
    split = [
        [0.02, 0.05, 0.12, 0.08, 0.02],
        [0.03, 0.12, 0.46, 0.28, 0.05],
        [0.05, 0.22, 0.95, 0.36, 0.08],
        [0.03, 0.42, 0.52, 0.14, 0.04],
        [0.02, 0.30, 0.26, 0.08, 0.02],
    ]
    cell = 56
    gap = 46
    total = 2 * 5 * cell + gap
    left_x = x0 + (x1 - x0 - total) // 2 + 4
    gy = y0 + 138
    grids = [
        full.draw_small_tower_grid(draw, (left_x, gy), cell, compact),
        full.draw_small_tower_grid(draw, (left_x + 5 * cell + gap, gy), cell, split),
    ]
    for grid in grids:
        gx0, gy0, _, _ = grid
        draw.rectangle((gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell), outline=(*full.SPHENIX_BLUE, 210), width=3)
        draw.rectangle((gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 3 * cell), outline=(*full.PHOTON_DARK, 220), width=3)
        draw.ellipse((gx0 + 2.5 * cell - 8, gy0 + 2.5 * cell - 8, gx0 + 2.5 * cell + 8, gy0 + 2.5 * cell + 8), fill=(0, 0, 0, 255))
    for label, grid in (("compact local", grids[0]), ("stretched / split", grids[1])):
        gx0, _, gx1, _ = grid
        lf = full.font(full.TIMES_ITALIC, 27)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y0 + 430), label, font=lf, fill=full.MUTED)
    draw_e3x2_e3x5_eq(draw, (x0 + 54, y0 + 528), 30)
    full.draw_wrapped(
        draw,
        "Photon-like: compact energy stays confined to the narrow strip.",
        (x0 + 54, y0 + 678),
        x1 - x0 - 96,
        full.font(full.TIMES_ITALIC, 30),
        fill=full.BLUE,
        line_gap=7,
    )


def draw_clean_full_synthesis() -> Image.Image:
    img = new_base(
        2,
        subtitle="A photon-like EMCal cluster has a concentrated core, narrow shoulders, and little evidence of elongation or splitting.",
    )
    boxes = panel_boxes()
    draw_core_summary_card(img, boxes[0])
    draw_shoulders_summary_card(img, boxes[1])
    draw_stretch_summary_card(img, boxes[2])
    bottom_takeaway(
        img,
        "Together, these variables ask whether the EMCal energy looks like one compact photon, rather than a broad, elongated, or split decay-like shower.",
    )
    full.draw_hp2026_identity_footer(img)
    return img.convert("RGB")


def variant_e_expand_collapse(frame: int) -> Image.Image:
    if frame == 2:
        return draw_clean_full_synthesis()

    subtitle = [
        "Start with one readable shower-shape question: is the energy concentrated in the local core?",
        "Keep the core check as context, then add the seed-excluded shoulder-width question.",
    ][frame]
    img = new_base(frame, subtitle=subtitle)
    if frame == 0:
        draw_large_core_focus_panel(img, (132, 324, W - 132, 1166))
    else:
        draw_core_summary_card(img, (132, 324, 885, 1166))
        draw_large_shoulders_panel(img, (935, 324, W - 132, 1166))
    full.draw_hp2026_identity_footer(img)
    return img.convert("RGB")


VARIANTS = [
    ("A", "Opaque placeholders", "Minimal replacement: same three-card geometry; future cards are present but informationally absent.", variant_a),
    ("B", "Current / completed / upcoming", "Most guided delivery: current card stays detailed; covered cards collapse to one-line memory traces.", variant_b),
    ("C", "Single active card + progress rail", "Cleanest one-idea-at-a-time reveal; strongest redesign and least future-reading risk.", variant_c),
]


def save_variant_frames() -> dict[str, list[Path]]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths: dict[str, list[Path]] = {}
    for key, _, _, builder in VARIANTS:
        paths[key] = []
        for frame in range(3):
            img = builder(frame).convert("RGB")
            path = OUT_DIR / f"variant_{key}_frame{frame + 1}_slide{8 + frame:02d}.png"
            img.save(path, "PNG")
            paths[key].append(path)
    return paths


def make_contact_sheet(paths: dict[str, list[Path]]) -> Path:
    thumb_w, thumb_h = 768, 432
    left_label_w = 410
    top_label_h = 88
    row_h = thumb_h + 88
    sheet = Image.new("RGB", (left_label_w + 3 * thumb_w, top_label_h + len(VARIANTS) * row_h), (246, 249, 252))
    draw = ImageDraw.Draw(sheet, "RGBA")
    title_font = full.font(full.TIMES_BOLD, 36)
    sub_font = full.font(full.TIMES, 24)
    small_font = full.font(full.TIMES_ITALIC, 21)
    col_font = full.font(full.TIMES_BOLD, 28)

    draw.text((26, 22), "Slides 8-10 no-blur strategy comparison", font=title_font, fill=full.INK)
    for col, label in enumerate(("Slide 8 / first click", "Slide 9 / second click", "Slide 10 / full reveal")):
        x = left_label_w + col * thumb_w
        tw, _ = full.text_box(draw, label, col_font)
        draw.text((x + (thumb_w - tw) / 2, 30), label, font=col_font, fill=full.BLUE)

    for row, (key, title, desc, _) in enumerate(VARIANTS):
        y = top_label_h + row * row_h
        draw.rounded_rectangle((18, y + 20, left_label_w - 22, y + row_h - 20), radius=14, fill=(255, 255, 255, 255), outline=(218, 226, 235, 255), width=2)
        draw.text((42, y + 52), f"Variant {key}", font=full.font(full.TIMES_BOLD, 34), fill=full.INK)
        draw.text((42, y + 104), title, font=sub_font, fill=full.BLUE)
        full.draw_wrapped(draw, desc, (42, y + 148), left_label_w - 84, small_font, fill=full.MUTED, line_gap=5)
        verdict = {
            "A": "Recommended default",
            "B": "Best guided narration",
            "C": "Cleanest but bigger redesign",
        }[key]
        draw.rounded_rectangle((42, y + row_h - 82, left_label_w - 58, y + row_h - 40), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
        draw.text((62, y + row_h - 75), verdict, font=full.font(full.TIMES_BOLD, 22), fill=full.TEAL if key == "A" else full.BLUE)

        for col, path in enumerate(paths[key]):
            thumb = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
            x = left_label_w + col * thumb_w
            sheet.paste(thumb, (x, y + 48))
            draw.rectangle((x, y + 48, x + thumb_w - 1, y + 48 + thumb_h - 1), outline=(204, 216, 228, 255), width=2)

    contact = OUT_DIR / "slides8_10_no_blur_strategy_comparison.png"
    sheet.save(contact, "PNG")
    return contact


def save_simple_opaque_option() -> tuple[list[Path], Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths = []
    for frame in range(3):
        img = variant_d_simple_opaque(frame)
        path = OUT_DIR / f"variant_D_simple_opaque_frame{frame + 1}_slide{8 + frame:02d}.png"
        img.save(path, "PNG")
        paths.append(path)

    thumb_w, thumb_h = 768, 432
    label_h = 82
    left_w = 460
    sheet = Image.new("RGB", (left_w + 3 * thumb_w, label_h + thumb_h + 40), (246, 249, 252))
    draw = ImageDraw.Draw(sheet, "RGBA")
    draw.rounded_rectangle((18, 22, left_w - 24, label_h + thumb_h + 18), radius=14, fill=(255, 255, 255, 255), outline=(218, 226, 235, 255), width=2)
    draw.text((42, 54), "Variant D", font=full.font(full.TIMES_BOLD, 38), fill=full.INK)
    draw.text((42, 112), "Simple opaque cover", font=full.font(full.TIMES, 28), fill=full.BLUE)
    full.draw_wrapped(
        draw,
        "Closest to the current blurred build: keep the same slide language and geometry; cover unrevealed cards with an opaque clean panel.",
        (42, 164),
        left_w - 88,
        full.font(full.TIMES_ITALIC, 24),
        fill=full.MUTED,
        line_gap=5,
    )
    draw.rounded_rectangle((42, label_h + thumb_h - 54, left_w - 58, label_h + thumb_h - 12), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
    draw.text((62, label_h + thumb_h - 47), "Lowest-risk replacement", font=full.font(full.TIMES_BOLD, 22), fill=full.TEAL)
    for col, (path, label) in enumerate(zip(paths, ("Slide 8 / first click", "Slide 9 / second click", "Slide 10 / full reveal"))):
        x = left_w + col * thumb_w
        lf = full.font(full.TIMES_BOLD, 28)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((x + (thumb_w - tw) / 2, 30), label, font=lf, fill=full.BLUE)
        thumb = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(thumb, (x, label_h))
        draw.rectangle((x, label_h, x + thumb_w - 1, label_h + thumb_h - 1), outline=(204, 216, 228, 255), width=2)
    contact = OUT_DIR / "variant_D_simple_opaque_cover_contact_sheet.png"
    sheet.save(contact, "PNG")
    return paths, contact


def save_expand_collapse_option() -> tuple[list[Path], Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths = []
    for frame in range(3):
        img = variant_e_expand_collapse(frame)
        path = OUT_DIR / f"variant_E_expand_collapse_frame{frame + 1}_slide{8 + frame:02d}.png"
        img.save(path, "PNG")
        paths.append(path)

    thumb_w, thumb_h = 768, 432
    label_h = 82
    left_w = 470
    sheet = Image.new("RGB", (left_w + 3 * thumb_w, label_h + thumb_h + 40), (246, 249, 252))
    draw = ImageDraw.Draw(sheet, "RGBA")
    draw.rounded_rectangle((18, 22, left_w - 24, label_h + thumb_h + 18), radius=14, fill=(255, 255, 255, 255), outline=(218, 226, 235, 255), width=2)
    draw.text((42, 54), "Variant E", font=full.font(full.TIMES_BOLD, 38), fill=full.INK)
    draw.text((42, 112), "Expand → collapse → synthesize", font=full.font(full.TIMES, 27), fill=full.BLUE)
    full.draw_wrapped(
        draw,
        "Uses space for one large readable active idea, then compresses covered content into context before the final three-card synthesis.",
        (42, 164),
        left_w - 88,
        full.font(full.TIMES_ITALIC, 24),
        fill=full.MUTED,
        line_gap=5,
    )
    draw.rounded_rectangle((42, label_h + thumb_h - 54, left_w - 58, label_h + thumb_h - 12), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
    draw.text((62, label_h + thumb_h - 47), "Best readability-first option", font=full.font(full.TIMES_BOLD, 22), fill=full.TEAL)
    for col, (path, label) in enumerate(zip(paths, ("Slide 8 / large core", "Slide 9 / core summary + shoulders", "Slide 10 / full synthesis"))):
        x = left_w + col * thumb_w
        lf = full.font(full.TIMES_BOLD, 28)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((x + (thumb_w - tw) / 2, 30), label, font=lf, fill=full.BLUE)
        thumb = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(thumb, (x, label_h))
        draw.rectangle((x, label_h, x + thumb_w - 1, label_h + thumb_h - 1), outline=(204, 216, 228, 255), width=2)
    contact = OUT_DIR / "variant_E_expand_collapse_contact_sheet.png"
    sheet.save(contact, "PNG")
    return paths, contact


def main() -> None:
    paths = save_variant_frames()
    contact = make_contact_sheet(paths)
    simple_paths, simple_contact = save_simple_opaque_option()
    expand_paths, expand_contact = save_expand_collapse_option()
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "source_generator": str((ROOT / "scripts/slides/hp2026/fulltalk/make_hp2026_fulltalk_candidates.py").relative_to(ROOT)),
        "source_functions": [
            "build_slide05_shower_shape_variables()",
            "slide05_focus_core_compactness()",
            "slide05_focus_shoulders()",
            "slide04()",
            "soften_focus_region() is the current blur/fade helper being replaced in these prototypes",
        ],
        "design_goal": "Replace blurred/faded unrevealed shower-shape cards with staged disclosure treatments where future content is structurally present but unreadable.",
        "variants": [
            {
                "key": key,
                "title": title,
                "description": desc,
                "outputs": [str(p.relative_to(ROOT)) for p in paths[key]],
            }
            for key, title, desc, _ in VARIANTS
        ],
        "contact_sheet": str(contact.relative_to(ROOT)),
        "simple_opaque_cover_contact_sheet": str(simple_contact.relative_to(ROOT)),
        "simple_opaque_cover_outputs": [str(p.relative_to(ROOT)) for p in simple_paths],
        "expand_collapse_contact_sheet": str(expand_contact.relative_to(ROOT)),
        "expand_collapse_outputs": [str(p.relative_to(ROOT)) for p in expand_paths],
        "notes": [
            "No Google Slides deck was changed.",
            "The primary full-talk generator was imported, not edited.",
            "These frames are local strategy prototypes for Justin to choose a no-blur pattern before propagating across HPslides_v1.",
        ],
    }
    manifest_path = OUT_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(contact)
    print(simple_contact)
    print(expand_contact)
    print(manifest_path)
    for key in sorted(paths):
        for path in paths[key]:
            print(path)


if __name__ == "__main__":
    main()
