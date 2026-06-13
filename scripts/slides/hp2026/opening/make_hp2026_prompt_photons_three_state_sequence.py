#!/usr/bin/env python3
"""Render a local prompt-photon opening sequence.

This is a source-native PNG prototype only. It does not mutate Google Slides.
The current sequence is designed as two conference-readable frames after
removing the heavy-ion/R_AA color-neutral detour:

1. p+p production channels and decay backgrounds
2. isolation definition and isolated/non-isolated intuition
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from functools import lru_cache
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_opening_motivation_slide as hp


OUT_DIR = (
    hp.ROOT
    / "outputs"
    / "manual-20260611-slide5_6_prompt_photons_pp_isolation"
)

CURRENT_ACCEPTED_BASELINE = (
    hp.ROOT
    / "outputs"
    / "manual-20260607-slide5_prompt_photons_refine"
    / "build_states_v71_iso_relationship_callout"
    / "slide5_state01_full_active_rhs_swapped.png"
)

FOOTER_TOP = 1326
BODY_TOP = 272
BODY_BOTTOM = 1290
TITLE_DIVIDER_Y = 232


def spring_points(
    start: tuple[float, float],
    end: tuple[float, float],
    amp: float,
    loops: float,
    steps: int = 140,
) -> list[tuple[float, float]]:
    sx, sy = start
    ex, ey = end
    vx, vy = ex - sx, ey - sy
    length = math.hypot(vx, vy)
    if length == 0:
        return [start]
    ux, uy = vx / length, vy / length
    px, py = -uy, ux
    pts: list[tuple[float, float]] = []
    for i in range(steps):
        t = i / (steps - 1)
        phase = math.tau * loops * t
        along = t * length + 0.9 * amp * math.sin(phase)
        off = amp * math.cos(phase)
        pts.append((sx + ux * along + px * off, sy + uy * along + py * off))
    return pts


def arrow_line(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    pos: float = 0.68,
    width: int = 3,
    head: int = 10,
) -> None:
    sx, sy = start
    ex, ey = end
    draw.line((sx, sy, ex, ey), fill=(*hp.INK, 222), width=width)
    hx, hy = sx + (ex - sx) * pos, sy + (ey - sy) * pos
    angle = math.atan2(ey - sy, ex - sx)
    left = (hx - head * math.cos(angle - 0.55), hy - head * math.sin(angle - 0.55))
    right = (hx - head * math.cos(angle + 0.55), hy - head * math.sin(angle + 0.55))
    draw.polygon([(hx, hy), left, right], fill=(*hp.INK, 222))


def draw_feynman(
    draw: ImageDraw.ImageDraw,
    cx: int,
    cy: int,
    scale: float,
    mode: str,
) -> None:
    def photon(start: tuple[float, float], end: tuple[float, float]) -> None:
        pts = hp.feynman_points(start, end, 5.3 * scale, 4.2, 86)
        hp.draw_polyline(draw, pts, (*hp.PHOTON_DARK, 235), max(3, round(4 * scale)))

    def gluon(start: tuple[float, float], end: tuple[float, float], loops: float) -> None:
        pts = spring_points(start, end, 6.0 * scale, loops)
        hp.draw_polyline(draw, pts, (*hp.TEAL, 235), max(3, round(4 * scale)))

    node_top = (cx, cy - 42 * scale)
    node_bot = (cx, cy + 64 * scale)
    if mode == "fragmentation":
        node_top = (cx, cy - 4 * scale)
        node_bot = (cx, cy + 86 * scale)
        arrow_line(draw, (cx - 108 * scale, cy - 88 * scale), node_top, pos=0.66, width=max(3, round(4 * scale)))
        arrow_line(draw, node_top, (cx + 112 * scale, cy - 82 * scale), pos=0.74, width=max(3, round(4 * scale)))
        photon((cx + 46 * scale, cy - 38 * scale), (cx + 116 * scale, cy + 8 * scale))
        gluon(node_top, node_bot, 6.6)
        arrow_line(draw, (cx - 108 * scale, cy + 150 * scale), node_bot, pos=0.58, width=max(3, round(4 * scale)))
        arrow_line(draw, node_bot, (cx + 112 * scale, cy + 150 * scale), pos=0.74, width=max(3, round(4 * scale)))
        labels = [
            ("q", cx - 132 * scale, cy - 110 * scale),
            ("q", cx + 124 * scale, cy - 110 * scale),
            ("γ", cx + 122 * scale, cy + 8 * scale),
            ("q", cx - 132 * scale, cy + 116 * scale),
            ("q", cx + 124 * scale, cy + 116 * scale),
        ]
    else:
        arrow_line(draw, (cx - 108 * scale, cy - 112 * scale), node_top, pos=0.66, width=max(3, round(4 * scale)))
        photon(node_top, (cx + 112 * scale, cy - 112 * scale))
        arrow_line(draw, node_top, node_bot, pos=0.72, width=max(3, round(4 * scale)))
        if mode == "compton":
            gluon((cx - 112 * scale, cy + 128 * scale), node_bot, 7.4)
            arrow_line(draw, node_bot, (cx + 112 * scale, cy + 128 * scale), pos=0.78, width=max(3, round(4 * scale)))
            labels = [
                ("q", cx - 132 * scale, cy - 134 * scale),
                ("γ", cx + 116 * scale, cy - 132 * scale),
                ("g", cx - 128 * scale, cy + 112 * scale),
                ("q", cx + 124 * scale, cy + 112 * scale),
            ]
        else:
            arrow_line(draw, (cx - 112 * scale, cy + 128 * scale), node_bot, pos=0.55, width=max(3, round(4 * scale)))
            gluon(node_bot, (cx + 112 * scale, cy + 128 * scale), 7.3)
            labels = [
                ("q", cx - 132 * scale, cy - 134 * scale),
                ("γ", cx + 116 * scale, cy - 132 * scale),
                ("qbar", cx - 132 * scale, cy + 112 * scale),
                ("g", cx + 116 * scale, cy + 112 * scale),
            ]

    for node in (node_top, node_bot):
        draw.ellipse(
            (
                node[0] - 7 * scale,
                node[1] - 7 * scale,
                node[0] + 7 * scale,
                node[1] + 7 * scale,
            ),
            fill=hp.INK,
        )
    label_font = hp.font(hp.TIMES_ITALIC, max(35, round(44 * scale)))
    for txt, lx, ly in labels:
        color = hp.PHOTON_DARK if txt == "γ" else hp.TEAL if txt == "g" else hp.INK
        label_text = "q" if txt == "qbar" else txt
        for dx, dy in ((-2, 0), (2, 0), (0, -2), (0, 2)):
            draw.text((lx + dx, ly + dy), label_text, font=label_font, fill=(255, 255, 255, 235))
        draw.text((lx, ly), label_text, font=label_font, fill=color)
        if txt == "qbar":
            tw, _ = hp.text_box(draw, "q", label_font)
            draw.line((lx, ly + 4 * scale, lx + tw, ly + 4 * scale), fill=(*color, 240), width=max(2, round(2 * scale)))


def box_to_list(box: tuple[float, float, float, float]) -> list[float]:
    return [round(v, 3) for v in box]


def text_bbox(
    draw: ImageDraw.ImageDraw,
    text: str,
    font_obj: object,
    x: float,
    y: float,
    *,
    rich: bool = False,
) -> tuple[float, float, float, float]:
    if rich:
        w, h = hp.rich_text_box(draw, text, font_obj)
    else:
        w, h = hp.text_box(draw, text, font_obj)
    return (x, y, x + w, y + h)


def add_text_node(
    nodes: list[dict[str, object]],
    *,
    name: str,
    text: str,
    font_px: int,
    bbox: tuple[float, float, float, float],
    role: str = "audience",
) -> None:
    nodes.append(
        {
            "name": name,
            "kind": "text",
            "role": role,
            "text": text,
            "font_px": font_px,
            "bbox": box_to_list(bbox),
        }
    )


@lru_cache(maxsize=None)
def feynman_ink_bbox(scale_key: int, mode: str) -> tuple[int, int, int, int]:
    """Return the alpha ink bounds for a labeled Feynman diagram.

    `scale_key` is `round(scale * 1000)` so the cached measurement survives
    float noise while still tracking the visible scale.
    """

    scale = scale_key / 1000
    tile_size = 620
    center = tile_size // 2
    tile = Image.new("RGBA", (tile_size, tile_size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(tile, "RGBA")
    draw_feynman(draw, center, center, scale, mode)
    bbox = tile.getchannel("A").getbbox()
    if bbox is None:
        return (center, center, center, center)
    return bbox


def feynman_visual_offset(scale: float, mode: str) -> tuple[float, float]:
    """Offset the nominal vertex center so the full visible ink is centered."""

    tile_size = 620
    center = tile_size / 2
    x0, y0, x1, y1 = feynman_ink_bbox(round(scale * 1000), mode)
    return (center - (x0 + x1) / 2, center - (y0 + y1) / 2)


def feynman_visual_box(cx: float, cy: float, scale: float, mode: str) -> tuple[float, float, float, float]:
    """Visible ink box on the slide after visual-envelope centering."""

    tile_size = 620
    center = tile_size / 2
    x0, y0, x1, y1 = feynman_ink_bbox(round(scale * 1000), mode)
    dx, dy = feynman_visual_offset(scale, mode)
    return (cx + dx + x0 - center, cy + dy + y0 - center, cx + dx + x1 - center, cy + dy + y1 - center)


def draw_feynman_centered(
    draw: ImageDraw.ImageDraw,
    cx: float,
    cy: float,
    scale: float,
    mode: str,
) -> None:
    dx, dy = feynman_visual_offset(scale, mode)
    draw_feynman(draw, round(cx + dx), round(cy + dy), scale, mode)


def draw_chrome(
    layout_nodes: list[dict[str, object]] | None = None,
    *,
    frame_id: str = "",
    title_text: str,
) -> tuple[Image.Image, ImageDraw.ImageDraw]:
    img = Image.new("RGBA", (hp.W, hp.H), (*hp.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, hp.W, hp.H), fill=(*hp.SOFT_BG, 255))
    draw.rectangle((0, 0, hp.W, 22), fill=(*hp.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, hp.W, 30), fill=(*hp.PHOTON, 255))

    logo_path = hp.TITLE_ASSET_DIR / "sphenix-logo-white-bg_0.png"
    if logo_path.exists():
        logo = hp.fit(hp.crop_visible(Image.open(logo_path).convert("RGBA"), white_threshold=252), 366, 107)
        img.alpha_composite(logo, (2432 - logo.width, 58))

    title_font_px = 86
    title_font = hp.font(hp.TIMES_BOLD, title_font_px)
    title_xy = (132, 76)
    draw.text(title_xy, title_text, font=title_font, fill=hp.INK)
    if layout_nodes is not None:
        add_text_node(
            layout_nodes,
            name=f"{frame_id} slide title",
            text=title_text,
            font_px=title_font_px,
            bbox=text_bbox(draw, title_text, title_font, *title_xy),
        )
    draw.line((132, TITLE_DIVIDER_Y, hp.W - 132, TITLE_DIVIDER_Y), fill=(221, 226, 232, 255), width=3)
    return img, draw


def apply_common_footer(img: Image.Image) -> None:
    """Paste one identical footer strip into every sequence frame."""
    if CURRENT_ACCEPTED_BASELINE.exists():
        footer_source = Image.open(CURRENT_ACCEPTED_BASELINE).convert("RGBA")
        if footer_source.size != (hp.W, hp.H):
            footer_source = footer_source.resize((hp.W, hp.H), Image.Resampling.LANCZOS)
        footer = footer_source.crop((0, FOOTER_TOP, hp.W, hp.H))
        img.alpha_composite(footer, (0, FOOTER_TOP))
        return
    hp.draw_recreated_footer(img)


def draw_card(
    draw: ImageDraw.ImageDraw,
    card: tuple[int, int, int, int],
    stripe: tuple[int, int, int],
    *,
    stripe_w: int = 14,
) -> None:
    draw.rounded_rectangle(card, radius=8, fill=(255, 255, 255, 248), outline=(222, 229, 236, 255), width=2)
    # Draw the accent strip explicitly rather than with a narrow rounded
    # rectangle. PIL's narrow rounded caps can leave a one-pixel white notch at
    # the top/bottom apex after scaling into Slides.
    x0, y0, _, y1 = card
    r = max(1, stripe_w // 2)
    color = (*stripe, 255)
    draw.ellipse((x0, y0, x0 + stripe_w, y0 + 2 * r), fill=color)
    draw.rectangle((x0, y0 + r, x0 + stripe_w, y1 - r), fill=color)
    draw.ellipse((x0, y1 - 2 * r, x0 + stripe_w, y1), fill=color)


def centered_label(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    text: str,
    color: tuple[int, int, int],
    fill: tuple[int, int, int],
    y_offset: int,
    size: int,
) -> None:
    x0, y0, x1, _ = box
    fnt = hp.font(hp.TIMES_BOLD, size)
    tw, th = hp.text_box(draw, text, fnt)
    lx0 = round((x0 + x1 - tw) / 2) - 22
    label_box = (lx0, y0 + y_offset, lx0 + 44 + tw, y0 + y_offset + 56)
    draw.rounded_rectangle(label_box, radius=20, fill=(*fill, 255), outline=(*color, 210), width=2)
    draw.text((label_box[0] + 18, label_box[1] + (56 - th) / 2 - 3), text, font=fnt, fill=color)


def draw_centered_segments(
    draw: ImageDraw.ImageDraw,
    segments: list[tuple[str, object, tuple[int, int, int], int]],
    y: int,
    x0: int,
    x1: int,
) -> tuple[float, float, float, float]:
    total = sum(hp.rich_text_box(draw, text, fnt)[0] + gap for text, fnt, _, gap in segments)
    x = (x0 + x1 - total) / 2
    start_x = x
    max_h = 0
    for text, fnt, fill, gap in segments:
        hp.draw_rich_text(draw, (round(x), y), text, fnt, fill)
        tw, th = hp.rich_text_box(draw, text, fnt)
        max_h = max(max_h, th)
        x += tw + gap
    return (start_x, y, x, y + max_h)


def draw_cone_heading(
    draw: ImageDraw.ImageDraw,
    *,
    text: str,
    box: tuple[int, int, int, int],
    y: int,
    font_px: int,
    fill: tuple[int, int, int],
    layout_nodes: list[dict[str, object]] | None = None,
    frame_id: str = "",
) -> tuple[float, float, float, float]:
    """Draw a large label above a cone sketch, centered on the sketch box."""

    fnt = hp.font(hp.TIMES_BOLD, font_px)
    tw, th = hp.text_box(draw, text, fnt)
    x = round((box[0] + box[2] - tw) / 2)
    draw.text((x, y), text, font=fnt, fill=fill)
    bbox = (x, y, x + tw, y + th)
    if layout_nodes is not None:
        add_text_node(
            layout_nodes,
            name=f"{frame_id} {text} cone label",
            text=text,
            font_px=font_px,
            bbox=bbox,
        )
    return bbox


def draw_production_card(
    img: Image.Image,
    card: tuple[int, int, int, int],
    *,
    expanded: bool,
    focused: bool = False,
    show_subtitle: bool = True,
    layout_nodes: list[dict[str, object]] | None = None,
    frame_id: str = "",
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    draw_card(draw, card, (197, 64, 48), stripe_w=14)

    if expanded:
        title_font_px = 76
        title_font = hp.font(hp.TIMES_BOLD, title_font_px)
        title_text = "Production channels"
        title_xy = (card[0] + 56, card[1] + 54)
        draw.text(title_xy, title_text, font=title_font, fill=hp.INK)
        if layout_nodes is not None:
            add_text_node(layout_nodes, name=f"{frame_id} production card title", text=title_text, font_px=title_font_px, bbox=text_bbox(draw, title_text, title_font, *title_xy))
        # Invisible-guide layout:
        #   - center the red production figure on the usable visual band
        #   - center blue/green sub-boxes on the red-box midline
        #   - center each Feynman diagram on equal visual cells
        # This keeps the frame artistically symmetric without drawing guides.
        visual_band_top = card[1] + 220
        visual_band_bottom = card[1] + 684
        visual_mid_x = (card[0] + card[2]) / 2
        visual_mid_y = (visual_band_top + visual_band_bottom) / 2
        prompt_w = 1650
        prompt_h = 420
        prompt_box = (
            round(visual_mid_x - prompt_w / 2),
            round(visual_mid_y - prompt_h / 2),
            round(visual_mid_x + prompt_w / 2),
            round(visual_mid_y + prompt_h / 2),
        )
        inner_top = prompt_box[1] + 100
        inner_bottom = prompt_box[3] - 58
        inner_mid_y = (inner_top + inner_bottom) / 2
        direct_w = 914
        frag_w = 490
        gap = 94
        total_inner_w = direct_w + gap + frag_w
        inner_x0 = round(visual_mid_x - total_inner_w / 2)
        direct_box = (inner_x0, inner_top, inner_x0 + direct_w, inner_bottom)
        frag_box = (inner_x0 + direct_w + gap, inner_top, inner_x0 + direct_w + gap + frag_w, inner_bottom)
        prompt_label_size = 40
        channel_label_size = 37
        feynman_scale = 0.84
        direct_center_y = inner_mid_y
        frag_center_y = inner_mid_y
        fey_items = [
            (((direct_box[0] * 3 + direct_box[2]) / 4, direct_center_y), "Compton scattering", hp.INK, "compton"),
            (((direct_box[0] + direct_box[2] * 3) / 4, direct_center_y), "Annihilation", hp.INK, "annihilation"),
            (((frag_box[0] + frag_box[2]) / 2, frag_center_y), "Fragmentation radiation", hp.TEAL, "fragmentation"),
        ]
        process_label_band = (inner_bottom, prompt_box[3])
        eq_y = card[1] + 722
        bg_y = card[1] + 824
        eq_font_px = 58
        bg_font_px = 58
        bg_decay_font_px = 66
        eq_font = hp.font(hp.TIMES_BOLD, eq_font_px)
        bg_label = hp.font(hp.TIMES_BOLD, bg_font_px)
        bg_body = hp.font(hp.TIMES, bg_font_px)
        bg_decay = hp.font(hp.TIMES_BOLD, bg_decay_font_px)
    else:
        title_font_px = 72 if not show_subtitle else 46
        title_font = hp.font(hp.TIMES_BOLD, title_font_px)
        title_text = "Production channels"
        title_xy = (card[0] + 48, card[1] + (36 if not show_subtitle else 44))
        draw.text(title_xy, title_text, font=title_font, fill=hp.INK)
        if layout_nodes is not None:
            add_text_node(layout_nodes, name=f"{frame_id} production card title", text=title_text, font_px=title_font_px, bbox=text_bbox(draw, title_text, title_font, *title_xy))
        if show_subtitle:
            subtitle_font_px = 35
            subtitle_font = hp.font(hp.TIMES_ITALIC, subtitle_font_px)
            subtitle_text = "Prompt photons include direct and fragmentation; decay photons are backgrounds."
            subtitle_xy = (card[0] + 50, card[1] + 96)
            draw.text(subtitle_xy, subtitle_text, font=subtitle_font, fill=hp.MUTED)
            if layout_nodes is not None:
                add_text_node(layout_nodes, name=f"{frame_id} production card subtitle", text=subtitle_text, font_px=subtitle_font_px, bbox=text_bbox(draw, subtitle_text, subtitle_font, *subtitle_xy))
        if focused:
            # Tighter "artist centered" geometry for the full-context frame:
            # the two direct diagrams sit on quarter-cell centers of the blue
            # box, while the fragmentation diagram sits on the green-box center.
            prompt_box = (card[0] + 120, card[1] + 190, card[2] - 116, card[1] + 590)
            direct_box = (card[0] + 164, card[1] + 270, card[0] + 744, card[1] + 520)
            frag_box = (card[0] + 840, card[1] + 270, card[2] - 168, card[1] + 520)
        else:
            prompt_box = (card[0] + 76, card[1] + 190, card[2] - 56, card[1] + 590)
            direct_box = (card[0] + 110, card[1] + 270, card[0] + 710, card[1] + 520)
            frag_box = (card[0] + 778, card[1] + 270, card[2] - 80, card[1] + 520)
        prompt_label_size = 35
        channel_label_size = 35
        feynman_scale = 0.74
        direct_center_y = (direct_box[1] + direct_box[3]) / 2
        frag_center_y = (frag_box[1] + frag_box[3]) / 2
        fey_items = [
            (((direct_box[0] * 3 + direct_box[2]) / 4, direct_center_y), "Compton scattering", hp.INK, "compton"),
            (((direct_box[0] + direct_box[2] * 3) / 4, direct_center_y), "Annihilation", hp.INK, "annihilation"),
            (((frag_box[0] + frag_box[2]) / 2, frag_center_y), "Fragmentation radiation", hp.TEAL, "fragmentation"),
        ]
        process_label_band = (direct_box[3], prompt_box[3])
        eq_y = card[3] - 292
        bg_y = card[3] - 186
        eq_font_px = 47
        bg_font_px = 45
        bg_decay_font_px = 54
        eq_font = hp.font(hp.TIMES_BOLD, eq_font_px)
        bg_label = hp.font(hp.TIMES_BOLD, bg_font_px)
        bg_body = hp.font(hp.TIMES, bg_font_px)
        bg_decay = hp.font(hp.TIMES_BOLD, bg_decay_font_px)

    draw.rounded_rectangle(prompt_box, radius=30, fill=(255, 247, 244, 120), outline=(197, 64, 48, 210), width=3)
    centered_label(draw, prompt_box, "Prompt photon production", (197, 64, 48), (255, 247, 244), -38, prompt_label_size)
    draw.rounded_rectangle(direct_box, radius=20, fill=(237, 247, 254, 230), outline=(*hp.SPHENIX_BLUE, 220), width=3)
    draw.rounded_rectangle(frag_box, radius=20, fill=(236, 247, 243, 230), outline=(*hp.TEAL, 220), width=3)
    centered_label(draw, direct_box, "Direct photons", hp.SPHENIX_BLUE, (237, 247, 254), -60, channel_label_size)
    centered_label(draw, frag_box, "Fragmentation photons", hp.TEAL, (236, 247, 243), -60, channel_label_size)
    if layout_nodes is not None:
        layout_nodes.extend(
            [
                {"name": f"{frame_id} production card", "kind": "box", "bbox": box_to_list(card)},
                {"name": f"{frame_id} prompt production box", "kind": "box", "bbox": box_to_list(prompt_box)},
                {"name": f"{frame_id} direct photon box", "kind": "box", "bbox": box_to_list(direct_box)},
                {"name": f"{frame_id} fragmentation photon box", "kind": "box", "bbox": box_to_list(frag_box)},
            ]
        )
        for text, size, box, y_offset in (
            ("Prompt photon production", prompt_label_size, prompt_box, -38),
            ("Direct photons", channel_label_size, direct_box, -60),
            ("Fragmentation photons", channel_label_size, frag_box, -60),
        ):
            fnt = hp.font(hp.TIMES_BOLD, size)
            tw, th = hp.text_box(draw, text, fnt)
            lx0 = round((box[0] + box[2] - tw) / 2) - 22
            label_box = (lx0, box[1] + y_offset, lx0 + 44 + tw, box[1] + y_offset + 56)
            add_text_node(
                layout_nodes,
                name=f"{frame_id} {text} bubble label",
                text=text,
                font_px=size,
                bbox=(label_box[0] + 18, label_box[1] + (56 - th) / 2 - 3, label_box[0] + 18 + tw, label_box[1] + (56 - th) / 2 - 3 + th),
            )

    for (cx, cy), label, color, mode in fey_items:
        draw_feynman_centered(draw, cx, cy, feynman_scale, mode)
        label_font_px = 36 if expanded else 35
        fnt = hp.font(hp.TIMES_BOLD, label_font_px)
        tw, th = hp.text_box(draw, label, fnt)
        label_y = round((process_label_band[0] + process_label_band[1] - th) / 2)
        draw.text((round(cx - tw / 2), label_y), label, font=fnt, fill=color)
        if layout_nodes is not None:
            if mode == "compton":
                cell_box = (direct_box[0], direct_box[1], (direct_box[0] + direct_box[2]) / 2, direct_box[3])
            elif mode == "annihilation":
                cell_box = ((direct_box[0] + direct_box[2]) / 2, direct_box[1], direct_box[2], direct_box[3])
            else:
                cell_box = frag_box
            caption_band = (cell_box[0], process_label_band[0], cell_box[2], process_label_band[1])
            layout_nodes.extend(
                [
                    {
                        "name": f"{frame_id} {mode} feynman cell",
                        "kind": "box",
                        "bbox": box_to_list(cell_box),
                    },
                    {
                        "name": f"{frame_id} {mode} feynman ink bbox",
                        "kind": "diagram_ink",
                        "mode": mode,
                        "bbox": box_to_list(feynman_visual_box(cx, cy, feynman_scale, mode)),
                        "parent": f"{frame_id} {mode} feynman cell",
                    },
                    {
                        "name": f"{frame_id} {mode} process-label band",
                        "kind": "box",
                        "bbox": box_to_list(caption_band),
                    },
                ]
            )
            add_text_node(
                layout_nodes,
                name=f"{frame_id} {mode} process label",
                text=label,
                font_px=label_font_px,
                bbox=(round(cx - tw / 2), label_y, round(cx - tw / 2) + tw, label_y + th),
            )

    eq_bbox = draw_centered_segments(
        draw,
        [
            ("prompt photons", eq_font, (197, 64, 48), 30 if expanded else 10),
            ("=", eq_font, hp.INK, 30 if expanded else 10),
            ("direct photons", eq_font, hp.SPHENIX_BLUE, 32 if expanded else 10),
            ("+", eq_font, hp.INK, 32 if expanded else 10),
            ("fragmentation photons", eq_font, hp.TEAL, 0),
        ],
        eq_y,
        card[0] + 70,
        card[2] - 70,
    )
    bg_label_bbox = draw_centered_segments(
        draw,
        [
            ("dominant backgrounds", bg_label, (38, 44, 52), 18 if expanded else 12),
            ("= meson-decay photons:", bg_body, hp.MUTED, 0),
        ],
        bg_y,
        card[0] + 70,
        card[2] - 70,
    )
    bg_decay_bbox = draw_centered_segments(
        draw,
        [
            ("π^0 → γγ", bg_decay, (38, 44, 52), 56 if expanded else 36),
            ("η → γγ", bg_decay, (38, 44, 52), 0),
        ],
        bg_y + (84 if expanded else 76),
        card[0] + 70,
        card[2] - 70,
    )
    if layout_nodes is not None:
        add_text_node(
            layout_nodes,
            name=f"{frame_id} production prompt summary",
            text="prompt photons = direct photons + fragmentation photons",
            font_px=eq_font_px,
            bbox=eq_bbox,
        )
        add_text_node(
            layout_nodes,
            name=f"{frame_id} production background label",
            text="meson-decay backgrounds: photons from",
            font_px=bg_font_px,
            bbox=bg_label_bbox,
        )
        add_text_node(
            layout_nodes,
            name=f"{frame_id} production decay line",
            text="π^0 → γγ   η → γγ",
            font_px=bg_decay_font_px,
            bbox=bg_decay_bbox,
        )


def draw_cone_card(
    img: Image.Image,
    box: tuple[int, int, int, int],
    *,
    busy: bool,
    show_label: bool = True,
    label_font_px: int = 24,
    rim_y: int | None = None,
) -> None:
    x0, y0, x1, y1 = box
    s = 4
    bw, bh = (x1 - x0) * s, (y1 - y0) * s
    tile = Image.new("RGBA", (bw, bh), (255, 255, 255, 0))
    td = ImageDraw.Draw(tile, "RGBA")
    cx = bw // 2
    color = (197, 64, 48) if busy else hp.PHOTON_DARK
    title_color = (197, 64, 48) if busy else hp.BLUE
    card_fill = (255, 248, 246, 255) if busy else (244, 249, 252, 255)
    cone_fill = (255, 247, 244, 178) if busy else (244, 250, 253, 208)
    label = "non-isolated" if busy else "isolated"

    def sc_pt(pt: tuple[float, float]) -> tuple[float, float]:
        return (pt[0] * s, pt[1] * s)

    def local_polyline(points: list[tuple[float, float]], fill: tuple[int, int, int, int], width: float) -> None:
        td.line([sc_pt(p) for p in points], fill=fill, width=max(1, round(width * s)), joint="curve")

    def local_arrow(start: tuple[float, float], end: tuple[float, float], *, width: float = 2.8, head: float = 8.5) -> None:
        sx, sy = sc_pt(start)
        ex, ey = sc_pt(end)
        td.line((sx, sy, ex, ey), fill=(255, 255, 255, 224), width=round((width + 2.4) * s))
        td.line((sx, sy, ex, ey), fill=(*hp.INK, 235), width=round(width * s))
        angle = math.atan2(ey - sy, ex - sx)
        h = head * s
        left = (ex - h * math.cos(angle - 0.58), ey - h * math.sin(angle - 0.58))
        right = (ex - h * math.cos(angle + 0.58), ey - h * math.sin(angle + 0.58))
        td.polygon([(ex, ey), left, right], fill=(*hp.INK, 235))

    td.rounded_rectangle((0, 0, bw - 1, bh - 1), radius=10 * s, fill=card_fill, outline=(222, 229, 236, 255), width=2 * s)
    if show_label:
        title_font = hp.font(hp.TIMES_BOLD, label_font_px * s)
        probe = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
        tw, _th = hp.text_box(probe, label, title_font)
        td.text(((bw - tw) / 2, 14 * s), label, font=title_font, fill=title_color)

    rim_y = rim_y if rim_y is not None else (108 if show_label else 78)
    local_w = x1 - x0
    local_h = y1 - y0
    apex = (local_w / 2, local_h - 48)
    rim_w = min(70, local_w / 2 - 26)
    td.polygon([sc_pt(apex), sc_pt((apex[0] - rim_w, rim_y)), sc_pt((apex[0] + rim_w, rim_y))], fill=cone_fill)
    td.line((sc_pt(apex), sc_pt((apex[0] - rim_w, rim_y))), fill=(*color, 230), width=4 * s)
    td.line((sc_pt(apex), sc_pt((apex[0] + rim_w, rim_y))), fill=(*color, 230), width=4 * s)
    td.arc(tuple(v * s for v in (apex[0] - rim_w, rim_y - 15, apex[0] + rim_w, rim_y + 15)), 0, 180, fill=(*color, 240), width=4 * s)
    td.arc(tuple(v * s for v in (apex[0] - rim_w, rim_y - 15, apex[0] + rim_w, rim_y + 15)), 180, 360, fill=(*color, 120), width=2 * s)

    photon_start = (apex[0] - 8, rim_y - 48)
    photon_end = (apex[0], apex[1] - 10)
    pts = hp.feynman_points(photon_start, photon_end, 7.3, 4.8, 128)
    local_polyline(pts, (*color, 246), 5.2)
    # Keep the photon arrowhead visually seated on the squiggle rather than drifting right.
    tip = (photon_start[0] - 2.6, photon_start[1] - 6.0)
    td.polygon(
        [sc_pt(tip), sc_pt((tip[0] - 10.0, tip[1] + 16.0)), sc_pt((tip[0] + 7.0, tip[1] + 13.0))],
        fill=(*color, 246),
    )
    gamma_font = hp.font(hp.TIMES_BOLD, 36 * s)
    gamma_xy = sc_pt((tip[0] + 19, tip[1] - 1))
    td.text(gamma_xy, "γ", font=gamma_font, fill=color)

    if busy:
        origin = (apex[0], apex[1] - 8)
        for start, end in (
            (origin, (apex[0] - 55, apex[1] - 52)),
            (origin, (apex[0] - 32, apex[1] - 128)),
            (origin, (apex[0] - 8, apex[1] - 88)),
            (origin, (apex[0] + 34, apex[1] - 126)),
            (origin, (apex[0] + 58, apex[1] - 58)),
        ):
            local_arrow(start, end)

    tile = tile.resize((x1 - x0, y1 - y0), Image.Resampling.LANCZOS)
    img.alpha_composite(tile, (x0, y0))


def draw_isolation_focus_card(
    img: Image.Image,
    card: tuple[int, int, int, int],
    *,
    show_subtitle: bool = True,
    title_font_px: int = 56,
    subtitle_font_px: int = 42,
    layout_nodes: list[dict[str, object]] | None = None,
    frame_id: str = "",
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    draw_card(draw, card, hp.PHOTON, stripe_w=14)
    left = card[0] + 52
    right = card[2] - 58
    title_text = "Isolation definition"
    title_font = hp.font(hp.TIMES_BOLD, title_font_px)
    title_xy = (left, card[1] + 48)
    draw.text(title_xy, title_text, font=title_font, fill=hp.INK)
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} isolation expanded title", text=title_text, font_px=title_font_px, bbox=text_bbox(draw, title_text, title_font, *title_xy))
    if show_subtitle:
        subtitle_text = "Selects photons with little nearby activity."
        subtitle_font = hp.font(hp.TIMES_ITALIC, subtitle_font_px)
        subtitle_xy = (left, card[1] + 116)
        draw.text(subtitle_xy, subtitle_text, font=subtitle_font, fill=hp.MUTED)
        if layout_nodes is not None:
            add_text_node(layout_nodes, name=f"{frame_id} isolation expanded subtitle", text=subtitle_text, font_px=subtitle_font_px, bbox=text_bbox(draw, subtitle_text, subtitle_font, *subtitle_xy))

    eq_box = (left, card[1] + 194, right, card[1] + 344)
    draw.rounded_rectangle(eq_box, radius=9, fill=(255, 251, 239, 255), outline=(238, 203, 128, 255), width=3)
    eq_font_px = 66
    eq_font = hp.font(hp.TIMES_BOLD, eq_font_px)
    eq_text = "E_T^iso = Σ E_T^tower - E_T^candidate"
    eq_w, eq_h = hp.rich_text_box(draw, eq_text, eq_font)
    eq_optical_y_offset = 8
    eq_xy = (
        round((eq_box[0] + eq_box[2] - eq_w) / 2),
        round((eq_box[1] + eq_box[3] - eq_h) / 2) + eq_optical_y_offset,
    )
    hp.draw_rich_text(
        draw,
        eq_xy,
        eq_text,
        eq_font,
        hp.INK,
    )
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} isolation expanded equation", text=eq_text, font_px=eq_font_px, bbox=text_bbox(draw, eq_text, eq_font, *eq_xy, rich=True))

    cue_box = (left + 38, card[1] + 362, right - 38, card[1] + 432)
    cue_text = "small E_T^iso  →  low nearby activity"
    cue_font_px = 44
    cue_font = hp.font(hp.TIMES_BOLD, cue_font_px)
    cue_w, cue_h = hp.rich_text_box(draw, cue_text, cue_font)
    hp.draw_rich_text(
        draw,
        (round((cue_box[0] + cue_box[2] - cue_w) / 2), round((cue_box[1] + cue_box[3] - cue_h) / 2) - 2),
        cue_text,
        cue_font,
        (58, 66, 76),
    )
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} isolation expanded cue", text=cue_text, font_px=cue_font_px, bbox=text_bbox(draw, cue_text, cue_font, round((cue_box[0] + cue_box[2] - cue_w) / 2), round((cue_box[1] + cue_box[3] - cue_h) / 2) - 2, rich=True))

    take_font_px = 42
    take_font = hp.font(hp.TIMES_BOLD, take_font_px)
    line1 = "Suppresses nearby"
    line2 = "jet activity."
    l1w, l1h = hp.text_box(draw, line1, take_font)
    l2w, _ = hp.text_box(draw, line2, take_font)
    take_band = (left + 30, card[1] + 440, right - 30, card[1] + 548)
    y_take = round((take_band[1] + take_band[3] - (l1h * 2 + 6)) / 2)
    draw.text((round((card[0] + card[2] - l1w) / 2), y_take), line1, font=take_font, fill=hp.INK)
    draw.text((round((card[0] + card[2] - l2w) / 2), y_take + l1h + 6), line2, font=take_font, fill=hp.INK)
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} isolation expanded takeaway line 1", text=line1, font_px=take_font_px, bbox=text_bbox(draw, line1, take_font, round((card[0] + card[2] - l1w) / 2), y_take))
        add_text_node(layout_nodes, name=f"{frame_id} isolation expanded takeaway line 2", text=line2, font_px=take_font_px, bbox=text_bbox(draw, line2, take_font, round((card[0] + card[2] - l2w) / 2), y_take + l1h + 6))

    cone_top = card[1] + 642
    cone_bottom = card[3] - 38
    cone_w = 322
    gap = 66
    total = cone_w * 2 + gap
    start_x = round((card[0] + card[2] - total) / 2)
    isolated_box = (start_x, cone_top, start_x + cone_w, cone_bottom)
    nonisolated_box = (start_x + cone_w + gap, cone_top, start_x + cone_w * 2 + gap, cone_bottom)
    label_y = card[1] + 596
    draw_cone_heading(
        draw,
        text="isolated",
        box=isolated_box,
        y=label_y,
        font_px=42,
        fill=hp.BLUE,
        layout_nodes=layout_nodes,
        frame_id=frame_id,
    )
    draw_cone_heading(
        draw,
        text="non-isolated",
        box=nonisolated_box,
        y=label_y,
        font_px=42,
        fill=(197, 64, 48),
        layout_nodes=layout_nodes,
        frame_id=frame_id,
    )
    draw_cone_card(img, isolated_box, busy=False, show_label=False, rim_y=78)
    draw_cone_card(img, nonisolated_box, busy=True, show_label=False, rim_y=78)


def draw_isolation_compact_card(
    img: Image.Image,
    card: tuple[int, int, int, int],
    *,
    layout_nodes: list[dict[str, object]] | None = None,
    frame_id: str = "",
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    draw_card(draw, card, hp.PHOTON, stripe_w=14)

    left = card[0] + 52
    title_text = "Isolation definition"
    title_font_px = 54
    title_font = hp.font(hp.TIMES_BOLD, title_font_px)
    title_xy = (left, card[1] + 18)
    draw.text(title_xy, title_text, font=title_font, fill=hp.INK)
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} isolation compact title", text=title_text, font_px=title_font_px, bbox=text_bbox(draw, title_text, title_font, *title_xy))

    divider_x = card[2] - 390
    eq_box = (left + 14, card[1] + 94, divider_x - 36, card[1] + 186)
    draw.rounded_rectangle(eq_box, radius=8, fill=(255, 251, 239, 255), outline=(238, 203, 128, 255), width=2)
    eq_text = "E_T^iso = Σ E_T^tower - E_T^candidate"
    eq_font_px = 38
    eq_font = hp.font(hp.TIMES_BOLD, eq_font_px)
    eq_w, eq_h = hp.rich_text_box(draw, eq_text, eq_font)
    eq_xy = (round((eq_box[0] + eq_box[2] - eq_w) / 2), round((eq_box[1] + eq_box[3] - eq_h) / 2))
    hp.draw_rich_text(draw, eq_xy, eq_text, eq_font, hp.INK)
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} isolation compact equation", text=eq_text, font_px=eq_font_px, bbox=text_bbox(draw, eq_text, eq_font, *eq_xy, rich=True))

    cue_text = "small E_T^iso  →  low nearby activity"
    cue_font_px = 36
    cue_font = hp.font(hp.TIMES_BOLD, cue_font_px)
    cue_w, cue_h = hp.rich_text_box(draw, cue_text, cue_font)
    cue_box = (left + 8, card[1] + 198, divider_x - 30, card[1] + 246)
    cue_xy = (round((cue_box[0] + cue_box[2] - cue_w) / 2), round((cue_box[1] + cue_box[3] - cue_h) / 2) - 1)
    hp.draw_rich_text(draw, cue_xy, cue_text, cue_font, hp.MUTED)
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} isolation compact cue", text=cue_text, font_px=cue_font_px, bbox=text_bbox(draw, cue_text, cue_font, *cue_xy, rich=True))

    line1 = "Suppresses nearby"
    line2 = "jet activity."
    take_font_px = 38
    take_font = hp.font(hp.TIMES_BOLD, take_font_px)
    l1w, l1h = hp.text_box(draw, line1, take_font)
    l2w, _ = hp.text_box(draw, line2, take_font)
    text_band = (left + 8, card[1] + 244, divider_x - 30, card[1] + 344)
    total_h = l1h * 2 + 4
    y_take = round((text_band[1] + text_band[3] - total_h) / 2)
    x1 = round((text_band[0] + text_band[2] - l1w) / 2)
    x2 = round((text_band[0] + text_band[2] - l2w) / 2)
    draw.text((x1, y_take), line1, font=take_font, fill=hp.INK)
    draw.text((x2, y_take + l1h + 4), line2, font=take_font, fill=hp.INK)
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} isolation compact takeaway line 1", text=line1, font_px=take_font_px, bbox=text_bbox(draw, line1, take_font, x1, y_take))
        add_text_node(layout_nodes, name=f"{frame_id} isolation compact takeaway line 2", text=line2, font_px=take_font_px, bbox=text_bbox(draw, line2, take_font, x2, y_take + l1h + 4))

    draw.line((divider_x, card[1] + 84, divider_x, card[3] - 28), fill=(222, 229, 236, 255), width=2)
    cone_top = card[1] + 112
    cone_bottom = card[3] - 22
    cone_w = 154
    gap = 34
    start_x = divider_x + 28
    isolated_box = (start_x, cone_top, start_x + cone_w, cone_bottom)
    nonisolated_box = (start_x + cone_w + gap, cone_top, start_x + cone_w * 2 + gap, cone_bottom)
    label_y = card[1] + 72
    draw_cone_heading(
        draw,
        text="isolated",
        box=isolated_box,
        y=label_y,
        font_px=35,
        fill=hp.BLUE,
        layout_nodes=layout_nodes,
        frame_id=frame_id,
    )
    draw_cone_heading(
        draw,
        text="non-isolated",
        box=nonisolated_box,
        y=label_y,
        font_px=35,
        fill=(197, 64, 48),
        layout_nodes=layout_nodes,
        frame_id=frame_id,
    )
    draw_cone_card(img, isolated_box, busy=False, show_label=False, rim_y=62)
    draw_cone_card(img, nonisolated_box, busy=True, show_label=False, rim_y=62)


def draw_color_neutral_compact_card(
    img: Image.Image,
    card: tuple[int, int, int, int],
    *,
    layout_nodes: list[dict[str, object]] | None = None,
    frame_id: str = "",
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    draw_card(draw, card, hp.SPHENIX_BLUE, stripe_w=14)
    left = card[0] + 52
    title_text = "Color-neutral behavior"
    title_font_px = 70
    title_font = hp.font(hp.TIMES_BOLD, title_font_px)
    title_xy = (left, card[1] + 32)
    draw.text(title_xy, title_text, font=title_font, fill=hp.INK)
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} color-neutral title", text=title_text, font_px=title_font_px, bbox=text_bbox(draw, title_text, title_font, *title_xy))

    plot_box = (left + 42, card[1] + 142, card[2] - 44, card[1] + 418)
    candidates = [
        hp.ASSET_DIR / "direct_gamma_raa_uncropped_user_20260607.png",
        hp.ASSET_DIR / "direct_gamma_raa_user_constructed_prl109_fig3_backup_slide14.png",
        hp.ASSET_DIR / "phenix_direct_gamma_raa.png",
        hp.ASSET_DIR / "direct_gamma_raa_attached.png",
        hp.ASSET_DIR / "direct_photon_raa.png",
    ]
    for candidate in candidates:
        if candidate.exists():
            hp.paste_fit(img, hp.open_rgba(candidate), plot_box)
            break

    take_line1 = "R_AA near unity:"
    take_line2 = "the photon calibrates the hard scale."
    take_font_px = 42
    take_font = hp.font(hp.TIMES_BOLD, take_font_px)
    take_band = (left + 34, card[1] + 426, card[2] - 34, card[3] - 28)
    l1w, l1h = hp.rich_text_box(draw, take_line1, take_font)
    l2w, l2h = hp.rich_text_box(draw, take_line2, take_font)
    total_h = l1h + l2h + 6
    y0 = round((take_band[1] + take_band[3] - total_h) / 2)
    line1_xy = (round((take_band[0] + take_band[2] - l1w) / 2), y0)
    line2_xy = (round((take_band[0] + take_band[2] - l2w) / 2), y0 + l1h + 6)
    hp.draw_rich_text(draw, line1_xy, take_line1, take_font, hp.INK)
    hp.draw_rich_text(draw, line2_xy, take_line2, take_font, hp.INK)
    if layout_nodes is not None:
        add_text_node(layout_nodes, name=f"{frame_id} color-neutral takeaway line 1", text=take_line1, font_px=take_font_px, bbox=text_bbox(draw, take_line1, take_font, *line1_xy, rich=True))
        add_text_node(layout_nodes, name=f"{frame_id} color-neutral takeaway line 2", text=take_line2, font_px=take_font_px, bbox=text_bbox(draw, take_line2, take_font, *line2_xy, rich=True))


def render_frame1(layout_nodes: list[dict[str, object]] | None = None) -> Path:
    img, _draw = draw_chrome(
        layout_nodes,
        frame_id="slide05",
        title_text="Prompt photons in p+p: production and backgrounds",
    )
    draw_production_card(img, (132, BODY_TOP, 2428, BODY_BOTTOM), expanded=True, layout_nodes=layout_nodes, frame_id="slide05")
    apply_common_footer(img)
    path = OUT_DIR / "slide05_state01_production_expanded.png"
    img.convert("RGB").save(path, "PNG")
    return path


def render_frame2(layout_nodes: list[dict[str, object]] | None = None) -> Path:
    img, _draw = draw_chrome(
        layout_nodes,
        frame_id="slide06",
        title_text="Isolation selects the clean prompt-photon sample",
    )
    draw_production_card(
        img,
        (132, BODY_TOP, 1448, BODY_BOTTOM),
        expanded=False,
        focused=True,
        show_subtitle=False,
        layout_nodes=layout_nodes,
        frame_id="slide06",
    )
    draw_isolation_focus_card(
        img,
        (1496, BODY_TOP, 2428, BODY_BOTTOM),
        show_subtitle=False,
        title_font_px=72,
        subtitle_font_px=40,
        layout_nodes=layout_nodes,
        frame_id="slide06",
    )
    apply_common_footer(img)
    path = OUT_DIR / "slide06_state02_isolation_expanded.png"
    img.convert("RGB").save(path, "PNG")
    return path


def render_frame3(layout_nodes: list[dict[str, object]] | None = None) -> Path:
    path = OUT_DIR / "slide07_state03_current_full_context.png"
    current, _draw = draw_chrome(
        layout_nodes,
        frame_id="slide07",
        title_text="Prompt photons as color-neutral hard probes",
    )
    draw_production_card(
        current,
        (132, BODY_TOP, 1448, BODY_BOTTOM),
        expanded=False,
        focused=True,
        show_subtitle=False,
        layout_nodes=layout_nodes,
        frame_id="slide07",
    )
    draw_isolation_compact_card(
        current,
        (1496, BODY_TOP, 2428, 720),
        layout_nodes=layout_nodes,
        frame_id="slide07",
    )
    draw_color_neutral_compact_card(
        current,
        (1496, 752, 2428, BODY_BOTTOM),
        layout_nodes=layout_nodes,
        frame_id="slide07",
    )
    apply_common_footer(current)
    current.convert("RGB").save(path, "PNG")
    return path


def render_sequence() -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    layout_nodes: list[dict[str, object]] = []
    paths = [render_frame1(layout_nodes), render_frame2(layout_nodes)]
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "source_generator": str(Path(__file__).resolve().relative_to(hp.ROOT)),
        "baseline_generator": "scripts/slides/hp2026/opening/make_hp2026_opening_motivation_slide.py",
        "accepted_current_baseline": str(CURRENT_ACCEPTED_BASELINE.relative_to(hp.ROOT)),
        "hp2026_main_header": {
            "deck": "hp2026_main_talk",
            "title_font_size": 86,
            "subtitle_font_size": None,
            "title_xy": [132, 76],
            "subtitle_xy": None,
            "divider_y": TITLE_DIVIDER_Y,
        },
        "outputs": [str(p.relative_to(hp.ROOT)) for p in paths],
        "sequence_intent": [
            "Frame 1 expands the Production channels card across the main content width and frames prompt photons in p+p as direct plus fragmentation photons, with neutral-meson decay photons as backgrounds.",
            "Frame 2 restores the Production card to a compact left context panel and expands the Isolation definition card on the right.",
        ],
        "notes": [
            "No Google Slides edits were made.",
            "Header matches the Slide 8-10 title-only contract; footer, logo placement, Times typography, card stripe language, and 2560x1440 canvas are preserved.",
            "The same source-native Feynman/cone drawing style is used for the new expanded states.",
            "The previous R_AA/color-neutral behavior frame is intentionally not emitted for this sequence.",
            "The previous campaign audit has hard-coded y=330 card positions and is not used for this y=272 body expansion.",
        ],
    }
    with (OUT_DIR / "manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    with (OUT_DIR / "layout_nodes.json").open("w", encoding="utf-8") as f:
        json.dump(
            {
                "generated_at": manifest["generated_at"],
                "source_generator": manifest["source_generator"],
                "minimum_audience_font_px": 35,
                "nodes": layout_nodes,
            },
            f,
            indent=2,
        )
        f.write("\n")
    return paths


if __name__ == "__main__":
    for p in render_sequence():
        print(p)
