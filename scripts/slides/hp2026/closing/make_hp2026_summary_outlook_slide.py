#!/usr/bin/env python3
"""Render the HP2026 final main-talk summary/outlook slide."""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
FULLTALK_DIR = ROOT / "scripts/slides/hp2026/fulltalk"
sys.path.insert(0, str(FULLTALK_DIR))

import make_hp2026_fulltalk_candidates as ft  # noqa: E402


OUTPUT_DIR = ROOT / "outputs/manual-20260608-hp2026-closing-slide"
PNG_NATIVE = OUTPUT_DIR / "hp2026_slide21_summary_outlook_closing_2560.png"
PNG_2X = OUTPUT_DIR / "hp2026_slide21_summary_outlook_closing_5120.png"
MANIFEST = OUTPUT_DIR / "hp2026_slide21_summary_outlook_closing_manifest.json"

HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": None,
    "title_xy": [132, 76],
    "subtitle_xy": None,
    "divider_y": 232,
}


def draw_centered_text(
    draw: ImageDraw.ImageDraw,
    text: str,
    y: int,
    fnt,
    fill: tuple[int, int, int],
    x0: int = 0,
    x1: int = ft.W,
) -> tuple[int, int, int, int]:
    tw, th = ft.text_box(draw, text, fnt)
    x = x0 + (x1 - x0 - tw) // 2
    draw.text((x, y), text, font=fnt, fill=fill)
    return (x, y, x + tw, y + th)


def draw_slide_number(base: Image.Image, number: str = "21") -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    fnt = ft.font(ft.TIMES_BOLD, 42)
    tw, th = ft.text_box(draw, number, fnt)
    draw.text((ft.W - 46 - tw, 1372 - th // 2), number, font=fnt, fill=(0, 0, 0))


def draw_flow_arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int]) -> None:
    ft.draw_arrow(draw, start, end, fill=(137, 156, 178), width=5)


def draw_hp2026_main_header(draw: ImageDraw.ImageDraw, title: str, subtitle: str | None = None) -> None:
    draw.rectangle((0, 0, ft.W, 22), fill=ft.SPHENIX_BLUE)
    draw.rectangle((0, 22, ft.W, 30), fill=ft.PHOTON)
    draw.text(tuple(HP2026_MAIN_HEADER["title_xy"]), title, font=ft.font(ft.TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]), fill=ft.INK)
    if subtitle and HP2026_MAIN_HEADER["subtitle_xy"] and HP2026_MAIN_HEADER["subtitle_font_size"]:
        draw.text(
            tuple(HP2026_MAIN_HEADER["subtitle_xy"]),
            subtitle,
            font=ft.font(ft.TIMES_ITALIC, HP2026_MAIN_HEADER["subtitle_font_size"]),
            fill=ft.BLUE,
        )
    y = HP2026_MAIN_HEADER["divider_y"]
    draw.line((132, y, ft.W - 132, y), fill=(221, 226, 232), width=3)


def wrapped_lines(draw: ImageDraw.ImageDraw, text: str, fnt, max_width: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    cur = ""
    for word in words:
        trial = word if not cur else f"{cur} {word}"
        if ft.text_box(draw, trial, fnt)[0] <= max_width:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines


def draw_wrapped_centered_y(
    draw: ImageDraw.ImageDraw,
    text: str,
    box: tuple[int, int, int, int],
    fnt,
    *,
    fill: tuple[int, int, int],
    line_gap: int = 7,
) -> None:
    lines = wrapped_lines(draw, text, fnt, box[2] - box[0])
    heights = [ft.text_box(draw, line, fnt)[1] for line in lines]
    total = sum(heights) + line_gap * max(0, len(lines) - 1)
    y = box[1] + (box[3] - box[1] - total) // 2
    for line, height in zip(lines, heights):
        draw.text((box[0], y), line, font=fnt, fill=fill)
        y += height + line_gap


def draw_outlook_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    accent: tuple[int, int, int],
    step: str,
    title: str,
    bullets: list[str],
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rectangle((x0 + 2, y0 + 2, x0 + 20, y1 - 2), fill=(*accent, 255))

    badge = (x0 + 46, y0 + 42, x0 + 108, y0 + 104)
    draw.ellipse(badge, fill=(248, 251, 253, 255), outline=(*accent, 230), width=3)
    draw_centered_text(draw, step, y0 + 57, ft.font(ft.TIMES_BOLD, 31), accent, badge[0], badge[2])

    title_x = x0 + 128
    title_font = ft.font(ft.TIMES_BOLD, 43)
    title_h = ft.text_box(draw, title, title_font)[1]
    title_y = badge[1] + (badge[3] - badge[1] - title_h) // 2 - 1
    ft.draw_wrapped(
        draw,
        title,
        (title_x, title_y),
        x1 - title_x - 40,
        title_font,
        fill=ft.INK,
        line_gap=7,
    )

    draw.line((x0 + 46, y0 + 142, x1 - 42, y0 + 142), fill=(221, 228, 236), width=2)

    row_top = y0 + 176
    row_bottom = y1 - 70
    row_h = (row_bottom - row_top) / max(1, len(bullets))
    body_font = ft.font(ft.TIMES, 41)
    for idx, bullet in enumerate(bullets):
        cy0 = int(row_top + idx * row_h)
        cy1 = int(row_top + (idx + 1) * row_h)
        marker_y = (cy0 + cy1) // 2
        draw.rounded_rectangle((x0 + 54, marker_y - 12, x0 + 78, marker_y + 12), radius=6, fill=(*accent, 235))
        draw_wrapped_centered_y(
            draw,
            bullet,
            (x0 + 102, cy0, x1 - 52, cy1),
            body_font,
            fill=ft.MUTED,
            line_gap=8,
        )


def render() -> Image.Image:
    img = Image.new("RGBA", (ft.W, ft.H), (*ft.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw_hp2026_main_header(draw, "Summary and outlook")
    ft.add_top_right_sphenix_logo_like_slide2(img)

    cards = [
        (
            (132, 298, 844, 1084),
            ft.PHOTON_DARK,
            "1",
            "Extend the p+p reference",
            ["Additional p+p luminosity", "Improved precision", "Broader kinematic reach"],
        ),
        (
            (924, 298, 1636, 1084),
            ft.SPHENIX_BLUE,
            "2",
            "Move into A+A systems",
            ["Train and apply the photon-ID BDT in Au+Au", "Compare isolated photons to the p+p reference", "Test medium effects in nuclear collisions"],
        ),
        (
            (1716, 298, 2428, 1084),
            ft.TEAL,
            "3",
            "Probe the QGP with γ+jet",
            ["Use the photon-tagged recoil-jet topology", "Calibrate the hard scale with the photon", "Study medium modification of the recoil jet"],
        ),
    ]

    mid_y = 691
    draw_flow_arrow(draw, (866, mid_y), (900, mid_y))
    draw_flow_arrow(draw, (1658, mid_y), (1692, mid_y))

    for card in cards:
        draw_outlook_card(img, *card)

    draw.line((256, 1126, ft.W - 256, 1126), fill=(224, 230, 237), width=2)
    draw_centered_text(draw, "Thank you", 1182, ft.font(ft.TIMES_BOLD, 78), ft.INK)

    ft.draw_hp2026_identity_footer(img)
    return img


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    img = render()
    img.convert("RGB").save(PNG_NATIVE, "PNG")
    PNG_NATIVE.with_suffix(".header.json").write_text(
        json.dumps({"hp2026_main_header": HP2026_MAIN_HEADER}, indent=2) + "\n",
        encoding="utf-8",
    )
    img.convert("RGB").resize((ft.W * 2, ft.H * 2), Image.Resampling.LANCZOS).save(PNG_2X, "PNG")
    MANIFEST.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().isoformat(timespec="seconds"),
                "script": str(Path(__file__).resolve()),
                "outputs": {"native_2560": str(PNG_NATIVE), "retina_5120": str(PNG_2X)},
                "deck_target": "HPslides_v1, final main-talk slide immediately before Backup",
                "content": {
                    "title": "Summary and outlook",
                    "subtitle": "First p+p isolated prompt photons establish the RHIC reference for A+A and gamma+jet.",
                    "panels": [
                        "Extend the p+p reference",
                        "Move into A+A systems",
                        "Probe the QGP with gamma+jet",
                    ],
                    "close": "Thank you",
                },
                "hp2026_main_header": HP2026_MAIN_HEADER,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG_2X)


if __name__ == "__main__":
    main()
