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
    ft.draw_arrow(draw, start, end, fill=(151, 169, 188), width=4)


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
    draw.rounded_rectangle((x0 + 18, y0 + 28, x0 + 30, y1 - 28), radius=6, fill=(*accent, 255))

    badge = (x0 + 54, y0 + 34, x0 + 108, y0 + 88)
    draw.ellipse(badge, fill=(248, 251, 253, 255), outline=(*accent, 230), width=3)
    draw_centered_text(draw, step, y0 + 46, ft.font(ft.TIMES_BOLD, 25), accent, badge[0], badge[2])

    title_x = x0 + 132
    ft.draw_wrapped(
        draw,
        title,
        (title_x, y0 + 36),
        x1 - title_x - 40,
        ft.font(ft.TIMES_BOLD, 35),
        fill=ft.INK,
        line_gap=5,
    )

    draw.line((x0 + 54, y0 + 122, x1 - 42, y0 + 122), fill=(221, 228, 236), width=2)

    y = y0 + 166
    for bullet in bullets:
        draw.rounded_rectangle((x0 + 74, y + 8, x0 + 92, y + 26), radius=4, fill=(*accent, 235))
        y = ft.draw_wrapped(
            draw,
            bullet,
            (x0 + 112, y),
            x1 - x0 - 162,
            ft.font(ft.TIMES, 30),
            fill=ft.MUTED,
            line_gap=7,
        )
        y += 38


def render() -> Image.Image:
    subtitle = (
        "The first sPHENIX p+p isolated prompt-photon cross section agrees with pQCD "
        "and anchors the future photon and γ+jet program at RHIC."
    )
    img = ft.base_slide("Summary and outlook", subtitle)
    ft.add_top_right_sphenix_logo_like_slide2(img)
    draw = ImageDraw.Draw(img, "RGBA")

    cards = [
        (
            (132, 408, 824, 938),
            ft.PHOTON_DARK,
            "1",
            "Extend the p+p reference",
            ["Additional p+p luminosity", "Improved precision and kinematic reach"],
        ),
        (
            (934, 408, 1626, 938),
            ft.SPHENIX_BLUE,
            "2",
            "Move into A+A systems",
            ["Train and apply the photon-ID BDT in Au+Au", "Compare isolated-photon results to the p+p reference"],
        ),
        (
            (1736, 408, 2428, 938),
            ft.TEAL,
            "3",
            "Probe the QGP with γ+jet",
            ["Use the photon-tagged recoil-jet topology", "Study medium modification with a calibrated hard scale"],
        ),
    ]

    mid_y = 673
    draw_flow_arrow(draw, (846, mid_y), (910, mid_y))
    draw_flow_arrow(draw, (1648, mid_y), (1712, mid_y))

    for card in cards:
        draw_outlook_card(img, *card)

    draw.line((256, 1018, ft.W - 256, 1018), fill=(224, 230, 237), width=2)
    draw_centered_text(draw, "Thank you", 1112, ft.font(ft.TIMES_BOLD, 64), ft.INK)

    ft.draw_hp2026_identity_footer(img)
    draw_slide_number(img, "21")
    return img


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    img = render()
    img.convert("RGB").save(PNG_NATIVE, "PNG")
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
                    "subtitle": "The first sPHENIX p+p isolated prompt-photon cross section agrees with pQCD and anchors the future photon and γ+jet program at RHIC.",
                    "panels": [
                        "Extend the p+p reference",
                        "Move into A+A systems",
                        "Probe the QGP with gamma+jet",
                    ],
                    "close": "Thank you",
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG_2X)


if __name__ == "__main__":
    main()
