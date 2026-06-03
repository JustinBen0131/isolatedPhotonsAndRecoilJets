#!/usr/bin/env python3
"""Render HP2026 title-slide PNG candidates.

This intentionally avoids generative text/logo rendering. Logos are loaded from
verified local assets, and all talk text is drawn deterministically.
"""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


W, H = 2560, 1440

TITLE = "sPHENIX measurement of isolated prompt photon production in p+p collisions at 200 GeV"
TITLE_LINES = [
    "sPHENIX measurement of isolated",
    "prompt photon production in",
    "p+p collisions at 200 GeV",
]
AUTHOR = "Justin Bennett"
AFFILIATION = "University of Illinois Urbana-Champaign"
COLLAB = "on behalf of the sPHENIX Collaboration"
EVENT = "Hard Probes 2026"
VENUE = "Vanderbilt University, Nashville, TN"
DATE = "June 24, 2026"
DOE_LINE_1 = "U.S. Department of Energy"
DOE_LINE_2 = "Office of Science"


ROOT = next(p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").exists())
DEFAULT_WORKSPACE = ROOT / "outputs/manual-20260601-hp2026-title/presentations/hp2026-title-slide"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
ARIAL = FONT_DIR / "Arial.ttf"
ARIAL_BOLD = FONT_DIR / "Arial Bold.ttf"

INK = (18, 22, 28)
MUTED = (71, 79, 91)
SOFT = (244, 246, 248)
BLUE = (19, 41, 75)
SPHENIX_BLUE = (30, 143, 214)
ILLINI_ORANGE = (255, 95, 5)
BNL_TEAL = (16, 92, 120)
SPHENIX_MAROON = (179, 14, 57)
PHOTON = (246, 185, 35)


def font(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(path), size)


def text_size(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def draw_lines(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    lines: list[str],
    fnt: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    line_gap: int,
) -> int:
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        _, h = text_size(draw, line, fnt)
        y += h + line_gap
    return y


def open_rgba(path: Path) -> Image.Image:
    return Image.open(path).convert("RGBA")


def crop_visible(img: Image.Image, white_threshold: int = 246) -> Image.Image:
    rgba = img.convert("RGBA")
    px = rgba.load()
    min_x, min_y, max_x, max_y = rgba.width, rgba.height, -1, -1
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = px[x, y]
            visible = a > 8 and not (r >= white_threshold and g >= white_threshold and b >= white_threshold)
            if visible:
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                max_x = max(max_x, x)
                max_y = max(max_y, y)
    if max_x < min_x:
        return rgba
    pad = 8
    return rgba.crop(
        (
            max(0, min_x - pad),
            max(0, min_y - pad),
            min(rgba.width, max_x + pad + 1),
            min(rgba.height, max_y + pad + 1),
        )
    )


def whiten_to_alpha(img: Image.Image, threshold: int = 248) -> Image.Image:
    rgba = img.convert("RGBA")
    data = []
    for r, g, b, a in rgba.getdata():
        if r >= threshold and g >= threshold and b >= threshold:
            data.append((r, g, b, 0))
        else:
            data.append((r, g, b, a))
    rgba.putdata(data)
    return rgba


def fit(img: Image.Image, max_w: int, max_h: int) -> Image.Image:
    scale = min(max_w / img.width, max_h / img.height)
    size = (max(1, int(img.width * scale)), max(1, int(img.height * scale)))
    return img.resize(size, Image.Resampling.LANCZOS)


def paste_fit(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int], anchor: str = "center") -> None:
    x0, y0, x1, y1 = box
    fitted = fit(img, x1 - x0, y1 - y0)
    if anchor == "left":
        x = x0
    elif anchor == "right":
        x = x1 - fitted.width
    else:
        x = x0 + ((x1 - x0) - fitted.width) // 2
    y = y0 + ((y1 - y0) - fitted.height) // 2
    base.alpha_composite(fitted, (x, y))


def load_assets(asset_dir: Path) -> dict[str, Image.Image]:
    hp = whiten_to_alpha(crop_visible(open_rgba(asset_dir / "hp2026_indico_logo.png")))
    sphenix = crop_visible(open_rgba(asset_dir / "sPHENIX_logo_maroon_transparent.png"))
    sphenix_blue_path = asset_dir / "sphenix-logo-white-bg_0.png"
    if sphenix_blue_path.exists():
        sphenix_blue = crop_visible(open_rgba(sphenix_blue_path), white_threshold=252)
    else:
        sphenix_blue = sphenix
    bnl_source = asset_dir / "bnl-logo-2021-squarewrap.svg.png"
    bnl = whiten_to_alpha(crop_visible(open_rgba(bnl_source), white_threshold=250))
    illinois = crop_visible(open_rgba(asset_dir / "illinois_logo_fullcolor_rgb.png"))
    doe_logo_path = asset_dir / "doe_logo_from_yeonju_title.png"
    doe = crop_visible(open_rgba(doe_logo_path), white_threshold=252) if doe_logo_path.exists() else None
    overview = open_rgba(asset_dir / "hp2026_overview_image.png")
    photon_accent_path = asset_dir.parent / "output/photon_detector_motif_best/hp2026_title_accent_best_photon_hits_detector.png"
    photon_accent = open_rgba(photon_accent_path) if photon_accent_path.exists() else None
    assets = {
        "hp": hp,
        "sphenix": sphenix,
        "sphenix_blue": sphenix_blue,
        "bnl": bnl,
        "illinois": illinois,
        "overview": overview,
    }
    if doe is not None:
        assets["doe"] = doe
    if photon_accent is not None:
        assets["photon_accent"] = photon_accent
    return assets


def background() -> Image.Image:
    img = Image.new("RGBA", (W, H), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)
    draw.rectangle((0, 0, W, H), fill=(252, 253, 254))
    return img


def add_footer_branding(base: Image.Image, assets: dict[str, Image.Image], compact: bool = False) -> None:
    draw = ImageDraw.Draw(base)
    y = 1268 if not compact else 1284
    draw.line((132, y - 24, W - 132, y - 24), fill=(221, 226, 232), width=3)
    paste_fit(base, assets["hp"], (132, y, 330, y + 94), "left")
    paste_fit(base, assets["sphenix"], (390, y + 6, 720, y + 94), "left")
    paste_fit(base, assets["bnl"], (785, y + 2, 1195, y + 94), "left")
    paste_fit(base, assets["illinois"], (1270, y + 5, 1336, y + 92), "left")
    draw.text((1354, y + 17), AFFILIATION, font=font(TIMES, 36), fill=BLUE)
    draw.text((1960, y + 17), EVENT, font=font(TIMES_BOLD, 38), fill=INK)
    draw.text((1960, y + 60), f"{VENUE} | {DATE}", font=font(TIMES, 28), fill=MUTED)


def add_footer_branding_refined(base: Image.Image, assets: dict[str, Image.Image], variant: str = "full") -> None:
    draw = ImageDraw.Draw(base)
    y = 1270
    draw.line((132, y - 22, W - 132, y - 22), fill=(221, 226, 232), width=3)
    paste_fit(base, assets["sphenix_blue"], (132, y - 2, 500, y + 98), "left")
    paste_fit(base, assets["bnl"], (602, y - 1, 1018, y + 96), "left")
    paste_fit(base, assets["illinois"], (1138, y + 4, 1205, y + 92), "left")
    draw.text((1228, y + 18), AFFILIATION, font=font(TIMES, 34), fill=BLUE)
    if variant == "full":
        paste_fit(base, assets["hp"], (1940, y + 2, 2098, y + 94), "left")
        draw.text((2122, y + 13), EVENT, font=font(TIMES_BOLD, 35), fill=INK)
        draw.text((2122, y + 53), DATE, font=font(TIMES, 28), fill=MUTED)


def add_footer_branding_with_doe_text(base: Image.Image, assets: dict[str, Image.Image]) -> None:
    draw = ImageDraw.Draw(base)
    y = 1268
    draw.line((132, y - 22, W - 132, y - 22), fill=(221, 226, 232), width=3)

    paste_fit(base, assets["sphenix_blue"], (132, y - 2, 458, y + 98), "left")
    paste_fit(base, assets["bnl"], (596, y - 1, 1006, y + 98), "left")
    paste_fit(base, assets["illinois"], (1156, y + 4, 1222, y + 92), "left")
    draw.text((1246, y + 18), AFFILIATION, font=font(TIMES, 32), fill=BLUE)

    draw.line((1842, y + 6, 1842, y + 88), fill=(223, 229, 236), width=2)
    draw.text((1904, y + 14), DOE_LINE_1, font=font(TIMES_BOLD, 31), fill=INK)
    draw.text((1904, y + 54), DOE_LINE_2, font=font(TIMES, 28), fill=MUTED)


def add_footer_branding_logo_rail(base: Image.Image, assets: dict[str, Image.Image]) -> None:
    draw = ImageDraw.Draw(base)
    rail_top = 1248
    logo_top = 1264
    logo_bottom = 1368
    draw.line((132, rail_top, W - 132, rail_top), fill=(221, 226, 232), width=3)

    paste_fit(base, assets["sphenix_blue"], (404, logo_top + 2, 732, logo_bottom - 2), "center")
    paste_fit(base, assets["bnl"], (870, logo_top - 2, 1215, logo_bottom + 2), "center")
    paste_fit(base, assets["illinois"], (1356, logo_top - 2, 1450, logo_bottom + 2), "center")
    if "doe" in assets:
        paste_fit(base, assets["doe"], (1584, logo_top - 4, 2076, logo_bottom + 4), "center")
    else:
        draw.text((1584, logo_top + 18), DOE_LINE_1, font=font(TIMES_BOLD, 34), fill=INK)
        draw.text((1584, logo_top + 60), DOE_LINE_2, font=font(TIMES, 28), fill=MUTED)

    for x in (804, 1288, 1518):
        draw.line((x, logo_top + 14, x, logo_bottom - 14), fill=(224, 229, 235), width=2)


def add_footer_branding_logo_rail_no_illinois(base: Image.Image, assets: dict[str, Image.Image]) -> None:
    draw = ImageDraw.Draw(base)
    rail_top = 1248
    logo_top = 1264
    logo_bottom = 1368
    draw.line((132, rail_top, W - 132, rail_top), fill=(221, 226, 232), width=3)

    paste_fit(base, assets["sphenix_blue"], (500, logo_top + 2, 840, logo_bottom - 2), "center")
    paste_fit(base, assets["bnl"], (970, logo_top - 2, 1328, logo_bottom + 2), "center")
    if "doe" in assets:
        paste_fit(base, assets["doe"], (1484, logo_top - 4, 1994, logo_bottom + 4), "center")
    else:
        draw.text((1484, logo_top + 18), DOE_LINE_1, font=font(TIMES_BOLD, 34), fill=INK)
        draw.text((1484, logo_top + 60), DOE_LINE_2, font=font(TIMES, 28), fill=MUTED)

    for x in (906, 1404):
        draw.line((x, logo_top + 14, x, logo_bottom - 14), fill=(224, 229, 235), width=2)


def add_footer_branding_logo_rail_final(base: Image.Image, assets: dict[str, Image.Image]) -> None:
    draw = ImageDraw.Draw(base)
    rail_top = 1238
    logo_top = 1246
    logo_bottom = 1398
    draw.line((132, rail_top, W - 132, rail_top), fill=(224, 229, 235), width=2)

    paste_fit(base, assets["sphenix_blue"], (190, logo_top, 760, logo_bottom), "center")
    paste_fit(base, assets["bnl"], (950, logo_top + 6, 1550, logo_bottom - 6), "center")
    if "doe" in assets:
        paste_fit(base, assets["doe"], (1752, logo_top - 2, 2388, logo_bottom + 2), "center")
    else:
        draw.text((1800, logo_top + 26), DOE_LINE_1, font=font(TIMES_BOLD, 39), fill=INK)
        draw.text((1800, logo_top + 74), DOE_LINE_2, font=font(TIMES, 32), fill=MUTED)

    for x in (820, 1686):
        draw.line((x, logo_top + 16, x, logo_bottom - 16), fill=(232, 236, 240), width=2)


def add_photo_panel(
    base: Image.Image,
    assets: dict[str, Image.Image],
    box: tuple[int, int, int, int],
    tint_alpha: int = 48,
    bottom_fade: int = 260,
) -> None:
    draw = ImageDraw.Draw(base)
    x0, y0, x1, y1 = box
    w, h = x1 - x0, y1 - y0
    right = crop_cover(assets["overview"], w, h).convert("RGBA")
    right = Image.alpha_composite(right, Image.new("RGBA", right.size, (255, 255, 255, tint_alpha)))
    mask = Image.new("L", right.size, 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, w, h), radius=36, fill=255)
    if bottom_fade:
        shade = Image.new("RGBA", (w, h), (0, 0, 0, 0))
        sd = ImageDraw.Draw(shade)
        for i in range(bottom_fade):
            alpha = int(118 * (i / bottom_fade) ** 1.55)
            sd.line((0, h - 1 - i, w, h - 1 - i), fill=(255, 255, 255, alpha), width=1)
        right = Image.alpha_composite(right, shade)
    base.paste(right, (x0, y0), mask)
    draw.rounded_rectangle(box, radius=36, outline=(221, 226, 232), width=3)


def draw_event_lockup(
    base: Image.Image,
    assets: dict[str, Image.Image],
    box: tuple[int, int, int, int],
    *,
    logo_side: str = "left",
    fill: tuple[int, int, int, int] = (255, 255, 255, 232),
    outline: tuple[int, int, int] = (225, 230, 236),
) -> None:
    draw = ImageDraw.Draw(base)
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=22, fill=fill, outline=outline, width=2)
    if logo_side == "right":
        paste_fit(base, assets["hp"], (x1 - 198, y0 + 20, x1 - 26, y1 - 20), "right")
        tx = x0 + 30
    else:
        paste_fit(base, assets["hp"], (x0 + 24, y0 + 18, x0 + 188, y1 - 18), "left")
        tx = x0 + 212
    draw.text((tx, y0 + 25), EVENT, font=font(TIMES_BOLD, 42), fill=INK)
    draw.text((tx, y0 + 76), VENUE, font=font(TIMES, 29), fill=MUTED)
    draw.text((tx, y0 + 113), DATE, font=font(TIMES, 29), fill=MUTED)


def draw_blue_sphenix_chip(
    base: Image.Image,
    assets: dict[str, Image.Image],
    box: tuple[int, int, int, int],
    *,
    fill: tuple[int, int, int, int] = (255, 255, 255, 224),
) -> None:
    draw = ImageDraw.Draw(base)
    draw.rounded_rectangle(box, radius=20, fill=fill, outline=(226, 231, 237), width=2)
    paste_fit(base, assets["sphenix_blue"], (box[0] + 22, box[1] + 10, box[2] - 22, box[3] - 10), "left")


def title_block_refined(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 82) -> int:
    draw.text((x, y - 56), "sPHENIX at RHIC", font=font(TIMES_ITALIC, 38), fill=SPHENIX_BLUE)
    end_y = draw_lines(draw, (x, y), TITLE_LINES, font(TIMES_BOLD, size), INK, 18)
    draw.text((x, end_y + 36), AUTHOR, font=font(TIMES_BOLD, 48), fill=BLUE)
    draw.text((x, end_y + 92), AFFILIATION, font=font(TIMES, 36), fill=MUTED)
    draw.text((x, end_y + 136), COLLAB, font=font(TIMES_ITALIC, 32), fill=MUTED)
    return end_y + 178


def title_block_refined_with_illinois_mark(
    base: Image.Image,
    draw: ImageDraw.ImageDraw,
    assets: dict[str, Image.Image],
    x: int,
    y: int,
    size: int = 82,
) -> int:
    draw.text((x, y - 56), "sPHENIX at RHIC", font=font(TIMES_ITALIC, 38), fill=SPHENIX_BLUE)
    end_y = draw_lines(draw, (x, y), TITLE_LINES, font(TIMES_BOLD, size), INK, 18)
    draw.text((x, end_y + 36), AUTHOR, font=font(TIMES_BOLD, 48), fill=BLUE)

    affiliation_y = end_y + 92
    affiliation_font = font(TIMES, 36)
    draw.text((x, affiliation_y), AFFILIATION, font=affiliation_font, fill=MUTED)
    affiliation_w, _ = text_size(draw, AFFILIATION, affiliation_font)
    paste_fit(base, assets["illinois"], (x + affiliation_w + 18, affiliation_y - 8, x + affiliation_w + 58, affiliation_y + 44), "left")

    draw.text((x, end_y + 136), COLLAB, font=font(TIMES_ITALIC, 32), fill=MUTED)
    return end_y + 178


def normalize_white_background(img: Image.Image, bg: tuple[int, int, int] = (252, 253, 254)) -> Image.Image:
    rgba = img.convert("RGBA")
    data = []
    for r, g, b, a in rgba.getdata():
        if a > 0 and r >= 246 and g >= 246 and b >= 246:
            data.append((*bg, a))
        else:
            data.append((r, g, b, a))
    rgba.putdata(data)
    return rgba


def add_photon_title_accent(base: Image.Image, assets: dict[str, Image.Image]) -> None:
    if "photon_accent" not in assets:
        return
    accent = normalize_white_background(assets["photon_accent"])
    accent = fit(accent, 1232, 300)
    x = 142
    y = 792
    base.alpha_composite(accent, (x, y))


def cubic_points(
    start: tuple[float, float],
    c1: tuple[float, float],
    c2: tuple[float, float],
    end: tuple[float, float],
    steps: int = 160,
) -> list[tuple[float, float]]:
    points = []
    for i in range(steps):
        t = i / (steps - 1)
        mt = 1.0 - t
        x = mt**3 * start[0] + 3 * mt**2 * t * c1[0] + 3 * mt * t**2 * c2[0] + t**3 * end[0]
        y = mt**3 * start[1] + 3 * mt**2 * t * c1[1] + 3 * mt * t**2 * c2[1] + t**3 * end[1]
        points.append((x, y))
    return points


def draw_polyline(draw: ImageDraw.ImageDraw, points: list[tuple[float, float]], fill: tuple[int, int, int, int], width: int) -> None:
    draw.line([(int(x), int(y)) for x, y in points], fill=fill, width=width, joint="curve")


def add_integrated_photon_detector_accent(base: Image.Image) -> None:
    overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    center_y = 906
    for i, (rx, alpha, width) in enumerate(((176, 28, 4), (250, 24, 3), (332, 20, 3), (428, 16, 2))):
        box = (1060 - rx, center_y - rx * 0.56, 1060 + rx, center_y + rx * 0.56)
        draw.arc(box, start=292, end=68, fill=(30, 143, 214, alpha), width=width)
        if i < 3:
            box2 = (520 - rx, center_y - rx * 0.58, 520 + rx, center_y + rx * 0.58)
            draw.arc(box2, start=112, end=248, fill=(30, 143, 214, max(10, alpha - 8)), width=max(1, width - 1))

    start = (198, center_y + 2)
    end = (1128, center_y - 8)
    path = cubic_points(start, (432, center_y - 24), (730, center_y + 30), end, 190)
    wave = []
    for idx, (x, y) in enumerate(path):
        t = idx / (len(path) - 1)
        envelope = math.sin(math.pi * t) ** 0.45
        wave.append((x, y + 17 * envelope * math.sin(16 * math.pi * t)))

    glow = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    glow_draw = ImageDraw.Draw(glow)
    draw_polyline(glow_draw, wave, (248, 194, 67, 78), 18)
    draw_polyline(glow_draw, wave, (250, 210, 100, 56), 34)
    overlay = Image.alpha_composite(overlay, glow.filter(ImageFilter.GaussianBlur(4.0)))
    draw = ImageDraw.Draw(overlay)
    draw_polyline(draw, wave, (237, 171, 41, 150), 6)
    draw_polyline(draw, wave, (255, 222, 126, 210), 2)

    source = (start[0], start[1])
    for r, alpha in ((30, 34), (18, 70), (7, 210)):
        draw.ellipse((source[0] - r, source[1] - r, source[0] + r, source[1] + r), fill=(250, 197, 70, alpha))
    for k in range(12):
        angle = 2 * math.pi * k / 12
        r0, r1 = 18, 46
        x0 = source[0] + r0 * math.cos(angle)
        y0 = source[1] + r0 * math.sin(angle)
        x1 = source[0] + r1 * math.cos(angle)
        y1 = source[1] + r1 * math.sin(angle)
        draw.line((x0, y0, x1, y1), fill=(250, 197, 70, 34), width=2)

    face_x, face_y = 1122, 792
    tile_w, tile_h = 42, 34
    rows, cols = 6, 5
    for c in range(cols):
        for r in range(rows):
            x = face_x + c * 39 + r * 8
            y = face_y + r * 40 - c * 7
            dx = 11 + c * 2
            near_impact = abs(r - 3) + abs(c - 1) <= 2
            fill = (229, 237, 247, 118)
            if near_impact:
                fill = (250, 184, 48, 122 + 18 * max(0, 3 - abs(r - 3) - abs(c - 1)))
            outline = (255, 255, 255, 148)
            poly = [(x, y), (x + tile_w, y - dx), (x + tile_w + 10, y + tile_h), (x + 8, y + tile_h + dx)]
            draw.polygon(poly, fill=fill, outline=outline)

    impact = (1143, center_y - 10)
    for r, alpha in ((64, 30), (42, 58), (22, 120), (7, 230)):
        draw.ellipse((impact[0] - r, impact[1] - r, impact[0] + r, impact[1] + r), fill=(252, 196, 59, alpha))
    for k in range(18):
        angle = 2 * math.pi * k / 18
        length = 54 if k % 3 else 72
        draw.line(
            (
                impact[0],
                impact[1],
                impact[0] + length * math.cos(angle),
                impact[1] + length * math.sin(angle),
            ),
            fill=(250, 205, 83, 46),
            width=2,
        )

    detector_edge = [(1328, 785), (1378, 812), (1392, 1032), (1322, 1060)]
    draw.line(detector_edge, fill=(30, 143, 214, 26), width=4, joint="curve")
    draw.line([(1306, 810), (1352, 836), (1362, 1008), (1300, 1034)], fill=(30, 143, 214, 18), width=3)
    base.alpha_composite(overlay)


def add_photon_arc(base: Image.Image, start: tuple[int, int], end: tuple[int, int], color=PHOTON) -> None:
    overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    sx, sy = start
    ex, ey = end
    for i in range(7):
        t = i / 6
        x = sx * (1 - t) + ex * t
        y = sy * (1 - t) + ey * t - 120 * math.sin(math.pi * t)
        r = int(20 + 14 * math.sin(math.pi * t))
        alpha = int(44 + 38 * math.sin(math.pi * t))
        draw.ellipse((x - r, y - r, x + r, y + r), fill=(*color, alpha))
    points = []
    for i in range(80):
        t = i / 79
        x = sx * (1 - t) + ex * t
        y = sy * (1 - t) + ey * t - 120 * math.sin(math.pi * t)
        points.append((x, y))
    draw.line(points, fill=(*color, 118), width=7, joint="curve")
    base.alpha_composite(overlay.filter(ImageFilter.GaussianBlur(0.4)))


def add_emcal_grid(base: Image.Image, origin: tuple[int, int], cols: int, rows: int, scale: int, alpha: int) -> None:
    overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    ox, oy = origin
    colors = [
        (231, 236, 243, alpha),
        (246, 185, 35, alpha + 8),
        (255, 95, 5, alpha),
        (22, 115, 150, alpha - 6),
    ]
    for r in range(rows):
        for c in range(cols):
            x = ox + c * scale + r * int(scale * 0.32)
            y = oy + r * int(scale * 0.74)
            jitter = ((c * 17 + r * 11) % 4) - 1
            col = colors[(c + 2 * r) % len(colors)]
            box = (x, y + jitter, x + int(scale * 0.72), y + int(scale * 0.56) + jitter)
            draw.rounded_rectangle(box, radius=6, fill=col, outline=(255, 255, 255, min(160, alpha + 50)), width=2)
    base.alpha_composite(overlay)


def add_analysis_path(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base)
    items = [
        (["EMCal", "clusters"], 930),
        (["isolation"], 1190),
        (["purity"], 1430),
        (["unfolding"], 1660),
    ]
    y = 846
    card_w = 196
    card_h = 106
    for i, (label_lines, x) in enumerate(items):
        fill = (255, 255, 255, 230)
        outline = (226, 231, 237, 255)
        accent = [PHOTON, ILLINI_ORANGE, SPHENIX_MAROON, BNL_TEAL][i]
        draw.rounded_rectangle((x, y, x + card_w, y + card_h), radius=26, fill=fill, outline=outline, width=2)
        draw.ellipse((x + 18, y + 24, x + 48, y + 54), fill=accent)
        label_y = y + 58 if len(label_lines) == 1 else y + 51
        for line in label_lines:
            tw, th = text_size(draw, line, font(TIMES_BOLD, 25))
            draw.text((x + card_w / 2 - tw / 2, label_y), line, font=font(TIMES_BOLD, 25), fill=INK)
            label_y += th + 4
        if i < len(items) - 1:
            next_x = items[i + 1][1]
            draw.line((x + card_w + 10, y + 52, next_x - 12, y + 52), fill=(155, 164, 178), width=3)
            draw.polygon(
                [
                    (next_x - 12, y + 52),
                    (next_x - 28, y + 43),
                    (next_x - 28, y + 61),
                ],
                fill=(155, 164, 178),
            )


def add_cross_section_motif(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base)
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=22, fill=(255, 255, 255, 222), outline=(225, 230, 236), width=2)
    ax0, ay0, ax1, ay1 = x0 + 44, y1 - 46, x1 - 28, y0 + 42
    draw.line((ax0, ay0, ax1, ay0), fill=(92, 101, 115), width=3)
    draw.line((ax0, ay0, ax0, ay1), fill=(92, 101, 115), width=3)
    points = [(ax0 + 4, ay0 - 8), (ax0 + 60, ay0 - 36), (ax0 + 116, ay0 - 82), (ax0 + 172, ay0 - 132), (ax0 + 228, ay0 - 188)]
    draw.line(points, fill=SPHENIX_MAROON, width=6)
    for x, y in points:
        draw.ellipse((x - 7, y - 7, x + 7, y + 7), fill=SPHENIX_MAROON)
    draw.text((x0 + 38, y0 + 18), "isolated prompt photon", font=font(TIMES_BOLD, 26), fill=INK)
    draw.text((x0 + 38, y0 + 52), "cross section", font=font(TIMES_ITALIC, 24), fill=MUTED)


def title_block(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 92) -> int:
    draw.text((x, y - 58), "sPHENIX at RHIC", font=font(TIMES_ITALIC, 38), fill=SPHENIX_MAROON)
    end_y = draw_lines(draw, (x, y), TITLE_LINES, font(TIMES_BOLD, size), INK, 18)
    draw.text((x, end_y + 36), AUTHOR, font=font(TIMES_BOLD, 48), fill=BLUE)
    draw.text((x, end_y + 92), AFFILIATION, font=font(TIMES, 36), fill=MUTED)
    draw.text((x, end_y + 136), COLLAB, font=font(TIMES_ITALIC, 32), fill=MUTED)
    return end_y + 178


def candidate_a(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_emcal_grid(base, (1510, 210), 9, 7, 112, 42)
    add_photon_arc(base, (1380, 800), (2230, 310), PHOTON)
    draw.rounded_rectangle((128, 112, W - 128, H - 126), radius=0, outline=(230, 234, 239), width=3)
    paste_fit(base, assets["hp"], (131, 120, 365, 250), "left")
    paste_fit(base, assets["sphenix"], (W - 620, 126, W - 132, 248), "right")
    title_block(draw, 185, 382, 94)
    draw.text((185, 1040), EVENT, font=font(TIMES_BOLD, 42), fill=INK)
    draw.text((185, 1092), f"{VENUE} | {DATE}", font=font(TIMES, 34), fill=MUTED)
    add_footer_branding(base, assets)
    base.convert("RGB").save(path, quality=95)


def candidate_b(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_emcal_grid(base, (105, 820), 8, 5, 96, 45)
    add_photon_arc(base, (650, 830), (1775, 470), PHOTON)
    title_block(draw, 160, 218, 84)
    add_analysis_path(base)
    add_cross_section_motif(base, (1888, 602, 2350, 954))
    draw.text((1815, 504), "From isolated EMCal clusters", font=font(TIMES_ITALIC, 31), fill=MUTED)
    draw.text((1815, 544), "to a corrected photon cross section.", font=font(TIMES_ITALIC, 31), fill=MUTED)
    paste_fit(base, assets["hp"], (2020, 122, 2320, 258), "right")
    paste_fit(base, assets["sphenix"], (1610, 126, 1950, 248), "right")
    add_footer_branding(base, assets)
    base.convert("RGB").save(path, quality=95)


def crop_cover(img: Image.Image, box_w: int, box_h: int) -> Image.Image:
    src_ratio = img.width / img.height
    dst_ratio = box_w / box_h
    if src_ratio > dst_ratio:
        new_w = int(img.height * dst_ratio)
        x0 = (img.width - new_w) // 2
        img = img.crop((x0, 0, x0 + new_w, img.height))
    else:
        new_h = int(img.width / dst_ratio)
        y0 = (img.height - new_h) // 2
        img = img.crop((0, y0, img.width, y0 + new_h))
    return img.resize((box_w, box_h), Image.Resampling.LANCZOS)


def candidate_c(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    right = crop_cover(assets["overview"], 910, 1040).convert("RGBA")
    tint = Image.new("RGBA", right.size, (255, 255, 255, 58))
    right = Image.alpha_composite(right, tint)
    mask = Image.new("L", right.size, 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, right.width, right.height), radius=36, fill=255)
    base.paste(right, (1538, 186), mask)
    draw.rounded_rectangle((1538, 186, 2448, 1226), radius=36, outline=(221, 226, 232), width=3)
    shade = Image.new("RGBA", (910, 1040), (0, 0, 0, 0))
    sd = ImageDraw.Draw(shade)
    sd.rectangle((0, 0, 910, 1040), fill=(255, 255, 255, 0))
    for i in range(260):
        alpha = int(105 * (i / 260) ** 1.6)
        sd.line((0, 1039 - i, 910, 1039 - i), fill=(255, 255, 255, alpha), width=1)
    base.paste(Image.alpha_composite(right, shade), (1538, 186), mask)
    paste_fit(base, assets["hp"], (2180, 222, 2390, 338), "right")
    add_photon_arc(base, (1060, 1045), (1850, 610), PHOTON)
    title_block(draw, 160, 280, 82)
    draw.text((160, 1060), EVENT, font=font(TIMES_BOLD, 42), fill=INK)
    draw.text((160, 1112), f"{VENUE} | {DATE}", font=font(TIMES, 34), fill=MUTED)
    paste_fit(base, assets["sphenix"], (1560, 1044, 1940, 1140), "left")
    add_footer_branding(base, assets)
    base.convert("RGB").save(path, quality=95)


def candidate_c1_left_lockup(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_photo_panel(base, assets, (1538, 186, 2448, 1226), tint_alpha=42, bottom_fade=280)
    title_block_refined(draw, 160, 280, 82)
    draw_event_lockup(base, assets, (160, 940, 960, 1096), logo_side="left")
    draw_blue_sphenix_chip(base, assets, (1580, 1048, 1964, 1148))
    add_footer_branding_refined(base, assets, variant="logos")
    base.convert("RGB").save(path, quality=95)


def candidate_c2_photo_lockup(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_photo_panel(base, assets, (1538, 186, 2448, 1226), tint_alpha=44, bottom_fade=260)
    draw_event_lockup(base, assets, (1584, 232, 2388, 388), logo_side="right", fill=(255, 255, 255, 226))
    title_block_refined(draw, 160, 280, 82)
    draw_blue_sphenix_chip(base, assets, (1584, 1048, 1978, 1148))
    add_footer_branding_refined(base, assets, variant="logos")
    base.convert("RGB").save(path, quality=95)


def candidate_c3_bottom_caption(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_photo_panel(base, assets, (1514, 164, 2456, 1218), tint_alpha=36, bottom_fade=330)
    title_block_refined(draw, 160, 274, 82)
    draw.rectangle((1518, 1012, 2452, 1214), fill=(255, 255, 255, 212))
    paste_fit(base, assets["hp"], (1556, 1044, 1728, 1154), "left")
    draw.text((1764, 1044), EVENT, font=font(TIMES_BOLD, 43), fill=INK)
    draw.text((1764, 1096), VENUE, font=font(TIMES, 30), fill=MUTED)
    draw.text((1764, 1136), DATE, font=font(TIMES, 30), fill=MUTED)
    add_footer_branding_refined(base, assets, variant="logos")
    base.convert("RGB").save(path, quality=95)


def candidate_c3_clean_vanderbilt_doe(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_photo_panel(base, assets, (1514, 164, 2456, 1218), tint_alpha=24, bottom_fade=0)
    title_block_refined(draw, 160, 274, 82)

    draw.rectangle((1518, 1012, 2452, 1214), fill=(255, 255, 255, 248))
    draw.line((1518, 1012, 2452, 1012), fill=(221, 226, 232), width=2)
    paste_fit(base, assets["hp"], (1556, 1044, 1728, 1154), "left")
    draw.text((1764, 1044), EVENT, font=font(TIMES_BOLD, 43), fill=INK)
    draw.text((1764, 1096), VENUE, font=font(TIMES, 30), fill=MUTED)
    draw.text((1764, 1136), DATE, font=font(TIMES, 30), fill=MUTED)

    add_footer_branding_with_doe_text(base, assets)
    base.convert("RGB").save(path, quality=95)


def candidate_c3_professional_logo_rail(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_photo_panel(base, assets, (1514, 164, 2456, 1218), tint_alpha=24, bottom_fade=0)
    title_block_refined(draw, 160, 274, 82)

    draw.rectangle((1518, 1012, 2452, 1214), fill=(255, 255, 255, 248))
    draw.line((1518, 1012, 2452, 1012), fill=(221, 226, 232), width=2)
    paste_fit(base, assets["hp"], (1556, 1044, 1728, 1154), "left")
    draw.text((1764, 1044), EVENT, font=font(TIMES_BOLD, 43), fill=INK)
    draw.text((1764, 1096), VENUE, font=font(TIMES, 30), fill=MUTED)
    draw.text((1764, 1136), DATE, font=font(TIMES, 30), fill=MUTED)

    add_footer_branding_logo_rail(base, assets)
    base.convert("RGB").save(path, quality=95)


def candidate_c3_integrated_photon_accent(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_photo_panel(base, assets, (1514, 164, 2456, 1218), tint_alpha=24, bottom_fade=0)
    title_block_refined_with_illinois_mark(base, draw, assets, 160, 274, 82)
    add_photon_title_accent(base, assets)

    draw.rectangle((1518, 1012, 2452, 1214), fill=(255, 255, 255, 248))
    draw.line((1518, 1012, 2452, 1012), fill=(221, 226, 232), width=2)
    paste_fit(base, assets["hp"], (1556, 1044, 1728, 1154), "left")
    draw.text((1764, 1044), EVENT, font=font(TIMES_BOLD, 43), fill=INK)
    draw.text((1764, 1096), VENUE, font=font(TIMES, 30), fill=MUTED)
    draw.text((1764, 1136), DATE, font=font(TIMES, 30), fill=MUTED)

    add_footer_branding_logo_rail_no_illinois(base, assets)
    base.convert("RGB").save(path, quality=95)


def draw_event_caption_band_final(base: Image.Image, assets: dict[str, Image.Image]) -> None:
    draw = ImageDraw.Draw(base)
    band = (1518, 1008, 2452, 1214)
    draw.rectangle(band, fill=(255, 255, 255, 250))
    draw.line((band[0], band[1], band[2], band[1]), fill=(222, 227, 233), width=2)
    paste_fit(base, assets["hp"], (1560, 1042, 1730, 1154), "left")
    draw.text((1772, 1039), EVENT, font=font(TIMES_BOLD, 44), fill=INK)
    draw.text((1772, 1093), VENUE, font=font(TIMES, 30), fill=MUTED)
    draw.text((1772, 1134), DATE, font=font(TIMES, 30), fill=MUTED)
    draw.rounded_rectangle((1514, 164, 2456, 1218), radius=36, outline=(221, 226, 232), width=3)


def candidate_c4_final_institutional_photon_polish(path: Path, assets: dict[str, Image.Image]) -> None:
    base = background()
    draw = ImageDraw.Draw(base)
    add_photo_panel(base, assets, (1514, 164, 2456, 1218), tint_alpha=22, bottom_fade=0)
    title_block_refined_with_illinois_mark(base, draw, assets, 160, 274, 82)
    add_integrated_photon_detector_accent(base)
    draw_event_caption_band_final(base, assets)
    add_footer_branding_logo_rail_final(base, assets)
    base.convert("RGB").save(path, quality=95)


def write_manifest(out_dir: Path, asset_dir: Path, files: list[Path]) -> None:
    manifest = {
        "title": TITLE,
        "size_px": [W, H],
        "outputs": [str(p) for p in files],
        "assets": {
            "hp2026_logo": {
                "path": str(asset_dir / "hp2026_indico_logo.png"),
                "source": "https://indico.cern.ch/event/1428985/logo-134929478.png",
            },
            "hp2026_overview_image": {
                "path": str(asset_dir / "hp2026_overview_image.png"),
                "source": "https://indico.cern.ch/event/1428985/attachments/3252768/5813862/image.png",
            },
            "sphenix_logo": {
                "path": str(asset_dir / "sPHENIX_logo_maroon_transparent.png"),
                "source": "/Users/patsfan753/Downloads/sPHENIX_logo_maroon_transparent.png",
            },
            "sphenix_blue_logo": {
                "path": str(asset_dir / "sphenix-logo-white-bg_0.png"),
                "source": "https://www.sphenix.bnl.gov/sites/default/files/sphenix-logo-white-bg_0.png",
            },
            "bnl_logo": {
                "path": str(asset_dir / "bnl-logo-2021.svg"),
                "source": "https://www.bnl.gov/assets/global/images/bnl-logo-2021.svg",
            },
            "illinois_block_i": {
                "path": str(asset_dir / "illinois_logo_fullcolor_rgb.png"),
                "source": "https://illinois.edu/wp-content/uploads/2025/10/Illinois_logo_fullcolor_rgb.png",
            },
            "doe_logo_reference_crop": {
                "path": str(asset_dir / "doe_logo_from_yeonju_title.png"),
                "source": "/Users/patsfan753/Desktop/ThesisAnalysis/usefulDocs/20260506_DIS_YeonjuGo.pdf page 1",
            },
            "doe_footer_svg_reference": {
                "path": str(asset_dir / "doe_footer_logo.svg"),
                "source": "https://www.energy.gov/themes/custom/energy_gov/svg/footer-logo.svg",
            },
            "photon_detector_title_accent": {
                "path": str(asset_dir.parent / "output/photon_detector_motif_best/hp2026_title_accent_best_photon_hits_detector.png"),
                "source": "built-in image_gen output copied and toned locally under output/photon_detector_motif_best",
            },
        },
        "web_references": {
            "sphenix_logo_page": {
                "path": str(asset_dir / "sphenix_logo_page.html"),
                "source": "https://www.sphenix.bnl.gov/logo",
            },
            "hp2024_public_contributions": {
                "path": str(asset_dir.parent / "webrefs/hp2024_contributions.html"),
                "source": "https://indico.cern.ch/event/1339555/contributions/",
            },
            "hp2020_public_contributions": {
                "path": str(asset_dir.parent / "webrefs/hp2020_contributions.html"),
                "source": "https://indico.cern.ch/event/751767/contributions/",
            },
        },
        "notes": [
            "No Google Slides mutation.",
            "No slide number or visible provenance footer baked into candidates.",
            "BNL SVG is also stored as a source asset; a local PNG render is used for placement.",
            "C-refinement variants use the official blue sPHENIX logo asset from sphenix.bnl.gov; the original C baseline is retained for comparison.",
            "The yellow photon arc/circles from the original C baseline are intentionally omitted from the C-refinement variants.",
            "Public HP reference pages supported the layout choice to group event identity, venue, and date as a coherent lockup instead of repeating conference text in multiple places.",
            "The C3 clean Vanderbilt/DOE-footer refinement includes the DOE Office of Science as a text label, not an unverified logo asset.",
            "The professional logo-rail refinement uses the DOE wordmark cropped from Yeonju Go's DIS 2026 title-slide reference, rather than AI-generated or text-recreated branding.",
            "The UIUC footer affiliation text is intentionally removed; the footer row uses the Block I logo only.",
            "The integrated photon-accent refinement places the photon-hit image in the left-side white space under the author/collaboration block and moves the Illinois Block I mark beside the UIUC affiliation line, leaving the footer as sPHENIX/BNL/DOE only.",
            "The C4 final polish pass uses a native-drawn photon/detector accent instead of pasting the standalone motif, keeps the official HP2026/Vanderbilt panel, and keeps the footer to sPHENIX/BNL/DOE.",
        ],
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, default=DEFAULT_WORKSPACE)
    parser.add_argument(
        "--only",
        choices=[
            "all",
            "c3-clean-vanderbilt-doe",
            "c3-professional-logo-rail",
            "c3-integrated-photon-accent",
            "c4-final-institutional-photon-polish",
        ],
        default="all",
        help="Render only one named final candidate when requested.",
    )
    args = parser.parse_args()
    asset_dir = args.workspace / "assets"
    out_dir = args.workspace / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    assets = load_assets(asset_dir)
    candidates = [
        ("all", out_dir / "hp2026_title_A_public_institutional_clean.png", candidate_a),
        ("all", out_dir / "hp2026_title_B_photon_story_motif.png", candidate_b),
        ("all", out_dir / "hp2026_title_C_yeonju_inspired_public_layout.png", candidate_c),
        ("all", out_dir / "hp2026_title_C1_left_event_lockup_blue_sphenix.png", candidate_c1_left_lockup),
        ("all", out_dir / "hp2026_title_C2_photo_event_lockup_blue_sphenix.png", candidate_c2_photo_lockup),
        ("all", out_dir / "hp2026_title_C3_photo_caption_blue_sphenix.png", candidate_c3_bottom_caption),
        (
            "c3-clean-vanderbilt-doe",
            out_dir / "hp2026_title_C3_clean_vanderbilt_doe_footer.png",
            candidate_c3_clean_vanderbilt_doe,
        ),
        (
            "c3-professional-logo-rail",
            out_dir / "hp2026_title_C3_professional_logo_rail.png",
            candidate_c3_professional_logo_rail,
        ),
        (
            "c3-integrated-photon-accent",
            out_dir / "hp2026_title_C3_integrated_photon_accent.png",
            candidate_c3_integrated_photon_accent,
        ),
        (
            "c4-final-institutional-photon-polish",
            out_dir / "hp2026_title_C4_final_institutional_photon_polish.png",
            candidate_c4_final_institutional_photon_polish,
        ),
    ]
    selected = candidates if args.only == "all" else [c for c in candidates if c[0] == args.only]
    files = [path for _, path, _ in selected]
    for _, path, render in selected:
        render(path, assets)
    write_manifest(out_dir, asset_dir, files)
    for path in files:
        print(path)


if __name__ == "__main__":
    main()
