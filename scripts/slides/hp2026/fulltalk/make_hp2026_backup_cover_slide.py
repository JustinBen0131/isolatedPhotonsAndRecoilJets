#!/usr/bin/env python3
"""Render a clean local HP2026 backup-divider cover slide.

This is a local PNG candidate only. It uses the collaboration-provided
lowercase "bACKUP" logo from the sPHENIX Public Material Repository screenshot
as the visual source. The slide is intentionally a visual divider, not another
main-talk content slide.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
SCRIPT_DIR = ROOT / "scripts/slides/hp2026/fulltalk"
sys.path.insert(0, str(SCRIPT_DIR))

import make_hp2026_fulltalk_candidates as ft  # noqa: E402


OUTDIR = ROOT / "outputs/manual-20260610-backup-cover-slide"
ASSET_DIR = OUTDIR / "assets"
PNG = OUTDIR / "slide21_backup_cover_lowercase_logo.png"
SCRIPT_MD = OUTDIR / "slide21_backup_cover_lowercase_logo_script.md"
MANIFEST = OUTDIR / "slide21_backup_cover_lowercase_logo_manifest.json"
HEADER_JSON = PNG.with_suffix(".header.json")

SOURCE_LOGO = ASSET_DIR / "lowercase_logo_box.png"
OFFICIAL_SPHENIX_LOGO = (
    ROOT
    / "outputs/manual-20260601-hp2026-title/presentations/hp2026-title-slide/assets/sPHENIX-logo.pdf.png"
)

TITLE = "Backup"
SUBTITLE = "Additional material for questions and technical details."

HEADER_SPEC = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": 56,
    "title_xy": [132, 76],
    "subtitle_xy": [136, 190],
    "divider_y": 286,
}


def draw_slide_shell() -> Image.Image:
    img = Image.new("RGBA", (ft.W, ft.H), (*ft.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, ft.W, ft.H), fill=(*ft.SOFT_BG, 255))
    draw.rectangle((0, 0, ft.W, 22), fill=(*ft.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, ft.W, 30), fill=(*ft.PHOTON, 255))
    draw.text(tuple(HEADER_SPEC["title_xy"]), TITLE, font=ft.font(ft.TIMES_BOLD, HEADER_SPEC["title_font_size"]), fill=ft.INK)
    draw.text(
        tuple(HEADER_SPEC["subtitle_xy"]),
        SUBTITLE,
        font=ft.font(ft.TIMES_ITALIC, HEADER_SPEC["subtitle_font_size"]),
        fill=ft.MUTED,
    )
    y = HEADER_SPEC["divider_y"]
    draw.line((132, y, ft.W - 132, y), fill=(221, 226, 232, 255), width=3)
    ft.add_top_right_sphenix_logo_like_slide2(img)
    return img


def draw_divider_shell() -> Image.Image:
    """Plain canvas for the backup divider mark."""

    img = Image.new("RGBA", (ft.W, ft.H), (*ft.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, ft.W, ft.H), fill=(255, 255, 255, 255))
    return img


def remove_table_border_and_crop(logo: Image.Image) -> Image.Image:
    """Remove the MediaWiki table border while keeping the real logo content."""

    rgba = logo.convert("RGBA")
    # Crop a few pixels from all edges to remove the gray wiki thumbnail frame.
    rgba = rgba.crop((4, 4, rgba.width - 4, rgba.height - 4))
    return ft.crop_visible(rgba, white_threshold=252, pad=8)


def median_blue_color(img: Image.Image) -> tuple[int, int, int]:
    rgba = img.convert("RGBA")
    blue_pixels = []
    for r, g, b, a in rgba.getdata():
        if a > 0 and b > 120 and g > 80 and b > r + 25 and g > r + 15:
            blue_pixels.append((r, g, b))
    if not blue_pixels:
        return ft.SPHENIX_BLUE
    rs, gs, bs = zip(*blue_pixels)
    return (int(sorted(rs)[len(rs) // 2]), int(sorted(gs)[len(gs) // 2]), int(sorted(bs)[len(bs) // 2]))


def blue_bbox(img: Image.Image) -> tuple[int, int, int, int]:
    rgba = img.convert("RGBA")
    xs: list[int] = []
    ys: list[int] = []
    pix = rgba.load()
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = pix[x, y]
            if a > 0 and b > 120 and g > 80 and b > r + 25 and g > r + 15:
                xs.append(x)
                ys.append(y)
    if not xs:
        return (0, 0, rgba.width, rgba.height)
    return (min(xs), min(ys), max(xs) + 1, max(ys) + 1)


def repair_clipped_thumbnail_logo(logo: Image.Image) -> Image.Image:
    """Repair the public-wiki thumbnail crop before enlarging it.

    The visible wiki thumbnail clips the blue sPHENIX circle at the top/bottom.
    This rebuilds the logo on a slightly taller white canvas, draws a matching
    blue contour behind the original crop, and then converts the result to clean
    blue/black layers over transparency. It avoids the flat edge without
    inventing a new visual style.
    """

    src = remove_table_border_and_crop(logo).convert("RGBA")
    blue = median_blue_color(src)

    # Build clean source masks from the thumbnail.
    blue_src = Image.new("L", src.size, 0)
    black_src = Image.new("L", src.size, 0)
    bp = blue_src.load()
    kp = black_src.load()
    sp = src.load()
    for y in range(src.height):
        for x in range(src.width):
            r, g, b, a = sp[x, y]
            if a == 0:
                continue
            if b > 115 and g > 75 and b > r + 18 and g > r + 8:
                bp[x, y] = max(120, min(255, int((b - r + 65) * 1.2)))
            elif r < 132 and g < 132 and b < 132:
                kp[x, y] = max(150, min(255, int((255 - max(r, g, b)) * 1.4)))

    top_extra = 16
    bottom_extra = 24
    side_pad = 8
    w, h = src.size
    canvas_size = (w + 2 * side_pad, h + top_extra + bottom_extra)
    blue_mask = Image.new("L", canvas_size, 0)
    black_mask = Image.new("L", canvas_size, 0)
    blue_mask.paste(blue_src, (side_pad, top_extra))
    black_mask.paste(black_src, (side_pad, top_extra))

    # Mirror the real edge bands from the official thumbnail. This removes the
    # flat crop without inventing a separate cap or filling the logo internals.
    top_band = blue_src.crop((0, 0, w, top_extra)).transpose(Image.Transpose.FLIP_TOP_BOTTOM)
    bottom_band = blue_src.crop((0, h - bottom_extra, w, h)).transpose(Image.Transpose.FLIP_TOP_BOTTOM)
    blue_mask.paste(top_band, (side_pad, 0), top_band)
    blue_mask.paste(bottom_band, (side_pad, top_extra + h), bottom_band)

    # Clean stair-step artifacts from the low-resolution thumbnail while keeping
    # the official wordmark proportions intact.
    scale = 8
    blue_hi = blue_mask.resize((canvas_size[0] * scale, canvas_size[1] * scale), Image.Resampling.LANCZOS)
    black_hi = black_mask.resize((canvas_size[0] * scale, canvas_size[1] * scale), Image.Resampling.LANCZOS)
    blue_hi = blue_hi.filter(ImageFilter.GaussianBlur(0.35))
    black_hi = black_hi.filter(ImageFilter.GaussianBlur(0.18))

    out_size = blue_hi.size
    out = Image.new("RGBA", out_size, (0, 0, 0, 0))
    blue_layer = Image.new("RGBA", out_size, (*blue, 255))
    black_layer = Image.new("RGBA", out_size, (8, 10, 12, 255))
    out.alpha_composite(Image.composite(blue_layer, Image.new("RGBA", out_size, (0, 0, 0, 0)), blue_hi))
    out.alpha_composite(Image.composite(black_layer, Image.new("RGBA", out_size, (0, 0, 0, 0)), black_hi))

    cleaned = ft.crop_visible(out, white_threshold=252, pad=60)
    debug_path = ASSET_DIR / "lowercase_logo_repaired_clean.png"
    cleaned.save(debug_path)
    return cleaned


def mask_from_threshold(img: Image.Image, kind: str) -> Image.Image:
    rgba = img.convert("RGBA")
    mask = Image.new("L", rgba.size, 0)
    mp = mask.load()
    pix = rgba.load()
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = pix[x, y]
            if a == 0:
                continue
            if kind == "blue" and b > 115 and g > 75 and b > r + 18 and g > r + 8:
                mp[x, y] = 255
            elif kind == "black" and r < 86 and g < 86 and b < 86:
                mp[x, y] = 255
    return mask


def compose_from_official_symbol() -> Image.Image:
    """Recompose the backup mark from higher-resolution source layers.

    The wiki thumbnail gives the exact lowercase-b backup word proportions but
    clips the sPHENIX circle. The title-slide asset has a clean high-resolution
    sPHENIX blue symbol. This function uses the thumbnail for text placement
    and the high-resolution asset for the circle, avoiding flat cropped edges.
    """

    if not SOURCE_LOGO.exists():
        raise FileNotFoundError(f"missing lowercase backup logo crop: {SOURCE_LOGO}")
    if not OFFICIAL_SPHENIX_LOGO.exists():
        raise FileNotFoundError(f"missing sPHENIX logo asset: {OFFICIAL_SPHENIX_LOGO}")

    thumb = Image.open(SOURCE_LOGO).convert("RGBA")
    thumb = thumb.crop((4, 4, thumb.width - 4, thumb.height - 4))
    blue = median_blue_color(thumb)
    bx0, by0, bx1, by1 = blue_bbox(thumb)
    black_mask = mask_from_threshold(thumb, "black")

    # Extract only the white K from the official thumbnail placement. The K
    # lives inside the blue disk; the surrounding white canvas is excluded by
    # the tight source-coordinate box.
    k_mask = Image.new("L", thumb.size, 0)
    kp = k_mask.load()
    pix = thumb.load()
    k_box = (
        max(0, bx0 + int((bx1 - bx0) * 0.08)),
        max(0, by0 + int((by1 - by0) * 0.33)),
        min(thumb.width, bx0 + int((bx1 - bx0) * 0.48)),
        min(thumb.height, by0 + int((by1 - by0) * 0.67)),
    )
    for y in range(k_box[1], k_box[3]):
        for x in range(k_box[0], k_box[2]):
            r, g, b, a = pix[x, y]
            if a > 0 and min(r, g, b) > 232:
                kp[x, y] = 255

    # Logical frame extends the blue symbol beyond the cropped thumbnail.
    top_ext = 10
    bottom_ext = 14
    side_pad = 4
    logical = Image.new(
        "RGBA",
        (thumb.width + 2 * side_pad, thumb.height + top_ext + bottom_ext),
        (0, 0, 0, 0),
    )
    shift = (side_pad, top_ext)

    official = Image.open(OFFICIAL_SPHENIX_LOGO).convert("RGBA")
    official_blue_mask = mask_from_threshold(official, "blue")
    obox = official_blue_mask.getbbox()
    if not obox:
        raise RuntimeError("could not locate blue sPHENIX symbol in official asset")
    symbol = Image.new("RGBA", (obox[2] - obox[0], obox[3] - obox[1]), (0, 0, 0, 0))
    symbol_mask = official_blue_mask.crop(obox)
    symbol.alpha_composite(Image.composite(Image.new("RGBA", symbol.size, (*blue, 255)), Image.new("RGBA", symbol.size, (0, 0, 0, 0)), symbol_mask))

    # Fill the old white N cutout before overlaying the backup K. This removes
    # remnants of the sPHENIX word while retaining the official blue contour.
    fill = Image.new("RGBA", symbol.size, (0, 0, 0, 0))
    fd = ImageDraw.Draw(fill, "RGBA")
    fd.rectangle(
        (
            int(symbol.width * 0.07),
            int(symbol.height * 0.39),
            int(symbol.width * 0.30),
            int(symbol.height * 0.62),
        ),
        fill=(*blue, 255),
    )
    symbol.alpha_composite(fill)

    target_w = bx1 - bx0 + 10
    target_h = by1 - by0 + top_ext + bottom_ext
    symbol = symbol.resize((target_w, target_h), Image.Resampling.LANCZOS)
    logical.alpha_composite(symbol, (bx0 + side_pad - 2, by0))

    scale = 8
    out = logical.resize((logical.width * scale, logical.height * scale), Image.Resampling.LANCZOS)
    black_hi = black_mask.resize((thumb.width * scale, thumb.height * scale), Image.Resampling.LANCZOS)
    k_hi = k_mask.resize((thumb.width * scale, thumb.height * scale), Image.Resampling.LANCZOS)
    text_layer = Image.new("RGBA", out.size, (0, 0, 0, 0))
    black_layer = Image.new("RGBA", black_hi.size, (8, 10, 12, 255))
    white_layer = Image.new("RGBA", k_hi.size, (255, 255, 255, 255))
    text_layer.alpha_composite(Image.composite(black_layer, Image.new("RGBA", black_hi.size, (0, 0, 0, 0)), black_hi), (side_pad * scale, top_ext * scale))
    text_layer.alpha_composite(Image.composite(white_layer, Image.new("RGBA", k_hi.size, (0, 0, 0, 0)), k_hi), (side_pad * scale, top_ext * scale))
    out.alpha_composite(text_layer)
    out = out.filter(ImageFilter.UnsharpMask(radius=0.8, percent=55, threshold=1))
    cleaned = ft.crop_visible(out, white_threshold=252, pad=80)
    (ASSET_DIR / "lowercase_logo_recomposed_official_symbol.png").parent.mkdir(parents=True, exist_ok=True)
    cleaned.save(ASSET_DIR / "lowercase_logo_recomposed_official_symbol.png")
    return cleaned


def blue_symbol_from_official() -> tuple[Image.Image, tuple[int, int, int]]:
    official = Image.open(OFFICIAL_SPHENIX_LOGO).convert("RGBA")
    blue = median_blue_color(official)
    mask = mask_from_threshold(official, "blue")
    obox = mask.getbbox()
    if not obox:
        raise RuntimeError("could not locate blue sPHENIX symbol in official asset")
    mask = mask.crop(obox)
    symbol = Image.new("RGBA", mask.size, (0, 0, 0, 0))
    symbol.alpha_composite(Image.composite(Image.new("RGBA", mask.size, (*blue, 255)), Image.new("RGBA", mask.size, (0, 0, 0, 0)), mask))
    # Remove the original white N gap so the backup K is the only white letter
    # in the blue sector.
    draw = ImageDraw.Draw(symbol, "RGBA")
    draw.rectangle(
        (
            int(symbol.width * 0.05),
            int(symbol.height * 0.37),
            int(symbol.width * 0.31),
            int(symbol.height * 0.63),
        ),
        fill=(*blue, 255),
    )
    return symbol, blue


def compose_pristine_symbol_text() -> Image.Image:
    """Compose a high-resolution clean backup divider mark.

    This sacrifices the low-resolution screenshot letter mask in favor of
    crisp Python-rendered lettering, while preserving the official lowercase-b
    backup spelling and the official sPHENIX symbol geometry.
    """

    symbol, blue = blue_symbol_from_official()
    canvas = Image.new("RGBA", (2400, 1220), (0, 0, 0, 0))
    symbol = symbol.resize((1020, 1020), Image.Resampling.LANCZOS)
    sym_x, sym_y = 1040, 85
    canvas.alpha_composite(symbol, (sym_x, sym_y))

    # DIN Alternate has the clean, engineered geometry closest to the public
    # backup lockup available on this machine. Letter placement is manual so
    # the white K lands in the blue part and UP lands in the center aperture.
    font_path = Path("/System/Library/Fonts/Supplemental/DIN Alternate Bold.ttf")
    if not font_path.exists():
        font_path = Path("/System/Library/Fonts/SFNSMono.ttf")
    font = ImageFont.truetype(str(font_path), 250)
    draw = ImageDraw.Draw(canvas, "RGBA")
    y = 480
    positions = [
        ("b", 250, (8, 10, 12, 255)),
        ("A", 450, (8, 10, 12, 255)),
        ("C", 650, (8, 10, 12, 255)),
        ("K", 940, (255, 255, 255, 255)),
        ("U", 1248, (8, 10, 12, 255)),
        ("P", 1460, (8, 10, 12, 255)),
    ]
    for ch, x, fill in positions:
        draw.text((x, y), ch, font=font, fill=fill)

    cleaned = ft.crop_visible(canvas, white_threshold=252, pad=60)
    cleaned.save(ASSET_DIR / "lowercase_logo_pristine_recomposed_text.png")
    return cleaned.filter(ImageFilter.UnsharpMask(radius=0.6, percent=45, threshold=1))


def cleaned_logo_asset() -> Image.Image:
    return compose_pristine_symbol_text()


def draw_center_panel(img: Image.Image) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    panel = (320, 374, 2240, 1164)
    ft.shadow(img, panel, radius=16)
    draw.rounded_rectangle(panel, radius=18, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)

    # A very quiet collaboration-blue accent keeps the divider slide in the
    # deck language without making the asset look like a generated graphic.
    accent = (panel[0] + 58, panel[1] + 78, panel[0] + 70, panel[3] - 78)
    draw.rounded_rectangle(accent, radius=6, fill=(*ft.SPHENIX_BLUE, 255))

    logo = cleaned_logo_asset()
    logo_box = (560, 480, 2020, 900)
    fitted = ft.fit(logo, logo_box[2] - logo_box[0], logo_box[3] - logo_box[1])
    x = logo_box[0] + (logo_box[2] - logo_box[0] - fitted.width) // 2
    y = logo_box[1] + (logo_box[3] - logo_box[1] - fitted.height) // 2
    img.alpha_composite(fitted, (x, y))

    rule_y = 968
    draw.line((600, rule_y, 1960, rule_y), fill=(222, 229, 237, 255), width=3)
    label = "Backup material"
    label_font = ft.font(ft.TIMES_BOLD, 58)
    tw, th = ft.text_box(draw, label, label_font)
    draw.text(((ft.W - tw) / 2, 1010), label, font=label_font, fill=ft.INK)
    note = "Detailed checks, supporting plots, and technical derivations."
    note_font = ft.font(ft.TIMES_ITALIC, 37)
    nw, nh = ft.text_box(draw, note, note_font)
    draw.text(((ft.W - nw) / 2, 1084), note, font=note_font, fill=ft.MUTED)


def draw_large_backup_mark(img: Image.Image) -> None:
    logo = cleaned_logo_asset()
    # Make the collaboration logo essentially the whole slide. No title,
    # footer, card, shadow, or extra identity chrome.
    logo_box = (110, 80, 2450, 1360)
    fitted = ft.fit(logo, logo_box[2] - logo_box[0], logo_box[3] - logo_box[1])
    x = logo_box[0] + (logo_box[2] - logo_box[0] - fitted.width) // 2
    y = logo_box[1] + (logo_box[3] - logo_box[1] - fitted.height) // 2
    img.alpha_composite(fitted, (x, y))


def write_script() -> None:
    SCRIPT_MD.write_text(
        """# HP2026 Slide 21 Script - Backup

This is the divider into backup material.

If questions come up, I have the supporting checks, additional plots, and technical details here.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    img = draw_divider_shell()
    draw_large_backup_mark(img)
    img.convert("RGB").save(PNG, "PNG")
    HEADER_JSON.write_text(
        json.dumps(
            {
                "divider_slide": {
                    "deck": "hp2026_backup_divider",
                    "exceptions": "No main-talk title/subtitle/footer/chrome; this is a plain enlarged backup-symbol divider.",
                    "canvas": [ft.W, ft.H],
                }
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    write_script()
    MANIFEST.write_text(
        json.dumps(
            {
                "created": datetime.now().isoformat(timespec="seconds"),
                "script": str(Path(__file__).resolve()),
                "output_png": str(PNG),
                "speaker_script": str(SCRIPT_MD),
                "source_logo_crop": str(SOURCE_LOGO),
                "source_page": "https://wiki.sphenix.bnl.gov/index.php?title=Public_Material_Repository",
                "source_note": "Lowercase-b sPHENIX backup slide logo from the Public Material Repository; local crop was taken from Justin's authenticated wiki screenshot because CLI access redirects to SDCC login.",
                "deck_mutation": "none; local Slide 21 replacement PNG candidate only",
                "divider_design": "plain white slide with enlarged lowercase backup symbol only",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG)


if __name__ == "__main__":
    main()
