#!/usr/bin/env python3
"""
Slide 4: Run 24 p+p dataset used for this measurement.

Local PNG generator for HP2026 full-talk slide 4.

This version removes the run-coordinator luminosity plot and the old
"key distinction" bookkeeping card.  The slide now defines only the analysis
sample and luminosity needed for the 15-minute talk narrative.

Output:
  outputs/manual-20260607-run24pp-lumi-asset/
      hp2026_slide04_analysis_dataset_context_v13_lint_expanded.png
      hp2026_slide04_analysis_dataset_context_v14_cards_expanded.png
      hp2026_slide04_analysis_dataset_context_v15_cleaner_cards.png
      hp2026_slide04_analysis_dataset_context_v16_lumi_card_aligned.png
      hp2026_slide04_analysis_dataset_context_v17_no_sample_row.png
      hp2026_slide04_analysis_dataset_context_v18_large_row_labels.png
      hp2026_slide04_analysis_dataset_context_v19_no_rhs_carry_note.png
      hp2026_slide04_analysis_dataset_context_v20_logo_matched.png
"""
from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont

W, H = 2560, 1440

HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": None,
    "title_xy": [132, 76],
    "subtitle_xy": None,
    "divider_y": 232,
}

ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
OUTPUT_DIR   = ROOT / "outputs/manual-20260607-run24pp-lumi-asset"
ASSET_DIR    = ROOT / "outputs/manual-20260601-hp2026-title/presentations/hp2026-title-slide/assets"

FONT_DIR   = Path("/System/Library/Fonts/Supplemental")
TIMES      = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITAL = FONT_DIR / "Times New Roman Italic.ttf"

# ── shared HP2026 palette ────────────────────────────────────────────────────
INK          = (18,  22,  28)
BLUE         = (19,  41,  75)
SPHENIX_BLUE = (30, 143, 214)
PHOTON       = (245, 181,  34)
PHOTON_DARK  = (196, 124,  10)
TEAL_SOFT    = (102, 173, 190)
MUTED        = (71,  79,  91)
LIGHT_MUTED  = (117, 126, 138)
SOFT_BG      = (252, 253, 254)
PANEL_EDGE   = (218, 226, 235)
FOOTER_RULE  = (217, 225, 233)
RED          = (192,  32,  32)

# ── FIX: amber card palette (replaces blue-tinted fill on the right card) ───
AMBER_FILL   = (255, 250, 237)   # warm pale amber
AMBER_BORDER = (238, 218, 176)   # warm amber outline


# ── font helpers ─────────────────────────────────────────────────────────────

def fnt(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(path), size)


def tbox(draw: ImageDraw.ImageDraw, text: str, f: ImageFont.ImageFont) -> tuple[int, int]:
    b = draw.textbbox((0, 0), text, font=f)
    return b[2] - b[0], b[3] - b[1]


def _fsize(f: ImageFont.ImageFont) -> int:
    return int(getattr(f, "size", 24))


def _script_fnt(f: ImageFont.ImageFont) -> ImageFont.FreeTypeFont:
    size = max(10, round(_fsize(f) * 0.58))
    return fnt(TIMES_BOLD if getattr(f, "path", "") == str(TIMES_BOLD) else TIMES, size)


def _next_run(text: str, idx: int) -> tuple[str, int]:
    idx += 1
    if idx >= len(text):
        return "", idx
    if text[idx] == "{":
        end = text.find("}", idx + 1)
        if end != -1:
            return text[idx + 1: end], end + 1
    start = idx
    if text[idx] in "+-":
        idx += 1
    while idx < len(text) and text[idx].isalnum():
        idx += 1
    if idx == start:
        idx += 1
    return text[start:idx], idx


def rich(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    f: ImageFont.ImageFont,
    fill: tuple[int, int, int],
) -> int:
    """Render text with ^ superscript markers; return x-advance."""
    sf = _script_fnt(f)
    x, y = xy
    i = 0
    while i < len(text):
        if text[i] in "^_":
            marker = text[i]
            run, i = _next_run(text, i)
            dy = -round(_fsize(f) * 0.22) if marker == "^" else round(_fsize(f) * 0.36)
            draw.text((x, y + dy), run, font=sf, fill=fill)
            x += tbox(draw, run, sf)[0]
        else:
            start = i
            while i < len(text) and text[i] not in "^_":
                i += 1
            run = text[start:i]
            draw.text((x, y), run, font=f, fill=fill)
            x += tbox(draw, run, f)[0]
    return x


def wrapped(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_w: int,
    f: ImageFont.ImageFont,
    fill: tuple[int, int, int] = INK,
    gap: int = 8,
) -> int:
    words = text.split()
    lines: list[str] = []
    cur = ""
    for word in words:
        trial = word if not cur else f"{cur} {word}"
        if tbox(draw, trial, f)[0] <= max_w:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=f, fill=fill)
        y += tbox(draw, line, f)[1] + gap
    return y


def _wrap_lines(
    draw: ImageDraw.ImageDraw,
    text: str,
    max_w: int,
    f: ImageFont.ImageFont,
) -> list[str]:
    words = text.split()
    lines: list[str] = []
    cur = ""
    for word in words:
        trial = word if not cur else f"{cur} {word}"
        if tbox(draw, trial, f)[0] <= max_w:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines


def wrapped_vcenter(
    draw: ImageDraw.ImageDraw,
    text: str,
    box: tuple[int, int, int, int],
    f: ImageFont.ImageFont,
    fill: tuple[int, int, int] = INK,
    gap: int = 7,
) -> None:
    """Draw wrapped text vertically centered inside a bounded box."""
    lines = _wrap_lines(draw, text, box[2] - box[0], f)
    if not lines:
        return
    heights = [tbox(draw, line, f)[1] for line in lines]
    total_h = sum(heights) + gap * (len(lines) - 1)
    y = box[1] + ((box[3] - box[1]) - total_h) // 2
    for line, h in zip(lines, heights):
        draw.text((box[0], y), line, font=f, fill=fill)
        y += h + gap


# ── image helpers ─────────────────────────────────────────────────────────────

def fit_img(img: Image.Image, max_w: int, max_h: int) -> Image.Image:
    s = min(max_w / img.width, max_h / img.height)
    return img.resize((max(1, int(img.width * s)), max(1, int(img.height * s))), Image.Resampling.LANCZOS)


def crop_vis(img: Image.Image, thr: int = 248) -> Image.Image:
    rgba = img.convert("RGBA")
    px = rgba.load()
    mnx, mny, mxx, mxy = rgba.width, rgba.height, -1, -1
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = px[x, y]
            if a > 8 and not (r >= thr and g >= thr and b >= thr):
                mnx = min(mnx, x); mny = min(mny, y)
                mxx = max(mxx, x); mxy = max(mxy, y)
    if mxx < mnx:
        return rgba
    p = 8
    return rgba.crop((max(0, mnx - p), max(0, mny - p),
                      min(rgba.width, mxx + p + 1), min(rgba.height, mxy + p + 1)))


def w2t(img: Image.Image, thr: int = 248) -> Image.Image:
    """White-to-transparent for logos."""
    rgba = img.convert("RGBA")
    out = Image.new("RGBA", rgba.size, (0, 0, 0, 0))
    src = rgba.load(); dst = out.load()
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = src[x, y]
            if a == 0:
                continue
            w = min(r, g, b)
            if w >= thr:
                continue
            alpha = a if w <= thr - 18 else int(a * (thr - w) / 18)
            dst[x, y] = (r, g, b, alpha)
    return out


def drop_shadow(base: Image.Image, box: tuple[int, int, int, int], r: int = 12) -> None:
    layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer, "RGBA")
    d.rounded_rectangle((box[0] + 8, box[1] + 10, box[2] + 8, box[3] + 10),
                         radius=r, fill=(30, 42, 58, 30))
    base.alpha_composite(layer.filter(ImageFilter.GaussianBlur(12)))


def round_panel(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    fill: tuple[int, int, int] = (255, 255, 255),
    edge: tuple[int, int, int] = PANEL_EDGE,
    r: int = 12,
) -> None:
    draw.rounded_rectangle(box, radius=r, fill=(*fill, 255), outline=(*edge, 255), width=2)


def side_accent(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    color: tuple[int, int, int],
    width: int = 14,
) -> None:
    """Draw a solid left-edge accent without panel-border seams."""
    r = width // 2
    x0 = box[0] - 1
    x1 = box[0] + width
    y0 = box[1]
    y1 = box[3]
    draw.rectangle((x0, y0 + r, x1, y1 - r), fill=(*color, 255))
    draw.ellipse((x0, y0, x1, y0 + width), fill=(*color, 255))
    draw.ellipse((x0, y1 - width, x1, y1), fill=(*color, 255))


# ── slide sections ────────────────────────────────────────────────────────────

def _draw_header(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))
    draw.text((132, 76), "Run 24 p+p dataset used for this measurement",
              font=fnt(TIMES_BOLD, 86), fill=INK)
    draw.line((132, 232, W - 132, 232), fill=(221, 226, 232), width=3)

    logo_path = ASSET_DIR / "sphenix-logo-white-bg_0.png"
    if logo_path.exists():
        # Match the logo geometry used on the following prompt-photon and
        # shower-shape slides exactly: same crop, fit box, top y, and right x.
        logo = fit_img(crop_vis(Image.open(logo_path).convert("RGBA"), thr=252), 366, 107)
        base.alpha_composite(logo, (2432 - logo.width, 58))


def _draw_text_centered(
    draw: ImageDraw.ImageDraw,
    text: str,
    center: tuple[int, int],
    f: ImageFont.ImageFont,
    fill: tuple[int, int, int] = INK,
) -> None:
    w, h = tbox(draw, text, f)
    draw.text((center[0] - w // 2, center[1] - h // 2), text, font=f, fill=fill)


def _draw_badge(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    fill: tuple[int, int, int],
    outline: tuple[int, int, int],
    text_fill: tuple[int, int, int] = BLUE,
) -> tuple[int, int, int, int]:
    f = fnt(TIMES_BOLD, 31)
    w, h = tbox(draw, text, f)
    box = (xy[0], xy[1], xy[0] + w + 34, xy[1] + h + 18)
    draw.rounded_rectangle(box, radius=12, fill=(*fill, 255), outline=(*outline, 255), width=2)
    draw.text((xy[0] + 17, xy[1] + 8), text, font=f, fill=text_fill)
    return box


def _draw_large_row_label(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    text: str,
) -> None:
    f = fnt(TIMES_BOLD, 39)
    tw, th = tbox(draw, text, f)
    draw.rounded_rectangle(box, radius=14, fill=(243, 248, 252, 255), outline=(210, 224, 236, 255), width=2)
    draw.text((box[0] + (box[2] - box[0] - tw) // 2, box[1] + (box[3] - box[1] - th) // 2 - 1),
              text, font=f, fill=BLUE)


def _draw_dataset_table_card(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    box = (132, 272, 2428, 1278)
    # Projector-friendly: keep a soft lift, but avoid the heavy laptop-screen
    # card shadow that can muddy on a room display.
    layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(layer, "RGBA")
    sd.rounded_rectangle((box[0] + 5, box[1] + 6, box[2] + 5, box[3] + 6),
                         radius=12, fill=(30, 42, 58, 14))
    base.alpha_composite(layer.filter(ImageFilter.GaussianBlur(9)))
    round_panel(draw, box)
    side_accent(draw, box, PHOTON, width=14)

    x0 = box[0] + 56
    x1 = box[2] - 54
    # Table geometry: one highlighted analysis row plus three forward samples.
    # The outer card spine carries the accent; the table itself stays calm.
    table = (x0, box[1] + 54, x1, box[3] - 56)
    col_sample = table[0] + 44
    col_lumi = table[0] + 820
    col_role = table[0] + 1320
    header_y = table[1]
    draw.rounded_rectangle(table, radius=14, fill=(255, 255, 255, 255), outline=(220, 228, 236, 255), width=2)
    draw.rectangle((table[0], table[1], table[2], table[1] + 94), fill=(244, 248, 251, 255))
    for label, x in (("Dataset", col_sample), ("Integrated luminosity", col_lumi), ("Role in this program", col_role)):
        draw.text((x, header_y + 24), label, font=fnt(TIMES_BOLD, 43), fill=BLUE)
    for vx in (col_lumi - 44, col_role - 44):
        draw.line((vx, table[1] + 14, vx, table[3] - 14), fill=(224, 231, 238, 255), width=2)

    rows = [
        {
            "sample": "Run 24 p+p\n√s = 200 GeV",
            "luminosity": "L = 64.4 pb^{-1}",
            "role": "analyzed sample\nfor this measurement",
            "fill": (255, 255, 255),
            "h": 270,
            "lumi_size": 72,
        },
        {
            "sample": "2026 p+p",
            "luminosity": "17 pb^{-1}",
            "role": "future p+p statistics",
            "fill": (248, 251, 253),
            "h": 178,
            "lumi_size": 50,
        },
        {
            "sample": "2025 Au+Au",
            "luminosity": "6.6 nb^{-1}",
            "role": "move the analysis into A+A",
            "fill": (255, 255, 255),
            "h": 178,
            "lumi_size": 50,
        },
        {
            "sample": "2026 O+O",
            "luminosity": "23.6 nb^{-1}",
            "role": "smaller-system comparison",
            "fill": (250, 250, 252),
            "h": 178,
            "lumi_size": 50,
        },
    ]
    y = table[1] + 94
    for idx, row in enumerate(rows):
        ry0, ry1 = y, y + row["h"]
        draw.rectangle((table[0], ry0, table[2], ry1), fill=(*row["fill"], 255))
        if idx:
            draw.line((table[0] + 28, ry0, table[2] - 28, ry0), fill=(224, 231, 238, 255), width=2)

        sample_lines = str(row["sample"]).split("\n")
        sample_font = fnt(TIMES_BOLD, 50 if idx == 0 else 46)
        y_sample = ry0 + (64 if idx == 0 else 61)
        for line in sample_lines:
            draw.text((col_sample, y_sample), line, font=sample_font, fill=INK)
            y_sample += 58
        rich(draw, (col_lumi, ry0 + (80 if idx == 0 else 60)), str(row["luminosity"]),
             fnt(TIMES_BOLD, int(row["lumi_size"])), RED if idx == 0 else BLUE)
        role_font = fnt(TIMES, 46 if idx == 0 else 43)
        wrapped_vcenter(draw, str(row["role"]), (col_role, ry0 + 22, table[2] - 50, ry1 - 22), role_font, fill=INK, gap=6)
        y = ry1


def _draw_footer(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft = 1326
    draw.rectangle((0, ft, W, H), fill=(*SOFT_BG, 255))
    draw.line((0, ft, W, ft), fill=(*FOOTER_RULE, 255), width=2)
    cy = 1384
    illinois = ASSET_DIR / "illinois_logo_fullcolor_rgb.png"
    if illinois.exists():
        il = fit_img(crop_vis(Image.open(illinois).convert("RGBA"), thr=252), 54, 62)
        base.alpha_composite(il, (30, cy - il.height // 2))
    draw.text((104, cy - 17), "Justin Bennett", font=fnt(TIMES, 31), fill=(43, 49, 57))
    center = "Hard Probes 2026 / June 24, 2026"
    cf = fnt(TIMES_BOLD, 31)
    cw, ch = tbox(draw, center, cf)
    hp_path = ASSET_DIR / "hp2026_indico_logo.png"
    hp = None
    if hp_path.exists():
        hp = fit_img(w2t(Image.open(hp_path).convert("RGBA"), thr=246), 116, 60)
    gw = (hp.width + 20 if hp else 0) + cw
    gx = (W - gw) // 2
    if hp:
        base.alpha_composite(hp, (gx, cy - hp.height // 2))
        gx += hp.width + 20
    draw.text((gx, cy - ch // 2 - 1), center, font=cf, fill=(43, 49, 57))


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    _draw_header(img)
    _draw_dataset_table_card(img)
    _draw_footer(img)
    out = OUTPUT_DIR / "hp2026_slide04_dataset_table_v24_projector_neutral_rows.png"
    img.convert("RGB").save(out, "PNG")
    out.with_suffix(".header.json").write_text(
        json.dumps({"hp2026_main_header": HP2026_MAIN_HEADER}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Saved: {out}")
    return out


if __name__ == "__main__":
    main()
