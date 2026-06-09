#!/usr/bin/env python3
"""
Slide 4: Run 24 p+p dataset used for this measurement.

Local PNG generator for HP2026 full-talk slide 4.

This version removes the run-coordinator luminosity plot and the old
"key distinction" bookkeeping card.  The slide now defines only the analysis
sample and luminosity needed for the 15-minute talk narrative.

Output:
  outputs/manual-20260607-run24pp-lumi-asset/
      hp2026_slide04_analysis_dataset_context_v12_story.png
"""
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont

W, H = 2560, 1440

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


# ── slide sections ────────────────────────────────────────────────────────────

def _draw_header(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))
    draw.text((132, 84), "Run 24 p+p dataset used for this measurement",
              font=fnt(TIMES_BOLD, 74), fill=INK)
    rich(draw, (136, 184),
         "Define only the analyzed p+p sample carried into the isolated prompt-photon cross section.",
         fnt(TIMES_ITAL, 36), BLUE)
    draw.line((132, 278, W - 132, 278), fill=(221, 226, 232), width=3)

    logo_path = ASSET_DIR / "sphenix-logo-white-bg_0.png"
    if logo_path.exists():
        logo = fit_img(crop_vis(Image.open(logo_path).convert("RGBA"), thr=252), 316, 92)
        base.alpha_composite(logo, (2444 - logo.width, 58 + (92 - logo.height) // 2))


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


def _draw_sample_definition_card(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    box = (132, 318, 1548, 1136)
    drop_shadow(base, box)
    round_panel(draw, box)
    draw.rounded_rectangle((box[0], box[1], box[0] + 14, box[3]), radius=6, fill=(*PHOTON, 255))

    x = box[0] + 54
    y = box[1] + 42
    draw.text((x, y), "Analysis sample used in this talk", font=fnt(TIMES_BOLD, 50), fill=INK)
    y += 72
    draw.text((x, y), "Run 24 p+p collisions at √s = 200 GeV", font=fnt(TIMES, 39), fill=MUTED)

    # Centerpiece luminosity statement.
    cx = (box[0] + box[2]) // 2
    draw.rounded_rectangle((x, 506, box[2] - 54, 764), radius=18,
                           fill=(255, 250, 237, 255), outline=(*AMBER_BORDER, 255), width=2)
    _draw_text_centered(draw, "PPG12 analysis luminosity", (cx, 556), fnt(TIMES_BOLD, 38), BLUE)
    rich(draw, (x + 314, 610), "L = 64.4 pb^{-1}", fnt(TIMES_BOLD, 118), RED)

    # Three digestible rows: sample, object, use.
    row_y = 828
    rows = [
        ("Sample", "analyzed Run 24 p+p data used for this measurement"),
        ("Object", "isolated prompt photons reconstructed in sPHENIX"),
        ("Use", "normalization for the reported p+p cross section"),
    ]
    for label, body in rows:
        _draw_badge(draw, (x, row_y), label, (243, 248, 252), (210, 224, 236), BLUE)
        draw.text((x + 190, row_y + 7), body, font=fnt(TIMES, 34), fill=INK)
        row_y += 72


def _draw_luminosity_definition_card(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    box = (1600, 318, 2428, 1136)
    drop_shadow(base, box)
    round_panel(draw, box)
    draw.rounded_rectangle((box[0], box[1], box[0] + 14, box[3]), radius=6, fill=(*SPHENIX_BLUE, 255))
    x0 = box[0] + 46
    y = box[1] + 42
    draw.text((x0, y), "How the luminosity is defined", font=fnt(TIMES_BOLD, 48), fill=INK)
    y += 72
    wrapped(draw,
            "Use the calibrated analysis luminosity, not a run-summary luminosity plot.",
            (x0, y), box[2] - x0 - 42, fnt(TIMES_ITAL, 34), fill=MUTED, gap=7)

    # Compact method equation.
    method = (x0, 548, box[2] - 46, 706)
    draw.rounded_rectangle(method, radius=14, fill=(244, 249, 252, 255),
                           outline=(205, 224, 238, 255), width=2)
    rich(draw, (method[0] + 72, method[1] + 48),
         "L_int = N_MB^analyzed / σ_MBD^Vernier",
         fnt(TIMES_BOLD, 45), BLUE)

    rows = [
        ("σ_MBD^Vernier", "MBD minimum-bias trigger cross section measured in a Vernier scan"),
        ("N_MB^analyzed", "minimum-bias-triggered events in the analyzed sample"),
    ]
    y = 756
    for label, body in rows:
        draw.rounded_rectangle((x0, y, box[2] - 46, y + 100), radius=12,
                               fill=(255, 255, 255, 255), outline=(222, 229, 236, 255), width=2)
        rich(draw, (x0 + 28, y + 25), label, fnt(TIMES_BOLD, 34), BLUE)
        wrapped(draw, body, (x0 + 265, y + 19), box[2] - x0 - 315,
                fnt(TIMES, 31), fill=INK, gap=5)
        y += 120

    draw.rounded_rectangle((x0, 1018, box[2] - 46, 1086), radius=12,
                           fill=(255, 250, 237, 255), outline=(*AMBER_BORDER, 255), width=2)
    draw.text((x0 + 28, 1033), "Carried forward as the cross-section normalization.",
              font=fnt(TIMES_BOLD, 31), fill=PHOTON_DARK)


def _draw_future_extensions(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    box = (132, 1162, 2428, 1278)
    draw.rounded_rectangle(box, radius=12, fill=(248, 251, 253, 255),
                           outline=(220, 228, 236, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 12, box[3]), radius=6, fill=(*SPHENIX_BLUE, 255))

    draw.text((box[0] + 42, box[1] + 20), "Future extensions already collected",
              font=fnt(TIMES_BOLD, 34), fill=INK)
    draw.text((box[0] + 42, box[1] + 62), "Not used in this measurement.",
              font=fnt(TIMES, 28), fill=MUTED)

    datasets = [
        ("2026 p+p", "17", "pb^{-1}", PHOTON, (255, 250, 237)),
        ("2025 Au+Au", "6.6", "nb^{-1}", TEAL_SOFT, (241, 250, 251)),
        ("2026 O+O", "23.6", "nb^{-1}", LIGHT_MUTED, (247, 248, 250)),
    ]
    x = box[0] + 710
    y = box[1] + 30
    for sys_label, value, unit, dot_col, fill_col in datasets:
        chip_w = 420
        draw.rounded_rectangle((x, y, x + chip_w, y + 58), radius=22,
                               fill=(*fill_col, 255), outline=(210, 220, 230, 255), width=2)
        draw.rounded_rectangle((x + 18, y + 13, x + 28, y + 45), radius=5, fill=(*dot_col, 255))
        draw.text((x + 48, y + 14), sys_label, font=fnt(TIMES_BOLD, 30), fill=INK)
        rich(draw, (x + 245, y + 14), f"{value} {unit}", fnt(TIMES_BOLD, 28), BLUE)
        x += chip_w + 42


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
    _draw_sample_definition_card(img)
    _draw_luminosity_definition_card(img)
    _draw_future_extensions(img)
    _draw_footer(img)
    out = OUTPUT_DIR / "hp2026_slide04_analysis_dataset_context_v12_story.png"
    img.convert("RGB").save(out, "PNG")
    print(f"Saved: {out}")
    return out


if __name__ == "__main__":
    main()
