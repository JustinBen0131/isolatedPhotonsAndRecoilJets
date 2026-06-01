#!/usr/bin/env python3
"""Compose the photon+jet stitching side-by-side slide candidate."""

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

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
REFERENCE_ASSETS = BASE / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_fullsim_20260521_1811/reference_ppg12_figures"
CURRENT_ASSETS = (
    BASE
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_exactStitchPhoton0p5_20260526_1307/slide_assets"
)

PPG12_REF = REFERENCE_ASSETS / "ppg12_analysis_note_combine.png"
CURRENT_RHS = CURRENT_ASSETS / "pp_currentIAN_photon_truth_stitch_slide8_side_by_side_insitu_contract.png"
OUT = CURRENT_ASSETS / "pp_exactstitch_slide8_photon_side_by_side_current_rhs.png"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT_REG = FONT_DIR / "Times New Roman.ttf"
FONT_BOLD = FONT_DIR / "Times New Roman Bold.ttf"

INK = "#0F172A"
MUTED = "#475569"
PANEL = "#F8FAFC"
BORDER = "#CBD5E1"
BOTTOM_FILL = "#FFF7ED"
BOTTOM_BORDER = "#FED7AA"
BOTTOM_LABEL = "#B91C1C"
CANVAS = (2560, 1440)


def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT_BOLD if bold else FONT_REG), size=size)


def crop_nonwhite(img: Image.Image, *, pad: int = 12) -> Image.Image:
    rgb = img.convert("RGB")
    pix = rgb.load()
    w, h = rgb.size

    # Some screen-grab references carry black capture bars outside the actual
    # plot. Drop only border-like strips; keep black axes, text, and data.
    for x in range(w):
        black = 0
        for y in range(h):
            r, g, b = pix[x, y]
            if r < 8 and g < 8 and b < 8:
                black += 1
        if black / max(h, 1) > 0.65:
            for y in range(h):
                pix[x, y] = (255, 255, 255)

    for y in range(h):
        black = 0
        for x in range(w):
            r, g, b = pix[x, y]
            if r < 8 and g < 8 and b < 8:
                black += 1
        if black / max(w, 1) > 0.65:
            for x in range(w):
                pix[x, y] = (255, 255, 255)

    xs: list[int] = []
    ys: list[int] = []
    for y in range(h):
        for x in range(w):
            r, g, b = pix[x, y]
            # Keep black axes/text/markers and colored points; drop plain white
            # margins and the black screenshot strip above the PPG12 reference.
            if max(r, g, b) < 245 and not (r < 8 and g < 8 and b < 8 and y < 40):
                xs.append(x)
                ys.append(y)
    if not xs:
        return rgb
    return rgb.crop(
        (
            max(min(xs) - pad, 0),
            max(min(ys) - pad, 0),
            min(max(xs) + pad + 1, w),
            min(max(ys) + pad + 1, h),
        )
    )


def fit_image(img: Image.Image, size: tuple[int, int]) -> Image.Image:
    target_w, target_h = size
    work = img.convert("RGBA")
    scale = min(target_w / work.width, target_h / work.height)
    resized = work.resize((int(work.width * scale), int(work.height * scale)), Image.Resampling.LANCZOS)
    canvas = Image.new("RGBA", size, "white")
    canvas.alpha_composite(resized, ((target_w - resized.width) // 2, (target_h - resized.height) // 2))
    return canvas.convert("RGB")


def main() -> None:
    slide = Image.new("RGB", CANVAS, "white")
    draw = ImageDraw.Draw(slide)

    ref = fit_image(crop_nonwhite(Image.open(PPG12_REF), pad=10), (1120, 940))
    rhs = fit_image(crop_nonwhite(Image.open(CURRENT_RHS), pad=10), (1120, 940))

    draw.text(
        (56, 38),
        "Photon+Jet Stitching Comparison to PPG12",
        font=font(72, bold=True),
        fill=INK,
    )
    draw.text(
        (58, 126),
        "Comparison of stitched photon+jet spectrum to PPG12 Fig. 5.",
        font=font(34),
        fill=MUTED,
    )

    panel_top = 166
    panel_bottom = 1370
    panel_w = 1205
    panel_h = panel_bottom - panel_top
    gutter = 30
    left_x = 58
    right_x = left_x + panel_w + gutter
    image_size = (panel_w - 28, panel_h - 88)

    for x0, title, plot in (
        (left_x, "PPG12 IAN Figure 5", ref),
        (right_x, "This Analysis Output", rhs),
    ):
        draw.rounded_rectangle(
            (x0, panel_top, x0 + panel_w, panel_bottom),
            radius=12,
            fill=PANEL,
            outline=BORDER,
            width=2,
        )
        draw.text((x0 + 30, panel_top + 22), title, font=font(36, bold=True), fill=INK)
        paste_x = x0 + 14 + (image_size[0] - plot.width) // 2
        paste_y = panel_top + 72 + (image_size[1] - plot.height) // 2
        slide.paste(plot, (paste_x, paste_y))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    slide.save(OUT, quality=95)
    print(OUT)


if __name__ == "__main__":
    main()
