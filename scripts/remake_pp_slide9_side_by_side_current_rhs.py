#!/usr/bin/env python3
"""Replace the RHS plot on the existing slide-9 side-by-side PNG."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_exactStitch_20260523_0915"
ASSETS = BASE / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN / "slide_assets"

TEMPLATE = ASSETS / "pp_exactstitch_slide9_jet_side_by_side_decision.png"
CURRENT_RHS = ASSETS / "pp_currentIAN_inclusive_jet_truth_stitch_slide9_side_by_side_insitu_contract.png"
OUT = ASSETS / "pp_exactstitch_slide9_jet_side_by_side_current_rhs.png"

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


def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT_BOLD if bold else FONT_REG), size=size)


def crop_white(img: Image.Image, pad: int = 12) -> Image.Image:
    rgb = img.convert("RGB")
    pix = rgb.load()
    w, h = rgb.size
    xs: list[int] = []
    ys: list[int] = []
    for y in range(h):
        for x in range(w):
            if min(pix[x, y]) < 247:
                xs.append(x)
                ys.append(y)
    if not xs:
        return img
    return img.crop(
        (
            max(min(xs) - pad, 0),
            max(min(ys) - pad, 0),
            min(max(xs) + pad, w),
            min(max(ys) + pad, h),
        )
    )


def fit_image(img: Image.Image, size: tuple[int, int], *, trim_white: bool = True) -> Image.Image:
    target_w, target_h = size
    work = (crop_white(img) if trim_white else img).convert("RGBA")
    scale = min(target_w / work.width, target_h / work.height)
    resized = work.resize((int(work.width * scale), int(work.height * scale)), Image.Resampling.LANCZOS)
    canvas = Image.new("RGBA", size, "white")
    canvas.alpha_composite(resized, ((target_w - resized.width) // 2, (target_h - resized.height) // 2))
    return canvas.convert("RGB")


def main() -> None:
    template = Image.open(TEMPLATE).convert("RGB")
    ppg12_ref = template.crop((101, 206, 1215, 1220))
    slide = Image.new("RGB", (2560, 1440), "white")

    draw = ImageDraw.Draw(slide)
    panel_w = 1160
    panel_top = 142
    panel_bottom = 1330
    panel_h = panel_bottom - panel_top
    left_x = 78
    right_x = left_x + panel_w + 80
    image_y = 203
    image_box = (1114, panel_h - 100)
    ref = fit_image(ppg12_ref, image_box, trim_white=False)
    rhs = fit_image(Image.open(CURRENT_RHS), image_box, trim_white=False)

    draw.text(
        (74, 28),
        "Inclusive-Jet Stitching Comparison to PPG12",
        font=font(66, bold=True),
        fill=INK,
    )
    draw.text(
        (76, 104),
        "RHS uses the same wiki-matched in-situ inclusive-jet plot shown on the previous slide.",
        font=font(28),
        fill=MUTED,
    )

    for x0, label, panel_img in (
        (left_x, "PPG12 IAN Figure 6", ref),
        (right_x, "This Analysis Output", rhs),
    ):
        draw.rounded_rectangle(
            (x0, panel_top, x0 + panel_w, panel_bottom),
            radius=20,
            fill=PANEL,
            outline=BORDER,
            width=2,
        )
        draw.text((x0 + 26, panel_top + 18), label, font=font(34, bold=True), fill=INK)
        slide.paste(panel_img, (x0 + 23, image_y))

    draw.rounded_rectangle(
        (70, 1348, 2490, 1420),
        radius=14,
        fill=BOTTOM_FILL,
        outline=BOTTOM_BORDER,
        width=2,
    )
    draw.text(
        (98, 1366),
        "Note: PPG12 used jet8 = 1.15e7 pb (wiki: 1.3013e7 pb); comparison is in backup when using matched PPG12 weight instead of wiki.",
        font=font(30),
        fill="#9A3412",
    )
    slide.save(OUT, quality=95)
    print(OUT)


if __name__ == "__main__":
    main()
