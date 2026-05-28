#!/usr/bin/env python3
"""Compose a slide-8-style photon+jet stitching comparison placeholder."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OLD_CAMPAIGN = "ppg12_basev3E_currentIAN_fullsim_20260521_1811"
ACTIVE_CAMPAIGN = "ppg12_basev3E_currentIAN_inSituStitch_20260522_0143"

OLD_BASE = REPO / "dataOutput/ppPhotonMLPipeline" / OLD_CAMPAIGN
ACTIVE_BASE = REPO / "dataOutput/ppPhotonMLPipeline" / ACTIVE_CAMPAIGN
OUT_DIR = ACTIVE_BASE / "slide_assets"

PPG12_REF = OLD_BASE / "reference_ppg12_figures/ppg12_analysis_note_combine.png"
CURRENT_PLACEHOLDER = OLD_BASE / "slide_assets/pp_currentIAN_photon_truth_stitch_ppg12_exact_style.png"
OUT = OUT_DIR / "pp_currentIAN_photon_truth_stitch_slide8_comparison_placeholder.png"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT_REG = FONT_DIR / "Times New Roman.ttf"
FONT_BOLD = FONT_DIR / "Times New Roman Bold.ttf"

INK = "#0F172A"
MUTED = "#475569"
PANEL = "#F8FAFC"
BORDER = "#CBD5E1"
YELLOW = "#FEF3C7"
YELLOW_BORDER = "#FCD34D"
YELLOW_TEXT = "#78350F"
BLUE = "#EAF2FF"
BLUE_BORDER = "#BFDBFE"
BLUE_TEXT = "#1E3A8A"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT_BOLD if bold else FONT_REG), size=size)


def crop_nonwhite(img: Image.Image, pad: int = 8) -> Image.Image:
    rgb = img.convert("RGB")
    pix = rgb.load()
    w, h = rgb.size
    xs: list[int] = []
    ys: list[int] = []
    for y in range(h):
        for x in range(w):
            r, g, b = pix[x, y]
            if min(r, g, b) < 248:
                xs.append(x)
                ys.append(y)
    if not xs:
        return img
    return img.crop(
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
    new_size = (round(work.width * scale), round(work.height * scale))
    resized = work.resize(new_size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGBA", size, "white")
    canvas.alpha_composite(resized, ((target_w - new_size[0]) // 2, (target_h - new_size[1]) // 2))
    return canvas


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    *,
    fill: str,
    font_obj: ImageFont.FreeTypeFont,
    max_width: int,
    line_gap: int = 4,
) -> None:
    words = text.split()
    lines: list[str] = []
    line = ""
    for word in words:
        test = f"{line} {word}".strip()
        if draw.textbbox((0, 0), test, font=font_obj)[2] <= max_width or not line:
            line = test
        else:
            lines.append(line)
            line = word
    if line:
        lines.append(line)

    x, y = xy
    line_h = draw.textbbox((0, 0), "Ag", font=font_obj)[3]
    for wrapped in lines:
        draw.text((x, y), wrapped, font=font_obj, fill=fill)
        y += line_h + line_gap


def make_slide() -> None:
    ref = crop_nonwhite(Image.open(PPG12_REF), pad=6)
    current = crop_nonwhite(Image.open(CURRENT_PLACEHOLDER), pad=6)

    slide = Image.new("RGB", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)

    draw.text(
        (72, 42),
        "Photon+jet stitching should reproduce the PPG12 signal-MC spectrum",
        font=font(54, True),
        fill=INK,
    )
    draw.text(
        (74, 104),
        "This is the pp signal-sample cross-section and ownership-window check before trusting the baseV3E BDT comparison.",
        font=font(28),
        fill=MUTED,
    )

    draw.rounded_rectangle((78, 146, 1246, 282), radius=18, fill=PANEL, outline=BORDER, width=2)
    draw.text((106, 160), "What is checked", font=font(31, True), fill=INK)
    draw_wrapped(
        draw,
        (106, 202),
        "photon5, photon10, and photon20 are cross-section weighted and stitched by leading truth photon E_T: 0-14, 14-22, and >=22 GeV.",
        fill="#334155",
        font_obj=font(25),
        max_width=1090,
        line_gap=3,
    )

    draw.rounded_rectangle((1302, 146, 2482, 282), radius=18, fill=YELLOW, outline=YELLOW_BORDER, width=2)
    draw.text((1330, 160), "Current status", font=font(31, True), fill=YELLOW_TEXT)
    draw_wrapped(
        draw,
        (1330, 202),
        "Right panel is a temporary pre-rerun comparison output. It will be replaced by in-situ RecoilJets histograms from the active rerun.",
        fill=YELLOW_TEXT,
        font_obj=font(25),
        max_width=1100,
        line_gap=3,
    )

    panel_y = 354
    panel_h = 884
    panel_w = 1162
    gap = 78
    left_x = 78
    right_x = left_x + panel_w + gap

    for x0, label, chip_text in [
        (left_x, "PPG12 IAN Figure 5", "reference target"),
        (right_x, "Current comparison output", "placeholder to replace"),
    ]:
        draw.rounded_rectangle(
            (x0, panel_y - 64, x0 + panel_w, panel_y + panel_h + 50),
            radius=20,
            fill=PANEL,
            outline=BORDER,
            width=2,
        )
        draw.text((x0 + 28, panel_y - 52), label, font=font(34, True), fill=INK)
        chip_w = draw.textbbox((0, 0), chip_text, font=font(21, True))[2] + 34
        draw.rounded_rectangle(
            (x0 + panel_w - chip_w - 26, panel_y - 51, x0 + panel_w - 26, panel_y - 16),
            radius=10,
            fill=BLUE if x0 == left_x else YELLOW,
            outline=BLUE_BORDER if x0 == left_x else YELLOW_BORDER,
            width=1,
        )
        draw.text(
            (x0 + panel_w - chip_w - 9, panel_y - 47),
            chip_text,
            font=font(21, True),
            fill=BLUE_TEXT if x0 == left_x else YELLOW_TEXT,
        )

    image_box = (panel_w - 54, panel_h - 18)
    slide.paste(fit_image(ref, image_box).convert("RGB"), (left_x + 27, panel_y + 2))
    slide.paste(fit_image(current, image_box).convert("RGB"), (right_x + 27, panel_y + 2))

    footer_y0 = 1310
    draw.rounded_rectangle((78, footer_y0, 2482, 1384), radius=18, fill="#F1F5F9", outline="#E2E8F0", width=1)
    footer_label = "Replacement criterion:"
    footer_font = font(27, True)
    draw.text((110, footer_y0 + 17), footer_label, font=footer_font, fill=INK)
    label_w = draw.textbbox((0, 0), footer_label, font=footer_font)[2]
    draw_wrapped(
        draw,
        (110 + label_w + 22, footer_y0 + 18),
        "the right panel must come from SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_* in the active in-situ run, not from PPG12 output files.",
        fill=MUTED,
        font_obj=font(25),
        max_width=2140,
        line_gap=0,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    slide.save(OUT, quality=95)
    print(OUT)


if __name__ == "__main__":
    make_slide()
