#!/usr/bin/env python3
"""Compose PPG12 IAN Fig. 6 and our reproduction side by side."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_fullsim_20260521_1811"
ASSET_DIR = REPO / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN / "slide_assets"
OUT = ASSET_DIR / "pp_currentIAN_inclusive_jet_truth_stitch_slide6_style.png"
BACKUP = ASSET_DIR / "pp_currentIAN_inclusive_jet_truth_stitch_slide6_style_single_plot_backup.png"
OUR = ASSET_DIR / "pp_currentIAN_inclusive_jet_truth_stitch_ppg12_identical_style.png"
PPG12_SCREENSHOT = Path(
    "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/NSIRD_screencaptureui_G8Wh3l/Screenshot 2026-05-22 at 12.46.59 AM.png"
)

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT_REG = FONT_DIR / "Times New Roman.ttf"
FONT_BOLD = FONT_DIR / "Times New Roman Bold.ttf"

INK = "#0F172A"
MUTED = "#475569"
PANEL = "#F8FAFC"
BORDER = "#CBD5E1"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT_BOLD if bold else FONT_REG), size=size)


def crop_white(img: Image.Image, pad: int = 16) -> Image.Image:
    rgb = img.convert("RGB")
    pix = rgb.load()
    w, h = rgb.size
    xs: list[int] = []
    ys: list[int] = []
    for y in range(h):
        for x in range(w):
            r, g, b = pix[x, y]
            if min(r, g, b) < 247:
                xs.append(x)
                ys.append(y)
    if not xs:
        return img
    box = (
        max(min(xs) - pad, 0),
        max(min(ys) - pad, 0),
        min(max(xs) + pad, w),
        min(max(ys) + pad, h),
    )
    return img.crop(box)


def fit_image(img: Image.Image, size: tuple[int, int]) -> Image.Image:
    target_w, target_h = size
    work = img.convert("RGBA")
    scale = min(target_w / work.width, target_h / work.height)
    new_size = (int(work.width * scale), int(work.height * scale))
    resized = work.resize(new_size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGBA", size, "white")
    canvas.alpha_composite(resized, ((target_w - new_size[0]) // 2, (target_h - new_size[1]) // 2))
    return canvas


def main() -> None:
    if OUT.exists() and not BACKUP.exists():
        BACKUP.write_bytes(OUT.read_bytes())

    ppg12 = crop_white(Image.open(PPG12_SCREENSHOT), pad=12)
    ours = crop_white(Image.open(OUR), pad=12)

    slide = Image.new("RGB", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)

    draw.text(
        (74, 44),
        "Inclusive-jet stitching reproduces the PPG12 Figure 6 structure",
        font=font(56, True),
        fill=INK,
    )
    draw.text(
        (76, 106),
        "Same current-IAN jet windows and axis ranges; compare the visible smoothness in the top spectrum and the MC/fit residual pattern.",
        font=font(28),
        fill=MUTED,
    )

    draw.rounded_rectangle((78, 146, 1260, 236), radius=18, fill="#F8FAFC", outline="#CBD5E1", width=2)
    draw.text((106, 158), "Stitching regions", font=font(29, True), fill=INK)
    draw.text(
        (106, 198),
        "jet8: 9-14  |  jet12: 14-21  |  jet20: 21-32  |  jet30: 32-42  |  jet40: >=42 GeV",
        font=font(26),
        fill="#334155",
    )
    draw.rounded_rectangle((1320, 146, 2482, 236), radius=18, fill="#FFF7ED", outline="#FED7AA", width=2)
    draw.text((1348, 158), "Weights and comparison", font=font(29, True), fill="#9A3412")
    draw.text(
        (1348, 198),
        "PPG12 truth-spectrum weights; same counts axis and MC/fit ratio.",
        font=font(26),
        fill="#9A3412",
    )

    panel_y = 300
    panel_h = 985
    panel_w = 1160
    gap = 80
    left_x = 78
    right_x = left_x + panel_w + gap

    for x0, label in [(left_x, "PPG12 IAN Figure 6"), (right_x, "This Analysis Output")]:
        draw.rounded_rectangle(
            (x0, panel_y - 58, x0 + panel_w, panel_y + panel_h + 48),
            radius=20,
            fill=PANEL,
            outline=BORDER,
            width=2,
        )
        draw.text((x0 + 26, panel_y - 46), label, font=font(34, True), fill=INK)

    image_box = (panel_w - 46, panel_h - 18)
    ppg12_fit = fit_image(ppg12, image_box)
    ours_fit = fit_image(ours, image_box)
    slide.paste(ppg12_fit.convert("RGB"), (left_x + 23, panel_y + 2))
    slide.paste(ours_fit.convert("RGB"), (right_x + 23, panel_y + 2))

    slide.save(OUT, quality=95)
    print(OUT)
    print(BACKUP)


if __name__ == "__main__":
    main()
