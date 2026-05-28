#!/usr/bin/env python3
"""Compose PPG12 IAN stitching figures against this-analysis reproductions."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_fullsim_20260521_1811"
BASE = REPO / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN
ASSET_DIR = BASE / "slide_assets"
REF_DIR = BASE / "reference_ppg12_figures"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT_REG = FONT_DIR / "Times New Roman.ttf"
FONT_BOLD = FONT_DIR / "Times New Roman Bold.ttf"

INK = "#0F172A"
MUTED = "#475569"
PANEL = "#F8FAFC"
BORDER = "#CBD5E1"
ORANGE_FILL = "#FFF7ED"
ORANGE_BORDER = "#FED7AA"
ORANGE_TEXT = "#9A3412"


@dataclass(frozen=True)
class SlideSpec:
    title: str
    subtitle: str
    left_label: str
    right_label: str
    box1_title: str
    box1_text: str
    box2_title: str
    box2_text: str
    reference_png: Path
    analysis_png: Path
    out_png: Path


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT_BOLD if bold else FONT_REG), size=size)


def crop_white(img: Image.Image, pad: int = 10) -> Image.Image:
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
        min(max(xs) + pad + 1, w),
        min(max(ys) + pad + 1, h),
    )
    return img.crop(box)


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
    _, _, _, line_h = draw.textbbox((0, 0), "Ag", font=font_obj)
    for line in lines:
        draw.text((x, y), line, font=font_obj, fill=fill)
        y += line_h + line_gap


def make_slide(spec: SlideSpec) -> None:
    ppg12 = crop_white(Image.open(spec.reference_png), pad=10)
    ours = crop_white(Image.open(spec.analysis_png), pad=10)

    slide = Image.new("RGB", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)

    draw.text((72, 42), spec.title, font=font(56, True), fill=INK)
    draw.text((74, 106), spec.subtitle, font=font(28), fill=MUTED)

    draw.rounded_rectangle((78, 146, 1248, 266), radius=18, fill=PANEL, outline=BORDER, width=2)
    draw.text((106, 158), spec.box1_title, font=font(30, True), fill=INK)
    draw_wrapped(
        draw,
        (106, 196),
        spec.box1_text,
        fill="#334155",
        font_obj=font(25),
        max_width=1098,
        line_gap=2,
    )

    draw.rounded_rectangle((1300, 146, 2482, 266), radius=18, fill=ORANGE_FILL, outline=ORANGE_BORDER, width=2)
    draw.text((1328, 158), spec.box2_title, font=font(30, True), fill=ORANGE_TEXT)
    draw_wrapped(
        draw,
        (1328, 196),
        spec.box2_text,
        fill=ORANGE_TEXT,
        font_obj=font(25),
        max_width=1104,
        line_gap=2,
    )

    panel_y = 328
    panel_h = 940
    panel_w = 1160
    gap = 80
    left_x = 78
    right_x = left_x + panel_w + gap

    for x0, label in [(left_x, spec.left_label), (right_x, spec.right_label)]:
        draw.rounded_rectangle(
            (x0, panel_y - 58, x0 + panel_w, panel_y + panel_h + 48),
            radius=20,
            fill=PANEL,
            outline=BORDER,
            width=2,
        )
        draw.text((x0 + 26, panel_y - 46), label, font=font(34, True), fill=INK)

    image_box = (panel_w - 46, panel_h - 16)
    slide.paste(fit_image(ppg12, image_box).convert("RGB"), (left_x + 23, panel_y + 2))
    slide.paste(fit_image(ours, image_box).convert("RGB"), (right_x + 23, panel_y + 2))

    spec.out_png.parent.mkdir(parents=True, exist_ok=True)
    slide.save(spec.out_png, quality=95)
    print(spec.out_png)


def main() -> None:
    specs = [
        SlideSpec(
            title="Photon+jet stitching now uses the exact PPG12 Figure 5 source",
            subtitle="The comparison isolates remaining differences to source path, binning, ownership windows, and ROOT fit behavior.",
            left_label="PPG12 IAN Figure 5",
            right_label="This Analysis Output",
            box1_title="Stitching regions",
            box1_text="photon5: 0-14  |  photon10: 14-22  |  photon20: >=22 GeV leading truth photon",
            box2_title="Weights and fit",
            box2_text="PPG12 input: photon_max_pT_uncut.root, 0.5 GeV bins, cross-section weighted; ROOT modified-power fit over 10-36 GeV.",
            reference_png=Path(
                "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/"
                "NSIRD_screencaptureui_KDLKg5/Screenshot 2026-05-22 at 1.07.00 AM.png"
            ),
            analysis_png=ASSET_DIR / "pp_currentIAN_photon_truth_stitch_ppg12_exact_style.png",
            out_png=ASSET_DIR / "pp_currentIAN_photon_truth_stitch_slide5_side_by_side.png",
        ),
        SlideSpec(
            title="Inclusive-jet stitching now uses the PPG12 Figure 6 histogram path",
            subtitle="The comparison checks the smooth top spectrum and the MC/fit residual structure after matching PPG12 binning.",
            left_label="PPG12 IAN Figure 6",
            right_label="This Analysis Output",
            box1_title="Stitching regions",
            box1_text="jet8: 9-14  |  jet12: 14-21  |  jet20: 21-32  |  jet30: 32-42  |  jet40: >=42 GeV leading truth jet",
            box2_title="Weights and fit",
            box2_text="PPG12 input: MC_efficiency_jet*_bdt_nom.root, h_max_truth_jet_pT, Rebin(10) to 1 GeV; ROOT modified-power fit over 10-50 GeV.",
            reference_png=Path(
                "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/"
                "NSIRD_screencaptureui_G8Wh3l/Screenshot 2026-05-22 at 12.46.59 AM.png"
            ),
            analysis_png=ASSET_DIR / "pp_currentIAN_inclusive_jet_truth_stitch_ppg12_exact_style.png",
            out_png=ASSET_DIR / "pp_currentIAN_inclusive_jet_truth_stitch_slide6_style.png",
        ),
    ]

    for spec in specs:
        make_slide(spec)


if __name__ == "__main__":
    main()
