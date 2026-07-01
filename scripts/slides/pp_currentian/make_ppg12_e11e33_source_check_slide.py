#!/usr/bin/env python3
"""Build a full-slide PNG explaining the PPG12 E11/E33 source check."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[3]
BASE = REPO / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
FIG_DIR = BASE / "shower_shape_reference_validation/fig13_e11_e33"
OUT_DIR = FIG_DIR / "slide_candidate"

IAN_CROP = (
    REPO
    / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611"
    / "comparison_slide/ian_extract/clean_ref_panels/ppg12_fig13_e11_to_e33_clean.png"
)
SDCC_REPLOT = FIG_DIR / "ppg12_sdcc_root_replot_fig13_e11_to_e33_clean_no_overlap.png"
CURRENT_OVERLAY = FIG_DIR / "ppg12_sdcc_vs_current_fullpp_e11_e33_data_overlay.png"

OUT_PNG = OUT_DIR / "ppg12_e11e33_source_check_slide.png"
OUT_SCRIPT = OUT_DIR / "ppg12_e11e33_source_check_slide_speaker_script.md"
OUT_MANIFEST = OUT_DIR / "ppg12_e11e33_source_check_slide_manifest.json"


W, H = 2560, 1440
MARGIN_X = 110
TITLE_Y = 60


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    candidates: list[Path] = []
    if bold and italic:
        candidates += [
            Path("/Library/Fonts/Times New Roman Bold Italic.ttf"),
            Path("/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf"),
        ]
    elif bold:
        candidates += [
            Path("/Library/Fonts/Times New Roman Bold.ttf"),
            Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"),
        ]
    elif italic:
        candidates += [
            Path("/Library/Fonts/Times New Roman Italic.ttf"),
            Path("/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"),
        ]
    else:
        candidates += [
            Path("/Library/Fonts/Times New Roman.ttf"),
            Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf"),
        ]
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size)
    fallback = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
    if Path(fallback).exists():
        return ImageFont.truetype(fallback, size)
    return ImageFont.load_default()


F_TITLE = font(54, bold=True)
F_SUB = font(29)
F_SECTION = font(30, bold=True)
F_BODY = font(25)
F_SMALL = font(23)
F_TINY = font(19)
F_BOLD = font(27, bold=True)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_width: int,
    fnt: ImageFont.FreeTypeFont,
    *,
    fill: str = "#1F2937",
    line_spacing: int = 8,
) -> int:
    words = text.split()
    lines: list[str] = []
    cur = ""
    for word in words:
        trial = f"{cur} {word}".strip()
        if draw.textbbox((0, 0), trial, font=fnt)[2] <= max_width or not cur:
            cur = trial
        else:
            lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        bbox = draw.textbbox((x, y), line, font=fnt)
        y += (bbox[3] - bbox[1]) + line_spacing
    return y


def rounded_rect(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    fill: str = "white",
    outline: str = "#CBD5E1",
    width: int = 3,
    radius: int = 18,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def paste_fit(
    canvas: Image.Image,
    image_path: Path,
    box: tuple[int, int, int, int],
    *,
    bg: str = "white",
) -> tuple[int, int, int, int]:
    img = Image.open(image_path).convert("RGB")
    x0, y0, x1, y1 = box
    bw, bh = x1 - x0, y1 - y0
    scale = min(bw / img.width, bh / img.height)
    nw, nh = int(img.width * scale), int(img.height * scale)
    panel = Image.new("RGB", (bw, bh), bg)
    panel.paste(img.resize((nw, nh), Image.Resampling.LANCZOS), ((bw - nw) // 2, (bh - nh) // 2))
    canvas.paste(panel, (x0, y0))
    return (x0 + (bw - nw) // 2, y0 + (bh - nh) // 2, x0 + (bw + nw) // 2, y0 + (bh + nh) // 2)


def bullet(draw: ImageDraw.ImageDraw, x: int, y: int, text: str, max_width: int) -> int:
    draw.polygon([(x, y + 9), (x, y + 29), (x + 17, y + 19)], fill="#1F5F9F")
    return draw_wrapped(draw, text, (x + 34, y), max_width - 34, F_BODY, line_spacing=7)


def build() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)

    draw.text((MARGIN_X, TITLE_Y), "E11/E33 pp QA: reference validated, current output now tracks PPG12", font=F_TITLE, fill="#111827")
    draw.text(
        (MARGIN_X, TITLE_Y + 70),
        "First validate the PPG12 SDCC reference against the IAN, then compare current PhotonClusterBuilder output.",
        font=F_SUB,
        fill="#374151",
    )

    left_x0, left_x1 = MARGIN_X, 1070
    right_x0, right_x1 = 1125, W - MARGIN_X
    top_y, bottom_y = 185, 1072
    callout_y0, callout_y1 = 1110, H - 64

    # Left reference card
    rounded_rect(draw, (left_x0, top_y, left_x1, bottom_y), outline="#B8C5D6", width=3, radius=20)
    draw.text((left_x0 + 34, top_y + 24), "Reference check: SDCC reproduces the IAN", font=F_SECTION, fill="#111827")
    draw.text((left_x0 + 34, top_y + 62), "Same PPG12 Fig. 13 E11/E33 source, shown two ways.", font=F_SMALL, fill="#4B5563")

    ian_box = (left_x0 + 40, top_y + 126, left_x1 - 40, top_y + 468)
    sdcc_box = (left_x0 + 40, top_y + 536, left_x1 - 40, bottom_y - 38)
    rounded_rect(draw, (ian_box[0] - 12, ian_box[1] - 34, ian_box[2] + 12, ian_box[3] + 16), outline="#E5E7EB", width=2, radius=12)
    draw.text((ian_box[0], ian_box[1] - 30), "IAN crop", font=F_SMALL, fill="#1F2937")
    paste_fit(canvas, IAN_CROP, ian_box)
    rounded_rect(draw, (sdcc_box[0] - 12, sdcc_box[1] - 34, sdcc_box[2] + 12, sdcc_box[3] + 16), outline="#E5E7EB", width=2, radius=12)
    draw.text((sdcc_box[0], sdcc_box[1] - 30), "Replotted from Shuhang SDCC ROOT", font=F_SMALL, fill="#1F2937")
    paste_fit(canvas, SDCC_REPLOT, sdcc_box)

    # Right current output card
    rounded_rect(draw, (right_x0, top_y, right_x1, bottom_y), outline="#8AB3DE", width=4, radius=20)
    draw.text((right_x0 + 34, top_y + 24), "Current pp output: PhotonClusterBuilder E11/E33", font=F_SECTION, fill="#111827")
    draw.text((right_x0 + 34, top_y + 62), "Compact PhotonClusterBuilder E11/E33 histograms vs the validated PPG12 SDCC data.", font=F_SMALL, fill="#4B5563")
    overlay_box = (right_x0 + 34, top_y + 100, right_x1 - 34, bottom_y - 34)
    paste_fit(canvas, CURRENT_OVERLAY, overlay_box)

    # Bottom explanation
    rounded_rect(draw, (MARGIN_X, callout_y0, W - MARGIN_X, callout_y1), fill="#FFF8DF", outline="#E1B84D", width=4, radius=18)
    draw.text((MARGIN_X + 34, callout_y0 + 28), "Why the earlier bump appeared", font=F_SECTION, fill="#111827")
    y = callout_y0 + 78
    y = bullet(
        draw,
        MARGIN_X + 42,
        y,
        "The bump came from an older THE42 table-QA diagnostic source: different ROOT file, different h2d/h1d histogram family, and different binning/stage provenance.",
        1030,
    )
    y = bullet(
        draw,
        MARGIN_X + 42,
        y + 12,
        "The current full-pp compact E11/E33 histograms are filled from PhotonClusterBuilder shower-shape parameters used by the photon-ID path and track the PPG12 SDCC reference closely.",
        1030,
    )
    draw.line((1255, callout_y0 + 34, 1255, callout_y1 - 34), fill="#E1B84D", width=3)
    draw.text((1290, callout_y0 + 42), "Clean takeaway", font=F_SECTION, fill="#111827")
    draw_wrapped(
        draw,
        "The reference source is validated. The old high-E11/E33 tail was a stale diagnostic-source mismatch, not the current pp shower-shape baseline. The remaining differences are small enough to study as real residuals.",
        (1290, callout_y0 + 92),
        1115,
        F_BODY,
        fill="#1F2937",
        line_spacing=9,
    )

    canvas.save(OUT_PNG)

    script = textwrap.dedent(
        f"""\
        # Speaker Script

        This slide separates two things that had been mixed together.

        On the left, the PPG12 reference is checked first. The top image is the IAN crop, and the lower image is the same Fig. 13 E11/E33 data replotted directly from Shuhang's SDCC ROOT output. That establishes that the reference data source is the right one.

        On the right, I compare that validated SDCC reference against the current full pp RecoilJets compact E11/E33 histogram. This compact histogram is filled from PhotonClusterBuilder shower-shape parameters, so it is the relevant check for the shower-shape information entering our photon-ID path.

        The earlier bump was not a current-pp physics feature. It came from an older THE42 table-QA diagnostic source: a different ROOT file and a different h2d/h1d histogram family with different binning and stage provenance. The current full pp compact source removes that artifact and tracks the PPG12 SDCC reference closely. The remaining differences are small and are now the right object to study.
        """
    )
    OUT_SCRIPT.write_text(script)

    manifest = {
        "png": str(OUT_PNG),
        "speaker_script": str(OUT_SCRIPT),
        "inputs": {
            "ian_crop": str(IAN_CROP),
            "sdcc_replot": str(SDCC_REPLOT),
            "current_overlay": str(CURRENT_OVERLAY),
        },
        "claim": (
            "PPG12 SDCC Fig.13 source is validated against the IAN; current full-pp compact "
            "PhotonClusterBuilder E11/E33 tracks that validated source. Earlier bump came from "
            "an old THE42 table-QA diagnostic source mismatch."
        ),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    print(OUT_PNG)
    print(OUT_SCRIPT)
    print(OUT_MANIFEST)


if __name__ == "__main__":
    build()
