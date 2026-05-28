#!/usr/bin/env python3
"""Compose the pp exact-stitch photon and jet spectrum QA plots into one slide."""

from __future__ import annotations

from pathlib import Path
import textwrap

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppPhotonMLPipeline/"
    "ppg12_basev3E_currentIAN_exactStitch_20260523_0915"
)
ASSET_DIR = ROOT / "slide_assets"
PHOTON = ASSET_DIR / "pp_currentIAN_photon_truth_stitch_slide8_side_by_side_insitu_contract.png"
JET = ASSET_DIR / "pp_currentIAN_inclusive_jet_truth_stitch_slide9_side_by_side_insitu_contract.png"
OUT = ASSET_DIR / "pp_exactstitch_slide8_9_spectrum_proof_1x2.png"

W, H = 2400, 1350


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/System/Library/Fonts/Times.ttc",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def contain(img: Image.Image, box_w: int, box_h: int) -> Image.Image:
    scale = min(box_w / img.width, box_h / img.height)
    new_size = (int(img.width * scale), int(img.height * scale))
    return img.resize(new_size, Image.Resampling.LANCZOS)


def main() -> None:
    photon = Image.open(PHOTON).convert("RGB")
    jet = Image.open(JET).convert("RGB")

    canvas = Image.new("RGB", (W, H), "#fbfbf8")
    draw = ImageDraw.Draw(canvas)

    f_title = font(48, bold=True)
    f_sub = font(25)
    f_panel = font(28, bold=True)
    f_note = font(22)

    draw.text((90, 62), "pp in-situ stitching: spectrum + fit proof", font=f_title, fill="#111827")
    draw.text(
        (90, 125),
        "This Analysis Output from Justin/RecoilJets exact-stitch histograms; PPG12/Shuhang files are not used as the plotted spectra",
        font=f_sub,
        fill="#475569",
    )

    left_box = (85, 205, 1135, 1190)
    right_box = (1265, 205, 2315, 1190)
    for box in (left_box, right_box):
        draw.rounded_rectangle(box, radius=24, fill="white", outline="#cfd6e3", width=2)

    draw.text((125, 230), "Photon+jet: leading truth photon E_T", font=f_panel, fill="#111827")
    draw.text((1305, 230), "Inclusive jet: leading R=0.4 truth jet p_T", font=f_panel, fill="#111827")

    p_img = contain(photon, 930, 845)
    j_img = contain(jet, 930, 845)
    canvas.paste(p_img, (left_box[0] + (left_box[2] - left_box[0] - p_img.width) // 2, 300))
    canvas.paste(j_img, (right_box[0] + (right_box[2] - right_box[0] - j_img.width) // 2, 300))

    note = (
        "Acceptance logic: spectra are weighted by xsec / Nevt / bin width, stitched in non-overlapping generator windows, "
        "then checked against a smooth modified-power-law reference. The ratio panels are the visual merge test; the boundary table is the numerical audit."
    )
    draw.rounded_rectangle((90, 1220, 2310, 1300), radius=20, fill="#f8fafc", outline="#cfd6e3", width=2)
    draw.multiline_text((125, 1234), "\n".join(textwrap.wrap(note, width=145)), font=font(19), fill="#334155", spacing=5)

    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    canvas.save(OUT)
    print(OUT)


if __name__ == "__main__":
    main()
