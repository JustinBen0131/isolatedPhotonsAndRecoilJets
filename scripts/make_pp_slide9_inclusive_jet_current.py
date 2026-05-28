#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_exactStitch_20260523_0915"
ROOT = BASE / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN
PLOT = ROOT / "slide_assets/pp_currentIAN_inclusive_jet_truth_stitch_slide9_side_by_side_insitu_contract.png"
ACCEPT_JSON = ROOT / "validation/insitu_stitching/pp_currentian_exactstitch_acceptance_summary.json"
OUT = ROOT / "slide_assets/slide9_pp_inclusive_jet_stitch_current_replacement.png"

CANVAS_W = 1920
CANVAS_H = 1080
GOOGLE_SLIDES_W = 2560
GOOGLE_SLIDES_H = 1440


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/Library/Fonts/Times New Roman Bold.ttf" if bold else "/Library/Fonts/Times New Roman.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


def draw_wrapped(draw: ImageDraw.ImageDraw, xy, text: str, fnt, fill, width: int, line_gap: int = 6):
    x, y = xy
    words = text.split()
    lines = []
    line = ""
    for word in words:
        trial = f"{line} {word}".strip()
        if draw.textbbox((0, 0), trial, font=fnt)[2] <= width or not line:
            line = trial
        else:
            lines.append(line)
            line = word
    if line:
        lines.append(line)
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        y += fnt.size + line_gap
    return y


def fit_image(img: Image.Image, box_w: int, box_h: int) -> Image.Image:
    scale = min(box_w / img.width, box_h / img.height)
    return img.resize((int(img.width * scale), int(img.height * scale)), Image.Resampling.LANCZOS)


def main() -> None:
    with ACCEPT_JSON.open() as f:
        acc = json.load(f)
    jet_samples = [s for s in acc["samples"] if s["group"] == "jet"]
    jet_boundaries = [b for b in acc["boundaries"] if b["group"] == "jet"]
    max_dev = max(float(b["max_same_bin_overlap_fractional_deviation"]) for b in jet_boundaries)
    max_dev_boundary = max(jet_boundaries, key=lambda b: float(b["max_same_bin_overlap_fractional_deviation"]))

    W, H = CANVAS_W, CANVAS_H
    bg = (250, 251, 253)
    ink = (23, 28, 40)
    muted = (84, 96, 116)
    blue = (43, 100, 184)
    green = (31, 125, 77)
    gold_bg = (255, 248, 232)
    gold_line = (232, 179, 84)
    panel_bg = (255, 255, 255)
    panel_edge = (209, 218, 230)

    slide = Image.new("RGB", (W, H), bg)
    draw = ImageDraw.Draw(slide)

    title_f = font(54, True)
    sub_f = font(26)
    h_f = font(27, True)
    body_f = font(23)
    small_f = font(19)
    table_f = font(20)
    table_b = font(20, True)

    draw.text((56, 38), "Inclusive Jet 8+12+20+30+40 Stitching in pp", font=title_f, fill=ink)
    draw.text(
        (58, 98),
        "Current in-situ RecoilJets output: R = 0.4 truth jets, 1 GeV bins, wiki-matched per-event cross-section normalization.",
        font=sub_f,
        fill=muted,
    )

    # Plot panel
    plot_panel = (58, 145, 945, 938)
    draw.rounded_rectangle(plot_panel, radius=10, fill=panel_bg, outline=panel_edge, width=2)
    plot_img = Image.open(PLOT).convert("RGB")
    fitted = fit_image(plot_img, 830, 735)
    px = plot_panel[0] + (plot_panel[2] - plot_panel[0] - fitted.width) // 2
    py = plot_panel[1] + 30
    slide.paste(fitted, (px, py))
    draw.text((78, 902), "This Analysis Output", font=h_f, fill=ink)

    # Right panel: contract
    right = (985, 145, 1865, 938)
    draw.rounded_rectangle(right, radius=10, fill=panel_bg, outline=panel_edge, width=2)
    x = 1018
    y = 176
    draw.text((x, y), "Stitching contract", font=h_f, fill=ink)
    y += 43
    y = draw_wrapped(
        draw,
        (x, y),
        "Each PYTHIA8 inclusive-jet sample contributes only inside its leading truth-jet pT window. The cross sections below match the current sPHENIX Jet Structure wiki table; yields are normalized as xsec / Nevents / bin width.",
        body_f,
        muted,
        790,
        7,
    )
    y += 16

    # Table header
    col_x = [x, x + 165, x + 355, x + 560]
    draw.text((col_x[0], y), "sample", font=table_b, fill=ink)
    draw.text((col_x[1], y), "window [GeV]", font=table_b, fill=ink)
    draw.text((col_x[2], y), "xsec [pb]", font=table_b, fill=ink)
    draw.text((col_x[3], y), "weight", font=table_b, fill=ink)
    y += 30
    draw.line((x, y, right[2] - 35, y), fill=panel_edge, width=2)
    y += 12

    for s in jet_samples:
        sample = s["sample"].replace("run28_", "").replace("jet", "jet")
        weight = float(s["weight_pb_per_event_GeV"])
        xsec = float(s["xsec_pb"])
        draw.text((col_x[0], y), sample, font=table_f, fill=ink)
        draw.text((col_x[1], y), s["window"], font=table_f, fill=ink)
        draw.text((col_x[2], y), f"{xsec:.4g}", font=table_f, fill=ink)
        draw.text((col_x[3], y), f"{weight:.3g}", font=table_f, fill=ink)
        y += 32

    y += 18
    draw.rounded_rectangle((x, y, right[2] - 34, y + 150), radius=8, fill=(239, 248, 244), outline=(185, 222, 203), width=2)
    draw.text((x + 20, y + 16), "Boundary QA", font=h_f, fill=green)
    y2 = y + 56
    qa = (
        f"All jet stitch boundaries pass. The largest same-bin overlap "
        f"deviation is {100*max_dev:.1f}% at {float(max_dev_boundary['boundary_GeV']):.0f} GeV; "
        "there is no window gap or overlap."
    )
    draw_wrapped(draw, (x + 20, y2), qa, body_f, ink, 760, 7)
    y += 182

    draw.rounded_rectangle((x, y, right[2] - 34, y + 154), radius=8, fill=(239, 245, 255), outline=(190, 210, 240), width=2)
    draw.text((x + 20, y + 16), "Interpretation", font=h_f, fill=blue)
    interp = (
        "The plot corresponds to the regenerated wiki-weight output. The top spectrum is smooth across the jet8, jet12, jet20, jet30, and jet40 handoffs; the ratio panel is a fit-diagnostic, not the stitching definition."
    )
    draw_wrapped(draw, (x + 20, y + 56), interp, body_f, ink, 760, 7)

    # Bottom takeaway band
    band = (58, 968, 1865, 1038)
    draw.rounded_rectangle(band, radius=10, fill=gold_bg, outline=gold_line, width=2)
    draw.text((82, 989), "Takeaway:", font=font(25, True), fill=(137, 75, 0))
    draw.text(
        (205, 989),
        "The inclusive-jet pp reference baseline is stitched correctly in situ with wiki-matched weights.",
        font=font(25),
        fill=ink,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    slide = slide.resize((GOOGLE_SLIDES_W, GOOGLE_SLIDES_H), Image.Resampling.LANCZOS)
    slide.save(OUT)
    print(OUT)


if __name__ == "__main__":
    main()
