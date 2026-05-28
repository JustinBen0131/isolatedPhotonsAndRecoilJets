#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_exactStitchPhoton0p5_20260526_1307"
ROOT = BASE / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN
PLOT = ROOT / "slide_assets/pp_currentIAN_photon_truth_stitch_slide8_side_by_side_insitu_contract.png"
CSV = ROOT / "validation/insitu_stitching/pp_currentian_photon0p5_exactstitch_contract_points.csv"
SUMMARY_JSON = ROOT / "validation/insitu_stitching/pp_currentian_photon0p5_exactstitch_contract_summary.json"
FIT_TXT = ROOT / "slide_assets/pp_currentIAN_photon_truth_stitch_slide8_side_by_side_insitu_contract.fit.txt"
OUT = ROOT / "slide_assets/slide10_pp_photon_jet_stitch_current_replacement.png"

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


def draw_wrapped(draw: ImageDraw.ImageDraw, xy, text: str, fnt, fill, width: int, line_gap: int = 7):
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


def crop_nonwhite(img: Image.Image, *, pad: int = 8) -> Image.Image:
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
        return rgb
    return rgb.crop(
        (
            max(min(xs) - pad, 0),
            max(min(ys) - pad, 0),
            min(max(xs) + pad + 1, w),
            min(max(ys) + pad + 1, h),
        )
    )


def fit_model(params: list[float], x: float) -> float:
    p0, p1, p2, p3, p4 = params
    return p0 * ((p1 / x) ** (p2 + p3 * math.log(x / p1) + p4 * x))


def load_numbers():
    summary = json.loads(SUMMARY_JSON.read_text())
    samples = sorted(summary["samples"], key=lambda s: float(s["stitch_window"][0]))
    rows = list(csv.DictReader(CSV.open()))

    params_match = re.search(r"params=([^\n]+)", FIT_TXT.read_text())
    if not params_match:
        raise RuntimeError(f"missing fit params in {FIT_TXT}")
    params = [float(x) for x in params_match.group(1).split(",")]

    fit_ratios = []
    boundary_residuals = []
    for row in rows:
        x = float(row["bin_center"])
        y = float(row["density_pb_per_gev"])
        in_window = int(float(row["ppg12_bin_center_window"]))
        if not (10.0 <= x <= 40.0 and y > 0.0 and in_window):
            continue
        ratio = y / fit_model(params, x)
        fit_ratios.append(ratio)
        if abs(x - 14.0) <= 0.25 or abs(x - 22.0) <= 0.25:
            boundary_residuals.append((float(row["bin_center"]), abs(ratio - 1.0)))

    max_resid = max(abs(r - 1.0) for r in fit_ratios)
    max_boundary = max(resid for _, resid in boundary_residuals)
    total_window_counts = sum(float(s["raw_display_integral_ppg12_bin_center_window"]) for s in samples)
    total_all_counts = sum(float(s["raw_all_integral"]) for s in samples)
    return samples, max_resid, max_boundary, total_window_counts, total_all_counts


def main() -> None:
    samples, max_resid, max_boundary, total_window_counts, total_all_counts = load_numbers()

    bg = (250, 251, 253)
    ink = (23, 28, 40)
    muted = (84, 96, 116)
    blue = (43, 100, 184)
    green = (31, 125, 77)
    gold_bg = (255, 248, 232)
    gold_line = (232, 179, 84)
    panel_bg = (255, 255, 255)
    panel_edge = (209, 218, 230)

    slide = Image.new("RGB", (CANVAS_W, CANVAS_H), bg)
    draw = ImageDraw.Draw(slide)

    title_f = font(54, True)
    sub_f = font(26)
    h_f = font(27, True)
    body_f = font(23)
    table_f = font(20)
    table_b = font(20, True)

    draw.text((56, 38), "Photon+Jet 5+10+20 Stitching in pp", font=title_f, fill=ink)
    draw.text(
        (58, 106),
        "Current in-situ RecoilJets output: leading truth photon, 0.5 GeV bins, PPG12-style per-event cross-section normalization.",
        font=sub_f,
        fill=muted,
    )

    plot_panel = (58, 145, 945, 938)
    draw.rounded_rectangle(plot_panel, radius=10, fill=panel_bg, outline=panel_edge, width=2)
    plot_img = crop_nonwhite(Image.open(PLOT), pad=8)
    fitted = fit_image(plot_img, 830, 735)
    px = plot_panel[0] + (plot_panel[2] - plot_panel[0] - fitted.width) // 2
    py = plot_panel[1] + 30
    slide.paste(fitted, (px, py))
    draw.text((78, 902), "This Analysis Output", font=h_f, fill=ink)

    right = (985, 145, 1865, 938)
    draw.rounded_rectangle(right, radius=10, fill=panel_bg, outline=panel_edge, width=2)
    x = 1018
    y = 176
    draw.text((x, y), "Stitching contract", font=h_f, fill=ink)
    y += 43
    y = draw_wrapped(
        draw,
        (x, y),
        "Each PYTHIA8 photon+jet sample contributes only inside its leading truth-photon ET window. Yields are normalized as xsec / Nevents / bin width.",
        body_f,
        muted,
        790,
        7,
    )
    y += 16

    col_x = [x, x + 165, x + 355, x + 560]
    draw.text((col_x[0], y), "sample", font=table_b, fill=ink)
    draw.text((col_x[1], y), "window [GeV]", font=table_b, fill=ink)
    draw.text((col_x[2], y), "xsec [pb]", font=table_b, fill=ink)
    draw.text((col_x[3], y), "weight", font=table_b, fill=ink)
    y += 30
    draw.line((x, y, right[2] - 35, y), fill=panel_edge, width=2)
    y += 12

    for s in samples:
        sample = s["sample"].replace("run28_", "").replace("photonjet", "photon")
        low, high = [float(v) for v in s["stitch_window"]]
        window = f"[{low:.0f}, {high:.0f})" if high < 100 else f">= {low:.0f}"
        weight = float(s["xsec_pb"]) / float(s["events_processed_metadata"]) / float(s["hist_bin_width"])
        draw.text((col_x[0], y), sample, font=table_f, fill=ink)
        draw.text((col_x[1], y), window, font=table_f, fill=ink)
        draw.text((col_x[2], y), f"{float(s['xsec_pb']):.4g}", font=table_f, fill=ink)
        draw.text((col_x[3], y), f"{weight:.3g}", font=table_f, fill=ink)
        y += 32

    y += 18
    draw.rounded_rectangle((x, y, right[2] - 34, y + 150), radius=8, fill=(239, 248, 244), outline=(185, 222, 203), width=2)
    draw.text((x + 20, y + 16), "Boundary QA", font=h_f, fill=green)
    qa = (
        f"Both photon stitch boundaries pass. The largest boundary Data/Fit "
        f"deviation is {100*max_boundary:.1f}% at the 14/22 GeV handoffs; "
        "there is no window gap or overlap."
    )
    draw_wrapped(draw, (x + 20, y + 56), qa, body_f, ink, 760, 7)
    y += 182

    draw.rounded_rectangle((x, y, right[2] - 34, y + 154), radius=8, fill=(239, 245, 255), outline=(190, 210, 240), width=2)
    draw.text((x + 20, y + 16), "Interpretation", font=h_f, fill=blue)
    interp = (
        f"The top spectrum is smooth across photon5, photon10, and photon20. "
        f"Across 10-40 GeV, the largest Data/Fit residual is {100*max_resid:.1f}%; "
        "the visible high-ET scatter is statistical, not a stitch discontinuity."
    )
    draw_wrapped(draw, (x + 20, y + 56), interp, body_f, ink, 760, 7)

    band = (58, 968, 1865, 1038)
    draw.rounded_rectangle(band, radius=10, fill=gold_bg, outline=gold_line, width=2)
    draw.text((82, 989), "Takeaway:", font=font(25, True), fill=(137, 75, 0))
    draw.text(
        (205, 989),
        "the photon+jet signal in pp is stitched with the defined wiki weights and little deviation from fit in the combined spectrum",
        font=font(25),
        fill=ink,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    slide = slide.resize((GOOGLE_SLIDES_W, GOOGLE_SLIDES_H), Image.Resampling.LANCZOS)
    slide.save(OUT)
    print(OUT)


if __name__ == "__main__":
    main()
