#!/usr/bin/env python3
"""Build the pp current-IAN E_T/eta reweighting QA slide in the slide-9 style."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw

import make_ppg12_exact_reweighting_verification_slide as base


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "ppg12_basev3E_currentIAN_fullsim_20260521_1811"
INPUT_DIR = (
    REPO
    / "dataOutput/ppPhotonMLPipeline"
    / CAMPAIGN
    / "validation/ppg12_exact_reweight_closure"
)
OUT_DIR = REPO / "dataOutput/ppPhotonMLPipeline" / CAMPAIGN / "slide_assets"
OUT = OUT_DIR / "pp_currentIAN_et_eta_reweighting_verification_slide9_exact.png"


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    title_font = base.font(70, True)
    sub_font = base.font(43, True)
    x, y = 78, 48
    draw.text((x, y), "E", font=title_font, fill=base.INK)
    x += base.text_width(draw, "E", title_font) + 2
    draw.text((x, y + 39), "T", font=sub_font, fill=base.INK)
    x += base.text_width(draw, "T", sub_font) + 4
    draw.text((x, y), "/eta-reweighting in pp", font=title_font, fill=base.INK)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    base.BASE = INPUT_DIR

    meta = json.loads((INPUT_DIR / "ppg12_exact_reweighting_metadata.json").read_text())
    weighting = meta["weighting"]
    et_data = base.load_binned("cluster_Et")
    eta_data = base.load_binned("cluster_Eta")

    sum0 = float(weighting["sum_weight_class0"])
    sum1 = float(weighting["sum_weight_class1"])
    imbalance = abs(sum1 - sum0) / ((sum1 + sum0) / 2.0) * 100.0
    et_ratio = [s["weighted"] / b["weighted"] for s, b in zip(et_data[1], et_data[0]) if b["weighted"] > 0]
    eta_ratio = [s["weighted"] / b["weighted"] for s, b in zip(eta_data[1], eta_data[0]) if b["weighted"] > 0]
    et_dev = max(abs(r - 1.0) for r in et_ratio) * 100.0
    eta_dev = max(abs(r - 1.0) for r in eta_ratio) * 100.0

    slide = Image.new("RGBA", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)

    draw_title(draw)
    base.draw_math_text(
        draw,
        (80, 134),
        "PPG12-style training weights are applied before pp baseV3E training; physics and stitching weights stay out of the training.",
        base.font(32),
        fill=base.MUTED,
    )

    base.draw_chip(draw, (80, 214, 575, 306), "1", "equal total signal and background weight")
    base.draw_chip(draw, (600, 214, 1095, 306), "2", "flatten eta separately for each truth class")
    base.draw_chip(draw, (1120, 214, 1615, 306), "3", "flatten E_T separately for each truth class")
    draw.rounded_rectangle((1660, 206, 2480, 316), radius=28, fill="#FFF7ED", outline="#FED7AA", width=2)
    draw.text((1692, 228), "Do not mix training weights with physics weights", font=base.font(31, True), fill="#9A3412")
    draw.text((1692, 274), "No event, cross-section, stitching, vertex, or centrality weights.", font=base.font(26), fill="#9A3412")

    base.draw_density_plot(
        draw,
        (80, 360, 1215, 760),
        et_data,
        "raw",
        "Before weighting: E_T populations are not comparable",
        (5, 35),
        "cluster E_T [GeV]",
        "density",
        (0.0, 0.22),
        legend_position="mid-right",
    )
    base.draw_ratio_plot(
        draw,
        (1280, 360, 2480, 760),
        et_data,
        "After weighting: E_T signal/background ratio",
        (5, 35),
        "cluster E_T [GeV]",
        et_dev,
    )
    base.draw_density_plot(
        draw,
        (80, 820, 1215, 1190),
        eta_data,
        "raw",
        "Before weighting: eta populations are also adjusted",
        (-0.7, 0.7),
        "cluster eta",
        "density",
        (0.0, 0.92),
        legend_position="bottom-right",
    )
    base.draw_ratio_plot(
        draw,
        (1280, 820, 2480, 1190),
        eta_data,
        "After weighting: eta signal/background ratio",
        (-0.7, 0.7),
        "cluster eta",
        eta_dev,
    )

    base.draw_metric(draw, (80, 1238, 560, 1380), f"{imbalance:.2f}%", "residual total S/B weight imbalance", base.GREEN)
    base.draw_metric(draw, (590, 1238, 1070, 1380), f"{et_dev:.1f}%", "largest E_T closure deviation over 5-35 GeV", base.PURPLE)
    base.draw_metric(draw, (1100, 1238, 1580, 1380), f"{eta_dev:.1f}%", "largest eta closure deviation over |eta| < 0.7", base.PURPLE)

    draw.rounded_rectangle((1610, 1238, 2480, 1380), radius=24, fill="#ECFDF5", outline="#A7F3D0", width=2)
    draw.text((1644, 1264), "Main readout", font=base.font(31, True), fill="#065F46")
    base.draw_wrapped(
        draw,
        "After reweighting, the pp BDT comparison is not dominated by the raw E_T/eta population difference.",
        (1644, 1312),
        780,
        base.font(25),
        fill="#064E3B",
        line_gap=3,
    )

    slide.convert("RGB").save(OUT, quality=95)
    print(OUT)


if __name__ == "__main__":
    main()
