#!/usr/bin/env python3
"""Compose a 16:9 before/after slide for the pp jet8 xsec update."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


W, H = 2560, 1440


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    paths = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/System/Library/Fonts/Times.ttc",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
    ]
    for path in paths:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


def fit(img: Image.Image, box_w: int, box_h: int) -> Image.Image:
    scale = min(box_w / img.width, box_h / img.height)
    return img.resize((int(img.width * scale), int(img.height * scale)), Image.Resampling.LANCZOS)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old", required=True, type=Path)
    parser.add_argument("--new", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    old = Image.open(args.old).convert("RGB")
    new = Image.open(args.new).convert("RGB")

    canvas = Image.new("RGB", (W, H), "#fbfbf8")
    draw = ImageDraw.Draw(canvas)
    title = font(66, True)
    sub = font(33)
    label = font(34, True)
    body = font(29)

    draw.text((80, 45), "Inclusive-Jet Stitching: Jet8 Weight Check", font=title, fill="#111827")
    draw.text(
        (82, 122),
        "All other pp PYTHIA weights match the wiki. This test changes only jet8 and keeps the same in-situ RecoilJets histograms.",
        font=sub,
        fill="#475569",
    )

    note = (80, 178, 2480, 250)
    draw.rounded_rectangle(note, radius=16, fill="#eef7f0", outline="#95cfa5", width=3)
    draw.text((112, 194), "Checked difference:", font=font(30, True), fill="#166534")
    draw.text(
        (405, 194),
        "jet8 old local value = 1.15e7 pb; wiki value = 1.3013e7 pb. RHS uses the wiki value.",
        font=font(30),
        fill="#111827",
    )

    boxes = [(80, 282, 1238, 1252), (1322, 282, 2480, 1252)]
    for box in boxes:
        draw.rounded_rectangle(box, radius=18, fill="white", outline="#cfd6e3", width=3)

    draw.text((122, 305), "Before: PPG12 repo jet8 xsec", font=label, fill="#111827")
    draw.text((1364, 305), "After: wiki jet8 xsec", font=label, fill="#111827")
    draw.text((122, 350), "jet8 = 1.15e7 pb, matching PPG12 CrossSectionWeights.h", font=body, fill="#6b7280")
    draw.text((1364, 350), "jet8 = 1.3013e7 pb", font=body, fill="#166534")

    old_fit = fit(old, 1060, 830)
    new_fit = fit(new, 1060, 830)
    canvas.paste(old_fit, (boxes[0][0] + (boxes[0][2] - boxes[0][0] - old_fit.width) // 2, 400))
    canvas.paste(new_fit, (boxes[1][0] + (boxes[1][2] - boxes[1][0] - new_fit.width) // 2, 400))

    band = (80, 1290, 2480, 1382)
    draw.rounded_rectangle(band, radius=16, fill="#fff8e8", outline="#e8b354", width=3)
    draw.text((118, 1312), "Takeaway:", font=font(31, True), fill="#8a4b00")
    draw.text(
        (285, 1312),
        "The wiki-matched jet8 weight improves the low-pT handoff test without rerunning Condor.",
        font=font(31),
        fill="#111827",
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(args.out)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
