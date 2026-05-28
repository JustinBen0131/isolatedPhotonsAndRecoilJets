#!/usr/bin/env python3
"""Plot pp baseline-BDT signal/background score separation."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[1]
IN_CSV = REPO / "dataOutput/ppPhotonMLPipeline/ppg12_fix4_score_overlays/score_separation/pp_baseline_bdt_score_histograms_compact.csv"
OUT_DIR = IN_CSV.parent
OUT_INTEGRATED = OUT_DIR / "pp_baseline_bdt_score_separation_integratedEt.png"
OUT_GRID = OUT_DIR / "pp_baseline_bdt_score_separation_8EtBins_4x2.png"
OUT_GRID_SLIDE31 = OUT_DIR / "pp_baseline_bdt_score_separation_8EtBins_4x2_slide31_body.png"
OUT_SUMMARY = OUT_DIR / "pp_baseline_bdt_score_separation_summary.json"

FONT_REG = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_BOLD_ITALIC = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf"

SIGNAL = (31, 119, 180)
BACKGROUND = (230, 126, 34)
AXIS = (60, 60, 60)
GRID = (230, 230, 230)
TEXT = (25, 25, 25)

ET_ORDER = ["6-8", "8-10", "10-12", "12-15", "15-18", "18-22", "22-28", "28-35"]


def font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size)


def load_rows():
    by_scope: dict[str, dict[str, list[dict[str, float]]]] = defaultdict(lambda: defaultdict(list))
    with IN_CSV.open() as f:
        for row in csv.DictReader(f):
            item = {
                "bin_low": float(row["bin_low"]),
                "bin_high": float(row["bin_high"]),
                "density": float(row["density"]),
                "count": int(row["count"]),
                "entries": int(row["entries"]),
                "auc": float(row["auc"]),
            }
            by_scope[row["scope"]][row["class"]].append(item)
    for cls_map in by_scope.values():
        for rows in cls_map.values():
            rows.sort(key=lambda r: r["bin_low"])
    return by_scope


def step_points(rows, x0, y0, w, h, ymax):
    pts = []
    for r in rows:
        x_left = x0 + r["bin_low"] * w
        x_right = x0 + r["bin_high"] * w
        y = y0 + h - min(r["density"], ymax) / ymax * h
        if not pts:
            pts.append((x_left, y0 + h))
            pts.append((x_left, y))
        else:
            pts.append((x_left, pts[-1][1]))
            pts.append((x_left, y))
        pts.append((x_right, y))
    if rows:
        pts.append((x0 + rows[-1]["bin_high"] * w, y0 + h))
    return [(int(x), int(y)) for x, y in pts]


def draw_rotated_label(img, text, xy, size):
    canvas = Image.new("RGBA", (520, 80), (255, 255, 255, 0))
    d = ImageDraw.Draw(canvas)
    d.text((0, 0), text, font=font(FONT_REG, size), fill=TEXT)
    rot = canvas.rotate(90, expand=True)
    img.paste(rot, xy, rot)


def draw_panel(draw, img, x0, y0, w, h, scope, data, title, show_ylabel=True, show_xlabel=True, ymax=None, small=False):
    tick_font = font(FONT_REG, 22 if small else 28)
    label_font = font(FONT_REG, 26 if small else 34)
    title_font = font(FONT_BOLD, 28 if small else 42)
    note_font = font(FONT_REG, 21 if small else 27)

    if ymax is None:
        ymax = max(r["density"] for cls in ("signal", "background") for r in data[cls]) * 1.08
        ymax = max(1.0, ymax)

    # Frame and grid.
    for i in range(6):
        frac = i / 5
        x = x0 + int(frac * w)
        y = y0 + h - int(frac * h)
        draw.line([(x, y0), (x, y0 + h)], fill=GRID, width=1)
        draw.line([(x0, y), (x0 + w, y)], fill=GRID, width=1)
    draw.rectangle([x0, y0, x0 + w, y0 + h], outline=AXIS, width=2)

    # Histograms.
    sig_pts = step_points(data["signal"], x0, y0, w, h, ymax)
    bkg_pts = step_points(data["background"], x0, y0, w, h, ymax)
    draw.line(bkg_pts, fill=BACKGROUND, width=4 if not small else 3, joint="curve")
    draw.line(sig_pts, fill=SIGNAL, width=4 if not small else 3, joint="curve")

    # Ticks.
    for i in range(6):
        val = i / 5
        x = x0 + int(val * w)
        tick = f"{val:.1f}"
        tw = draw.textlength(tick, font=tick_font)
        draw.text((x - tw / 2, y0 + h + 8), tick, font=tick_font, fill=TEXT)
    y_ticks = [0.0, ymax / 2.0, ymax]
    for val in y_ticks:
        y = y0 + h - int(val / ymax * h)
        lab = f"{val:.0f}" if ymax > 10 else f"{val:.1f}"
        tw = draw.textlength(lab, font=tick_font)
        draw.text((x0 - tw - 10, y - 12), lab, font=tick_font, fill=TEXT)

    tw = draw.textlength(title, font=title_font)
    draw.text((x0 + (w - tw) / 2, y0 - (42 if small else 58)), title, font=title_font, fill=TEXT)
    if show_xlabel:
        xlabel = "BDT score"
        tw = draw.textlength(xlabel, font=label_font)
        draw.text((x0 + (w - tw) / 2, y0 + h + (46 if small else 58)), xlabel, font=label_font, fill=TEXT)
    if show_ylabel:
        draw_rotated_label(img, "Area-normalized density", (x0 - (112 if small else 140), y0 + h // 2 - 250), 28 if small else 34)

    auc = data["signal"][0]["auc"]
    ns = data["signal"][0]["entries"]
    nb = data["background"][0]["entries"]
    box_w, box_h = (225 if small else 300), (75 if small else 92)
    bx, by = x0 + 18, y0 + 18
    draw.rounded_rectangle([bx, by, bx + box_w, by + box_h], radius=6, fill=(255, 255, 255), outline=(210, 210, 210), width=1)
    draw.text((bx + 12, by + 10), f"AUC {auc:.3f}", font=font(FONT_BOLD, 25 if small else 32), fill=TEXT)
    draw.text((bx + 12, by + (40 if small else 52)), f"S {ns:,}   B {nb:,}", font=note_font, fill=(70, 70, 70))


def draw_legend(draw, x, y, font_size=30):
    f = font(FONT_REG, font_size)
    draw.line([(x, y + 14), (x + 70, y + 14)], fill=SIGNAL, width=5)
    draw.text((x + 88, y), "Signal", font=f, fill=TEXT)
    x2 = x + 245
    draw.line([(x2, y + 14), (x2 + 70, y + 14)], fill=BACKGROUND, width=5)
    draw.text((x2 + 88, y), "Background", font=f, fill=TEXT)


def make_integrated(data):
    W, H = 1500, 1050
    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)
    draw.text((80, 44), "pp baseline BDT score separation", font=font(FONT_BOLD, 52), fill=TEXT)
    draw.text((80, 106), "Integrated over photon ET; curves area-normalized separately", font=font(FONT_REG, 31), fill=(70, 70, 70))
    sfont = font(FONT_BOLD_ITALIC, 31)
    draw.text((80, 158), "sPHENIX", font=sfont, fill=TEXT)
    draw.text((80 + int(draw.textlength("sPHENIX", font=sfont)) + 12, 161), "Internal", font=font(FONT_REG, 28), fill=TEXT)
    draw.text((80, 198), "pp baseline BDT: 9-feature PPG12 base-v1E, no isolation input", font=font(FONT_REG, 27), fill=(45, 45, 45))
    draw_panel(draw, img, 170, 295, 1120, 560, "integrated", data["integrated"], "", True, True, ymax=32.0)
    draw_legend(draw, 805, 245, 31)
    img.save(OUT_INTEGRATED)


def make_grid(data):
    W, H = 2300, 1550
    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)
    draw.text((70, 38), "pp baseline BDT score separation by photon ET", font=font(FONT_BOLD, 54), fill=TEXT)
    draw.text((70, 104), "4x2 table of area-normalized signal/background score densities", font=font(FONT_REG, 31), fill=(70, 70, 70))
    sfont = font(FONT_BOLD_ITALIC, 31)
    draw.text((70, 151), "sPHENIX", font=sfont, fill=TEXT)
    draw.text((70 + int(draw.textlength("sPHENIX", font=sfont)) + 12, 154), "Internal", font=font(FONT_REG, 28), fill=TEXT)
    draw.text((70, 190), "pp baseline BDT: 9-feature PPG12 base-v1E, no isolation input", font=font(FONT_REG, 27), fill=(45, 45, 45))
    draw_legend(draw, 1565, 142, 31)

    cols, rows = 4, 2
    panel_w, panel_h = 440, 365
    x_start, y_start = 210, 330
    x_gap, y_gap = 95, 185
    for idx, scope in enumerate(ET_ORDER):
        r = idx // cols
        c = idx % cols
        x = x_start + c * (panel_w + x_gap)
        y = y_start + r * (panel_h + y_gap)
        title = f"{scope} GeV"
        draw_panel(draw, img, x, y, panel_w, panel_h, scope, data[scope], title, show_ylabel=(c == 0), show_xlabel=(r == rows - 1), ymax=32.0, small=True)
    img.save(OUT_GRID)


def make_grid_slide31_body(data):
    W, H = 2700, 1240
    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)

    sfont = font(FONT_BOLD_ITALIC, 32)
    draw.text((70, 34), "sPHENIX", font=sfont, fill=TEXT)
    draw.text((70 + int(draw.textlength("sPHENIX", font=sfont)) + 12, 37), "Internal", font=font(FONT_REG, 29), fill=TEXT)
    draw.text(
        (70, 82),
        "pp baseline BDT: 9-feature PPG12 base-v1E, no isolation input; signal/background normalized separately",
        font=font(FONT_REG, 29),
        fill=(45, 45, 45),
    )
    draw_legend(draw, 1840, 52, 31)
    draw_rotated_label(img, "Area-normalized density", (78, 530), 32)

    cols, rows = 4, 2
    panel_w, panel_h = 500, 390
    x_start, y_start = 220, 245
    x_gap, y_gap = 120, 170
    for idx, scope in enumerate(ET_ORDER):
        r = idx // cols
        c = idx % cols
        x = x_start + c * (panel_w + x_gap)
        y = y_start + r * (panel_h + y_gap)
        title = f"{scope} GeV"
        draw_panel(
            draw,
            img,
            x,
            y,
            panel_w,
            panel_h,
            scope,
            data[scope],
            title,
            show_ylabel=False,
            show_xlabel=(r == rows - 1),
            ymax=32.0,
            small=True,
        )
    img.save(OUT_GRID_SLIDE31)


def main():
    data = load_rows()
    make_integrated(data)
    make_grid(data)
    make_grid_slide31_body(data)
    summary = {
        "input_csv": str(IN_CSV),
        "model": "ppg12_base_v1E_bdt_noIso",
        "feature_note": "9-feature PPG12 base-v1E pp baseline, no isolation inputs; not the 11-feature AuAu baseV3E model.",
        "normalization": "signal and background histograms are each area-normalized separately within each ET scope",
        "et_bins": ET_ORDER,
        "outputs": {
            "integrated": str(OUT_INTEGRATED),
            "et_grid": str(OUT_GRID),
            "et_grid_slide31_body": str(OUT_GRID_SLIDE31),
        },
        "auc_by_scope": {scope: data[scope]["signal"][0]["auc"] for scope in ["integrated"] + ET_ORDER},
        "entries_by_scope": {
            scope: {
                "signal": data[scope]["signal"][0]["entries"],
                "background": data[scope]["background"][0]["entries"],
            }
            for scope in ["integrated"] + ET_ORDER
        },
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(OUT_INTEGRATED)
    print(OUT_GRID)
    print(OUT_GRID_SLIDE31)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
