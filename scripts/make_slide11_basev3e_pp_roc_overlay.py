#!/usr/bin/env python3
"""Make slide-11 ROC overlay with Au+Au baseV3E curves and pp baseline."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[1]
BASE_ROC = REPO / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/slideReady/basev3e_roc/baseBDT_v3E_withCentrality_roc_coarse3_centrality.csv"
PP_ROC = REPO / "dataOutput/ppPhotonMLPipeline/ppg12_fix4_score_overlays/pp_baseline_bdt_roc_points.csv"
OUT_DIR = BASE_ROC.parent
OUT_PNG = OUT_DIR / "baseBDT_v3E_withCentrality_roc_coarse3_plus_ppBaseline.png"
OUT_CSV = OUT_DIR / "baseBDT_v3E_withCentrality_roc_coarse3_plus_ppBaseline.csv"

FONT_REG = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_ITALIC = "/System/Library/Fonts/Supplemental/Arial Italic.ttf"
FONT_BOLD_ITALIC = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf"


def font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size)


def read_base_curves() -> dict[str, list[dict[str, float]]]:
    curves: dict[str, list[dict[str, float]]] = defaultdict(list)
    with BASE_ROC.open() as f:
        for row in csv.DictReader(f):
            curves[row["centrality_bin"]].append(
                {
                    "fpr": float(row["false_positive_rate"]),
                    "tpr": float(row["true_positive_rate"]),
                    "auc": float(row["auc_source"]),
                    "threshold": float(row["threshold"]),
                }
            )
    for rows in curves.values():
        rows.sort(key=lambda r: (r["fpr"], r["tpr"]))
    return dict(curves)


def read_pp_curve() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with PP_ROC.open() as f:
        for row in csv.DictReader(f):
            if row["model"] != "score_ppg12_base_v1E_bdt_noIso":
                continue
            rows.append(
                {
                    "fpr": float(row["fpr"]),
                    "tpr": float(row["tpr"]),
                    "auc": float(row["auc"]),
                    "threshold": float(row["threshold"]) if row["threshold"] != "inf" else float("inf"),
                }
            )
    rows.sort(key=lambda r: (r["fpr"], r["tpr"]))
    # Keep the high-stat pp curve smooth but avoid an oversized slide-side CSV.
    if len(rows) > 1000:
        step = max(1, len(rows) // 1000)
        sampled = rows[::step]
        if sampled[-1] != rows[-1]:
            sampled.append(rows[-1])
        rows = sampled
    return rows


def draw_dashed_line(draw: ImageDraw.ImageDraw, points, fill, width=4, dash=18, gap=12):
    cycle = dash + gap
    phase = 0.0
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        dx = x1 - x0
        dy = y1 - y0
        dist = (dx * dx + dy * dy) ** 0.5
        if dist <= 0:
            continue
        pos = 0.0
        while pos < dist:
            in_dash = phase < dash
            remaining = (dash - phase) if in_dash else (cycle - phase)
            end = min(dist, pos + remaining)
            if in_dash:
                a = pos / dist
                b = end / dist
                draw.line(
                    [(x0 + dx * a, y0 + dy * a), (x0 + dx * b, y0 + dy * b)],
                    fill=fill,
                    width=width,
                )
            phase = (phase + (end - pos)) % cycle
            pos = end


def polyline(draw: ImageDraw.ImageDraw, points, fill, width=5, dashed=False):
    if dashed:
        draw_dashed_line(draw, points, fill, width=width)
    else:
        draw.line(points, fill=fill, width=width, joint="curve")


def main() -> None:
    scale = 2
    width, height = 1936 * scale, 1540 * scale
    left, right, top, bottom = 210 * scale, 135 * scale, 170 * scale, 185 * scale
    plot_w = width - left - right
    plot_h = height - top - bottom

    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)

    title_font = font(FONT_BOLD, 56 * scale)
    label_font = font(FONT_REG, 38 * scale)
    tick_font = font(FONT_REG, 30 * scale)
    small_font = font(FONT_REG, 27 * scale)
    small_bold = font(FONT_BOLD, 27 * scale)
    italic_bold = font(FONT_BOLD_ITALIC, 31 * scale)
    legend_font = font(FONT_REG, 29 * scale)

    def xy(fpr: float, tpr: float) -> tuple[int, int]:
        x = left + int(max(0.0, min(1.0, fpr)) * plot_w)
        y = top + int((1.0 - max(0.0, min(1.0, tpr))) * plot_h)
        return x, y

    # Grid and frame.
    grid = (232, 232, 232)
    axis = (80, 80, 80)
    for i in range(11):
        frac = i / 10
        x = left + int(frac * plot_w)
        y = top + int(frac * plot_h)
        draw.line([(x, top), (x, top + plot_h)], fill=grid, width=1 * scale)
        draw.line([(left, y), (left + plot_w, y)], fill=grid, width=1 * scale)
    draw.rectangle([left, top, left + plot_w, top + plot_h], outline=axis, width=3 * scale)

    # Ticks and labels.
    for i in range(6):
        val = i / 5
        x = left + int(val * plot_w)
        y = top + plot_h - int(val * plot_h)
        tick = f"{val:.1f}"
        draw.line([(x, top + plot_h), (x, top + plot_h + 14 * scale)], fill=axis, width=2 * scale)
        tw = draw.textlength(tick, font=tick_font)
        draw.text((x - tw / 2, top + plot_h + 22 * scale), tick, fill=(25, 25, 25), font=tick_font)
        draw.line([(left - 14 * scale, y), (left, y)], fill=axis, width=2 * scale)
        tw = draw.textlength(tick, font=tick_font)
        draw.text((left - 28 * scale - tw, y - 17 * scale), tick, fill=(25, 25, 25), font=tick_font)

    title = "Base v3E + centrality BDT ROC with pp baseline"
    tw = draw.textlength(title, font=title_font)
    draw.text(((width - tw) / 2, 48 * scale), title, fill=(25, 25, 25), font=title_font)

    xlabel = "Background efficiency"
    tw = draw.textlength(xlabel, font=label_font)
    draw.text((left + (plot_w - tw) / 2, height - 80 * scale), xlabel, fill=(20, 20, 20), font=label_font)

    ylabel = "Signal efficiency"
    y_img = Image.new("RGBA", (520 * scale, 70 * scale), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_img)
    y_draw.text((0, 0), ylabel, fill=(20, 20, 20), font=label_font)
    y_rot = y_img.rotate(90, expand=True)
    img.paste(y_rot, (38 * scale, top + plot_h // 2 - y_rot.height // 2), y_rot)

    # Curves.
    colors = {
        "0-20": (31, 119, 180),
        "20-50": (214, 39, 40),
        "50-80": (44, 160, 44),
    }
    labels = {
        "0-20": "Au+Au 0-20%",
        "20-50": "Au+Au 20-50%",
        "50-80": "Au+Au 50-80%",
    }
    base = read_base_curves()
    pp = read_pp_curve()

    # Random baseline is light dotted so pp dashed remains visually distinct.
    draw_dashed_line(draw, [xy(0, 0), xy(1, 1)], fill=(170, 170, 170), width=2 * scale, dash=8 * scale, gap=10 * scale)

    legend_items = []
    for cent in ["0-20", "20-50", "50-80"]:
        pts = [xy(r["fpr"], r["tpr"]) for r in base[cent]]
        polyline(draw, pts, colors[cent], width=6 * scale)
        legend_items.append((labels[cent], base[cent][0]["auc"], colors[cent], False))

    pp_pts = [xy(r["fpr"], r["tpr"]) for r in pp]
    pp_color = (95, 47, 135)
    polyline(draw, pp_pts, pp_color, width=6 * scale, dashed=True)
    legend_items.append(("pp baseline BDT", pp[0]["auc"], pp_color, True))
    legend_items.append(("Random", None, (170, 170, 170), True))

    # In-plot provenance label.
    lx, ly = left + int(0.58 * plot_w), top + int(0.55 * plot_h)
    draw.text((lx, ly), "sPHENIX", fill=(20, 20, 20), font=italic_bold)
    sx = lx + int(draw.textlength("sPHENIX", font=italic_bold)) + 9 * scale
    draw.text((sx, ly + 3 * scale), "Internal", fill=(20, 20, 20), font=small_font)
    note_lines = [
        "Au+Au embedded validation",
        "15 < photon ET < 35 GeV",
        "Photon12+20 and Jet12+20+30",
        "pp dashed: baseline pp BDT",
    ]
    yy = ly + 40 * scale
    for line in note_lines:
        draw.text((lx, yy), line, fill=(45, 45, 45), font=small_font)
        yy += 34 * scale

    # Legend.
    leg_x, leg_y = left + int(0.62 * plot_w), top + int(0.76 * plot_h)
    leg_w, leg_h = 575 * scale, 205 * scale
    draw.rectangle([leg_x, leg_y, leg_x + leg_w, leg_y + leg_h], fill=(255, 255, 255), outline=(210, 210, 210), width=2 * scale)
    yy = leg_y + 18 * scale
    for label, auc, color, dashed in legend_items:
        x0 = leg_x + 22 * scale
        x1 = leg_x + 108 * scale
        y0 = yy + 18 * scale
        if dashed:
            draw_dashed_line(draw, [(x0, y0), (x1, y0)], fill=color, width=5 * scale, dash=15 * scale, gap=9 * scale)
        else:
            draw.line([(x0, y0), (x1, y0)], fill=color, width=5 * scale)
        text = f"{label}  AUC = {auc:.3f}" if auc is not None else label
        draw.text((leg_x + 125 * scale, yy), text, fill=(25, 25, 25), font=legend_font)
        yy += 36 * scale

    # Write combined CSV provenance.
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["source", "group", "fpr", "tpr", "auc", "threshold"])
        for cent, rows in base.items():
            for r in rows:
                writer.writerow(["auau_basev3e_centrality", cent, r["fpr"], r["tpr"], r["auc"], r["threshold"]])
        for r in pp:
            writer.writerow(["pp_baseline_bdt_noIso", "pp", r["fpr"], r["tpr"], r["auc"], r["threshold"]])

    img = img.resize((width // scale, height // scale), Image.Resampling.LANCZOS)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img.save(OUT_PNG)
    print(OUT_PNG)
    print(OUT_CSV)


if __name__ == "__main__":
    main()
