#!/usr/bin/env python3
"""Build a slide-sized PNG summarizing PPG12-exact E_T/eta reweighting QA."""

from __future__ import annotations

import csv
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
BASE = REPO / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/slideReady/ppg12_exact_reweight_bdt_remote"
OUT_DIR = REPO / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/slideReady/ppg12_exact_reweight_bdt"
OUT = OUT_DIR / "ppg12_exact_et_eta_reweighting_verification_slide.png"
SUBTITLE = "PPG12-style training weights are applied before routed E_T training; physics and stitching weights stay out of the training."
CONTEXT_LABEL = ""
ET_XMIN = 5.0
ET_XMAX = 35.0

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT_REG = FONT_DIR / "Times New Roman.ttf"
FONT_BOLD = FONT_DIR / "Times New Roman Bold.ttf"

INK = "#0F172A"
MUTED = "#475569"
GRID = "#E5E7EB"
BLUE = "#2563EB"
ORANGE = "#EA580C"
GREEN = "#059669"
PURPLE = "#7C3AED"
SIGNAL_RED = "#DC2626"
BACKGROUND_BLUE = BLUE


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT_BOLD if bold else FONT_REG), size=size)


def text_width(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont) -> int:
    if "E_T" in text:
        width = 0
        remaining = text
        while "E_T" in remaining:
            before, remaining = remaining.split("E_T", 1)
            width += draw.textbbox((0, 0), before, font=fnt)[2]
            sub_fnt = font(max(12, int(fnt.size * 0.62)), fnt.path == str(FONT_BOLD) if hasattr(fnt, "path") else False)
            width += draw.textbbox((0, 0), "E", font=fnt)[2] + draw.textbbox((0, 0), "T", font=sub_fnt)[2] + 4
        width += draw.textbbox((0, 0), remaining, font=fnt)[2]
        return int(width)
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0]


def draw_math_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    text: str,
    fnt: ImageFont.FreeTypeFont,
    fill: str = INK,
) -> None:
    x, y = xy
    if "E_T" not in text:
        draw.text((x, y), text, font=fnt, fill=fill)
        return
    remaining = text
    while "E_T" in remaining:
        before, remaining = remaining.split("E_T", 1)
        if before:
            draw.text((x, y), before, font=fnt, fill=fill)
            x += text_width(draw, before, fnt)
        sub_fnt = font(max(12, int(fnt.size * 0.62)))
        draw.text((x, y), "E", font=fnt, fill=fill)
        x += text_width(draw, "E", fnt) + 1
        draw.text((x, y + fnt.size * 0.48), "T", font=sub_fnt, fill=fill)
        x += text_width(draw, "T", sub_fnt) + 4
    if remaining:
        draw.text((x, y), remaining, font=fnt, fill=fill)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_width: int,
    fnt: ImageFont.FreeTypeFont,
    fill: str = INK,
    line_gap: int = 8,
) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = (current + " " + word).strip()
        if text_width(draw, trial, fnt) <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    x, y = xy
    for line in lines:
        draw_math_text(draw, (x, y), line, fnt, fill=fill)
        y += fnt.size + line_gap
    return y


def load_binned(axis: str) -> dict[int, list[dict[str, float]]]:
    fine_csv = BASE / "ppg12_exact_fine_reweight_histograms.csv"
    if fine_csv.exists() and fine_csv.stat().st_size > 0:
        out: dict[int, list[dict[str, float]]] = {0: [], 1: []}
        with fine_csv.open() as f:
            for row in csv.DictReader(f):
                if row["axis"] != axis:
                    continue
                cls = int(row["class"])
                out[cls].append(
                    {
                        "lo": float(row["bin_low"]),
                        "hi": float(row["bin_high"]),
                        "raw": float(row["raw_density"]),
                        "weighted": float(row["weighted_density"]),
                    }
                )
        for cls in (0, 1):
            out[cls].sort(key=lambda r: r["lo"])
        return out

    grouped: dict[tuple[int, float, float], dict[str, float]] = defaultdict(lambda: {"n": 0.0, "w": 0.0})
    with (BASE / "ppg12_exact_sample_inventory_binned.csv").open() as f:
        for row in csv.DictReader(f):
            if row["axis"] != axis:
                continue
            key = (int(row["class"]), float(row["bin_low"]), float(row["bin_high"]))
            grouped[key]["n"] += float(row["n_rows"])
            grouped[key]["w"] += float(row["sum_ppg12_exact_weight"])
    out: dict[int, list[dict[str, float]]] = {0: [], 1: []}
    for cls in (0, 1):
        total_n = sum(v["n"] for (c, _, _), v in grouped.items() if c == cls)
        total_w = sum(v["w"] for (c, _, _), v in grouped.items() if c == cls)
        for (c, lo, hi), v in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1])):
            if c != cls:
                continue
            width = hi - lo
            out[cls].append(
                {
                    "lo": lo,
                    "hi": hi,
                    "raw": v["n"] / total_n / width if total_n > 0 and width > 0 else 0.0,
                    "weighted": v["w"] / total_w / width if total_w > 0 and width > 0 else 0.0,
                }
            )
    return out


def panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill: str = "white") -> None:
    draw.rounded_rectangle(box, radius=22, fill=fill, outline="#CBD5E1", width=2)


def nice_ticks(vmin: float, vmax: float, n: int = 4) -> list[float]:
    return [vmin + i * (vmax - vmin) / n for i in range(n + 1)]


def centers(rows: list[dict[str, float]]) -> list[float]:
    return [0.5 * (r["lo"] + r["hi"]) for r in rows]


def max_density(data: dict[int, list[dict[str, float]]], field: str, pad: float = 1.15) -> float:
    return max(max(r[field] for r in data[0]), max(r[field] for r in data[1])) * pad


def draw_axes(
    draw: ImageDraw.ImageDraw,
    plot: tuple[int, int, int, int],
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    xlabel: str,
    ylabel: str,
    xticks: list[float],
    yticks: list[float],
) -> tuple[callable, callable]:
    left, top, right, bottom = plot
    f_tick = font(17)
    f_axis = font(21)

    def sx(x: float) -> float:
        return left + (x - xlim[0]) / (xlim[1] - xlim[0]) * (right - left)

    def sy(y: float) -> float:
        return bottom - (y - ylim[0]) / (ylim[1] - ylim[0]) * (bottom - top)

    for y in yticks:
        yy = sy(y)
        draw.line((left, yy, right, yy), fill=GRID, width=2)
        label = f"{y:.2f}" if ymax_small(ylim) else f"{y:.1f}"
        draw.text((left - 64, yy - 11), label, font=f_tick, fill=MUTED)
    for x in xticks:
        xx = sx(x)
        draw.line((xx, top, xx, bottom), fill="#F1F5F9", width=1)
        if abs(x) < 2:
            label = f"{x:.1f}"
        else:
            label = f"{x:.0f}"
        draw.text((xx - text_width(draw, label, f_tick) / 2, bottom + 11), label, font=f_tick, fill=MUTED)
    draw.rectangle((left, top, right, bottom), outline=INK, width=2)
    draw_math_text(draw, ((left + right) // 2 - text_width(draw, xlabel, f_axis) // 2, bottom + 40), xlabel, f_axis, fill=INK)
    if ylabel:
        draw.text((left, top - 32), ylabel, font=f_axis, fill=INK)
    return sx, sy


def ymax_small(ylim: tuple[float, float]) -> bool:
    return ylim[1] <= 1.2


def draw_density_plot(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    data: dict[int, list[dict[str, float]]],
    field: str,
    title: str,
    xlim: tuple[float, float],
    xlabel: str,
    ylabel: str,
    ylim: tuple[float, float] | None = None,
    legend_position: str = "top-right",
) -> None:
    panel(draw, box)
    x0, y0, x1, y1 = box
    draw_math_text(draw, (x0 + 24, y0 + 18), title, font(31, True), fill=INK)
    left, top, right, bottom = x0 + 92, y0 + 78, x1 - 30, y1 - 76
    if ylim is None:
        ylim = (0.0, max_density(data, field))
    xticks = nice_ticks(*xlim, n=4)
    yticks = nice_ticks(*ylim, n=4)
    sx, sy = draw_axes(draw, (left, top, right, bottom), xlim, ylim, xlabel, ylabel, xticks, yticks)

    for cls, color, label in [(1, SIGNAL_RED, "signal"), (0, BACKGROUND_BLUE, "background")]:
        rows = data[cls]
        pts: list[tuple[float, float]] = []
        for r in rows:
            pts.extend([(sx(r["lo"]), sy(r[field])), (sx(r["hi"]), sy(r[field]))])
        draw.line(pts, fill=color, width=5)
        for a, b in zip(rows[:-1], rows[1:]):
            xx = sx(a["hi"])
            draw.line((xx, sy(a[field]), xx, sy(b[field])), fill=color, width=5)

    if legend_position == "mid-right":
        lx, ly = right - 220, top + 80
    elif legend_position == "bottom-right":
        lx, ly = right - 220, bottom - 72
    else:
        lx, ly = right - 220, top + 18
    draw.rounded_rectangle((lx - 18, ly - 22, lx + 210, ly + 66), radius=14, fill=(255, 255, 255, 230), outline="#E2E8F0", width=1)
    draw.line((lx, ly, lx + 56, ly), fill=SIGNAL_RED, width=5)
    draw.text((lx + 68, ly - 13), "signal", font=font(21), fill=INK)
    draw.line((lx, ly + 34, lx + 56, ly + 34), fill=BACKGROUND_BLUE, width=5)
    draw.text((lx + 68, ly + 21), "background", font=font(21), fill=INK)


def draw_ratio_plot(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    data: dict[int, list[dict[str, float]]],
    title: str,
    xlim: tuple[float, float],
    xlabel: str,
    max_dev_percent: float,
) -> None:
    panel(draw, box)
    x0, y0, x1, y1 = box
    draw_math_text(draw, (x0 + 24, y0 + 18), title, font(31, True), fill=INK)
    left, top, right, bottom = x0 + 92, y0 + 78, x1 - 30, y1 - 76
    ylim = (0.88, 1.12)
    xticks = nice_ticks(*xlim, n=4)
    yticks = [0.90, 0.95, 1.00, 1.05, 1.10]
    sx, sy = draw_axes(draw, (left, top, right, bottom), xlim, ylim, xlabel, "S/B", xticks, yticks)

    # Pale closure band helps the audience read the target without treating it as a fit.
    draw.rectangle((left, sy(1.05), right, sy(0.95)), fill="#DCFCE7")
    for y in yticks:
        yy = sy(y)
        draw.line((left, yy, right, yy), fill=GRID if abs(y - 1.0) > 1e-6 else INK, width=2 if abs(y - 1.0) > 1e-6 else 3)
    draw.rectangle((left, top, right, bottom), outline=INK, width=2)

    pts: list[tuple[float, float] | None] = []
    for sig, bkg in zip(data[1], data[0]):
        x = 0.5 * (sig["lo"] + sig["hi"])
        if bkg["weighted"] <= 0 or sig["weighted"] <= 0:
            pts.append(None)
            continue
        ratio = sig["weighted"] / bkg["weighted"]
        if not math.isfinite(ratio):
            pts.append(None)
            continue
        pts.append((sx(x), sy(ratio)))

    segment: list[tuple[float, float]] = []
    for pt in pts + [None]:
        if pt is None:
            if len(segment) >= 2:
                draw.line(segment, fill=PURPLE, width=5)
            segment = []
        else:
            segment.append(pt)

    for pt in pts:
        if pt is None:
            continue
        x, y = pt
        draw.ellipse((x - 8, y - 8, x + 8, y + 8), fill=PURPLE, outline="white", width=2)

    badge = f"max |S/B - 1| = {max_dev_percent:.1f}%"
    bw = text_width(draw, badge, font(24, True)) + 42
    draw.rounded_rectangle((right - bw, top + 14, right - 16, top + 62), radius=16, fill="#F5F3FF", outline="#DDD6FE", width=2)
    draw.text((right - bw + 21, top + 25), badge, font=font(24, True), fill="#5B21B6")
    draw.text((left + 20, top + 18), "green band: +/-5%", font=font(20), fill="#166534")


def draw_chip(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], number: str, text: str) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=28, fill="#EFF6FF", outline="#BFDBFE", width=2)
    draw.ellipse((x0 + 22, y0 + 21, x0 + 66, y0 + 65), fill=BLUE)
    draw.text((x0 + 44 - text_width(draw, number, font(24, True)) / 2, y0 + 28), number, font=font(24, True), fill="white")
    draw_wrapped(draw, text, (x0 + 84, y0 + 18), x1 - x0 - 106, font(25, True), fill="#1E3A8A", line_gap=2)


def draw_metric(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], value: str, label: str, color: str) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=24, fill="#F8FAFC", outline="#CBD5E1", width=2)
    draw.text((x0 + 24, y0 + 22), value, font=font(43, True), fill=color)
    draw_wrapped(draw, label, (x0 + 24, y0 + 82), x1 - x0 - 48, font(21), fill="#334155", line_gap=2)


def draw_title_with_et(draw: ImageDraw.ImageDraw) -> None:
    title_font = font(70, True)
    sub_font = font(43, True)
    x, y = 78, 48
    draw.text((x, y), "E", font=title_font, fill=INK)
    x += text_width(draw, "E", title_font) + 2
    draw.text((x, y + 39), "T", font=sub_font, fill=INK)
    x += text_width(draw, "T", sub_font) + 4
    draw.text((x, y), "/eta-reweighting Au+Au, Embedded Sim", font=title_font, fill=INK)


def main() -> None:
    global BASE, OUT_DIR, OUT, SUBTITLE, CONTEXT_LABEL, ET_XMIN, ET_XMAX
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE, help="Directory with ppg12_exact_* metadata/CSV inputs.")
    parser.add_argument("--out", type=Path, default=OUT, help="Output slide-sized PNG path.")
    parser.add_argument("--subtitle", default=SUBTITLE, help="Subtitle text below the title.")
    parser.add_argument("--context-label", default=CONTEXT_LABEL, help="Optional context label drawn near the title.")
    parser.add_argument("--et-xmin", type=float, default=ET_XMIN, help="Lower E_T axis bound for the E_T panels.")
    parser.add_argument("--et-xmax", type=float, default=ET_XMAX, help="Upper E_T axis bound for the E_T panels.")
    args = parser.parse_args()
    BASE = args.base
    OUT = args.out
    OUT_DIR = OUT.parent
    SUBTITLE = args.subtitle
    CONTEXT_LABEL = args.context_label
    ET_XMIN = args.et_xmin
    ET_XMAX = args.et_xmax

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = json.loads((BASE / "ppg12_exact_reweighting_metadata.json").read_text())
    weighting = meta["weighting"]
    et_data = load_binned("cluster_Et")
    eta_data = load_binned("cluster_Eta")

    sum0 = float(weighting["sum_weight_class0"])
    sum1 = float(weighting["sum_weight_class1"])
    imbalance = abs(sum1 - sum0) / ((sum1 + sum0) / 2.0) * 100.0
    et_ratio = [s["weighted"] / b["weighted"] for s, b in zip(et_data[1], et_data[0])]
    eta_ratio = [s["weighted"] / b["weighted"] for s, b in zip(eta_data[1], eta_data[0])]
    et_dev = max(abs(r - 1.0) for r in et_ratio) * 100.0
    eta_dev = max(abs(r - 1.0) for r in eta_ratio) * 100.0

    slide = Image.new("RGBA", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)

    draw_title_with_et(draw)
    draw_math_text(
        draw,
        (80, 134),
        SUBTITLE,
        font(32),
        fill=MUTED,
    )
    if CONTEXT_LABEL:
        label_font = font(28, True)
        label_w = text_width(draw, CONTEXT_LABEL, label_font) + 50
        draw.rounded_rectangle((2480 - label_w, 62, 2480, 114), radius=18, fill="#F8FAFC", outline="#CBD5E1", width=2)
        draw.text((2505 - label_w, 75), CONTEXT_LABEL, font=label_font, fill="#334155")

    draw_chip(draw, (80, 214, 575, 306), "1", "equal total signal and background weight")
    draw_chip(draw, (600, 214, 1095, 306), "2", "flatten eta separately for each truth class")
    draw_chip(draw, (1120, 214, 1615, 306), "3", "flatten E_T separately for each truth class")
    draw.rounded_rectangle((1660, 206, 2480, 316), radius=28, fill="#FFF7ED", outline="#FED7AA", width=2)
    draw.text((1692, 228), "Do not mix training weights with physics weights", font=font(31, True), fill="#9A3412")
    draw.text((1692, 274), "No event, cross-section, stitching, vertex, or centrality weights.", font=font(26), fill="#9A3412")

    draw_density_plot(
        draw,
        (80, 360, 1215, 760),
        et_data,
        "raw",
        "Before weighting: E_T populations are not comparable",
        (ET_XMIN, ET_XMAX),
        "cluster E_T [GeV]",
        "density",
        (0.0, 0.22),
        legend_position="mid-right",
    )
    draw_ratio_plot(
        draw,
        (1280, 360, 2480, 760),
        et_data,
        "After weighting: E_T signal/background ratio",
        (ET_XMIN, ET_XMAX),
        "cluster E_T [GeV]",
        et_dev,
    )
    draw_density_plot(
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
    draw_ratio_plot(
        draw,
        (1280, 820, 2480, 1190),
        eta_data,
        "After weighting: eta signal/background ratio",
        (-0.7, 0.7),
        "cluster eta",
        eta_dev,
    )

    draw_metric(draw, (80, 1238, 560, 1380), f"{imbalance:.2f}%", "residual total S/B weight imbalance", GREEN)
    et_range_label = f"largest E_T closure deviation over {ET_XMIN:g}-{ET_XMAX:g} GeV"
    draw_metric(draw, (590, 1238, 1070, 1380), f"{et_dev:.1f}%", et_range_label, PURPLE)
    draw_metric(draw, (1100, 1238, 1580, 1380), f"{eta_dev:.1f}%", "largest eta closure deviation over |eta| < 0.7", PURPLE)

    draw.rounded_rectangle((1610, 1238, 2480, 1380), radius=24, fill="#ECFDF5", outline="#A7F3D0", width=2)
    draw.text((1644, 1264), "Main readout", font=font(31, True), fill="#065F46")
    draw_wrapped(
        draw,
        "After reweighting, the BDT comparison is not dominated by the raw E_T/eta population difference.",
        (1644, 1312),
        780,
        font(25),
        fill="#064E3B",
        line_gap=3,
    )

    slide.convert("RGB").save(OUT, quality=95)
    print(OUT)


if __name__ == "__main__":
    main()
