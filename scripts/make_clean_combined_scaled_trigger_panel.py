#!/usr/bin/env python3
"""Render a combined 1x2 scaled-trigger QA panel from summed bin CSV."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau/scaledTriggerRunByRunQA/"
    "scaled_trigger_run_by_run_20260521_202504"
)
IN_CSV = BASE / "clean_efficiency_turnon_tables" / "combined_clean63_scaled_trigger_bins.csv"
OUT_PNG = BASE / "clean_efficiency_turnon_tables" / "combined_clean63_scaled_trigger_1x2.png"
OUT_SUMMARY = BASE / "clean_efficiency_turnon_tables" / "combined_clean63_scaled_trigger_summary.txt"

FONT_REG = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"
FONT_ITALIC = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"


def font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size=size)


W, H = 2400, 1400
INK = (18, 24, 32)
MUTED = (82, 92, 105)
GRID = (221, 226, 232)
AXIS = (38, 44, 51)
BLACK = (0, 0, 0)
BLUE = (0, 34, 220)
RED = (220, 0, 0)
GREEN = (22, 105, 67)
GRAY = (92, 98, 108)
FIT_X_MAX = 16.0
FIT_ERROR_FLOOR = 0.04

F_TITLE = font(FONT_BOLD, 58)
F_PANEL = font(FONT_BOLD, 34)
F_TEXT = font(FONT_REG, 28)
F_SMALL = font(FONT_REG, 23)
F_TINY = font(FONT_REG, 20)
F_LABEL = font(FONT_REG, 27)
F_SPHX = font(FONT_BOLD, 27)
F_ITAL = font(FONT_ITALIC, 27)
F_SUBSCRIPT = font(FONT_ITALIC, 18)


class TurnOnFit:
    def __init__(self, floor: float, plateau: float, x50: float, width: float, x90: float, used_points: int, loss: float) -> None:
        self.floor = floor
        self.plateau = plateau
        self.x50 = x50
        self.width = width
        self.x90 = x90
        self.used_points = used_points
        self.loss = loss

    def value(self, x: float) -> float:
        arg = max(min(-(x - self.x50) / self.width, 60.0), -60.0)
        return self.floor + (self.plateau - self.floor) / (1.0 + math.exp(arg))


def load_bins() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with IN_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["BIN"] != "BIN":
                continue
            rows.append(
                {
                    "bin": float(row["bin"]),
                    "center": float(row["center"]),
                    "width": float(row["width"]),
                    "mbd": float(row["mbd"]),
                    "p10": float(row["p10"]),
                    "p12": float(row["p12"]),
                }
            )
    if not rows:
        raise RuntimeError(f"No BIN rows found in {IN_CSV}")
    return rows


def integral(rows: list[dict[str, float]], key: str, lo: float, hi: float) -> float:
    return sum(row[key] for row in rows if lo <= row["center"] < hi)


def ratio_err(num: float, den: float) -> tuple[float, float]:
    if den <= 0:
        return math.nan, math.nan
    ratio = num / den
    if num <= 0:
        return ratio, math.sqrt(max(num, 1.0)) / den
    return ratio, ratio * math.sqrt(1.0 / num + 1.0 / den)


def fit_turnon(rows: list[dict[str, float]], key: str) -> TurnOnFit:
    # Free-plateau robust sigmoid. Fit the dense turn-on and plateau shoulder,
    # but leave sparse high-E tail points visible as data rather than using
    # them to drive the guide curve.
    data: list[tuple[float, float, float]] = []
    for row in rows:
        if row["mbd"] <= 0 or not (1.0 <= row["center"] <= FIT_X_MAX):
            continue
        value, err = ratio_err(row[key], row["mbd"])
        if math.isfinite(value):
            data.append((row["center"], value, max(err, FIT_ERROR_FLOOR)))

    def model(params: list[float], xval: float) -> float:
        floor, plateau, x50, width = params
        arg = max(min(-(xval - x50) / width, 60.0), -60.0)
        return floor + (plateau - floor) / (1.0 + math.exp(arg))

    def loss(params: list[float]) -> float:
        floor, plateau, x50, width = params
        if not (0.0 <= floor <= 0.08 and 0.82 <= plateau <= 1.08 and 3.0 <= x50 <= 9.5 and 0.35 <= width <= 4.0):
            return 1.0e99
        if plateau <= floor + 0.20:
            return 1.0e99
        total = 0.0
        for xval, yval, err in data:
            residual = (model(params, xval) - yval) / err
            total += 2.0 * (math.sqrt(1.0 + residual * residual) - 1.0)
        return total

    starts = [
        [floor, plateau, x50, width]
        for floor in [0.0, 0.004, 0.01, 0.02]
        for plateau in [0.88, 0.91, 0.94, 0.97, 1.00]
        for x50 in [5.0, 5.5, 6.0, 6.5, 7.0]
        for width in [0.7, 1.0, 1.3, 1.7, 2.2]
    ]
    best_params = starts[0]
    best_loss = loss(best_params)
    for params in starts[1:]:
        current = loss(params)
        if current < best_loss:
            best_loss = current
            best_params = params

    steps = [0.010, 0.010, 0.10, 0.10]
    for _ in range(90):
        improved = False
        for idx, step in enumerate(steps):
            for sign in [-1.0, 1.0]:
                trial = best_params[:]
                trial[idx] += sign * step
                current = loss(trial)
                if current < best_loss:
                    best_loss = current
                    best_params = trial
                    improved = True
        if not improved:
            steps = [step * 0.70 for step in steps]

    floor, plateau, x50, width = best_params
    x90 = x50 + width * math.log(0.9 / 0.1)
    return TurnOnFit(floor, plateau, x50, width, x90, len(data), best_loss)


def text_size(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def draw_collision_label(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    prefix = "Internal  Au+Au, "
    draw.text((x, y), prefix, font=F_ITAL, fill=INK)
    px, _ = text_size(draw, prefix, F_ITAL)
    sx = x + px
    draw.text((sx, y), "√s", font=F_ITAL, fill=INK)
    root_w, _ = text_size(draw, "√s", F_ITAL)
    draw.text((sx + root_w - 1, y + 16), "NN", font=F_SUBSCRIPT, fill=INK)
    draw.text((sx + root_w + 24, y), "=200 GeV", font=F_ITAL, fill=INK)


class Panel:
    def __init__(self, x0: int, y0: int, x1: int, y1: int, *, xlo: float, xhi: float, ylo: float, yhi: float, logy: bool = False) -> None:
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.logy = logy

    def x(self, value: float) -> int:
        frac = (value - self.xlo) / (self.xhi - self.xlo)
        return int(round(self.x0 + frac * (self.x1 - self.x0)))

    def y(self, value: float) -> int:
        if self.logy:
            value = max(value, self.ylo)
            frac = (math.log10(value) - math.log10(self.ylo)) / (math.log10(self.yhi) - math.log10(self.ylo))
        else:
            frac = (value - self.ylo) / (self.yhi - self.ylo)
        return int(round(self.y1 - frac * (self.y1 - self.y0)))


def draw_axes(
    img: Image.Image,
    draw: ImageDraw.ImageDraw,
    panel: Panel,
    *,
    xlabel: str,
    ylabel: str,
    x_ticks: list[float],
    y_ticks: list[float],
    y_labels: list[str] | None = None,
) -> None:
    draw.rectangle((panel.x0, panel.y0, panel.x1, panel.y1), outline=AXIS, width=3)
    for tick in x_ticks:
        x = panel.x(tick)
        draw.line((x, panel.y1, x, panel.y1 + 10), fill=AXIS, width=2)
        label = f"{tick:g}"
        tw, th = text_size(draw, label, F_SMALL)
        draw.text((x - tw // 2, panel.y1 + 14), label, font=F_SMALL, fill=INK)
        if panel.xlo < tick < panel.xhi:
            draw.line((x, panel.y0, x, panel.y1), fill=GRID, width=1)
    labels = y_labels if y_labels is not None else [f"{tick:g}" for tick in y_ticks]
    for tick, label in zip(y_ticks, labels):
        y = panel.y(tick)
        draw.line((panel.x0 - 10, y, panel.x0, y), fill=AXIS, width=2)
        tw, th = text_size(draw, label, F_SMALL)
        draw.text((panel.x0 - 16 - tw, y - th // 2), label, font=F_SMALL, fill=INK)
        if panel.ylo < tick < panel.yhi:
            draw.line((panel.x0, y, panel.x1, y), fill=GRID, width=1)
    tw, th = text_size(draw, xlabel, F_LABEL)
    draw.text(((panel.x0 + panel.x1 - tw) // 2, panel.y1 + 62), xlabel, font=F_LABEL, fill=INK)
    label_w, label_h = text_size(draw, ylabel, F_LABEL)
    label_img = Image.new("RGBA", (label_w + 10, label_h + 10), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label_img)
    label_draw.text((5, 5), ylabel, font=F_LABEL, fill=INK)
    label_img = label_img.rotate(90, expand=True)
    img.paste(label_img, (panel.x0 - 132, (panel.y0 + panel.y1 - label_img.height) // 2), label_img)


def step_points(rows: list[dict[str, float]], key: str, panel: Panel) -> list[tuple[int, int]]:
    points: list[tuple[int, int]] = []
    for idx, row in enumerate(rows):
        lo = row["center"] - row["width"] / 2.0
        hi = row["center"] + row["width"] / 2.0
        if hi < panel.xlo or lo > panel.xhi:
            continue
        y = panel.y(max(row[key], panel.ylo))
        xlo = panel.x(max(lo, panel.xlo))
        xhi = panel.x(min(hi, panel.xhi))
        if not points:
            points.append((xlo, y))
        else:
            prev_x, prev_y = points[-1]
            if prev_x != xlo:
                points.append((xlo, prev_y))
            points.append((xlo, y))
        points.append((xhi, y))
    return points


def draw_marker(draw: ImageDraw.ImageDraw, x: int, y: int, color: tuple[int, int, int], open_marker: bool = False) -> None:
    r = 5
    if open_marker:
        draw.ellipse((x - r, y - r, x + r, y + r), outline=color, width=3)
    else:
        draw.ellipse((x - r, y - r, x + r, y + r), fill=color, outline=color)


def draw_efficiency(draw: ImageDraw.ImageDraw, panel: Panel, rows: list[dict[str, float]], key: str, color: tuple[int, int, int], *, open_marker: bool = False) -> None:
    for row in rows:
        center = row["center"]
        if center < panel.xlo or center > panel.xhi:
            continue
        den = row["mbd"]
        num = row[key]
        if den <= 0:
            continue
        value, err = ratio_err(num, den)
        if math.isnan(value):
            continue
        x = panel.x(center)
        y = panel.y(value)
        xlo = panel.x(center - row["width"] / 2.0)
        xhi = panel.x(center + row["width"] / 2.0)
        ylo = panel.y(max(panel.ylo, value - err))
        yhi = panel.y(min(panel.yhi, value + err))
        draw.line((xlo, y, xhi, y), fill=color, width=2)
        draw.line((x, ylo, x, yhi), fill=color, width=2)
        draw.line((x - 5, ylo, x + 5, ylo), fill=color, width=2)
        draw.line((x - 5, yhi, x + 5, yhi), fill=color, width=2)
        draw_marker(draw, x, y, color, open_marker=open_marker)


def draw_fit_curve(draw: ImageDraw.ImageDraw, panel: Panel, fit: TurnOnFit, color: tuple[int, int, int]) -> None:
    points: list[tuple[int, int]] = []
    n_steps = 190
    for idx in range(n_steps + 1):
        xval = panel.xlo + idx * (panel.xhi - panel.xlo) / n_steps
        yval = min(max(fit.value(xval), panel.ylo), panel.yhi)
        points.append((panel.x(xval), panel.y(yval)))
    draw.line(points, fill=color, width=4)


def draw_dashed_vertical(draw: ImageDraw.ImageDraw, panel: Panel, xval: float, ytop: float, color: tuple[int, int, int]) -> None:
    x = panel.x(xval)
    y0 = panel.y(0.0)
    y1 = panel.y(min(max(ytop, panel.ylo), panel.yhi))
    step = 9
    yy = y0
    while yy > y1:
        draw.line((x, yy, x, max(y1, yy - 5)), fill=color, width=2)
        yy -= step


def draw_legend(draw: ImageDraw.ImageDraw, x: int, y: int, items: list[tuple[str, tuple[int, int, int], str]]) -> None:
    for idx, (label, color, style) in enumerate(items):
        yy = y + idx * 38
        if style == "step":
            draw.line((x, yy + 13, x + 46, yy + 13), fill=color, width=5)
        elif style == "open":
            draw_marker(draw, x + 22, yy + 13, color, open_marker=True)
        else:
            draw_marker(draw, x + 22, yy + 13, color, open_marker=False)
        draw.text((x + 58, yy), label, font=F_SMALL, fill=INK)


def main() -> int:
    rows = load_bins()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)

    tail_mbd = integral(rows, "mbd", 15.0, 20.0)
    tail_p10 = integral(rows, "p10", 15.0, 20.0)
    tail_p12 = integral(rows, "p12", 15.0, 20.0)
    low_mbd = integral(rows, "mbd", 1.0, 3.0)
    low_p10 = integral(rows, "p10", 1.0, 3.0)
    low_p12 = integral(rows, "p12", 1.0, 3.0)
    mid_mbd = integral(rows, "mbd", 6.0, 9.0)
    mid_p10 = integral(rows, "p10", 6.0, 9.0)
    mid_p12 = integral(rows, "p12", 6.0, 9.0)

    r10_tail = tail_p10 / tail_mbd
    r12_tail = tail_p12 / tail_mbd
    r10_low = low_p10 / low_mbd
    r12_low = low_p12 / low_mbd
    r10_mid = mid_p10 / mid_mbd
    r12_mid = mid_p12 / mid_mbd
    fit10 = fit_turnon(rows, "p10")
    fit12 = fit_turnon(rows, "p12")

    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)

    draw.text((72, 46), "Combined clean full-efficiency scaled-trigger QA", font=F_TITLE, fill=INK)
    draw.text(
        (76, 115),
        "63 selected runs summed bin-by-bin before forming Trigger/MBD; RHS fit is a free robust sigmoid over the dense turn-on region.",
        font=F_TEXT,
        fill=MUTED,
    )
    draw.text((76, 160), "sPHENIX", font=F_SPHX, fill=INK)
    sx, _ = text_size(draw, "sPHENIX", F_SPHX)
    draw_collision_label(draw, 82 + sx, 160)

    left = Panel(205, 310, 1088, 1090, xlo=1.0, xhi=20.0, ylo=1.0e2, yhi=1.0e9, logy=True)
    right = Panel(1388, 310, 2268, 1090, xlo=1.0, xhi=20.0, ylo=0.0, yhi=2.1)

    draw.text((left.x0, 255), "Summed max-cluster energy overlay", font=F_PANEL, fill=INK)
    draw_axes(
        img,
        draw,
        left,
        xlabel="Max cluster energy [GeV], Eclus > 1 GeV",
        ylabel="Live/scaled-corrected counts",
        x_ticks=[2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
        y_ticks=[1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9],
        y_labels=[r"10^2", r"10^3", r"10^4", r"10^5", r"10^6", r"10^7", r"10^8", r"10^9"],
    )
    for key, color in [("mbd", BLACK), ("p10", BLUE), ("p12", RED)]:
        pts = step_points(rows, key, left)
        if len(pts) >= 2:
            draw.line(pts, fill=color, width=4)
    draw.text((left.x0 + 28, left.y0 + 32), "MBD N&S >= 2, |vz| < 150 cm", font=F_SMALL, fill=INK)
    draw_legend(
        draw,
        left.x0 + 490,
        left.y0 + 68,
        [
            ("MBD reference", BLACK, "step"),
            ("Photon 10, scaled", BLUE, "step"),
            ("Photon 12, scaled", RED, "step"),
        ],
    )

    draw.text((right.x0, 255), "Combined turn-on efficiency", font=F_PANEL, fill=INK)
    draw_axes(
        img,
        draw,
        right,
        xlabel="Max cluster energy [GeV], Eclus > 1 GeV",
        ylabel="Trigger / MBD efficiency",
        x_ticks=[2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
        y_ticks=[0.0, 0.5, 1.0, 1.5, 2.0],
    )
    unity_y = right.y(1.0)
    draw.line((right.x0, unity_y, right.x1, unity_y), fill=(175, 180, 188), width=4)
    draw_fit_curve(draw, right, fit10, BLUE)
    draw_fit_curve(draw, right, fit12, RED)
    draw_dashed_vertical(draw, right, fit10.x90, fit10.value(fit10.x90), BLUE)
    draw_dashed_vertical(draw, right, fit12.x90, fit12.value(fit12.x90), RED)
    draw_efficiency(draw, right, rows, "p10", BLUE)
    draw_efficiency(draw, right, rows, "p12", RED, open_marker=True)
    draw.text((right.x0 + 24, right.y0 + 22), "sPHENIX", font=F_SPHX, fill=INK)
    rhs_sx, _ = text_size(draw, "sPHENIX", F_SPHX)
    draw.text((right.x0 + 30 + rhs_sx, right.y0 + 22), "Internal", font=F_ITAL, fill=INK)
    draw.text(
        (right.x0 + 24, right.y0 + 62),
        f"Photon 10 plateau = {fit10.plateau:.3f}; 90% = {fit10.x90:.2f} GeV",
        font=F_SMALL,
        fill=BLUE,
    )
    draw.text(
        (right.x0 + 24, right.y0 + 94),
        f"Photon 12 plateau = {fit12.plateau:.3f}; 90% = {fit12.x90:.2f} GeV",
        font=F_SMALL,
        fill=RED,
    )
    draw.text(
        (right.x(13.4), right.y(0.76)),
        "Fit: free robust sigmoid",
        font=F_SMALL,
        fill=INK,
    )
    draw.text(
        (right.x(13.4), right.y(0.68)),
        "1<Emax<16 GeV",
        font=F_SMALL,
        fill=INK,
    )
    draw.text(
        (right.x(13.4), right.y(0.60)),
        "tail points shown, not fitted",
        font=F_SMALL,
        fill=INK,
    )
    draw_legend(
        draw,
        right.x0 + 500,
        right.y0 + 575,
        [
            (f"Photon 10 / MBD, 90%={fit10.x90:.2f}", BLUE, "closed"),
            (f"Photon 12 / MBD, 90%={fit12.x90:.2f}", RED, "open"),
        ],
    )

    stats = (
        f"Clean group: 63/620 runs = {63 / 620 * 100:.1f}%    "
        f"low 1-3 GeV P10/P12={r10_low:.4f}/{r12_low:.4f}    "
        f"mid 6-9 GeV P10/P12={r10_mid:.3f}/{r12_mid:.3f}    "
        f"tail 15-20 GeV P10/P12={r10_tail:.3f}/{r12_tail:.3f}"
    )
    draw.rounded_rectangle((72, 1270, W - 72, 1350), radius=10, fill=(235, 247, 240), outline=(189, 215, 199), width=2)
    draw.text((96, 1291), stats, font=F_TEXT, fill=INK)

    img.save(OUT_PNG, quality=95)
    with OUT_SUMMARY.open("w") as f:
        f.write("Combined clean full-efficiency scaled-trigger QA\n")
        f.write(f"Input bins: {IN_CSV}\n")
        f.write("Run group: clean_full_turn_on\n")
        f.write("Runs: 63 of 620 (10.2%)\n")
        f.write(f"MBD tail 15-20 GeV: {tail_mbd:.6g}\n")
        f.write(f"Photon10/MBD tail 15-20 GeV: {r10_tail:.6f}\n")
        f.write(f"Photon12/MBD tail 15-20 GeV: {r12_tail:.6f}\n")
        f.write(f"Photon10/MBD low 1-3 GeV: {r10_low:.6f}\n")
        f.write(f"Photon12/MBD low 1-3 GeV: {r12_low:.6f}\n")
        f.write(f"Photon10/MBD mid 6-9 GeV: {r10_mid:.6f}\n")
        f.write(f"Photon12/MBD mid 6-9 GeV: {r12_mid:.6f}\n")
        f.write(f"Photon10 fit floor: {fit10.floor:.6f}\n")
        f.write("Fit method: free robust sigmoid; floor/plateau/x50/width fitted over 1-16 GeV; sparse high-E tail points shown but not fitted\n")
        f.write(f"Fit x window GeV: 1.000-{FIT_X_MAX:.3f}\n")
        f.write(f"Fit error floor: {FIT_ERROR_FLOOR:.6f}\n")
        f.write(f"Photon10 fit plateau: {fit10.plateau:.6f}\n")
        f.write(f"Photon10 fit x50 GeV: {fit10.x50:.6f}\n")
        f.write(f"Photon10 fit width GeV: {fit10.width:.6f}\n")
        f.write(f"Photon10 fit x90 GeV: {fit10.x90:.6f}\n")
        f.write(f"Photon10 fit points used: {fit10.used_points}\n")
        f.write(f"Photon10 fit robust loss: {fit10.loss:.6f}\n")
        f.write(f"Photon12 fit floor: {fit12.floor:.6f}\n")
        f.write(f"Photon12 fit plateau: {fit12.plateau:.6f}\n")
        f.write(f"Photon12 fit x50 GeV: {fit12.x50:.6f}\n")
        f.write(f"Photon12 fit width GeV: {fit12.width:.6f}\n")
        f.write(f"Photon12 fit x90 GeV: {fit12.x90:.6f}\n")
        f.write(f"Photon12 fit points used: {fit12.used_points}\n")
        f.write(f"Photon12 fit robust loss: {fit12.loss:.6f}\n")
    print(OUT_PNG)
    print(OUT_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
