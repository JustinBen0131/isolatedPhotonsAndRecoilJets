#!/usr/bin/env python3
"""Build a compact pp exact-stitch acceptance slide from local QA artifacts."""

from __future__ import annotations

import csv
import json
import math
import textwrap
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppPhotonMLPipeline/"
    "ppg12_basev3E_currentIAN_exactStitch_20260523_0915"
)
QA_DIR = ROOT / "validation" / "insitu_stitching"
SLIDE_DIR = ROOT / "slide_assets"

CONTRACT_JSON = QA_DIR / "pp_currentian_exactstitch_contract_summary.json"
BOUNDARY_CSV = QA_DIR / "pp_currentian_exactstitch_boundary_continuity_qa.csv"
OUT_CSV = QA_DIR / "pp_currentian_exactstitch_acceptance_summary.csv"
OUT_JSON = QA_DIR / "pp_currentian_exactstitch_acceptance_summary.json"
OUT_PNG = SLIDE_DIR / "pp_exactstitch_slide10_acceptance_summary.png"

W, H = 2400, 1350


def font(size: int, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    candidates = []
    if bold and italic:
        candidates += [
            "/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf",
            "/System/Library/Fonts/Times.ttc",
        ]
    elif bold:
        candidates += [
            "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
            "/System/Library/Fonts/Times.ttc",
        ]
    elif italic:
        candidates += [
            "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",
            "/System/Library/Fonts/Times.ttc",
        ]
    else:
        candidates += [
            "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
            "/System/Library/Fonts/Times.ttc",
        ]
    candidates.append("/System/Library/Fonts/Supplemental/Arial.ttf")
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


F_TITLE = font(48, bold=True)
F_SUB = font(24)
F_HEAD = font(21, bold=True)
F_CELL = font(20)
F_SMALL = font(18)
F_TINY = font(16)


def fmt_num(x: float, sig: int = 4) -> str:
    if x == 0:
        return "0"
    ax = abs(x)
    if ax >= 1e5 or ax < 1e-2:
        return f"{x:.{sig - 1}e}"
    if ax >= 1000:
        return f"{x:,.0f}"
    if ax >= 100:
        return f"{x:.1f}"
    if ax >= 10:
        return f"{x:.2f}"
    return f"{x:.3f}"


def load_contract() -> dict:
    with CONTRACT_JSON.open() as f:
        return json.load(f)


def load_boundary_rows() -> list[dict]:
    rows: list[dict] = []
    with BOUNDARY_CSV.open() as f:
        for row in csv.DictReader(f):
            parsed = dict(row)
            for key in [
                "boundary_GeV",
                "bin_width_GeV",
                "display_adjacent_left_over_right",
                "max_same_bin_overlap_fractional_deviation",
            ]:
                parsed[key] = float(parsed[key])
            parsed["pass_boundary"] = str(parsed["pass_boundary"]).lower() == "true"
            rows.append(parsed)
    return rows


def sample_rows(contract: dict) -> list[dict]:
    rows = []
    group_order = {"photon": 0, "jet": 1}
    for s in sorted(contract["samples"], key=lambda r: (group_order.get(r["group"], 9), r["stitch_window"][0])):
        weight = float(s["xsec_pb"]) / float(s["events_processed_metadata"]) / float(s["hist_bin_width"])
        rows.append(
            {
                "group": s["group"],
                "sample": s["sample"].replace("run28_", ""),
                "window": f"[{fmt_num(float(s['stitch_window'][0]))}, {fmt_num(float(s['stitch_window'][1]))})",
                "bin_width_GeV": float(s["hist_bin_width"]),
                "xsec_pb": float(s["xsec_pb"]),
                "events": float(s["events_processed_metadata"]),
                "weight_pb_per_event_GeV": weight,
            }
        )
    return rows


def write_acceptance_summary(contract: dict, boundaries: list[dict], samples: list[dict]) -> None:
    max_dev = max(float(r["max_same_bin_overlap_fractional_deviation"]) for r in boundaries)
    summary = {
        "schema": "PP_CURRENTIAN_EXACTSTITCH_ACCEPTANCE_SUMMARY_V1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_contract_json": str(CONTRACT_JSON),
        "source_boundary_csv": str(BOUNDARY_CSV),
        "normalization": "density = raw_all_events * xsec_pb / events_processed_metadata / bin_width",
        "stitch_variable": {
            "photon": "max truth photon Et from Justin/RecoilJets in-situ histograms",
            "inclusive_jet": "max R=0.4 truth jet pT from Justin/RecoilJets in-situ histograms",
        },
        "overall_pass": all(r["pass_boundary"] for r in boundaries),
        "max_same_bin_overlap_fractional_deviation": max_dev,
        "samples": samples,
        "boundaries": boundaries,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n")

    with OUT_CSV.open("w", newline="") as f:
        fieldnames = [
            "group",
            "sample",
            "window",
            "bin_width_GeV",
            "xsec_pb",
            "events",
            "weight_pb_per_event_GeV",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(samples)


def rounded(draw: ImageDraw.ImageDraw, xy, fill, outline="#d9dee7", radius=20, width=2):
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy, s: str, fnt, fill="#111827", anchor=None, spacing=4):
    draw.multiline_text(xy, s, font=fnt, fill=fill, anchor=anchor, spacing=spacing)


def fit_text(s: str, max_chars: int) -> str:
    return "\n".join(textwrap.wrap(s, width=max_chars, break_long_words=False))


def draw_table(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    h: int,
    headers: list[str],
    rows: list[list[str]],
    col_fracs: list[float],
    header_fill="#eef3fb",
    row_font=F_CELL,
    head_font=F_HEAD,
):
    rounded(draw, (x, y, x + w, y + h), "white", "#cfd6e3", radius=18, width=2)
    n_rows = len(rows) + 1
    row_h = h / n_rows
    col_x = [x]
    for frac in col_fracs[:-1]:
        col_x.append(col_x[-1] + int(round(w * frac)))
    col_w = [int(round(w * frac)) for frac in col_fracs]

    draw.rounded_rectangle((x, y, x + w, int(y + row_h)), radius=18, fill=header_fill, outline="#cfd6e3", width=1)
    draw.rectangle((x, int(y + row_h * 0.55), x + w, int(y + row_h)), fill=header_fill)

    for i, htxt in enumerate(headers):
        text(draw, (col_x[i] + 12, y + row_h / 2), htxt, head_font, "#1e2a36", anchor="lm")

    for ridx, row in enumerate(rows):
        yy = y + row_h * (ridx + 1)
        if ridx % 2 == 1:
            draw.rectangle((x + 1, int(yy), x + w - 1, int(yy + row_h)), fill="#f8fafc")
        for i, value in enumerate(row):
            text(draw, (col_x[i] + 12, yy + row_h / 2), value, row_font, "#25313d", anchor="lm")

    for cx in col_x[1:]:
        draw.line((cx, y, cx, y + h), fill="#e2e8f0", width=1)
    for i in range(n_rows + 1):
        yy = int(y + i * row_h)
        draw.line((x, yy, x + w, yy), fill="#e2e8f0", width=1)


def make_slide(contract: dict, boundaries: list[dict], samples: list[dict]) -> None:
    img = Image.new("RGB", (W, H), "#fbfbf8")
    draw = ImageDraw.Draw(img)

    text(draw, (108, 72), "pp exact-stitch acceptance check", F_TITLE, "#111827")
    text(
        draw,
        (108, 135),
        "This Analysis Output: Justin/RecoilJets in-situ histograms, explicit xsec/Nevt/bin-width weights, non-overlapping stitch windows",
        F_SUB,
        "#475569",
    )

    pass_all = all(r["pass_boundary"] for r in boundaries)
    rounded(draw, (1835, 72, 2295, 172), "#e8f8ee" if pass_all else "#fff0f0", "#55b878" if pass_all else "#d94c4c", 20, 3)
    text(draw, (1875, 105), "Boundary QA: PASS" if pass_all else "Boundary QA: FAIL", font(29, bold=True), "#166534")
    max_dev = max(100.0 * float(r["max_same_bin_overlap_fractional_deviation"]) for r in boundaries)
    text(draw, (1875, 142), f"max same-bin deviation = {max_dev:.1f}%", F_SMALL, "#166534")

    sample_table = []
    for s in samples:
        sample_table.append(
            [
                s["sample"].replace("photonjet", "pho").replace("jet", "jet"),
                s["window"],
                fmt_num(s["bin_width_GeV"]),
                fmt_num(s["xsec_pb"]),
                fmt_num(s["events"]),
                fmt_num(s["weight_pb_per_event_GeV"], 3),
            ]
        )

    draw_table(
        draw,
        108,
        250,
        2184,
        500,
        ["sample", "stitch window [GeV]", "bin", "xsec [pb]", "Nevt", "weight"],
        sample_table,
        [0.14, 0.22, 0.08, 0.16, 0.17, 0.23],
        row_font=font(20),
        head_font=font(21, bold=True),
    )

    interp = (
        "Interpretation: the final stitched spectra are rate-normalized by sample cross section and generated-event count. "
        "The boundary check compares adjacent samples at the same truth-pT bin before window zeroing; factor-level mistakes would appear here."
    )
    text(draw, (108, 790), fit_text(interp, 132), font(21), "#334155", spacing=7)

    boundary_table = []
    labels = []
    values = []
    groups = []
    for b in boundaries:
        pair = (
            b["left_sample"].replace("run28_", "").replace("photonjet", "pho")
            + " -> "
            + b["right_sample"].replace("run28_", "").replace("photonjet", "pho")
        )
        dev_pct = 100.0 * float(b["max_same_bin_overlap_fractional_deviation"])
        boundary_table.append(
            [
                b["group"],
                fmt_num(float(b["boundary_GeV"])),
                pair,
                f"{float(b['display_adjacent_left_over_right']):.3f}",
                f"{dev_pct:.2f}%",
                "pass" if b["pass_boundary"] else "fail",
            ]
        )
        labels.append(f"{b['group']} {fmt_num(float(b['boundary_GeV']))}")
        values.append(dev_pct)
        groups.append(b["group"])

    draw_table(
        draw,
        108,
        930,
        1370,
        330,
        ["group", "edge", "adjacent samples", "display L/R", "overlap dev", "QA"],
        boundary_table,
        [0.12, 0.10, 0.33, 0.15, 0.17, 0.13],
        row_font=font(18),
        head_font=font(19, bold=True),
    )

    text(draw, (1595, 914), "Adjacent-sample agreement at stitch edges", font(22, bold=True), "#111827")
    chart_x, chart_y, chart_w, chart_h = 1600, 965, 610, 245
    draw.line((chart_x, chart_y + chart_h, chart_x + chart_w, chart_y + chart_h), fill="#334155", width=2)
    draw.line((chart_x, chart_y, chart_x, chart_y + chart_h), fill="#334155", width=2)
    threshold_x = chart_x + int(chart_w * 15.0 / 16.0)
    draw.line((threshold_x, chart_y, threshold_x, chart_y + chart_h), fill="#334155", width=2)
    text(draw, (threshold_x - 5, chart_y + chart_h + 8), "15% gate", F_TINY, "#334155", anchor="ra")
    bar_h = 24
    gap = 14
    for i, (label, value, group) in enumerate(zip(labels, values, groups)):
        yy = chart_y + 10 + i * (bar_h + gap)
        color = "#3b82f6" if group == "photon" else "#f97316"
        draw.rectangle((chart_x, yy, chart_x + int(chart_w * value / 16.0), yy + bar_h), fill=color)
        text(draw, (chart_x - 14, yy + bar_h / 2), label, F_TINY, "#334155", anchor="rm")
        text(draw, (chart_x + int(chart_w * value / 16.0) + 8, yy + bar_h / 2), f"{value:.1f}%", F_TINY, "#334155", anchor="lm")
    text(draw, (chart_x, chart_y + chart_h + 8), "0", F_TINY, "#334155")
    text(draw, (chart_x + chart_w, chart_y + chart_h + 8), "16%", F_TINY, "#334155", anchor="ra")

    rounded(draw, (1588, 748, 2292, 895), "#f8fafc", "#cfd6e3", 18, 2)
    text(draw, (1625, 778), "Acceptance rule", font(23, bold=True), "#111827")
    text(
        draw,
        (1625, 815),
        "1. no stitch-window gaps or overlaps\n2. adjacent samples agree within 15% at the edge\n3. displayed spectrum falls smoothly across the boundary",
        font(18),
        "#334155",
        spacing=5,
    )

    SLIDE_DIR.mkdir(parents=True, exist_ok=True)
    img.save(OUT_PNG)


def main() -> None:
    contract = load_contract()
    boundaries = load_boundary_rows()
    samples = sample_rows(contract)
    write_acceptance_summary(contract, boundaries, samples)
    make_slide(contract, boundaries, samples)
    print(OUT_CSV)
    print(OUT_JSON)
    print(OUT_PNG)


if __name__ == "__main__":
    main()
