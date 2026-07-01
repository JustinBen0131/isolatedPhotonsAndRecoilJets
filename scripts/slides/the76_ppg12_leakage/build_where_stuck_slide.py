#!/usr/bin/env python3
"""Build the THE-76 PPG12 leakage parity status slide as a 2560x1440 PNG."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
    "ppg12_photon_yield_v1_data_20260620/purity_fig29_comparison/"
    "ppg12_ratio_diagnostic/leakage_components"
)
CSV_PATH = BASE / "current_vs_ppg12_fig29_leakage_components_bcd.csv"
OUT_DIR = BASE / "slide_candidates"
PNG_PATH = OUT_DIR / "the76_where_stuck_ppg12_leakage_bcd_slide_20260629.png"
SCRIPT_PATH = OUT_DIR / "the76_where_stuck_ppg12_leakage_bcd_slide_20260629_speaker.md"
MANIFEST_PATH = OUT_DIR / "the76_where_stuck_ppg12_leakage_bcd_slide_20260629_manifest.json"
AUDIT_PATH = OUT_DIR / "the76_where_stuck_ppg12_leakage_bcd_slide_20260629_slide_audit.json"
LAYOUT_PATH = OUT_DIR / "the76_where_stuck_ppg12_leakage_bcd_slide_20260629_layout_nodes.json"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"


def font(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(path), size=size)


F_TITLE = font(TIMES_BOLD, 68)
F_HEADER = font(TIMES_BOLD, 40)
F_BODY = font(TIMES, 36)
F_BODY_SMALL = font(TIMES, 37)
F_PANEL = font(TIMES_BOLD, 36)
F_AXIS = font(TIMES, 31)
F_AXIS_SMALL = font(TIMES, 28)
F_LEGEND = font(TIMES, 29)
F_FOOT = font(TIMES, 37)


def read_rows() -> list[dict[str, float]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def sha256(path: Path) -> str:
    return subprocess.check_output(["shasum", "-a", "256", str(path)], text=True).split()[0]


def draw_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font_obj: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int] = (0, 0, 0),
    *,
    anchor: str = "la",
) -> tuple[int, int, int, int]:
    draw.text(xy, text, font=font_obj, fill=fill, anchor=anchor)
    return draw.textbbox(xy, text, font=font_obj, anchor=anchor)


def bold_label_line(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    label: str,
    body: str,
    *,
    size: int = 36,
) -> list[dict]:
    x, y = xy
    f_bold = font(TIMES_BOLD, size)
    f_regular = font(TIMES, size)
    bbox_label = draw.textbbox((x, y), label, font=f_bold, anchor="la")
    draw.text((x, y), label, font=f_bold, fill=(0, 0, 0), anchor="la")
    draw.text((bbox_label[2] + 8, y), body, font=f_regular, fill=(0, 0, 0), anchor="la")
    bbox_body = draw.textbbox((bbox_label[2] + 8, y), body, font=f_regular, anchor="la")
    return [
        {
            "kind": "text",
            "role": "audience",
            "name": label.strip(),
            "text": label + " " + body,
            "bbox": [x, bbox_label[1], bbox_body[2], bbox_body[3]],
            "font_px": size,
            "text_runs": [
                {"text": label, "bold": True},
                {"text": " " + body, "bold": False},
            ],
        }
    ]


def map_xy(
    x: float,
    y: float,
    plot: tuple[int, int, int, int],
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
) -> tuple[float, float]:
    x0, y0, x1, y1 = plot
    px = x0 + (x - xmin) / (xmax - xmin) * (x1 - x0)
    py = y1 - (y - ymin) / (ymax - ymin) * (y1 - y0)
    return px, py


def draw_marker(
    draw: ImageDraw.ImageDraw,
    x: float,
    y: float,
    *,
    color: tuple[int, int, int],
    shape: str,
    size: int = 9,
    open_marker: bool = True,
) -> None:
    if shape == "circle":
        box = [x - size, y - size, x + size, y + size]
        draw.ellipse(box, outline=color, fill=None if open_marker else color, width=2)
    elif shape == "square":
        box = [x - size, y - size, x + size, y + size]
        draw.rectangle(box, outline=color, fill=None if open_marker else color, width=2)
    elif shape == "triangle":
        pts = [(x, y - size * 1.1), (x - size * 1.05, y + size * 0.9), (x + size * 1.05, y + size * 0.9)]
        draw.line([pts[0], pts[1], pts[2], pts[0]], fill=color, width=2)


def draw_axis_label_rotated(
    img: Image.Image,
    text: str,
    center: tuple[int, int],
    font_obj: ImageFont.FreeTypeFont,
) -> None:
    tmp = Image.new("RGBA", (420, 60), (255, 255, 255, 0))
    d = ImageDraw.Draw(tmp)
    bbox = d.textbbox((0, 0), text, font=font_obj)
    d.text(((420 - (bbox[2] - bbox[0])) // 2, (60 - (bbox[3] - bbox[1])) // 2 - bbox[1]), text, font=font_obj, fill=(0, 0, 0, 255))
    rot = tmp.rotate(90, expand=True)
    img.alpha_composite(rot, (center[0] - rot.width // 2, center[1] - rot.height // 2))


def draw_panel(
    img: Image.Image,
    draw: ImageDraw.ImageDraw,
    rows: list[dict[str, float]],
    rect: tuple[int, int, int, int],
    title: str,
    current_key: str,
    ppg_key: str | None,
    ratio_key: str | None,
    *,
    ymin: float,
    ymax: float,
    y_ticks: list[float],
    color: tuple[int, int, int],
    shape: str,
    show_y_label: bool,
    y_label: str,
    show_x_label: bool,
    is_ratio: bool = False,
) -> None:
    x0, y0, x1, y1 = rect
    draw.rectangle(rect, outline=(0, 0, 0), width=2)
    for xt in [10, 15, 20, 25, 30, 35]:
        px, _ = map_xy(xt, ymin, rect, 9.5, 36.5, ymin, ymax)
        draw.line([(px, y1), (px, y1 + 12)], fill=(0, 0, 0), width=2)
        draw.line([(px, y0), (px, y0 - 10)], fill=(0, 0, 0), width=1)
        draw.text((px, y1 + 17), f"{xt}", font=F_AXIS_SMALL, fill=(0, 0, 0), anchor="ma")
    for yt in y_ticks:
        _, py = map_xy(10, yt, rect, 9.5, 36.5, ymin, ymax)
        draw.line([(x0 - 12, py), (x0, py)], fill=(0, 0, 0), width=2)
        draw.line([(x1, py), (x1 + 10, py)], fill=(0, 0, 0), width=1)
        label = f"{yt:.2f}" if ymax <= 0.2 else f"{yt:.1f}"
        draw.text((x0 - 18, py), label, font=F_AXIS_SMALL, fill=(0, 0, 0), anchor="rm")
    if is_ratio:
        px0, py = map_xy(9.5, 1.0, rect, 9.5, 36.5, ymin, ymax)
        px1, _ = map_xy(36.5, 1.0, rect, 9.5, 36.5, ymin, ymax)
        for xs in range(int(px0), int(px1), 18):
            draw.line([(xs, py), (xs + 9, py)], fill=(120, 120, 120), width=2)

    if title:
        draw.text((x0 - 72, y0 - 48), title, font=F_PANEL, fill=(0, 0, 0), anchor="la")
    if show_y_label:
        draw_axis_label_rotated(img, y_label, (x0 - 78, (y0 + y1) // 2), F_AXIS)
    if show_x_label:
        draw.text(((x0 + x1) // 2, y1 + 55), "reco cluster ET [GeV]", font=F_AXIS, fill=(0, 0, 0), anchor="ma")

    if ppg_key:
        for row in rows:
            px, py = map_xy(row["pt_center"], row[ppg_key], rect, 9.5, 36.5, ymin, ymax)
            draw_marker(draw, px, py, color=(0, 0, 0), shape="circle", size=8, open_marker=False)
    key = ratio_key if ratio_key else current_key
    for row in rows:
        px, py = map_xy(row["pt_center"], row[key], rect, 9.5, 36.5, ymin, ymax)
        draw_marker(draw, px, py, color=color, shape=shape, size=8, open_marker=True)


def build_slide() -> None:
    rows = read_rows()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    img = Image.new("RGBA", (2560, 1440), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)
    nodes: list[dict] = []

    bbox = draw_text(draw, (116, 56), "Where we are stuck: PPG12 leakage parity", F_TITLE)
    nodes.append({"kind": "text", "role": "title", "name": "slide title", "text": "Where we are stuck: PPG12 leakage parity", "bbox": list(bbox), "title_anchor": True, "font_px": 68})
    blue_fill, blue_line = (239, 246, 255), (82, 125, 190)
    gold_fill, gold_line = (255, 249, 237), (213, 145, 38)
    left_box = (116, 138, 1285, 296)
    right_box = (1335, 138, 2428, 296)
    draw.rectangle(left_box, fill=blue_fill, outline=blue_line, width=2)
    draw.rectangle(right_box, fill=gold_fill, outline=gold_line, width=2)
    nodes += bold_label_line(draw, (162, 170), "Remaining issue:", "no leakage region tracks PPG12 cleanly yet.", size=37)
    nodes += bold_label_line(draw, (162, 228), "Goal:", "reproduce PPG12 pp yield inside PhotonClusterBuilder.", size=37)
    nodes += bold_label_line(draw, (1380, 170), "PPG12 quirks:", "G4/double, MBD vertex, topo iso.", size=37)
    nodes += bold_label_line(draw, (1380, 228), "Still suspect:", "BDT/candidate semantics; high-ET binning.", size=37)
    nodes += bold_label_line(draw, (116, 330), "Black points:", "direct SDCC PPG12 ROOT extract; 32-36 bin is not visibly exposed in IAN Fig. 28.", size=37)

    blue, green, orange = (31, 119, 180), (44, 160, 44), (217, 95, 2)
    col_x = [245, 978, 1765]
    top_y0, top_y1 = 455, 850
    bot_y0, bot_y1 = 912, 1218
    w = 565
    panels = [(x, top_y0, x + w, top_y1) for x in col_x] + [(x, bot_y0, x + w, bot_y1) for x in col_x]

    draw_panel(img, draw, rows, panels[0], "Region B", "current_f_b", "ppg12_f_b", None, ymin=0, ymax=0.165, y_ticks=[0, 0.05, 0.10, 0.15], color=blue, shape="circle", show_y_label=True, y_label="leakage / A", show_x_label=False)
    draw_panel(img, draw, rows, panels[1], "Region C", "current_f_c", "ppg12_f_c", None, ymin=0, ymax=0.90, y_ticks=[0, 0.2, 0.4, 0.6, 0.8], color=green, shape="square", show_y_label=False, y_label="", show_x_label=False)
    draw_panel(img, draw, rows, panels[2], "Region D", "current_f_d", "ppg12_f_d", None, ymin=0, ymax=0.095, y_ticks=[0, 0.02, 0.04, 0.06, 0.08], color=orange, shape="triangle", show_y_label=False, y_label="", show_x_label=False)
    draw_panel(img, draw, rows, panels[3], "", "", None, "current_f_b_over_ppg12", ymin=0.30, ymax=1.10, y_ticks=[0.4, 0.6, 0.8, 1.0], color=blue, shape="circle", show_y_label=True, y_label="current / PPG12", show_x_label=True, is_ratio=True)
    draw_panel(img, draw, rows, panels[4], "", "", None, "current_f_c_over_ppg12", ymin=0.30, ymax=1.10, y_ticks=[0.4, 0.6, 0.8, 1.0], color=green, shape="square", show_y_label=False, y_label="", show_x_label=True, is_ratio=True)
    draw_panel(img, draw, rows, panels[5], "", "", None, "current_f_d_over_ppg12", ymin=0.30, ymax=1.10, y_ticks=[0.4, 0.6, 0.8, 1.0], color=orange, shape="triangle", show_y_label=False, y_label="", show_x_label=True, is_ratio=True)

    # Legends in the upper-left open space of each top panel.
    for rect, ccol, shp in [(panels[0], blue, "circle"), (panels[1], green, "square"), (panels[2], orange, "triangle")]:
        lx, ly = rect[0] + 42, rect[1] + 28
        draw_marker(draw, lx, ly, color=(0, 0, 0), shape="circle", size=7, open_marker=False)
        draw.text((lx + 28, ly - 16), "PPG12 reference", font=F_LEGEND, fill=(0, 0, 0))
        draw_marker(draw, lx, ly + 40, color=ccol, shape=shp, size=7, open_marker=True)
        draw.text((lx + 28, ly + 24), "Current rerun", font=F_LEGEND, fill=(0, 0, 0))

    x, y = 116, 1344
    read_label = "Readout:"
    read_body = "do not over-interpret the visible IAN panel; compare exact ROOT bins, then run same-cluster stage parity."
    bbox_label = draw.textbbox((x, y), read_label, font=font(TIMES_BOLD, 37), anchor="la")
    draw.text((x, y), read_label, font=font(TIMES_BOLD, 37), fill=(100, 100, 100), anchor="la")
    draw.text((bbox_label[2] + 8, y), " " + read_body, font=F_FOOT, fill=(100, 100, 100), anchor="la")
    bbox_body = draw.textbbox((bbox_label[2] + 8, y), " " + read_body, font=F_FOOT, anchor="la")
    nodes.append({
        "kind": "text",
        "role": "audience",
        "name": "readout",
        "text": read_label + " " + read_body,
        "bbox": [x, bbox_label[1], bbox_body[2], bbox_body[3]],
        "font_px": 37,
        "text_runs": [
            {"text": read_label, "bold": True},
            {"text": " " + read_body, "bold": False},
        ],
    })

    img.convert("RGB").save(PNG_PATH)
    LAYOUT_PATH.write_text(json.dumps({"title_axis_x": 116, "nodes": nodes}, indent=2), encoding="utf-8")

    SCRIPT_PATH.write_text(
        """# THE-76 Slide Script - Where We Are Stuck In PPG12 Leakage Parity

Here I want to be clear about where this comparison is stuck. This is not a clean PPG12 match yet, and the disagreement becomes worse at higher reconstructed cluster energy.

The left column is region B, the middle column is region C, and the right column is region D. The upper panels show the leakage fraction relative to region A, and the lower panels show the current result divided by the PPG12 reference.

The black points are extracted directly from Shuhang's SDCC PPG12 ROOT files. The exact ROOT histograms contain a 32 to 36 GeV bin, but the visible current IAN Figure 28 panel does not expose that high bin cleanly. So the correct comparison must be made against the ROOT bin contents, not by eye from the embedded PDF figure.

The important progress is that we have learned several concrete PPG12 mechanics that matter. The component mixture, the double-sample G4 reconstruction path, the MBD vertex choices, the topo-isolation settings, and the PPG12 BDT sideband logic all had to be treated carefully to move the result closer.

But the remaining issue is broader than just one leakage region. None of the regions track PPG12 cleanly across the full ROOT-bin range, and the high-ET behavior must be treated carefully because the last bin is not visually obvious in the IAN panel.

The purpose of this exercise is not to rerun Shuhang's code as a black box. The purpose is to understand exactly what PPG12 did and transfer the relevant pp photon-yield logic into our PhotonClusterBuilder baseline, so the pp reference is unified, controlled, and less convoluted. The next useful step is therefore same-cluster PPG12-versus-RecoilJets stage parity, especially the candidate tagging and score-sideband semantics before isolation.
""",
        encoding="utf-8",
    )

    png_hash = sha256(PNG_PATH)
    csv_hash = sha256(CSV_PATH)
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "artifact": "THE-76 where stuck PPG12 leakage parity slide candidate",
                "generated_at_local": datetime.now().astimezone().isoformat(timespec="seconds"),
                "png": str(PNG_PATH),
                "speaker_script": str(SCRIPT_PATH),
                "generator": str(Path(__file__).resolve()),
                "input_csv": str(CSV_PATH),
                "layout_nodes": str(LAYOUT_PATH),
                "audit_json": str(AUDIT_PATH),
                "png_sha256": png_hash,
                "input_csv_sha256": csv_hash,
                "dimensions_px": [2560, 1440],
                "font_policy": "Times New Roman used for full-slide title, callouts, labels, and axis text",
                "status": "local full-slide PNG candidate; no Google Slides mutation performed",
                "on_slide_claim": "No leakage region cleanly tracks the PPG12 ROOT-bin reference yet; the visible IAN Fig. 28 panel does not expose the 32-36 GeV bin cleanly, so exact ROOT-bin provenance is required before interpreting high-ET disagreement",
                "source_context": "THE-76 pp-SIM signal-only global-MBD component-mix leakage diagnostic",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG_PATH)


if __name__ == "__main__":
    build_slide()
