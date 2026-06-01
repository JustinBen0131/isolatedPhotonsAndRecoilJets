#!/usr/bin/env python3
"""Build RHS-efficiency-only review sheets for scaled-trigger categories."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import csv
import math
import textwrap
from collections import Counter
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau/scaledTriggerRunByRunQA/"
    "scaled_trigger_run_by_run_20260521_202504"
)
CSV_PATH = BASE / "run_shape_categories.csv"
OUT_DIR = BASE / "efficiency_category_review"

FONT_REG = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"
W, H = 3840, 2160

CATEGORIES = [
    {
        "key": "clean_full_turn_on",
        "label": "Clean full efficiency turn-on",
        "rule": "low ~0, turn-on rises, tail in [0.97, 1.08]",
        "check": "This is the accepted efficient shape. Both P10/MBD and P12/MBD should approach 1 only after the turn-on.",
        "color": (28, 111, 70),
        "fill": (229, 246, 236),
    },
    {
        "key": "missing_photon_trigger_data",
        "label": "Missing / sparse photon-trigger efficiency",
        "rule": "trigger numerator is absent or too sparse for a meaningful curve",
        "check": "Look for mostly empty RHS panels, ratios near zero, or very sparse trigger entries.",
        "color": (153, 45, 45),
        "fill": (255, 246, 244),
    },
    {
        "key": "flat_unity_from_threshold",
        "label": "Flat 100% from threshold",
        "rule": "low-E efficiency is already ~1, so there is no physical turn-on",
        "check": "These should not be counted as efficient turn-ons, even though the tail is unity.",
        "color": (171, 91, 23),
        "fill": (255, 249, 238),
    },
    {
        "key": "early_turn_on_outlier",
        "label": "Early-on tail-unity outlier",
        "rule": "tail is near unity, but low-E efficiency starts too high",
        "check": "71260-like: not flat at 1, but already significantly above zero before the expected turn-on.",
        "color": (184, 72, 42),
        "fill": (255, 242, 235),
    },
    {
        "key": "tail_inefficient",
        "label": "Severely inefficient high-E tail",
        "rule": "15-20 GeV tail efficiency < 0.80",
        "check": "Usually has a real turn-on shape but does not reach full efficiency.",
        "color": (128, 57, 140),
        "fill": (247, 240, 250),
    },
    {
        "key": "tail_overscaled",
        "label": "Severely overscaled high-E tail",
        "rule": "15-20 GeV tail efficiency > 1.25",
        "check": "The high-E points overshoot unity too strongly to be treated as a clean efficiency plateau.",
        "color": (40, 101, 154),
        "fill": (236, 246, 255),
    },
    {
        "key": "tail_intermediate_deviation",
        "label": "Moderate tail deviation",
        "rule": "tail outside [0.92, 1.15], excluding severe failures",
        "check": "These are the main threshold-tuning boundary between failed and maybe-usable runs.",
        "color": (81, 87, 99),
        "fill": (244, 245, 247),
    },
    {
        "key": "mild_tail_deviation",
        "label": "Near-unity but not strict clean",
        "rule": "tail near unity but fails the strict [0.97, 1.08] clean window",
        "check": "These are close in the tail; human review should decide whether this bucket is too strict.",
        "color": (81, 87, 99),
        "fill": (244, 245, 247),
    },
    {
        "key": "low_stat_tail",
        "label": "Low-stat high-E tail",
        "rule": "MBD tail < 100 counts, so tail efficiency is noisy/inconclusive",
        "check": "Do not over-interpret the high-E efficiency because the denominator is small.",
        "color": (81, 87, 99),
        "fill": (244, 245, 247),
    },
]


def font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size=size)


F_TITLE = font(FONT_BOLD, 72)
F_SUB = font(FONT_REG, 35)
F_H1 = font(FONT_BOLD, 42)
F_H2 = font(FONT_BOLD, 33)
F_BODY = font(FONT_REG, 29)
F_BODY_B = font(FONT_BOLD, 29)
F_SMALL = font(FONT_REG, 23)
F_SMALL_B = font(FONT_BOLD, 24)
F_TINY = font(FONT_REG, 20)

INK = (20, 26, 33)
MUTED = (83, 91, 100)
LINE = (202, 207, 214)
LIGHT = (231, 234, 238)


def load_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with CSV_PATH.open() as f:
        for row in csv.DictReader(f):
            parsed: dict[str, object] = dict(row)
            for key in [
                "run",
                "mbd_tail",
                "r10_tail",
                "r12_tail",
                "r10_low",
                "r12_low",
                "r10_mid",
                "r12_mid",
                "r10_total",
                "r12_total",
            ]:
                parsed[key] = int(row[key]) if key == "run" else float(row[key])
            parsed["png_path"] = BASE / str(row["png"])
            rows.append(parsed)
    return rows


def sorted_rows(rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(rows, key=lambda row: (float(row["mbd_tail"]), int(row["run"])), reverse=True)


def category_info(key: str) -> dict[str, object]:
    for cat in CATEGORIES:
        if cat["key"] == key:
            return cat
    raise KeyError(key)


def eff_label(row: dict[str, object]) -> str:
    return (
        f"Run {row['run']}  low={float(row['r10_low']):.3f}/{float(row['r12_low']):.3f}  "
        f"mid={float(row['r10_mid']):.2f}/{float(row['r12_mid']):.2f}  "
        f"tail={float(row['r10_tail']):.2f}/{float(row['r12_tail']):.2f}"
    )


def efficiency_crop(row: dict[str, object], max_size: tuple[int, int]) -> Image.Image:
    path = Path(row["png_path"])  # type: ignore[arg-type]
    im = Image.open(path).convert("RGB")
    # Keep only the RHS efficiency panel plus its title/legend. The crop is
    # intentionally generous because ROOT margins vary slightly across panels.
    x0 = int(im.width * 0.50)
    crop = im.crop((x0, 20, im.width - 8, im.height - 8))
    crop.thumbnail(max_size, Image.Resampling.LANCZOS)
    return crop


def draw_wrapped(draw: ImageDraw.ImageDraw, text: str, x: int, y: int, width: int, font_obj, fill=INK, line_h=33) -> int:
    for idx, line in enumerate(textwrap.wrap(text, width=width)):
        draw.text((x, y + idx * line_h), line, font=font_obj, fill=fill)
    return y + max(1, len(textwrap.wrap(text, width=width))) * line_h


def draw_eff_cell(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    row: dict[str, object],
    x: int,
    y: int,
    cell_w: int,
    cell_h: int,
    color: tuple[int, int, int],
    fill: tuple[int, int, int],
) -> None:
    draw.rounded_rectangle((x, y, x + cell_w, y + cell_h), radius=8, outline=color, width=3, fill=(254, 254, 254))
    draw.rectangle((x + 2, y + 2, x + cell_w - 2, y + 50), fill=fill)
    draw.text((x + 14, y + 13), eff_label(row), font=F_SMALL_B, fill=color)
    thumb = efficiency_crop(row, (cell_w - 32, cell_h - 80))
    canvas.paste(thumb, (x + (cell_w - thumb.width) // 2, y + 62 + ((cell_h - 82) - thumb.height) // 2))


def make_overview(rows: list[dict[str, object]], counts: Counter[str]) -> Path:
    out = OUT_DIR / "00_efficiency_category_overview.png"
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    d.text((86, 58), "Scaled-trigger categories by actual efficiency behavior", font=F_TITLE, fill=INK)
    d.text((88, 150), "Every example below is the RHS Trigger/MBD panel only; categories are defined by low/mid/tail efficiency ratios.", font=F_SUB, fill=MUTED)
    d.line((88, 220, W - 88, 220), fill=LINE, width=3)

    card_w, card_h = 1190, 565
    gap_x, gap_y = 56, 42
    top = 282
    for idx, cat in enumerate(CATEGORIES):
        key = str(cat["key"])
        color = cat["color"]  # type: ignore[assignment]
        fill = cat["fill"]  # type: ignore[assignment]
        examples = sorted_rows(row for row in rows if row["category"] == key)[:3]
        x = 88 + (idx % 3) * (card_w + gap_x)
        y = top + (idx // 3) * (card_h + gap_y)
        d.rounded_rectangle((x, y, x + card_w, y + card_h), radius=12, fill=fill, outline=color, width=3)
        d.text((x + 24, y + 18), f"{cat['label']}: {counts[key]}", font=F_H1, fill=color)
        y2 = draw_wrapped(d, f"Rule: {cat['rule']}", x + 24, y + 76, 65, F_SMALL, MUTED, 29)
        y2 = draw_wrapped(d, f"Check: {cat['check']}", x + 24, y2 + 8, 65, F_SMALL, INK, 29)
        for j, row in enumerate(examples):
            thumb = efficiency_crop(row, (350, 280))
            px = x + 24 + j * 375
            py = y + 248
            img.paste(thumb, (px, py))
            d.text((px, py + thumb.height + 8), f"Run {row['run']}", font=F_TINY, fill=INK)
            d.text((px, py + thumb.height + 31), f"tail {float(row['r10_tail']):.2f}/{float(row['r12_tail']):.2f}", font=F_TINY, fill=MUTED)
    img.save(out, quality=95)
    return out


def make_boundary(rows: list[dict[str, object]], counts: Counter[str]) -> Path:
    out = OUT_DIR / "00_efficiency_boundary_tuning.png"
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    d.text((86, 58), "Efficiency boundary checks for iterative tuning", font=F_TITLE, fill=INK)
    d.text((88, 150), "Use this sheet to decide which category boundaries are too strict, too loose, or mislabeled.", font=F_SUB, fill=MUTED)
    d.line((88, 220, W - 88, 220), fill=LINE, width=3)

    panels = [
        ("clean_full_turn_on", "Accepted efficient"),
        ("mild_tail_deviation", "Near-unity tail, not strict"),
        ("tail_intermediate_deviation", "Moderate tail deviation"),
        ("tail_inefficient", "Severe low tail"),
        ("tail_overscaled", "Severe high tail"),
        ("missing_photon_trigger_data", "Missing/sparse efficiency"),
        ("flat_unity_from_threshold", "Flat from threshold"),
        ("early_turn_on_outlier", "Early-on before threshold"),
    ]
    cell_w, cell_h = 1780, 430
    gap_x, gap_y = 84, 42
    top = 284
    for idx, (key, title) in enumerate(panels):
        cat = category_info(key)
        color = cat["color"]  # type: ignore[assignment]
        fill = cat["fill"]  # type: ignore[assignment]
        x = 88 + (idx % 2) * (cell_w + gap_x)
        y = top + (idx // 2) * (cell_h + gap_y)
        d.rounded_rectangle((x, y, x + cell_w, y + cell_h), radius=10, fill=fill, outline=color, width=3)
        d.text((x + 24, y + 16), f"{title}: {counts[key]}", font=F_H2, fill=color)
        examples = sorted_rows(row for row in rows if row["category"] == key)[:4]
        for j, row in enumerate(examples):
            thumb = efficiency_crop(row, (395, 285))
            px = x + 24 + j * 430
            py = y + 82
            img.paste(thumb, (px, py))
            d.text((px, py + thumb.height + 8), f"Run {row['run']} tail {float(row['r10_tail']):.2f}/{float(row['r12_tail']):.2f}", font=F_TINY, fill=INK)
    img.save(out, quality=95)
    return out


def make_category_pages(rows: list[dict[str, object]], counts: Counter[str]) -> list[Path]:
    paths: list[Path] = []
    per_page = 12
    cols, rows_per_page = 4, 3
    left, top = 70, 342
    gap_x, gap_y = 28, 34
    cell_w = (W - 2 * left - (cols - 1) * gap_x) // cols
    cell_h = (H - top - 92 - (rows_per_page - 1) * gap_y) // rows_per_page

    for cat in CATEGORIES:
        key = str(cat["key"])
        color = cat["color"]  # type: ignore[assignment]
        fill = cat["fill"]  # type: ignore[assignment]
        cat_rows = sorted_rows(row for row in rows if row["category"] == key)
        total_pages = max(1, math.ceil(len(cat_rows) / per_page))
        for page in range(total_pages):
            out = OUT_DIR / f"{key}_eff_page_{page + 1:02d}_of_{total_pages:02d}.png"
            img = Image.new("RGB", (W, H), "white")
            d = ImageDraw.Draw(img)
            d.text((86, 58), f"{cat['label']} ({counts[key]} runs)", font=F_TITLE, fill=color)
            d.text((88, 145), f"Efficiency rule: {cat['rule']}  Page {page + 1} of {total_pages}.", font=F_SUB, fill=MUTED)
            d.line((88, 220, W - 88, 220), fill=LINE, width=3)
            draw_wrapped(d, f"What to check: {cat['check']}", 88, 245, 150, F_BODY, INK, 34)
            page_rows = cat_rows[page * per_page : (page + 1) * per_page]
            for idx, row in enumerate(page_rows):
                col = idx % cols
                rpage = idx // cols
                x = left + col * (cell_w + gap_x)
                y = top + rpage * (cell_h + gap_y)
                draw_eff_cell(img, d, row, x, y, cell_w, cell_h, color, fill)
            img.save(out, quality=95)
            paths.append(out)
    return paths


def write_index(rows: list[dict[str, object]], counts: Counter[str], overview: Path, boundary: Path, pages: list[Path]) -> Path:
    out = OUT_DIR / "README.md"
    by_key: dict[str, list[Path]] = {str(cat["key"]): [] for cat in CATEGORIES}
    for path in pages:
        for cat in CATEGORIES:
            key = str(cat["key"])
            if path.name.startswith(key):
                by_key[key].append(path)
                break

    with out.open("w") as f:
        f.write("# Efficiency-only category review pack\n\n")
        f.write("This pack intentionally crops every run to the RHS `Trigger / MBD` efficiency panel.\n\n")
        f.write(f"Input CSV: `{CSV_PATH}`\n\n")
        f.write(f"Overview: `{overview}`\n\n")
        f.write(f"Boundary tuning sheet: `{boundary}`\n\n")
        f.write("| Category | Runs | Share | Efficiency rule | Pages |\n")
        f.write("| --- | ---: | ---: | --- | --- |\n")
        for cat in CATEGORIES:
            key = str(cat["key"])
            page_links = ", ".join(f"`{p.name}`" for p in by_key[key])
            f.write(
                f"| {cat['label']} | {counts[key]} | {counts[key] / len(rows) * 100:.1f}% | "
                f"{cat['rule']} | {page_links} |\n"
            )
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    missing = [str(row["png_path"]) for row in rows if not Path(row["png_path"]).exists()]
    if missing:
        raise FileNotFoundError("Missing PNGs:\\n" + "\\n".join(missing[:20]))
    counts = Counter(str(row["category"]) for row in rows)
    overview = make_overview(rows, counts)
    boundary = make_boundary(rows, counts)
    pages = make_category_pages(rows, counts)
    index = write_index(rows, counts, overview, boundary, pages)
    print(f"overview={overview}")
    print(f"boundary={boundary}")
    print(f"pages={len(pages)}")
    print(f"index={index}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
