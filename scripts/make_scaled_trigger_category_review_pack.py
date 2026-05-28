#!/usr/bin/env python3
"""Build visual review sheets for scaled-trigger run-shape categories.

The output is intentionally review-facing rather than slide-facing: one
overview image, one boundary/tuning image, and paginated contact sheets that
show every run in each current category.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau/scaledTriggerRunByRunQA/"
    "scaled_trigger_run_by_run_20260521_202504"
)
CSV_PATH = BASE / "run_shape_categories.csv"
PNG_DIR = BASE / "png_all"
OUT_DIR = BASE / "category_review"

FONT_REG = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"

W, H = 3840, 2160

CATEGORIES = [
    {
        "key": "clean_full_turn_on",
        "label": "Clean full turn-on",
        "short": "efficient",
        "rule": "Tail unity and low-E starts near zero.",
        "what_to_check": "Curves should begin near zero, turn on smoothly, and reach ~1 in the 15-20 GeV tail.",
        "color": (28, 111, 70),
        "fill": (229, 246, 236),
    },
    {
        "key": "missing_photon_trigger_data",
        "label": "Missing or sparse photon-trigger data",
        "short": "missing/sparse",
        "rule": "Photon-trigger spectra absent/nearly absent relative to MBD.",
        "what_to_check": "Look for empty blue/red spectra, near-zero turn-on points, or too few trigger entries to interpret.",
        "color": (153, 45, 45),
        "fill": (255, 246, 244),
    },
    {
        "key": "flat_unity_from_threshold",
        "label": "Flat unity from threshold",
        "short": "flat unity",
        "rule": "Trigger/MBD already ~1 at low E and stays flat.",
        "what_to_check": "This is not an efficiency turn-on: the ratio is already unity before the photon-trigger threshold.",
        "color": (171, 91, 23),
        "fill": (255, 249, 238),
    },
    {
        "key": "early_turn_on_outlier",
        "label": "Early turn-on outlier",
        "short": "early-on",
        "rule": "Tail unity but low-E ratios are already high.",
        "what_to_check": "71260-like behavior: curves eventually reach unity, but the low-energy start is not near zero.",
        "color": (184, 72, 42),
        "fill": (255, 242, 235),
    },
    {
        "key": "tail_inefficient",
        "label": "Tail inefficient",
        "short": "low tail",
        "rule": "15-20 GeV trigger/MBD ratio below 0.80.",
        "what_to_check": "Usually a real-looking turn-on that never reaches full high-E efficiency.",
        "color": (128, 57, 140),
        "fill": (247, 240, 250),
    },
    {
        "key": "tail_overscaled",
        "label": "Tail overscaled",
        "short": "high tail",
        "rule": "15-20 GeV trigger/MBD ratio above 1.25.",
        "what_to_check": "High-E points overshoot MBD; often visually clear in the RHS overlay tail.",
        "color": (40, 101, 154),
        "fill": (236, 246, 255),
    },
    {
        "key": "tail_intermediate_deviation",
        "label": "Intermediate tail deviation",
        "short": "8-25% off",
        "rule": "Tail outside [0.92, 1.15] but not severe.",
        "what_to_check": "Borderline between useful-but-not-full and severe tail failure.",
        "color": (81, 87, 99),
        "fill": (244, 245, 247),
    },
    {
        "key": "mild_tail_deviation",
        "label": "Mild tail deviation",
        "short": "near tail",
        "rule": "Tail is near unity but fails strict [0.97, 1.08].",
        "what_to_check": "These may be acceptable or not, depending on how strict we want the clean category.",
        "color": (81, 87, 99),
        "fill": (244, 245, 247),
    },
    {
        "key": "low_stat_tail",
        "label": "Low-stat tail",
        "short": "low stats",
        "rule": "MBD tail has fewer than 100 counts.",
        "what_to_check": "Tail behavior is too noisy to classify robustly; inspect separately from shape pathologies.",
        "color": (81, 87, 99),
        "fill": (244, 245, 247),
    },
]


def font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size=size)


F_TITLE = font(FONT_BOLD, 74)
F_SUB = font(FONT_REG, 36)
F_H1 = font(FONT_BOLD, 46)
F_H2 = font(FONT_BOLD, 34)
F_BODY = font(FONT_REG, 29)
F_BODY_B = font(FONT_BOLD, 29)
F_SMALL = font(FONT_REG, 24)
F_SMALL_B = font(FONT_BOLD, 25)
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
                "r10_total",
                "r12_total",
            ]:
                parsed[key] = int(row[key]) if key == "run" else float(row[key])
            parsed["png_path"] = BASE / str(row["png"])
            rows.append(parsed)
    return rows


def run_label(row: dict[str, object]) -> str:
    return (
        f"Run {row['run']}  tail={float(row['r10_tail']):.2f}/{float(row['r12_tail']):.2f}  "
        f"low={float(row['r10_low']):.3f}/{float(row['r12_low']):.3f}"
    )


def open_thumb(row: dict[str, object], max_size: tuple[int, int]) -> Image.Image:
    path = Path(row["png_path"])  # type: ignore[arg-type]
    im = Image.open(path).convert("RGB")
    im = im.crop((12, 8, im.width - 10, im.height - 6))
    im.thumbnail(max_size, Image.Resampling.LANCZOS)
    return im


def wrap(draw: ImageDraw.ImageDraw, text: str, x: int, y: int, width_chars: int, font_obj, fill=INK, line_h=34) -> int:
    import textwrap

    lines = textwrap.wrap(text, width=width_chars)
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * line_h), line, font=font_obj, fill=fill)
    return y + max(1, len(lines)) * line_h


def category_info(key: str) -> dict[str, object]:
    for category in CATEGORIES:
        if category["key"] == key:
            return category
    raise KeyError(key)


def sorted_rows(rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(rows, key=lambda row: (float(row["mbd_tail"]), int(row["run"])), reverse=True)


def draw_cell(
    img: Image.Image,
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
    draw.rectangle((x + 2, y + 2, x + cell_w - 2, y + 48), fill=fill)
    draw.text((x + 14, y + 11), run_label(row), font=F_SMALL_B, fill=color)
    thumb = open_thumb(row, (cell_w - 28, cell_h - 76))
    img.paste(thumb, (x + (cell_w - thumb.width) // 2, y + 58 + ((cell_h - 76) - thumb.height) // 2))


def make_overview(rows: list[dict[str, object]], counts: Counter[str]) -> Path:
    out = OUT_DIR / "00_category_overview.png"
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    d.text((86, 62), "Scaled-trigger run categories: current qualitative map", font=F_TITLE, fill=INK)
    d.text((88, 154), "Counts are mutually exclusive; each category shows two high-stat examples for fast visual feedback.", font=F_SUB, fill=MUTED)
    d.line((88, 224, W - 88, 224), fill=LINE, width=3)

    card_w, card_h = 1190, 560
    gap_x, gap_y = 56, 44
    top = 286
    for idx, cat in enumerate(CATEGORIES):
        x = 88 + (idx % 3) * (card_w + gap_x)
        y = top + (idx // 3) * (card_h + gap_y)
        key = str(cat["key"])
        color = cat["color"]  # type: ignore[assignment]
        fill = cat["fill"]  # type: ignore[assignment]
        examples = sorted_rows(row for row in rows if row["category"] == key)[:2]
        d.rounded_rectangle((x, y, x + card_w, y + card_h), radius=12, fill=fill, outline=color, width=3)
        d.text((x + 28, y + 22), f"{cat['label']}: {counts[key]}", font=F_H1, fill=color)
        y2 = wrap(d, str(cat["rule"]), x + 28, y + 84, 64, F_BODY, MUTED, 35)
        d.text((x + 28, y2 + 10), "Check:", font=F_SMALL_B, fill=INK)
        wrap(d, str(cat["what_to_check"]), x + 106, y2 + 10, 60, F_SMALL, MUTED, 29)
        for j, row in enumerate(examples):
            thumb = open_thumb(row, (500, 235))
            px = x + 28 + j * 560
            py = y + 236
            img.paste(thumb, (px, py))
            d.text((px, py + thumb.height + 10), run_label(row), font=F_TINY, fill=INK)
    img.save(out, quality=95)
    return out


def make_boundary_sheet(rows: list[dict[str, object]], counts: Counter[str]) -> Path:
    out = OUT_DIR / "00_boundary_tuning_examples.png"
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    d.text((86, 62), "Boundary checks for category tuning", font=F_TITLE, fill=INK)
    d.text((88, 154), "These are the places where human feedback is most likely to change a threshold or split/merge a category.", font=F_SUB, fill=MUTED)
    d.line((88, 224, W - 88, 224), fill=LINE, width=3)

    panels = [
        ("clean_full_turn_on", "High-stat clean examples"),
        ("mild_tail_deviation", "Near-unity but not strict clean"),
        ("tail_intermediate_deviation", "Moderate tail deviation"),
        ("tail_inefficient", "Severe inefficient tail"),
        ("tail_overscaled", "Severe overscaled tail"),
        ("missing_photon_trigger_data", "Missing/sparse trigger data"),
        ("flat_unity_from_threshold", "Flat unity pathology"),
        ("early_turn_on_outlier", "Early turn-on pathology"),
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
        examples = sorted_rows(row for row in rows if row["category"] == key)[:3]
        for j, row in enumerate(examples):
            thumb = open_thumb(row, (535, 260))
            px = x + 24 + j * 575
            py = y + 82
            img.paste(thumb, (px, py))
            d.text((px, py + thumb.height + 8), run_label(row), font=F_TINY, fill=INK)
    img.save(out, quality=95)
    return out


def make_category_pages(rows: list[dict[str, object]], counts: Counter[str]) -> list[Path]:
    paths: list[Path] = []
    per_page = 16
    cols, rows_per_page = 4, 4
    left, top = 70, 312
    gap_x, gap_y = 28, 28
    cell_w = (W - 2 * left - (cols - 1) * gap_x) // cols
    cell_h = (H - top - 82 - (rows_per_page - 1) * gap_y) // rows_per_page

    for cat in CATEGORIES:
        key = str(cat["key"])
        color = cat["color"]  # type: ignore[assignment]
        fill = cat["fill"]  # type: ignore[assignment]
        cat_rows = sorted_rows(row for row in rows if row["category"] == key)
        total_pages = max(1, math.ceil(len(cat_rows) / per_page))
        slug = key
        for page in range(total_pages):
            out = OUT_DIR / f"{slug}_page_{page + 1:02d}_of_{total_pages:02d}.png"
            img = Image.new("RGB", (W, H), "white")
            d = ImageDraw.Draw(img)
            d.text((86, 58), f"{cat['label']} ({counts[key]} runs)", font=F_TITLE, fill=color)
            d.text((88, 146), f"Rule: {cat['rule']}  Page {page + 1} of {total_pages}.", font=F_SUB, fill=MUTED)
            d.line((88, 224, W - 88, 224), fill=LINE, width=3)
            wrap(d, f"What to check: {cat['what_to_check']}", 88, 246, 150, F_BODY, INK, 34)
            page_rows = cat_rows[page * per_page : (page + 1) * per_page]
            for idx, row in enumerate(page_rows):
                col = idx % cols
                rpage = idx // cols
                x = left + col * (cell_w + gap_x)
                y = top + rpage * (cell_h + gap_y)
                draw_cell(img, d, row, x, y, cell_w, cell_h, color, fill)
            img.save(out, quality=95)
            paths.append(out)
    return paths


def write_index(rows: list[dict[str, object]], counts: Counter[str], overview: Path, boundary: Path, pages: list[Path]) -> Path:
    out = OUT_DIR / "README.md"
    page_by_key: dict[str, list[Path]] = {str(cat["key"]): [] for cat in CATEGORIES}
    for path in pages:
        for cat in CATEGORIES:
            if path.name.startswith(str(cat["key"])):
                page_by_key[str(cat["key"])].append(path)
                break

    with out.open("w") as f:
        f.write("# Scaled-trigger category review pack\\n\\n")
        f.write(f"Input CSV: `{CSV_PATH}`\\n\\n")
        f.write(f"Overview: `{overview}`\\n\\n")
        f.write(f"Boundary tuning sheet: `{boundary}`\\n\\n")
        f.write("| Category | Runs | Share | Rule | Pages |\\n")
        f.write("| --- | ---: | ---: | --- | --- |\\n")
        total = len(rows)
        for cat in CATEGORIES:
            key = str(cat["key"])
            page_links = ", ".join(f"`{p.name}`" for p in page_by_key[key])
            f.write(f"| {cat['label']} | {counts[key]} | {counts[key] / total * 100:.1f}% | {cat['rule']} | {page_links} |\\n")
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    missing_pngs = [str(row["png_path"]) for row in rows if not Path(row["png_path"]).exists()]
    if missing_pngs:
        raise FileNotFoundError("Missing PNGs:\\n" + "\\n".join(missing_pngs[:20]))
    counts = Counter(str(row["category"]) for row in rows)

    overview = make_overview(rows, counts)
    boundary = make_boundary_sheet(rows, counts)
    pages = make_category_pages(rows, counts)
    index = write_index(rows, counts, overview, boundary, pages)

    print(f"overview={overview}")
    print(f"boundary={boundary}")
    print(f"pages={len(pages)}")
    print(f"index={index}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
