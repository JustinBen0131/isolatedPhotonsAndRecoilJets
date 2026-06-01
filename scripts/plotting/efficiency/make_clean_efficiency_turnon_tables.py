#!/usr/bin/env python3
"""Make turn-on efficiency contact sheets for clean scaled-trigger runs."""

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
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau/scaledTriggerRunByRunQA/"
    "scaled_trigger_run_by_run_20260521_202504"
)
CSV_PATH = BASE / "run_shape_categories.csv"
OUT_DIR = BASE / "clean_efficiency_turnon_tables"

FONT_REG = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"


def font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size=size)


# Large enough that 6x6 remains readable in Preview/Slides.
W, H = 4800, 3600
F_TITLE = font(FONT_BOLD, 76)
F_SUB = font(FONT_REG, 38)
F_RUN = font(FONT_BOLD, 34)
F_METRIC = font(FONT_REG, 25)
F_SMALL = font(FONT_REG, 24)

W_SLIDE, H_SLIDE = 7680, 4320
F_SLIDE_TITLE = font(FONT_BOLD, 122)
F_SLIDE_SUB = font(FONT_REG, 58)
F_SLIDE_RUN = font(FONT_BOLD, 70)
F_SLIDE_METRIC = font(FONT_REG, 42)
F_SLIDE_SMALL = font(FONT_REG, 38)

INK = (20, 26, 33)
MUTED = (83, 91, 100)
GREEN = (28, 111, 70)
FILL = (229, 246, 236)
LINE = (190, 198, 207)
LIGHT = (245, 248, 246)


def load_clean_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with CSV_PATH.open() as f:
        for row in csv.DictReader(f):
            if row["category"] != "clean_full_turn_on":
                continue
            parsed: dict[str, object] = dict(row)
            for key in [
                "run",
                "mbd_tail",
                "r10_low",
                "r12_low",
                "r10_mid",
                "r12_mid",
                "r10_tail",
                "r12_tail",
            ]:
                parsed[key] = int(row[key]) if key == "run" else float(row[key])
            parsed["png_path"] = BASE / row["png"]
            rows.append(parsed)
    return sorted(rows, key=lambda r: int(r["run"]))


def crop_efficiency(row: dict[str, object], max_size: tuple[int, int]) -> Image.Image:
    im = Image.open(Path(row["png_path"])).convert("RGB")  # type: ignore[arg-type]
    # RHS panel crop. Keep title, axes, unity line, legend, and tail labels.
    crop = im.crop((int(im.width * 0.50), 18, im.width - 8, im.height - 8))
    crop.thumbnail(max_size, Image.Resampling.LANCZOS)
    return crop


def crop_efficiency_tight(row: dict[str, object], max_size: tuple[int, int]) -> Image.Image:
    im = Image.open(Path(row["png_path"])).convert("RGB")  # type: ignore[arg-type]
    # Tight RHS plot crop for dense one-slide contact sheets. The run header
    # carries the run/tail info, so keep the axes and curves dominant.
    crop = im.crop((int(im.width * 0.525), 106, im.width - 14, im.height - 78))
    crop.thumbnail(max_size, Image.Resampling.LANCZOS)
    return crop


def draw_cell(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    row: dict[str, object],
    x: int,
    y: int,
    cell_w: int,
    cell_h: int,
) -> None:
    draw.rounded_rectangle((x, y, x + cell_w, y + cell_h), radius=10, outline=LINE, width=3, fill=(253, 254, 253))
    draw.rectangle((x + 2, y + 2, x + cell_w - 2, y + 76), fill=FILL)
    run = int(row["run"])
    draw.text((x + 16, y + 12), f"Run {run}", font=F_RUN, fill=GREEN)
    draw.text(
        (x + 16, y + 47),
        (
            f"low {float(row['r10_low']):.3f}/{float(row['r12_low']):.3f}   "
            f"mid {float(row['r10_mid']):.2f}/{float(row['r12_mid']):.2f}   "
            f"tail {float(row['r10_tail']):.3f}/{float(row['r12_tail']):.3f}"
        ),
        font=F_METRIC,
        fill=INK,
    )
    thumb = crop_efficiency(row, (cell_w - 30, cell_h - 98))
    canvas.paste(thumb, (x + (cell_w - thumb.width) // 2, y + 88 + ((cell_h - 102) - thumb.height) // 2))


def draw_slide_cell(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    row: dict[str, object],
    x: int,
    y: int,
    cell_w: int,
    cell_h: int,
) -> None:
    draw.rounded_rectangle((x, y, x + cell_w, y + cell_h), radius=8, outline=LINE, width=3, fill=(253, 254, 253))
    header_h = 74
    draw.rectangle((x + 2, y + 2, x + cell_w - 2, y + header_h), fill=FILL)
    run = int(row["run"])
    draw.text((x + 18, y + 5), f"{run}", font=F_SLIDE_RUN, fill=GREEN)
    draw.text(
        (x + 290, y + 22),
        f"tail {float(row['r10_tail']):.3f}/{float(row['r12_tail']):.3f}",
        font=F_SLIDE_METRIC,
        fill=INK,
    )
    thumb = crop_efficiency_tight(row, (cell_w - 28, cell_h - header_h - 18))
    canvas.paste(thumb, (x + (cell_w - thumb.width) // 2, y + header_h + 8 + ((cell_h - header_h - 18) - thumb.height) // 2))


def draw_slide_summary_cell(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    cell_w: int,
    cell_h: int,
    total_clean: int,
) -> None:
    draw.rounded_rectangle((x, y, x + cell_w, y + cell_h), radius=8, outline=(156, 190, 168), width=4, fill=(235, 247, 240))
    draw.text((x + 28, y + 34), "Clean group", font=F_SLIDE_RUN, fill=GREEN)
    draw.text((x + 28, y + 128), f"{total_clean}/620 runs", font=font(FONT_BOLD, 84), fill=INK)
    draw.text((x + 28, y + 222), f"{total_clean / 620 * 100:.1f}% of sample", font=F_SLIDE_SUB, fill=INK)
    draw.text((x + 28, y + 315), "All panels are RHS", font=F_SLIDE_SMALL, fill=MUTED)
    draw.text((x + 28, y + 365), "Trigger/MBD turn-ons", font=F_SLIDE_SMALL, fill=MUTED)


def make_one_slide_8x8(rows: list[dict[str, object]]) -> Path:
    out = OUT_DIR / "clean_full_efficiency_turnons_8x8_one_slide.png"
    img = Image.new("RGB", (W_SLIDE, H_SLIDE), "white")
    draw = ImageDraw.Draw(img)

    draw.text((86, 54), f"Clean full-efficiency turn-on runs ({len(rows)} total)", font=F_SLIDE_TITLE, fill=INK)
    draw.text(
        (90, 184),
        "8x8 one-slide contact sheet; runs sorted increasing. Cell header gives Photon10/Photon12 tail ratios over MBD.",
        font=F_SLIDE_SUB,
        fill=MUTED,
    )
    draw.line((86, 278, W_SLIDE - 86, 278), fill=LINE, width=5)

    cols, grid_rows = 8, 8
    left, top = 72, 330
    gap_x, gap_y = 12, 12
    bottom = 72
    cell_w = (W_SLIDE - 2 * left - (cols - 1) * gap_x) // cols
    cell_h = (H_SLIDE - top - bottom - (grid_rows - 1) * gap_y) // grid_rows

    for idx, row in enumerate(rows):
        col = idx % cols
        grid_row = idx // cols
        x = left + col * (cell_w + gap_x)
        y = top + grid_row * (cell_h + gap_y)
        draw_slide_cell(img, draw, row, x, y, cell_w, cell_h)

    summary_idx = len(rows)
    col = summary_idx % cols
    grid_row = summary_idx // cols
    x = left + col * (cell_w + gap_x)
    y = top + grid_row * (cell_h + gap_y)
    draw_slide_summary_cell(draw, x, y, cell_w, cell_h, len(rows))

    img.save(out, quality=95)
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_clean_rows()
    missing = [str(row["png_path"]) for row in rows if not Path(row["png_path"]).exists()]
    if missing:
        raise FileNotFoundError("Missing PNGs:\n" + "\n".join(missing))

    per_page = 36
    pages = math.ceil(len(rows) / per_page)
    cols, grid_rows = 6, 6
    left, top = 58, 292
    gap_x, gap_y = 18, 18
    cell_w = (W - 2 * left - (cols - 1) * gap_x) // cols
    cell_h = (H - top - 84 - (grid_rows - 1) * gap_y) // grid_rows

    outputs: list[Path] = []
    slide_output = make_one_slide_8x8(rows)
    outputs.append(slide_output)
    for page in range(pages):
        start = page * per_page
        page_rows = rows[start : start + per_page]
        out = OUT_DIR / f"clean_full_efficiency_turnons_6x6_page_{page + 1:02d}_of_{pages:02d}.png"
        img = Image.new("RGB", (W, H), "white")
        draw = ImageDraw.Draw(img)
        run_span = f"{int(page_rows[0]['run'])}-{int(page_rows[-1]['run'])}" if page_rows else ""
        draw.text((72, 54), f"Clean full-efficiency turn-on runs ({len(rows)} total)", font=F_TITLE, fill=INK)
        draw.text(
            (74, 144),
            f"Page {page + 1} of {pages}; runs sorted increasing ({run_span}). Each cell shows RHS Trigger/MBD efficiency only.",
            font=F_SUB,
            fill=MUTED,
        )
        draw.line((72, 222, W - 72, 222), fill=LINE, width=3)
        draw.text(
            (74, 238),
            "Cell header: low=1-3 GeV, mid=6-9 GeV, tail=15-20 GeV ratios for Photon10/Photon12 over MBD.",
            font=F_SMALL,
            fill=MUTED,
        )

        for idx, row in enumerate(page_rows):
            col = idx % cols
            grid_row = idx // cols
            x = left + col * (cell_w + gap_x)
            y = top + grid_row * (cell_h + gap_y)
            draw_cell(img, draw, row, x, y, cell_w, cell_h)
        img.save(out, quality=95)
        outputs.append(out)

    index = OUT_DIR / "README.md"
    with index.open("w") as f:
        f.write("# Clean full-efficiency turn-on tables\n\n")
        f.write(f"Input CSV: `{CSV_PATH}`\n\n")
        f.write("Runs are sorted in increasing run number. Each PNG shows RHS `Trigger / MBD` efficiency panels only.\n\n")
        f.write(f"- `{slide_output.name}`: one-slide 8x8 contact sheet with the final cell used for summary statistics.\n")
        for out in outputs:
            if out == slide_output:
                continue
            f.write(f"- `{out.name}`\n")
    print(index)
    for out in outputs:
        print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
