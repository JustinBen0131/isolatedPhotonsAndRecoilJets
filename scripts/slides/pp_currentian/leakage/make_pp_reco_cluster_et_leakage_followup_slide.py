#!/usr/bin/env python3
"""Build a slide-4-style pp reco-cluster ET leakage follow-up slide."""

from __future__ import annotations

import csv
import json
import math
import argparse
from io import BytesIO
from pathlib import Path
from textwrap import wrap

import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[4]
DEFAULT_BASE = REPO / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_recoEt50TruthWindow_20260605"
DEFAULT_CSV_PATH = DEFAULT_BASE / "validation/reco_cluster_et_leakage/pp_reco_cluster_et_leakage_components.csv"
DEFAULT_SUMMARY_PATH = DEFAULT_BASE / "validation/reco_cluster_et_leakage/pp_reco_cluster_et_leakage_summary.json"
DEFAULT_OUT = REPO / "dataOutput/slides/wp_gammajets_6_1_26/pp_reco_cluster_et_leakage_20260605_truthwindow"

W, H = 2560, 1440
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
TIMES_BOLD_ITALIC = FONT_DIR / "Times New Roman Bold Italic.ttf"

INK = (18, 24, 38)
MUTED = (78, 86, 101)
LINE = (208, 216, 228)
JET8 = (128, 128, 128)
JET12 = (33, 76, 195)
JET20 = (238, 114, 18)
JET30 = (204, 43, 168)
JET40 = (35, 150, 55)
BACKGROUND_ACCENT = (92, 78, 152)
BACKGROUND_SOFT = (244, 241, 250)
BACKGROUND_LINE = (204, 198, 230)
SIGNAL_ACCENT = (26, 111, 116)
SIGNAL_SOFT = (235, 248, 248)
SIGNAL_LINE = (174, 214, 216)
RESULT_ACCENT = (42, 91, 117)
RESULT_SOFT = (241, 246, 249)

SAMPLES = ("jet8", "jet12", "jet20", "jet30", "jet40")
COLORS = {
    "jet8": JET8,
    "jet12": JET12,
    "jet20": JET20,
    "jet30": JET30,
    "jet40": JET40,
    "Sum": INK,
}
MARKERS = {"jet8": "o", "jet12": "s", "jet20": "^", "jet30": "v", "jet40": "D", "Sum": "D"}
WINDOWS = {
    "jet8": "<14",
    "jet12": "14-21",
    "jet20": "21-32",
    "jet30": "32-42",
    "jet40": ">=42",
}


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    path = TIMES
    if bold and italic:
        path = TIMES_BOLD_ITALIC
    elif bold:
        path = TIMES_BOLD
    elif italic:
        path = TIMES_ITALIC
    try:
        return ImageFont.truetype(str(path), size)
    except OSError:
        return ImageFont.load_default()


F = {
    "title": font(58, bold=True),
    "subtitle": font(31),
    "panel_title": font(50, bold=True),
    "panel_subtitle": font(31),
    "card_title": font(43, bold=True),
    "card_body": font(35),
    "card_body_bold": font(35, bold=True),
    "body": font(30),
    "body_bold": font(30, bold=True),
    "small": font(26),
    "small_bold": font(26, bold=True),
}


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill, outline=LINE, width=2, radius=16) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], msg: str, fnt, fill=INK, anchor=None) -> None:
    draw.text(xy, msg, font=fnt, fill=fill, anchor=anchor)


def wrapped(draw: ImageDraw.ImageDraw, xy: tuple[int, int], msg: str, fnt, fill=INK, width_chars=44, spacing=8) -> int:
    x, y = xy
    for para in msg.split("\n"):
        for line in wrap(para, width=width_chars) or [""]:
            draw.text((x, y), line, font=fnt, fill=fill)
            y += fnt.size + spacing
    return y


def rich_line(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    parts: list[tuple[str, ImageFont.FreeTypeFont, tuple[int, int, int]]],
) -> int:
    cursor = x
    max_h = 0
    for msg, fnt, fill in parts:
        bbox = draw.textbbox((cursor, y), msg, font=fnt)
        draw.text((cursor, y), msg, font=fnt, fill=fill)
        cursor += bbox[2] - bbox[0]
        max_h = max(max_h, bbox[3] - bbox[1])
    return y + max_h


def paste_contained(canvas: Image.Image, img: Image.Image, box: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    bw, bh = x1 - x0, y1 - y0
    scale = min(bw / img.width, bh / img.height)
    nw, nh = int(img.width * scale), int(img.height * scale)
    resized = img.resize((nw, nh), Image.Resampling.LANCZOS)
    bg = Image.new("RGB", (bw, bh), "white")
    ox, oy = (bw - nw) // 2, (bh - nh) // 2
    bg.paste(resized, (ox, oy))
    canvas.paste(bg, (x0, y0))
    return (x0 + ox, y0 + oy, x0 + ox + nw, y0 + oy + nh)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV_PATH)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY_PATH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    return parser.parse_args()


def load_rows(csv_path: Path) -> list[dict]:
    rows: list[dict] = []
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            out = dict(row)
            for key in ("x_low", "x_high", "x_center", "x_err_low", "x_err_high", "weighted_entries", "weighted_error"):
                out[key] = float(out[key])
            rows.append(out)
    return rows


def pivot_rows(rows: list[dict]) -> list[dict]:
    by_bin: dict[tuple[float, float, float], dict] = {}
    for row in rows:
        key = (row["x_low"], row["x_high"], row["x_center"])
        target = by_bin.setdefault(
            key,
            {
                "x_low": row["x_low"],
                "x_high": row["x_high"],
                "x_center": row["x_center"],
                "x_err_low": row["x_err_low"],
                "x_err_high": row["x_err_high"],
            },
        )
        sample = row["sample"]
        target[f"weighted_entries_{sample}"] = row["weighted_entries"]
        target[f"weighted_error_{sample}"] = row["weighted_error"]

    out: list[dict] = []
    for _, row in sorted(by_bin.items()):
        sum_y = 0.0
        err2 = 0.0
        for sample in SAMPLES:
            y = float(row.get(f"weighted_entries_{sample}", 0.0))
            ey = float(row.get(f"weighted_error_{sample}", 0.0))
            row[f"weighted_entries_{sample}"] = y
            row[f"weighted_error_{sample}"] = ey
            sum_y += y
            err2 += ey * ey
        row["weighted_entries_Sum"] = sum_y
        row["weighted_error_Sum"] = math.sqrt(err2)
        for sample in SAMPLES:
            row[f"fraction_{sample}"] = row[f"weighted_entries_{sample}"] / sum_y if sum_y > 0 else 0.0
        out.append(row)
    return out


def pct(value: float) -> str:
    return "n/a" if not math.isfinite(value) else f"{100.0 * value:.1f}%"


def selection_label(summary: dict) -> str:
    candidate_filter = summary.get("candidate_filter") or "preselected"
    if candidate_filter == "abcd-reference":
        fixed = summary.get("abcd_fixed_iso_gev")
        try:
            fixed_s = f"{float(fixed):g}"
        except (TypeError, ValueError):
            fixed_s = "2"
        cap = " + caps" if summary.get("sample_cap_enforced") else ""
        return f"Truth-window + ABCD iso{fixed_s}{cap}"
    return str(summary.get("truth_window_selection") or "PPG12 truth-window pass")


def title_for_summary(summary: dict) -> str:
    if summary.get("candidate_filter") == "abcd-reference":
        if summary.get("sample_cap_enforced"):
            return "pp reco-cluster energy leakage check with ABCD rows and sample caps"
        return "pp reco-cluster energy leakage check with ABCD-matched rows"
    return "pp reco-cluster energy leakage check after truth-window rerun"


def readout_lines(summary: dict, max_selected_et: float) -> list[tuple[str, str]]:
    if summary.get("candidate_filter") == "abcd-reference":
        scope = (
            "matches ABCD rows with sample caps."
            if summary.get("sample_cap_enforced")
            else "matches the embedded ABCD rows."
        )
        return [
            ("Complete: ", "all five samples have 400 ROOTs."),
            ("Scope: ", scope),
            ("Reach: ", f"selected support reaches {max_selected_et:.2f} GeV."),
        ]
    return [
        ("Complete: ", "all five samples have 400 ROOTs."),
        ("Scope: ", "all truth-window preselected tree rows."),
        ("Reach: ", f"selected support reaches {max_selected_et:.2f} GeV."),
    ]


def range_fraction(pivot: list[dict], lo: float, hi: float, sample: str) -> float:
    total = 0.0
    sample_total = 0.0
    for row in pivot:
        if row["x_low"] >= lo - 1.0e-9 and row["x_high"] <= hi + 1.0e-9:
            total += row["weighted_entries_Sum"]
            sample_total += row[f"weighted_entries_{sample}"]
    return sample_total / total if total > 0 else float("nan")


def render_plot(pivot: list[dict]) -> Image.Image:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(12.0, 8.25), dpi=120, facecolor="white")
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[3.3, 1.15],
        left=0.105,
        right=0.965,
        bottom=0.115,
        top=0.885,
        hspace=0.055,
    )
    ax = fig.add_subplot(gs[0])
    fax = fig.add_subplot(gs[1], sharex=ax)

    display = [row for row in pivot if 12.0 <= row["x_center"] <= 50.0]
    xs = np.array([r["x_center"] for r in display], dtype=float)
    xerr = np.array([[r["x_err_low"] for r in display], [r["x_err_high"] for r in display]], dtype=float)

    for sample in SAMPLES + ("Sum",):
        y = np.array([r[f"weighted_entries_{sample}"] for r in display], dtype=float)
        ey = np.array([r[f"weighted_error_{sample}"] for r in display], dtype=float)
        mask = y > 0.0
        color = np.array(COLORS[sample]) / 255.0
        label = sample if sample != "Sum" else "Weighted sum"
        ax.errorbar(
            xs[mask],
            y[mask],
            yerr=ey[mask],
            xerr=xerr[:, mask],
            fmt=MARKERS[sample],
            markersize=5.0 if sample != "Sum" else 5.8,
            markerfacecolor=color if sample != "Sum" else "white",
            markeredgecolor=color,
            markeredgewidth=1.0,
            linestyle="none",
            color=color,
            ecolor=color,
            elinewidth=0.95,
            capsize=1.6,
            label=label,
            zorder=5 if sample == "Sum" else 4,
        )

    for sample in SAMPLES:
        color = np.array(COLORS[sample]) / 255.0
        vals = np.array([r[f"fraction_{sample}"] for r in display], dtype=float)
        sum_vals = np.array([r["weighted_entries_Sum"] for r in display], dtype=float)
        mask = sum_vals > 0.0
        fax.plot(xs[mask], vals[mask], marker=MARKERS[sample], markersize=5.2, linestyle="none", color=color, label=sample, zorder=3)

    ax.set_yscale("log")
    ax.set_xlim(12.0, 50.0)
    positive = [r[f"weighted_entries_{sample}"] for r in display for sample in (*SAMPLES, "Sum") if r[f"weighted_entries_{sample}"] > 0.0]
    ax.set_ylim(max(1.0e-5, min(positive) * 0.45), max(positive) * 2.4)
    ax.set_ylabel("Weighted entries / 1 GeV bin", fontsize=15.0)
    ax.tick_params(labelbottom=False, labelsize=12.5, length=6)
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(which="major", color="0.84", linestyle=":", linewidth=0.75)
    ax.grid(which="minor", color="0.91", linestyle=":", linewidth=0.45)

    fax.set_ylim(-0.035, 1.05)
    fax.set_yticks([0.0, 0.25, 0.50, 0.75, 1.0])
    fax.set_ylabel("Fraction of\nweighted sum", fontsize=13.5)
    fax.set_xlabel(r"Reco photon-cluster $E_T$ ($p_T^\gamma$) [GeV]", fontsize=15.0, loc="right")
    fax.tick_params(labelsize=12.5, length=6)
    fax.grid(which="major", color="0.84", linestyle=":", linewidth=0.75)

    ax.text(
        0.965,
        0.980,
        r"$\it{\bf{sPHENIX}}$ Internal" + "\nPYTHIA8 pp inclusive jet, 200 GeV",
        transform=ax.transAxes,
        fontsize=14.0,
        ha="right",
        va="top",
        linespacing=1.05,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 2.0},
        zorder=20,
    )
    ax.text(
        0.040,
        0.105,
        "jet8: <14 GeV      jet12: 14-21 GeV      jet20: 21-32 GeV\njet30: 32-42 GeV      jet40: >=42 GeV",
        transform=ax.transAxes,
        fontsize=15.8,
        ha="left",
        va="bottom",
        color=np.array(INK) / 255.0,
        linespacing=1.18,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "0.82", "alpha": 0.90, "linewidth": 0.8},
        zorder=20,
    )
    legend = ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.958, 0.785),
        frameon=False,
        fontsize=13.0,
        handlelength=1.8,
        labelspacing=0.25,
        borderaxespad=0.2,
        markerscale=1.0,
    )
    for item in legend.get_texts():
        item.set_fontsize(13.0)

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120, facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def add_plot_panel(canvas: Image.Image, draw: ImageDraw.ImageDraw, plot: Image.Image) -> None:
    box = (72, 184, 1654, 1348)
    x0, y0, x1, y1 = box
    rounded(draw, box, (255, 255, 255), outline=BACKGROUND_LINE, width=3, radius=18)
    draw.rounded_rectangle((x0 + 24, y0 + 22, x0 + 37, y0 + 100), radius=6, fill=BACKGROUND_ACCENT)
    text(draw, (x0 + 68, y0 + 14), "pp inclusive-jet source composition", F["panel_title"], INK)
    text(draw, (x0 + 70, y0 + 72), "Weighted reco-cluster energy in 1 GeV bins, with sample fractions below", F["panel_subtitle"], MUTED)
    paste_contained(canvas, plot, (x0 + 50, y0 + 120, x1 - 50, y1 - 34))


def sample_chip(draw: ImageDraw.ImageDraw, x: int, y: int, name: str, window: str, color: tuple[int, int, int]) -> None:
    rounded(draw, (x, y, x + 304, y + 58), (255, 255, 255), outline=(224, 230, 239), radius=10)
    draw.ellipse((x + 20, y + 21, x + 38, y + 39), fill=color)
    text(draw, (x + 52, y + 4), name, F["small_bold"], color)
    text(draw, (x + 156, y + 7), window, F["small"], INK)


def add_info_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title_msg: str,
    accent: tuple[int, int, int],
    fill: tuple[int, int, int],
    outline: tuple[int, int, int],
    lines: list[tuple[str, str]],
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, fill, outline=outline, width=3, radius=18)
    draw.rounded_rectangle((x0 + 20, y0 + 24, x0 + 31, y1 - 24), radius=5, fill=accent)
    text(draw, (x0 + 54, y0 + 22), title_msg, F["card_title"], INK)
    y = y0 + 82
    for label, body in lines:
        y = rich_line(draw, x0 + 54, y, [(label, F["card_body_bold"], accent), (body, F["card_body"], MUTED)])
        y += 12


def add_fraction_card(draw: ImageDraw.ImageDraw, pivot: list[dict]) -> None:
    box = (1710, 800, 2488, 1098)
    x0, y0, x1, y1 = box
    rounded(draw, box, RESULT_SOFT, outline=(190, 208, 218), width=3, radius=18)
    draw.rounded_rectangle((x0 + 20, y0 + 24, x0 + 31, y1 - 24), radius=5, fill=RESULT_ACCENT)
    text(draw, (x0 + 54, y0 + 22), "Composition readout", F["card_title"], INK)
    entries = [
        ("12-15 GeV", [("jet8", range_fraction(pivot, 12, 15, "jet8"), JET8), ("jet12", range_fraction(pivot, 12, 15, "jet12"), JET12)]),
        ("15-23 GeV", [("jet12", range_fraction(pivot, 15, 23, "jet12"), JET12), ("jet20", range_fraction(pivot, 15, 23, "jet20"), JET20)]),
        ("23-35 GeV", [("jet20", range_fraction(pivot, 23, 35, "jet20"), JET20), ("jet30", range_fraction(pivot, 23, 35, "jet30"), JET30)]),
        ("35-40 GeV", [("jet30", range_fraction(pivot, 35, 40, "jet30"), JET30), ("jet40", range_fraction(pivot, 35, 40, "jet40"), JET40)]),
    ]
    y = y0 + 82
    for label, pieces in entries:
        text(draw, (x0 + 54, y), label, F["small_bold"], INK)
        cursor = x0 + 218
        for name, value, color in pieces:
            rounded(draw, (cursor, y - 4, cursor + 220, y + 40), (255, 255, 255), outline=(224, 230, 239), radius=9)
            draw.ellipse((cursor + 14, y + 10, cursor + 28, y + 24), fill=color)
            text(draw, (cursor + 38, y + 2), f"{name} {pct(value)}", F["small_bold"], color)
            cursor += 232
        y += 50


def build(csv_path: Path, summary_path: Path, outdir: Path) -> tuple[Path, Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(summary_path.read_text())
    rows = load_rows(csv_path)
    pivot = pivot_rows(rows)
    plot = render_plot(pivot)
    max_selected_et = max(
        float(sample.get("cluster_et_max") or 0.0)
        for sample in summary.get("samples", [])
    )

    plot_png = outdir / "standardized_pp_reco_cluster_et_leakage_plot.png"
    plot.save(plot_png)

    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    title = title_for_summary(summary)
    text(draw, (72, 46), title, F["title"], INK)
    subtitle = "The diagnostic asks which pp inclusive-jet source populates the reconstructed photon-cluster-energy tail after weighting."
    text(draw, (76, 126), subtitle, F["subtitle"], MUTED)

    add_plot_panel(canvas, draw, plot)
    add_info_card(
        draw,
        (1710, 184, 2488, 414),
        "What is plotted",
        BACKGROUND_ACCENT,
        BACKGROUND_SOFT,
        BACKGROUND_LINE,
        [
            ("Object: ", "Reco photon-cluster energy"),
            ("Binning: ", "1 GeV bins, 12-50 GeV axis"),
            ("Selection: ", selection_label(summary)),
        ],
    )
    add_info_card(
        draw,
        (1710, 444, 2488, 770),
        "Main readout",
        SIGNAL_ACCENT,
        SIGNAL_SOFT,
        SIGNAL_LINE,
        readout_lines(summary, max_selected_et),
    )
    add_fraction_card(draw, pivot)

    rounded(draw, (1710, 1128, 2488, 1398), (255, 255, 255), outline=BACKGROUND_LINE, width=3, radius=18)
    text(draw, (1744, 1150), "Truth-jet windows", F["card_title"], INK)
    sample_chip(draw, 1744, 1210, "jet8", "<14", JET8)
    sample_chip(draw, 2078, 1210, "jet12", "14-21", JET12)
    sample_chip(draw, 1744, 1272, "jet20", "21-32", JET20)
    sample_chip(draw, 2078, 1272, "jet30", "32-42", JET30)
    sample_chip(draw, 1744, 1334, "jet40", ">=42", JET40)

    png = outdir / "pp_reco_cluster_et_leakage_followup_slide.png"
    script = outdir / "pp_reco_cluster_et_leakage_followup_slide_script.md"
    manifest = outdir / "pp_reco_cluster_et_leakage_followup_slide_manifest.json"
    canvas.save(png)

    if summary.get("candidate_filter") == "abcd-reference":
        script_text = f"""# WP GammaJets pp Follow-up Script - Reco-cluster Energy Leakage

This is the pp version of the reco-cluster energy leakage diagnostic. It uses the same visual organization as the embedded slide: the plot on the left shows the weighted source composition, and the cards on the right summarize what the audience should take from it.

The object plotted here is reconstructed photon-cluster energy from the pp inclusive-jet Monte Carlo. The selection is {selection_label(summary)}. That matters because the embedded counterpart is built from ABCD histograms, so this version removes the earlier mismatch where pp used all preselected tree rows.

The points are weighted with the current pp inclusive cross sections, and the lower panel shows the fraction of the weighted sum coming from each source sample in each 1 GeV bin. The diagnostic now reaches {max_selected_et:.2f} GeV in selected cluster energy, so the 40-50 GeV region is populated in the same candidate scope used for the embedded comparison.

The main point is that this is now an apples-to-apples bookkeeping check against the embedded ABCD counterpart. It is not a new cross-section estimate or a new production pass; it is a corrected compact extraction from the recovered canonical-tree pp outputs.
"""
    else:
        script_text = f"""# WP GammaJets pp Follow-up Script - Reco-cluster Energy Leakage

This is the pp version of the reco-cluster energy leakage diagnostic. It uses the same visual organization as the embedded slide: the plot on the left shows the weighted source composition, and the cards on the right summarize what the audience should take from it.

The object plotted here is reconstructed photon-cluster energy from the pp inclusive-jet Monte Carlo. The selection is {selection_label(summary)}. The points are weighted with the current pp inclusive cross sections, and the lower panel shows the fraction of the weighted sum coming from each source sample in each 1 GeV bin.

The main point is that this preselected-row view is useful for checking the production output, but it is not the same candidate scope as the embedded ABCD histogram counterpart. The selected support reaches {max_selected_et:.2f} GeV in this extraction.

This slide is not a new cross-section estimate; it is a rendering of the extracted pp leakage CSV and summary sidecars.
"""
    script.write_text(script_text, encoding="utf-8")
    manifest.write_text(
        json.dumps(
            {
                "png": str(png),
                "script": str(script),
                "standardized_plot_png": str(plot_png),
                "google_slides_mutation": False,
                "working_deck_style_reference": {
                    "presentation_id": "1D_m5vPM3eJIBI0nMz9RxAzZcjr15KQaG2uLOTjRTXOY",
                    "slide_object_id": "g3e831441f4b_0_261",
                    "image_title_readback": "reco_cluster_et_leakage_followup_slide.png",
                },
                "pp_reference_slide": {
                    "presentation_id": "1GAoEcN9UGOkUVxTgT0m1vRhghl3O1jd28LwSHg9klBs",
                    "slide_object_id": "g3a89b5423f1_0_192",
                    "image_title_readback": "inclusive_jet1234_reco_cluster_et_leakage_slide7_followup_12to50.png",
                },
                "source_summary": str(summary_path),
                "source_csv": str(csv_path),
                "source_schema": summary.get("schema"),
                "truth_window_selection": summary.get("truth_window_selection"),
                "candidate_filter": summary.get("candidate_filter"),
                "candidate_filter_description": summary.get("candidate_filter_description"),
                "abcd_fixed_iso_gev": summary.get("abcd_fixed_iso_gev"),
                "sample_cap_enforced": summary.get("sample_cap_enforced"),
                "weight_formula": summary.get("weight_formula"),
                "max_truth_window_selected_cluster_et": max_selected_et,
                "qa": summary.get("qa"),
                "composition_ranges": [
                    {
                        "range": label,
                        "fractions": {sample: range_fraction(pivot, lo, hi, sample) for sample in SAMPLES},
                    }
                    for label, lo, hi in [
                        ("12-15", 12, 15),
                        ("15-23", 15, 23),
                        ("23-35", 23, 35),
                        ("35-40", 35, 40),
                        ("40-50", 40, 50),
                    ]
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return png, script, manifest


if __name__ == "__main__":
    args = parse_args()
    for path in build(args.csv, args.summary, args.outdir):
        print(path)
