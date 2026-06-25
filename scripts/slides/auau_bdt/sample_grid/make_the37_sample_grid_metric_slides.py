#!/usr/bin/env python3
"""Render THE-37/THE-43 sample-grid and depth-sweep slide candidates."""

from __future__ import annotations

import argparse
import json
import math
import textwrap
import sys as _codex_sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle

_THIS_FILE = Path(__file__).resolve()
_SCRIPTS_DIR = next((p for p in _THIS_FILE.parents if p.name == "scripts"), _THIS_FILE.parent)
if str(_SCRIPTS_DIR) not in _codex_sys.path:
    _codex_sys.path.append(str(_SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


DEFAULT_OUTDIR = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/"
    "auauTightBDTValidation/THE37_sample_grid_slides_20260621"
)
DEFAULT_SNAPSHOT = DEFAULT_OUTDIR / "the37_sample_grid_registry_snapshot_20260621.json"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
FONT_FAMILY = "Times New Roman"

COLORS = {
    "ink": "#111827",
    "muted": "#5b6472",
    "line": "#cbd5e1",
    "panel": "#f8fafc",
    "signal": "#7c3aed",
    "background": "#d97706",
    "good": "#e8f5ee",
    "watch": "#fff1d6",
    "risk": "#fde8e8",
    "blue": "#2563ad",
    "orange": "#d97706",
    "red": "#c94f43",
}


def configure_fonts() -> None:
    for font_path in TIMES_FONTS:
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
    plt.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.serif": [FONT_FAMILY],
            "mathtext.fontset": "custom",
            "mathtext.rm": FONT_FAMILY,
            "mathtext.it": f"{FONT_FAMILY}:italic",
            "mathtext.bf": f"{FONT_FAMILY}:bold",
            "axes.unicode_minus": False,
        }
    )


def fmt(value: float | int | None, digits: int = 3, missing: str = "--") -> str:
    if value is None:
        return missing
    try:
        value = float(value)
    except (TypeError, ValueError):
        return missing
    if not math.isfinite(value):
        return missing
    return f"{value:.{digits}f}"


def fmt_rows(value: float | int | None) -> str:
    if value is None:
        return "--"
    value = int(value)
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.0f}k"
    return str(value)


def add_round_box(fig, x: float, y: float, w: float, h: float, *, face: str, edge: str, lw: float = 1.2) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle="round,pad=0.008,rounding_size=0.012",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=-10,
        )
    )


def gap_face(auc_gap: float | None) -> str:
    if auc_gap is None:
        return COLORS["panel"]
    a = abs(float(auc_gap))
    if a < 0.01:
        return COLORS["good"]
    if a < 0.04:
        return COLORS["watch"]
    return COLORS["risk"]


def draw_header(fig, title: str, kicker: str, *, title_size: float = 31.0) -> None:
    fig.text(0.045, 0.94, title, ha="left", va="top", fontsize=title_size, fontweight="bold", color=COLORS["ink"])
    fig.text(0.047, 0.888, kicker, ha="left", va="top", fontsize=16.5, color=COLORS["muted"])


def draw_metric_cell(fig, item: dict, x: float, y: float, w: float, h: float, *, title: str, compact: bool = False) -> None:
    auc_gap = item.get("auc_gap")
    add_round_box(fig, x, y, w, h, face=gap_face(auc_gap), edge="#c7d2e2")
    title_size = 15.2 if compact else 16.5
    if title:
        fig.text(x + 0.018, y + h - 0.026, title, ha="left", va="top", fontsize=title_size, fontweight="bold", color=COLORS["ink"])
    fig.text(
        x + w - 0.018,
        y + h - 0.026,
        f"N={fmt_rows(item.get('n_rows'))}",
        ha="right",
        va="top",
        fontsize=12.5,
        color=COLORS["muted"],
    )

    rows = [
        ("Holdout AUC", fmt(item.get("holdout_auc"), 3)),
        ("Holdout logloss", fmt(item.get("holdout_logloss"), 3)),
        ("AUC gap", fmt(item.get("auc_gap"), 3)),
        ("Logloss gap", fmt(item.get("logloss_gap"), 3)),
    ]
    left = x + (0.018 if compact else 0.022)
    val_x = x + w - (0.020 if compact else 0.026)
    y0 = y + h - (0.058 if compact else 0.070)
    dy = h * (0.172 if compact else 0.150)
    for idx, (label, value) in enumerate(rows):
        yy = y0 - idx * dy
        color = COLORS["red"] if "gap" in label.lower() and value not in {"--", "0.000"} else COLORS["ink"]
        fig.text(left, yy, label, ha="left", va="center", fontsize=11.8 if compact else 13.2, color=COLORS["muted"])
        fig.text(val_x, yy, value, ha="right", va="center", fontsize=13.4 if compact else 15.0, fontweight="bold", color=color)

    if not compact:
        fig.text(
            left,
            y + 0.020,
            f"train AUC {fmt(item.get('train_auc'), 3)}   train logloss {fmt(item.get('train_logloss'), 3)}",
            ha="left",
            va="bottom",
            fontsize=10.8,
            color=COLORS["muted"],
        )


def summarize_grid(items: list[dict]) -> str:
    best = max(items, key=lambda x: x.get("holdout_auc") if x.get("holdout_auc") is not None else -999)
    calm = min(items, key=lambda x: abs(x.get("auc_gap") or 999))
    risk = max(items, key=lambda x: abs(x.get("auc_gap") or 0))
    return (
        f"Best holdout AUC: {best['signal_label']} / {best['background_label']} = {fmt(best.get('holdout_auc'), 3)}.  "
        f"Smallest |AUC gap|: {calm['signal_label']} / {calm['background_label']} = {fmt(calm.get('auc_gap'), 3)}.  "
        f"Largest |AUC gap|: {risk['signal_label']} / {risk['background_label']} = {fmt(risk.get('auc_gap'), 3)}."
    )


def readout_items(items: list[dict]) -> list[tuple[str, str]]:
    best = max(items, key=lambda x: x.get("holdout_auc") if x.get("holdout_auc") is not None else -999)
    calm = min(items, key=lambda x: abs(x.get("auc_gap") or 999))
    risk = max(items, key=lambda x: abs(x.get("auc_gap") or 0))
    return [
        ("Best holdout AUC", f"{best['signal_label']} / {best['background_label']} = {fmt(best.get('holdout_auc'), 3)}"),
        ("Smallest AUC gap", f"{calm['signal_label']} / {calm['background_label']} = {fmt(calm.get('auc_gap'), 3)}"),
        ("Largest AUC gap", f"{risk['signal_label']} / {risk['background_label']} = {fmt(risk.get('auc_gap'), 3)}"),
    ]


def draw_readout_cards(fig, items: list[dict]) -> None:
    labels = readout_items(items)
    xs = [0.060, 0.365, 0.670]
    colors = [("#edf2fb", "#9bb7df"), ("#e8f5ee", "#98c7aa"), ("#fff2d6", "#dfb85d")]
    for x, (head, body), (face, edge) in zip(xs, labels, colors):
        add_round_box(fig, x, 0.060, 0.270, 0.095, face=face, edge=edge)
        fig.text(x + 0.018, 0.122, head, ha="left", va="center", fontsize=16.0, fontweight="bold", color=COLORS["ink"])
        wrapped = "\n".join(textwrap.wrap(body, width=35))
        fig.text(x + 0.018, 0.087, wrapped, ha="left", va="center", fontsize=11.5, color=COLORS["ink"], linespacing=0.9)


def render_grid_slide(items: list[dict], *, system: str, outdir: Path) -> Path:
    configure_fonts()
    if system == "AuAu":
        row_order = ["12", "20"]
        col_order = ["12_20", "12_20_30", "12_20_30_40"]
        title = "AuAu sample-grid controls: signal and inclusive-jet coverage"
        kicker = "Rows vary the embedded photon signal sample; columns vary the embedded inclusive-jet background ladder."
        output = outdir / "the37_auau_sample_grid_metric_comparison_slide.png"
        row_title = "Signal training sample"
        col_title = "Background training sample"
    else:
        row_order = ["5", "5_10", "5_10_20"]
        col_order = ["8_12", "8_12_20", "8_12_20_30", "8_12_20_30_40"]
        title = "pp sample-grid controls: currentIAN baseV3E"
        kicker = "Rows vary the photon+jet signal ladder; columns vary the inclusive-jet background ladder."
        output = outdir / "the37_pp_sample_grid_metric_comparison_slide.png"
        row_title = "Signal training sample"
        col_title = "Background training sample"

    by_key = {(str(i["tag_signal"]), str(i["tag_background"])): i for i in items}
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    draw_header(fig, title, kicker, title_size=30.0 if system == "AuAu" else 29.0)

    grid_left = 0.165
    grid_right = 0.955
    grid_bottom = 0.190
    grid_top = 0.770
    gap_x = 0.018
    gap_y = 0.030
    n_rows = len(row_order)
    n_cols = len(col_order)
    cell_w = (grid_right - grid_left - (n_cols - 1) * gap_x) / n_cols
    cell_h = (grid_top - grid_bottom - (n_rows - 1) * gap_y) / n_rows

    fig.text(grid_left, 0.812, col_title, ha="left", va="center", fontsize=14.5, fontweight="bold", color=COLORS["muted"])
    fig.text(0.045, grid_top + 0.010, row_title, ha="left", va="bottom", fontsize=13.5, fontweight="bold", color=COLORS["muted"])
    for c, col in enumerate(col_order):
        x = grid_left + c * (cell_w + gap_x)
        label = by_key[(row_order[0], col)]["background_label"]
        label = label.replace("Jet", "Jet ")
        fig.text(x + cell_w / 2, 0.785, label, ha="center", va="center", fontsize=13.5 if system == "pp" else 15.5, fontweight="bold", color=COLORS["ink"])

    for r, row in enumerate(row_order):
        y = grid_top - (r + 1) * cell_h - r * gap_y
        label = by_key[(row, col_order[0])]["signal_label"]
        add_round_box(fig, 0.040, y, 0.105, cell_h, face="#f5f3ff", edge="#c4b5fd")
        label = label.replace("Photon", "Photon ")
        fig.text(0.092, y + cell_h / 2, label, ha="center", va="center", fontsize=13.5 if system == "pp" else 16.5, fontweight="bold", color=COLORS["ink"], wrap=True)
        for c, col in enumerate(col_order):
            x = grid_left + c * (cell_w + gap_x)
            title_text = "" if system == "pp" else by_key[(row, col)]["background_label"]
            draw_metric_cell(fig, by_key[(row, col)], x, y, cell_w, cell_h, title=title_text, compact=(system == "pp"))

    fig.text(0.955, 0.812, "Cell shade = |train-holdout AUC gap| warning", ha="right", va="center", fontsize=12.5, color=COLORS["muted"])
    draw_readout_cards(fig, items)
    fig.savefig(output, dpi=SLIDE_DPI)
    plt.close(fig)
    return output


def padded(values: np.ndarray, frac: float = 0.14, min_pad: float = 0.001) -> tuple[float, float]:
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    if math.isclose(lo, hi):
        return lo - min_pad, hi + min_pad
    pad = max(min_pad, (hi - lo) * frac)
    return lo - pad, hi + pad


def style_ax(ax) -> None:
    ax.grid(True, color="#e5e9f0", linewidth=1.0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#9aa7b8")
    ax.tick_params(labelsize=12.5, colors=COLORS["muted"])


def depth_regions(ax, depths: np.ndarray) -> None:
    ax.axvspan(float(depths.min()) - 0.45, 3.5, color="#edf6f0", zorder=0)
    ax.axvspan(3.5, 4.5, color="#edf2fb", zorder=0)
    ax.axvspan(4.5, float(depths.max()) + 0.45, color="#fff2d6", zorder=0)
    ax.axvline(4, color=COLORS["ink"], linewidth=1.6, linestyle=(0, (4, 3)), zorder=1)


def render_depth_slide(items: list[dict], outdir: Path) -> Path:
    configure_fonts()
    items = sorted(items, key=lambda i: i["depth"])
    depths = np.array([int(i["depth"]) for i in items], dtype=float)
    train_auc = np.array([float(i["train_auc"]) for i in items])
    holdout_auc = np.array([float(i["holdout_auc"]) for i in items])
    train_ll = np.array([float(i["train_logloss"]) for i in items])
    holdout_ll = np.array([float(i["holdout_logloss"]) for i in items])
    auc_gap = np.array([float(i["auc_gap"]) for i in items])
    ll_gap = np.array([float(i["logloss_gap"]) for i in items])

    output = outdir / "the43_depth_2to10_overfit_metric_sweep_slide.png"
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    draw_header(
        fig,
        "Tree-depth sweep: holdout gain versus capacity gap",
        "Corrected THE8 baseline, 90/10 row split. The overfit metric is the train-holdout gap.",
    )

    add_round_box(fig, 0.055, 0.570, 0.420, 0.270, face=COLORS["panel"], edge="#d3dce9")
    add_round_box(fig, 0.525, 0.570, 0.420, 0.270, face=COLORS["panel"], edge="#d3dce9")
    add_round_box(fig, 0.055, 0.245, 0.890, 0.235, face=COLORS["panel"], edge="#d3dce9")

    ax_auc = fig.add_axes([0.095, 0.610, 0.350, 0.175])
    ax_ll = fig.add_axes([0.565, 0.610, 0.350, 0.175])
    ax_gap = fig.add_axes([0.105, 0.292, 0.800, 0.128])
    for ax in (ax_auc, ax_ll, ax_gap):
        depth_regions(ax, depths)

    ax_auc.plot(depths, train_auc, color=COLORS["orange"], marker="o", linewidth=2.8, markersize=7.5, label="Train")
    ax_auc.plot(depths, holdout_auc, color=COLORS["blue"], marker="o", linewidth=2.8, markersize=7.5, label="Holdout")
    ax_auc.set_title("AUC: train versus holdout", fontsize=17.0, fontweight="bold")
    ax_auc.set_xlabel("XGBoost max depth", fontsize=13.5)
    ax_auc.set_ylabel("AUC", fontsize=13.5)
    ax_auc.set_xticks(depths)
    ax_auc.set_ylim(*padded(np.r_[train_auc, holdout_auc], min_pad=0.001))
    ax_auc.legend(frameon=False, fontsize=12.5, loc="lower right")
    style_ax(ax_auc)

    ax_ll.plot(depths, train_ll, color=COLORS["orange"], marker="o", linewidth=2.8, markersize=7.5, label="Train")
    ax_ll.plot(depths, holdout_ll, color=COLORS["blue"], marker="o", linewidth=2.8, markersize=7.5, label="Holdout")
    ax_ll.set_title("Logloss: calibration check", fontsize=17.0, fontweight="bold")
    ax_ll.set_xlabel("XGBoost max depth", fontsize=13.5)
    ax_ll.set_ylabel("Logloss", fontsize=13.5)
    ax_ll.set_xticks(depths)
    ax_ll.set_ylim(*padded(np.r_[train_ll, holdout_ll], min_pad=0.001))
    ax_ll.legend(frameon=False, fontsize=12.5, loc="upper right")
    style_ax(ax_ll)

    width = 0.34
    ax_gap.axhline(0, color="#94a3b8", linewidth=1.1)
    ax_gap.bar(depths - width / 2, auc_gap, width=width, color=COLORS["blue"], label="AUC gap")
    ax_gap.bar(depths + width / 2, ll_gap, width=width, color=COLORS["red"], label="Logloss gap")
    ax_gap.set_title("Train-holdout gaps are the capacity warning", fontsize=17.0, fontweight="bold")
    ax_gap.set_xlabel("XGBoost max depth", fontsize=13.5)
    ax_gap.set_ylabel("Gap size", fontsize=13.5)
    ax_gap.set_xticks(depths)
    ax_gap.legend(frameon=False, fontsize=12.5, loc="upper center", ncol=2, bbox_to_anchor=(0.73, 1.30))
    style_ax(ax_gap)
    lim = max(0.004, float(np.nanmax(np.abs(np.r_[auc_gap, ll_gap]))) * 1.35)
    ax_gap.set_ylim(-lim, lim)

    best = max(items, key=lambda i: i["holdout_auc"])
    calm = min(items, key=lambda i: abs(i["auc_gap"]))
    risk = max(items, key=lambda i: abs(i["auc_gap"]))
    cards = [
        ("Best holdout AUC", f"depth {best['depth']} = {fmt(best['holdout_auc'], 3)}", "#edf2fb", "#9bb7df"),
        ("Smallest AUC gap", f"depth {calm['depth']} = {fmt(calm['auc_gap'], 4)}", "#e8f5ee", "#98c7aa"),
        ("Largest AUC gap", f"depth {risk['depth']} = {fmt(risk['auc_gap'], 4)}", "#fff2d6", "#dfb85d"),
    ]
    x0s = [0.075, 0.365, 0.655]
    for (head, body, face, edge), x in zip(cards, x0s):
        add_round_box(fig, x, 0.070, 0.255, 0.105, face=face, edge=edge)
        fig.text(x + 0.020, 0.137, head, ha="left", va="center", fontsize=18.0, fontweight="bold", color=COLORS["ink"])
        fig.text(x + 0.020, 0.095, body, ha="left", va="center", fontsize=16.0, color=COLORS["ink"])
    fig.text(0.945, 0.034, "Depth 4 is the corrected baseline reference.", ha="right", va="bottom", fontsize=12.0, color=COLORS["muted"])
    fig.savefig(output, dpi=SLIDE_DPI)
    plt.close(fig)
    return output


def write_notes(outdir: Path, outputs: dict[str, Path], payload: dict) -> Path:
    note = outdir / "the37_the43_sample_grid_slide_speaker_notes.md"
    lines = [
        "# THE-37 / THE-43 sample-grid slide notes",
        "",
        "These are local PNG candidates only. No Google Slides deck was mutated.",
        "",
        "## How to say the sample-grid slides",
        "Rows change the signal training sample and columns change the inclusive-jet background training ladder. Each cell reports the exact registry train/holdout diagnostic for that trained lane. The color is deliberately not a physics result; it is a visual warning for the size of the train-holdout AUC gap.",
        "",
        "For pp, this isolates sample-composition behavior without the AuAu underlying event. For AuAu, it shows how the corrected truth-default training responds when we vary signal coverage and inclusive-jet background coverage.",
        "",
        "## How to say the depth slide",
        "The metric Justin cares about for overfitting is the train-holdout gap. Holdout AUC tells us whether the model generalizes better, but a growing train-holdout AUC or logloss gap is the capacity warning. Depth 4 is the corrected baseline reference.",
        "",
        "## Outputs",
    ]
    for key, path in outputs.items():
        lines.append(f"- {key}: `{path}`")
    lines.append("")
    lines.append(f"Snapshot: `{DEFAULT_SNAPSHOT}`")
    lines.append(f"Remote snapshot time: {payload.get('generated_remote_time')}")
    note.write_text("\n".join(lines) + "\n")
    return note


def render_all(snapshot: Path, outdir: Path) -> dict[str, Path]:
    payload = json.loads(snapshot.read_text())
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "auau_grid": render_grid_slide(payload["auau"], system="AuAu", outdir=outdir),
        "pp_grid": render_grid_slide(payload["pp"], system="pp", outdir=outdir),
        "depth_sweep": render_depth_slide(payload["depth"], outdir=outdir),
    }
    notes = write_notes(outdir, outputs, payload)
    manifest = outdir / "the37_the43_sample_grid_metric_slides_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE37_THE43_SAMPLE_GRID_SLIDE_MANIFEST_V1",
                "snapshot": str(snapshot),
                "outputs": {k: str(v) for k, v in outputs.items()},
                "speaker_notes": str(notes),
                "deck_mutated": False,
                "slide_numbers_baked_in": False,
                "font_family": FONT_FAMILY,
                "font_files": [str(p) for p in TIMES_FONTS if p.exists()],
                "size_px": [2560, 1440],
                "dpi": SLIDE_DPI,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    outputs["speaker_notes"] = notes
    outputs["manifest"] = manifest
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    outputs = render_all(args.snapshot, args.outdir)
    for path in outputs.values():
        print(path)


if __name__ == "__main__":
    main()
