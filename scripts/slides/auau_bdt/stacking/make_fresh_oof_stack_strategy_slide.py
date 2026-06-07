#!/usr/bin/env python3
"""Render audience-facing slides explaining the fresh OOF stack strategy."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


W, H = 2560, 1440
DPI = 200
OUTDIR = Path("dataOutput/auauTightBDTValidation/fresh_oof_stack_strategy_20260604")

COLORS = {
    "ink": "#171717",
    "muted": "#5f6368",
    "line": "#d6d9de",
    "soft_gray": "#eef1f4",
    "fold": "#eaf3ef",
    "fold_edge": "#9fc6b1",
    "heldout": "#f4dd9a",
    "heldout_edge": "#c89b18",
    "test": "#ded8ee",
    "test_edge": "#7a67a8",
    "stack": "#e8f0f7",
    "stack_edge": "#6e9abf",
    "green": "#2e7d5b",
    "purple": "#6b4ea3",
    "gold": "#b27b00",
    "dark_blue": "#2f5f83",
}


def setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def fig_ax():
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    return fig, ax


def add_text(ax, x, y, text, size=20, weight="normal", color=None, ha="left", va="top", linespacing=1.18):
    return ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=size,
        fontweight=weight,
        color=color or COLORS["ink"],
        ha=ha,
        va=va,
        linespacing=linespacing,
    )


def box(ax, x, y, w, h, fc, ec=None, lw=1.4, radius=0.012, zorder=1):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        transform=ax.transAxes,
        linewidth=lw,
        edgecolor=ec or fc,
        facecolor=fc,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, x1, y1, x2, y2, color=None, lw=2.4):
    patch = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        transform=ax.transAxes,
        arrowstyle="-|>",
        mutation_scale=18,
        linewidth=lw,
        color=color or COLORS["muted"],
        shrinkA=4,
        shrinkB=4,
        zorder=3,
    )
    ax.add_patch(patch)
    return patch


def save_png_pdf(fig, stem: str) -> tuple[Path, Path]:
    png = OUTDIR / f"{stem}.png"
    pdf = OUTDIR / f"{stem}.pdf"
    fig.savefig(png, dpi=DPI)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def draw_header(ax, title: str, subtitle: str) -> None:
    add_text(ax, 0.045, 0.935, title, 31, "bold")
    add_text(ax, 0.045, 0.884, subtitle, 18.5, color=COLORS["muted"])


def draw_event_key_band(ax, x: float, y: float, w: float, text: str) -> None:
    box(ax, x, y, w, 0.058, COLORS["soft_gray"], COLORS["line"], lw=1.0, radius=0.007)
    add_text(ax, x + 0.020, y + 0.037, text, 15.5, "bold", va="center", color=COLORS["ink"])


def draw_setup_band(ax) -> None:
    box(ax, 0.070, 0.735, 0.880, 0.095, COLORS["soft_gray"], COLORS["line"], lw=1.0, radius=0.007)
    add_text(ax, 0.095, 0.797, "Simple holdout stacking:", 18.5, "bold", va="center", color=COLORS["ink"])
    add_text(
        ax,
        0.350,
        0.797,
        "64% base training | 16% stack training | 20% locked test",
        17.0,
        color=COLORS["ink"],
        va="center",
    )
    add_text(
        ax,
        0.095,
        0.764,
        "Core rule:",
        15.8,
        "bold",
        color=COLORS["muted"],
        va="center",
    )
    add_text(
        ax,
        0.195,
        0.764,
        "the stacker trains only on base scores from events the base models never saw.",
        15.8,
        color=COLORS["muted"],
        va="center",
    )


def draw_flow_card(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    label: str,
    title: str,
    body: str,
    fc: str,
    ec: str,
    accent: str,
) -> None:
    box(ax, x, y, w, h, fc, ec, lw=1.5, radius=0.010)
    add_text(ax, x + 0.024, y + h - 0.037, label, 19, "bold", color=accent)
    add_text(ax, x + 0.024, y + h - 0.087, title, 22, "bold", color=COLORS["ink"])
    add_text(ax, x + 0.024, y + h - 0.150, body, 16.0, color=COLORS["ink"], linespacing=1.25)


def draw_callout(ax, x: float, y: float, w: float, h: float, title: str, body: str, fc: str, ec: str, accent: str) -> None:
    box(ax, x, y, w, h, fc, ec, lw=1.35, radius=0.010)
    add_text(ax, x + 0.025, y + h - 0.032, title, 21, "bold", color=accent)
    add_text(ax, x + 0.025, y + h - 0.080, body, 17.0, color=COLORS["ink"], linespacing=1.24)


def draw_simple_holdout_scaffold(ax) -> None:
    draw_header(
        ax,
        "The stacker needs honest base-model scores",
        "The stacker is another trained model, so its training rows must be scored by base models that did not see those events.",
    )
    draw_setup_band(ax)

    y, h = 0.395, 0.285
    w = 0.265
    x1, x2, x3 = 0.070, 0.370, 0.670
    draw_flow_card(
        ax,
        x1,
        y,
        w,
        h,
        "Set A: 64%",
        "Train base models",
        "BDT and MLP learn the\nsignal/background boundary\nusing Set A events only.",
        COLORS["fold"],
        COLORS["fold_edge"],
        COLORS["green"],
    )
    draw_flow_card(
        ax,
        x2,
        y,
        w,
        h,
        "Set B: 16%",
        "Train stacker",
        "Stacker trains on honest\nBDT + MLP scores.\nContext variant adds $E_T$\n(+ centrality for AuAu).",
        COLORS["stack"],
        COLORS["stack_edge"],
        COLORS["dark_blue"],
    )
    draw_flow_card(
        ax,
        x3,
        y,
        w,
        h,
        "Set C: 20%",
        "Locked test",
        "No training choices use\nthese events. They are\nscored only after models\nare fixed.",
        COLORS["test"],
        COLORS["test_edge"],
        COLORS["purple"],
    )
    arrow(ax, x1 + w + 0.010, y + h * 0.55, x2 - 0.012, y + h * 0.55, COLORS["muted"])
    arrow(ax, x2 + w + 0.010, y + h * 0.55, x3 - 0.012, y + h * 0.55, COLORS["muted"])

    draw_callout(
        ax,
        0.070,
        0.180,
        0.405,
        0.160,
        "Why this is valid",
        "Set-B rows are scored by models\ntrained only on Set A.",
        "#eef7f2",
        COLORS["fold_edge"],
        COLORS["green"],
    )
    draw_callout(
        ax,
        0.525,
        0.180,
        0.405,
        0.160,
        "Cost of the simple version",
        "Only 16% of events train the simple\nstacker, motivating the fold strategy.",
        "#fff2c2",
        "#d8bf5f",
        COLORS["gold"],
    )

    box(ax, 0.070, 0.050, 0.880, 0.105, "white", COLORS["line"], lw=1.2, radius=0.009)
    add_text(ax, 0.100, 0.108, "Next: out-of-fold (OOF)", 19, "bold", color=COLORS["dark_blue"], va="center")
    add_text(
        ax,
        0.345,
        0.108,
        "OOF means each event is scored by a model that did not train on it.\nRepeating that rule across folds gives honest stack rows for the full 80% development region.",
        14.6,
        color=COLORS["ink"],
        va="center",
        linespacing=1.20,
    )


def draw_oof_partition_bar(ax):
    x, y, w, h = 0.070, 0.715, 0.880, 0.090
    add_text(ax, x, y + 0.132, "Event-level split", 20.5, "bold")
    add_text(
        ax,
        x + 0.215,
        y + 0.132,
        "keep candidates from the same event in the same split",
        14.2,
        "bold",
        color=COLORS["muted"],
    )
    ax.add_patch(Rectangle((x, y), w * 0.80, h, transform=ax.transAxes, facecolor=COLORS["fold"], edgecolor=COLORS["fold_edge"], linewidth=1.5))
    ax.add_patch(Rectangle((x + w * 0.80, y), w * 0.20, h, transform=ax.transAxes, facecolor=COLORS["test"], edgecolor=COLORS["test_edge"], linewidth=1.5))
    add_text(ax, x + w * 0.40, y + h / 2, "80% development region: five OOF folds", 18.5, "bold", ha="center", va="center", color=COLORS["green"])
    add_text(ax, x + w * 0.90, y + h / 2, "20% locked test\nnever trains models", 15.4, "bold", ha="center", va="center", color=COLORS["purple"], linespacing=1.08)


def draw_oof_fold_matrix(ax):
    x0, y0 = 0.125, 0.360
    w, h = 0.420, 0.285
    add_text(ax, x0 - 0.005, y0 + h + 0.060, "Rotate which fold is held out", 21.5, "bold")
    add_text(ax, x0 - 0.005, y0 + h + 0.027, "Green folds train BDT+MLP; the gold fold is scored but not trained on.", 15.2, color=COLORS["muted"])
    col_w = w / 5.0
    row_h = h / 5.0
    for row in range(5):
        y = y0 + h - (row + 1) * row_h
        add_text(ax, x0 - 0.025, y + row_h / 2, f"pass {row + 1}", 14.5, "bold", ha="right", va="center", color=COLORS["ink"])
        for col in range(5):
            x = x0 + col * col_w
            held = row == col
            fc = COLORS["heldout"] if held else COLORS["fold"]
            ec = COLORS["heldout_edge"] if held else COLORS["fold_edge"]
            ax.add_patch(Rectangle((x, y), col_w - 0.003, row_h - 0.004, transform=ax.transAxes, facecolor=fc, edgecolor=ec, linewidth=1.05))
            if held:
                add_text(ax, x + col_w / 2, y + row_h / 2, f"score\nF{row}", 12.6, "bold", ha="center", va="center", color=COLORS["gold"], linespacing=0.95)
    for col in range(5):
        add_text(ax, x0 + col_w * (col + 0.5), y0 - 0.016, f"F{col}", 14.5, "bold", ha="center", va="top", color=COLORS["muted"])


def draw_oof_stack_table(ax):
    x, y, w, h = 0.612, 0.360, 0.338, 0.320
    box(ax, x, y, w, h, COLORS["stack"], COLORS["stack_edge"], lw=1.5, radius=0.010)
    add_text(ax, x + 0.024, y + h - 0.035, "Why do five passes help?", 20.0, "bold", color=COLORS["dark_blue"])
    add_text(
        ax,
        x + 0.024,
        y + h - 0.080,
        "Each pass unlocks one fold for the stacker.",
        13.2,
        color=COLORS["muted"],
        linespacing=1.12,
    )
    add_text(ax, x + 0.043, y + h - 0.128, "Simple holdout:", 14.4, "bold", va="center", color=COLORS["ink"])
    add_text(ax, x + 0.176, y + h - 0.128, "one 16% gold piece", 12.8, va="center", color=COLORS["muted"])

    add_text(ax, x + 0.043, y + h - 0.176, "Out-of-fold:", 14.4, "bold", va="center", color=COLORS["ink"])
    add_text(ax, x + 0.176, y + h - 0.176, "rotate the gold piece", 12.8, va="center", color=COLORS["muted"])

    tile_y = y + 0.078
    tile_h = 0.055
    tile_w = 0.036
    gap = 0.008
    start_x = x + 0.052
    for idx in range(5):
        tx = start_x + idx * (tile_w + gap)
        ax.add_patch(
            Rectangle(
                (tx, tile_y),
                tile_w,
                tile_h,
                transform=ax.transAxes,
                facecolor=COLORS["heldout"],
                edgecolor=COLORS["heldout_edge"],
                linewidth=1.2,
                zorder=2,
            )
        )
        add_text(ax, tx + tile_w / 2, tile_y + tile_h / 2, "16%", 10.9, "bold", ha="center", va="center", color=COLORS["gold"])
        if idx < 4:
            add_text(ax, tx + tile_w + gap / 2, tile_y + tile_h / 2, "+", 12.5, "bold", ha="center", va="center", color=COLORS["muted"])
    add_text(ax, start_x + 5 * (tile_w + gap) + 0.010, tile_y + tile_h / 2, "= 80%", 16.5, "bold", ha="left", va="center", color=COLORS["green"])
    add_text(ax, x + 0.043, y + 0.045, "Stitch the five honest pieces together.", 13.2, color=COLORS["ink"], va="center")


def draw_oof_locked_test(ax):
    x, y, w, h = 0.612, 0.175, 0.338, 0.120
    box(ax, x, y, w, h, COLORS["test"], COLORS["test_edge"], lw=1.5, radius=0.010)
    add_text(ax, x + 0.030, y + 0.078, "Locked-test evaluation", 20, "bold", color=COLORS["purple"], va="center")
    add_text(
        ax,
        x + 0.030,
        y + 0.038,
        "The separate 20% is scored only\nafter training choices are fixed.",
        13.4,
        va="center",
        linespacing=1.10,
    )


def draw_analogy_mapping_band(ax):
    box(ax, 0.070, 0.760, 0.880, 0.082, COLORS["soft_gray"], COLORS["line"], lw=1.0, radius=0.008)
    add_text(
        ax,
        0.095,
        0.803,
        "Analogy:",
        17.0,
        "bold",
        va="center",
        color=COLORS["ink"],
    )
    mapping = [
        ("events", "students"),
        ("BDT + MLP", "teachers"),
        ("scores", "grades"),
        ("stacker", "judge"),
    ]
    start_x = 0.205
    card_w = 0.158
    gap = 0.025
    for idx, (analysis_word, analogy_word) in enumerate(mapping):
        x = start_x + idx * (card_w + gap)
        box(ax, x, 0.778, card_w, 0.045, "white", COLORS["line"], lw=1.0, radius=0.006)
        add_text(ax, x + card_w / 2, 0.806, analysis_word, 12.4, "bold", ha="center", va="center")
        add_text(ax, x + card_w / 2, 0.784, analogy_word, 11.5, ha="center", va="center", color=COLORS["muted"])


def draw_student_tiles(ax, x, y, w, h, graded_idx, labels=True, alpha_train=1.0):
    tile_w = w / 5.0
    for idx in range(5):
        tx = x + idx * tile_w
        is_graded = idx == graded_idx
        fc = COLORS["heldout"] if is_graded else COLORS["fold"]
        ec = COLORS["heldout_edge"] if is_graded else COLORS["fold_edge"]
        rect = Rectangle(
            (tx, y),
            tile_w - 0.004,
            h,
            transform=ax.transAxes,
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.2,
            alpha=1.0 if is_graded else alpha_train,
        )
        ax.add_patch(rect)
        if labels:
            text = "graded" if is_graded else "taught"
            color = COLORS["gold"] if is_graded else COLORS["green"]
            add_text(ax, tx + tile_w / 2, y + h / 2, text, 11.8, "bold", ha="center", va="center", color=color)


def draw_simple_analogy_panel(ax):
    x, y, w, h = 0.070, 0.225, 0.405, 0.325
    box(ax, x, y, w, h, "white", COLORS["line"], lw=1.3, radius=0.010)
    add_text(ax, x + 0.027, y + h - 0.039, "Simple honest holdout", 20.0, "bold", color=COLORS["ink"])
    add_text(ax, x + 0.027, y + h - 0.078, "One 16% group gets graded for the stacker.", 13.7, color=COLORS["muted"])

    add_text(ax, x + 0.027, y + h - 0.113, "Teacher studies 64%", 15.0, "bold", color=COLORS["green"])
    add_text(ax, x + 0.270, y + h - 0.113, "grades 16%", 15.0, "bold", color=COLORS["gold"])
    draw_student_tiles(ax, x + 0.037, y + h - 0.183, w - 0.074, 0.052, graded_idx=4)

    arrow(ax, x + w * 0.50, y + h - 0.208, x + w * 0.50, y + 0.145, COLORS["muted"], lw=1.9)
    box(ax, x + 0.058, y + 0.033, w - 0.116, 0.104, COLORS["stack"], COLORS["stack_edge"], lw=1.25, radius=0.008)
    add_text(ax, x + w / 2, y + 0.109, "Judge learns from", 13.8, "bold", ha="center", va="center", color=COLORS["dark_blue"])
    add_text(ax, x + w / 2, y + 0.075, "only one honest graded slice", 13.2, ha="center", va="center")
    add_text(ax, x + w / 2, y + 0.046, "16% stack-training coverage", 14.5, "bold", ha="center", va="center", color=COLORS["gold"])


def draw_oof_analogy_panel(ax):
    x, y, w, h = 0.525, 0.225, 0.425, 0.325
    box(ax, x, y, w, h, COLORS["stack"], COLORS["stack_edge"], lw=1.4, radius=0.010)
    add_text(ax, x + 0.027, y + h - 0.039, "Out-of-fold rotation", 20.0, "bold", color=COLORS["dark_blue"])
    add_text(ax, x + 0.027, y + h - 0.077, "Same 64% teaches, but the 16% graded slice rotates.", 13.4, color=COLORS["muted"])

    grid_x = x + 0.040
    grid_y = y + 0.068
    grid_w = w - 0.080
    row_h = 0.027
    row_gap = 0.004
    for row in range(5):
        ry = grid_y + (4 - row) * (row_h + row_gap)
        add_text(ax, grid_x - 0.012, ry + row_h / 2, f"{row + 1}", 12.2, "bold", ha="right", va="center", color=COLORS["muted"])
        draw_student_tiles(ax, grid_x, ry, grid_w, row_h, graded_idx=row, labels=False, alpha_train=0.85)
        col_w = grid_w / 5.0
        add_text(
            ax,
            grid_x + col_w * (row + 0.5),
            ry + row_h / 2,
            "grade",
            9.2,
            "bold",
            ha="center",
            va="center",
            color=COLORS["gold"],
        )

    arrow(ax, x + w * 0.50, y + 0.073, x + w * 0.50, y + 0.054, COLORS["muted"], lw=1.9)
    box(ax, x + 0.055, y + 0.001, w - 0.110, 0.052, "white", COLORS["stack_edge"], lw=1.15, radius=0.008)
    add_text(ax, x + w / 2, y + 0.035, "five honest slices stitch together", 12.3, "bold", ha="center", va="center", color=COLORS["dark_blue"])
    add_text(ax, x + w / 2, y + 0.016, "16% + 16% + 16% + 16% + 16% = 80%", 11.2, "bold", ha="center", va="center", color=COLORS["green"])


def draw_analogy_split_bar(ax):
    x, y, w, h = 0.070, 0.580, 0.880, 0.058
    add_text(ax, x, y + h + 0.028, "Same event split", 15.3, "bold", va="center")
    add_text(ax, x + 0.175, y + h + 0.028, "100% events = 80% development + 20% locked test", 13.2, "bold", color=COLORS["muted"], va="center")
    dev_w = w * 0.80
    fold_w = dev_w / 5.0
    for idx in range(5):
        fx = x + idx * fold_w
        ax.add_patch(
            Rectangle(
                (fx, y),
                fold_w - 0.003,
                h,
                transform=ax.transAxes,
                facecolor=COLORS["fold"],
                edgecolor=COLORS["fold_edge"],
                linewidth=1.1,
            )
        )
        add_text(ax, fx + fold_w / 2, y + h / 2, f"F{idx}\n16%", 10.3, "bold", ha="center", va="center", color=COLORS["green"], linespacing=0.95)
    ax.add_patch(
        Rectangle(
            (x + dev_w, y),
            w * 0.20,
            h,
            transform=ax.transAxes,
            facecolor=COLORS["test"],
            edgecolor=COLORS["test_edge"],
            linewidth=1.2,
        )
    )
    add_text(ax, x + dev_w + w * 0.10, y + h / 2, "locked test\n20%", 10.8, "bold", ha="center", va="center", color=COLORS["purple"], linespacing=0.95)


def draw_stack_analogy_slide(ax) -> None:
    add_text(ax, 0.045, 0.935, "OOF lets every event get an honest score", 31, "bold")
    add_text(
        ax,
        0.045,
        0.884,
        "Like grading: every group gets scored by a teacher that did not teach that group.",
        17.4,
        color=COLORS["muted"],
    )
    draw_analogy_mapping_band(ax)

    box(ax, 0.070, 0.700, 0.880, 0.042, "#fff2c2", "#d8bf5f", lw=1.1, radius=0.007)
    add_text(ax, 0.095, 0.721, "Rule:", 14.5, "bold", va="center", color=COLORS["gold"])
    add_text(
        ax,
        0.160,
        0.721,
        "the stacker can learn only from scores on events the base models did not train on.",
        13.7,
        va="center",
        color=COLORS["ink"],
    )

    draw_analogy_split_bar(ax)
    draw_simple_analogy_panel(ax)
    draw_oof_analogy_panel(ax)
    arrow(ax, 0.482, 0.407, 0.518, 0.407, COLORS["muted"], lw=2.2)

    box(ax, 0.070, 0.055, 0.880, 0.125, "#eef7f2", COLORS["fold_edge"], lw=1.25, radius=0.010)
    add_text(ax, 0.095, 0.118, "Plain-English takeaway", 17.8, "bold", color=COLORS["green"], va="center")
    add_text(
        ax,
        0.330,
        0.140,
        "Simple holdout: the stacker sees one honest 16% slice.",
        15.2,
        "bold",
        va="center",
    )
    add_text(
        ax,
        0.330,
        0.101,
        "OOF: repeat the same 64% teach / 16% grade rule five times.",
        14.8,
        va="center",
        color=COLORS["ink"],
    )
    add_text(
        ax,
        0.330,
        0.072,
        "After stitching the five graded slices, the stacker sees honest scores for the full 80%.",
        14.8,
        va="center",
        color=COLORS["ink"],
    )


def draw_oof_data_efficient_contract(ax) -> None:
    add_text(ax, 0.045, 0.935, "Out-of-fold makes honest stacking data-efficient", 31, "bold")
    add_text(
        ax,
        0.045,
        0.876,
        "OOF: base models train on four folds, score the fifth, then stitch honest scores for the stacker.",
        16.3,
        color=COLORS["muted"],
    )
    draw_oof_partition_bar(ax)
    draw_oof_fold_matrix(ax)
    draw_oof_stack_table(ax)
    draw_oof_locked_test(ax)
    arrow(ax, 0.558, 0.505, 0.612, 0.508, COLORS["muted"])
    arrow(ax, 0.781, 0.360, 0.781, 0.297, COLORS["muted"])

    box(ax, 0.125, 0.175, 0.420, 0.120, "white", COLORS["line"], lw=1.2, radius=0.009)
    add_text(ax, 0.150, 0.258, "Honesty preserved in every row", 18.3, "bold", color=COLORS["dark_blue"], va="center")
    add_text(ax, 0.150, 0.216, "For any fold, BDT/MLP scores come from models\ntrained on the other four folds.", 14.2, color=COLORS["ink"], va="center", linespacing=1.20)

    box(ax, 0.070, 0.055, 0.880, 0.095, "#fff2c2", "#d8bf5f")
    add_text(ax, 0.095, 0.112, "Takeaway", 20, "bold", color=COLORS["gold"], va="center")
    add_text(
        ax,
        0.215,
        0.124,
        "OOF rotates the simple held-out scoring rule across all five folds.",
        16.0,
        "bold",
        color=COLORS["ink"],
        va="center",
    )
    add_text(
        ax,
        0.215,
        0.084,
        "Stitching the five honest 16% pieces gives full 80% stack coverage; the locked test stays untouched.",
        15.6,
        "normal",
        color=COLORS["ink"],
        va="center",
    )


def draw_partition_bar(ax):
    x, y, w, h = 0.070, 0.690, 0.330, 0.070
    box(ax, x - 0.012, y - 0.055, w + 0.024, h + 0.115, "white", COLORS["line"])
    add_text(ax, x, y + 0.105, "Event-level split", 21, "bold")
    add_text(ax, x, y - 0.012, "Event key: source_sample / run / evt", 13.5, color=COLORS["muted"])
    ax.add_patch(Rectangle((x, y), w * 0.80, h, transform=ax.transAxes, facecolor=COLORS["fold"], edgecolor=COLORS["fold_edge"], linewidth=1.4))
    ax.add_patch(Rectangle((x + w * 0.80, y), w * 0.20, h, transform=ax.transAxes, facecolor=COLORS["test"], edgecolor=COLORS["test_edge"], linewidth=1.4))
    add_text(ax, x + w * 0.40, y + 0.036, "80% stack-train region", 16, "bold", ha="center", va="center", color=COLORS["green"])
    add_text(ax, x + w * 0.90, y + 0.036, "20% locked\ntest", 14, "bold", ha="center", va="center", color=COLORS["purple"], linespacing=0.95)


def draw_fold_matrix(ax):
    x0, y0 = 0.070, 0.260
    w, h = 0.465, 0.330
    box(ax, x0 - 0.012, y0 - 0.055, w + 0.024, h + 0.105, "white", COLORS["line"])
    add_text(ax, x0, y0 + h + 0.065, "Five out-of-fold passes", 23, "bold")
    add_text(ax, x0, y0 + h + 0.028, "Train BDT+MLP on four folds, then score the one held-out fold.", 15, color=COLORS["muted"])
    col_w = w / 5.0
    row_h = h / 5.0
    for row in range(5):
        y = y0 + h - (row + 1) * row_h
        add_text(ax, x0 - 0.025, y + row_h / 2, f"Fold {row}", 15, "bold", ha="right", va="center")
        for col in range(5):
            x = x0 + col * col_w
            held = row == col
            fc = COLORS["heldout"] if held else COLORS["fold"]
            ec = COLORS["heldout_edge"] if held else COLORS["fold_edge"]
            ax.add_patch(Rectangle((x, y), col_w - 0.003, row_h - 0.004, transform=ax.transAxes, facecolor=fc, edgecolor=ec, linewidth=1.0))
            if held:
                add_text(ax, x + col_w / 2, y + row_h / 2, "score", 14, "bold", ha="center", va="center", color=COLORS["gold"])
    for col in range(5):
        add_text(ax, x0 + col_w * (col + 0.5), y0 - 0.014, f"F{col}", 14, "bold", ha="center", va="top", color=COLORS["muted"])
    add_text(ax, x0 + w + 0.018, y0 + h * 0.63, "green =\ntrain", 15, "bold", color=COLORS["green"], linespacing=1.05)
    add_text(ax, x0 + w + 0.018, y0 + h * 0.34, "gold =\nscore", 15, "bold", color=COLORS["gold"], linespacing=1.05)


def draw_stack_panel(ax):
    x, y, w, h = 0.615, 0.310, 0.335, 0.475
    box(ax, x, y, w, h, COLORS["stack"], COLORS["stack_edge"])
    add_text(ax, x + 0.025, y + h - 0.035, "Stack training uses honest scores", 22, "bold", color=COLORS["dark_blue"])
    add_text(
        ax,
        x + 0.025,
        y + h - 0.090,
        "For each event in the 80% region:",
        16,
        "bold",
    )
    rows = [
        ("BDT score", "held-out score for this event", "#d9e9f5"),
        ("MLP score", "held-out score for this event", "#e5ddf4"),
        ("truth label", "is_signal is the target", "#f7e8b5"),
    ]
    yy = y + h - 0.155
    for label, body, fc in rows:
        box(ax, x + 0.030, yy - 0.056, w - 0.060, 0.060, fc, COLORS["line"], lw=1.0, radius=0.006)
        add_text(ax, x + 0.045, yy - 0.010, label, 16, "bold", va="top")
        add_text(ax, x + 0.150, yy - 0.010, body, 15, color=COLORS["ink"], va="top")
        yy -= 0.075
    box(ax, x + 0.052, y + 0.018, w - 0.104, 0.082, "white", COLORS["stack_edge"], lw=1.1)
    add_text(ax, x + w / 2, y + 0.066, "Train logistic / GBM / MLP stackers", 16.5, "bold", ha="center", va="center", color=COLORS["dark_blue"])
    add_text(ax, x + w / 2, y + 0.038, "inputs = BDT score + MLP score", 14, ha="center", va="center", color=COLORS["muted"])


def draw_test_panel(ax):
    x, y, w, h = 0.615, 0.180, 0.335, 0.105
    box(ax, x, y, w, h, COLORS["test"], COLORS["test_edge"])
    add_text(ax, x + 0.025, y + h - 0.027, "Locked-test evaluation", 20, "bold", color=COLORS["purple"])
    add_text(ax, x + 0.025, y + h - 0.064, "Final BDT+MLP train on 80% non-test rows.\nThen they score the untouched 20% test rows.", 13.8, linespacing=1.03)


def draw_bottom_band(ax):
    box(ax, 0.070, 0.050, 0.880, 0.095, "#fff2c2", "#d8bf5f")
    add_text(ax, 0.095, 0.122, "Leakage control", 20, "bold", color=COLORS["gold"])
    add_text(
        ax,
        0.270,
        0.122,
        "Stack rows use BDT/MLP scores from models that excluded that event.",
        17.5,
        "bold",
        va="top",
    )
    add_text(
        ax,
        0.270,
        0.083,
        "Cost: 5 fold BDTs + 5 fold MLPs + final BDT/MLP per domain before stack metrics are available.",
        14.2,
        color=COLORS["muted"],
        va="top",
    )


def render_original_slide() -> tuple[Path, Path, Path]:
    setup()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, ax = fig_ax()
    add_text(ax, 0.045, 0.935, "The stack is trained on out-of-fold base scores", 31, "bold")
    add_text(
        ax,
        0.045,
        0.884,
        "BDT and MLP scores used for stack training come from models that did not train on that event.",
        18.5,
        color=COLORS["muted"],
    )

    draw_partition_bar(ax)
    draw_fold_matrix(ax)
    draw_stack_panel(ax)
    draw_test_panel(ax)
    draw_bottom_band(ax)

    arrow(ax, 0.410, 0.710, 0.615, 0.695, COLORS["muted"])
    arrow(ax, 0.552, 0.430, 0.615, 0.500, COLORS["muted"])
    arrow(ax, 0.782, 0.310, 0.782, 0.285, COLORS["muted"])

    box(ax, 0.420, 0.645, 0.145, 0.078, COLORS["soft_gray"], COLORS["line"])
    add_text(ax, 0.4925, 0.687, "repeat for\nall 5 folds", 16, "bold", ha="center", va="center", linespacing=1.05)

    png = OUTDIR / "fresh_oof_stack_strategy_contract_v1.png"
    script = OUTDIR / "fresh_oof_stack_strategy_contract_v1.speaker.md"
    manifest = OUTDIR / "fresh_oof_stack_strategy_contract_v1.json"
    fig.savefig(png, dpi=DPI)
    plt.close(fig)

    script.write_text(
        "# Fresh OOF Stack Strategy Script - Leakage-Controlled Stacking\n\n"
        "At a high level, this slide is explaining why the stack is not trained on circular information.\n\n"
        "First, I split the rows at the event level using the source sample, run, and event number. "
        "Twenty percent of events are locked away as the final test set. The remaining eighty percent is the stack-training region.\n\n"
        "Inside that eighty percent, I make five folds. For each fold, I train the base BDT and the base MLP on the other four folds, "
        "and I only use those models to score the held-out fold. That means every row used to train the stack has a BDT score and an MLP score "
        "from models that did not train on that event.\n\n"
        "The stack then learns from those two honest base scores, together with the truth label. "
        "For the final locked-test comparison, I train final base models only on the non-test region and evaluate the whole model chain on the untouched test events.\n\n"
        "The cost is that this is computationally heavier than a quick stack on existing score caches: each domain has to train five fold BDTs, five fold MLPs, and then final base models before the stack metrics are available.\n"
    )
    manifest.write_text(
        json.dumps(
            {
                "schema": "RJ_FRESH_OOF_STACK_STRATEGY_SLIDE_V1",
                "png": str(png),
                "speaker_script": str(script),
                "dimensions_px": [W, H],
                "training_class_definition": "signal = is_signal == 1; background = is_signal == 0",
                "partition_contract": "event key = source_sample/run/evt; locked test = 20%; OOF stack train region = remaining 80% split into 5 folds",
                "base_score_contract": "fold k base scores are produced by BDT/MLP models trained on folds excluding k",
                "locked_test_contract": "final BDT/MLP models train on non-test events only, then score locked-test rows for final comparison",
                "source_code": str(Path(__file__)),
                "audience_visible_caveat": "Cost: five fold BDTs, five fold MLPs, and final base models per domain before stack metrics are available.",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return png, script, manifest


def render_simple_holdout_slide() -> tuple[Path, Path, Path, Path]:
    setup()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, ax = fig_ax()
    draw_simple_holdout_scaffold(ax)

    stem = "fresh_stack_honest_scores_simple_holdout_scaffold_v1"
    png, pdf = save_png_pdf(fig, stem)
    script = OUTDIR / f"{stem}.speaker.md"
    manifest = OUTDIR / f"{stem}.json"
    script.write_text(
        "# Honest Base-Score Scaffold Script\n\n"
        "Before introducing folds, I would first explain the simple version of honest stacking.\n\n"
        "The stacker is not supposed to learn from base-model scores on events that those same base models trained on. "
        "In the simple held-out version, sixty-four percent of the events train the base BDT and MLP. "
        "Those fixed base models then score a separate sixteen percent block, and that block trains the stacker. "
        "The final twenty percent remains locked away for the final metric comparison.\n\n"
        "The stackers being tested use the BDT and MLP scores as the extra learned inputs. "
        "The score-only stack uses just those two scores, while the context stack also keeps the candidate kinematic context used by the campaign, such as cluster transverse energy and centrality for AuAu.\n\n"
        "This is valid because the stacker sees base scores from events excluded from base-model training. "
        "The reason we do not stop here is that it is data-inefficient: separate event blocks are consumed for base training, stack training, and final testing.\n"
    )
    manifest.write_text(
        json.dumps(
            {
                "schema": "RJ_FRESH_STACK_HONEST_SCORES_SCAFFOLD_V1",
                "png": str(png),
                "pdf": str(pdf),
                "speaker_script": str(script),
                "dimensions_px": [W, H],
                "purpose": "Introduce the simple leakage-controlled stacking idea before the OOF fold implementation.",
                "training_class_definition": "signal = is_signal == 1; background = is_signal == 0; source_sample is provenance and part of the event key, not the supervised class label",
                "event_key": "source_sample/run/evt",
                "visual_contract": "Simple holdout uses 64% base training, 16% stack training, and 20% locked test. Set A trains BDT+MLP; Set B is scored by those base models and trains the stacker; Set C is the locked test set.",
                "stack_input_contract": "score_only stack inputs = bdt_score, mlp_score; score_context stack inputs = bdt_score, mlp_score, cluster_Et, plus centrality for AuAu.",
                "source_code": str(Path(__file__)),
                "provenance_source_slide": str(OUTDIR / "fresh_oof_stack_strategy_contract_v1.png"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return png, pdf, script, manifest


def render_oof_contract_slide() -> tuple[Path, Path, Path, Path]:
    setup()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, ax = fig_ax()
    draw_oof_data_efficient_contract(ax)

    stem = "fresh_stack_oof_data_efficient_contract_v1"
    png, pdf = save_png_pdf(fig, stem)
    script = OUTDIR / f"{stem}.speaker.md"
    manifest = OUTDIR / f"{stem}.json"
    script.write_text(
        "# OOF Data-Efficient Stack Strategy Script\n\n"
        "After the simple honest-stacking scaffold, this slide explains the actual campaign implementation.\n\n"
        "OOF means out-of-fold: an event is scored by a model that did not train on that event. "
        "The event-level split first locks away twenty percent of events for the final test set. "
        "Inside the remaining eighty percent, the campaign makes five folds. For each pass, the base BDT and MLP train on four folds and score only the held-out fold.\n\n"
        "The change from simple holdout is coverage, not the leakage rule. "
        "Instead of using one sixteen percent block to train the stacker, OOF rotates the same held-out scoring rule across all five folds. "
        "Each pass gives one fold honest base-model scores; after five passes, five times sixteen percent becomes the full eighty percent development region. "
        "When those five held-out score blocks are stitched together, every development event has an honest BDT score and an honest MLP score. "
        "The locked test events remain separate and are only used after all training choices are fixed.\n"
    )
    manifest.write_text(
        json.dumps(
            {
                "schema": "RJ_FRESH_STACK_OOF_DATA_EFFICIENT_CONTRACT_V1",
                "png": str(png),
                "pdf": str(pdf),
                "speaker_script": str(script),
                "dimensions_px": [W, H],
                "oof_definition": "out-of-fold = an event is scored by a model that did not train on that event",
                "partition_contract": "event key = source_sample/run/evt; locked test = 20%; OOF development region = remaining 80% split into 5 folds",
                "base_score_contract": "fold k base scores are produced by BDT/MLP models trained on folds excluding k",
                "training_class_definition": "signal = is_signal == 1; background = is_signal == 0; source_sample is provenance and part of the event key, not the supervised class label",
                "locked_test_contract": "locked-test events are not used for stack or base training; they are scored only for final metrics",
                "source_code": str(Path(__file__)),
                "provenance_source_slide": str(OUTDIR / "fresh_oof_stack_strategy_contract_v1.png"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return png, pdf, script, manifest


def render_teacher_grader_analogy_slide() -> tuple[Path, Path, Path, Path]:
    setup()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, ax = fig_ax()
    draw_stack_analogy_slide(ax)

    stem = "fresh_stack_teacher_grader_oof_analogy_v2"
    png, pdf = save_png_pdf(fig, stem)
    script = OUTDIR / f"{stem}.speaker.md"
    manifest = OUTDIR / f"{stem}.json"
    script.write_text(
        "# Teacher-Grader Analogy Script - Why OOF Recovers 80 Percent\n\n"
        "The simplest way I think about the out-of-fold strategy is like grading.\n\n"
        "The events are like students, the BDT and MLP are like teachers, their scores are like grades, and the stacker is like a judge learning from those grades. "
        "The key rule is that a fair grade should come from a teacher who did not teach that student. In the ML language, the stacker should train only on BDT and MLP scores from models that did not train on that event.\n\n"
        "In the simple holdout version, four groups teach the base models and one separate group gets graded. That is honest, but only that one group can train the stacker, so the stacker only gets one sixteen-percent slice of the development data.\n\n"
        "Out-of-fold keeps the same honesty rule but rotates it. Each pass changes which group is held out and graded. After five passes, every group has a fair grade from a model that did not train on it. Then those five fair-grade blocks are stitched together.\n\n"
        "So OOF is not relaxing the leakage control. It is just rotating the held-out group, which gives the stacker honest BDT and MLP scores for the full eighty-percent development region while keeping the locked test untouched.\n"
    )
    manifest.write_text(
        json.dumps(
            {
                "schema": "RJ_FRESH_STACK_TEACHER_GRADER_ANALOGY_V2",
                "png": str(png),
                "pdf": str(pdf),
                "speaker_script": str(script),
                "dimensions_px": [W, H],
                "purpose": "Audience-facing teacher/grader analogy explaining the numeric split that lets OOF give honest stack-training scores for the full 80% development region.",
                "analogy_map": {
                    "events": "students",
                    "BDT_MLP": "teachers",
                    "scores": "grades",
                    "stacker": "judge",
                },
                "core_rule": "the stacker trains only on base scores from events excluded from base-model training",
                "event_split_message": "100% events = 80% development region plus 20% locked test; the 80% development region is five 16% folds",
                "simple_holdout_message": "one 64% base-training block scores one 16% stack-training block",
                "oof_message": "each pass uses the same 64% teach and 16% grade rule, rotating the graded fold across all five folds so five honest 16% blocks stitch into 80% stack coverage",
                "locked_test_contract": "locked-test events remain separate and are not used for training choices",
                "source_code": str(Path(__file__)),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return png, pdf, script, manifest


def write_split_note(
    simple_paths: tuple[Path, Path, Path, Path],
    oof_paths: tuple[Path, Path, Path, Path],
    analogy_paths: tuple[Path, Path, Path, Path],
) -> Path:
    note = OUTDIR / "fresh_stack_oof_scaffold_split_v1.md"
    note.write_text(
        "# Fresh OOF stack strategy scaffold split v1\n\n"
        "## Source\n"
        "- Preserved source slide: `fresh_oof_stack_strategy_contract_v1.png`\n"
        "- Generator updated: `scripts/slides/auau_bdt/stacking/make_fresh_oof_stack_strategy_slide.py`\n\n"
        "## New recommended opening sequence\n"
        f"1. `{analogy_paths[0].name}` gives the teacher/grader analogy for why rotation recovers 80% without breaking honesty.\n"
        f"2. `{simple_paths[0].name}` introduces the simple honest stacking rule without folds.\n"
        f"3. `{oof_paths[0].name}` shows the five-fold OOF implementation as the data-efficient extension of that rule.\n\n"
        "## What changed\n"
        "- Split the overloaded one-slide explanation into a simple-holdout scaffold followed by the actual OOF campaign contract.\n"
        "- Kept the current slide style: white background, Times New Roman, rounded semantic boxes, green train region, gold held-out scoring, purple locked test, blue stack table, and yellow takeaway band.\n"
        "- Made the locked test visually separate from stack training and stated that source_sample is provenance/event-key material, not the supervised class label.\n"
        "- Defined OOF before using the acronym and removed the compute-cost message from the opening scaffold so the audience first sees the leakage-control idea.\n\n"
        "## Files\n"
        f"- Analogy PNG/PDF/JSON/script: `{analogy_paths[0].name}`, `{analogy_paths[1].name}`, `{analogy_paths[3].name}`, `{analogy_paths[2].name}`\n"
        f"- Simple scaffold PNG/PDF/JSON/script: `{simple_paths[0].name}`, `{simple_paths[1].name}`, `{simple_paths[3].name}`, `{simple_paths[2].name}`\n"
        f"- OOF contract PNG/PDF/JSON/script: `{oof_paths[0].name}`, `{oof_paths[1].name}`, `{oof_paths[3].name}`, `{oof_paths[2].name}`\n"
    )
    return note


def render() -> dict[str, str]:
    original_png, original_script, original_manifest = render_original_slide()
    analogy_paths = render_teacher_grader_analogy_slide()
    simple_paths = render_simple_holdout_slide()
    oof_paths = render_oof_contract_slide()
    note = write_split_note(simple_paths, oof_paths, analogy_paths)
    return {
        "original_png": str(original_png),
        "original_speaker_script": str(original_script),
        "original_manifest": str(original_manifest),
        "analogy_png": str(analogy_paths[0]),
        "analogy_pdf": str(analogy_paths[1]),
        "analogy_speaker_script": str(analogy_paths[2]),
        "analogy_manifest": str(analogy_paths[3]),
        "simple_png": str(simple_paths[0]),
        "simple_pdf": str(simple_paths[1]),
        "simple_speaker_script": str(simple_paths[2]),
        "simple_manifest": str(simple_paths[3]),
        "oof_png": str(oof_paths[0]),
        "oof_pdf": str(oof_paths[1]),
        "oof_speaker_script": str(oof_paths[2]),
        "oof_manifest": str(oof_paths[3]),
        "split_note": str(note),
    }


def main() -> None:
    print(json.dumps(render(), sort_keys=True))


if __name__ == "__main__":
    main()
