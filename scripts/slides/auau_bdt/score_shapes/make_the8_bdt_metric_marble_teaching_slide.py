#!/usr/bin/env python3
"""Build a pedagogical toy slide explaining BDT metric behavior."""

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
import json
import math
from collections import Counter
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUTDIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "fixed_sample_controls/bdt_metric_marble_teaching_20260611"
)
OUTPNG = OUTDIR / "bdt_metric_marble_teaching_slide.png"
OUTNOTE = OUTDIR / "bdt_metric_marble_teaching_note.md"
OUTSCRIPT = OUTDIR / "bdt_metric_marble_teaching_speaker_script.md"
OUTCSV = OUTDIR / "bdt_metric_marble_teaching_toy_metrics.csv"
OUTMANIFEST = OUTDIR / "bdt_metric_marble_teaching_manifest.json"

NAVY = "#101828"
INK = "#1D2939"
MUTED = "#475467"
BORDER = "#D5DCE8"
GRID = "#E5EAF2"
RED = "#C5392F"
BLUE = "#1F77B4"
PURPLE = "#7651A6"
GOLD = "#C28400"
GREEN = "#228A4D"
SOFT_BLUE = "#F4F8FF"
SOFT_GOLD = "#FFF8E5"
SOFT_RED = "#FFF5F3"
SOFT_GRAY = "#F8FAFC"


TOY = {
    "old": {
        "title": "OLD MODEL",
        "red": [4.0, 5.0, 6.0, 6.0, 7.0],
        "blue": [1.0, 2.0, 3.5, 4.5, 8.5],
    },
    "new": {
        "title": "NEW MODEL",
        "red": [6.0, 7.0, 7.0, 8.0, 8.0],
        "blue": [0.5, 1.0, 2.0, 3.0, 8.5],
    },
}


def add_box(
    fig: plt.Figure,
    x: float,
    y: float,
    w: float,
    h: float,
    face: str,
    *,
    edge: str = BORDER,
    lw: float = 1.2,
    radius: float = 0.012,
    zorder: int = -2,
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def add_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    weight: str = "normal",
    color: str = INK,
    ha: str = "left",
    va: str = "top",
    linespacing: float = 1.15,
) -> None:
    fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
        family="Times New Roman",
    )


def auc(red: list[float], blue: list[float]) -> float:
    total = 0.0
    for r in red:
        for b in blue:
            if r > b:
                total += 1.0
            elif r == b:
                total += 0.5
    return total / (len(red) * len(blue))


def logloss(red: list[float], blue: list[float]) -> float:
    loss = 0.0
    n = 0
    for score in red:
        p = min(max(score / 10.0, 1e-6), 1.0 - 1e-6)
        loss -= math.log(p)
        n += 1
    for score in blue:
        p = min(max(score / 10.0, 1e-6), 1.0 - 1e-6)
        loss -= math.log(1.0 - p)
        n += 1
    return loss / n


def median(values: list[float]) -> float:
    vals = sorted(values)
    return vals[len(vals) // 2]


def wp80_threshold(red: list[float]) -> float:
    vals = sorted(red)
    return vals[1]


def wp80_fake(red: list[float], blue: list[float]) -> float:
    cut = wp80_threshold(red)
    return sum(1 for score in blue if score >= cut) / len(blue)


def metrics(red: list[float], blue: list[float]) -> dict[str, float]:
    return {
        "auc": auc(red, blue),
        "logloss": logloss(red, blue),
        "median_gap": median(red) - median(blue),
        "wp80_cut": wp80_threshold(red),
        "wp80_fake": wp80_fake(red, blue),
    }


def jittered(scores: list[float], y_base: float) -> list[tuple[float, float, float]]:
    counts = Counter(scores)
    used: Counter[float] = Counter()
    out: list[tuple[float, float, float]] = []
    for score in scores:
        idx = used[score]
        used[score] += 1
        if counts[score] == 1:
            dx = 0.0
        else:
            dx = (idx - (counts[score] - 1) / 2.0) * 0.22
        out.append((score + dx, y_base, score))
    return out


def score_label(score: float) -> str:
    return str(int(score)) if score == int(score) else f"{score:.1f}"


def draw_metric_value_row(fig: plt.Figure, x: float, y: float, w: float, m: dict[str, float]) -> None:
    parts = [
        ("AUC", f"{m['auc']:.2f}", BLUE),
        ("logloss", f"{m['logloss']:.2f}", PURPLE),
        ("gap", f"{m['median_gap']:.1f}", GOLD),
        ("WP80 fake", f"{100*m['wp80_fake']:.0f}%", RED),
    ]
    positions = [0.16, 0.39, 0.61, 0.84]
    for i, (label, value, color) in enumerate(parts):
        xx = x + w * positions[i]
        add_text(fig, xx, y, label, size=10.6, weight="bold", color=color, ha="center")
        add_text(fig, xx, y - 0.023, value, size=13.0, weight="bold", color=NAVY, ha="center")


def draw_ruler_panel(
    fig: plt.Figure,
    *,
    key: str,
    x: float,
    y: float,
    w: float,
    h: float,
    face: str,
    annotation: str,
) -> dict[str, float]:
    sample = TOY[key]
    red = sample["red"]
    blue = sample["blue"]
    m = metrics(red, blue)
    add_box(fig, x, y, w, h, face, edge=BORDER, lw=1.35, radius=0.016)
    add_text(fig, x + 0.024, y + h - 0.034, sample["title"], size=19.0, weight="bold", color=NAVY)
    add_text(fig, x + 0.024, y + h - 0.072, annotation, size=12.5, color=MUTED, linespacing=1.05)
    draw_metric_value_row(fig, x + 0.042, y + 0.073, w - 0.084, m)

    ax = fig.add_axes([x + 0.056, y + 0.154, w - 0.112, 0.130])
    ax.set_xlim(-0.2, 10.2)
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([])
    ax.set_xticks([0, 5, 10])
    ax.set_xticklabels(["0", "5", "10"], fontsize=11.5, fontfamily="Times New Roman")
    ax.set_xlabel("BDT-like score", fontsize=11.5, fontfamily="Times New Roman", labelpad=5, color=MUTED)
    ax.tick_params(axis="x", direction="inout", length=5, width=1.0, colors=INK)
    for spine in ("top", "left", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_position(("data", 0.50))
    ax.spines["bottom"].set_color("#667085")
    ax.spines["bottom"].set_linewidth(1.3)
    for tick in range(0, 11):
        ax.plot([tick, tick], [0.47, 0.53], color="#98A2B3", lw=0.8, zorder=0)

    cut = m["wp80_cut"]
    ax.axvline(cut, color=NAVY, lw=1.5, ls=(0, (5, 4)), alpha=0.85, zorder=1)
    ax.text(
        cut,
        0.95,
        f"WP80 cut = {score_label(cut)}",
        ha="center",
        va="top",
        fontsize=10.7,
        fontweight="bold",
        color=NAVY,
        fontfamily="Times New Roman",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.6},
    )

    for plot_score, yy, true_score in jittered(red, 0.68):
        ax.scatter([plot_score], [yy], s=315, color=RED, edgecolor="white", linewidth=1.6, zorder=3)
        ax.text(plot_score, yy, score_label(true_score), color="white", ha="center", va="center", fontsize=8.2, fontweight="bold")
    for plot_score, yy, true_score in jittered(blue, 0.32):
        ax.scatter([plot_score], [yy], s=315, color=BLUE, edgecolor="white", linewidth=1.6, zorder=3)
        ax.text(plot_score, yy, score_label(true_score), color="white", ha="center", va="center", fontsize=8.2, fontweight="bold")

    ax.text(-0.10, 0.68, "signal", color=RED, ha="right", va="center", fontsize=11.0, fontweight="bold", fontfamily="Times New Roman")
    ax.text(-0.10, 0.32, "background", color=BLUE, ha="right", va="center", fontsize=11.0, fontweight="bold", fontfamily="Times New Roman")

    ax.annotate(
        "hard blue tail\nstill signal-like",
        xy=(blue[-1], 0.32),
        xytext=(8.35, 0.06),
        arrowprops={"arrowstyle": "->", "color": BLUE, "lw": 1.2},
        fontsize=10.4,
        fontweight="bold",
        color=BLUE,
        fontfamily="Times New Roman",
        ha="center",
    )
    return m


def draw_explain_card(
    fig: plt.Figure,
    *,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    sees: str,
    body: str,
    accent: str,
    face: str,
) -> None:
    add_box(fig, x, y, w, h, face, edge=BORDER, lw=1.0, radius=0.012)
    fig.patches.append(
        Rectangle((x + 0.012, y + 0.020), 0.006, h - 0.040, transform=fig.transFigure, facecolor=accent, edgecolor="none", zorder=-1)
    )
    add_text(fig, x + 0.030, y + h - 0.030, title, size=16.2, weight="bold", color=NAVY)
    add_text(fig, x + 0.030, y + h - 0.070, sees, size=12.3, weight="bold", color=accent)
    add_text(fig, x + 0.030, y + h - 0.107, fill(body, 44), size=12.5, color=INK, linespacing=1.12)


def write_outputs(old_m: dict[str, float], new_m: dict[str, float]) -> None:
    rows = []
    for key, m in (("old", old_m), ("new", new_m)):
        rows.append(
            {
                "model": key,
                "red_scores": " ".join(score_label(v) for v in TOY[key]["red"]),
                "blue_scores": " ".join(score_label(v) for v in TOY[key]["blue"]),
                "auc": f"{m['auc']:.6f}",
                "logloss": f"{m['logloss']:.6f}",
                "median_gap": f"{m['median_gap']:.6f}",
                "wp80_cut": f"{m['wp80_cut']:.6f}",
                "wp80_fake_rate": f"{m['wp80_fake']:.6f}",
            }
        )
    with OUTCSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    note = f"""# BDT Metric Marble Teaching Slide

## Toy scores

Old model:
- Signal/red scores: {rows[0]['red_scores']}
- Background/blue scores: {rows[0]['blue_scores']}

New model:
- Signal/red scores: {rows[1]['red_scores']}
- Background/blue scores: {rows[1]['blue_scores']}

## Exact toy metrics

| Model | AUC | Logloss | Median score gap | WP80 cut | WP80 fake |
| --- | ---: | ---: | ---: | ---: | ---: |
| Old | {old_m['auc']:.2f} | {old_m['logloss']:.2f} | {old_m['median_gap']:.1f} | {old_m['wp80_cut']:.1f} | {100*old_m['wp80_fake']:.0f}% |
| New | {new_m['auc']:.2f} | {new_m['logloss']:.2f} | {new_m['median_gap']:.1f} | {new_m['wp80_cut']:.1f} | {100*new_m['wp80_fake']:.0f}% |

## Why this satisfies the pedagogical constraints

- AUC changes only from {old_m['auc']:.2f} to {new_m['auc']:.2f}: the old model already ranks most red marbles above most blue marbles, and the new model only fixes one additional red-blue pair.
- Logloss improves from {old_m['logloss']:.2f} to {new_m['logloss']:.2f}: red scores move higher and easy blue scores move lower, so the model is more confident on the bulk.
- Median gap doubles from {old_m['median_gap']:.1f} to {new_m['median_gap']:.1f}: the class centers visibly separate.
- WP80 fake stays {100*old_m['wp80_fake']:.0f}% to {100*new_m['wp80_fake']:.0f}%: the same single hard blue marble remains beyond the 80%-signal cut.

Real motivating case printed on the slide uses the existing slide/source values: fixed Jet12+20 validation, Jet12+20-trained BDT to Jet12+20+30+40-trained BDT; AUC 0.830 to 0.835, logloss 0.562 to 0.433, median gap 0.27 to 0.42, and WP80 fake about 31% with little change.
"""
    OUTNOTE.write_text(note)

    script = """# Speaker Script: Why AUC Can Barely Move

This is a toy version of what is happening in the BDT matrix.

The same five signal marbles and five background marbles are scored by two classifiers. In the old model, most red marbles are already to the right of most blue marbles, so the ranking is mostly correct. That is what AUC sees.

In the new model, the easy backgrounds move left and the signal marbles move right. Visually, the red and blue groups separate much more cleanly. Logloss and median BDT gap respond to that because they care about confidence and score distance, not only pair ordering.

But the hard blue tail is still there. The WP80 cut moves with the signal distribution, and the same one blue marble is still beyond the cut. That is why the WP80 fake rate can barely move even though the bulk score split becomes much cleaner.

So the key lesson is: a small AUC change does not mean the visible split is meaningless. It can mean the ranking was already mostly correct, while the new model is making the score scale more decisive.
"""
    OUTSCRIPT.write_text(script)

    OUTMANIFEST.write_text(
        json.dumps(
            {
                "schema": "BDT_METRIC_MARBLE_TEACHING_SLIDE_V1",
                "status": "READY",
                "slide_png": str(OUTPNG),
                "note": str(OUTNOTE),
                "speaker_script": str(OUTSCRIPT),
                "toy_metrics_csv": str(OUTCSV),
                "toy_scores": TOY,
                "toy_metrics": {"old": old_m, "new": new_m},
                "real_motivating_case": {
                    "validation": "Jet12+20 fixed validation, 0-20% centrality",
                    "training_comparison": "Jet12+20 trained BDT -> Jet12+20+30+40 trained BDT",
                    "auc": "0.830 -> 0.835",
                    "logloss": "0.562 -> 0.433",
                    "median_gap": "0.27 -> 0.42",
                    "wp80_fake": "about 31%, little change",
                },
                "google_slides_mutated": False,
            },
            indent=2,
        )
        + "\n"
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    old_m = metrics(TOY["old"]["red"], TOY["old"]["blue"])
    new_m = metrics(TOY["new"]["red"], TOY["new"]["blue"])

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "#FFFFFF",
            "savefig.facecolor": "#FFFFFF",
        }
    )
    fig = plt.figure(figsize=slide_figsize(SLIDE_DPI), dpi=SLIDE_DPI)
    fig.subplots_adjust(0, 0, 1, 1)

    add_text(fig, 0.055, 0.948, "Why the BDT split can look much cleaner while AUC barely changes", size=25.8, weight="bold", color=NAVY)
    add_text(
        fig,
        0.056,
        0.902,
        "Toy red-vs-blue marbles: the ranking changes little, but the score scale becomes more decisive.",
        size=14.5,
        color=MUTED,
    )

    old_metrics = draw_ruler_panel(
        fig,
        key="old",
        x=0.055,
        y=0.505,
        w=0.430,
        h=0.345,
        face="#FFFFFF",
        annotation="Mostly ordered, but not decisive",
    )
    new_metrics = draw_ruler_panel(
        fig,
        key="new",
        x=0.515,
        y=0.505,
        w=0.430,
        h=0.345,
        face="#FFFFFF",
        annotation="Bulk scores pull apart; tail remains",
    )

    draw_explain_card(
        fig,
        x=0.055,
        y=0.220,
        w=0.285,
        h=0.190,
        title="AUC sees ranking",
        sees="red-blue pair ordering",
        body=f"Most pairs are already ordered, so AUC gains only one pair: {old_metrics['auc']:.2f} to {new_metrics['auc']:.2f}.",
        accent=BLUE,
        face=SOFT_BLUE,
    )
    draw_explain_card(
        fig,
        x=0.357,
        y=0.220,
        w=0.285,
        h=0.190,
        title="Logloss and gap see decisiveness",
        sees="confidence and center separation",
        body=f"Reds move higher and easy blues move lower: gap {old_metrics['median_gap']:.1f} to {new_metrics['median_gap']:.1f}, logloss {old_metrics['logloss']:.2f} to {new_metrics['logloss']:.2f}.",
        accent=PURPLE,
        face="#F7F2FF",
    )
    draw_explain_card(
        fig,
        x=0.659,
        y=0.220,
        w=0.285,
        h=0.190,
        title="WP80 fake sees the tail",
        sees="blue beyond the 80%-signal cut",
        body=f"The same hard blue marble still passes the cut, so fake rate stays {100*old_metrics['wp80_fake']:.0f}% to {100*new_metrics['wp80_fake']:.0f}%.",
        accent=RED,
        face=SOFT_RED,
    )

    add_box(fig, 0.055, 0.065, 0.890, 0.105, SOFT_GOLD, edge="#E7C46A", lw=1.15, radius=0.012)
    add_text(fig, 0.075, 0.143, "Real motivating case: fixed Jet12+20 validation", size=14.6, weight="bold", color=NAVY)
    add_text(
        fig,
        0.075,
        0.111,
        "Jet12+20 training to Jet12+20+30+40 training: AUC 0.830 to 0.835; logloss 0.562 to 0.433; median gap 0.27 to 0.42; WP80 fake stays about 31%.",
        size=12.7,
        color=INK,
    )
    add_text(
        fig,
        0.075,
        0.083,
        "Lesson: small AUC/fake-rate movement can coexist with a real cleanup of the score scale.",
        size=13.0,
        weight="bold",
        color="#7A4B00",
    )

    fig.savefig(OUTPNG, dpi=SLIDE_DPI)
    plt.close(fig)
    write_outputs(old_m, new_m)
    print(OUTPNG)
    print(OUTNOTE)
    print(OUTSCRIPT)
    print(OUTCSV)
    print(OUTMANIFEST)


if __name__ == "__main__":
    main()
