#!/usr/bin/env python3
"""Build the THE-79 training-window problem slide candidate."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX, slide_figsize  # noqa: E402


OUT_DIR = REPO / "dataOutput/auauBDTModelSelection/THE79_training_window_problem_20260627"
OUT_PNG = OUT_DIR / "the79_training_window_problem_slide.png"
OUT_CSV = OUT_DIR / "the79_training_window_problem_metrics.csv"
OUT_MANIFEST = OUT_DIR / "the79_training_window_problem_manifest.json"
OUT_SCRIPT = OUT_DIR / "the79_training_window_problem_speaker_script.md"
OUT_LAYOUT = OUT_DIR / "the79_training_window_problem_layout_nodes.json"
OUT_AUDIT = OUT_DIR / "the79_training_window_problem_slide_audit.json"

WINDOW_CSV = (
    REPO
    / "dataOutput/auauTightBDTValidation/model_window_comparisons/"
    / "20260625_0_20_score_separation/auau_bdt_score_separation_0_20_pt5to40_vs_pt15to35.csv"
)
PT5_SCORE_DIR = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_cent3x3_20260510_221419/score_caches"
PT15_SCORE_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE57_baseline_validation_20260615/"
    / "a2_nocut_score_caches/score_caches"
)
PT5_REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_cent3x3_20260510_221419"
PT15_REPORT = REPO / "dataOutput/auauTightBDTValidation/THE57_baseline_validation_20260615/a2_nocut"

PT_BINS = [(6.0, 10.0), (10.0, 15.0), (15.0, 20.0), (20.0, 25.0), (25.0, 35.0)]
PT_LABELS = ["6-10", "10-15", "15-20", "20-25", "25-35"]
CENT_RANGE = (0.0, 20.0)

PT5_KEY = "score_centAsFeat3x3_pt5to40"
PT15_KEY = "score_centAsFeatBase3x3_pt15to35"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]

INK = "#0F172A"
MUTED = "#475569"
GRID = "#D9E2EF"
EDGE = "#C9D4E5"
BLUE = "#1F77B4"
ORANGE = "#E66100"
RED = "#B91C1C"
GREEN = "#15803D"
PANEL = "#FBFDFF"
SOFT_BLUE = "#EEF5FC"


@dataclass(frozen=True)
class ModelSpec:
    key: str
    label: str
    score_dir: Path
    score_key: str
    color: str
    product: str
    report: Path


MODELS = [
    ModelSpec(
        key="pt5to40",
        label="5-40 training",
        score_dir=PT5_SCORE_DIR,
        score_key=PT5_KEY,
        color=BLUE,
        product="centAsFeat3x3_pt5to40",
        report=PT5_REPORT,
    ),
    ModelSpec(
        key="pt15to35",
        label="15-35 default",
        score_dir=PT15_SCORE_DIR,
        score_key=PT15_KEY,
        color=ORANGE,
        product="centAsFeatBase3x3_pt15to35",
        report=PT15_REPORT,
    ),
]


def setup_style() -> None:
    for font in TIMES_FONTS:
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": PANEL,
            "axes.edgecolor": EDGE,
            "axes.linewidth": 1.0,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.unicode_minus": False,
        }
    )


def font_px(points: float) -> float:
    return points * SLIDE_DPI / 72.0


def text_bbox_px(fig: plt.Figure, artist) -> list[float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(renderer=fig.canvas.get_renderer())
    scale_x = fig.bbox.width / SLIDE_WIDTH_PX
    scale_y = fig.bbox.height / SLIDE_HEIGHT_PX
    height = fig.bbox.height
    return [
        float(bbox.x0 / scale_x),
        float((height - bbox.y1) / scale_y),
        float(bbox.x1 / scale_x),
        float((height - bbox.y0) / scale_y),
    ]


def add_text(
    fig: plt.Figure,
    nodes: list[dict],
    name: str,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    role: str = "audience",
    weight: str = "normal",
    color: str = INK,
    ha: str = "left",
    va: str = "center",
    linespacing: float = 1.0,
    **flags,
) -> None:
    artist = fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
    )
    nodes.append(
        {
            "name": name,
            "kind": "text",
            "role": role,
            "text": text,
            "font_px": font_px(size),
            "bbox": text_bbox_px(fig, artist),
            "title_anchor": role == "title",
            **flags,
        }
    )


def rounded_box(
    fig: plt.Figure,
    nodes: list[dict],
    name: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str,
    edge: str,
    lw: float = 1.3,
    radius: float = 0.006,
    role: str = "card",
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.006,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=3,
        )
    )
    nodes.append(
        {
            "name": name,
            "kind": "rect",
            "role": role,
            "bbox": [
                x * SLIDE_WIDTH_PX,
                (1.0 - y - h) * SLIDE_HEIGHT_PX,
                (x + w) * SLIDE_WIDTH_PX,
                (1.0 - y) * SLIDE_HEIGHT_PX,
            ],
        }
    )


def auc_rank(y_true: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(score, dtype="float64")
    finite = np.isfinite(s)
    y = y[finite]
    s = s[finite]
    n_sig = int(y.sum())
    n_bkg = int((~y).sum())
    if n_sig == 0 or n_bkg == 0:
        return math.nan
    order = np.argsort(s, kind="mergesort")
    sorted_s = s[order]
    ranks = np.empty(len(s), dtype="float64")
    i = 0
    while i < len(s):
        j = i + 1
        while j < len(s) and sorted_s[j] == sorted_s[i]:
            j += 1
        ranks[order[i:j]] = (i + 1 + j) / 2.0
        i = j
    return float((ranks[y].sum() - n_sig * (n_sig + 1) / 2.0) / (n_sig * n_bkg))


def wp80_background_fake_rate(y_true: np.ndarray, score: np.ndarray) -> tuple[float, float]:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(score, dtype="float64")
    finite = np.isfinite(s)
    y = y[finite]
    s = s[finite]
    sig = s[y]
    bkg = s[~y]
    if len(sig) == 0 or len(bkg) == 0:
        return math.nan, math.nan
    threshold = float(np.quantile(sig, 0.20, method="lower"))
    fake = float(np.mean(bkg >= threshold))
    return threshold, fake


def score_separation(y_true: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(score, dtype="float64")
    finite = np.isfinite(s)
    y = y[finite]
    s = s[finite]
    sig = s[y]
    bkg = s[~y]
    if len(sig) == 0 or len(bkg) == 0:
        return math.nan
    pooled = math.sqrt(0.5 * (float(np.var(sig)) + float(np.var(bkg))))
    if pooled <= 0:
        return math.nan
    return float((float(np.mean(sig)) - float(np.mean(bkg))) / pooled)


def binary_logloss(y_true: np.ndarray, score: np.ndarray, *, balanced: bool = False) -> float:
    y = np.asarray(y_true, dtype="float64")
    s = np.asarray(score, dtype="float64")
    finite = np.isfinite(y) & np.isfinite(s)
    y = y[finite]
    s = s[finite]
    if len(y) == 0:
        return math.nan
    eps = 1.0e-6
    p = np.clip(s, eps, 1.0 - eps)
    loss = -(y * np.log(p) + (1.0 - y) * np.log(1.0 - p))
    if not balanced:
        return float(np.mean(loss))
    n_sig = float(np.sum(y == 1.0))
    n_bkg = float(np.sum(y == 0.0))
    if n_sig == 0.0 or n_bkg == 0.0:
        return math.nan
    weights = np.where(y == 1.0, 0.5 / n_sig, 0.5 / n_bkg)
    return float(np.sum(weights * loss))


def load_model_arrays(spec: ModelSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ys: list[np.ndarray] = []
    ets: list[np.ndarray] = []
    scores: list[np.ndarray] = []
    files = sorted(spec.score_dir.glob("score_cache_*.npz"))
    if not files:
        raise FileNotFoundError(f"No score caches found under {spec.score_dir}")
    for path in files:
        with np.load(path, allow_pickle=False) as cache:
            y = np.asarray(cache["is_signal"], dtype=bool)
            et = np.asarray(cache["cluster_Et"], dtype="float64")
            cent = np.asarray(cache["centrality"], dtype="float64")
            score = np.asarray(cache[spec.score_key], dtype="float64")
        mask = (
            (cent >= CENT_RANGE[0])
            & (cent < CENT_RANGE[1])
            & np.isfinite(et)
            & np.isfinite(score)
        )
        if mask.any():
            ys.append(y[mask])
            ets.append(et[mask])
            scores.append(score[mask])
    if not ys:
        return np.array([], dtype=bool), np.array([], dtype="float64"), np.array([], dtype="float64")
    return np.concatenate(ys), np.concatenate(ets), np.concatenate(scores)


def compute_metrics() -> list[dict]:
    rows: list[dict] = []
    source_counts = {}
    for spec in MODELS:
        y, et, score = load_model_arrays(spec)
        source_counts[spec.key] = int(len(score))
        for (lo, hi), label in zip(PT_BINS, PT_LABELS):
            mask = (et >= lo) & (et < hi) & np.isfinite(score)
            yy = y[mask]
            ss = score[mask]
            threshold, fake = wp80_background_fake_rate(yy, ss)
            row = {
                "model_key": spec.key,
                "model_label": spec.label,
                "product": spec.product,
                "pt_bin": label,
                "pt_low_GeV": lo,
                "pt_high_GeV": hi,
                "entries": int(len(yy)),
                "signal_entries": int(yy.sum()) if len(yy) else 0,
                "background_entries": int((~yy).sum()) if len(yy) else 0,
                "auc": auc_rank(yy, ss) if len(yy) else math.nan,
                "wp80_threshold": threshold,
                "wp80_background_pass_rate": fake,
                "score_separation": score_separation(yy, ss),
                "score_logloss": binary_logloss(yy, ss),
                "balanced_score_logloss": binary_logloss(yy, ss, balanced=True),
                "score_min": float(np.nanmin(ss)) if len(ss) else math.nan,
                "score_max": float(np.nanmax(ss)) if len(ss) else math.nan,
            }
            rows.append(row)

    lookup = {(row["model_key"], row["pt_bin"]): row for row in rows}
    for label in PT_LABELS:
        r5 = lookup[("pt5to40", label)]
        r15 = lookup[("pt15to35", label)]
        if math.isfinite(float(r5["auc"])) and math.isfinite(float(r15["auc"])):
            r5["delta_auc_vs_15to35"] = float(r5["auc"]) - float(r15["auc"])
            r5["delta_wp80_background_pass_rate_vs_15to35_pp"] = (
                float(r5["wp80_background_pass_rate"]) - float(r15["wp80_background_pass_rate"])
            ) * 100.0
            r5["delta_score_separation_vs_15to35"] = float(r5["score_separation"]) - float(r15["score_separation"])
            r5["delta_score_logloss_vs_15to35"] = float(r5["score_logloss"]) - float(r15["score_logloss"])
            r5["delta_balanced_score_logloss_vs_15to35"] = (
                float(r5["balanced_score_logloss"]) - float(r15["balanced_score_logloss"])
            )
        else:
            r5["delta_auc_vs_15to35"] = math.nan
            r5["delta_wp80_background_pass_rate_vs_15to35_pp"] = math.nan
            r5["delta_score_separation_vs_15to35"] = math.nan
            r5["delta_score_logloss_vs_15to35"] = math.nan
            r5["delta_balanced_score_logloss_vs_15to35"] = math.nan
        r15["delta_auc_vs_15to35"] = 0.0 if math.isfinite(float(r15["auc"])) else math.nan
        r15["delta_wp80_background_pass_rate_vs_15to35_pp"] = 0.0 if math.isfinite(float(r15["auc"])) else math.nan
        r15["delta_score_separation_vs_15to35"] = 0.0 if math.isfinite(float(r15["auc"])) else math.nan
        r15["delta_score_logloss_vs_15to35"] = 0.0 if math.isfinite(float(r15["score_logloss"])) else math.nan
        r15["delta_balanced_score_logloss_vs_15to35"] = (
            0.0 if math.isfinite(float(r15["balanced_score_logloss"])) else math.nan
        )

    return rows


def write_metrics(rows: list[dict]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = [
        "model_key",
        "model_label",
        "product",
        "pt_bin",
        "pt_low_GeV",
        "pt_high_GeV",
        "entries",
        "signal_entries",
        "background_entries",
        "auc",
        "wp80_threshold",
        "wp80_background_pass_rate",
        "score_separation",
        "score_logloss",
        "balanced_score_logloss",
        "score_min",
        "score_max",
        "delta_auc_vs_15to35",
        "delta_wp80_background_pass_rate_vs_15to35_pp",
        "delta_score_separation_vs_15to35",
        "delta_score_logloss_vs_15to35",
        "delta_balanced_score_logloss_vs_15to35",
    ]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def metric_vectors(rows: list[dict]) -> dict[str, np.ndarray]:
    lookup = {(row["model_key"], row["pt_bin"]): row for row in rows}

    def values(model: str, field: str) -> np.ndarray:
        return np.asarray([float(lookup[(model, label)][field]) for label in PT_LABELS], dtype="float64")

    return {
        "auc5": values("pt5to40", "auc"),
        "auc15": values("pt15to35", "auc"),
        "wp5": values("pt5to40", "wp80_background_pass_rate") * 100.0,
        "wp15": values("pt15to35", "wp80_background_pass_rate") * 100.0,
        "dprime5": values("pt5to40", "score_separation"),
        "dprime15": values("pt15to35", "score_separation"),
        "logloss5": values("pt5to40", "balanced_score_logloss"),
        "logloss15": values("pt15to35", "balanced_score_logloss"),
    }


def bar_label(ax: plt.Axes, x: float, y: float, text: str, *, color: str, dy: float = 0.0, size: float = 12.5) -> None:
    ax.text(x, y + dy, text, ha="center", va="bottom" if dy >= 0 else "top", fontsize=size, color=color, fontweight="bold")


def decorate_axis(ax: plt.Axes, *, title: str, y_label: str | None = None, x_label: str | None = None) -> None:
    ax.set_title(title, loc="left", fontsize=17.5, fontweight="bold", color=INK, pad=8)
    if y_label:
        ax.set_ylabel(y_label, fontsize=14.0, color=INK, labelpad=8)
    if x_label:
        ax.set_xlabel(x_label, fontsize=16.0, fontweight="bold", color=INK, labelpad=8)
    ax.grid(axis="y", color=GRID, linewidth=1.1, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(EDGE)
    ax.tick_params(labelsize=13.5, colors=INK)


def draw_auc_panel(fig: plt.Figure, vectors: dict[str, np.ndarray]) -> None:
    ax = fig.add_axes([0.083, 0.445, 0.835, 0.225])
    x = np.arange(len(PT_LABELS))
    for idx in range(len(PT_LABELS)):
        if idx % 2 == 0:
            ax.axvspan(idx - 0.5, idx + 0.5, color=SOFT_BLUE, alpha=0.55, zorder=0)
    ax.plot(
        x,
        vectors["auc5"],
        color=BLUE,
        marker="o",
        markersize=9.0,
        linewidth=3.0,
        zorder=4,
    )
    finite15 = np.isfinite(vectors["auc15"])
    ax.plot(
        x[finite15],
        vectors["auc15"][finite15],
        color=ORANGE,
        marker="s",
        markersize=8.5,
        linewidth=3.0,
        zorder=5,
    )
    for idx, value in enumerate(vectors["auc5"]):
        ax.text(idx, value + 0.012, f"{value:.3f}", ha="center", va="bottom", fontsize=12.5, color=BLUE, fontweight="bold")
    for idx, value in enumerate(vectors["auc15"]):
        if np.isfinite(value):
            ax.text(idx, value + 0.012, f"{value:.3f}", ha="center", va="bottom", fontsize=12.5, color=ORANGE, fontweight="bold")
    ax.axhline(0.5, color="#94A3B8", linewidth=1.1, linestyle=(0, (4, 4)), zorder=1)
    ax.set_xticks(x, [])
    ax.set_ylim(0.50, 0.90)
    decorate_axis(ax, title="ROC AUC by cluster $E_T$ bin", y_label="AUC")
    ax.legend(
        handles=[
            Line2D([0], [0], color=BLUE, marker="o", markersize=8.5, linewidth=3.0, label="5-40 training"),
            Line2D([0], [0], color=ORANGE, marker="s", markersize=8.0, linewidth=3.0, label="15-35 default"),
        ],
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.985, 1.18),
        fontsize=13.5,
        ncol=2,
        handlelength=1.1,
        columnspacing=1.1,
    )


def draw_logloss_panel(fig: plt.Figure, vectors: dict[str, np.ndarray]) -> None:
    ax = fig.add_axes([0.083, 0.165, 0.835, 0.195])
    x = np.arange(len(PT_LABELS))
    for idx in range(len(PT_LABELS)):
        if idx % 2 == 0:
            ax.axvspan(idx - 0.5, idx + 0.5, color=SOFT_BLUE, alpha=0.55, zorder=0)
    ax.plot(
        x,
        vectors["logloss5"],
        color=BLUE,
        marker="o",
        markersize=9.0,
        linewidth=3.0,
        zorder=4,
    )
    finite15 = np.isfinite(vectors["logloss15"])
    ax.plot(
        x[finite15],
        vectors["logloss15"][finite15],
        color=ORANGE,
        marker="s",
        markersize=8.5,
        linewidth=3.0,
        zorder=5,
    )
    for idx, value in enumerate(vectors["logloss5"]):
        if np.isfinite(value):
            ax.text(idx, value + 0.035, f"{value:.2f}", ha="center", va="bottom", fontsize=12.5, color=BLUE, fontweight="bold")
    for idx, value in enumerate(vectors["logloss15"]):
        if np.isfinite(value):
            ax.text(idx, value - 0.035, f"{value:.2f}", ha="center", va="top", fontsize=12.5, color=ORANGE, fontweight="bold")
    ax.set_xticks(x, PT_LABELS, fontsize=14.5, fontweight="bold")
    finite_values = np.concatenate([vectors["logloss5"][np.isfinite(vectors["logloss5"])], vectors["logloss15"][np.isfinite(vectors["logloss15"])]])
    y_min = max(0.0, float(np.nanmin(finite_values)) - 0.12)
    y_max = float(np.nanmax(finite_values)) + 0.18
    ax.set_ylim(y_min, y_max)
    decorate_axis(
        ax,
        title="Balanced score logloss by cluster $E_T$ bin",
        y_label="logloss",
        x_label="cluster $E_T$ bin (GeV)",
    )


def add_header(fig: plt.Figure, nodes: list[dict]) -> None:
    add_text(
        fig,
        nodes,
        "title",
        0.050,
        0.943,
        r"Broad 5-40 training loses high-$E_T$ separation",
        size=30.0,
        role="title",
        weight="bold",
        va="top",
    )
    rounded_box(fig, nodes, "problem_callout_box", 0.058, 0.758, 0.560, 0.105, face="white", edge="#EF4444", lw=1.6)
    add_text(
        fig,
        nodes,
        "problem_callout_lead",
        0.076,
        0.836,
        "SIMPLE WINDOW TEST",
        size=15.4,
        role="audience",
        weight="bold",
        color="#B91C1C",
    )
    add_text(
        fig,
        nodes,
        "problem_callout_scope",
        0.076,
        0.806,
        "0-20% centrality; embedded signal vs inclusive background",
        size=15.0,
        role="audience",
        weight="bold",
        color="#B91C1C",
    )
    add_text(
        fig,
        nodes,
        "problem_callout_body",
        0.076,
        0.778,
        r"5-40 covers low $E_T$; high-$E_T$ AUC flattens and logloss rises.",
        size=15.0,
        role="audience",
        weight="bold",
        color="#B91C1C",
    )

def add_readout_card(fig: plt.Figure, nodes: list[dict], rows: list[dict]) -> None:
    lookup = {(row["model_key"], row["pt_bin"]): row for row in rows}
    r5 = lookup[("pt5to40", "25-35")]
    r15 = lookup[("pt15to35", "25-35")]
    auc_text = f"{float(r15['auc']):.3f} → {float(r5['auc']):.3f}"
    logloss_text = f"{float(r15['balanced_score_logloss']):.2f} → {float(r5['balanced_score_logloss']):.2f}"
    rounded_box(fig, nodes, "readout_card", 0.650, 0.758, 0.268, 0.105, face="#FFFFFF", edge="#CBD5E1", lw=1.2)
    add_text(
        fig,
        nodes,
        "readout_card_title",
        0.668,
        0.833,
        "25-35 GeV readout",
        size=15.0,
        role="audience",
        color=INK,
        weight="bold",
        title_band_exception=True,
    )
    add_text(
        fig,
        nodes,
        "readout_card_auc",
        0.668,
        0.801,
        f"AUC {auc_text}",
        size=16.2,
        role="audience",
        color=RED,
        weight="bold",
        title_band_exception=True,
    )
    add_text(
        fig,
        nodes,
        "readout_card_wp80",
        0.668,
        0.771,
        f"balanced logloss {logloss_text}",
        size=16.2,
        role="audience",
        color=RED,
        weight="bold",
        title_band_exception=True,
    )


def add_conclusion(fig: plt.Figure, nodes: list[dict]) -> None:
    add_text(
        fig,
        nodes,
        "conclusion",
        0.500,
        0.055,
        "Use routed or binned 5-40 models; one broad global model is not enough.",
        size=17.0,
        role="audience",
        color=INK,
        weight="bold",
        ha="center",
    )


def write_manifest(rows: list[dict], audit_status: str | None) -> None:
    payload = {
        "schema": "THE79_TRAINING_WINDOW_PROBLEM_SLIDE_V1",
        "png": str(OUT_PNG),
        "metrics_csv": str(OUT_CSV),
        "speaker_script": str(OUT_SCRIPT),
        "layout_nodes": str(OUT_LAYOUT),
        "audit_json": str(OUT_AUDIT),
        "audit_status": audit_status,
        "visible_scope": {
            "centrality": "0-20%",
            "comparison": "5-40 training versus 15-35 default/reference training",
            "signal": "embedded photon signal",
            "background": "embedded inclusive-jet background",
            "finite_scores_only": True,
        },
        "source_inputs": {
            "existing_window_comparison_csv": str(WINDOW_CSV),
            "pt5to40_score_caches": str(PT5_SCORE_DIR),
            "pt15to35_score_caches": str(PT15_SCORE_DIR),
            "pt5to40_report": str(PT5_REPORT),
            "pt15to35_report": str(PT15_REPORT),
        },
        "model_products": {
            "pt5to40": "centAsFeat3x3_pt5to40",
            "pt15to35": "centAsFeatBase3x3_pt15to35",
        },
        "metrics": {
            "auc": "rank-based ROC AUC computed from local score caches in 0-20 centrality and pT bin",
            "wp80_background_pass_rate": "background fraction passing threshold chosen for 80% signal efficiency in the same 0-20 centrality and pT bin",
            "score_separation": "mean signal-background score gap divided by pooled RMS; written to CSV, not shown on canvas",
            "score_logloss": "raw binary logloss computed from clipped score-cache outputs in [1e-6, 1 - 1e-6]",
            "balanced_score_logloss": "same clipped binary logloss with signal and background each contributing half the bin weight; lower is better",
        },
        "rows": rows,
    }
    OUT_MANIFEST.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_script(rows: list[dict]) -> None:
    lookup = {(row["model_key"], row["pt_bin"]): row for row in rows}
    high = ["15-20", "20-25", "25-35"]
    auc_drops = [
        1000.0 * (float(lookup[("pt5to40", pt)]["auc"]) - float(lookup[("pt15to35", pt)]["auc"]))
        for pt in high
    ]
    wp_increases = [
        100.0
        * (
            float(lookup[("pt5to40", pt)]["wp80_background_pass_rate"])
            - float(lookup[("pt15to35", pt)]["wp80_background_pass_rate"])
        )
        for pt in high
    ]
    logloss_changes = [
        float(lookup[("pt5to40", pt)]["balanced_score_logloss"])
        - float(lookup[("pt15to35", pt)]["balanced_score_logloss"])
        for pt in high
    ]
    OUT_SCRIPT.write_text(
        "\n".join(
            [
                "# Training-window problem slide",
                "",
                "This slide is the clean motivation for the expanded BDT campaign.",
                "The comparison is deliberately simple: 0-20% centrality, embedded photon signal against embedded inclusive-jet background, using local finite-score validation rows.",
                "",
                "The blue model is the broad 5-40 GeV training window. It gives finite scores below 15 GeV, which is the reason we were tempted by it.",
                "But in the high-ET overlap, the AUC trace flattens downward relative to the 15-35 default reference.",
                "",
                f"Across 15-20, 20-25, and 25-35 GeV, the AUC changes are {auc_drops[0]:.0f}, {auc_drops[1]:.0f}, and {auc_drops[2]:.0f} in units of 1e-3.",
                f"The balanced score-logloss changes are {logloss_changes[0]:+.2f}, {logloss_changes[1]:+.2f}, and {logloss_changes[2]:+.2f}; lower logloss is better.",
                f"The WP80 background pass-rate increases by {wp_increases[0]:.1f}, {wp_increases[1]:.1f}, and {wp_increases[2]:.1f} percentage points.",
                "",
                "The conclusion is not that 5-40 is impossible. The conclusion is that one broad global BDT is the wrong 5-40 solution.",
                "This is the direct motivation for testing routed or binned models: keep the low-ET coverage while recovering high-ET separation.",
                "",
                "The logloss panel uses the cached score outputs directly after clipping to [1e-6, 1 - 1e-6], with signal and background balanced inside each bin. It is a score-calibration diagnostic, while AUC remains the primary ranking metric.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def run_audit() -> str:
    cmd = [
        sys.executable,
        str(SCRIPTS / "slides/common/post_render_slide_audit.py"),
        "--png",
        str(OUT_PNG),
        "--layout-nodes",
        str(OUT_LAYOUT),
        "--output",
        str(OUT_AUDIT),
    ]
    proc = subprocess.run(cmd, cwd=REPO, text=True, capture_output=True)
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise SystemExit(proc.returncode)
    return "passed"


def build_slide(rows: list[dict]) -> list[dict]:
    setup_style()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[dict] = [
        {
            "name": "slide_canvas",
            "kind": "canvas",
            "role": "canvas",
            "bbox": [0, 0, SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
            "width_px": SLIDE_WIDTH_PX,
            "height_px": SLIDE_HEIGHT_PX,
        },
        {
            "name": "title_axis",
            "kind": "axis",
            "role": "layout",
            "title_axis_x": 0.050 * SLIDE_WIDTH_PX,
        },
    ]
    add_header(fig, nodes)
    add_readout_card(fig, nodes, rows)
    vectors = metric_vectors(rows)
    draw_auc_panel(fig, vectors)
    draw_logloss_panel(fig, vectors)
    add_conclusion(fig, nodes)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=SLIDE_DPI)
    plt.close(fig)

    layout_payload = {
        "nodes": nodes,
        "title_axis_x": 0.050 * SLIDE_WIDTH_PX,
    }
    OUT_LAYOUT.write_text(json.dumps(layout_payload, indent=2), encoding="utf-8")
    return nodes


def main() -> None:
    rows = compute_metrics()
    write_metrics(rows)
    build_slide(rows)
    audit_status = run_audit()
    write_manifest(rows, audit_status)
    write_script(rows)
    print(OUT_PNG)
    print(OUT_SCRIPT)
    print(OUT_CSV)
    print(OUT_MANIFEST)
    print(OUT_AUDIT)


if __name__ == "__main__":
    main()
